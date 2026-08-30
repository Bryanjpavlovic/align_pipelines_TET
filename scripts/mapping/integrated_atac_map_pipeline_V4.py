#!/usr/bin/env python3
"""
Integrated Mapping Pipeline for 10X Multiome ATAC-seq Data - V3

Changes from V2:
  - Symlinks now include run tags: LibName_S##_L##__RunTag__R#_001.fastq.gz
    Run numbering is per-library (Run001, Run002...), ordered by folder age
  - Collision-safe: no symlink can ever overwrite another, script aborts if it would
  - Manifest file records original path -> symlinked name for every file
  - --libraries flag to process a subset by number (e.g. --libraries 3 5 12 34)
  - --consolidate mode called from sbatch for run-tagged symlinking
  - All SLURM time limits set to 24 days (576 hours)

Usage:
  python3 integrated_atac_map_pipeline_V3.py \
    --input-dirs /mnt/beegfs/reads/10X_ATAC_multiome/BP12957 \
                 /mnt/beegfs/reads/10X_ATAC_multiome/BP14536B \
    --output-base /mnt/beegfs/tetmultiome_atac/ \
    --libraries 3 5 12 34 \
    --submit-jobs
"""

import argparse
import glob
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path
from datetime import datetime
from collections import defaultdict
import re


RELEASE = "2026-08-30-v9-align-repo-migration"
ALIGN_PIPELINES_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_ATAC_WORKFLOW = str(ALIGN_PIPELINES_ROOT / "workflows" / "align_atac.nf")
DEFAULT_ATAC_LIB_PREFIX = "Tet_2025_Multiome-ATAC_"


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def get_file_size(filepath):
    """Get file size in bytes, following symlinks."""
    try:
        return os.path.getsize(filepath)
    except OSError:
        return 0


def find_fastq_triplets(input_dir):
    """Find all R1/R2/R3 FASTQ file triplets in the input directory."""
    r3_patterns = [
        os.path.join(input_dir, "*_R3_*.fastq.gz"),
        os.path.join(input_dir, "*_R3_*.fq.gz"),
    ]
    
    triplets = []
    for pattern in r3_patterns:
        for r3_file in glob.glob(pattern):
            r1_file = r3_file.replace('_R3_', '_R1_')
            r2_file = r3_file.replace('_R3_', '_R2_')
            if os.path.exists(r1_file) and os.path.exists(r2_file):
                triplets.append((r1_file, r2_file, r3_file))
    
    return sorted(list(set(triplets)))


def extract_library_name(fastq_file, prefix=DEFAULT_ATAC_LIB_PREFIX):
    """Extract a canonical library name from legacy or annotated FASTQs."""
    basename = os.path.basename(fastq_file)
    match = re.match(r'(.+?)_S\d+_L\d+', basename)
    if match:
        candidate = match.group(1)
    else:
        candidate = re.sub(r'_R[123]_.*', '', basename)
    if prefix:
        prefixed = re.match(rf'({re.escape(prefix)}\d+)(?:_|$)', candidate)
        if prefixed:
            return prefixed.group(1)
    return re.sub(r'(?:_L\d+)+$', '', candidate)


def insert_run_tag(basename, run_tag):
    """
    Insert __RunTag between _S##_L## and _R#_ parts of a FASTQ filename.
    
    Tet_2025_Multiome-ATAC_3_S3_L002_R3_001.fastq.gz
    -> Tet_2025_Multiome-ATAC_3_S3_L002__Run001_R3_001.fastq.gz
    
    The NF regex captures everything before _R3 as match[2], so match[2]
    will be: Tet_2025_Multiome-ATAC_3_S3_L002__Run001
    Then strip_run_tag removes __Run001 to get the original prefix back.
    """
    m = re.match(r'^(.+_S\d+_L\d+)(_R[123]_.+)$', basename)
    if m:
        return f"{m.group(1)}__{run_tag}{m.group(2)}"
    # Fallback: insert before _R#_
    m = re.match(r'^(.+?)(_R[123]_.+)$', basename)
    if m:
        return f"{m.group(1)}__{run_tag}{m.group(2)}"
    # Last resort
    return f"{basename}__{run_tag}"


def get_folder_sort_key(folder_path):
    """Sort key for input folders: by modification time (oldest first), name as tiebreaker."""
    try:
        mtime = os.path.getmtime(folder_path)
    except OSError:
        mtime = 0
    return (mtime, os.path.basename(folder_path))


def extract_flowcell_lane(fastq_path):
    """
    Extract flowcell ID and lane number from the first read header of a FASTQ.
    
    Read header format: @InstrumentID:RunNum:FlowcellID:Lane:...
    e.g. @A01535:265:HW5MYDRX3:2:1101:5327:1000 1:N:0:ACGTATCA
    
    Returns (flowcell, lane_num) e.g. ('HW5MYDRX3', '2')
    Returns ('unknown', '0') if extraction fails.
    """
    import gzip
    try:
        with gzip.open(fastq_path, 'rt') as f:
            header = f.readline().strip()
        fields = header.split(':')
        if len(fields) >= 4:
            flowcell = fields[2]
            lane_num = fields[3]
            return (flowcell, lane_num)
    except Exception as e:
        print(f"    WARNING: Could not extract flowcell/lane from {fastq_path}: {e}")
    return ('unknown', '0')


def extract_sample_lane(fastq_basename):
    """
    Extract S# and L### from a FASTQ basename.
    e.g. Tet_2025_Multiome-ATAC_3_S3_L002_R1_001.fastq.gz -> ('S3', 'L002')
    """
    m = re.search(r'_(S\d+)_(L\d+)', fastq_basename)
    if m:
        return (m.group(1), m.group(2))
    return ('S0', 'L000')


def library_nums_to_names(lib_nums, prefix=DEFAULT_ATAC_LIB_PREFIX):
    """Convert library numbers like [3, 5, 12] to full names."""
    return {f"{prefix}{n}" for n in lib_nums}


# ============================================================================
# CONSOLIDATION
# ============================================================================

def consolidate_inputs(input_dirs, consolidated_dir, project_dir,
                       filter_lib_names=None, lib_prefix=DEFAULT_ATAC_LIB_PREFIX):
    """
    Create run-tagged symlinks for all FASTQ triplets.
    
    Every file gets a run tag inserted: __Run001__, __Run002__, etc.
    Run numbering is per-library, ordered by source folder age (oldest first),
    with name sort as tiebreaker for same-age folders.
    
    Reuses only an existing symlink to the exact same source and aborts if any
    target would collide or if unrelated FASTQs are already staged.
    Writes manifest and run assignment files for full traceability.
    """
    os.makedirs(consolidated_dir, exist_ok=True)
    os.makedirs(project_dir, exist_ok=True)
    
    # -------------------------------------------------------------------------
    # Step 1: Scan all input dirs, group triplets by library
    # -------------------------------------------------------------------------
    # Sort input dirs by age (oldest first), name as tiebreaker
    sorted_dirs = sorted(input_dirs, key=get_folder_sort_key)
    
    # lib_name -> [(folder_path, r1, r2, r3), ...]
    lib_folder_triplets = defaultdict(list)
    
    for input_dir in sorted_dirs:
        if not os.path.isdir(input_dir):
            print(f"  WARNING: {input_dir} not found, skipping")
            continue
        
        run_id = os.path.basename(input_dir)
        triplets = find_fastq_triplets(input_dir)
        
        if not triplets:
            print(f"  WARNING: No ATAC triplets found in {run_id}")
            continue
        
        print(f"  Scanned {run_id}: {len(triplets)} triplets")
        
        for r1, r2, r3 in triplets:
            lib_name = extract_library_name(r1, lib_prefix)
            if filter_lib_names is not None and lib_name not in filter_lib_names:
                continue
            lib_folder_triplets[lib_name].append((input_dir, r1, r2, r3))
    
    if not lib_folder_triplets:
        print("  ERROR: No libraries found to consolidate!")
        sys.exit(1)
    
    # -------------------------------------------------------------------------
    # Step 2: Assign per-library run numbers
    # -------------------------------------------------------------------------
    # For each library, unique folders in order of appearance (age-sorted).
    # Each folder gets Run001, Run002, etc for that library.
    
    # lib_name -> {folder_path: "Run001"}
    lib_run_tags = {}
    
    for lib_name, entries in lib_folder_triplets.items():
        seen_folders = []
        for folder_path, _, _, _ in entries:
            if folder_path not in seen_folders:
                seen_folders.append(folder_path)
        
        run_map = {}
        for idx, folder_path in enumerate(seen_folders):
            run_map[folder_path] = f"Run{idx + 1:03d}"
        
        lib_run_tags[lib_name] = run_map
    
    # -------------------------------------------------------------------------
    # Step 3: Create symlinks with collision checking + extract RG metadata
    # -------------------------------------------------------------------------
    symlink_registry = {}  # symlink_basename -> original_full_path
    manifest = []
    rg_metadata = []  # (fnbase, library, bp_id, sample_idx, lane, flowcell, flowcell_lane)
    total_triplets = 0
    
    print(f"\n  Creating run-tagged symlinks and extracting RG metadata...")
    
    for lib_name in sorted(lib_folder_triplets.keys()):
        entries = lib_folder_triplets[lib_name]
        run_map = lib_run_tags[lib_name]
        
        lib_triplet_count = 0
        lib_runs_seen = []
        
        for folder_path, r1, r2, r3 in entries:
            run_tag = run_map[folder_path]
            folder_name = os.path.basename(folder_path)
            
            if folder_name not in [x[1] for x in lib_runs_seen]:
                lib_runs_seen.append((run_tag, folder_name))
            
            # Extract flowcell and lane from R1 header (once per triplet)
            flowcell, lane_num = extract_flowcell_lane(r1)
            sample_idx, lane_id = extract_sample_lane(os.path.basename(r1))
            
            # Determine fnbase from the R3 symlink name (what NF will see as match[2])
            r3_orig_base = os.path.basename(r3)
            r3_new_base = insert_run_tag(r3_orig_base, run_tag)
            # fnbase = everything before _R3 in the symlink name
            fnbase_match = re.match(r'^(.+?)_R3', r3_new_base)
            if fnbase_match:
                fnbase = fnbase_match.group(1)
            else:
                fnbase = r3_new_base
            
            # Build the full RG string that NF will use in minimap2 -R
            # Single \t in the file — Groovy interpolation doesn't re-escape,
            # so minimap2 receives literal \t and converts to tab internally
            rg_id = f"{lib_name}_{folder_name}_{sample_idx}_{lane_id}"
            pu = f"{flowcell}.{lane_num}"
            bs = chr(92)  # single backslash
            sep = bs + "t"  # \t in the file
            rg_string = f"@RG{sep}ID:{rg_id}{sep}SM:{lib_name}{sep}PL:Illumina{sep}PU:{pu}{sep}DS:{folder_name}"
            
            rg_metadata.append((fnbase, lib_name, folder_name, sample_idx, lane_id,
                               flowcell, lane_num, pu, rg_id, rg_string))
            
            for src in [r1, r2, r3]:
                orig_base = os.path.basename(src)
                new_base = insert_run_tag(orig_base, run_tag)
                
                # COLLISION CHECK: registry
                if new_base in symlink_registry:
                    existing_src = symlink_registry[new_base]
                    print(f"\n  !! FATAL COLLISION DETECTED !!")
                    print(f"     Symlink name:  {new_base}")
                    print(f"     Already from:  {existing_src}")
                    print(f"     New source:    {src}")
                    print(f"     Aborting to prevent data loss.")
                    sys.exit(1)
                
                dst = os.path.join(consolidated_dir, new_base)
                
                # COLLISION CHECK: filesystem.  A failed downstream Nextflow
                # run may be safely resumed: reuse only an existing symlink
                # that resolves to this exact source.  Everything else is a
                # hard failure, so no file can be silently replaced.
                if os.path.lexists(dst):
                    if os.path.islink(dst) and os.path.realpath(dst) == os.path.realpath(src):
                        pass
                    else:
                        print(f"\n  !! FATAL: Target already exists on disk !!")
                        print(f"     {dst}")
                        print(f"     Existing target is not the expected source: {src}")
                        print(f"     Aborting to prevent data loss.")
                        sys.exit(1)
                else:
                    os.symlink(src, dst)
                symlink_registry[new_base] = src
                manifest.append((src, new_base, lib_name, run_tag, folder_name))
            
            lib_triplet_count += 1
            total_triplets += 1
        
        run_str = ", ".join(f"{t}={f}" for t, f in lib_runs_seen)
        print(f"    {lib_name}: {lib_triplet_count} triplets, {len(run_map)} run(s) [{run_str}]")
    
    # -------------------------------------------------------------------------
    # Step 4: Reject stale inputs and verify counts
    # -------------------------------------------------------------------------
    unexpected_existing = sorted(
        path for path in glob.glob(os.path.join(consolidated_dir, '*.fastq.gz'))
        if os.path.basename(path) not in symlink_registry
    )
    if unexpected_existing:
        print("\n  !! FATAL: Unrelated FASTQs already exist in mapping_input !!")
        for path in unexpected_existing[:5]:
            print(f"     {path}")
        print("     Refusing to mix staged runs; use a new orchestrator --run-name.")
        sys.exit(1)

    num_r1 = len(glob.glob(os.path.join(consolidated_dir, '*_R1_*.fastq.gz')))
    num_r2 = len(glob.glob(os.path.join(consolidated_dir, '*_R2_*.fastq.gz')))
    num_r3 = len(glob.glob(os.path.join(consolidated_dir, '*_R3_*.fastq.gz')))
    
    print(f"\n  Consolidation summary:")
    print(f"    Libraries:      {len(lib_folder_triplets)}")
    print(f"    Total triplets: {total_triplets}")
    print(f"    R1 files:       {num_r1}")
    print(f"    R2 files:       {num_r2}")
    print(f"    R3 files:       {num_r3}")
    print(f"    Total symlinks: {len(symlink_registry)}")
    
    if num_r1 != num_r3 or num_r2 != num_r3:
        print("    ERROR: Mismatch between R1, R2, and R3 file counts!")
        sys.exit(1)
    
    # -------------------------------------------------------------------------
    # Step 5: Write manifest, run assignments, and RG metadata
    # -------------------------------------------------------------------------
    manifest_path = os.path.join(project_dir, 'symlink_manifest.tsv')
    with open(manifest_path, 'w') as f:
        f.write("original_path\tsymlink_name\tlibrary\trun_tag\tsource_folder\n")
        for orig, sym, lib, tag, folder in manifest:
            f.write(f"{orig}\t{sym}\t{lib}\t{tag}\t{folder}\n")
    print(f"    Manifest:        {manifest_path}")
    
    run_summary_path = os.path.join(project_dir, 'run_assignments.tsv')
    with open(run_summary_path, 'w') as f:
        f.write("library\trun_tag\tfolder_name\tfolder_path\n")
        for lib_name in sorted(lib_run_tags.keys()):
            run_map = lib_run_tags[lib_name]
            for folder_path, run_tag in run_map.items():
                f.write(f"{lib_name}\t{run_tag}\t{os.path.basename(folder_path)}\t{folder_path}\n")
    print(f"    Run assignments: {run_summary_path}")
    
    rg_metadata_path = os.path.join(project_dir, 'rg_metadata.tsv')
    with open(rg_metadata_path, 'w') as f:
        f.write("fnbase\tlibrary\tbp_id\tsample_idx\tlane\tflowcell\tlane_num\tpu\trg_id\trg_string\n")
        for entry in rg_metadata:
            f.write('\t'.join(str(x) for x in entry) + '\n')
    print(f"    RG metadata:     {rg_metadata_path}")
    
    return total_triplets


# ============================================================================
# LIBS FILE GENERATION
# ============================================================================

def generate_libs_file(libraries, output_base):
    """Generate libs.txt file with library names."""
    libs_file = os.path.join(output_base, 'libs.txt')
    with open(libs_file, 'w') as f:
        for lib in sorted(libraries):
            f.write(f"{lib}\n")
    return libs_file


# ============================================================================
# SBATCH GENERATION
# ============================================================================

def generate_mapping_script(output_base, input_dirs, libs_file, atac_ref,
                           rna_whitelist, atac_whitelist, memgb, threads,
                           num_chunks, workflow_file, filter_lib_nums=None,
                           lib_prefix=DEFAULT_ATAC_LIB_PREFIX):
    """Generate SLURM script for ATAC mapping via Nextflow."""
    
    project_dir = os.path.join(output_base, 'mapping_project')
    os.makedirs(project_dir, exist_ok=True)
    
    script_name = os.path.join(output_base, 'run_atac_mapping.sbatch')
    consolidated_dir = os.path.join(output_base, 'mapping_input')
    output_dir = os.path.join(output_base, 'mapping_output')
    script_path = os.path.abspath(__file__)
    workflow_file = os.path.abspath(workflow_file)
    with open(workflow_file, 'rb') as workflow_handle:
        workflow_sha256 = hashlib.sha256(workflow_handle.read()).hexdigest()
    input_dirs_str = ' '.join([f'"{d}"' for d in input_dirs])
    
    # Build the consolidate command
    consolidate_cmd = f'python3 {script_path} --consolidate \\\n'
    consolidate_cmd += f'    --input-dirs {input_dirs_str} \\\n'
    consolidate_cmd += f'    --output-base {output_base} \\\n'
    consolidate_cmd += f'    --lib-prefix {lib_prefix}'
    if filter_lib_nums:
        lib_nums_str = ' '.join(str(n) for n in filter_lib_nums)
        consolidate_cmd += f' \\\n    --libraries {lib_nums_str}'
    
    script_content = f'''#!/bin/bash
#SBATCH --job-name=atac_mapping
#SBATCH --output={project_dir}/mapping_%j.out
#SBATCH --error={project_dir}/mapping_%j.err
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=576:00:00
#SBATCH --exclude=squirtle
#SBATCH --chdir={project_dir}

# ==============================================================================
# MULTIOME ATAC-SEQ MAPPING PIPELINE - V3
# ==============================================================================

set -euo pipefail

echo "========================================================================"
echo "MULTIOME ATAC-SEQ MAPPING PIPELINE - V3"
echo "========================================================================"
echo "Started at: $(date)"
echo "Running on: $(hostname)"
echo ""

mkdir -p {project_dir}/logs
mkdir -p {consolidated_dir}
mkdir -p {output_dir}

module purge
module load miniforge/3 nextflow/latest align_pipelines/bjp

echo "ALIGN_PIPELINES_HOME: ${{ALIGN_PIPELINES_HOME:-NOT SET}}"

# The RG-aware ATAC workflow is owned by the align_pipelines/bjp module.
WORKFLOW_SHA256={workflow_sha256!r}
if [[ -z "${{ALIGN_PIPELINES_HOME:-}}" || ! -f "${{ALIGN_PIPELINES_HOME}}/align_pipelines.nf" ]]; then
    echo "ERROR: align_pipelines/bjp did not provide ALIGN_PIPELINES_HOME"
    exit 1
fi
MODULE_WORKFLOW="${{ALIGN_PIPELINES_HOME}}/workflows/align_atac.nf"
if [[ ! -f "${{MODULE_WORKFLOW}}" ]]; then
    echo "ERROR: align_pipelines/bjp ATAC workflow is missing: ${{MODULE_WORKFLOW}}"
    exit 1
fi
if [[ "$(sha256sum "${{MODULE_WORKFLOW}}" | awk '{{print $1}}')" != "${{WORKFLOW_SHA256}}" ]]; then
    echo "WARNING: installed align_pipelines/bjp align_atac.nf differs from the requested workflow; continuing" >&2
fi
if ! command -v atac_fq_preprocess &> /dev/null; then
    echo "ERROR: atac_fq_preprocess not found in PATH"
    exit 1
fi

cd {project_dir}

# ==============================================================================
# STEP 1: CONSOLIDATE INPUT FILES (collision-safe, run-tagged)
# ==============================================================================
echo ""
echo "Step 1: Consolidating input files with run tags..."
echo ""

{consolidate_cmd}

CONSOLIDATE_EXIT=$?
if [ $CONSOLIDATE_EXIT -ne 0 ]; then
    echo "ERROR: Consolidation failed with exit code $CONSOLIDATE_EXIT"
    exit $CONSOLIDATE_EXIT
fi

NUM_R1=$(ls -1 {consolidated_dir}/*_R1_*.fastq.gz 2>/dev/null | wc -l)
NUM_R2=$(ls -1 {consolidated_dir}/*_R2_*.fastq.gz 2>/dev/null | wc -l)
NUM_R3=$(ls -1 {consolidated_dir}/*_R3_*.fastq.gz 2>/dev/null | wc -l)

echo ""
echo "R1: $NUM_R1 | R2: $NUM_R2 | R3: $NUM_R3"

if [ $NUM_R1 -eq 0 ] || [ $NUM_R2 -eq 0 ] || [ $NUM_R3 -eq 0 ]; then
    echo "ERROR: Missing ATAC-seq triplet files!"
    exit 1
fi
if [ $NUM_R1 -ne $NUM_R3 ] || [ $NUM_R2 -ne $NUM_R3 ]; then
    echo "ERROR: Mismatch between R1, R2, and R3 file counts!"
    exit 1
fi

echo ""

# ==============================================================================
# STEP 2: CREATE NEXTFLOW CONFIGURATION
# ==============================================================================
echo "Step 2: Creating Nextflow configuration..."

cat > params_atac.yml << EOF
libs: '{libs_file}'
output_directory: '{output_dir}'
memgb: '{memgb}'
threads: {threads}
num_chunks: {num_chunks}
atac_dir: '{consolidated_dir}'
atac_ref: '{atac_ref}'
rna_whitelist: '{rna_whitelist}'
atac_whitelist: '{atac_whitelist}'
rg_metadata: '{project_dir}/rg_metadata.tsv'
multiome: true
EOF

echo "params_atac.yml:"
cat params_atac.yml
echo ""

cat > nextflow.config << 'NFCONFIG'
process {{
    executor = 'slurm'
    queue = 'compute'
    clusterOptions = '--exclude=squirtle'
    beforeScript = 'module purge && module load miniforge/3 nextflow/latest align_pipelines/bjp && module load htslib/1.20 samtools/1.20 minimap2/2.28 && command -v minimap2 >/dev/null && command -v samtools >/dev/null && command -v bgzip >/dev/null && command -v tabix >/dev/null && command -v sinto >/dev/null'
    
    cpus = 4
    memory = '32 GB'
    time = '576 h'
    errorStrategy = 'finish'
    
    withName: 'preproc_atac_files_multiome' {{
        cpus = 4
        memory = '32 GB'
        time = '576 h'
        errorStrategy = 'retry'
        maxRetries = 3
    }}
    withName: 'split_reads_atac' {{
        cpus = 4
        memory = '32 GB'
        time = '576 h'
        errorStrategy = 'retry'
        maxRetries = 3
    }}
    withName: 'align_atac_files' {{
        cpus = {threads}
        memory = '{memgb} GB'
        time = '576 h'
        errorStrategy = 'retry'
        maxRetries = 3
    }}
    withName: 'cat_atac_bams' {{
        cpus = {threads}
        memory = '{memgb} GB'
        time = '576 h'
        errorStrategy = 'retry'
        maxRetries = 3
    }}
    withName: 'atac_mkdup' {{
        cpus = {threads}
        memory = '100 GB'
        time = '576 h'
        errorStrategy = 'retry'
        maxRetries = 3
    }}
    withName: 'atac_namesort' {{
        cpus = {threads}
        memory = '100 GB'
        time = '576 h'
        errorStrategy = 'retry'
        maxRetries = 3
    }}
    withName: 'atac_fragments' {{
        cpus = {threads}
        memory = '100 GB'
        time = '576 h'
        errorStrategy = 'retry'
        maxRetries = 3
    }}
}}

workDir = "./work"
resume = true

executor {{
    name = 'slurm'
    queueSize = 30
    submitRateLimit = '10/1min'
    pollInterval = '30 sec'
    exitReadTimeout = '10 min'
}}
NFCONFIG

echo "nextflow.config created"
echo ""

# ==============================================================================
# STEP 3: RUN ATAC MAPPING PIPELINE
# ==============================================================================
echo "========================================================================"
echo "Step 3: Running ATAC alignment pipeline"
echo "========================================================================"
echo "Working directory: {project_dir}"
echo "Output directory: {output_dir}"
echo ""
echo "Processing libraries:"
cat {libs_file}
echo ""

rm -rf .nextflow/cache/*/db/LOCK 2>/dev/null || true

nextflow run "${{ALIGN_PIPELINES_HOME}}/align_pipelines.nf" \\
    -params-file params_atac.yml \\
    -resume \\
    -with-report report.html \\
    -with-trace trace.txt \\
    -with-timeline timeline.html \\
    -ansi-log false

NF_EXIT=$?

if [ $NF_EXIT -eq 0 ]; then
    echo ""
    echo "========================================================================"
    echo "SUCCESS: Pipeline completed at $(date)"
    echo "========================================================================"
    echo "Results: {output_dir}"
else
    echo ""
    echo "========================================================================"
    echo "ERROR: Pipeline failed at $(date) with exit code $NF_EXIT"
    echo "========================================================================"
    echo "To resume: resubmit this staged wrapper; it uses ${{ALIGN_PIPELINES_HOME}}/align_pipelines.nf"
    exit $NF_EXIT
fi
'''
    
    with open(script_name, 'w') as f:
        f.write(script_content)
    
    os.chmod(script_name, 0o755)
    return script_name


# ============================================================================
# MAIN
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Integrated ATAC-seq mapping pipeline for 10X Multiome data - V3',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument('--input-dirs', '-i', required=True, nargs='+',
                       help='Input directories containing ATAC FASTQ triplets')
    parser.add_argument('--output-base', '-o', required=True,
                       help='Base output directory')
    parser.add_argument('--libraries', nargs='+', type=int, default=None,
                       help='Library numbers to process (e.g. --libraries 3 5 12 34)')
    parser.add_argument('--lib-prefix', default=DEFAULT_ATAC_LIB_PREFIX,
                       help='Prefix before the numeric ATAC library ID')
    parser.add_argument('--atac-ref',
                       default='/mnt/beegfs/genomes_annotations/ancestral_genomes/litterbox/human_chimp_bonobo/human_chimp_bonobo.mm2',
                       help='Path to minimap2 index')
    parser.add_argument('--rna-whitelist',
                       default='/mnt/beegfs/genomes_annotations/white_lists_adapters/cellranger_10x/RNA-737K-arc-v1.txt.gz',
                       help='Path to RNA barcode whitelist')
    parser.add_argument('--atac-whitelist',
                       default='/mnt/beegfs/genomes_annotations/white_lists_adapters/cellranger_10x/ATAC-737K-arc-v1.txt.gz',
                       help='Path to ATAC barcode whitelist')
    parser.add_argument('--filter-libs', '-f', default=None,
                       help='Path to file with library names to process')
    parser.add_argument('--workflow-file',
                       default=DEFAULT_ATAC_WORKFLOW,
                       help='ATAC workflow in this align_pipelines repository/package')
    parser.add_argument('--memgb', default='100',
                       help='Memory in GB for alignment jobs (default: 100)')
    parser.add_argument('--threads', type=int, default=8,
                       help='Number of threads for mapping (default: 8)')
    parser.add_argument('--num-chunks', type=int, default=4,
                       help='Number of chunks for parallel alignment (default: 4)')
    parser.add_argument('--submit-jobs', action='store_true',
                       help='Automatically submit jobs to SLURM')
    parser.add_argument('--test-run', action='store_true',
                       help='Test mode: only process smallest library')
    parser.add_argument('--consolidate', action='store_true',
                       help='Run consolidation only (called by sbatch)')
    
    args = parser.parse_args()
    args.output_base = args.output_base.rstrip('/')
    
    # Build filter set from --libraries and/or --filter-libs
    filter_lib_names = None
    if args.libraries:
        filter_lib_names = library_nums_to_names(args.libraries, args.lib_prefix)
    if args.filter_libs:
        from_file = set()
        with open(args.filter_libs, 'r') as f:
            from_file = set(line.strip() for line in f if line.strip())
        if filter_lib_names:
            filter_lib_names = filter_lib_names | from_file
        else:
            filter_lib_names = from_file
    
    # Handle --consolidate mode (called from within sbatch)
    if args.consolidate:
        consolidated_dir = os.path.join(args.output_base, 'mapping_input')
        project_dir = os.path.join(args.output_base, 'mapping_project')
        print(f"Consolidating FASTQs from {len(args.input_dirs)} runs into {consolidated_dir}")
        if filter_lib_names:
            print(f"Filtering to libraries: {sorted(filter_lib_names)}")
        consolidate_inputs(
            args.input_dirs,
            consolidated_dir,
            project_dir,
            filter_lib_names,
            args.lib_prefix,
        )
        return 0
    
    # Validate inputs
    for input_dir in args.input_dirs:
        if not os.path.isdir(input_dir):
            print(f"ERROR: Input directory not found: {input_dir}")
            return 1
    
    if not os.path.isfile(args.atac_ref):
        print(f"ERROR: ATAC reference not found: {args.atac_ref}")
        return 1
    if not os.path.isfile(args.rna_whitelist):
        print(f"ERROR: RNA whitelist not found: {args.rna_whitelist}")
        return 1
    if not os.path.isfile(args.atac_whitelist):
        print(f"ERROR: ATAC whitelist not found: {args.atac_whitelist}")
        return 1
    if not os.path.isfile(args.workflow_file):
        print(f"ERROR: ATAC Nextflow workflow not found: {args.workflow_file}")
        return 1
    
    os.makedirs(args.output_base, exist_ok=True)
    
    print("\n" + "="*70)
    print("INTEGRATED ATAC-SEQ MAPPING PIPELINE - V3")
    print("="*70)
    print(f"\nInput directories: {len(args.input_dirs)}")
    for d in args.input_dirs:
        print(f"  - {d}")
    print(f"\nOutput base: {args.output_base}")
    if args.libraries:
        print(f"Libraries: {args.libraries}")
    
    # Identify libraries
    print("\n" + "="*70)
    print("STEP 1: IDENTIFYING LIBRARIES")
    print("="*70)
    
    print("\nScanning input directories for ATAC triplets...")
    all_libraries = set()
    library_sizes = defaultdict(int)
    library_triplets = defaultdict(list)
    total_triplets = 0
    
    for input_dir in args.input_dirs:
        triplets = find_fastq_triplets(input_dir)
        run_id = os.path.basename(input_dir)
        
        if not triplets:
            print(f"\n  WARNING: No ATAC triplets found in {run_id}")
            continue
        
        print(f"\n  {run_id}: {len(triplets)} triplets found")
        
        for r1, r2, r3 in triplets:
            lib_name = extract_library_name(r1, args.lib_prefix)
            if filter_lib_names is None or lib_name in filter_lib_names:
                all_libraries.add(lib_name)
                total_triplets += 1
                triplet_size = get_file_size(r1) + get_file_size(r2) + get_file_size(r3)
                library_sizes[lib_name] += triplet_size
                library_triplets[lib_name].append((r1, r2, r3))
    
    libraries = sorted(list(all_libraries))
    
    if len(libraries) == 0:
        print("\n  ERROR: No libraries found to process!")
        return 1
    
    print(f"\n  Found {len(libraries)} unique libraries ({total_triplets} total triplets)")
    
    # Test-run mode
    if args.test_run:
        print("\n" + "-"*70)
        print("TEST-RUN MODE: Selecting smallest library")
        smallest_lib = min(library_sizes.keys(), key=lambda x: library_sizes[x])
        smallest_size = library_sizes[smallest_lib]
        size_str = f"{smallest_size/1e9:.2f} GB" if smallest_size > 1e9 else f"{smallest_size/1e6:.2f} MB"
        print(f"  Smallest library: {smallest_lib} ({size_str})")
        libraries = [smallest_lib]
        print("-"*70)
    
    # Generate libs file
    print("\n" + "="*70)
    print("STEP 2: GENERATING CONFIGURATION")
    print("="*70)
    
    libs_file = generate_libs_file(libraries, args.output_base)
    print(f"\n  Created: {libs_file}")
    
    print(f"\nLibraries to be mapped:")
    for lib in sorted(libraries):
        print(f"  - {lib}")
    
    consolidated_dir = os.path.join(args.output_base, 'mapping_input')
    os.makedirs(consolidated_dir, exist_ok=True)
    
    # Generate mapping script
    print("\n" + "="*70)
    print("STEP 3: GENERATING MAPPING SCRIPT")
    print("="*70)
    
    mapping_script = generate_mapping_script(
        args.output_base, args.input_dirs, libs_file,
        args.atac_ref, args.rna_whitelist, args.atac_whitelist,
        args.memgb, args.threads, args.num_chunks, args.workflow_file,
        filter_lib_nums=args.libraries, lib_prefix=args.lib_prefix
    )
    print(f"\n  Generated: {mapping_script}")
    
    if args.submit_jobs:
        result = subprocess.run(['sbatch', mapping_script],
                              capture_output=True, text=True)
        if result.returncode == 0:
            job_id = result.stdout.strip().split()[-1]
            print(f"  Submitted: Job ID {job_id}")
        else:
            print(f"  Failed to submit: {result.stderr}")
    else:
        print(f"  To submit: sbatch {mapping_script}")
    
    print("\n" + "="*70)
    print("PIPELINE SETUP COMPLETE")
    print("="*70)
    print(f"\nOutput directories:")
    print(f"  - Consolidated input: {args.output_base}/mapping_input/")
    print(f"  - Nextflow work: {args.output_base}/mapping_project/")
    print(f"  - Final results: {args.output_base}/mapping_output/")
    
    if args.submit_jobs:
        print(f"\nMonitor with: squeue -u $USER")
    
    print("\n" + "="*70 + "\n")
    return 0


if __name__ == '__main__':
    sys.exit(main())
