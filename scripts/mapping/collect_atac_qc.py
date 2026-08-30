#!/usr/bin/env python3
"""
ATAC-seq QC Statistics Collection - V3

Collects QC metrics from ATAC-seq BAM and fragment files with proper filtering:
- Uses RNA-seq filtered barcodes as the valid cell list
- Analyzes with and without mitochondrial reads
- Reports metrics for all 4 combinations
- FIXED: Samples fragment sizes throughout entire file for valid_bc slices

Usage:
  python collect_atac_qc.py --generate-test-sbatch --library 1
  python collect_atac_qc.py --generate-all-sbatch --submit
  python collect_atac_qc.py --run --library 1 --ramdisk /dev/shm/atac_qc_12345
"""

import os
import sys
import re
import json
import subprocess
import gzip
import argparse
from pathlib import Path
from collections import defaultdict
import numpy as np

# ============================================================================
# CONFIGURATION
# ============================================================================

ATAC_BASE_PATH = "/mnt/beegfs/tetmultiome_atac/mapping_output"
RNA_BASE_PATH = "/mnt/beegfs/tetmultiome_rna_mapped/mapping_output"
STATS_DIR = "/mnt/beegfs/tetmultiome_atac/mapping_output/atac_qc_stats"
SCRIPT_PATH = str(Path(__file__).resolve())

MITO_CHROM = "chrM"

# Number of fragment sizes to sample PER SLICE
# We sample throughout the file, not just first N lines
FRAG_SIZE_SAMPLE_TARGET = 500000

# ============================================================================
# SBATCH TEMPLATE
# ============================================================================

SBATCH_TEMPLATE = """#!/bin/bash
#SBATCH --job-name=atac_qc_{lib_num}
#SBATCH --output={stats_dir}/logs/atac_qc_lib{lib_num}_%j.out
#SBATCH --error={stats_dir}/logs/atac_qc_lib{lib_num}_%j.err
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task={cpus}
#SBATCH --mem={mem}
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --exclude=squirtle
{exclusive}

echo "=========================================="
echo "ATAC QC Collection: Library {lib_num}"
echo "Started at: $(date)"
echo "Running on: $(hostname)"
echo "=========================================="

RAMDISK=/dev/shm/atac_qc_${{SLURM_JOB_ID}}
mkdir -p ${{RAMDISK}}

cleanup() {{
    echo "Cleaning up RAM disk..."
    rm -rf ${{RAMDISK}}
    echo "Finished at: $(date)"
}}
trap cleanup EXIT

module purge
module load miniforge/3
module load htslib/1.20
module load samtools/1.20
module load genomics-base

ATAC_DIR="{atac_base}/{atac_lib_name}"
RNA_DIR="{rna_base}/{rna_lib_name}"

BAM_ORIG="${{ATAC_DIR}}/atac.bam"
BAI_ORIG="${{ATAC_DIR}}/atac.bam.bai"
FRAG_ORIG="${{ATAC_DIR}}/atac_fragments.tsv.gz"
BARCODES_ORIG="${{RNA_DIR}}/filtered/barcodes.tsv.gz"

echo "Checking input files..."
for f in "$BAM_ORIG" "$FRAG_ORIG" "$BARCODES_ORIG"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: File not found: $f"
        exit 1
    fi
done

echo "Copying files to RAM disk..."
time cp $BAM_ORIG ${{RAMDISK}}/atac.bam
[ -f "$BAI_ORIG" ] && cp $BAI_ORIG ${{RAMDISK}}/atac.bam.bai
time cp $FRAG_ORIG ${{RAMDISK}}/atac_fragments.tsv.gz
cp $BARCODES_ORIG ${{RAMDISK}}/barcodes.tsv.gz

echo "Running QC collection..."
python3 {script_path} --run --library {lib_num} --ramdisk ${{RAMDISK}} --threads $SLURM_CPUS_PER_TASK

echo "COMPLETED: Library {lib_num}"
"""


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def get_lib_names(lib_num):
    atac_name = f"Tet_2025_Multiome-ATAC_{lib_num}"
    rna_name = f"Tet_2025_Multiome-RNA_{lib_num}"
    return atac_name, rna_name


def get_all_library_numbers():
    base = Path(ATAC_BASE_PATH)
    nums = []
    for d in base.iterdir():
        if d.is_dir() and d.name.startswith('Tet_2025_Multiome-ATAC_'):
            match = re.search(r'ATAC_(\d+)$', d.name)
            if match:
                nums.append(int(match.group(1)))
    return sorted(nums)


def get_smallest_library():
    base = Path(ATAC_BASE_PATH)
    smallest_num = None
    smallest_size = float('inf')
    for lib_num in get_all_library_numbers():
        atac_name, _ = get_lib_names(lib_num)
        bam = base / atac_name / "atac.bam"
        if bam.exists():
            size = bam.stat().st_size
            if size < smallest_size:
                smallest_size = size
                smallest_num = lib_num
    return smallest_num, smallest_size


def load_rna_barcodes(barcode_file):
    barcodes = set()
    opener = gzip.open if str(barcode_file).endswith('.gz') else open
    with opener(barcode_file, 'rt') as f:
        for line in f:
            bc = line.strip()
            if bc:
                barcodes.add(bc)
    return barcodes


# ============================================================================
# BAM ANALYSIS
# ============================================================================

def analyze_bam(bam_path, valid_barcodes, threads=4):
    print(f"  Analyzing BAM file...")
    print(f"    Valid barcodes from RNA: {len(valid_barcodes):,}")
    
    barcode_reads_all = defaultdict(int)
    barcode_reads_no_mito = defaultdict(int)
    
    total_reads = 0
    mapped_reads = 0
    duplicate_reads = 0
    mito_reads = 0
    
    cmd = ['samtools', 'view', '-@', str(threads), str(bam_path)]
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    
    for i, line in enumerate(process.stdout):
        if i % 10000000 == 0 and i > 0:
            print(f"      Processed {i:,} reads...")
        
        fields = line.strip().split('\t')
        if len(fields) < 11:
            continue
        
        total_reads += 1
        flag = int(fields[1])
        chrom = fields[2]
        
        if not (flag & 4):
            mapped_reads += 1
        if flag & 1024:
            duplicate_reads += 1
        
        is_mito = (chrom == MITO_CHROM)
        if is_mito:
            mito_reads += 1
        
        cb = None
        for field in fields[11:]:
            if field.startswith('CB:Z:'):
                cb = field[5:]
                break
        
        if cb:
            barcode_reads_all[cb] += 1
            if not is_mito:
                barcode_reads_no_mito[cb] += 1
    
    process.wait()
    
    print(f"    Total reads: {total_reads:,}")
    print(f"    Mapped reads: {mapped_reads:,}")
    
    stats = {
        'total_reads': total_reads,
        'mapped_reads': mapped_reads,
        'duplicate_reads': duplicate_reads,
        'mito_reads': mito_reads,
        'mapping_rate': (mapped_reads / total_reads * 100) if total_reads > 0 else 0,
        'duplicate_rate': (duplicate_reads / total_reads * 100) if total_reads > 0 else 0,
        'mito_fraction': (mito_reads / mapped_reads * 100) if mapped_reads > 0 else 0,
        'total_barcodes_in_bam': len(barcode_reads_all),
    }
    
    slices = {
        'all_bc_all_reads': (set(barcode_reads_all.keys()), barcode_reads_all),
        'all_bc_no_mito': (set(barcode_reads_no_mito.keys()), barcode_reads_no_mito),
        'valid_bc_all_reads': (valid_barcodes, {bc: barcode_reads_all.get(bc, 0) for bc in valid_barcodes}),
        'valid_bc_no_mito': (valid_barcodes, {bc: barcode_reads_no_mito.get(bc, 0) for bc in valid_barcodes}),
    }
    
    for slice_name, (bc_set, reads_dict) in slices.items():
        read_counts = [reads_dict[bc] for bc in bc_set if reads_dict.get(bc, 0) > 0]
        if read_counts:
            stats[f'{slice_name}_num_cells'] = len(read_counts)
            stats[f'{slice_name}_total_reads'] = sum(read_counts)
            stats[f'{slice_name}_median_reads_per_cell'] = float(np.median(read_counts))
            stats[f'{slice_name}_mean_reads_per_cell'] = float(np.mean(read_counts))
            stats[f'{slice_name}_frac_reads'] = (sum(read_counts) / total_reads * 100) if total_reads > 0 else 0
        else:
            for suffix in ['_num_cells', '_total_reads', '_median_reads_per_cell', '_mean_reads_per_cell', '_frac_reads']:
                stats[f'{slice_name}{suffix}'] = 0
    
    return stats


# ============================================================================
# FRAGMENT FILE ANALYSIS - FIXED SAMPLING
# ============================================================================

def analyze_fragments(frag_path, valid_barcodes):
    """
    Analyze fragment file with proper sampling for all slices.
    
    KEY FIX: We now collect fragment sizes for each slice separately,
    sampling throughout the entire file rather than just first N lines.
    """
    print(f"  Analyzing fragment file...")
    print(f"    Valid barcodes from RNA: {len(valid_barcodes):,}")
    
    bc_frags_all = defaultdict(int)
    bc_frags_no_mito = defaultdict(int)
    
    # Collect fragment sizes for each slice separately
    # We'll use reservoir sampling to get a representative sample
    frag_sizes = {
        'all_bc_all_reads': [],
        'all_bc_no_mito': [],
        'valid_bc_all_reads': [],
        'valid_bc_no_mito': [],
    }
    
    total_fragments = 0
    mito_fragments = 0
    
    # Count totals for each slice (for reservoir sampling)
    slice_counts = {k: 0 for k in frag_sizes.keys()}
    
    opener = gzip.open if str(frag_path).endswith('.gz') else open
    
    # First pass: count fragments and do reservoir sampling
    np.random.seed(42)  # Reproducibility
    
    with opener(frag_path, 'rt') as f:
        for i, line in enumerate(f):
            if i % 10000000 == 0 and i > 0:
                print(f"      Processed {i:,} fragments...")
            
            parts = line.strip().split('\t')
            if len(parts) < 4:
                continue
            
            chrom, start, end, barcode = parts[0], int(parts[1]), int(parts[2]), parts[3]
            frag_size = end - start
            
            total_fragments += 1
            is_mito = (chrom == MITO_CHROM)
            is_valid = (barcode in valid_barcodes)
            
            if is_mito:
                mito_fragments += 1
            
            bc_frags_all[barcode] += 1
            if not is_mito:
                bc_frags_no_mito[barcode] += 1
            
            # Reservoir sampling for each slice
            def reservoir_sample(slice_key, should_include):
                if not should_include:
                    return
                slice_counts[slice_key] += 1
                n = slice_counts[slice_key]
                
                if len(frag_sizes[slice_key]) < FRAG_SIZE_SAMPLE_TARGET:
                    frag_sizes[slice_key].append(frag_size)
                else:
                    # Reservoir sampling: replace with probability target/n
                    j = np.random.randint(0, n)
                    if j < FRAG_SIZE_SAMPLE_TARGET:
                        frag_sizes[slice_key][j] = frag_size
            
            reservoir_sample('all_bc_all_reads', True)
            reservoir_sample('all_bc_no_mito', not is_mito)
            reservoir_sample('valid_bc_all_reads', is_valid)
            reservoir_sample('valid_bc_no_mito', is_valid and not is_mito)
    
    print(f"    Total fragments: {total_fragments:,}")
    print(f"    Mitochondrial fragments: {mito_fragments:,}")
    print(f"    Fragment sizes sampled per slice:")
    for k, v in frag_sizes.items():
        print(f"      {k}: {len(v):,}")
    
    stats = {
        'total_fragments': total_fragments,
        'mito_fragments': mito_fragments,
        'mito_fragment_fraction': (mito_fragments / total_fragments * 100) if total_fragments > 0 else 0,
        'total_barcodes_in_frags': len(bc_frags_all),
    }
    
    # Calculate per-slice statistics
    slice_data = {
        'all_bc_all_reads': (set(bc_frags_all.keys()), bc_frags_all),
        'all_bc_no_mito': (set(bc_frags_no_mito.keys()), bc_frags_no_mito),
        'valid_bc_all_reads': (valid_barcodes, {bc: bc_frags_all.get(bc, 0) for bc in valid_barcodes}),
        'valid_bc_no_mito': (valid_barcodes, {bc: bc_frags_no_mito.get(bc, 0) for bc in valid_barcodes}),
    }
    
    for slice_name, (bc_set, frags_dict) in slice_data.items():
        frag_counts = [frags_dict[bc] for bc in bc_set if frags_dict.get(bc, 0) > 0]
        sizes = frag_sizes[slice_name]
        
        if frag_counts:
            stats[f'{slice_name}_num_cells'] = len(frag_counts)
            stats[f'{slice_name}_total_frags'] = sum(frag_counts)
            stats[f'{slice_name}_median_frags_per_cell'] = float(np.median(frag_counts))
            stats[f'{slice_name}_mean_frags_per_cell'] = float(np.mean(frag_counts))
            stats[f'{slice_name}_frac_frags'] = (sum(frag_counts) / total_fragments * 100) if total_fragments > 0 else 0
        else:
            for suffix in ['_num_cells', '_total_frags', '_median_frags_per_cell', '_mean_frags_per_cell', '_frac_frags']:
                stats[f'{slice_name}{suffix}'] = 0
        
        if sizes:
            stats[f'{slice_name}_median_frag_size'] = float(np.median(sizes))
            stats[f'{slice_name}_mean_frag_size'] = float(np.mean(sizes))
            
            n_frags = len(sizes)
            stats[f'{slice_name}_nfr_fraction'] = sum(1 for s in sizes if s < 147) / n_frags * 100
            stats[f'{slice_name}_mono_nuc_fraction'] = sum(1 for s in sizes if 147 <= s < 294) / n_frags * 100
            stats[f'{slice_name}_di_nuc_fraction'] = sum(1 for s in sizes if 294 <= s < 441) / n_frags * 100
            
            hist, _ = np.histogram(sizes, bins=np.arange(0, 1001, 5))
            stats[f'{slice_name}_frag_size_hist'] = hist.tolist()
    
    stats['frag_size_bins'] = np.arange(0, 1001, 5)[:-1].tolist()
    
    return stats


# ============================================================================
# MAIN COLLECTION
# ============================================================================

def collect_library_stats(lib_num, ramdisk, threads=4):
    print(f"\n{'='*60}")
    print(f"Collecting stats for Library {lib_num}")
    print(f"{'='*60}")
    
    ramdisk = Path(ramdisk)
    bam_path = ramdisk / "atac.bam"
    frag_path = ramdisk / "atac_fragments.tsv.gz"
    barcode_path = ramdisk / "barcodes.tsv.gz"
    
    atac_name, rna_name = get_lib_names(lib_num)
    
    stats = {
        'library_number': lib_num,
        'atac_library': atac_name,
        'rna_library': rna_name,
    }
    
    print(f"\nLoading RNA-seq filtered barcodes...")
    valid_barcodes = load_rna_barcodes(barcode_path)
    stats['rna_filtered_cells'] = len(valid_barcodes)
    print(f"  Loaded {len(valid_barcodes):,} valid barcodes")
    
    if bam_path.exists():
        print(f"\nBAM Analysis:")
        bam_stats = analyze_bam(bam_path, valid_barcodes, threads)
        stats.update(bam_stats)
    
    if frag_path.exists():
        print(f"\nFragment Analysis:")
        frag_stats = analyze_fragments(frag_path, valid_barcodes)
        stats.update(frag_stats)
    
    # Summary
    prefix = 'valid_bc_no_mito'
    print(f"\n{'='*60}")
    print(f"SUMMARY for Library {lib_num}")
    print(f"{'='*60}")
    print(f"  RNA-filtered cells: {stats.get('rna_filtered_cells', 0):,}")
    print(f"  Mapping rate: {stats.get('mapping_rate', 0):.1f}%")
    print(f"  Duplicate rate: {stats.get('duplicate_rate', 0):.1f}%")
    print(f"  Mito fraction: {stats.get('mito_fraction', 0):.1f}%")
    print(f"  Valid cells (no mito):")
    print(f"    Cells with ATAC: {stats.get(f'{prefix}_num_cells', 0):,}")
    print(f"    Median reads/cell: {stats.get(f'{prefix}_median_reads_per_cell', 0):,.0f}")
    print(f"    Median frags/cell: {stats.get(f'{prefix}_median_frags_per_cell', 0):,.0f}")
    print(f"    Frac reads in cells: {stats.get(f'{prefix}_frac_reads', 0):.1f}%")
    
    return stats


def save_stats(stats, lib_num):
    stats_dir = Path(STATS_DIR)
    stats_dir.mkdir(parents=True, exist_ok=True)
    output_file = stats_dir / f"library_{lib_num}_stats.json"
    with open(output_file, 'w') as f:
        json.dump(stats, f, indent=2)
    print(f"\nSaved stats to: {output_file}")
    return output_file


# ============================================================================
# SBATCH GENERATION
# ============================================================================

# Libraries that need more time/resources
SLOW_LIBRARIES = {1, 2, 3, 4, 30, 32, 38}

def generate_sbatch_script(lib_num, test_mode=False):
    atac_name, rna_name = get_lib_names(lib_num)
    
    is_slow = lib_num in SLOW_LIBRARIES
    
    if test_mode:
        cpus = 96
        mem = "500G"
        time = "4:00:00"
        exclusive = "#SBATCH --exclusive"
    elif is_slow:
        cpus = 24
        mem = "300G"
        time = "12:00:00"
        exclusive = ""
    else:
        cpus = 16
        mem = "250G"
        time = "6:00:00"
        exclusive = ""
    
    logs_dir = Path(STATS_DIR) / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)
    sbatch_dir = Path(STATS_DIR) / "sbatch"
    sbatch_dir.mkdir(parents=True, exist_ok=True)
    
    # Update template with time
    template = SBATCH_TEMPLATE.replace("#SBATCH --time=6:00:00", f"#SBATCH --time={time}")
    
    script_content = template.format(
        lib_num=lib_num, atac_lib_name=atac_name, rna_lib_name=rna_name,
        atac_base=ATAC_BASE_PATH, rna_base=RNA_BASE_PATH,
        stats_dir=STATS_DIR, script_path=SCRIPT_PATH,
        cpus=cpus, mem=mem, exclusive=exclusive
    )
    
    script_path = sbatch_dir / f"collect_qc_lib{lib_num}.sbatch"
    with open(script_path, 'w') as f:
        f.write(script_content)
    os.chmod(script_path, 0o755)
    return script_path, is_slow


# ============================================================================
# MAIN
# ============================================================================

def main():
    parser = argparse.ArgumentParser(description='ATAC-seq QC Statistics Collection V3')
    parser.add_argument('--generate-test-sbatch', action='store_true')
    parser.add_argument('--generate-all-sbatch', action='store_true')
    parser.add_argument('--library', type=int)
    parser.add_argument('--submit', action='store_true')
    parser.add_argument('--list-libraries', action='store_true')
    parser.add_argument('--run', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('--ramdisk', type=str, help=argparse.SUPPRESS)
    parser.add_argument('--threads', type=int, default=4, help=argparse.SUPPRESS)
    
    args = parser.parse_args()
    
    if args.list_libraries:
        print(f"\nLibraries:\n")
        for lib_num in get_all_library_numbers():
            atac_name, rna_name = get_lib_names(lib_num)
            bam = Path(ATAC_BASE_PATH) / atac_name / "atac.bam"
            bc = Path(RNA_BASE_PATH) / rna_name / "filtered" / "barcodes.tsv.gz"
            bam_size = f"{bam.stat().st_size / (1024**3):.1f} GB" if bam.exists() else "missing"
            bc_status = "OK" if bc.exists() else "MISSING"
            print(f"  {lib_num:2d}: ATAC {bam_size:>10s} | RNA barcodes: {bc_status}")
        return
    
    if args.run:
        if args.library is None or not args.ramdisk:
            print("ERROR: --run requires --library and --ramdisk")
            sys.exit(1)
        stats = collect_library_stats(args.library, args.ramdisk, args.threads)
        save_stats(stats, args.library)
        return
    
    if args.generate_test_sbatch:
        lib_num = args.library if args.library else get_smallest_library()[0]
        if lib_num is None:
            print("ERROR: No libraries found")
            sys.exit(1)
        script, is_slow = generate_sbatch_script(lib_num, test_mode=True)
        print(f"Generated: {script}")
        if args.submit:
            result = subprocess.run(['sbatch', str(script)], capture_output=True, text=True)
            print(f"Submitted: {result.stdout.strip()}" if result.returncode == 0 else f"Failed: {result.stderr}")
        return
    
    if args.generate_all_sbatch:
        lib_nums = get_all_library_numbers()
        print(f"\nGenerating sbatch scripts for {len(lib_nums)} libraries...")
        print(f"Slow libraries (24 cores, 12hr): {sorted(SLOW_LIBRARIES)}")
        
        normal_scripts = []
        slow_scripts = []
        
        for lib_num in lib_nums:
            stats_file = Path(STATS_DIR) / f"library_{lib_num}_stats.json"
            if stats_file.exists():
                print(f"  Skipping Library {lib_num} (already done)")
                continue
            script, is_slow = generate_sbatch_script(lib_num, test_mode=False)
            if is_slow:
                slow_scripts.append(script)
                print(f"  Generated: {script.name} [SLOW - 24 cores, 12hr]")
            else:
                normal_scripts.append(script)
                print(f"  Generated: {script.name}")
        
        print(f"\nGenerated {len(normal_scripts)} normal + {len(slow_scripts)} slow scripts")
        
        if args.submit:
            print("\nSubmitting normal jobs first...")
            for script in normal_scripts:
                result = subprocess.run(['sbatch', str(script)], capture_output=True, text=True)
                print(f"  {script.name}: {result.stdout.strip()}" if result.returncode == 0 else f"  FAILED: {result.stderr}")
            
            print("\nSubmitting slow jobs (will queue behind normal jobs)...")
            for script in slow_scripts:
                result = subprocess.run(['sbatch', str(script)], capture_output=True, text=True)
                print(f"  {script.name}: {result.stdout.strip()}" if result.returncode == 0 else f"  FAILED: {result.stderr}")
        return
    
    parser.print_help()


if __name__ == "__main__":
    main()
