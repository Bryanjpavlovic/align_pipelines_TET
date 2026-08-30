#!/usr/bin/env python3
"""Generate collision-safe 10X RNA trimming and STARsolo SLURM jobs.

V10 completes the pieces that were split across the historical V9 branches:

* collision-safe ``__RunNNN`` FASTQ names and source-to-BAM read groups;
* automatic validation of standard/long-R2 versus 5' PE150 geometry;
* PE150 R1 handling that preserves CB16+UMI12, removes the internal 13 bp TSO,
  trims the retained R1 cDNA, and maps both cDNA mates;
* legitimate single-end R2 cutadapt info files for barcode-linked trimming QC;
* direct use of the custom ``align_pipelines/bjp`` module; a mismatch between
  the requested and installed workflow is reported but does not block a run.

The normal mode only generates scripts unless ``--submit-jobs`` is supplied.
The higher-level ``orchestrate_10x_mapping_qc.py`` normally handles submission.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import os
import re
import shlex
import subprocess
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from statistics import median
from typing import Iterable, Sequence


RELEASE = "2026-08-30-v17-align-repo-migration"
ALIGN_PIPELINES_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_RNA_WORKFLOW = str(ALIGN_PIPELINES_ROOT / "workflows" / "align_rna.nf")
DEFAULT_RNA3_LIB_PREFIX = "Tet_2025_Multiome-RNA_"
DEFAULT_RNA5_LIB_PREFIX = "Tet_2025_RNA_5P_"
CB_LENGTH = 16
UMI_LENGTH = 12
BARCODE_BLOCK_LENGTH = CB_LENGTH + UMI_LENGTH
FIVE_PRIME_TSO = "TTTCTTATATGGG"


class PipelineError(RuntimeError):
    """A user-facing pipeline configuration or data error."""


def q(value: os.PathLike[str] | str) -> str:
    return shlex.quote(str(value))


def open_text(path: Path, mode: str):
    return gzip.open(path, mode) if path.name.endswith(".gz") else path.open(mode)


def fastq_pairs(input_dir: Path) -> list[tuple[Path, Path]]:
    """Return complete R1/R2 pairs and fail if an R1 mate is missing."""
    r1_files: set[Path] = set()
    for pattern in ("*_R1_*.fastq.gz", "*_R1_*.fq.gz"):
        r1_files.update(input_dir.glob(pattern))
    pairs: list[tuple[Path, Path]] = []
    missing: list[Path] = []
    for r1 in sorted(r1_files):
        r2 = Path(str(r1).replace("_R1_", "_R2_"))
        if r2.is_file():
            # Preserve the named FASTQ path in provenance tables even when a
            # staging/test input is itself a symlink.  Collision checks still
            # compare resolved targets where identity matters.
            pairs.append((r1.absolute(), r2.absolute()))
        else:
            missing.append(r2)
    if missing:
        preview = ", ".join(str(path) for path in missing[:5])
        raise PipelineError(f"missing R2 mate(s): {preview}")
    return pairs


def sample_lengths(path: Path, n_reads: int = 100) -> list[int]:
    lengths: list[int] = []
    with open_text(path, "rt") as handle:
        while len(lengths) < n_reads:
            header = handle.readline()
            if not header:
                break
            sequence = handle.readline().rstrip("\r\n")
            plus = handle.readline()
            quality = handle.readline().rstrip("\r\n")
            if not header.startswith("@") or not plus.startswith("+"):
                raise PipelineError(f"invalid FASTQ record in {path}")
            if len(sequence) != len(quality):
                raise PipelineError(f"sequence/quality length mismatch in {path}")
            lengths.append(len(sequence))
    if not lengths:
        raise PipelineError(f"no FASTQ records found in {path}")
    return lengths


def detect_read_format(r1: Path, r2: Path) -> tuple[str, float, float]:
    """Classify one pair using median lengths from up to 100 records."""
    r1_length = float(median(sample_lengths(r1)))
    r2_length = float(median(sample_lengths(r2)))
    if 24 <= r1_length <= 35 and r2_length >= 80:
        return "long-r2", r1_length, r2_length
    if r1_length >= 100 and r2_length >= 100:
        return "pe150", r1_length, r2_length
    return "unknown", r1_length, r2_length


def validate_run_format(
    pairs: Sequence[tuple[Path, Path]], chemistry: str, requested: str
) -> str:
    detected: set[str] = set()
    details: list[str] = []
    for r1, r2 in pairs:
        fmt, r1_length, r2_length = detect_read_format(r1, r2)
        detected.add(fmt)
        details.append(f"{r1.name}: R1={r1_length:.0f}, R2={r2_length:.0f}, {fmt}")
    if "unknown" in detected or len(detected) != 1:
        raise PipelineError(
            "a run must contain one recognized read geometry; observed: "
            + "; ".join(details)
        )
    observed = next(iter(detected))
    if requested != "auto" and requested != observed:
        raise PipelineError(
            f"--read-format {requested} disagrees with detected {observed}: "
            + "; ".join(details)
        )
    if chemistry == "3prime" and observed != "long-r2":
        raise PipelineError(
            f"3prime currently requires a short barcode R1 plus cDNA R2; detected {observed}"
        )
    return observed


def library_name(path: Path, prefix: str | None = None) -> str:
    """Canonicalize legacy and bcl-convert-annotated library names."""
    match = re.match(r"(.+?)_S\d+_L\d+", path.name)
    if not match:
        raise PipelineError(f"cannot extract library from FASTQ name: {path.name}")
    candidate = match.group(1)
    if prefix:
        prefixed = re.match(rf"({re.escape(prefix)}\d+)(?:_|$)", candidate)
        if prefixed:
            return prefixed.group(1)
    return re.sub(r"(?:_L\d+)+$", "", candidate)


def library_number(name: str, prefix: str | None = None) -> int | None:
    if prefix:
        match = re.fullmatch(rf"{re.escape(prefix)}(\d+)", name)
        if match:
            return int(match.group(1))
    match = re.search(r"_(\d+)$", name)
    return int(match.group(1)) if match else None


def sample_key(r2: Path) -> str:
    name = r2.name.replace("_R2_", "_")
    return re.sub(r"\.(?:fastq|fq)\.gz$", "", name)


def insert_run_tag(name: str, run_tag: str) -> str:
    match = re.match(r"^(.+_S\d+_L\d+)(_R[12]_.+)$", name)
    if not match:
        raise PipelineError(f"FASTQ name lacks _S#_L###_R# structure: {name}")
    return f"{match.group(1)}__{run_tag}{match.group(2)}"


def extract_sample_lane(name: str) -> tuple[str, str]:
    match = re.search(r"_(S\d+)_(L\d+)", name)
    return (match.group(1), match.group(2)) if match else ("S0", "L000")


def extract_flowcell_lane(path: Path) -> tuple[str, str]:
    try:
        with open_text(path, "rt") as handle:
            header = handle.readline().strip()
        fields = header.split()[0].split(":")
        if len(fields) >= 4:
            return fields[2], fields[3]
    except (OSError, EOFError, UnicodeError):
        pass
    return "unknown", "0"


def r2_adapter_lines(chemistry: str, paired: bool) -> list[str]:
    front = "-G" if paired else "-g"
    back = "-A" if paired else "-a"
    if chemistry == "3prime":
        adapters = [
            (front, "TSO=^AAGCAGTGGTATCAACGCAGAGTACATGGG"),
            (back, "TSO_RC=CCCATGTACTCTGCGTTGATACCACTGCTT"),
        ]
    else:
        adapters = [
            (front, "TSO_5prime_marker=^TTTCTTATATGGG"),
            (front, "TSO_5prime_marker_RC=^CCCATATAAGAAA"),
        ]
    adapters.extend(
        [
            (back, "TruSeq_Read2_Universal=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"),
            (back, "TruSeq_Read1_RC=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"),
            (back, "P7_adapter=CAAGCAGAAGACGGCATACGAGAT"),
            (back, "P7_RC=ATCTCGTATGCCGTCTTCTGCTTG"),
            (back, "P5_adapter_RC=GTGTAGATCTCGGTGGTCGCCGTATCATT"),
            (back, "TruSeq_Read2_RC_Primer=TCTAGCCTTCTCGTGTGCAGACTTGAGGTCAGTG"),
        ]
    )
    if chemistry == "5prime":
        adapters.append(
            (back, "TruSeq_Read2_Primer=GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT")
        )
    adapters.append((back, "PolyG=G{20}"))
    return [f"  {flag} {q(spec)} \\" for flag, spec in adapters]


def pe150_r1_adapter_lines() -> list[str]:
    specs = [
        "TruSeq_Read2_Universal=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
        "TruSeq_Read1_RC=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
        "P7_adapter=CAAGCAGAAGACGGCATACGAGAT",
        "P7_RC=ATCTCGTATGCCGTCTTCTGCTTG",
        "P5_adapter_RC=GTGTAGATCTCGGTGGTCGCCGTATCATT",
        "TruSeq_Read2_RC_Primer=TCTAGCCTTCTCGTGTGCAGACTTGAGGTCAGTG",
        "TruSeq_Read2_Primer=GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT",
        "PolyG=G{20}",
    ]
    return [f"  -a {q(spec)} \\" for spec in specs]


def bash_array(name: str, values: Iterable[Path | str]) -> str:
    return name + "=(\n" + "".join(f"  {q(value)}\n" for value in values) + ")\n"


def generate_trim_script(
    input_dir: Path,
    output_base: Path,
    chemistry: str,
    read_format: str,
    no_info_file: bool,
    dry_run: bool,
    min_tso_match_fraction: float,
    max_cores: int | None,
    library_numbers: set[int] | None,
    filter_names: set[str] | None,
    lib_prefix: str,
) -> Path:
    pairs = [
        pair
        for pair in fastq_pairs(input_dir)
        if selected_library(
            library_name(pair[0], lib_prefix),
            library_numbers,
            filter_names,
            lib_prefix,
        )
    ]
    if not pairs:
        raise PipelineError(
            f"no selected R1/R2 FASTQ pairs found in {input_dir}"
        )
    observed = validate_run_format(pairs, chemistry, read_format)
    run_id = input_dir.name
    trim_dir = output_base / "trimming" / run_id
    log_dir = trim_dir / "logs"
    trim_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)
    suffix = "_noinfo" if no_info_file else ""
    script_path = output_base / f"trim_{chemistry}_{run_id}{suffix}.sh"
    samples = [sample_key(r2) for _, r2 in pairs]
    driver = Path(__file__).resolve()
    trim_cpus = min(8, max_cores) if max_cores is not None else 8
    array_limit = max(1, max_cores // trim_cpus) if max_cores is not None else None
    array_spec = f"0-{len(pairs) - 1}"
    if array_limit is not None:
        array_spec += f"%{array_limit}"

    r1_adapter_text = "\n".join(pe150_r1_adapter_lines()) + "\n" if observed == "pe150" else ""
    r2_adapter_text = "\n".join(r2_adapter_lines(chemistry, paired=True)) + "\n"
    r2_record_text = "\n".join(r2_adapter_lines(chemistry, paired=False)) + "\n"

    if observed == "pe150":
        prepare = f"""PREPARED_R1=\"${{TASK_TMP}}/prepared_R1.fastq.gz\"
TSO_INFO_TMP=\"${{TASK_TMP}}/R1_tso_info.tsv.gz\"
python3 {q(driver)} --prepare-5p-pe150-r1 \\
  --worker-input \"${{R1_FILE}}\" \\
  --worker-output \"${{PREPARED_R1}}\" \\
  --worker-tso-info \"${{TSO_INFO_TMP}}\" \\
  --min-pe150-tso-match-fraction {min_tso_match_fraction}
TRIM_INPUT_R1=\"${{PREPARED_R1}}\"
"""
        min_length = "53:25"
    else:
        prepare = 'TRIM_INPUT_R1="${R1_FILE}"\n'
        min_length = "25:25"

    info_vars = ""
    info_main_flag = ""
    info_sidecar = ""
    info_publish = ""
    info_guard = ""
    if not no_info_file:
        info_vars = """R2_INFO=\"${LOGDIR}/${SAMPLE}_R2_adapter_info.txt.gz\"
R2_RECORD_JSON=\"${LOGDIR}/${SAMPLE}_R2_record.cutadapt.json\"
R2_RECORD_REPORT=\"${LOGDIR}/${SAMPLE}_R2_record_report.txt\"
R2_INFO_TMP=\"${TASK_TMP}/R2_adapter_info.txt.gz\"
R2_RECORD_JSON_TMP=\"${TASK_TMP}/R2_record.cutadapt.json\"
R2_RECORD_REPORT_TMP=\"${TASK_TMP}/R2_record_report.txt\"
"""
        info_guard = ' && -s "${R2_INFO}" && -s "${R2_RECORD_JSON}" && -s "${R2_RECORD_REPORT}"'
        info_publish = """mv \"${R2_INFO_TMP}\" \"${R2_INFO}\"
mv \"${R2_RECORD_JSON_TMP}\" \"${R2_RECORD_JSON}\"
mv \"${R2_RECORD_REPORT_TMP}\" \"${R2_RECORD_REPORT}\"
"""
        if observed == "pe150":
            info_vars += """R1_INFO=\"${LOGDIR}/${SAMPLE}_R1_adapter_info.txt.gz\"
R1_TSO_INFO=\"${LOGDIR}/${SAMPLE}_R1_tso_info.tsv.gz\"
R1_INFO_TMP=\"${TASK_TMP}/R1_adapter_info.txt.gz\"
"""
            info_main_flag = '  --info-file="${R1_INFO_TMP}" \\\n'
            info_guard += ' && -s "${R1_INFO}" && -s "${R1_TSO_INFO}"'
            info_publish += """mv \"${R1_INFO_TMP}\" \"${R1_INFO}\"
mv \"${TSO_INFO_TMP}\" \"${R1_TSO_INFO}\"
"""
        info_sidecar = f"""cutadapt \\
{r2_record_text}  -O 10 \\
  --times 3 \\
  --info-file=\"${{R2_INFO_TMP}}\" \\
  --json=\"${{R2_RECORD_JSON_TMP}}\" \\
  --cores \"${{SLURM_CPUS_PER_TASK}}\" \\
  -o /dev/null \\
  \"${{R2_FILE}}\" > \"${{R2_RECORD_REPORT_TMP}}\" 2>&1
"""

    if dry_run:
        output_block = """CUTADAPT_R1=/dev/null
CUTADAPT_R2=/dev/null
"""
        publish_block = ""
        guard = ""
        collision_guard = ""
    else:
        output_block = """OUTPUT_R1=\"${OUTDIR}/$(basename \"${R1_FILE}\")\"
OUTPUT_R2=\"${OUTDIR}/$(basename \"${R2_FILE}\")\"
CUTADAPT_R1=\"${TASK_TMP}/trimmed_R1.fastq.gz\"
CUTADAPT_R2=\"${TASK_TMP}/trimmed_R2.fastq.gz\"
"""
        guard = f"""if [[ -s \"${{OUTPUT_R1}}\" && -s \"${{OUTPUT_R2}}\" && -s \"${{JSON_REPORT}}\" && -s \"${{REPORT}}\"{info_guard} ]]; then
  echo \"SKIP: completed trimming outputs already exist for ${{SAMPLE}}\"
  exit 0
fi
"""
        final_variables = ["OUTPUT_R1", "OUTPUT_R2", "JSON_REPORT", "REPORT"]
        if not no_info_file:
            final_variables.extend(["R2_INFO", "R2_RECORD_JSON", "R2_RECORD_REPORT"])
            if observed == "pe150":
                final_variables.extend(["R1_INFO", "R1_TSO_INFO"])
        final_array = " ".join(f'\"${{{name}}}\"' for name in final_variables)
        collision_guard = f"""for FINAL_PATH in {final_array}; do
  if [[ -e \"${{FINAL_PATH}}\" ]]; then
    echo \"ERROR: incomplete prior output exists; refusing to overwrite: ${{FINAL_PATH}}\" >&2
    echo \"Use a new orchestrator --run-name or inspect this staged run.\" >&2
    exit 1
  fi
done
"""
        publish_block = """mv \"${CUTADAPT_R1}\" \"${OUTPUT_R1}\"
mv \"${CUTADAPT_R2}\" \"${OUTPUT_R2}\"
"""

    script = f"""#!/bin/bash
#SBATCH --job-name=trim_{chemistry}_{run_id}
#SBATCH --output={log_dir}/trim_%A_%a.out
#SBATCH --error={log_dir}/trim_%A_%a.err
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={trim_cpus}
#SBATCH --mem=24G
#SBATCH --time=5-00:00:00
#SBATCH --array={array_spec}
#SBATCH --chdir={output_base}

set -euo pipefail
module purge
module load miniforge/3 genomics-base/latest

OUTDIR={q(trim_dir)}
LOGDIR={q(log_dir)}
TEMP_ROOT={q(output_base / 'tmp' / 'trim' / run_id)}
mkdir -p \"${{OUTDIR}}\" \"${{LOGDIR}}\" \"${{TEMP_ROOT}}\"
{bash_array('R1_FILES', [r1 for r1, _ in pairs])}{bash_array('R2_FILES', [r2 for _, r2 in pairs])}{bash_array('SAMPLES', samples)}
R1_FILE=\"${{R1_FILES[$SLURM_ARRAY_TASK_ID]}}\"
R2_FILE=\"${{R2_FILES[$SLURM_ARRAY_TASK_ID]}}\"
SAMPLE=\"${{SAMPLES[$SLURM_ARRAY_TASK_ID]}}\"
REPORT=\"${{LOGDIR}}/${{SAMPLE}}_report.txt\"
JSON_REPORT=\"${{LOGDIR}}/${{SAMPLE}}_cutadapt.json\"
TASK_TMP=$(mktemp -d \"${{TEMP_ROOT}}/tet_rna_trim_${{SLURM_JOB_ID:-manual}}_${{SLURM_ARRAY_TASK_ID}}_XXXXXX\")
trap 'rm -rf \"${{TASK_TMP}}\"' EXIT
MAIN_REPORT_TMP=\"${{TASK_TMP}}/cutadapt_report.txt\"
{info_vars}{output_block}{guard}{collision_guard}
for INPUT in \"${{R1_FILE}}\" \"${{R2_FILE}}\"; do
  [[ -s \"${{INPUT}}\" ]] || {{ echo \"ERROR: missing FASTQ: ${{INPUT}}\" >&2; exit 1; }}
done
command -v cutadapt >/dev/null || {{ echo \"ERROR: cutadapt not found\" >&2; exit 1; }}

echo \"Release: {RELEASE}\"
echo \"Chemistry: {chemistry}; geometry: {observed}\"
echo \"R1: ${{R1_FILE}}\"
echo \"R2: ${{R2_FILE}}\"
{prepare}
cutadapt \\
{r1_adapter_text}{r2_adapter_text}  -O 10 \\
  --minimum-length {min_length} \\
  --times 3 \\
  --pair-filter=any \\
{info_main_flag}  --json=\"${{TASK_TMP}}/cutadapt.json\" \\
  --cores \"${{SLURM_CPUS_PER_TASK}}\" \\
  -o \"${{CUTADAPT_R1}}\" \\
  -p \"${{CUTADAPT_R2}}\" \\
  \"${{TRIM_INPUT_R1}}\" \\
  \"${{R2_FILE}}\" > \"${{MAIN_REPORT_TMP}}\" 2>&1

{info_sidecar}mv \"${{TASK_TMP}}/cutadapt.json\" \"${{JSON_REPORT}}\"
mv \"${{MAIN_REPORT_TMP}}\" \"${{REPORT}}\"
{publish_block}{info_publish}
echo \"SUCCESS: {chemistry} {observed} trimming complete for ${{SAMPLE}}\"
"""
    script_path.write_text(script, encoding="utf-8")
    script_path.chmod(0o755)
    return script_path


def prepare_pe150_r1(
    input_path: Path,
    output_path: Path,
    info_path: Path,
    min_match_fraction: float,
) -> None:
    """Delete the fixed TSO block at R1 bases 29-41 while retaining both flanks.

    A conventional cutadapt 5' adapter cannot do this because it would also
    remove the 28 bp barcode/UMI prefix.  The fixed-position transform is
    therefore explicit and auditable.  Every record is logged.  The job fails
    if the expected TSO is not recognizable in enough reads.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    info_path.parent.mkdir(parents=True, exist_ok=True)
    total = 0
    recognized = 0
    with open_text(input_path, "rt") as source, gzip.open(output_path, "wt") as dest, gzip.open(
        info_path, "wt"
    ) as info:
        info.write("read_name\traw_length\ttso_start\ttso_length\terrors\tstatus\n")
        while True:
            header = source.readline()
            if not header:
                break
            sequence = source.readline().rstrip("\r\n")
            plus = source.readline()
            quality = source.readline().rstrip("\r\n")
            if not header.startswith("@") or not plus.startswith("+"):
                raise PipelineError(f"invalid FASTQ record in {input_path}")
            if len(sequence) != len(quality):
                raise PipelineError(f"sequence/quality mismatch in {input_path}")
            if len(sequence) < BARCODE_BLOCK_LENGTH + len(FIVE_PRIME_TSO):
                raise PipelineError(f"PE150 R1 is too short ({len(sequence)} bp): {input_path}")
            observed = sequence[
                BARCODE_BLOCK_LENGTH : BARCODE_BLOCK_LENGTH + len(FIVE_PRIME_TSO)
            ]
            errors = sum(a != b for a, b in zip(observed, FIVE_PRIME_TSO))
            status = "exact" if errors == 0 else "one_mismatch" if errors == 1 else "unrecognized"
            if errors <= 1:
                recognized += 1
            total += 1
            # The documented library structure fixes the TSO at bases 29-41.
            # Remove exactly that block even when a sequencing error obscures it.
            new_sequence = sequence[:BARCODE_BLOCK_LENGTH] + sequence[BARCODE_BLOCK_LENGTH + len(FIVE_PRIME_TSO) :]
            new_quality = quality[:BARCODE_BLOCK_LENGTH] + quality[BARCODE_BLOCK_LENGTH + len(FIVE_PRIME_TSO) :]
            dest.write(header)
            dest.write(new_sequence + "\n")
            dest.write(plus)
            dest.write(new_quality + "\n")
            read_name = header[1:].split()[0]
            info.write(
                f"{read_name}\t{len(sequence)}\t{BARCODE_BLOCK_LENGTH}\t{len(FIVE_PRIME_TSO)}\t{errors}\t{status}\n"
            )
    fraction = recognized / total if total else 0.0
    if not total:
        raise PipelineError(f"no records found in {input_path}")
    if fraction < min_match_fraction:
        raise PipelineError(
            f"only {recognized}/{total} ({fraction:.3%}) PE150 R1 reads matched the expected "
            f"TSO at base 29; required {min_match_fraction:.3%}. Refusing to map this file."
        )


def selected_library(
    name: str,
    library_numbers: set[int] | None,
    names: set[str] | None,
    lib_prefix: str | None = None,
) -> bool:
    """Apply numeric and explicit-name filters as a union, matching V9."""
    if library_numbers is None and names is None:
        return True
    numeric_match = (
        library_numbers is not None
        and library_number(name, lib_prefix) in library_numbers
    )
    name_match = names is not None and name in names
    return numeric_match or name_match


def folder_sort_key(path: Path) -> tuple[float, str]:
    try:
        return path.stat().st_mtime, path.name
    except OSError:
        return 0.0, path.name


def atomic_write(path: Path, text: str) -> None:
    temp = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    temp.write_text(text, encoding="utf-8")
    os.replace(temp, path)


def consolidate(
    input_dirs: Sequence[Path],
    output_base: Path,
    library_numbers: set[int] | None,
    filter_names: set[str] | None,
    lib_prefix: str,
) -> list[str]:
    mapping_input = output_base / "mapping_input"
    project = output_base / "mapping_project"
    mapping_input.mkdir(parents=True, exist_ok=True)
    project.mkdir(parents=True, exist_ok=True)
    grouped: dict[str, list[tuple[Path, Path, Path]]] = defaultdict(list)
    for input_dir in sorted(input_dirs, key=folder_sort_key):
        trim_dir = output_base / "trimming" / input_dir.name
        pairs = fastq_pairs(trim_dir) if trim_dir.is_dir() else []
        if not pairs:
            raise PipelineError(f"no complete trimmed pairs found in {trim_dir}")
        for r1, r2 in pairs:
            lib = library_name(r1, lib_prefix)
            if selected_library(lib, library_numbers, filter_names, lib_prefix):
                grouped[lib].append((input_dir, r1, r2))
    if not grouped:
        raise PipelineError("no libraries remained after consolidation filters")

    manifest_rows: list[str] = []
    assignment_rows: list[str] = []
    rg_rows: list[str] = []
    destinations: dict[str, Path] = {}
    for lib in sorted(grouped):
        source_folders: list[Path] = []
        for folder, _, _ in grouped[lib]:
            if folder not in source_folders:
                source_folders.append(folder)
        run_map = {folder: f"Run{index + 1:03d}" for index, folder in enumerate(source_folders)}
        for folder in source_folders:
            assignment_rows.append(f"{lib}\t{run_map[folder]}\t{folder.name}\t{folder}\n")
        for folder, r1, r2 in grouped[lib]:
            run_tag = run_map[folder]
            flowcell, flowcell_lane = extract_flowcell_lane(r1)
            sample_idx, lane_id = extract_sample_lane(r1.name)
            r2_new = insert_run_tag(r2.name, run_tag)
            fnbase = re.sub(r"_R[12]_.*$", "", r2_new)
            rg_id = f"{lib}_{folder.name}_{sample_idx}_{lane_id}"
            pu = f"{flowcell}.{flowcell_lane}"
            rg_string = (
                f"@RG\\tID:{rg_id}\\tSM:{lib}\\tPL:Illumina\\tPU:{pu}\\tDS:{folder.name}"
            )
            raw_r1 = (folder / r1.name).resolve(strict=False)
            raw_r2 = (folder / r2.name).resolve(strict=False)
            rg_rows.append(
                "\t".join(
                    [
                        fnbase,
                        lib,
                        folder.name,
                        sample_idx,
                        lane_id,
                        flowcell,
                        flowcell_lane,
                        pu,
                        rg_id,
                        rg_string,
                        str(raw_r1),
                        str(raw_r2),
                        str(r1),
                        str(r2),
                    ]
                )
                + "\n"
            )
            for source in (r1, r2):
                new_name = insert_run_tag(source.name, run_tag)
                if new_name in destinations and destinations[new_name] != source:
                    raise PipelineError(
                        f"run-tag collision: {new_name} maps to both {destinations[new_name]} and {source}"
                    )
                destinations[new_name] = source
                target = mapping_input / new_name
                if target.is_symlink():
                    if target.resolve() != source.resolve():
                        raise PipelineError(f"existing symlink points to a different source: {target}")
                elif target.exists():
                    raise PipelineError(f"refusing to overwrite non-symlink mapping input: {target}")
                else:
                    target.symlink_to(source)
                manifest_rows.append(
                    f"{lib}\t{run_tag}\t{folder.name}\t{folder / source.name}\t"
                    f"{source}\t{target}\t{new_name}\n"
                )

    existing_fastqs: set[Path] = set()
    for pattern in ("*_R[12]_*.fastq.gz", "*_R[12]_*.fq.gz"):
        existing_fastqs.update(mapping_input.glob(pattern))
    unexpected = sorted(path for path in existing_fastqs if path.name not in destinations)
    if unexpected:
        preview = ", ".join(str(path) for path in unexpected[:5])
        raise PipelineError(
            "unexpected pre-existing mapping-input FASTQ(s); refusing to mix staged runs: "
            + preview
        )

    atomic_write(
        project / "symlink_manifest.tsv",
        "library\trun_tag\tbp_id\traw_source_fastq\ttrimmed_fastq\t"
        "symlink_path\tsymlink_name\n"
        + "".join(manifest_rows),
    )
    atomic_write(
        project / "run_assignments.tsv",
        "library\trun_tag\tfolder_name\tfolder_path\n" + "".join(assignment_rows),
    )
    atomic_write(
        project / "rg_metadata.tsv",
        "fnbase\tlibrary\tbp_id\tsample_idx\tlane\tflowcell\tlane_num\tpu\t"
        "rg_id\trg_string\traw_r1\traw_r2\ttrimmed_r1\ttrimmed_r2\n"
        + "".join(rg_rows),
    )
    return sorted(grouped)


def workflow_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate_workflow_compatibility(path: Path, read_format: str) -> None:
    """Fail before submission when the selected RNA workflow cannot do the job."""
    text = path.read_text(encoding="utf-8")
    base_markers = ("params.rg_metadata", "--outSAMattrRGline")
    missing_base = [marker for marker in base_markers if marker not in text]
    if missing_base:
        raise PipelineError(
            f"RNA workflow {path} cannot attach source-FASTQ read groups; "
            "missing: " + ", ".join(missing_base)
        )

    pe150_markers = (
        "rna_geometry",
        "--soloBarcodeMate 1",
        "--clip5pNbases 28 0",
        '"${r1} ${r2}"',
    )
    if read_format == "pe150":
        missing_pe150 = [marker for marker in pe150_markers if marker not in text]
        if missing_pe150:
            raise PipelineError(
                "5' PE150 mapping requires the V10-compatible align_rna.nf; "
                f"{path} is missing: "
                + ", ".join(missing_pe150)
            )

    strict_marker = "Missing RNA read-group metadata for input unit(s)"
    if strict_marker not in text:
        print(
            "WARNING: align_rna.nf supports RG provenance but uses a generic "
            "fallback if an input unit is absent from rg_metadata.tsv. The "
            "generated metadata is expected to be complete; install the V10 "
            "workflow to make this fail closed.",
            file=sys.stderr,
        )


def generate_mapping_script(
    output_base: Path,
    input_dirs: Sequence[Path],
    libs_file: Path,
    chemistry: str,
    read_format: str,
    rna_ref: Path,
    whitelist: Path,
    workflow_file: Path,
    memgb: str,
    threads: int,
    library_numbers: set[int] | None,
    filter_file: Path | None,
    max_cores: int | None,
    lib_prefix: str,
) -> Path:
    project = output_base / "mapping_project"
    project.mkdir(parents=True, exist_ok=True)
    script_path = output_base / "run_mapping.sbatch"
    driver = Path(__file__).resolve()
    workflow_hash = workflow_sha256(workflow_file)
    nextflow_queue_size = (
        max(1, (max_cores - 2) // threads) if max_cores is not None else 50
    )
    consolidate_command = [
        "python3",
        str(driver),
        "--consolidate",
        "--input-dirs",
        *[str(path) for path in input_dirs],
        "--chemistry",
        chemistry,
        "--read-format",
        read_format,
        "--output-base",
        str(output_base),
        "--rna-ref",
        str(rna_ref),
        "--rna-whitelist",
        str(whitelist),
        "--workflow-file",
        str(workflow_file),
        "--lib-prefix",
        lib_prefix,
    ]
    if library_numbers:
        consolidate_command.extend(["--libraries", *[str(value) for value in sorted(library_numbers)]])
    if filter_file:
        consolidate_command.extend(["--filter-libs", str(filter_file)])
    consolidate_text = (" " + chr(92) + "\n  ").join(
        q(part) for part in consolidate_command
    )

    script = f"""#!/bin/bash
#SBATCH --job-name=rna_{read_format}_mapping
#SBATCH --output={project}/mapping_%j.out
#SBATCH --error={project}/mapping_%j.err
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=7-00:00:00
#SBATCH --exclude=squirtle
#SBATCH --chdir={project}

set -euo pipefail
module purge
module load miniforge/3 nextflow/latest align_pipelines/bjp

PROJECT={q(project)}
OUTPUT_BASE={q(output_base)}
WORKFLOW_SHA256={q(workflow_hash)}
mkdir -p \"${{PROJECT}}/logs\"
cd \"${{PROJECT}}\"

{consolidate_text}
CONSOLIDATED_DIR=\"${{OUTPUT_BASE}}/mapping_input\"
shopt -s nullglob
R1_FILES=(\"${{CONSOLIDATED_DIR}}\"/*_R1_*.fastq.gz)
R2_FILES=(\"${{CONSOLIDATED_DIR}}\"/*_R2_*.fastq.gz)
(( ${{#R1_FILES[@]}} > 0 && ${{#R1_FILES[@]}} == ${{#R2_FILES[@]}} )) || {{
  echo \"ERROR: consolidated R1/R2 counts are missing or unequal\" >&2
  exit 1
}}

[[ -n \"${{ALIGN_PIPELINES_HOME:-}}\" && -f \"${{ALIGN_PIPELINES_HOME}}/align_pipelines.nf\" ]] || {{
  echo \"ERROR: align_pipelines/bjp did not provide ALIGN_PIPELINES_HOME\" >&2
  exit 1
}}
MODULE_WORKFLOW=\"${{ALIGN_PIPELINES_HOME}}/workflows/align_rna.nf\"
[[ -f \"${{MODULE_WORKFLOW}}\" ]] || {{
  echo \"ERROR: align_pipelines/bjp RNA workflow is missing: ${{MODULE_WORKFLOW}}\" >&2
  exit 1
}}
if [[ \"$(sha256sum \"${{MODULE_WORKFLOW}}\" | awk '{{print $1}}')\" != \"${{WORKFLOW_SHA256}}\" ]]; then
  echo \"WARNING: installed align_pipelines/bjp align_rna.nf differs from the requested workflow; continuing\" >&2
fi
cat > params_rna.yml <<EOF
libs: '{libs_file}'
output_directory: '{output_base}/mapping_output'
memgb: '{memgb}'
threads: {threads}
rna_dir: '${{CONSOLIDATED_DIR}}'
rna_ref: '{rna_ref}'
rna_whitelist: '{whitelist}'
rg_metadata: '{project}/rg_metadata.tsv'
rna_geometry: '{read_format}'
EOF

cat > nextflow.config <<'NFCONFIG'
process {{
    executor = 'slurm'
    queue = 'compute'
    clusterOptions = '--exclude=squirtle'
    beforeScript = 'module purge && module load miniforge/3 nextflow/latest align_pipelines/bjp && module load htslib/1.20 samtools/1.20 star/2.7.11b && command -v STAR >/dev/null && command -v samtools >/dev/null'
    cpus = 4
    memory = '32 GB'
    time = '168 h'
    withName: 'map_rna' {{
        cpus = {threads}
        memory = '{memgb} GB'
        time = '168 h'
    }}
}}
workDir = './work'
resume = true
executor {{
    name = 'slurm'
    queueSize = {nextflow_queue_size}
    submitRateLimit = '20/1min'
    pollInterval = '30 sec'
}}
NFCONFIG

rm -f .nextflow/cache/*/db/LOCK 2>/dev/null || true
nextflow run \"${{ALIGN_PIPELINES_HOME}}/align_pipelines.nf\" \\
  -params-file params_rna.yml \\
  -resume \\
  -with-report report.html \\
  -with-trace trace.txt \\
  -with-timeline timeline.html \\
  -ansi-log false

echo \"Results: {output_base}/mapping_output\"
"""
    script_path.write_text(script, encoding="utf-8")
    script_path.chmod(0o755)
    return script_path


def read_filter_names(path: Path | None) -> set[str] | None:
    if path is None:
        return None
    return {line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Generate 10X RNA trim/map jobs with RG provenance and 5' format routing"
    )
    parser.add_argument("--input-dirs", nargs="+", default=None)
    parser.add_argument("--chemistry", choices=("3prime", "5prime"), default=None)
    parser.add_argument("--read-format", choices=("auto", "long-r2", "pe150"), default="auto")
    parser.add_argument("--output-base", default=None)
    parser.add_argument("--rna-ref", default=None)
    parser.add_argument("--rna-whitelist", default=None)
    parser.add_argument(
        "--workflow-file",
        default=DEFAULT_RNA_WORKFLOW,
        help=(
            "RNA workflow in this align_pipelines repository/package; default: "
            f"{DEFAULT_RNA_WORKFLOW}"
        ),
    )
    parser.add_argument("--filter-libs", default=None)
    parser.add_argument("--libraries", nargs="+", type=int, default=None)
    parser.add_argument(
        "--lib-prefix",
        default=None,
        help=(
            "Prefix before the numeric library ID. Defaults to "
            f"{DEFAULT_RNA3_LIB_PREFIX!r} for 3prime and "
            f"{DEFAULT_RNA5_LIB_PREFIX!r} for 5prime."
        ),
    )
    parser.add_argument("--memgb", default="80")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument(
        "--max-cores",
        type=int,
        default=None,
        help=(
            "Hard RNA core ceiling used to throttle trim arrays and concurrent "
            "Nextflow mapping workers"
        ),
    )
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--no-info-file", action="store_true")
    parser.add_argument("--skip-trimming", action="store_true")
    parser.add_argument("--skip-mapping", action="store_true")
    parser.add_argument("--submit-jobs", action="store_true")
    parser.add_argument("--consolidate", action="store_true")
    parser.add_argument("--min-pe150-tso-match-fraction", type=float, default=0.80)
    parser.add_argument("--prepare-5p-pe150-r1", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--worker-input", default=None, help=argparse.SUPPRESS)
    parser.add_argument("--worker-output", default=None, help=argparse.SUPPRESS)
    parser.add_argument("--worker-tso-info", default=None, help=argparse.SUPPRESS)
    return parser


def normal_main(args: argparse.Namespace) -> int:
    required = {
        "--input-dirs": args.input_dirs,
        "--chemistry": args.chemistry,
        "--output-base": args.output_base,
        "--rna-ref": args.rna_ref,
        "--rna-whitelist": args.rna_whitelist,
    }
    missing = [name for name, value in required.items() if not value]
    if missing:
        raise PipelineError("missing required argument(s): " + ", ".join(missing))
    if not 0 <= args.min_pe150_tso_match_fraction <= 1:
        raise PipelineError("--min-pe150-tso-match-fraction must be between 0 and 1")
    if args.max_cores is not None:
        if args.max_cores < 3:
            raise PipelineError("--max-cores must be at least 3")
        if not args.skip_mapping and args.threads + 2 > args.max_cores:
            raise PipelineError(
                "--max-cores must accommodate one mapping worker plus the "
                f"2-core Nextflow controller ({args.threads + 2} cores required)"
            )
    input_dirs = [Path(value).expanduser().resolve() for value in args.input_dirs]
    output_base = Path(args.output_base).expanduser().resolve(strict=False)
    rna_ref = Path(args.rna_ref).expanduser().resolve()
    whitelist = Path(args.rna_whitelist).expanduser().resolve()
    workflow = Path(args.workflow_file).expanduser().resolve()
    filter_file = Path(args.filter_libs).expanduser().resolve() if args.filter_libs else None
    for path in input_dirs:
        if not path.is_dir():
            raise PipelineError(f"input directory does not exist: {path}")
    if not rna_ref.is_dir():
        raise PipelineError(f"STAR reference directory does not exist: {rna_ref}")
    for label, path in (("RNA whitelist", whitelist), ("align_rna.nf", workflow)):
        if not path.is_file():
            raise PipelineError(f"{label} does not exist: {path}")
    if args.dry_run and not args.skip_mapping:
        raise PipelineError("--dry-run produces no trimmed FASTQs; combine it with --skip-mapping")
    output_base.mkdir(parents=True, exist_ok=True)
    library_numbers = set(args.libraries) if args.libraries else None
    filter_names = read_filter_names(filter_file)
    lib_prefix = args.lib_prefix or (
        DEFAULT_RNA3_LIB_PREFIX
        if args.chemistry != "5prime"
        else DEFAULT_RNA5_LIB_PREFIX
    )

    formats: set[str] = set()
    all_libraries: set[str] = set()
    active_input_dirs: list[Path] = []
    for input_dir in input_dirs:
        all_pairs = fastq_pairs(input_dir)
        if not all_pairs:
            raise PipelineError(f"no R1/R2 pairs found in {input_dir}")
        pairs = [
            pair
            for pair in all_pairs
            if selected_library(
                library_name(pair[0], lib_prefix),
                library_numbers,
                filter_names,
                lib_prefix,
            )
        ]
        if not pairs:
            print(f"Skipping {input_dir.name}: no selected libraries are present")
            continue
        active_input_dirs.append(input_dir)
        observed = validate_run_format(pairs, args.chemistry, args.read_format)
        formats.add(observed)
        for r1, _ in pairs:
            all_libraries.add(library_name(r1, lib_prefix))
    if not all_libraries:
        raise PipelineError("no libraries remained after filters")
    if len(formats) != 1:
        raise PipelineError(
            "one driver invocation must contain one read geometry; group long-r2 and PE150 runs separately"
        )
    read_format = next(iter(formats))
    validate_workflow_compatibility(workflow, read_format)
    libs_file = output_base / "libs.txt"
    atomic_write(libs_file, "".join(f"{lib}\n" for lib in sorted(all_libraries)))

    trim_scripts: list[Path] = []
    if not args.skip_trimming:
        for input_dir in active_input_dirs:
            trim_scripts.append(
                generate_trim_script(
                    input_dir,
                    output_base,
                    args.chemistry,
                    read_format,
                    args.no_info_file,
                    args.dry_run,
                    args.min_pe150_tso_match_fraction,
                    args.max_cores,
                    library_numbers,
                    filter_names,
                    lib_prefix,
                )
            )
    mapping_script: Path | None = None
    if not args.skip_mapping:
        mapping_script = generate_mapping_script(
            output_base,
            active_input_dirs,
            libs_file,
            args.chemistry,
            read_format,
            rna_ref,
            whitelist,
            workflow,
            args.memgb,
            args.threads,
            library_numbers,
            filter_file,
            args.max_cores,
            lib_prefix,
        )

    print(f"RNA driver release: {RELEASE}")
    print(f"Detected geometry: {read_format}")
    for script in trim_scripts:
        print(f"Generated trim job: {script}")
    if mapping_script:
        print(f"Generated mapping job: {mapping_script}")
    if args.submit_jobs:
        trim_ids: list[str] = []
        previous_trim_id: str | None = None
        for script in trim_scripts:
            command = ["sbatch", "--parsable"]
            if args.max_cores is not None and previous_trim_id is not None:
                command.extend(["--dependency", f"afterok:{previous_trim_id}"])
            command.append(str(script))
            result = subprocess.run(
                command, check=True, text=True, capture_output=True
            )
            previous_trim_id = result.stdout.strip().split(";")[0]
            trim_ids.append(previous_trim_id)
        if mapping_script:
            command = ["sbatch", "--parsable"]
            if trim_ids:
                command.extend(["--dependency", "afterok:" + ":".join(trim_ids)])
            command.append(str(mapping_script))
            subprocess.run(command, check=True)
    return 0


def main() -> int:
    args = build_parser().parse_args()
    try:
        if args.prepare_5p_pe150_r1:
            if not (args.worker_input and args.worker_output and args.worker_tso_info):
                raise PipelineError("PE150 worker arguments are incomplete")
            prepare_pe150_r1(
                Path(args.worker_input),
                Path(args.worker_output),
                Path(args.worker_tso_info),
                args.min_pe150_tso_match_fraction,
            )
            return 0
        if args.consolidate:
            required = (args.input_dirs, args.output_base)
            if not all(required):
                raise PipelineError("--consolidate requires --input-dirs and --output-base")
            lib_prefix = args.lib_prefix or (
                DEFAULT_RNA3_LIB_PREFIX
                if args.chemistry != "5prime"
                else DEFAULT_RNA5_LIB_PREFIX
            )
            consolidate(
                [Path(value).resolve() for value in args.input_dirs],
                Path(args.output_base).resolve(),
                set(args.libraries) if args.libraries else None,
                read_filter_names(Path(args.filter_libs).resolve() if args.filter_libs else None),
                lib_prefix,
            )
            return 0
        return normal_main(args)
    except (PipelineError, OSError, subprocess.CalledProcessError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
