#!/usr/bin/env python3
"""Stage, submit, and validate Tet2025 10X trimming/mapping/QC jobs.

This is a control layer around the historical RNA and ATAC drivers plus the
bundled provenance and 5' corrections.  It does not reimplement STARsolo,
minimap2, Nextflow, or the plotting code.  Its responsibilities are
deliberately narrow:

* resolve short BP run names below the known raw-data roots;
* reject accidental reuse of an output run directory;
* generate all child scripts below one run-specific staging directory;
* submit the SLURM dependency graph only when ``--submit`` is present;
* collect trimming, RNA mapping, and ATAC QC plots;
* validate durable outputs without publishing over production data;
* mark Nextflow work caches as manually removable, but never delete them.

Dry-run script generation is the default.  The same command with ``--submit``
launches the jobs.  Use ``--resume`` only to continue the same immutable run.

Repository and package layout
-----------------------------
Keep this file and the runtime helpers listed below together in
``scripts/mapping`` in the align_pipelines repository. Install that directory
into ``bin/mapping`` below the deployed align_pipelines package. The
orchestrator resolves helpers beside its own path and workflows from the
repository/package root, so no resource override is needed in either layout.

* integrated_RNA_trim_map_pipeline_V10.py
* integrated_atac_map_pipeline_V4.py
* trim_barcode_aggregator.py
* analyze_trimming_V15.py
* plot_mapping_stats_V5.py
* collect_atac_qc.py
* plot_atac_qc_v2.py
* run_rna_bam_evidence.py
* rna_evidence_common.py
* analyze_gained_cells.py
* analyze_depth_history.py
* join_trim_cell_metrics.py
* export_repooling_evidence.py
* plot_rna_evidence.py

The Nextflow workflows are owned by the same repository and remain below its
top-level ``workflows`` directory.

Supported chemistry status
--------------------------
* 10X 3' Multiome RNA: trimming + STARsolo + plots.
* 10X Multiome ATAC: mapping + RNA-barcode-aware QC + plots.
* 10X 5' long-R2 and PE150: auto-detected, trimmed, mapped with their correct
  STARsolo geometry, and plotted separately by default.
* Cross-format ``together`` mode is allowed only when the library sets are
  disjoint.  Overlapping libraries are blocked because a BAM concatenation
  does not re-deduplicate molecules across STARsolo invocations.

Example (the four newly sequenced 3' runs)
------------------------------------------
python3 orchestrate_10x_mapping_qc.py \
  --run-name rna3_library19_2026_08_22 \
  --rna3-runs BP16462 BP16480 BP16481 BP16482 \
  --libraries 19 \
  --max-cores 24 \
  --submit

Example (RNA + ATAC; BP names are resolved below their default roots)
---------------------------------------------------------------------
python3 orchestrate_10x_mapping_qc.py \
  --run-name multiome_refresh_2026_08_22 \
  --rna3-runs BP16462 BP16480 BP16481 BP16482 \
  --atac-runs BP_EXAMPLE_1 BP_EXAMPLE_2 \
  --submit

Example (all 5' runs; correct geometry is detected automatically)
------------------------------------------------------------------
python3 orchestrate_10x_mapping_qc.py \
  --run-name rna5_complete_2026_08_22 \
  --rna5-runs BP12954 BP12955 BP14541 BP14542 \
  --rna5-map-mode separate \
  --submit
"""

from __future__ import annotations

import argparse
import csv
import fcntl
import gzip
import hashlib
import importlib.util
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Sequence


RELEASE = "2026-09-15-v25-central-figures-flexible-run-dir"

DEFAULT_RNA3_RAW_ROOT = "/mnt/beegfs/reads/3P_Multiome_10XRNA"
DEFAULT_ATAC_RAW_ROOT = "/mnt/beegfs/reads/10X_ATAC_multiome"
DEFAULT_RNA5_RAW_ROOT = "/mnt/beegfs/reads/5P_10XRNA"
DEFAULT_STAGING_ROOT = "/mnt/beegfs/tetraploid_multiome_cis_trans"
DEFAULT_RNA3_LIB_PREFIX = "Tet_2025_Multiome-RNA_"
DEFAULT_RNA5_LIB_PREFIX = "Tet_2025_RNA_5P_"
DEFAULT_ATAC_LIB_PREFIX = "Tet_2025_Multiome-ATAC_"
MAPPING_SCRIPT_DIR = Path(__file__).resolve().parent
ALIGN_PIPELINES_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_RNA_WORKFLOW = str(ALIGN_PIPELINES_ROOT / "workflows" / "align_rna.nf")
DEFAULT_ATAC_WORKFLOW = str(ALIGN_PIPELINES_ROOT / "workflows" / "align_atac.nf")

DEFAULT_RNA_REF = (
    "/mnt/beegfs/genomes_annotations/ancestral_genomes/litterbox/"
    "human_chimp_bonobo/STARv2.7.11b"
)
DEFAULT_ATAC_REF = (
    "/mnt/beegfs/genomes_annotations/ancestral_genomes/litterbox/"
    "human_chimp_bonobo/human_chimp_bonobo.mm2"
)
DEFAULT_WHITELIST_ROOT = (
    "/mnt/beegfs/genomes_annotations/white_lists_adapters/cellranger_10x"
)
DEFAULT_RNA3_WHITELIST = f"{DEFAULT_WHITELIST_ROOT}/RNA-737K-arc-v1.txt.gz"
DEFAULT_RNA5_WHITELIST = DEFAULT_RNA3_WHITELIST
DEFAULT_ATAC_WHITELIST = f"{DEFAULT_WHITELIST_ROOT}/ATAC-737K-arc-v1.txt.gz"
DEFAULT_PRODUCTION_RNA = (
    "/mnt/beegfs/tetraploid_multiome_cis_trans/3P/mapping_output"
)
DEFAULT_PRODUCTION_ATAC = (
    "/mnt/beegfs/tetraploid_multiome_cis_trans/ATAC/mapping_output"
)
DEFAULT_RNA3_FIGURE_BASE = (
    "/mnt/beegfs/tetraploid_multiome_cis_trans/3P/figures"
)
DEFAULT_RNA5_FIGURE_BASE = (
    "/mnt/beegfs/tetraploid_multiome_cis_trans/5P/figures"
)
DEFAULT_ATAC_FIGURE_BASE = (
    "/mnt/beegfs/tetraploid_multiome_cis_trans/ATAC/figures"
)

ALL_STAGES = ("trim", "map", "qc", "plot", "validate")
DEFAULT_BAM_EVIDENCE_HASH_SEED = 1469598103934665603
RNA_BAM_EVIDENCE_STARSOLO_FEATURE = "GeneFull_Ex50pAS"
RNA_BAM_EVIDENCE_STARSOLO_UMI_FILTERING = "MultiGeneUMI_CR"
RNA_BAM_EVIDENCE_STARSOLO_UMI_DEDUP = "1MM_CR"
RNA_BAM_EVIDENCE_STARSOLO_MULTIMAPPERS = "EM"
RNA_BAM_EVIDENCE_HASH_ALGORITHM = "fnv1a64_seeded_v1"
RNA_BAM_EVIDENCE_ORDINARY_COUNTEDU_READ_DEFINITION = (
    "distinct_RG_QNAME_declared_RG_valid_CB_singleton_feature_GX_"
    "NH_unrestricted_UB_not_required_STARsolo_pre_UMI_filter_countedU"
)
RNA_BAM_EVIDENCE_ORDINARY_MOLECULE_DEFINITION = (
    "distinct_CB_GX_valid_STARsolo_corrected_UB_after_1MM_CR_and_"
    "MultiGeneUMI_CR_NH_unrestricted"
)
RNA_BAM_EVIDENCE_MULTIMAPPER_DEFINITION = (
    "STARsolo_multi_gene_EM_unavailable_from_standard_uppercase_GX_UB_"
    "BAM_tags_NH_is_not_EM_membership"
)
RNA_BAM_EVIDENCE_SUMMARY_UNIQUE_READ_METRIC = (
    "Unique Reads in Cells Mapped to GeneFull_Ex50pAS"
)
RNA_BAM_EVIDENCE_NH_GT1_UNIQUE_GENE_DEFINITION = (
    "subset_of_ordinary_countedU_reads_with_NH_gt1_and_singleton_feature_GX"
)
RNA_BAM_EVIDENCE_EM_AVAILABILITY = (
    "unavailable_from_standard_uppercase_GX_UB_BAM_tags"
)
class OrchestratorError(RuntimeError):
    """A user-facing orchestration error."""


@dataclass
class Resources:
    rna_driver: Path | None
    atac_driver: Path | None
    rna_workflow: Path | None
    atac_workflow: Path | None
    barcode_aggregator: Path | None
    trim_plotter: Path | None
    mapping_plotter: Path | None
    atac_collector: Path | None
    atac_plotter: Path | None
    bam_evidence_runner: Path | None
    bam_evidence_profiler: Path | None
    star_diagnostic_promoter: Path | None


@dataclass
class Inputs:
    rna3: list[Path] = field(default_factory=list)
    atac: list[Path] = field(default_factory=list)
    rna5: list[Path] = field(default_factory=list)


@dataclass
class Discovery:
    libraries: dict[str, list[str]] = field(default_factory=dict)
    library_numbers: dict[str, list[int]] = field(default_factory=dict)
    fastq_units: dict[str, int] = field(default_factory=dict)
    read_lengths: dict[str, list[dict[str, object]]] = field(default_factory=dict)
    fastq_files: dict[str, list[dict[str, object]]] = field(default_factory=dict)
    rna5_formats: dict[str, str] = field(default_factory=dict)
    rna5_libraries_by_format: dict[str, list[str]] = field(default_factory=dict)
    active_runs: dict[str, list[Path]] = field(default_factory=dict)


@dataclass
class JobSpec:
    label: str
    script: Path
    dependencies: list[str] = field(default_factory=list)
    job_id: str | None = None


def q(value: os.PathLike[str] | str) -> str:
    return shlex.quote(str(value))


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def atomic_text(path: Path, text: str, mode: int | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    temp.write_text(text, encoding="utf-8")
    if mode is not None:
        temp.chmod(mode)
    os.replace(temp, path)


def atomic_json(path: Path, payload: object) -> None:
    atomic_text(path, json.dumps(payload, indent=2, sort_keys=True) + "\n")


def acquire_run_lock(run_dir: Path):
    """Hold one prepare/generate/submit transaction per staged run."""
    run_dir.parent.mkdir(parents=True, exist_ok=True)
    lock_path = run_dir.parent / f".{run_dir.name}.orchestrator.lock"
    handle = lock_path.open("a+", encoding="utf-8")
    try:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError as exc:
        handle.close()
        raise OrchestratorError(
            f"another orchestrator process is active for {run_dir}; no jobs "
            "were generated or submitted"
        ) from exc
    return handle


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_stages(raw: Sequence[str]) -> set[str]:
    values: list[str] = []
    for item in raw:
        values.extend(part.strip().lower() for part in item.split(","))
    values = [value for value in values if value]
    if "all" in values:
        return set(ALL_STAGES)
    bad = sorted(set(values) - set(ALL_STAGES))
    if bad:
        raise OrchestratorError(
            f"unknown stage(s): {', '.join(bad)}; valid: all, "
            + ", ".join(ALL_STAGES)
        )
    return set(values)


def safe_run_name(value: str) -> str:
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,95}", value):
        raise OrchestratorError(
            "--run-name must be 1-96 characters and contain only letters, "
            "numbers, dot, underscore, or dash"
        )
    return value


def safe_nodelist(value: str) -> str:
    """Validate a SLURM host-list expression before embedding it in scripts."""
    if not re.fullmatch(r"[A-Za-z0-9_.\-\[\],]+", value):
        raise OrchestratorError(
            "--nodelist must be a SLURM host-list expression containing only "
            "letters, numbers, dot, underscore, dash, brackets, or commas"
        )
    if "squirtle" in value.lower():
        raise OrchestratorError(
            "--nodelist cannot include squirtle because the mapping workflows "
            "explicitly exclude that node"
        )
    return value


def resolve_runs(values: Sequence[str] | None, root: str) -> list[Path]:
    resolved: list[Path] = []
    for value in values or []:
        candidate = Path(value).expanduser()
        if not candidate.is_absolute() and len(candidate.parts) == 1:
            candidate = Path(root).expanduser() / candidate
        resolved.append(candidate.resolve(strict=False))
    return resolved


def ensure_unique_run_basenames(label: str, runs: Sequence[Path]) -> list[str]:
    names = [path.name for path in runs]
    duplicates = sorted({name for name in names if names.count(name) > 1})
    if duplicates:
        return [
            f"{label}: duplicate run-folder basename(s) would collide in the "
            f"child driver: {', '.join(duplicates)}"
        ]
    return []


def find_resource(
    explicit: str | None,
    resource_root: Path,
    relative_candidates: Sequence[str],
    exact_names: Sequence[str],
) -> Path | None:
    if explicit:
        return Path(explicit).expanduser().resolve(strict=False)
    for relative in relative_candidates:
        candidate = resource_root / relative
        if candidate.is_file():
            return candidate.resolve()
    for name in exact_names:
        matches = sorted(resource_root.rglob(name)) if resource_root.is_dir() else []
        if matches:
            return matches[0].resolve()
    return None


def find_bam_evidence_profiler(explicit: str | None) -> Path | None:
    """Resolve the compiled profiler without depending on a module-mutated PATH."""
    if explicit:
        return Path(explicit).expanduser().resolve(strict=False)
    candidates = (
        MAPPING_SCRIPT_DIR.parent / "rna_bam_evidence",
        ALIGN_PIPELINES_ROOT / "rna_bam_evidence",
    )
    for candidate in candidates:
        if candidate.is_file():
            return candidate.resolve()
    return None


def resolve_resources(args: argparse.Namespace) -> Resources:
    root = Path(args.resource_root).expanduser().resolve(strict=False)
    return Resources(
        rna_driver=find_resource(
            args.rna_driver,
            root,
            (
                "integrated_RNA_trim_map_pipeline_V10.py",
            ),
            ("integrated_RNA_trim_map_pipeline_V10.py",),
        ),
        atac_driver=find_resource(
            args.atac_driver,
            root,
            (
                "integrated_atac_map_pipeline_V4.py",
                "ATAC_Map/integrated_atac_map_pipeline_V4.py",
            ),
            ("integrated_atac_map_pipeline_V4.py",),
        ),
        rna_workflow=find_resource(
            args.rna_workflow,
            root,
            (),
            (),
        ),
        atac_workflow=find_resource(
            args.atac_workflow,
            root,
            (),
            (),
        ),
        barcode_aggregator=find_resource(
            args.barcode_aggregator,
            root,
            ("trim_barcode_aggregator.py", "Trimming_QC/trim_barcode_aggregator.py"),
            ("trim_barcode_aggregator.py",),
        ),
        trim_plotter=find_resource(
            args.trim_plotter,
            root,
            (
                "analyze_trimming_V15.py",
                "Trimming_QC/analyze_trimming_V15.py",
            ),
            ("analyze_trimming_V15.py", "analyze_trimming_V12.py"),
        ),
        mapping_plotter=find_resource(
            args.mapping_plotter,
            root,
            (
                "plot_mapping_stats_V5.py",
                "Mapping_QC/plot_mapping_stats_V5.py",
            ),
            ("plot_mapping_stats_V5.py",),
        ),
        atac_collector=find_resource(
            args.atac_collector,
            root,
            ("collect_atac_qc.py", "ATAC_QC/collect_atac_qc.py"),
            ("collect_atac_qc.py",),
        ),
        atac_plotter=find_resource(
            args.atac_plotter,
            root,
            ("plot_atac_qc_v2.py", "ATAC_QC/plot_atac_qc_v2.py"),
            ("plot_atac_qc_v2.py",),
        ),
        bam_evidence_runner=find_resource(
            args.rna3_bam_evidence_runner,
            root,
            ("run_rna_bam_evidence.py",),
            ("run_rna_bam_evidence.py",),
        ),
        bam_evidence_profiler=find_bam_evidence_profiler(
            args.rna3_bam_evidence_profiler
        ),
        star_diagnostic_promoter=find_resource(
            args.rna3_star_diagnostic_promoter,
            root,
            ("promote_rna_star_diagnostics.py",),
            ("promote_rna_star_diagnostics.py",),
        ),
    )


def fastq_pairs(run: Path) -> list[tuple[Path, Path]]:
    pairs: set[tuple[Path, Path]] = set()
    for pattern in ("*_R1_*.fastq.gz", "*_R1_*.fq.gz"):
        for r1 in run.glob(pattern):
            r2 = Path(str(r1).replace("_R1_", "_R2_"))
            if r2.is_file():
                pairs.add((r1, r2))
    return sorted(pairs)


def fastq_triplets(run: Path) -> list[tuple[Path, Path, Path]]:
    triplets: set[tuple[Path, Path, Path]] = set()
    for pattern in ("*_R1_*.fastq.gz", "*_R1_*.fq.gz"):
        for r1 in run.glob(pattern):
            r2 = Path(str(r1).replace("_R1_", "_R2_"))
            r3 = Path(str(r1).replace("_R1_", "_R3_"))
            if r2.is_file() and r3.is_file():
                triplets.add((r1, r2, r3))
    return sorted(triplets)


def library_from_fastq(path: Path, prefix: str | None = None) -> str:
    """Return a canonical library name from legacy or annotated FASTQ names.

    New bcl-convert names may insert annotations such as ``_L1_L2`` between
    the actual library ID and the ``_S#_L###`` sample/lane fields.  When a
    prefix is supplied, ``<prefix><number>`` is authoritative.  The generic
    fallback removes trailing ``_L<number>`` annotations before the sample
    token so explicit-name filtering remains usable for other prefixes.
    """
    match = re.match(r"(.+?)_S\d+_L\d+", path.name)
    candidate = match.group(1) if match else re.sub(r"_R[123]_.*$", "", path.name)
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


def first_sequence_length(path: Path) -> int:
    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        header = handle.readline()
        sequence = handle.readline().strip()
    if not header.startswith("@") or not sequence:
        raise OrchestratorError(f"could not read the first FASTQ record from {path}")
    return len(sequence)


def detect_five_prime_format(r1_length: int, r2_length: int) -> str:
    if r1_length <= 35 and r2_length >= 80:
        return "long-r2"
    if r1_length >= 100 and r2_length >= 100:
        return "pe150"
    return "unknown"


def discover(
    inputs: Inputs,
    selected_library_numbers: set[int] | None = None,
    library_prefixes: dict[str, str] | None = None,
) -> tuple[Discovery, list[str], list[str]]:
    result = Discovery()
    failures: list[str] = []
    warnings: list[str] = []

    modalities: list[tuple[str, Sequence[Path], bool]] = [
        ("rna3", inputs.rna3, False),
        ("atac", inputs.atac, True),
        ("rna5", inputs.rna5, False),
    ]
    prefixes = library_prefixes or {
        "rna3": DEFAULT_RNA3_LIB_PREFIX,
        "rna5": DEFAULT_RNA5_LIB_PREFIX,
        "atac": DEFAULT_ATAC_LIB_PREFIX,
    }
    expected_patterns = {
        modality: re.compile(rf"^{re.escape(prefix)}\d+$")
        for modality, prefix in prefixes.items()
    }

    for modality, runs, is_atac in modalities:
        prefix = prefixes[modality]
        libraries: set[str] = set()
        lengths: list[dict[str, object]] = []
        file_records: list[dict[str, object]] = []
        active_runs: list[Path] = []
        unit_count = 0
        failures.extend(ensure_unique_run_basenames(modality, runs))
        for run in runs:
            if not run.is_dir():
                failures.append(f"{modality}: input directory does not exist: {run}")
                continue
            all_units = fastq_triplets(run) if is_atac else fastq_pairs(run)
            if not all_units:
                kind = "R1/R2/R3 triplets" if is_atac else "R1/R2 pairs"
                failures.append(f"{modality}: no complete {kind} found in {run}")
                continue
            units = [
                unit
                for unit in all_units
                if selected_library_numbers is None
                or library_number(library_from_fastq(unit[0], prefix), prefix)
                in selected_library_numbers
            ]
            if not units:
                warnings.append(
                    f"{modality}: none of the selected libraries are present in "
                    f"{run.name}; this run will be skipped"
                )
                continue
            active_runs.append(run)
            unit_count += len(units)
            libraries.update(library_from_fastq(unit[0], prefix) for unit in units)
            for fastq in sorted({path for unit in units for path in unit}):
                try:
                    stat = fastq.stat()
                    file_records.append(
                        {
                            "path": str(fastq.resolve()),
                            "size": stat.st_size,
                            "mtime_ns": stat.st_mtime_ns,
                        }
                    )
                except OSError as exc:
                    failures.append(f"{modality}: could not stat FASTQ {fastq}: {exc}")
            try:
                r1_length = first_sequence_length(units[0][0])
                r2_length = first_sequence_length(units[0][1])
                record: dict[str, object] = {
                    "run": str(run),
                    "r1": r1_length,
                    "r2": r2_length,
                }
                if is_atac:
                    record["r3"] = first_sequence_length(units[0][2])
                lengths.append(record)
                if modality == "rna5":
                    unit_formats: set[str] = set()
                    for unit in units:
                        unit_r1 = first_sequence_length(unit[0])
                        unit_r2 = first_sequence_length(unit[1])
                        unit_format = detect_five_prime_format(unit_r1, unit_r2)
                        unit_formats.add(unit_format)
                        lib = library_from_fastq(unit[0], prefix)
                        current = set(result.rna5_libraries_by_format.get(unit_format, []))
                        current.add(lib)
                        result.rna5_libraries_by_format[unit_format] = sorted(current)
                    result.rna5_formats[str(run)] = (
                        next(iter(unit_formats)) if len(unit_formats) == 1 else "mixed"
                    )
            except (OSError, EOFError, OrchestratorError) as exc:
                failures.append(f"{modality}: FASTQ sampling failed in {run}: {exc}")

        unexpected = sorted(
            name for name in libraries if not expected_patterns[modality].match(name)
        )
        if unexpected:
            warnings.append(
                f"{modality}: nonstandard library name(s): {', '.join(unexpected[:8])}"
            )
        result.libraries[modality] = sorted(libraries)
        result.library_numbers[modality] = sorted(
            number
            for number in (library_number(name, prefix) for name in libraries)
            if number is not None
        )
        result.fastq_units[modality] = unit_count
        result.read_lengths[modality] = lengths
        result.fastq_files[modality] = sorted(
            file_records, key=lambda record: str(record["path"])
        )
        result.active_runs[modality] = active_runs
        if runs and selected_library_numbers is not None and not active_runs:
            failures.append(
                f"{modality}: none of the requested libraries "
                f"{sorted(selected_library_numbers)} were found in the supplied runs"
            )

    return result, failures, warnings


def resource_failures(
    inputs: Inputs, stages: set[str], resources: Resources, args: argparse.Namespace
) -> list[str]:
    failures: list[str] = []

    def require_file(label: str, path: Path | None) -> None:
        if path is None or not path.is_file():
            failures.append(f"required {label} was not found: {path or '[unresolved]'}")

    def require_path(label: str, value: str, directory: bool = False) -> None:
        path = Path(value)
        good = path.is_dir() if directory else path.is_file()
        if not good:
            kind = "directory" if directory else "file"
            failures.append(f"required {label} {kind} was not found: {path}")

    rna_work = bool(inputs.rna3 or inputs.rna5) and bool({"trim", "map"} & stages)
    if rna_work:
        require_file("RNA driver", resources.rna_driver)
        require_file("RNA Nextflow workflow", resources.rna_workflow)
        require_path("RNA reference", args.rna_ref, directory=True)
    if inputs.rna3 and rna_work:
        require_path("3' RNA whitelist", args.rna3_whitelist)
    if inputs.rna5 and rna_work:
        require_path("5' RNA whitelist", args.rna5_whitelist)

    if inputs.atac and "map" in stages:
        require_file("ATAC driver", resources.atac_driver)
        require_file("ATAC Nextflow workflow", resources.atac_workflow)
        require_path("ATAC minimap2 reference", args.atac_ref)
        require_path("RNA Multiome whitelist", args.rna3_whitelist)
        require_path("ATAC Multiome whitelist", args.atac_whitelist)

    if "plot" in stages:
        if inputs.rna3 or inputs.rna5:
            require_file("trimming plotter", resources.trim_plotter)
            require_file("mapping plotter", resources.mapping_plotter)
        if args.rna3_mapping_baseline_stats:
            require_path(
                "3' RNA mapping baseline statistics",
                args.rna3_mapping_baseline_stats,
            )
        if args.rna3_mapping_baseline_root:
            require_path(
                "3' RNA mapping baseline root",
                args.rna3_mapping_baseline_root,
                directory=True,
            )
        if args.rna3_mapping_baseline_work_root:
            require_path(
                "3' RNA mapping baseline Nextflow work root",
                args.rna3_mapping_baseline_work_root,
                directory=True,
            )
        if inputs.atac:
            require_file("ATAC plotter", resources.atac_plotter)
    if inputs.atac and "qc" in stages:
        require_file("ATAC QC collector", resources.atac_collector)
    if (inputs.rna3 or inputs.rna5) and "qc" in stages and not args.no_trim_info:
        require_file("barcode-trimming aggregator", resources.barcode_aggregator)
        if args.empty_drop_roster:
            require_path("explicit empty-drop roster", args.empty_drop_roster)
        if "trim" not in stages or "map" not in stages:
            failures.append(
                "RNA barcode-linked trimming QC requires both trim and map stages "
                "(or use --no-trim-info to skip that QC product)"
            )

    if inputs.rna3 and args.rna3_bam_evidence_from_bam:
        require_file("RNA BAM evidence runner", resources.bam_evidence_runner)
        require_file("compiled RNA BAM evidence profiler", resources.bam_evidence_profiler)
        if (
            resources.bam_evidence_profiler is not None
            and resources.bam_evidence_profiler.is_file()
            and not os.access(resources.bam_evidence_profiler, os.X_OK)
        ):
            failures.append(
                "compiled RNA BAM evidence profiler is not executable: "
                f"{resources.bam_evidence_profiler}"
            )
        require_file(
            "RNA STAR diagnostic promoter", resources.star_diagnostic_promoter
        )
        if args.rna3_bam_evidence_baseline_root:
            require_path(
                "RNA BAM evidence baseline root",
                args.rna3_bam_evidence_baseline_root,
                directory=True,
            )
        if args.rna3_bam_evidence_class_manifest:
            require_path(
                "RNA BAM evidence classification manifest",
                args.rna3_bam_evidence_class_manifest,
            )

    if args.submit and shutil.which("sbatch") is None:
        failures.append("--submit was requested but sbatch is not available in PATH")
    return failures


def protected_output_failures(run_dir: Path) -> list[str]:
    failures: list[str] = []
    candidate = run_dir.resolve(strict=False)
    for protected in (Path(DEFAULT_PRODUCTION_RNA), Path(DEFAULT_PRODUCTION_ATAC)):
        protected = protected.resolve(strict=False)
        if candidate == protected or protected in candidate.parents:
            failures.append(
                f"staging run directory may not be production mapping output: {candidate}"
            )
    return failures


def print_diagnosis(
    run_dir: Path,
    stages: set[str],
    inputs: Inputs,
    resources: Resources,
    discovery: Discovery,
    failures: Sequence[str],
    warnings: Sequence[str],
) -> None:
    print("=" * 72)
    print("10X MAPPING/QC ORCHESTRATOR DIAGNOSTICS")
    print(f"  Release: {RELEASE}")
    print(f"  Staging run: {run_dir}")
    print(f"  Stages: {', '.join(stage for stage in ALL_STAGES if stage in stages)}")
    print("  Production output overwrite: DISABLED")
    print("=" * 72)
    for modality, runs in (
        ("rna3", inputs.rna3),
        ("atac", inputs.atac),
        ("rna5", inputs.rna5),
    ):
        if not runs:
            continue
        print(f"\n{modality}: {len(runs)} run(s), "
              f"{discovery.fastq_units.get(modality, 0)} FASTQ unit(s), "
              f"{len(discovery.libraries.get(modality, []))} library/libraries")
        for record in discovery.read_lengths.get(modality, []):
            detail = f"R1={record['r1']} R2={record['r2']}"
            if "r3" in record:
                detail += f" R3={record['r3']}"
            if modality == "rna5":
                detail += f" format={discovery.rna5_formats.get(str(record['run']))}"
            print(f"  {Path(str(record['run'])).name}: {detail}")
        libs = discovery.libraries.get(modality, [])
        if libs:
            preview = ", ".join(libs[:8])
            if len(libs) > 8:
                preview += f", ... (+{len(libs) - 8})"
            print(f"  libraries: {preview}")

    print("\nResolved resources:")
    for name, value in resources.__dict__.items():
        print(f"  {name}: {value or '[not found]'}")
    if warnings:
        print("\nWarnings:")
        for warning in warnings:
            print(f"  WARNING: {warning}")
    if failures:
        print("\nBlocking problems:")
        for failure in failures:
            print(f"  ERROR: {failure}")
    else:
        print("\nPreflight: PASS")


def config_payload(
    args: argparse.Namespace,
    inputs: Inputs,
    resources: Resources,
    discovery: Discovery,
) -> dict[str, object]:
    return {
        "release": RELEASE,
        "run_name": args.run_name,
        "run_dir": args.run_dir,
        "inputs": {
            "rna3": [str(path) for path in inputs.rna3],
            "atac": [str(path) for path in inputs.atac],
            "rna5": [str(path) for path in inputs.rna5],
        },
        "references": {
            "rna_ref": str(Path(args.rna_ref).resolve(strict=False)),
            "atac_ref": str(Path(args.atac_ref).resolve(strict=False)),
            "rna3_whitelist": str(Path(args.rna3_whitelist).resolve(strict=False)),
            "rna5_whitelist": (
                str(Path(args.rna5_whitelist).resolve(strict=False))
                if args.rna5_whitelist
                else None
            ),
            "atac_whitelist": str(Path(args.atac_whitelist).resolve(strict=False)),
            "empty_drop_roster": (
                {
                    "path": str(Path(args.empty_drop_roster).resolve(strict=False)),
                    "sha256": file_sha256(Path(args.empty_drop_roster)),
                }
                if args.empty_drop_roster
                else None
            ),
        },
        "resources": {
            key: (
                {"path": str(value), "sha256": file_sha256(value)}
                if value and value.is_file()
                else None
            )
            for key, value in resources.__dict__.items()
        },
        "libraries": discovery.libraries,
        "input_fastqs": discovery.fastq_files,
        "rna5_formats": discovery.rna5_formats,
        "stages": [
            stage for stage in ALL_STAGES if stage in parse_stages(args.stages)
        ],
        "options": {
            "libraries": sorted(args.libraries) if args.libraries else None,
            "rna3_lib_prefix": args.rna3_lib_prefix,
            "rna5_lib_prefix": args.rna5_lib_prefix,
            "atac_lib_prefix": args.atac_lib_prefix,
            "rna_mem_gb": args.rna_mem_gb,
            "rna_threads": args.rna_threads,
            "max_cores": args.max_cores,
            "atac_mem_gb": args.atac_mem_gb,
            "atac_threads": args.atac_threads,
            "atac_chunks": args.atac_chunks,
            "no_trim_info": args.no_trim_info,
            "rna5_map_mode": args.rna5_map_mode,
            "min_pe150_tso_match_fraction": args.min_pe150_tso_match_fraction,
            "rna_barcode_base": args.rna_barcode_base,
            "atac_qc_copy_inputs": args.atac_qc_copy_inputs,
        },
        "scheduling": (
            {
                "nodelist": args.nodelist,
                "array_max_concurrent": args.array_max_concurrent,
            }
            if args.nodelist is not None or args.array_max_concurrent is not None
            else None
        ),
        "reporting": {
            "rna3_figure_root": args.rna3_figure_root,
            "rna5_figure_root": args.rna5_figure_root,
            "atac_figure_root": args.atac_figure_root,
            "rna3_mapping_baseline_stats": (
                {
                    "path": str(
                        Path(args.rna3_mapping_baseline_stats).resolve(strict=False)
                    ),
                    "sha256": file_sha256(Path(args.rna3_mapping_baseline_stats)),
                }
                if args.rna3_mapping_baseline_stats
                else None
            ),
            "rna3_mapping_baseline_label": args.rna3_mapping_baseline_label,
            "rna3_mapping_current_label": args.rna3_mapping_current_label,
            "rna3_mapping_baseline_root": args.rna3_mapping_baseline_root,
            "rna3_mapping_baseline_work_root": args.rna3_mapping_baseline_work_root,
            "rna3_bam_evidence_from_bam": args.rna3_bam_evidence_from_bam,
            "rna3_bam_evidence_cpus": args.rna3_bam_evidence_cpus,
            "rna3_bam_evidence_memory_gb": args.rna3_bam_evidence_memory_gb,
            "rna3_bam_evidence_max_concurrent": (
                args.rna3_bam_evidence_max_concurrent
            ),
            "rna3_bam_evidence_baseline_root": (
                args.rna3_bam_evidence_baseline_root
            ),
            "rna3_bam_evidence_source_order": (
                args.rna3_bam_evidence_source_order
            ),
            "rna3_bam_evidence_hash_bins": args.rna3_bam_evidence_hash_bins,
            "rna3_bam_evidence_class_manifest": (
                args.rna3_bam_evidence_class_manifest
            ),
            "rna3_bam_evidence_no_biological_classification": (
                args.rna3_bam_evidence_no_biological_classification
            ),
            "rna3_bam_evidence_reset_failed_run": (
                args.rna3_bam_evidence_reset_failed_run
            ),
        },
    }


def payload_hash(payload: dict[str, object]) -> str:
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def scientific_config(payload: dict[str, object]) -> dict[str, object]:
    """Return immutable mapping inputs; BAM evidence has its own frozen config."""
    return {
        key: value
        for key, value in payload.items()
        if key not in {"release", "resources", "scheduling", "reporting"}
    }


def runtime_config(payload: dict[str, object]) -> dict[str, object]:
    """Return warning-only software provenance and scheduling controls."""
    return {
        "release": payload.get("release"),
        "resources": payload.get("resources"),
        "scheduling": payload.get("scheduling"),
        "reporting": payload.get("reporting"),
    }


def stored_payload(stored: dict[str, object], current: dict[str, object]) -> dict[str, object]:
    """Recover the payload fields from either legacy or current run_config.json."""
    return {key: stored.get(key) for key in current}


def prepare_run_directory(
    run_dir: Path,
    resume: bool,
    payload: dict[str, object],
    requested_resubmits: Sequence[str] | None,
) -> None:
    config_path = run_dir / "control" / "run_config.json"
    if run_dir.exists() and not resume:
        raise OrchestratorError(
            f"run directory already exists: {run_dir}\n"
            "Choose a new --run-name, or pass --resume to continue this exact run."
        )
    if resume:
        if not config_path.is_file():
            raise OrchestratorError(
                f"--resume requires the existing immutable config: {config_path}"
            )
        old = json.loads(config_path.read_text(encoding="utf-8"))
        old_payload = stored_payload(old, payload)
        old_scientific_hash = payload_hash(scientific_config(old_payload))
        new_scientific_hash = payload_hash(scientific_config(payload))
        if old_scientific_hash != new_scientific_hash:
            raise OrchestratorError(
                "--resume configuration does not match the existing run. "
                "Use the original inputs/references/options or choose a new --run-name."
            )
        old_runtime_hash = payload_hash(runtime_config(old_payload))
        new_runtime_hash = payload_hash(runtime_config(payload))
        if old_runtime_hash != new_runtime_hash:
            print(
                "WARNING: runtime software or scheduling controls differ from the "
                "existing run; resume will continue.",
                file=sys.stderr,
            )
            print(
                f"  previous release: {old_payload.get('release')}",
                file=sys.stderr,
            )
            print(f"  current release:  {payload.get('release')}", file=sys.stderr)
            if old_payload.get("scheduling") != payload.get("scheduling"):
                print(
                    f"  previous scheduling: {old_payload.get('scheduling')}",
                    file=sys.stderr,
                )
                print(
                    f"  current scheduling:  {payload.get('scheduling')}",
                    file=sys.stderr,
                )
            history = list(old.get("runtime_update_history") or [])
            history.append(
                {
                    "accepted_utc": utc_now(),
                    "requested_resubmits": list(requested_resubmits or []),
                    "previous": runtime_config(old_payload),
                    "replacement": runtime_config(payload),
                    "previous_runtime_hash": old_runtime_hash,
                    "replacement_runtime_hash": new_runtime_hash,
                    "policy": "warn_and_continue",
                }
            )
            stored = dict(payload)
            stored["configuration_hash"] = new_scientific_hash
            stored["runtime_hash"] = new_runtime_hash
            stored["runtime_update_history"] = history
            stored["created_utc"] = old.get("created_utc", utc_now())
            stored["production_promotion"] = old.get(
                "production_promotion", "disabled"
            )
            atomic_json(config_path, stored)
    else:
        run_dir.mkdir(parents=True, exist_ok=False)
        stored = dict(payload)
        stored["configuration_hash"] = payload_hash(scientific_config(payload))
        stored["runtime_hash"] = payload_hash(runtime_config(payload))
        stored["runtime_update_history"] = []
        stored["created_utc"] = utc_now()
        stored["production_promotion"] = "disabled"
        atomic_json(config_path, stored)

    for relative in (
        "control/generation_logs",
        "control/slurm",
        "logs",
        "rna3",
        "rna5",
        "atac",
        "qc",
        "validation",
    ):
        (run_dir / relative).mkdir(parents=True, exist_ok=True)


def run_generator(command: Sequence[str], log_path: Path) -> None:
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    combined = (
        "$ " + " ".join(q(item) for item in command) + "\n\n"
        + result.stdout
        + ("\nSTDERR:\n" + result.stderr if result.stderr else "")
    )
    atomic_text(log_path, combined)
    if result.returncode != 0:
        raise OrchestratorError(
            f"child script generation failed (exit {result.returncode}); see {log_path}"
        )


def required_replace(text: str, old: str, new: str, context: str) -> str:
    if old not in text:
        raise OrchestratorError(
            f"cannot safely normalize generated {context}: expected text was not found"
        )
    return text.replace(old, new)


def normalize_rna_trim_script(path: Path) -> None:
    text = path.read_text(encoding="utf-8")
    for marker in (
        "module load miniforge/3 genomics-base/latest",
        "#SBATCH --chdir=",
        "SKIP: completed trimming outputs already exist",
        "--pair-filter=any",
    ):
        if marker not in text:
            raise OrchestratorError(
                f"generated RNA trimming script lacks required safety marker {marker!r}: {path}"
            )
    atomic_text(path, text, mode=0o755)


def add_success_guard(text: str, marker: Path, context: str) -> str:
    marker_q = q(marker)
    start = "set -euo pipefail\n"
    guarded = start + f"""

SUCCESS_MARKER={marker_q}
if [[ -s "$SUCCESS_MARKER" ]]; then
    echo "SKIP: {context} already completed for this staged run"
    exit 0
fi
"""
    text = required_replace(text, start, guarded, f"{context} success guard")
    success_line = "echo \"Results: "
    locations = [match.start() for match in re.finditer(re.escape(success_line), text)]
    if not locations:
        raise OrchestratorError(
            f"cannot safely normalize generated {context}: success branch was not found"
        )
    line_end = text.find("\n", locations[-1])
    text = text[: line_end + 1] + '    date -Is > "$SUCCESS_MARKER"\n' + text[line_end + 1 :]
    return text


def normalize_rna_mapping_script(path: Path, marker: Path) -> None:
    text = path.read_text(encoding="utf-8")
    for required in (
        "align_pipelines/bjp",
        "htslib/1.20 samtools/1.20 star/2.7.11b",
        "rg_metadata.tsv",
        "rna_geometry",
        "--lib-prefix",
        "#SBATCH --chdir=",
    ):
        if required not in text:
            raise OrchestratorError(
                f"generated RNA mapping script lacks {required!r}: {path}"
            )
    text = add_success_guard(text, marker, "RNA mapping")
    atomic_text(path, text, mode=0o755)


def normalize_atac_mapping_script(path: Path, marker: Path) -> None:
    text = path.read_text(encoding="utf-8")
    text = text.replace("#SBATCH --time=576:00:00", "#SBATCH --time=7-00:00:00")
    text = text.replace("time = '576 h'", "time = '168 h'")
    for required in (
        "align_pipelines/bjp",
        "rg_metadata.tsv",
        "htslib/1.20 samtools/1.20 minimap2/2.28",
        "--lib-prefix",
        "#SBATCH --chdir=",
    ):
        if required not in text:
            raise OrchestratorError(
                f"generated ATAC mapping script lacks {required!r}: {path}"
            )
    text = add_success_guard(text, marker, "ATAC mapping")
    atomic_text(path, text, mode=0o755)


def sbatch_text(
    job_name: str,
    log_root: Path,
    body: str,
    *,
    cpus: int,
    memory: str,
    walltime: str = "7-00:00:00",
    array: str | None = None,
) -> str:
    directives = [
        "#!/bin/bash",
        f"#SBATCH --job-name={job_name}",
        f"#SBATCH --output={log_root}/{job_name}_%A_%a.out" if array else
        f"#SBATCH --output={log_root}/{job_name}_%j.out",
        f"#SBATCH --error={log_root}/{job_name}_%A_%a.err" if array else
        f"#SBATCH --error={log_root}/{job_name}_%j.err",
        "#SBATCH --partition=compute",
        "#SBATCH --ntasks=1",
        f"#SBATCH --cpus-per-task={cpus}",
        f"#SBATCH --mem={memory}",
        f"#SBATCH --time={walltime}",
        f"#SBATCH --chdir={log_root.parent}",
    ]
    if array:
        directives.append(f"#SBATCH --array={array}")
    return "\n".join(directives) + "\n\nset -euo pipefail\n\n" + body.rstrip() + "\n"


def write_sbatch(path: Path, content: str) -> Path:
    atomic_text(path, content, mode=0o755)
    return path


def apply_slurm_controls(
    path: Path,
    nodelist: str | None,
    array_max_concurrent: int | None,
) -> None:
    """Apply user-requested placement and array limits to one generated job.

    Mapping scripts contain a Nextflow configuration inside their shell body.
    The outer ``#SBATCH --nodelist`` controls only the Nextflow controller, so
    the same placement expression is also appended to ``clusterOptions`` for
    every Nextflow worker submitted by that controller.
    """
    text = path.read_text(encoding="utf-8")
    original = text

    if nodelist is not None:
        nodelist_directive = f"#SBATCH --nodelist={nodelist}"
        existing_nodelist = re.search(
            r"^#SBATCH\s+(?:--nodelist(?:=|\s+)|-w\s+)[^\s]+\s*$",
            text,
            flags=re.MULTILINE,
        )
        if existing_nodelist is not None:
            text = (
                text[: existing_nodelist.start()]
                + nodelist_directive
                + text[existing_nodelist.end() :]
            )
        else:
            partition = re.search(
                r"^#SBATCH\s+--partition=[^\s]+\s*$", text, flags=re.MULTILINE
            )
            if partition is None:
                raise OrchestratorError(
                    f"cannot apply --nodelist; generated job lacks a partition "
                    f"directive: {path}"
                )
            text = (
                text[: partition.end()]
                + "\n"
                + nodelist_directive
                + text[partition.end() :]
            )

        def add_nextflow_nodelist(match: re.Match[str]) -> str:
            options = match.group(2)
            if re.search(r"(?:^|\s)--nodelist(?:=|\s+)\S+", options):
                options = re.sub(
                    r"(?:^|\s)--nodelist(?:=|\s+)\S+",
                    f" --nodelist={nodelist}",
                    options,
                    count=1,
                ).strip()
            else:
                options = f"{options} --nodelist={nodelist}".strip()
            return f"{match.group(1)}{options}{match.group(3)}"

        text = re.sub(
            r"(clusterOptions\s*=\s*')([^']*)(')",
            add_nextflow_nodelist,
            text,
        )

    if array_max_concurrent is not None:
        array_match = re.search(
            r"^#SBATCH\s+--array=([^\s]+)\s*$", text, flags=re.MULTILINE
        )
        if array_match is not None:
            array_value = array_match.group(1)
            array_range, separator, old_limit = array_value.partition("%")
            limit = array_max_concurrent
            if separator:
                try:
                    limit = min(limit, int(old_limit))
                except ValueError as exc:
                    raise OrchestratorError(
                        f"generated array has a non-numeric concurrency limit: {path}"
                    ) from exc
            replacement = f"#SBATCH --array={array_range}%{limit}"
            text = text[: array_match.start()] + replacement + text[array_match.end() :]

    if text != original:
        atomic_text(path, text, mode=0o755)


def enforce_core_ceiling(jobs: Sequence[JobSpec], max_cores: int) -> None:
    """Serialize outer jobs and throttle every SLURM array to a hard ceiling.

    RNA mapping is a two-core controller that submits Nextflow workers. The RNA
    driver separately limits that internal queue; serializing these outer jobs
    prevents plotting/QC allocations from overlapping it or each other.
    """
    seen: set[str] = set()
    previous: str | None = None
    for job in jobs:
        unavailable = sorted(set(job.dependencies) - seen)
        if unavailable:
            raise OrchestratorError(
                f"cannot core-throttle non-topological job {job.label}; "
                f"unavailable dependencies: {', '.join(unavailable)}"
            )
        text = job.script.read_text(encoding="utf-8")
        cpus_match = re.search(
            r"^#SBATCH\s+--cpus-per-task=(\d+)\s*$", text, flags=re.MULTILINE
        )
        if cpus_match is None:
            raise OrchestratorError(
                f"cannot enforce --max-cores; job lacks cpus-per-task: {job.script}"
            )
        cpus = int(cpus_match.group(1))
        if cpus > max_cores:
            raise OrchestratorError(
                f"job requests {cpus} cores above --max-cores {max_cores}: {job.script}"
            )
        array_match = re.search(
            r"^#SBATCH\s+--array=([^\s]+)\s*$", text, flags=re.MULTILINE
        )
        if array_match is not None:
            array_value = array_match.group(1)
            array_range, separator, requested_limit = array_value.partition("%")
            max_active = max(1, max_cores // cpus)
            if separator:
                try:
                    max_active = min(max_active, int(requested_limit))
                except ValueError as exc:
                    raise OrchestratorError(
                        f"generated array has a non-numeric concurrency limit: "
                        f"{job.script}"
                    ) from exc
            replacement = f"#SBATCH --array={array_range}%{max_active}"
            text = text[: array_match.start()] + replacement + text[array_match.end() :]
            atomic_text(job.script, text, mode=0o755)
        if previous is not None and previous not in job.dependencies:
            job.dependencies.append(previous)
        seen.add(job.label)
        previous = job.label


def generate_rna_jobs(
    modality: str,
    chemistry: str,
    read_format: str,
    runs: Sequence[Path],
    whitelist: str,
    root: Path,
    run_dir: Path,
    stages: set[str],
    args: argparse.Namespace,
    resources: Resources,
) -> list[JobSpec]:
    if not ({"trim", "map"} & stages):
        return []
    assert resources.rna_driver is not None
    assert resources.rna_workflow is not None
    command = [
        sys.executable,
        str(resources.rna_driver),
        "--input-dirs",
        *[str(path) for path in runs],
        "--chemistry",
        chemistry,
        "--read-format",
        read_format,
        "--output-base",
        str(root),
        "--rna-ref",
        args.rna_ref,
        "--rna-whitelist",
        whitelist,
        "--workflow-file",
        str(resources.rna_workflow),
        "--memgb",
        str(args.rna_mem_gb),
        "--threads",
        str(args.rna_threads),
        "--min-pe150-tso-match-fraction",
        str(args.min_pe150_tso_match_fraction),
        "--lib-prefix",
        args.rna3_lib_prefix if chemistry == "3prime" else args.rna5_lib_prefix,
    ]
    if args.max_cores is not None:
        command.extend(["--max-cores", str(args.max_cores)])
    if args.libraries:
        command.extend(["--libraries", *[str(value) for value in args.libraries]])
    if "trim" not in stages:
        command.append("--skip-trimming")
    if "map" not in stages:
        command.append("--skip-mapping")
    if args.no_trim_info:
        command.append("--no-info-file")
    run_generator(
        command, run_dir / "control" / "generation_logs" / f"{modality}_driver.log"
    )

    jobs: list[JobSpec] = []
    trim_labels: list[str] = []
    if "trim" in stages:
        suffix = "_noinfo" if args.no_trim_info else ""
        for run in runs:
            script = root / f"trim_{chemistry}_{run.name}{suffix}.sh"
            if not script.is_file():
                raise OrchestratorError(f"RNA driver did not generate expected script: {script}")
            normalize_rna_trim_script(script)
            label = f"{modality}_trim_{run.name}"
            trim_labels.append(label)
            jobs.append(JobSpec(label=label, script=script))

    if "map" in stages:
        script = root / "run_mapping.sbatch"
        if not script.is_file():
            raise OrchestratorError(f"RNA driver did not generate expected script: {script}")
        normalize_rna_mapping_script(
            script, root / "mapping_project" / "MAPPING_COMPLETE.ok"
        )
        jobs.append(
            JobSpec(
                label=f"{modality}_map",
                script=script,
                dependencies=trim_labels,
            )
        )
    return jobs


def generate_atac_mapping_job(
    runs: Sequence[Path],
    root: Path,
    run_dir: Path,
    args: argparse.Namespace,
    resources: Resources,
) -> JobSpec:
    assert resources.atac_driver is not None
    assert resources.atac_workflow is not None
    command = [
        sys.executable,
        str(resources.atac_driver),
        "--input-dirs",
        *[str(path) for path in runs],
        "--output-base",
        str(root),
        "--atac-ref",
        args.atac_ref,
        "--rna-whitelist",
        args.rna3_whitelist,
        "--atac-whitelist",
        args.atac_whitelist,
        "--workflow-file",
        str(resources.atac_workflow),
        "--memgb",
        str(args.atac_mem_gb),
        "--threads",
        str(args.atac_threads),
        "--num-chunks",
        str(args.atac_chunks),
        "--lib-prefix",
        args.atac_lib_prefix,
    ]
    if args.libraries:
        command.extend(["--libraries", *[str(value) for value in args.libraries]])
    run_generator(
        command, run_dir / "control" / "generation_logs" / "atac_driver.log"
    )
    script = root / "run_atac_mapping.sbatch"
    if not script.is_file():
        raise OrchestratorError(f"ATAC driver did not generate expected script: {script}")
    normalize_atac_mapping_script(
        script, root / "mapping_project" / "MAPPING_COMPLETE.ok"
    )
    return JobSpec(label="atac_map", script=script)


def trim_sample_key(r2: Path) -> str:
    name = r2.name.replace("_R2_", "_")
    return re.sub(r"\.(?:fastq|fq)\.gz$", "", name)


def expected_rna_rg_id(library: str, run_name: str, fastq: Path) -> str:
    match = re.search(r"_(S\d+)_(L\d+)", fastq.name)
    if match is None:
        raise OrchestratorError(
            f"FASTQ name lacks the S#/L### fields required for an RG ID: {fastq}"
        )
    return f"{library}_{run_name}_{match.group(1)}_{match.group(2)}"


def load_rg_metadata_by_raw_r1(root: Path) -> dict[Path, dict[str, str]]:
    """Index the mapping driver's authoritative source-FASTQ/RG metadata."""
    path = root / "mapping_project" / "rg_metadata.tsv"
    if not path.is_file():
        raise OrchestratorError(
            f"RNA driver did not generate source read-group metadata: {path}"
        )
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"raw_r1", "raw_r2", "rg_id", "bp_id"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise OrchestratorError(
                f"RG metadata is missing column(s) {', '.join(sorted(missing))}: {path}"
            )
        records: dict[Path, dict[str, str]] = {}
        for row in reader:
            raw_r1 = Path(row["raw_r1"]).resolve(strict=False)
            if raw_r1 in records:
                raise OrchestratorError(f"duplicate raw R1 in RG metadata: {raw_r1}")
            records[raw_r1] = row
    if not records:
        raise OrchestratorError(f"RG metadata contains no FASTQ records: {path}")
    return records


def generate_barcode_trim_jobs(
    modality: str,
    groups: Sequence[tuple[str, Sequence[Path], Path]],
    selected_libraries: Sequence[str],
    lib_prefix: str,
    whitelist: str,
    empty_drop_roster: str | None,
    run_dir: Path,
    resources: Resources,
    dependencies: Sequence[str],
    use_evidence_bridge: bool = False,
) -> list[JobSpec]:
    """Generate full-population trim joins plus labelled derived views."""
    assert resources.barcode_aggregator is not None
    control = run_dir / "control"
    qc_root = run_dir / "qc" / modality / "trim_by_barcode"
    per_library = qc_root / "per_library"
    manifest = control / f"{modality}_trim_barcode_manifest.tsv"
    rows: list[str] = []
    libraries: set[str] = set()
    selected = set(selected_libraries)
    for read_format, runs, root in groups:
        rg_metadata = root / "mapping_project" / "rg_metadata.tsv"
        rg_by_r1 = load_rg_metadata_by_raw_r1(root) if rg_metadata.is_file() else None
        for run in runs:
            log_dir = root / "trimming" / run.name / "logs"
            for r1, r2 in fastq_pairs(run):
                library = library_from_fastq(r1, lib_prefix)
                if library not in selected:
                    continue
                libraries.add(library)
                sample = trim_sample_key(r2)
                source_id = expected_rna_rg_id(library, run.name, r1)
                source_fastq = str(r2.resolve(strict=False))
                if rg_by_r1 is not None:
                    rg_record = rg_by_r1.get(r1.resolve(strict=False))
                    if rg_record is None:
                        raise OrchestratorError(
                            f"raw FASTQ has no matching RG metadata record: {r1}"
                        )
                    if rg_record["bp_id"] != run.name:
                        raise OrchestratorError(
                            f"RG metadata BP mismatch for {r1}: "
                            f"{rg_record['bp_id']} != {run.name}"
                        )
                    if rg_record["rg_id"] != source_id:
                        raise OrchestratorError(
                            f"RG metadata ID mismatch for {r1}: "
                            f"{rg_record['rg_id']} != {source_id}"
                        )
                    source_fastq = rg_record["raw_r2"]
                raw = root / "mapping_output" / library / "raw" / "barcodes.tsv.gz"
                filtered = root / "mapping_output" / library / "filtered" / "barcodes.tsv.gz"
                bam = root / "mapping_output" / library / "gex.bam"
                bridge = (
                    root
                    / "mapping_output"
                    / library
                    / "bam_evidence"
                    / "raw_to_corrected_barcode_counts.tsv.gz"
                )
                rows.append(
                    "\t".join(
                        [
                            library,
                            run.name,
                            source_id,
                            read_format,
                            "R2",
                            str(r1),
                            source_fastq,
                            str(log_dir / f"{sample}_R2_adapter_info.txt.gz"),
                            "",
                            str(raw),
                            str(filtered),
                            "" if use_evidence_bridge else str(bam),
                            str(rg_metadata),
                            str(bridge) if use_evidence_bridge else "",
                        ]
                    )
                    + "\n"
                )
                if read_format == "pe150":
                    rows.append(
                        "\t".join(
                            [
                                library,
                                run.name,
                                source_id,
                                read_format,
                                "R1_cDNA",
                                str(r1),
                                str(r1.resolve(strict=False)),
                                str(log_dir / f"{sample}_R1_adapter_info.txt.gz"),
                                str(log_dir / f"{sample}_R1_tso_info.tsv.gz"),
                                str(raw),
                                str(filtered),
                                "" if use_evidence_bridge else str(bam),
                                str(rg_metadata),
                                str(bridge) if use_evidence_bridge else "",
                            ]
                        )
                        + "\n"
                    )
    if not rows or not libraries:
        raise OrchestratorError(f"no RNA FASTQ pairs available for {modality} barcode QC")
    atomic_text(
        manifest,
        "library\trun_id\tsource_id\tread_format\tmate\tbarcode_fastq\t"
        "source_fastq\tinfo_file\ttso_info_file\traw_barcodes\tfiltered_barcodes\t"
        "bam\trg_metadata\tbarcode_correction_bridge\n"
        + "".join(rows),
    )
    library_file = control / f"{modality}_trim_barcode_libraries.txt"
    atomic_text(library_file, "".join(f"{library}\n" for library in sorted(libraries)))
    array_script = run_dir / "control" / "slurm" / f"{modality}_trim_by_barcode.sbatch"
    empty_drop_arg = (
        f" \\\n  --empty-drop-roster {q(empty_drop_roster)}"
        if empty_drop_roster
        else ""
    )
    array_body = f"""module purge
module load miniforge/3 genomics-base/latest
module load htslib/1.20 samtools/1.20
command -v samtools >/dev/null || {{ echo "ERROR: samtools/1.20 module did not provide samtools" >&2; exit 1; }}
mkdir -p {q(per_library)}
LIBRARY=$(sed -n "${{SLURM_ARRAY_TASK_ID}}p" {q(library_file)})
[[ -n "$LIBRARY" ]] || {{ echo "ERROR: empty library for array task" >&2; exit 1; }}
python3 {q(resources.barcode_aggregator)} process \\
  --manifest {q(manifest)} \\
  --library "$LIBRARY" \\
  --whitelist {q(whitelist)} \\
  --output-dir {q(per_library)}{empty_drop_arg}
"""
    write_sbatch(
        array_script,
        sbatch_text(
            f"{modality}_trim_bc",
            run_dir / "logs",
            array_body,
            cpus=2,
            memory="24G",
            walltime="3-00:00:00",
            array=f"1-{len(libraries)}",
        ),
    )
    array_label = f"{modality}_trim_by_barcode"
    merge_script = run_dir / "control" / "slurm" / f"{modality}_trim_by_barcode_merge.sbatch"
    marker = qc_root / "TRIM_BY_BARCODE_COMPLETE.ok"
    merge_body = f"""module purge
module load miniforge/3
if [[ -s {q(marker)} ]]; then
  echo "SKIP: {modality} trim-by-barcode merge already completed"
  exit 0
fi
mkdir -p {q(qc_root)}
python3 {q(resources.barcode_aggregator)} merge \\
  --per-library-dir {q(per_library)} \\
  --output-dir {q(qc_root)}
date -Is > {q(marker)}
"""
    write_sbatch(
        merge_script,
        sbatch_text(
            f"{modality}_trim_bc_merge",
            run_dir / "logs",
            merge_body,
            cpus=2,
            memory="8G",
            walltime="1-00:00:00",
        ),
    )
    return [
        JobSpec(array_label, array_script, dependencies=list(dependencies)),
        JobSpec(
            f"{modality}_trim_by_barcode_merge",
            merge_script,
            dependencies=[array_label],
        ),
    ]


def generate_rna5_together_view_job(
    groups: Sequence[tuple[str, Sequence[Path], Path]],
    discovery: Discovery,
    run_dir: Path,
    dependencies: Sequence[str],
) -> JobSpec:
    """Create a common output view when formats contain disjoint libraries."""
    view_root = run_dir / "rna5" / "mapping_output"
    manifest = run_dir / "control" / "rna5_together_view.tsv"
    rows: list[str] = []
    for read_format, _runs, root in groups:
        for library in discovery.rna5_libraries_by_format.get(read_format, []):
            rows.append(
                f"{read_format}\t{library}\t{root / 'mapping_output' / library}\t"
                f"{view_root / library}\n"
            )
    atomic_text(manifest, "format\tlibrary\tsource\ttarget\n" + "".join(rows))
    script = run_dir / "control" / "slurm" / "rna5_together_view.sbatch"
    marker = run_dir / "rna5" / "RNA5_TOGETHER_VIEW_COMPLETE.ok"
    body = f"""if [[ -s {q(marker)} ]]; then
  echo "SKIP: RNA5 together view already completed"
  exit 0
fi
mkdir -p {q(view_root)}
while IFS=$'\t' read -r FORMAT LIBRARY SOURCE TARGET; do
  [[ "$FORMAT" == "format" ]] && continue
  [[ -d "$SOURCE" ]] || {{ echo "ERROR: missing mapping output $SOURCE" >&2; exit 1; }}
  if [[ -L "$TARGET" ]]; then
    [[ "$(readlink -f "$TARGET")" == "$(readlink -f "$SOURCE")" ]] || {{
      echo "ERROR: existing view symlink has a different source: $TARGET" >&2
      exit 1
    }}
  elif [[ -e "$TARGET" ]]; then
    echo "ERROR: refusing to overwrite existing together-view target: $TARGET" >&2
    exit 1
  else
    ln -s "$SOURCE" "$TARGET"
  fi
done < {q(manifest)}
date -Is > {q(marker)}
"""
    write_sbatch(
        script,
        sbatch_text(
            "rna5_together_view",
            run_dir / "logs",
            body,
            cpus=1,
            memory="2G",
            walltime="2:00:00",
        ),
    )
    return JobSpec("rna5_together_view", script, dependencies=list(dependencies))


def generate_trim_plot_job(
    modality: str,
    root: Path,
    run_dir: Path,
    figure_root: Path,
    resources: Resources,
    dependencies: Sequence[str],
) -> JobSpec:
    assert resources.trim_plotter is not None
    script = run_dir / "control" / "slurm" / f"plot_{modality}_trimming.sbatch"
    work = figure_root / "trimming"
    report_dir = run_dir / "qc" / modality / "trimming" / "json_reports"
    marker = work / "TRIMMING_PLOTS_COMPLETE.ok"
    body = f"""if [[ -s {q(marker)} ]]; then
    echo "SKIP: {modality} trimming plots already completed"
    exit 0
fi
module purge
module load miniforge/3 genomics-base/latest
export MPLBACKEND=Agg
mkdir -p {q(work)} {q(report_dir)}
python3 {q(Path(__file__).resolve())} --_gather-trim-json \\
    --_trim-root {q(root / 'trimming')} \\
    --_report-dir {q(report_dir)}
cd {q(work)}
python3 {q(resources.trim_plotter)} --path {q(report_dir)}
date -Is > {q(marker)}
"""
    write_sbatch(
        script,
        sbatch_text(
            f"plot_{modality}_trim",
            run_dir / "logs",
            body,
            cpus=4,
            memory="24G",
            walltime="1-00:00:00",
        ),
    )
    return JobSpec(
        label=f"{modality}_trim_plot",
        script=script,
        dependencies=list(dependencies),
    )


def generate_mapping_plot_job(
    modality: str,
    root: Path,
    run_dir: Path,
    figure_root: Path,
    resources: Resources,
    dependencies: Sequence[str],
    baseline_stats: str | None = None,
    baseline_label: str = "Baseline",
    current_label: str = "Current",
    baseline_base_path: str | None = None,
    baseline_cell_reads_work_root: str | None = None,
) -> JobSpec:
    assert resources.mapping_plotter is not None
    script = run_dir / "control" / "slurm" / f"plot_{modality}_mapping.sbatch"
    out = figure_root / "mapping_QC"
    marker = out / "MAPPING_PLOTS_COMPLETE.ok"
    expected_outputs = [
        out / "processedstats.tsv",
        out / "quality_control_dashboard.png",
        out / "cell_quality_matrix.png",
        out / "outlier_detection_report.png",
        out / "reads_per_cell.tsv.gz",
        out / "reads_per_cell_distribution.png",
    ]
    comparison_args = ""
    if baseline_stats:
        expected_outputs.extend(
            [
                out / "mapping_stats_deltas.tsv",
                out / "mapping_delta_dashboard.png",
            ]
        )
        comparison_args = (
            f" \\\n    --baseline-stats {q(baseline_stats)}"
        )
    cell_read_args = (
        f" \\\n    --cell-read-distribution"
        f" \\\n    --cell-reads-work-root {q(root / 'mapping_project' / 'work')}"
        f" \\\n    --baseline-label {q(baseline_label)}"
        f" \\\n    --current-label {q(current_label)}"
    )
    if baseline_base_path:
        cell_read_args += f" \\\n    --baseline-base-path {q(baseline_base_path)}"
    if baseline_cell_reads_work_root:
        cell_read_args += (
            f" \\\n    --baseline-cell-reads-work-root {q(baseline_cell_reads_work_root)}"
        )
    output_checks = "\n".join(
        f"test -s {q(path)} || {{ echo {q('ERROR: missing mapping plot output: ' + str(path))} >&2; exit 1; }}"
        for path in expected_outputs
    )
    completion_condition = " && ".join(
        [f"[[ -s {q(marker)} ]]", *[f"[[ -s {q(path)} ]]" for path in expected_outputs]]
    )
    body = f"""if {completion_condition}; then
    echo "SKIP: {modality} mapping plots already completed"
    exit 0
fi
module purge
module load miniforge/3 genomics-base/latest
export MPLBACKEND=Agg
mkdir -p {q(out)}
python3 {q(resources.mapping_plotter)} \\
    --base-path {q(root / 'mapping_output')} \\
    --output-dir {q(out)} \\
    --plot-type all{comparison_args}{cell_read_args}
{output_checks}
date -Is > {q(marker)}
"""
    write_sbatch(
        script,
        sbatch_text(
            f"plot_{modality}_map",
            run_dir / "logs",
            body,
            cpus=4,
            memory="24G",
            walltime="1-00:00:00",
        ),
    )
    return JobSpec(
        label=f"{modality}_mapping_plot",
        script=script,
        dependencies=list(dependencies),
    )


def _existing_optional_product(directory: Path, relative: str) -> str:
    """Resolve an optional plain/gzip baseline product without guessing content."""
    candidate = directory / relative
    alternatives = [candidate]
    if candidate.suffix == ".gz":
        alternatives.append(candidate.with_suffix(""))
    else:
        alternatives.append(candidate.with_name(candidate.name + ".gz"))
    return str(next((path.resolve() for path in alternatives if path.is_file()), ""))


def _scientific_text_snapshot(path: Path, label: str) -> dict[str, object]:
    if not path.is_file() or path.stat().st_size == 0:
        raise OrchestratorError(f"{label} is missing or empty: {path}")
    raw = path.read_bytes()
    try:
        contents = raw.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise OrchestratorError(f"{label} is not UTF-8 text: {path}") from exc
    return {
        "path": str(path.resolve()),
        "sha256": hashlib.sha256(raw).hexdigest(),
        "contents": contents,
    }


def _frozen_bam_evidence_config(
    manifest_fields: Sequence[str],
    rows: Sequence[dict[str, object]],
    source_order_path: Path,
    source_order_text: str,
    source_order_provenance: str,
    source_values: Sequence[str],
    baseline_root: Path | None,
    class_manifest: str,
    hash_bins: int,
    profiler: Path,
) -> dict[str, object]:
    source_bytes = source_order_text.encode("utf-8")
    classification: dict[str, object]
    if class_manifest:
        classification = {
            "status": "manifest",
            **_scientific_text_snapshot(
                Path(class_manifest), "RNA BAM evidence class manifest"
            ),
        }
    else:
        classification = {
            "status": "explicit_unavailable",
            "path": "",
            "sha256": "",
            "contents": "",
        }
    payload: dict[str, object] = {
        "schema_version": 1,
        "immutability": "created_by_first_run_and_required_unchanged_on_resume",
        "selected_libraries": [str(row["library"]) for row in rows],
        "manifest_fields": list(manifest_fields),
        "manifest_rows": [
            {field: str(row[field]) for field in manifest_fields}
            for row in rows
        ],
        "resolved_current_inputs": {
            str(row["library"]): {
                field: str(row[field])
                for field in (
                    "library_dir", "bam", "bam_index", "summary",
                    "raw_barcodes", "filtered_barcodes", "raw_features",
                    "raw_matrix", "filtered_matrix", "rg_metadata",
                    "native_cell_reads",
                )
            }
            for row in rows
        },
        "baseline": {
            "root": str(baseline_root) if baseline_root else "",
            "resolved_inputs": {
                str(row["library"]): {
                    "old_raw_barcodes": str(row["old_raw_barcodes"]),
                    "old_filtered_barcodes": str(row["old_filtered_barcodes"]),
                }
                for row in rows
            },
        },
        "source_order": {
            "path": str(source_order_path.resolve()),
            "provenance": source_order_provenance,
            "values": list(source_values),
            "sha256": hashlib.sha256(source_bytes).hexdigest(),
            "contents": source_order_text,
        },
        "biological_classification": classification,
        "hash": {
            "algorithm": RNA_BAM_EVIDENCE_HASH_ALGORITHM,
            "seed": DEFAULT_BAM_EVIDENCE_HASH_SEED,
            "bins": hash_bins,
        },
        "profiler_executable": {
            "path": str(profiler.resolve()),
            "sha256": file_sha256(profiler),
            "bytes": profiler.stat().st_size,
        },
        "starsolo": {
            "feature": RNA_BAM_EVIDENCE_STARSOLO_FEATURE,
            "umi_filtering": RNA_BAM_EVIDENCE_STARSOLO_UMI_FILTERING,
            "umi_dedup": RNA_BAM_EVIDENCE_STARSOLO_UMI_DEDUP,
            "multimappers": RNA_BAM_EVIDENCE_STARSOLO_MULTIMAPPERS,
            "ordinary_countedU_read_definition": (
                RNA_BAM_EVIDENCE_ORDINARY_COUNTEDU_READ_DEFINITION
            ),
            "ordinary_molecule_definition": (
                RNA_BAM_EVIDENCE_ORDINARY_MOLECULE_DEFINITION
            ),
            "multimapper_definition": RNA_BAM_EVIDENCE_MULTIMAPPER_DEFINITION,
            "summary_unique_read_metric": (
                RNA_BAM_EVIDENCE_SUMMARY_UNIQUE_READ_METRIC
            ),
            "nh_gt1_unique_gene_countedU_definition": (
                RNA_BAM_EVIDENCE_NH_GT1_UNIQUE_GENE_DEFINITION
            ),
            "starsolo_EM_evidence_availability": (
                RNA_BAM_EVIDENCE_EM_AVAILABILITY
            ),
        },
    }
    payload["configuration_hash"] = payload_hash(payload)
    return payload


def _validate_or_create_bam_evidence_config(
    path: Path,
    expected: dict[str, object],
) -> None:
    if path.is_file():
        try:
            stored = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise OrchestratorError(
                f"could not read frozen RNA BAM evidence configuration: {path}"
            ) from exc
        if stored != expected:
            raise OrchestratorError(
                "RNA BAM evidence scientific configuration differs from the "
                "frozen configuration; no BAM-reader job was generated"
            )
        return
    atomic_json(path, expected)


def _validate_or_create_frozen_text(
    path: Path,
    expected: str,
    label: str,
) -> None:
    if path.is_file():
        try:
            actual = path.read_text(encoding="utf-8")
        except OSError as exc:
            raise OrchestratorError(f"could not read frozen {label}: {path}") from exc
        if actual != expected:
            raise OrchestratorError(
                f"{label} differs from its frozen copy; no BAM-reader job was generated"
            )
        return
    atomic_text(path, expected)


def _reset_failed_bam_evidence_run(
    run_dir: Path,
    rows: Sequence[dict[str, object]],
    scientific_config_path: Path,
    frozen_config: dict[str, object],
    source_order_path: Path,
    source_order_text: str,
    manifest: Path,
    manifest_text: str,
) -> None:
    """Archive an obsolete failed evidence run after proving its jobs terminal."""
    if not scientific_config_path.is_file():
        raise OrchestratorError(
            "--rna3-bam-evidence-reset-failed-run requires an existing frozen "
            "BAM-evidence configuration"
        )

    plan_path = run_dir / "control" / "job_plan.json"
    if not plan_path.is_file() or plan_path.stat().st_size == 0:
        raise OrchestratorError(
            "refusing failed-run reset because the prior submission ledger "
            f"is missing or empty: {plan_path}"
        )
    try:
        plan = json.loads(plan_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise OrchestratorError(
            f"cannot verify prior BAM-reader state from {plan_path}: {exc}"
        ) from exc
    if (
        not isinstance(plan, dict)
        or not isinstance(plan.get("jobs"), list)
        or not plan["jobs"]
        or any(not isinstance(job, dict) for job in plan["jobs"])
    ):
        raise OrchestratorError(
            f"cannot verify prior BAM-reader state from malformed or empty {plan_path}"
        )
    for job in plan["jobs"]:
        label = str(job.get("label", ""))
        raw_job_id = job.get("job_id")
        protected = (
            label.startswith("rna3_bam_evidence")
            or label == "rna3_star_diagnostics"
        )
        if not raw_job_id or not protected:
            continue
        job_id = str(raw_job_id)
        states = slurm_accounting_states(job_id)
        if not states:
            raise OrchestratorError(
                "refusing failed-run reset because SLURM accounting "
                f"could not prove RNA evidence job {job_id} ({label}) is "
                "terminal"
            )
        active = sorted(set(states).intersection(SLURM_ACTIVE_STATES))
        if active:
            raise OrchestratorError(
                "refusing failed-run reset while prior RNA evidence job "
                f"{job_id} ({label}) is active: {','.join(active)}"
            )
        unknown = sorted(
            set(states) - SLURM_FAILURE_STATES - {"COMPLETED"}
        )
        if unknown:
            raise OrchestratorError(
                "refusing failed-run reset because SLURM did not prove "
                f"RNA evidence job {job_id} ({label}) terminal: "
                f"{','.join(unknown)}"
            )

    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    old_hash = file_sha256(scientific_config_path)[:16]
    archive = (
        run_dir / "control" / "rna3_bam_evidence_failed_run_migrations"
        / f"{timestamp}_{old_hash}"
    )
    archive.mkdir(parents=True, exist_ok=False)
    for path in (
        scientific_config_path,
        source_order_path,
        manifest,
        run_dir / "control" / "rna3_bam_evidence_pilot_selection.json",
        run_dir / "control" / "rna3_bam_evidence_phase.json",
        run_dir / "control" / "rna3_bam_evidence_execution.json",
        plan_path,
    ):
        if path.is_file():
            shutil.copy2(path, archive / path.name)
    archived_scripts = archive / "slurm"
    old_scripts = list((run_dir / "control" / "slurm").glob(
        "rna3_bam_evidence*.sbatch"
    ))
    diagnostic_script = (
        run_dir / "control" / "slurm" / "rna3_star_diagnostic_promotion.sbatch"
    )
    if diagnostic_script.is_file():
        old_scripts.append(diagnostic_script)
    for path in sorted(set(old_scripts)):
        archived_scripts.mkdir(parents=True, exist_ok=True)
        shutil.move(str(path), archived_scripts / path.name)

    archived_evidence = archive / "evidence"
    for row in rows:
        output_dir = Path(str(row["output_dir"]))
        if output_dir.exists():
            destination = archived_evidence / str(row["library"])
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(output_dir), destination)
    gathered = run_dir / "qc" / "rna3" / "bam_evidence"
    if gathered.exists():
        shutil.move(str(gathered), archive / "gathered_bam_evidence")
    atomic_json(
        archive / "migration.json",
        {
            "created_utc": utc_now(),
            "reason": "explicit_obsolete_failed_run_semantic_reset",
            "prior_jobs_proven_terminal": True,
            "prior_evidence_moved_to_archive": True,
            "replacement_configuration_hash": frozen_config["configuration_hash"],
        },
    )
    if plan_path.is_file():
        # Preserve unrelated mapping job IDs, but make every BAM-evidence job
        # new for restore_prior_job_ids(). Reusing the failed array ID
        # would otherwise abort or, with expired accounting, falsely satisfy a
        # dependency after this explicit scientific-contract migration.
        invalidated = {
            str(job.get("label", ""))
            for job in plan.get("jobs", [])
            if (
                str(job.get("label", "")).startswith("rna3_bam_evidence")
                or str(job.get("label", "")) == "rna3_star_diagnostics"
            )
        }
        changed = True
        while changed:
            changed = False
            for job in plan.get("jobs", []):
                label = str(job.get("label", ""))
                dependencies = {
                    str(value) for value in job.get("dependencies", [])
                }
                if label not in invalidated and dependencies.intersection(invalidated):
                    invalidated.add(label)
                    changed = True
        retained_jobs = [
            job for job in plan.get("jobs", [])
            if str(job.get("label", "")) not in invalidated
        ]
        replacement_plan = dict(plan)
        replacement_plan["jobs"] = retained_jobs
        replacement_plan["bam_evidence_failed_run_reset_utc"] = utc_now()
        atomic_json(plan_path, replacement_plan)
    for obsolete in (
        run_dir / "control" / "rna3_bam_evidence_phase.json",
        run_dir / "control" / "rna3_bam_evidence_pilot_selection.json",
        run_dir / "control" / "rna3_bam_evidence_execution.json",
    ):
        obsolete.unlink(missing_ok=True)
    atomic_json(scientific_config_path, frozen_config)
    atomic_text(source_order_path, source_order_text)
    atomic_text(manifest, manifest_text)


def generate_rna_bam_evidence_jobs(
    libraries: Sequence[str],
    root: Path,
    run_dir: Path,
    args: argparse.Namespace,
    resources: Resources,
    dependencies: Sequence[str],
    inferred_source_order: Sequence[str],
) -> list[JobSpec]:
    """Create one bounded all-library evidence array and BAM-free gather."""
    if not libraries:
        raise OrchestratorError("no 3' RNA libraries were selected for BAM evidence")
    if len(libraries) >= 1000:
        raise OrchestratorError(
            "RNA BAM evidence requires fewer than 1000 selected library tasks"
        )
    assert resources.bam_evidence_runner is not None
    assert resources.bam_evidence_profiler is not None
    assert resources.star_diagnostic_promoter is not None
    mapping_output = root / "mapping_output"
    rg_metadata = root / "mapping_project" / "rg_metadata.tsv"
    source_values = list(
        args.rna3_bam_evidence_source_order or inferred_source_order
    )
    if not source_values or len(source_values) != len(set(source_values)):
        raise OrchestratorError(
            "RNA BAM evidence source order must be nonempty and contain no duplicates"
        )
    source_order_provenance = (
        "explicit --rna3-bam-evidence-source-order"
        if args.rna3_bam_evidence_source_order
        else "original user-provided --rna3-runs order"
    )
    control = run_dir / "control"
    source_order_path = control / "rna3_bam_evidence_source_order.tsv"
    source_order_text = "source_id\n" + "".join(
        f"{value}\n" for value in source_values
    )
    scientific_config_path = control / "rna3_bam_evidence_scientific_config.json"
    baseline_root = (
        Path(args.rna3_bam_evidence_baseline_root).resolve(strict=False)
        if args.rna3_bam_evidence_baseline_root
        else None
    )
    class_manifest = (
        str(Path(args.rna3_bam_evidence_class_manifest).resolve(strict=False))
        if args.rna3_bam_evidence_class_manifest
        else ""
    )
    rows: list[dict[str, object]] = []
    for library in libraries:
        library_dir = (mapping_output / library).resolve(strict=False)
        bam = library_dir / "gex.bam"
        baseline_library = baseline_root / library if baseline_root else None
        native_cell_reads = _existing_optional_product(
            library_dir, "CellReads.stats.gz"
        )
        rows.append(
            {
                "library": library,
                "library_dir": str(library_dir),
                "bam": str(bam),
                "bam_index": str(library_dir / "gex.bam.bai"),
                "summary": str(library_dir / "Summary.csv"),
                "raw_barcodes": str(library_dir / "raw" / "barcodes.tsv.gz"),
                "filtered_barcodes": str(
                    library_dir / "filtered" / "barcodes.tsv.gz"
                ),
                "raw_features": str(library_dir / "raw" / "features.tsv.gz"),
                "raw_matrix": str(library_dir / "raw" / "matrix.mtx.gz"),
                "filtered_matrix": str(
                    library_dir / "filtered" / "matrix.mtx.gz"
                ),
                "rg_metadata": str(rg_metadata.resolve(strict=False)),
                "source_order": str(source_order_path.resolve()),
                "scientific_config": str(scientific_config_path.resolve()),
                "source_order_provenance": source_order_provenance,
                "starsolo_feature": RNA_BAM_EVIDENCE_STARSOLO_FEATURE,
                "starsolo_umi_filtering": (
                    RNA_BAM_EVIDENCE_STARSOLO_UMI_FILTERING
                ),
                "starsolo_umi_dedup": RNA_BAM_EVIDENCE_STARSOLO_UMI_DEDUP,
                "starsolo_multimappers": (
                    RNA_BAM_EVIDENCE_STARSOLO_MULTIMAPPERS
                ),
                "biological_classification_intent": (
                    "manifest" if class_manifest else "explicit_unavailable"
                ),
                "old_raw_barcodes": (
                    _existing_optional_product(
                        baseline_library, "raw/barcodes.tsv.gz"
                    )
                    if baseline_library else ""
                ),
                "old_filtered_barcodes": (
                    _existing_optional_product(
                        baseline_library, "filtered/barcodes.tsv.gz"
                    )
                    if baseline_library else ""
                ),
                "class_manifest": class_manifest,
                "native_cell_reads": native_cell_reads,
                "output_dir": str(library_dir / "bam_evidence"),
                "bam_bytes": bam.stat().st_size if bam.is_file() else 0,
            }
        )

    def library_number_from_name(name: object) -> tuple[int, str]:
        match = re.search(r"(\d+)$", str(name))
        return (int(match.group(1)) if match else sys.maxsize, str(name))

    if scientific_config_path.is_file():
        try:
            frozen_order = json.loads(
                scientific_config_path.read_text(encoding="utf-8")
            ).get("selected_libraries", [])
        except (OSError, json.JSONDecodeError) as exc:
            raise OrchestratorError(
                f"could not read frozen RNA BAM evidence configuration: "
                f"{scientific_config_path}"
            ) from exc
        by_library = {str(row["library"]): row for row in rows}
        if (
            not isinstance(frozen_order, list)
            or len(frozen_order) != len(rows)
            or set(frozen_order) != set(by_library)
        ):
            raise OrchestratorError(
                "selected RNA BAM evidence libraries differ from the frozen "
                "configuration; no BAM-reader job was generated"
            )
        rows = [by_library[str(library)] for library in frozen_order]
    else:
        rows.sort(key=lambda row: library_number_from_name(row["library"]))

    manifest = control / "rna3_bam_evidence_manifest.tsv"
    manifest_fields = [
        "library", "library_dir", "bam", "bam_index", "summary",
        "raw_barcodes", "filtered_barcodes", "raw_features", "raw_matrix",
        "filtered_matrix", "rg_metadata", "source_order", "scientific_config",
        "source_order_provenance",
        "starsolo_feature", "starsolo_umi_filtering", "starsolo_umi_dedup",
        "starsolo_multimappers", "biological_classification_intent",
        "old_raw_barcodes", "old_filtered_barcodes", "class_manifest",
        "native_cell_reads", "output_dir",
    ]
    manifest_text = (
        "\t".join(manifest_fields) + "\n"
        + "".join(
            "\t".join(str(row[field]) for field in manifest_fields) + "\n"
            for row in rows
        )
    )
    frozen_config = _frozen_bam_evidence_config(
        manifest_fields,
        rows,
        source_order_path,
        source_order_text,
        source_order_provenance,
        source_values,
        baseline_root,
        class_manifest,
        args.rna3_bam_evidence_hash_bins,
        resources.bam_evidence_profiler,
    )
    # Freeze every value-affecting input before generating the one full array.
    if getattr(args, "rna3_bam_evidence_reset_failed_run", False):
        _reset_failed_bam_evidence_run(
            run_dir,
            rows,
            scientific_config_path,
            frozen_config,
            source_order_path,
            source_order_text,
            manifest,
            manifest_text,
        )
    else:
        _validate_or_create_bam_evidence_config(
            scientific_config_path, frozen_config
        )
        _validate_or_create_frozen_text(
            source_order_path,
            source_order_text,
            "RNA BAM evidence source order",
        )
        _validate_or_create_frozen_text(
            manifest, manifest_text, "RNA BAM evidence manifest"
        )
    reader_preflight = f"""module purge
module load miniforge/3 htslib/1.20 samtools/1.20
module list
command -v python3
command -v samtools
test -x {q(resources.bam_evidence_profiler)}
{q(resources.bam_evidence_profiler)} --version
"""
    # This is a total-process RLIMIT_AS cap in the C++ process, not a molecule
    # table allowance. Ten percent remains for the Python wrapper/cgroup margin.
    memory_limit = max(1.0, args.rna3_bam_evidence_memory_gb * 0.90)
    worker = (
        f"python3 {q(resources.bam_evidence_runner)} run "
        f"--manifest {q(manifest)} "
        f"--profiler {q(resources.bam_evidence_profiler)} "
        '--threads "$SLURM_CPUS_PER_TASK" '
        f"--hash-bins {args.rna3_bam_evidence_hash_bins} "
        f"--hash-seed {DEFAULT_BAM_EVIDENCE_HASH_SEED} "
        f"--max-memory-gb {memory_limit:.3f}"
    )
    logs = run_dir / "logs"
    slurm_dir = control / "slurm"

    diagnostic_audit = control / "rna3_star_diagnostic_promotion.json"
    diagnostic_marker = control / "STAR_DIAGNOSTIC_PROMOTION_COMPLETE.ok"
    diagnostic_script = slurm_dir / "rna3_star_diagnostic_promotion.sbatch"
    diagnostic_body = (
        "module purge\n"
        "module load miniforge/3\n"
        "module list\n"
        f"python3 {q(resources.star_diagnostic_promoter)} \\\n"
        f"    --manifest {q(manifest)} \\\n"
        f"    --work-root {q(root / 'mapping_project' / 'work')} \\\n"
        f"    --audit {q(diagnostic_audit)} \\\n"
        f"    --marker {q(diagnostic_marker)}\n"
        f"test -s {q(diagnostic_audit)}\n"
        f"test -s {q(diagnostic_marker)}\n"
    )
    write_sbatch(
        diagnostic_script,
        sbatch_text(
            "rna3_star_diagnostics",
            logs,
            diagnostic_body,
            cpus=1,
            memory="4G",
            walltime="12:00:00",
        ),
    )
    diagnostic_job = JobSpec(
        "rna3_star_diagnostics",
        diagnostic_script,
        dependencies=list(dependencies),
    )

    # One invocation submits every selected library as a bounded array and then
    # runs the BAM-free gather only after both the array and STAR diagnostics
    # complete successfully.
    production_cap = min(args.rna3_bam_evidence_max_concurrent, 3)
    array_script = slurm_dir / "rna3_bam_evidence.sbatch"
    write_sbatch(
        array_script,
        sbatch_text(
            "rna3_bam_evidence",
            logs,
            reader_preflight
            + worker
            + ' --row-index "$SLURM_ARRAY_TASK_ID"\n',
            cpus=args.rna3_bam_evidence_cpus,
            memory=f"{args.rna3_bam_evidence_memory_gb}G",
            array=f"0-{len(rows) - 1}%{production_cap}",
        ),
    )

    gather_output = run_dir / "qc" / "rna3" / "bam_evidence"
    gather_script = slurm_dir / "rna3_bam_evidence_gather.sbatch"
    gather_body = (
        "module purge\n"
        + "module load miniforge/3\n"
        + "module list\n"
        + f"python3 {q(resources.bam_evidence_runner)} gather \\\n"
        + f"    --manifest {q(manifest)} \\\n"
        + f"    --output-dir {q(gather_output)}\n"
        + f"test -s {q(gather_output / 'BAM_EVIDENCE_GATHER_COMPLETE.ok')}\n"
    )
    write_sbatch(
        gather_script,
        sbatch_text(
            "rna3_bam_evidence_gather",
            logs,
            gather_body,
            cpus=2,
            memory="16G",
        ),
    )
    atomic_json(
        control / "rna3_bam_evidence_execution.json",
        {
            "mode": "direct_full_run",
            "selected_libraries": [str(row["library"]) for row in rows],
            "array": f"0-{len(rows) - 1}%{production_cap}",
            "selected_max_concurrent": production_cap,
            "default_ceiling": 3,
            "automatic_gather_after_all_readers": True,
        },
    )
    return [
        diagnostic_job,
        JobSpec(
            "rna3_bam_evidence",
            array_script,
            dependencies=list(dependencies),
        ),
        JobSpec(
            "rna3_bam_evidence_gather",
            gather_script,
            dependencies=["rna3_star_diagnostics", "rna3_bam_evidence"],
        ),
    ]


def generate_atac_qc_job(
    libraries: Sequence[int],
    atac_base: Path,
    rna_base: Path,
    run_dir: Path,
    args: argparse.Namespace,
    resources: Resources,
    dependencies: Sequence[str],
) -> JobSpec:
    if not libraries:
        raise OrchestratorError("ATAC QC was selected, but no numeric ATAC libraries were found")
    assert resources.atac_collector is not None
    stats = run_dir / "qc" / "atac" / "stats"
    script = run_dir / "control" / "slurm" / "collect_atac_qc.sbatch"
    copy_or_link = (
        'cp "$BAM" "$TASK_TMP/atac.bam"\n'
        'cp "$FRAG" "$TASK_TMP/atac_fragments.tsv.gz"\n'
        'cp "$BARCODES" "$TASK_TMP/barcodes.tsv.gz"'
        if args.atac_qc_copy_inputs
        else
        'ln -s "$BAM" "$TASK_TMP/atac.bam"\n'
        'ln -s "$FRAG" "$TASK_TMP/atac_fragments.tsv.gz"\n'
        'ln -s "$BARCODES" "$TASK_TMP/barcodes.tsv.gz"'
    )
    library_array = " ".join(str(value) for value in libraries)
    body = f"""module purge
module load miniforge/3 genomics-base/latest
module load htslib/1.20
module load samtools/1.20
mkdir -p {q(stats)}
LIBRARIES=({library_array})
LIB="${{LIBRARIES[$SLURM_ARRAY_TASK_ID]}}"
ATAC_PREFIX={q(args.atac_lib_prefix)}
RNA_PREFIX={q(args.rna3_lib_prefix)}
OUT={q(stats)}/library_${{LIB}}_stats.json
if [[ -s "$OUT" ]]; then
    echo "SKIP: ATAC QC library $LIB already completed"
    exit 0
fi
ATAC_DIR={q(atac_base)}/${{ATAC_PREFIX}}${{LIB}}
RNA_DIR={q(rna_base)}/${{RNA_PREFIX}}${{LIB}}
BAM="$ATAC_DIR/atac.bam"
FRAG="$ATAC_DIR/atac_fragments.tsv.gz"
BARCODES="$RNA_DIR/filtered/barcodes.tsv.gz"
for INPUT in "$BAM" "$FRAG" "$BARCODES"; do
    if [[ ! -s "$INPUT" ]]; then
        echo "ERROR: missing required input: $INPUT" >&2
        exit 1
    fi
done
TASK_TMP="${{SLURM_TMPDIR:-/tmp}}/tet_atac_qc_${{SLURM_JOB_ID}}_${{SLURM_ARRAY_TASK_ID}}"
mkdir -p "$TASK_TMP"
cleanup() {{ rm -rf -- "$TASK_TMP"; }}
trap cleanup EXIT
{copy_or_link}
python3 {q(Path(__file__).resolve())} --_run-atac-qc \\
    --_collector {q(resources.atac_collector)} \\
    --_library "$LIB" \\
    --_ramdisk "$TASK_TMP" \\
    --_stats-dir {q(stats)} \\
    --_threads "$SLURM_CPUS_PER_TASK"
test -s "$OUT"
"""
    write_sbatch(
        script,
        sbatch_text(
            "collect_atac_qc",
            run_dir / "logs",
            body,
            cpus=args.atac_qc_cpus,
            memory=args.atac_qc_memory,
            array=f"0-{len(libraries) - 1}",
        ),
    )
    return JobSpec(
        label="atac_qc", script=script, dependencies=list(dependencies)
    )


def generate_atac_plot_job(
    run_dir: Path,
    figure_root: Path,
    resources: Resources,
    dependencies: Sequence[str],
) -> JobSpec:
    assert resources.atac_plotter is not None
    stats = run_dir / "qc" / "atac" / "stats"
    out = figure_root / "mapping_QC"
    marker = out / "ATAC_PLOTS_COMPLETE.ok"
    script = run_dir / "control" / "slurm" / "plot_atac_qc.sbatch"
    body = f"""if [[ -s {q(marker)} ]]; then
    echo "SKIP: ATAC plots already completed"
    exit 0
fi
module purge
module load miniforge/3 genomics-base/latest
export MPLBACKEND=Agg
mkdir -p {q(out)}
python3 {q(resources.atac_plotter)} --stats-dir {q(stats)} --output-dir {q(out)}
date -Is > {q(marker)}
"""
    write_sbatch(
        script,
        sbatch_text(
            "plot_atac_qc",
            run_dir / "logs",
            body,
            cpus=4,
            memory="24G",
            walltime="1-00:00:00",
        ),
    )
    return JobSpec(
        label="atac_plot", script=script, dependencies=list(dependencies)
    )


def rna_library_output_checks(
    modality: str,
    library: str,
    base: Path,
    bam_evidence: bool = False,
) -> list[tuple[str, str, str]]:
    """Return the durable per-library RNA products required before work cleanup."""
    relative_paths: tuple[str, ...] = (
        "gex.bam",
        "gex.bam.bai",
        "Barcodes.stats",
        "Features.stats",
        "Summary.csv",
        "UMIperCellSorted.txt",
        "STAR_Log.out",
        "STAR_Log.final.out",
        "STAR_SJ.out.tab.gz",
        "raw/barcodes.tsv.gz",
        "raw/features.tsv.gz",
        "raw/matrix.mtx.gz",
        "filtered/barcodes.tsv.gz",
        "filtered/features.tsv.gz",
        "filtered/matrix.mtx.gz",
    )
    if bam_evidence:
        relative_paths += (
            "bam_evidence/barcode_read_metrics.tsv.gz",
            "bam_evidence/barcode_rg_metrics.tsv.gz",
            "bam_evidence/molecule_source_hash_bins.tsv.gz",
            "bam_evidence/rg_summary.tsv",
            "bam_evidence/rg_contig_class_summary.tsv.gz",
            "bam_evidence/raw_to_corrected_barcode_counts.tsv.gz",
            "bam_evidence/CellReads.countedU.from_bam.tsv.gz",
            "bam_evidence/audit.json",
            "bam_evidence/BAM_EVIDENCE_COMPLETE.ok",
        )
    else:
        relative_paths += ("CellReads.stats.gz",)
    return [
        (modality, library, str(base / relative_path))
        for relative_path in relative_paths
    ]


def nextflow_project_output_checks(
    modality: str,
    project: Path,
) -> list[tuple[str, str, str]]:
    """Return durable execution metadata needed after Nextflow cache removal."""
    relative_paths = (
        "MAPPING_COMPLETE.ok",
        "rg_metadata.tsv",
        "symlink_manifest.tsv",
        "params_rna.yml" if modality.startswith("rna") else "params_atac.yml",
        "nextflow.config",
        "report.html",
        "trace.txt",
        "timeline.html",
    )
    return [
        (f"{modality}_provenance", relative_path, str(project / relative_path))
        for relative_path in relative_paths
    ]


def nextflow_work_cleanup_targets(
    inputs: Inputs,
    discovery: Discovery,
    run_dir: Path,
) -> list[tuple[str, str]]:
    """List cache directories that may be removed manually after validation."""
    targets: list[tuple[str, str]] = []
    if inputs.rna3:
        targets.append(("rna3", str(run_dir / "rna3" / "mapping_project" / "work")))
    if inputs.rna5:
        for read_format in ("long-r2", "pe150"):
            if discovery.rna5_libraries_by_format.get(read_format):
                targets.append(
                    (
                        f"rna5_{read_format}",
                        str(
                            run_dir
                            / "rna5"
                            / "formats"
                            / read_format
                            / "mapping_project"
                            / "work"
                        ),
                    )
                )
    if inputs.atac:
        targets.append(("atac", str(run_dir / "atac" / "mapping_project" / "work")))
    return targets


def validation_checks(
    inputs: Inputs,
    discovery: Discovery,
    run_dir: Path,
    stages: set[str],
    rna5_map_mode: str,
    include_trim_barcode_qc: bool,
    include_empty_drop_qc: bool,
    rna3_bam_evidence: bool,
) -> list[tuple[str, str, str]]:
    checks: list[tuple[str, str, str]] = []
    if inputs.rna3:
        for library in discovery.libraries.get("rna3", []):
            base = run_dir / "rna3" / "mapping_output" / library
            checks.extend(
                rna_library_output_checks(
                    "rna3",
                    library,
                    base,
                    bam_evidence=rna3_bam_evidence,
                )
            )
        checks.extend(
            nextflow_project_output_checks(
                "rna3", run_dir / "rna3" / "mapping_project"
            )
        )
        checks.append(
            ("rna3_provenance", "libs.txt", str(run_dir / "rna3" / "libs.txt"))
        )
        checks.append(
            (
                "rna3_provenance",
                "run_mapping.sbatch",
                str(run_dir / "rna3" / "run_mapping.sbatch"),
            )
        )
        if rna3_bam_evidence:
            evidence_gather = run_dir / "qc" / "rna3" / "bam_evidence"
            for relative in (
                "bam_evidence_inventory.tsv",
                "all_libraries_rg_summary.tsv.gz",
                "all_libraries_barcode_summary.tsv.gz",
                "all_libraries_source_yield.tsv",
                "project_audit.json",
                "BAM_EVIDENCE_GATHER_COMPLETE.ok",
            ):
                checks.append(
                    ("rna3_bam_evidence_gather", relative, str(evidence_gather / relative))
                )
        if include_trim_barcode_qc:
            checks.extend(
                [
                    ("rna3_trim_by_barcode", "all", str(run_dir / "qc" / "rna3" / "trim_by_barcode" / "all_libraries_trim_by_barcode.tsv.gz")),
                    ("rna3_trim_by_barcode", "adapters", str(run_dir / "qc" / "rna3" / "trim_by_barcode" / "all_libraries_trim_adapters_by_barcode.tsv.gz")),
                    ("rna3_trim_by_barcode", "by_source", str(run_dir / "qc" / "rna3" / "trim_by_barcode" / "all_libraries_trim_by_barcode_by_source.tsv.gz")),
                    ("rna3_trim_by_barcode", "filtered_cells", str(run_dir / "qc" / "rna3" / "trim_by_barcode" / "all_libraries_trim_by_barcode_filtered_cells.tsv.gz")),
                    ("rna3_trim_by_barcode", "observed_not_filtered", str(run_dir / "qc" / "rna3" / "trim_by_barcode" / "all_libraries_trim_by_barcode_observed_not_filtered.tsv.gz")),
                ]
            )
            if include_empty_drop_qc:
                checks.append(
                    ("rna3_trim_by_barcode", "empty_drops", str(run_dir / "qc" / "rna3" / "trim_by_barcode" / "all_libraries_trim_by_barcode_empty_drops.tsv.gz"))
                )
    if inputs.rna5:
        for read_format in ("long-r2", "pe150"):
            for library in discovery.rna5_libraries_by_format.get(read_format, []):
                if rna5_map_mode == "together":
                    base = run_dir / "rna5" / "mapping_output" / library
                else:
                    base = run_dir / "rna5" / "formats" / read_format / "mapping_output" / library
                checks.extend(
                    rna_library_output_checks(
                        f"rna5_{read_format}", library, base
                    )
                )
            format_project = run_dir / "rna5" / "formats" / read_format / "mapping_project"
            if discovery.rna5_libraries_by_format.get(read_format):
                checks.extend(
                    nextflow_project_output_checks(
                        f"rna5_{read_format}", format_project
                    )
                )
                checks.append(
                    (
                        f"rna5_{read_format}_provenance",
                        "libs.txt",
                        str(
                            run_dir
                            / "rna5"
                            / "formats"
                            / read_format
                            / "libs.txt"
                        ),
                    )
                )
                checks.append(
                    (
                        f"rna5_{read_format}_provenance",
                        "run_mapping.sbatch",
                        str(
                            run_dir
                            / "rna5"
                            / "formats"
                            / read_format
                            / "run_mapping.sbatch"
                        ),
                    )
                )
        if include_trim_barcode_qc:
            checks.extend(
                [
                    ("rna5_trim_by_barcode", "all", str(run_dir / "qc" / "rna5" / "trim_by_barcode" / "all_libraries_trim_by_barcode.tsv.gz")),
                    ("rna5_trim_by_barcode", "adapters", str(run_dir / "qc" / "rna5" / "trim_by_barcode" / "all_libraries_trim_adapters_by_barcode.tsv.gz")),
                    ("rna5_trim_by_barcode", "by_source", str(run_dir / "qc" / "rna5" / "trim_by_barcode" / "all_libraries_trim_by_barcode_by_source.tsv.gz")),
                    ("rna5_trim_by_barcode", "filtered_cells", str(run_dir / "qc" / "rna5" / "trim_by_barcode" / "all_libraries_trim_by_barcode_filtered_cells.tsv.gz")),
                    ("rna5_trim_by_barcode", "observed_not_filtered", str(run_dir / "qc" / "rna5" / "trim_by_barcode" / "all_libraries_trim_by_barcode_observed_not_filtered.tsv.gz")),
                ]
            )
            if include_empty_drop_qc:
                checks.append(
                    ("rna5_trim_by_barcode", "empty_drops", str(run_dir / "qc" / "rna5" / "trim_by_barcode" / "all_libraries_trim_by_barcode_empty_drops.tsv.gz"))
                )
    if inputs.atac:
        for library in discovery.libraries.get("atac", []):
            base = run_dir / "atac" / "mapping_output" / library
            checks.extend(
                [
                    ("atac", library, str(base / "atac.bam")),
                    ("atac", library, str(base / "atac.bam.bai")),
                    ("atac", library, str(base / "atac_namesort.bam")),
                    ("atac", library, str(base / "atac_fragments.tsv.gz")),
                    ("atac", library, str(base / "atac_fragments.tsv.gz.tbi")),
                ]
            )
        checks.extend(
            nextflow_project_output_checks(
                "atac", run_dir / "atac" / "mapping_project"
            )
        )
        checks.append(
            ("atac_provenance", "libs.txt", str(run_dir / "atac" / "libs.txt"))
        )
        checks.append(
            (
                "atac_provenance",
                "run_atac_mapping.sbatch",
                str(run_dir / "atac" / "run_atac_mapping.sbatch"),
            )
        )
        if "qc" in stages:
            for number in discovery.library_numbers.get("atac", []):
                checks.append(
                    (
                        "atac_qc",
                        str(number),
                        str(run_dir / "qc" / "atac" / "stats" / f"library_{number}_stats.json"),
                    )
                )
    return checks


def generate_validation_job(
    inputs: Inputs,
    discovery: Discovery,
    run_dir: Path,
    stages: set[str],
    dependencies: Sequence[str],
    rna5_map_mode: str,
    include_trim_barcode_qc: bool,
    include_empty_drop_qc: bool,
    rna3_bam_evidence: bool,
) -> JobSpec:
    checks = validation_checks(
        inputs,
        discovery,
        run_dir,
        stages,
        rna5_map_mode,
        include_trim_barcode_qc,
        include_empty_drop_qc,
        rna3_bam_evidence,
    )
    check_file = run_dir / "control" / "expected_outputs.tsv"
    atomic_text(
        check_file,
        "modality\tlibrary\tpath\n"
        + "".join(f"{modality}\t{library}\t{path}\n" for modality, library, path in checks),
    )
    cleanup_plan = run_dir / "control" / "work_cleanup_targets.tsv"
    cleanup_targets = nextflow_work_cleanup_targets(inputs, discovery, run_dir)
    atomic_text(
        cleanup_plan,
        "modality\tpath\taction\tconsequence\n"
        + "".join(
            f"{modality}\t{path}\tMANUAL_DELETE_ONLY_AFTER_WORK_CLEANUP_READY\t"
            "removes_nextflow_resume_cache_not_published_results\n"
            for modality, path in cleanup_targets
        ),
    )
    summary = run_dir / "validation" / "output_validation.tsv"
    marker = run_dir / "validation" / "RUN_COMPLETE.ok"
    cleanup_manifest = run_dir / "validation" / "work_cleanup_targets.tsv"
    cleanup_marker = run_dir / "validation" / "WORK_CLEANUP_READY.ok"
    script = run_dir / "control" / "slurm" / "validate_run.sbatch"
    body = f"""module purge
module load htslib/1.20
module load samtools/1.20
mkdir -p {q(run_dir / 'validation')}
TMP={q(summary)}.tmp.${{SLURM_JOB_ID}}
printf 'modality\\tlibrary\\tstatus\\tpath\\n' > "$TMP"
FAIL=0
while IFS=$'\\t' read -r MODALITY LIBRARY PATH_VALUE; do
    [[ "$MODALITY" == "modality" ]] && continue
    if [[ ! -s "$PATH_VALUE" ]]; then
        printf '%s\\t%s\\tMISSING\\t%s\\n' "$MODALITY" "$LIBRARY" "$PATH_VALUE" >> "$TMP"
        FAIL=1
        continue
    fi
    if [[ "$PATH_VALUE" == *.bam ]] && ! samtools quickcheck "$PATH_VALUE"; then
        printf '%s\\t%s\\tBAM_QUICKCHECK_FAILED\\t%s\\n' "$MODALITY" "$LIBRARY" "$PATH_VALUE" >> "$TMP"
        FAIL=1
        continue
    fi
    if [[ "$PATH_VALUE" == *.bam ]] && ! samtools view -H "$PATH_VALUE" | grep -q '^@RG'; then
        printf '%s\\t%s\\tBAM_RG_HEADER_MISSING\\t%s\\n' "$MODALITY" "$LIBRARY" "$PATH_VALUE" >> "$TMP"
        FAIL=1
        continue
    fi
    printf '%s\\t%s\\tOK\\t%s\\n' "$MODALITY" "$LIBRARY" "$PATH_VALUE" >> "$TMP"
done < {q(check_file)}
mv "$TMP" {q(summary)}
if [[ "$FAIL" -ne 0 ]]; then
    echo "Validation failed; see {summary}" >&2
    exit 1
fi
date -Is > {q(marker)}
cp {q(cleanup_plan)} {q(cleanup_manifest)}
{{
    printf 'validated_utc\\t'
    date -Is
    printf 'cleanup_manifest\\t%s\\n' {q(str(cleanup_manifest))}
    printf 'policy\\tmanual deletion only; published results are retained; Nextflow resume cache is lost\\n'
}} > {q(cleanup_marker)}
echo "Validation passed: {marker}"
echo "Manual work cleanup is safe: {cleanup_marker}"
"""
    write_sbatch(
        script,
        sbatch_text(
            "validate_10x_run",
            run_dir / "logs",
            body,
            cpus=2,
            memory="8G",
            walltime="1-00:00:00",
        ),
    )
    return JobSpec(
        label="validate", script=script, dependencies=list(dependencies)
    )


def submit_job(script: Path, dependency_ids: Sequence[str]) -> str:
    command = ["sbatch", "--parsable"]
    if dependency_ids:
        command.extend(["--dependency", "afterok:" + ":".join(dependency_ids)])
    command.append(str(script))
    last_detail = ""
    for attempt in (1, 2):
        try:
            result = subprocess.run(command, capture_output=True, text=True, check=False)
        except OSError as exc:
            result = None
            last_detail = str(exc)
        if result is not None and result.returncode == 0:
            job_id = result.stdout.strip().split(";", 1)[0]
            if job_id.isdigit():
                return job_id
            last_detail = f"unexpected sbatch response: {result.stdout.strip()!r}"
        elif result is not None:
            last_detail = result.stderr.strip() or result.stdout.strip()
        if attempt == 1:
            print(f"WARNING: sbatch failed for {script}: {last_detail}; retrying once")
            time.sleep(1)
    raise OrchestratorError(f"sbatch did not submit {script}: {last_detail}")


SLURM_ACTIVE_STATES = {
    "CONFIGURING",
    "COMPLETING",
    "PENDING",
    "REQUEUED",
    "REQUEUE_FED",
    "REQUEUE_HOLD",
    "RESIZING",
    "RUNNING",
    "SIGNALING",
    "STAGE_OUT",
    "STOPPED",
    "SUSPENDED",
}
SLURM_FAILURE_STATES = {
    "BOOT_FAIL",
    "CANCELLED",
    "DEADLINE",
    "FAILED",
    "NODE_FAIL",
    "OUT_OF_MEMORY",
    "PREEMPTED",
    "REVOKED",
    "SPECIAL_EXIT",
    "TIMEOUT",
}


def slurm_accounting_states(job_id: str) -> list[str]:
    """Return normalized parent/array states retained by SLURM accounting."""
    command = ["sacct", "-X", "-n", "-P", "-j", job_id, "-o", "State"]
    try:
        result = subprocess.run(command, capture_output=True, text=True, check=False)
    except OSError:
        return []
    if result.returncode != 0:
        return []
    states: list[str] = []
    for line in result.stdout.splitlines():
        state = line.strip().split("|", 1)[0].split(None, 1)[0].rstrip("+")
        if state:
            states.append(state)
    return states


def reused_job_dependency(job_id: str, label: str) -> str | None:
    """Decide whether a reused job ID still belongs in a dependency list.

    Completed jobs have already satisfied their dependency and must not be sent
    back to ``sbatch`` after they leave SLURM's active-controller memory. Active
    jobs remain dependencies. A known failed job is never silently reused.
    """
    states = slurm_accounting_states(job_id)
    # These jobs publish scientific/provenance products consumed by later
    # RNA-evidence phases.  When accounting no longer proves what happened,
    # require an explicit resubmission instead of silently releasing their
    # dependants.  The replacement job then recreates and validates its marker
    # and audit in the normal job body.
    requires_accounting_proof = (
        label.startswith("rna3_bam_evidence")
        or label == "rna3_star_diagnostics"
    )
    if states and all(state == "COMPLETED" for state in states):
        print(f"Prior submission already completed: {label} -> job {job_id}")
        return None

    failed = sorted({state for state in states if state in SLURM_FAILURE_STATES})
    if failed:
        raise OrchestratorError(
            f"cannot reuse {label} job {job_id}: SLURM state is "
            f"{','.join(failed)}. Add {label} to --resubmit-jobs."
        )

    if any(state in SLURM_ACTIVE_STATES for state in states):
        print(f"Prior submission is still active: {label} -> job {job_id}")
        return job_id

    if states:
        if requires_accounting_proof:
            raise OrchestratorError(
                f"cannot reuse {label} job {job_id}: SLURM returned only "
                f"unrecognized state(s) {','.join(sorted(set(states)))}. "
                f"Add {label} to --resubmit-jobs."
            )
        print(
            f"WARNING: unrecognized SLURM state for reused job {label} "
            f"({job_id}): {','.join(sorted(set(states)))}; treating its "
            "dependency as already satisfied"
        )
    else:
        if requires_accounting_proof:
            raise OrchestratorError(
                f"cannot reuse {label} job {job_id}: SLURM accounting no "
                f"longer proves its state. Add {label} to --resubmit-jobs."
            )
        print(
            f"WARNING: SLURM no longer retains accounting state for reused "
            f"job {label} ({job_id}); treating its dependency as already "
            "satisfied"
        )
    return None


def validate_completed_star_diagnostic_reuse(
    promoter: Path,
    run_dir: Path,
) -> None:
    """Use the promoter's authoritative read-only audit validator on reuse."""
    control = run_dir / "control"
    command = [
        sys.executable,
        str(promoter),
        "--validate-only",
        "--manifest", str(control / "rna3_bam_evidence_manifest.tsv"),
        "--work-root", str(run_dir / "rna3" / "mapping_project" / "work"),
        "--audit", str(control / "rna3_star_diagnostic_promotion.json"),
        "--marker", str(control / "STAR_DIAGNOSTIC_PROMOTION_COMPLETE.ok"),
    ]
    try:
        result = subprocess.run(
            command, capture_output=True, text=True, check=False
        )
    except OSError as exc:
        raise OrchestratorError(
            f"could not validate reused STAR diagnostics: {exc}"
        ) from exc
    if result.returncode != 0:
        detail = result.stderr.strip() or result.stdout.strip()
        raise OrchestratorError(
            "completed rna3_star_diagnostics job lacks a current valid "
            "manifest-bound audit/marker/output set; add "
            "rna3_star_diagnostics to --resubmit-jobs. "
            f"Validator detail: {detail}"
        )


def plan_payload(
    run_dir: Path, stages: set[str], jobs: Sequence[JobSpec], submitted: bool
) -> dict[str, object]:
    return {
        "release": RELEASE,
        "updated_utc": utc_now(),
        "run_directory": str(run_dir),
        "stages": [stage for stage in ALL_STAGES if stage in stages],
        "submitted": submitted,
        "production_promotion": "disabled",
        "jobs": [
            {
                "label": job.label,
                "script": str(job.script),
                "dependencies": job.dependencies,
                "job_id": job.job_id,
            }
            for job in jobs
        ],
    }


def print_plan(
    run_dir: Path,
    jobs: Sequence[JobSpec],
    submitted: bool,
    max_cores: int | None,
    nodelist: str | None,
    array_max_concurrent: int | None,
) -> None:
    print("\n" + "=" * 72)
    print("10X TRIMMING / MAPPING / QC PLAN")
    print(f"  Run directory: {run_dir}")
    print(f"  Mode: {'SUBMITTED' if submitted else 'DRY RUN (scripts only)'}")
    print("  Production mapping_output overwrite: DISABLED")
    if max_cores is not None:
        print(f"  Hard RNA core ceiling: {max_cores}")
    if nodelist is not None:
        print(f"  SLURM node list: {nodelist}")
    if array_max_concurrent is not None:
        print(f"  Maximum simultaneous array tasks: {array_max_concurrent}")
    print("=" * 72)
    for job in jobs:
        deps = ",".join(job.dependencies) if job.dependencies else "none"
        suffix = f" -> job {job.job_id}" if job.job_id else ""
        print(f"  {job.label}: {job.script} [after: {deps}]{suffix}")
    if not submitted:
        print("\nNo jobs were submitted. Re-run this exact command with --submit --resume")
        print("after inspecting the generated scripts.")
    else:
        print("\nMonitor with: squeue -u $USER")
    print(f"Validation marker: {run_dir / 'validation' / 'RUN_COMPLETE.ok'}")
    print(
        "Work-cleanup marker: "
        f"{run_dir / 'validation' / 'WORK_CLEANUP_READY.ok'}"
    )


def restore_prior_job_ids(
    plan_path: Path,
    jobs: Sequence[JobSpec],
    requested_resubmits: Sequence[str] | None,
    star_diagnostic_promoter: Path | None = None,
) -> None:
    """Reuse accepted SLURM IDs and clear explicitly retried dependency branches."""
    if not plan_path.is_file():
        if requested_resubmits:
            raise OrchestratorError(
                "--resubmit-jobs requires an existing job_plan.json with submitted jobs"
            )
        return
    try:
        prior = json.loads(plan_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise OrchestratorError(f"could not read existing job plan {plan_path}: {exc}")

    prior_ids = {
        item.get("label"): str(item.get("job_id"))
        for item in prior.get("jobs", [])
        if item.get("label") and item.get("job_id")
    }
    labels = {job.label for job in jobs}
    requested = set(requested_resubmits or [])
    if "all" in requested:
        requested = set(labels)
    unknown = sorted(requested - labels)
    if unknown:
        raise OrchestratorError(
            "--resubmit-jobs contains unknown label(s): " + ", ".join(unknown)
        )
    for label in sorted(requested):
        prior_id = prior_ids.get(label)
        protected = (
            label.startswith("rna3_bam_evidence")
            or label == "rna3_star_diagnostics"
        )
        if not prior_id or not protected:
            continue
        states = slurm_accounting_states(prior_id)
        active = sorted(set(states).intersection(SLURM_ACTIVE_STATES))
        unknown_states = sorted(
            set(states) - SLURM_ACTIVE_STATES - SLURM_FAILURE_STATES
            - {"COMPLETED"}
        )
        if active:
            raise OrchestratorError(
                f"cannot resubmit {label} while prior job {prior_id} is "
                f"active: {','.join(active)}"
            )
        if not states or unknown_states:
            detail = "no state" if not states else ",".join(unknown_states)
            raise OrchestratorError(
                f"cannot resubmit {label}: SLURM did not prove prior job "
                f"{prior_id} terminal ({detail})"
            )

    diagnostic_id = prior_ids.get("rna3_star_diagnostics")
    if (
        diagnostic_id
        and "rna3_star_diagnostics" in labels
        and "rna3_star_diagnostics" not in requested
    ):
        diagnostic_states = slurm_accounting_states(diagnostic_id)
        if diagnostic_states and all(
            state == "COMPLETED" for state in diagnostic_states
        ):
            if star_diagnostic_promoter is None:
                raise OrchestratorError(
                    "cannot validate completed rna3_star_diagnostics reuse: "
                    "promoter path is unavailable"
                )
            validate_completed_star_diagnostic_reuse(
                star_diagnostic_promoter, plan_path.parent.parent
            )

    # A requested retry or a newly introduced job invalidates every prior
    # downstream submission in this plan. This prevents an old completed plot
    # or validation job from being reused when a new reporting stage is added
    # during --resume.
    retry = set(requested)
    retry.update(label for label in labels if label not in prior_ids)
    changed = True
    while changed:
        changed = False
        for job in jobs:
            if job.label not in retry and retry.intersection(job.dependencies):
                retry.add(job.label)
                changed = True

    for job in jobs:
        if job.label not in retry and job.label in prior_ids:
            job.job_id = prior_ids[job.label]


def gather_trim_json(trim_root: Path, report_dir: Path) -> None:
    if not trim_root.is_dir():
        raise OrchestratorError(f"trimming directory does not exist: {trim_root}")
    report_dir.mkdir(parents=True, exist_ok=True)
    reports = sorted(trim_root.glob("*/logs/*_cutadapt.json"))
    if not reports:
        raise OrchestratorError(f"no cutadapt JSON reports found below {trim_root}")
    for source in reports:
        run_id = source.parent.parent.name
        suffix = "_cutadapt.json"
        stem = source.name[: -len(suffix)] if source.name.endswith(suffix) else source.stem
        destination = report_dir / f"{stem}__{run_id}{suffix}"
        if destination.is_symlink():
            if destination.resolve() == source.resolve():
                continue
            raise OrchestratorError(f"report-link collision: {destination}")
        if destination.exists():
            raise OrchestratorError(f"refusing to overwrite report path: {destination}")
        destination.symlink_to(source.resolve())
    print(f"Gathered {len(reports)} cutadapt JSON reports into {report_dir}")


def run_atac_qc_worker(args: argparse.Namespace) -> None:
    collector = Path(args._collector).resolve()
    spec = importlib.util.spec_from_file_location("tet_atac_qc_collector", collector)
    if spec is None or spec.loader is None:
        raise OrchestratorError(f"could not import ATAC QC collector: {collector}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    stats = module.collect_library_stats(args._library, args._ramdisk, args._threads)
    stats_dir = Path(args._stats_dir)
    stats_dir.mkdir(parents=True, exist_ok=True)
    destination = stats_dir / f"library_{args._library}_stats.json"
    temp = destination.with_name(f".{destination.name}.tmp.{os.getpid()}")
    with temp.open("w", encoding="utf-8") as handle:
        json.dump(stats, handle, indent=2)
        handle.write("\n")
    os.replace(temp, destination)
    print(f"Saved stats atomically to {destination}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Orchestrate staged 10X RNA/ATAC trimming, mapping, QC, and plots",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "BP names are expanded below the assay-specific raw roots. Absolute "
            "input paths are also accepted. Dry run is the default; --submit is explicit."
        ),
    )
    parser.add_argument(
        "--run-name",
        help=("Unique run label. Required unless --run-dir is supplied; with "
              "--run-dir it must match that directory's basename."),
    )
    parser.add_argument(
        "--run-dir", default=None,
        help=("Exact absolute output directory. This bypasses the historical "
              "<staging-root>/<run-name> construction."),
    )
    parser.add_argument("--rna3-runs", nargs="+", default=None, metavar="RUN")
    parser.add_argument("--atac-runs", nargs="+", default=None, metavar="RUN")
    parser.add_argument("--rna5-runs", nargs="+", default=None, metavar="RUN")
    parser.add_argument(
        "--libraries",
        nargs="+",
        type=int,
        default=None,
        metavar="N",
        help=(
            "Process only these library numbers across the supplied multiplexed "
            "run directories (for example: --libraries 19)"
        ),
    )
    parser.add_argument(
        "--rna3-lib-prefix",
        default=DEFAULT_RNA3_LIB_PREFIX,
        help="Prefix before the numeric 3' RNA library ID",
    )
    parser.add_argument(
        "--rna5-lib-prefix",
        default=DEFAULT_RNA5_LIB_PREFIX,
        help="Prefix before the numeric 5' RNA library ID",
    )
    parser.add_argument(
        "--atac-lib-prefix",
        default=DEFAULT_ATAC_LIB_PREFIX,
        help="Prefix before the numeric ATAC library ID",
    )
    parser.add_argument("--rna3-raw-root", default=DEFAULT_RNA3_RAW_ROOT)
    parser.add_argument("--atac-raw-root", default=DEFAULT_ATAC_RAW_ROOT)
    parser.add_argument("--rna5-raw-root", default=DEFAULT_RNA5_RAW_ROOT)
    parser.add_argument("--staging-root", default=DEFAULT_STAGING_ROOT)
    parser.add_argument(
        "--rna3-figure-root", default=None,
        help=("Exact physical root for 3P mapping figures; default: "
              "3P/figures/<run-name>/Mapping"),
    )
    parser.add_argument(
        "--rna5-figure-root", default=None,
        help=("Exact physical root for 5P mapping figures; default: "
              "5P/figures/<run-name>/Mapping"),
    )
    parser.add_argument(
        "--atac-figure-root", default=None,
        help=("Exact physical root for ATAC mapping figures; default: "
              "ATAC/figures/<run-name>/Mapping"),
    )
    parser.add_argument(
        "--stages",
        nargs="+",
        default=["all"],
        metavar="STAGE",
        help="all, or any of: trim map qc plot validate (comma-separated also accepted)",
    )
    parser.add_argument("--submit", action="store_true", help="Submit generated jobs to SLURM")
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Continue an existing run only when its immutable configuration matches",
    )
    parser.add_argument(
        "--resubmit-jobs",
        nargs="+",
        default=None,
        metavar="LABEL",
        help=(
            "With --resume --submit, resubmit named failed job labels (or all) "
            "and automatically resubmit their downstream dependents"
        ),
    )
    parser.add_argument("--diagnose-only", action="store_true", help="Read-only preflight report")

    parser.add_argument("--rna-ref", default=DEFAULT_RNA_REF)
    parser.add_argument("--atac-ref", default=DEFAULT_ATAC_REF)
    parser.add_argument("--rna3-whitelist", default=DEFAULT_RNA3_WHITELIST)
    parser.add_argument(
        "--rna5-whitelist",
        default=DEFAULT_RNA5_WHITELIST,
        help="5' RNA barcode whitelist (the documented default is RNA-737K-arc-v1)",
    )
    parser.add_argument("--atac-whitelist", default=DEFAULT_ATAC_WHITELIST)
    parser.add_argument(
        "--rna-barcode-base",
        default=None,
        help=(
            "RNA mapping_output used by ATAC QC. Default: this staged 3' run "
            "when selected, otherwise existing production RNA output."
        ),
    )

    parser.add_argument(
        "--resource-root",
        default=str(MAPPING_SCRIPT_DIR),
        help=(
            "Directory containing the bundled helper scripts; defaults to the "
            "orchestrator's own deployed directory"
        ),
    )
    parser.add_argument("--rna-driver", default=None)
    parser.add_argument("--atac-driver", default=None)
    parser.add_argument(
        "--rna-workflow",
        default=DEFAULT_RNA_WORKFLOW,
        help=(
            "RNA workflow in this align_pipelines repository/package; default: "
            f"{DEFAULT_RNA_WORKFLOW}"
        ),
    )
    parser.add_argument(
        "--atac-workflow",
        default=DEFAULT_ATAC_WORKFLOW,
        help=(
            "ATAC workflow in this align_pipelines repository/package; default: "
            f"{DEFAULT_ATAC_WORKFLOW}"
        ),
    )
    parser.add_argument("--barcode-aggregator", default=None)
    parser.add_argument("--trim-plotter", default=None)
    parser.add_argument("--mapping-plotter", default=None)
    parser.add_argument("--atac-collector", default=None)
    parser.add_argument("--atac-plotter", default=None)
    parser.add_argument(
        "--rna3-bam-evidence-runner",
        default=None,
        help="Override the bundled run_rna_bam_evidence.py helper",
    )
    parser.add_argument(
        "--rna3-bam-evidence-profiler",
        default=None,
        help=(
            "Override the compiled rna_bam_evidence executable; the default is "
            "the package bin/ sibling or repository-root build"
        ),
    )
    parser.add_argument(
        "--rna3-star-diagnostic-promoter",
        default=None,
        help="Override the bundled promote_rna_star_diagnostics.py helper",
    )

    parser.add_argument(
        "--rna3-mapping-baseline-stats",
        default=None,
        help=(
            "Prior RNA processedstats.tsv. When supplied, the RNA mapping plot "
            "job also writes mapping_delta_dashboard.png and mapping_stats_deltas.tsv."
        ),
    )
    parser.add_argument(
        "--rna3-mapping-baseline-label",
        default="Baseline",
        help="Display label for --rna3-mapping-baseline-stats",
    )
    parser.add_argument(
        "--rna3-mapping-current-label",
        default="Current",
        help="Display label for the current 3' RNA mapping statistics",
    )
    parser.add_argument(
        "--rna3-mapping-baseline-root",
        default=None,
        help=(
            "Optional prior RNA mapping_output directory. When supplied, the "
            "reads-per-cell figure compares native called-cell distributions and "
            "the per-cell change for shared library/barcode pairs."
        ),
    )
    parser.add_argument(
        "--rna3-mapping-baseline-work-root",
        default=None,
        help=(
            "Optional prior Nextflow work directory used to recover baseline "
            "CellReads.stats only when that prior STAR run generated it with "
            "--soloCellReadStats Standard"
        ),
    )
    parser.add_argument(
        "--rna3-bam-evidence-from-bam",
        action="store_true",
        help=(
            "Opt in to the compiled one-pass per-barcode/RG/molecule evidence "
            "profiler, bounded all-library array, and BAM-free gather"
        ),
    )
    parser.add_argument(
        "--rna3-cell-reads-backfill-from-bam",
        action="store_true",
        help=(
            "Deprecated one-release alias for --rna3-bam-evidence-from-bam; "
            "it invokes the compiled profiler and never the retired AWK worker"
        ),
    )
    parser.add_argument(
        "--rna3-bam-evidence-cpus",
        type=int,
        default=4,
        metavar="N",
        help="CPUs per BAM evidence task (default: 4)",
    )
    parser.add_argument(
        "--rna3-bam-evidence-memory-gb",
        type=int,
        default=48,
        metavar="N",
        help="SLURM memory in GiB per BAM evidence task (default: 48)",
    )
    parser.add_argument(
        "--rna3-bam-evidence-max-concurrent",
        type=int,
        default=3,
        metavar="N",
        help=(
            "Maximum simultaneous large BAM readers in the full array "
            "(default: 3)"
        ),
    )
    parser.add_argument(
        "--rna3-bam-evidence-reset-failed-run",
        dest="rna3_bam_evidence_reset_failed_run",
        action="store_true",
        help=(
            "After all prior RNA-evidence jobs are terminal, archive their "
            "outputs and obsolete frozen configuration before a clean full rerun"
        ),
    )
    parser.add_argument(
        "--rna3-bam-evidence-baseline-root",
        default=None,
        help="Optional historical mapping_output root for exact old barcode rosters",
    )
    parser.add_argument(
        "--rna3-bam-evidence-source-order",
        nargs="+",
        default=None,
        metavar="SOURCE",
        help=(
            "Explicit chronological bp_id order. Default: preserve the original "
            "--rna3-runs argument order"
        ),
    )
    parser.add_argument(
        "--rna3-bam-evidence-hash-bins",
        type=int,
        default=100,
        metavar="N",
        help="Nested stable minimum-QNAME-hash bins (default: 100)",
    )
    parser.add_argument(
        "--rna3-bam-evidence-class-manifest",
        default=None,
        help=(
            "TSV with validated contig and/or feature/GX biological classes"
        ),
    )
    parser.add_argument(
        "--rna3-bam-evidence-no-biological-classification",
        action="store_true",
        help=(
            "Explicitly acknowledge that mitochondrial/rRNA/species values will "
            "be unavailable. One of this flag or --rna3-bam-evidence-class-manifest "
            "is mandatory."
        ),
    )

    parser.add_argument("--rna-mem-gb", default="80")
    parser.add_argument("--rna-threads", type=int, default=8)
    parser.add_argument(
        "--max-cores",
        type=int,
        default=None,
        help=(
            "Hard total core ceiling for an RNA-only run. RNA jobs are serialized, "
            "SLURM arrays are throttled, and Nextflow mapping concurrency is capped."
        ),
    )
    parser.add_argument(
        "--nodelist",
        default=None,
        help=(
            "Restrict every generated SLURM job, including Nextflow mapping "
            "workers, to this SLURM host-list expression (for example: char "
            "or char,pika)"
        ),
    )
    parser.add_argument(
        "--array-max-concurrent",
        type=int,
        default=None,
        metavar="N",
        help=(
            "Maximum number of tasks from any generated SLURM array that may "
            "run simultaneously; combines with --max-cores using the stricter limit"
        ),
    )
    parser.add_argument("--atac-mem-gb", default="100")
    parser.add_argument("--atac-threads", type=int, default=8)
    parser.add_argument("--atac-chunks", type=int, default=4)
    parser.add_argument("--no-trim-info", action="store_true")
    parser.add_argument(
        "--empty-drop-roster",
        default=None,
        help=(
            "Optional TSV(.gz) with library and barcode/cell_barcode columns. "
            "It creates an explicit empty-drop trimming view; barcodes absent "
            "from STARsolo filtered output are never assumed to be empty drops."
        ),
    )
    parser.add_argument(
        "--rna5-map-mode",
        choices=("separate", "together"),
        default="separate",
        help=(
            "separate maps long-R2 and PE150 to format-specific outputs; together "
            "creates one output view only when their library sets do not overlap"
        ),
    )
    parser.add_argument(
        "--min-pe150-tso-match-fraction",
        type=float,
        default=0.80,
        help="Minimum fraction of PE150 R1 reads recognizing the expected TSO at base 29",
    )
    parser.add_argument("--atac-qc-cpus", type=int, default=16)
    parser.add_argument("--atac-qc-memory", default="250G")
    parser.add_argument(
        "--atac-qc-copy-inputs",
        action="store_true",
        help="Copy BAM/fragments/barcodes to task-local temp instead of symlinking",
    )

    # Internal workers used only by generated sbatch scripts.
    parser.add_argument("--_gather-trim-json", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--_trim-root", default=None, help=argparse.SUPPRESS)
    parser.add_argument("--_report-dir", default=None, help=argparse.SUPPRESS)
    parser.add_argument("--_run-atac-qc", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--_collector", default=None, help=argparse.SUPPRESS)
    parser.add_argument("--_library", type=int, default=None, help=argparse.SUPPRESS)
    parser.add_argument("--_ramdisk", default=None, help=argparse.SUPPRESS)
    parser.add_argument("--_stats-dir", default=None, help=argparse.SUPPRESS)
    parser.add_argument("--_threads", type=int, default=4, help=argparse.SUPPRESS)
    return parser


def orchestrate(args: argparse.Namespace) -> int:
    requested_run_dir = None
    if args.run_dir:
        requested_run_dir = Path(args.run_dir).expanduser()
        if not requested_run_dir.is_absolute():
            raise OrchestratorError("--run-dir must be an absolute path")
        requested_run_dir = requested_run_dir.resolve(strict=False)
        if args.run_name and args.run_name != requested_run_dir.name:
            raise OrchestratorError(
                "--run-name must match the basename of --run-dir")
        if not args.run_name:
            args.run_name = requested_run_dir.name
    elif not args.run_name:
        raise OrchestratorError("--run-name or --run-dir is required")
    if args.rna3_cell_reads_backfill_from_bam:
        print(
            "WARNING: --rna3-cell-reads-backfill-from-bam is deprecated; "
            "running --rna3-bam-evidence-from-bam with the compiled profiler",
            file=sys.stderr,
        )
        args.rna3_bam_evidence_from_bam = True
    args.run_name = safe_run_name(args.run_name)
    run_dir = (
        requested_run_dir
        if requested_run_dir is not None
        else Path(args.staging_root).expanduser().resolve(strict=False)
        / args.run_name
    )
    args.run_dir = str(run_dir)
    figure_defaults = (
        ("--rna3-figure-root", "rna3_figure_root",
         DEFAULT_RNA3_FIGURE_BASE),
        ("--rna5-figure-root", "rna5_figure_root",
         DEFAULT_RNA5_FIGURE_BASE),
        ("--atac-figure-root", "atac_figure_root",
         DEFAULT_ATAC_FIGURE_BASE),
    )
    for option, attribute, default_base in figure_defaults:
        supplied = getattr(args, attribute)
        value = Path(
            supplied or os.path.join(default_base, args.run_name, "Mapping")
        ).expanduser()
        if not value.is_absolute():
            raise OrchestratorError(f"{option} must be an absolute path")
        setattr(args, attribute, str(value.resolve(strict=False)))
    stages = parse_stages(args.stages)
    inputs = Inputs(
        rna3=resolve_runs(args.rna3_runs, args.rna3_raw_root),
        atac=resolve_runs(args.atac_runs, args.atac_raw_root),
        rna5=resolve_runs(args.rna5_runs, args.rna5_raw_root),
    )
    if not (inputs.rna3 or inputs.atac or inputs.rna5):
        raise OrchestratorError("select at least one of --rna3-runs, --atac-runs, --rna5-runs")
    if args.libraries and any(value < 1 for value in args.libraries):
        raise OrchestratorError("--libraries values must be positive integers")
    if args.libraries:
        args.libraries = sorted(set(args.libraries))
    if args.resubmit_jobs and not (args.resume and args.submit):
        raise OrchestratorError("--resubmit-jobs requires both --resume and --submit")
    if not 0 <= args.min_pe150_tso_match_fraction <= 1:
        raise OrchestratorError("--min-pe150-tso-match-fraction must be between 0 and 1")
    if args.max_cores is not None:
        if args.max_cores < 3:
            raise OrchestratorError("--max-cores must be at least 3")
        if not (inputs.rna3 or inputs.rna5):
            raise OrchestratorError("--max-cores currently requires an RNA input")
        if inputs.atac:
            raise OrchestratorError(
                "--max-cores currently guarantees RNA-only scheduling; run ATAC "
                "separately rather than claiming an unsafe combined ceiling"
            )
        if "map" in stages and args.rna_threads + 2 > args.max_cores:
            raise OrchestratorError(
                "--max-cores must accommodate one RNA mapping worker plus the "
                f"2-core Nextflow controller ({args.rna_threads + 2} cores required)"
            )
    if args.nodelist is not None:
        args.nodelist = safe_nodelist(args.nodelist)
    if args.array_max_concurrent is not None and args.array_max_concurrent < 1:
        raise OrchestratorError("--array-max-concurrent must be at least 1")
    if args.rna3_bam_evidence_cpus < 2:
        raise OrchestratorError("--rna3-bam-evidence-cpus must be at least 2")
    if args.rna3_bam_evidence_memory_gb < 4:
        raise OrchestratorError("--rna3-bam-evidence-memory-gb must be at least 4")
    if args.rna3_bam_evidence_max_concurrent < 1:
        raise OrchestratorError(
            "--rna3-bam-evidence-max-concurrent must be at least 1"
        )
    if args.rna3_bam_evidence_max_concurrent > 3:
        raise OrchestratorError(
            "--rna3-bam-evidence-max-concurrent may not exceed the audited ceiling of 3"
        )
    if args.rna3_bam_evidence_hash_bins < 1:
        raise OrchestratorError("--rna3-bam-evidence-hash-bins must be positive")
    if args.rna3_bam_evidence_from_bam:
        if not inputs.rna3:
            raise OrchestratorError(
                "--rna3-bam-evidence-from-bam requires a 3' RNA input"
            )
        if bool(args.rna3_bam_evidence_class_manifest) == bool(
            args.rna3_bam_evidence_no_biological_classification
        ):
            raise OrchestratorError(
                "RNA BAM evidence requires exactly one of "
                "--rna3-bam-evidence-class-manifest or "
                "--rna3-bam-evidence-no-biological-classification"
            )
        if args.rna3_bam_evidence_baseline_root:
            args.rna3_bam_evidence_baseline_root = str(
                Path(args.rna3_bam_evidence_baseline_root)
                .expanduser()
                .resolve(strict=False)
            )
        if args.rna3_bam_evidence_class_manifest:
            args.rna3_bam_evidence_class_manifest = str(
                Path(args.rna3_bam_evidence_class_manifest)
                .expanduser()
                .resolve(strict=False)
            )
    if args.rna3_mapping_baseline_stats:
        if not inputs.rna3:
            raise OrchestratorError(
                "--rna3-mapping-baseline-stats requires a 3' RNA input"
            )
        if "plot" not in stages:
            raise OrchestratorError(
                "--rna3-mapping-baseline-stats requires the plot stage"
            )
        args.rna3_mapping_baseline_stats = str(
            Path(args.rna3_mapping_baseline_stats)
            .expanduser()
            .resolve(strict=False)
        )
    if args.rna3_mapping_baseline_root:
        if not inputs.rna3:
            raise OrchestratorError(
                "--rna3-mapping-baseline-root requires a 3' RNA input"
            )
        if "plot" not in stages:
            raise OrchestratorError(
                "--rna3-mapping-baseline-root requires the plot stage"
            )
        args.rna3_mapping_baseline_root = str(
            Path(args.rna3_mapping_baseline_root)
            .expanduser()
            .resolve(strict=False)
        )
    if args.rna3_mapping_baseline_work_root:
        if not args.rna3_mapping_baseline_root:
            raise OrchestratorError(
                "--rna3-mapping-baseline-work-root requires "
                "--rna3-mapping-baseline-root"
            )
        args.rna3_mapping_baseline_work_root = str(
            Path(args.rna3_mapping_baseline_work_root)
            .expanduser()
            .resolve(strict=False)
        )
    if not args.rna3_mapping_baseline_label.strip():
        raise OrchestratorError("--rna3-mapping-baseline-label may not be empty")
    if not args.rna3_mapping_current_label.strip():
        raise OrchestratorError("--rna3-mapping-current-label may not be empty")

    resources = resolve_resources(args)
    discovery, failures, warnings = discover(
        inputs,
        set(args.libraries) if args.libraries else None,
        {
            "rna3": args.rna3_lib_prefix,
            "rna5": args.rna5_lib_prefix,
            "atac": args.atac_lib_prefix,
        },
    )
    failures.extend(protected_output_failures(run_dir))
    failures.extend(resource_failures(inputs, stages, resources, args))

    unsupported_5p = [
        f"{Path(run).name} ({fmt})"
        for run, fmt in discovery.rna5_formats.items()
        if fmt not in {"long-r2", "pe150"}
    ]
    if unsupported_5p and ({"trim", "map"} & stages):
        failures.append(
            "5' run contains an unknown or mixed read geometry: "
            + ", ".join(unsupported_5p)
            + ". Each input folder must contain exactly one recognized geometry."
        )
    if inputs.rna5 and args.rna5_map_mode == "together":
        long_libs = set(discovery.rna5_libraries_by_format.get("long-r2", []))
        pe_libs = set(discovery.rna5_libraries_by_format.get("pe150", []))
        overlap = sorted(long_libs & pe_libs)
        if overlap:
            failures.append(
                "--rna5-map-mode together is unsafe for libraries present in both formats: "
                + ", ".join(overlap)
                + ". STARsolo cannot re-deduplicate molecules by concatenating its BAMs; "
                "use the default separate mode until a validated UMI-level recount is available."
            )

    print_diagnosis(
        run_dir, stages, inputs, resources, discovery, failures, warnings
    )
    if args.diagnose_only:
        return 1 if failures else 0
    if failures:
        raise OrchestratorError(
            f"preflight found {len(failures)} blocking problem(s); no jobs were generated"
        )

    # Keep this handle live through the final sbatch/plan write. The advisory
    # lock prevents concurrent phase invocations from interleaving plans and
    # releasing readers against different continuation records.
    _run_lock = acquire_run_lock(run_dir)
    payload = config_payload(args, inputs, resources, discovery)
    prepare_run_directory(run_dir, args.resume, payload, args.resubmit_jobs)

    jobs: list[JobSpec] = []
    rna3_groups: list[tuple[str, Sequence[Path], Path]] = []
    rna5_groups: list[tuple[str, Sequence[Path], Path]] = []
    rna3_runs = discovery.active_runs.get("rna3", [])
    rna5_runs = discovery.active_runs.get("rna5", [])
    atac_runs = discovery.active_runs.get("atac", [])
    if rna3_runs:
        rna3_root = run_dir / "rna3"
        rna3_groups.append(("long-r2", rna3_runs, rna3_root))
        jobs.extend(
            generate_rna_jobs(
                "rna3",
                "3prime",
                "long-r2",
                rna3_runs,
                args.rna3_whitelist,
                rna3_root,
                run_dir,
                stages,
                args,
                resources,
            )
        )
    if rna5_runs:
        for read_format in ("long-r2", "pe150"):
            format_runs = [
                run
                for run in rna5_runs
                if discovery.rna5_formats.get(str(run)) == read_format
            ]
            if not format_runs:
                continue
            label = f"rna5_{read_format.replace('-', '_')}"
            root = run_dir / "rna5" / "formats" / read_format
            rna5_groups.append((read_format, format_runs, root))
            jobs.extend(
                generate_rna_jobs(
                    label,
                    "5prime",
                    read_format,
                    format_runs,
                    args.rna5_whitelist,
                    root,
                    run_dir,
                    stages,
                    args,
                    resources,
                )
            )
    if atac_runs and "map" in stages:
        jobs.append(
            generate_atac_mapping_job(
                atac_runs, run_dir / "atac", run_dir, args, resources
            )
        )

    labels = {job.label for job in jobs}
    if inputs.rna5 and args.rna5_map_mode == "together" and "map" in stages:
        map_dependencies = sorted(
            label
            for label in labels
            if label.startswith("rna5_") and label.endswith("_map")
        )
        jobs.append(
            generate_rna5_together_view_job(
                rna5_groups, discovery, run_dir, map_dependencies
            )
        )

    if args.rna3_bam_evidence_from_bam:
        current_labels = {job.label for job in jobs}
        jobs.extend(
            generate_rna_bam_evidence_jobs(
                discovery.libraries.get("rna3", []),
                run_dir / "rna3",
                run_dir,
                args,
                resources,
                ["rna3_map"] if "rna3_map" in current_labels else [],
                [path.name for path in inputs.rna3],
            )
        )

    if "qc" in stages and not args.no_trim_info:
        if rna3_groups:
            rna3_dependencies = sorted(
                job.label
                for job in jobs
                if (
                    job.label.startswith("rna3_trim_")
                    or job.label == "rna3_map"
                    or job.label == "rna3_bam_evidence_gather"
                )
            )
            jobs.extend(
                generate_barcode_trim_jobs(
                    "rna3",
                    rna3_groups,
                    discovery.libraries.get("rna3", []),
                    args.rna3_lib_prefix,
                    args.rna3_whitelist,
                    args.empty_drop_roster,
                    run_dir,
                    resources,
                    rna3_dependencies,
                    use_evidence_bridge=args.rna3_bam_evidence_from_bam,
                )
            )
        if rna5_groups:
            rna5_dependencies = sorted(
                job.label
                for job in jobs
                if (
                    job.label.startswith("rna5_")
                    and ("_trim_" in job.label or job.label.endswith("_map"))
                )
            )
            jobs.extend(
                generate_barcode_trim_jobs(
                    "rna5",
                    rna5_groups,
                    discovery.libraries.get("rna5", []),
                    args.rna5_lib_prefix,
                    args.rna5_whitelist,
                    args.empty_drop_roster,
                    run_dir,
                    resources,
                    rna5_dependencies,
                )
            )

    labels = {job.label for job in jobs}
    if "plot" in stages:
        if inputs.rna3:
            trim_dependencies = sorted(
                label for label in labels if label.startswith("rna3_trim_")
            )
            jobs.append(
                generate_trim_plot_job(
                    "rna3",
                    run_dir / "rna3",
                    run_dir,
                    Path(args.rna3_figure_root),
                    resources,
                    trim_dependencies,
                )
            )
            if "rna3_bam_evidence_gather" in labels:
                map_dependencies = ["rna3_bam_evidence_gather"]
            else:
                map_dependencies = ["rna3_map"] if "rna3_map" in labels else []
            jobs.append(
                generate_mapping_plot_job(
                    "rna3",
                    run_dir / "rna3",
                    run_dir,
                    Path(args.rna3_figure_root),
                    resources,
                    map_dependencies,
                    baseline_stats=args.rna3_mapping_baseline_stats,
                    baseline_label=args.rna3_mapping_baseline_label,
                    current_label=args.rna3_mapping_current_label,
                    baseline_base_path=args.rna3_mapping_baseline_root,
                    baseline_cell_reads_work_root=args.rna3_mapping_baseline_work_root,
                )
            )
        for read_format, _runs, root in rna5_groups:
            modality_label = f"rna5_{read_format.replace('-', '_')}"
            trim_dependencies = sorted(
                label for label in labels if label.startswith(f"{modality_label}_trim_")
            )
            jobs.append(
                generate_trim_plot_job(
                    modality_label, root, run_dir,
                    Path(args.rna5_figure_root) / "formats" / read_format,
                    resources, trim_dependencies
                )
            )
            if args.rna5_map_mode == "separate":
                map_label = f"{modality_label}_map"
                jobs.append(
                    generate_mapping_plot_job(
                        modality_label,
                        root,
                        run_dir,
                        Path(args.rna5_figure_root) / "formats" / read_format,
                        resources,
                        [map_label] if map_label in labels else [],
                    )
                )
        if rna5_groups and args.rna5_map_mode == "together":
            jobs.append(
                generate_mapping_plot_job(
                    "rna5",
                    run_dir / "rna5",
                    run_dir,
                    Path(args.rna5_figure_root) / "together",
                    resources,
                    ["rna5_together_view"] if "rna5_together_view" in labels else [],
                )
            )

    if inputs.atac and "qc" in stages:
        if args.rna_barcode_base:
            rna_barcode_base = Path(args.rna_barcode_base).resolve(strict=False)
        elif inputs.rna3:
            rna_barcode_base = run_dir / "rna3" / "mapping_output"
        else:
            rna_barcode_base = Path(DEFAULT_PRODUCTION_RNA)
        qc_dependencies = []
        if "atac_map" in {job.label for job in jobs}:
            qc_dependencies.append("atac_map")
        if rna_barcode_base == run_dir / "rna3" / "mapping_output" and "rna3_map" in {
            job.label for job in jobs
        }:
            qc_dependencies.append("rna3_map")
        jobs.append(
            generate_atac_qc_job(
                discovery.library_numbers.get("atac", []),
                run_dir / "atac" / "mapping_output",
                rna_barcode_base,
                run_dir,
                args,
                resources,
                qc_dependencies,
            )
        )
    if inputs.atac and "plot" in stages:
        plot_dependencies = ["atac_qc"] if "atac_qc" in {job.label for job in jobs} else []
        jobs.append(generate_atac_plot_job(
            run_dir, Path(args.atac_figure_root), resources,
            plot_dependencies))

    if "validate" in stages:
        current_labels = [job.label for job in jobs]
        jobs.append(
            generate_validation_job(
                inputs,
                discovery,
                run_dir,
                stages,
                current_labels,
                args.rna5_map_mode,
                bool("qc" in stages and not args.no_trim_info),
                bool(args.empty_drop_roster),
                args.rna3_bam_evidence_from_bam,
            )
        )

    if not jobs:
        raise OrchestratorError("the selected inputs/stages produced no jobs")

    for job in jobs:
        apply_slurm_controls(
            job.script,
            args.nodelist,
            args.array_max_concurrent,
        )

    if args.max_cores is not None:
        enforce_core_ceiling(jobs, args.max_cores)

    plan_path = run_dir / "control" / "job_plan.json"
    if args.resume:
        restore_prior_job_ids(
            plan_path,
            jobs,
            args.resubmit_jobs,
            resources.star_diagnostic_promoter,
        )
    atomic_json(
        plan_path,
        plan_payload(
            run_dir, stages, jobs, submitted=any(job.job_id for job in jobs)
        ),
    )

    if args.submit:
        submitted: dict[str, str] = {}
        dependency_ids: dict[str, str | None] = {}
        for job in jobs:
            missing = [label for label in job.dependencies if label not in submitted]
            if missing:
                raise OrchestratorError(
                    f"internal dependency order error for {job.label}: {', '.join(missing)}"
                )
            if job.job_id:
                submitted[job.label] = job.job_id
                print(f"Reusing prior submission: {job.label} -> job {job.job_id}")
                dependency_ids[job.label] = reused_job_dependency(
                    job.job_id, job.label
                )
                continue
            job.job_id = submit_job(
                job.script,
                [
                    dependency_id
                    for label in job.dependencies
                    if (dependency_id := dependency_ids[label]) is not None
                ],
            )
            submitted[job.label] = job.job_id
            dependency_ids[job.label] = job.job_id
            atomic_json(
                plan_path,
                plan_payload(run_dir, stages, jobs, submitted=True),
            )

    print_plan(
        run_dir,
        jobs,
        args.submit,
        args.max_cores,
        args.nodelist,
        args.array_max_concurrent,
    )
    return 0


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    try:
        if args._gather_trim_json:
            if not args._trim_root or not args._report_dir:
                raise OrchestratorError("internal gather worker requires trim/report paths")
            gather_trim_json(Path(args._trim_root), Path(args._report_dir))
            return 0
        if args._run_atac_qc:
            required = (args._collector, args._library, args._ramdisk, args._stats_dir)
            if any(value is None for value in required):
                raise OrchestratorError("internal ATAC QC worker arguments are incomplete")
            run_atac_qc_worker(args)
            return 0
        return orchestrate(args)
    except OrchestratorError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
