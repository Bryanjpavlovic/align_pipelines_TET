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
* validate expected files without publishing over production data.

Dry-run script generation is the default.  The same command with ``--submit``
launches the jobs.  Use ``--resume`` only to continue the same immutable run.

Repository and package layout
-----------------------------
Keep this file and the seven runtime helpers listed below together in
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


RELEASE = "2026-08-30-v14-align-repo-migration"

DEFAULT_RNA3_RAW_ROOT = "/mnt/beegfs/reads/3P_Multiome_10XRNA"
DEFAULT_ATAC_RAW_ROOT = "/mnt/beegfs/reads/10X_ATAC_multiome"
DEFAULT_RNA5_RAW_ROOT = "/mnt/beegfs/reads/5P_10XRNA"
DEFAULT_STAGING_ROOT = "/mnt/beegfs/tet2025_mapping_staging"
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
DEFAULT_PRODUCTION_RNA = "/mnt/beegfs/tetmultiome_rna_mapped/mapping_output"
DEFAULT_PRODUCTION_ATAC = "/mnt/beegfs/tetmultiome_atac/mapping_output"

ALL_STAGES = ("trim", "map", "qc", "plot", "validate")


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
    }


def payload_hash(payload: dict[str, object]) -> str:
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def scientific_config(payload: dict[str, object]) -> dict[str, object]:
    """Return the immutable data/reference/analysis portion of a run config."""
    return {
        key: value
        for key, value in payload.items()
        if key not in {"release", "resources", "scheduling"}
    }


def runtime_config(payload: dict[str, object]) -> dict[str, object]:
    """Return warning-only software provenance and scheduling controls."""
    return {
        "release": payload.get("release"),
        "resources": payload.get("resources"),
        "scheduling": payload.get("scheduling"),
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
                            str(bam),
                            str(rg_metadata),
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
                                str(bam),
                                str(rg_metadata),
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
        "bam\trg_metadata\n"
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
    resources: Resources,
    dependencies: Sequence[str],
) -> JobSpec:
    assert resources.trim_plotter is not None
    script = run_dir / "control" / "slurm" / f"plot_{modality}_trimming.sbatch"
    work = run_dir / "qc" / modality / "trimming"
    report_dir = work / "json_reports"
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
    resources: Resources,
    dependencies: Sequence[str],
) -> JobSpec:
    assert resources.mapping_plotter is not None
    script = run_dir / "control" / "slurm" / f"plot_{modality}_mapping.sbatch"
    out = run_dir / "qc" / modality / "mapping"
    marker = out / "MAPPING_PLOTS_COMPLETE.ok"
    body = f"""if [[ -s {q(marker)} ]]; then
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
    --plot-type all
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
    resources: Resources,
    dependencies: Sequence[str],
) -> JobSpec:
    assert resources.atac_plotter is not None
    stats = run_dir / "qc" / "atac" / "stats"
    out = run_dir / "qc" / "atac" / "plots"
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


def validation_checks(
    inputs: Inputs,
    discovery: Discovery,
    run_dir: Path,
    stages: set[str],
    rna5_map_mode: str,
    include_trim_barcode_qc: bool,
    include_empty_drop_qc: bool,
) -> list[tuple[str, str, str]]:
    checks: list[tuple[str, str, str]] = []
    if inputs.rna3:
        for library in discovery.libraries.get("rna3", []):
            base = run_dir / "rna3" / "mapping_output" / library
            checks.extend(
                [
                    ("rna3", library, str(base / "gex.bam")),
                    ("rna3", library, str(base / "gex.bam.bai")),
                    ("rna3", library, str(base / "Summary.csv")),
                    ("rna3", library, str(base / "raw" / "barcodes.tsv.gz")),
                    ("rna3", library, str(base / "filtered" / "barcodes.tsv.gz")),
                ]
            )
        checks.extend(
            [
                ("rna3_provenance", "metadata", str(run_dir / "rna3" / "mapping_project" / "rg_metadata.tsv")),
                ("rna3_provenance", "manifest", str(run_dir / "rna3" / "mapping_project" / "symlink_manifest.tsv")),
            ]
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
                    [
                        (f"rna5_{read_format}", library, str(base / "gex.bam")),
                        (f"rna5_{read_format}", library, str(base / "gex.bam.bai")),
                        (f"rna5_{read_format}", library, str(base / "Summary.csv")),
                        (f"rna5_{read_format}", library, str(base / "raw" / "barcodes.tsv.gz")),
                        (f"rna5_{read_format}", library, str(base / "filtered" / "barcodes.tsv.gz")),
                    ]
                )
            format_project = run_dir / "rna5" / "formats" / read_format / "mapping_project"
            if discovery.rna5_libraries_by_format.get(read_format):
                checks.extend(
                    [
                        (f"rna5_{read_format}_provenance", "metadata", str(format_project / "rg_metadata.tsv")),
                        (f"rna5_{read_format}_provenance", "manifest", str(format_project / "symlink_manifest.tsv")),
                    ]
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
                    ("atac", library, str(base / "atac_fragments.tsv.gz")),
                    ("atac", library, str(base / "atac_fragments.tsv.gz.tbi")),
                ]
            )
        checks.extend(
            [
                ("atac_provenance", "metadata", str(run_dir / "atac" / "mapping_project" / "rg_metadata.tsv")),
                ("atac_provenance", "manifest", str(run_dir / "atac" / "mapping_project" / "symlink_manifest.tsv")),
            ]
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
) -> JobSpec:
    checks = validation_checks(
        inputs,
        discovery,
        run_dir,
        stages,
        rna5_map_mode,
        include_trim_barcode_qc,
        include_empty_drop_qc,
    )
    check_file = run_dir / "control" / "expected_outputs.tsv"
    atomic_text(
        check_file,
        "modality\tlibrary\tpath\n"
        + "".join(f"{modality}\t{library}\t{path}\n" for modality, library, path in checks),
    )
    summary = run_dir / "validation" / "output_validation.tsv"
    marker = run_dir / "validation" / "RUN_COMPLETE.ok"
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
echo "Validation passed: {marker}"
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
        print(
            f"WARNING: unrecognized SLURM state for reused job {label} "
            f"({job_id}): {','.join(sorted(set(states)))}; treating its "
            "dependency as already satisfied"
        )
    else:
        print(
            f"WARNING: SLURM no longer retains accounting state for reused "
            f"job {label} ({job_id}); treating its dependency as already "
            "satisfied"
        )
    return None


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


def restore_prior_job_ids(
    plan_path: Path,
    jobs: Sequence[JobSpec],
    requested_resubmits: Sequence[str] | None,
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

    # A retried upstream job invalidates every prior downstream submission in
    # this plan; walk the dependency graph until the descendant set stabilizes.
    retry = set(requested)
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
    parser.add_argument("--run-name", help="Unique staged run name (required for normal operation)")
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
    if not args.run_name:
        raise OrchestratorError("--run-name is required")
    args.run_name = safe_run_name(args.run_name)
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

    run_dir = Path(args.staging_root).expanduser().resolve(strict=False) / args.run_name
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

    if "qc" in stages and not args.no_trim_info:
        if rna3_groups:
            rna3_dependencies = sorted(
                job.label
                for job in jobs
                if job.label.startswith("rna3_trim_") or job.label == "rna3_map"
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
                    resources,
                    trim_dependencies,
                )
            )
            map_dependencies = ["rna3_map"] if "rna3_map" in labels else []
            jobs.append(
                generate_mapping_plot_job(
                    "rna3",
                    run_dir / "rna3",
                    run_dir,
                    resources,
                    map_dependencies,
                )
            )
        for read_format, _runs, root in rna5_groups:
            modality_label = f"rna5_{read_format.replace('-', '_')}"
            trim_dependencies = sorted(
                label for label in labels if label.startswith(f"{modality_label}_trim_")
            )
            jobs.append(
                generate_trim_plot_job(
                    modality_label, root, run_dir, resources, trim_dependencies
                )
            )
            if args.rna5_map_mode == "separate":
                map_label = f"{modality_label}_map"
                jobs.append(
                    generate_mapping_plot_job(
                        modality_label,
                        root,
                        run_dir,
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
        jobs.append(generate_atac_plot_job(run_dir, resources, plot_dependencies))

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
        restore_prior_job_ids(plan_path, jobs, args.resubmit_jobs)
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
