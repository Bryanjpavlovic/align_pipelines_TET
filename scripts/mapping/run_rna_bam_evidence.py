#!/usr/bin/env python3
"""Run, validate, and gather one-pass RNA BAM evidence.

Python never decodes BAM records. The compiled profiler owns the sole BAM
traversal. Large profiler products are validated and reconciled as sorted
streams; only filtered-cell vectors and small RG/library summaries are retained.
"""

from __future__ import annotations

import argparse
import csv
import fcntl
import gzip
import hashlib
import io
import json
import os
import resource
import shutil
import subprocess
import sys
import tempfile
import time
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Iterator, TextIO


RELEASE = "2026-09-05-rna-bam-evidence-v7-direct-full-run"
HASH_ALGORITHM = "fnv1a64_seeded_v1"
DEFAULT_HASH_SEED = 1469598103934665603
STARSOLO_FEATURE = "GeneFull_Ex50pAS"
STARSOLO_UMI_FILTERING = "MultiGeneUMI_CR"
STARSOLO_UMI_DEDUP = "1MM_CR"
STARSOLO_MULTIMAPPERS = "EM"
SUMMARY_UNIQUE_READ_METRIC = (
    f"Unique Reads in Cells Mapped to {STARSOLO_FEATURE}"
)
ORDINARY_COUNTEDU_READ_DEFINITION = (
    "distinct_RG_QNAME_declared_RG_valid_CB_singleton_feature_GX_"
    "NH_unrestricted_UB_not_required_STARsolo_pre_UMI_filter_countedU"
)
ORDINARY_MOLECULE_DEFINITION = (
    "distinct_CB_GX_valid_STARsolo_corrected_UB_after_1MM_CR_and_"
    "MultiGeneUMI_CR_NH_unrestricted"
)
MULTIMAPPER_DEFINITION = (
    "STARsolo_multi_gene_EM_unavailable_from_standard_uppercase_GX_UB_"
    "BAM_tags_NH_is_not_EM_membership"
)
NH_GT1_UNIQUE_GENE_COUNTEDU_DEFINITION = (
    "subset_of_ordinary_countedU_reads_with_NH_gt1_and_singleton_feature_GX"
)
STARSOLO_EM_EVIDENCE_AVAILABILITY = (
    "unavailable_from_standard_uppercase_GX_UB_BAM_tags"
)
MIN_TWO_READER_OVERLAP_FRACTION = 0.5
GATHER_GZIP_LEVEL = 1

REQUIRED_MANIFEST_COLUMNS = {
    "library",
    "library_dir",
    "bam",
    "bam_index",
    "summary",
    "raw_barcodes",
    "filtered_barcodes",
    "raw_features",
    "raw_matrix",
    "filtered_matrix",
    "rg_metadata",
    "source_order",
    "scientific_config",
    "starsolo_feature",
    "starsolo_umi_filtering",
    "starsolo_umi_dedup",
    "starsolo_multimappers",
    "biological_classification_intent",
    "output_dir",
}
SIGNATURE_INPUT_COLUMNS = (
    "bam", "bam_index", "summary", "raw_barcodes", "filtered_barcodes",
    "raw_features", "raw_matrix", "filtered_matrix", "rg_metadata",
    "source_order", "old_raw_barcodes", "old_filtered_barcodes",
    "class_manifest", "native_cell_reads",
)
COMPRESSED_PRODUCTS = {
    "barcode_read_metrics.tsv": "barcode_read_metrics.tsv.gz",
    "barcode_rg_metrics.tsv": "barcode_rg_metrics.tsv.gz",
    "molecule_source_hash_bins.tsv": "molecule_source_hash_bins.tsv.gz",
    "rg_contig_class_summary.tsv": "rg_contig_class_summary.tsv.gz",
    "raw_to_corrected_barcode_counts.tsv": "raw_to_corrected_barcode_counts.tsv.gz",
}
PLAIN_PRODUCTS = {"rg_summary.tsv": "rg_summary.tsv"}
REQUIRED_PRODUCTS = tuple(COMPRESSED_PRODUCTS.values()) + tuple(PLAIN_PRODUCTS.values())

# Candidate molecules are not additive across RGs because one molecule may
# occur in several RGs. All fields below are record counts and are additive.
CORE_ADDITIVE_FIELDS = (
    "all_records",
    "primary_records",
    "primary_mapped_reads",
    "secondary_records",
    "supplementary_records",
    "qcfail_records",
    "bam_duplicate_flag_reads",
    "nh1_reads",
    "nh_gt1_reads",
    "unambiguous_gx_reads",
    "unique_gene_tagged_reads",
    "candidate_countedU_reads",
    "nh_gt1_unique_gene_countedU_reads",
    "current_cell_associated_reads",
    "raw_nonfiltered_droplet_reads",
    "biological_classified_reads",
    "classification_by_contig_reads",
    "classification_by_feature_reads",
    "mitochondrial_reads",
    "rrna_reads",
)
OPTIONAL_CLASS_COUNT_FIELDS = {"mitochondrial_reads", "rrna_reads"}

SCHEMAS = {
    "barcode_read_metrics.tsv.gz": {
        "CB",
        "current_raw_member",
        "current_filtered_member",
        "historical_raw_member",
        "historical_filtered_member",
        "barcode_category",
        "primary_mapped_reads",
        "nh1_reads",
        "nh_gt1_reads",
        "candidate_countedU_reads",
        "nh_gt1_unique_gene_countedU_reads",
        "biological_classified_reads",
        "mitochondrial_reads",
        "rrna_reads",
        "biological_classification_status",
        "candidate_matrix_molecules",
    },
    "barcode_rg_metrics.tsv.gz": {
        "CB",
        "RG",
        "source_id",
        "source_index",
        "primary_mapped_reads",
        "candidate_countedU_reads",
        "nh_gt1_unique_gene_countedU_reads",
        "biological_classification_status",
    },
    "molecule_source_hash_bins.tsv.gz": {
        "CB",
        "source_mask_hex",
        "nested_min_hash_bin",
        "n_hash_bins",
        "candidate_matrix_molecules",
    },
    "rg_summary.tsv": {
        "library",
        "RG",
        "source_id",
        "source_index",
        "all_records",
        "primary_mapped_reads",
        "candidate_matrix_molecules",
    },
    "rg_contig_class_summary.tsv.gz": {
        "library",
        "RG",
        "source_index",
        "contig_index",
        "contig",
        "contig_class_available",
        "species",
        "mitochondrial",
        "rrna",
        "reference_class",
        "biological_classification_status",
    },
    "raw_to_corrected_barcode_counts.tsv.gz": {
        "CR",
        "RG",
        "CB",
        "read_count",
        "within_rg_conflict",
        "cross_rg_conflict",
        "global_conflict",
    },
}
UNIQUE_KEYS = {
    "barcode_read_metrics.tsv.gz": ("CB",),
    "barcode_rg_metrics.tsv.gz": ("CB", "RG"),
    "molecule_source_hash_bins.tsv.gz": (
        "CB",
        "source_mask_hex",
        "nested_min_hash_bin",
    ),
    "rg_summary.tsv": ("library", "RG"),
    "rg_contig_class_summary.tsv.gz": ("library", "RG", "contig"),
    "raw_to_corrected_barcode_counts.tsv.gz": ("CR", "RG", "CB"),
}

GATHER_SCHEMAS = {
    "bam_evidence_inventory.tsv": {
        "library", "product", "path", "rows", "bytes", "sha256"
    },
    "all_libraries_rg_summary.tsv.gz": {
        "library", "RG", "source_id", "source_index", "all_records"
    },
    "all_libraries_source_yield.tsv": {
        "library", "source_id", "source_order_index", "primary_mapped_reads",
        "candidate_matrix_molecules",
        "marginal_candidate_molecules_first_observed",
        "biological_classification_status",
    },
    "all_libraries_barcode_summary.tsv.gz": {
        "library", "CB", "candidate_countedU_reads",
        "nh_gt1_unique_gene_countedU_reads", "candidate_matrix_molecules",
        "biological_classification_status",
    },
}
GATHER_KEYS = {
    "bam_evidence_inventory.tsv": ("library", "product"),
    "all_libraries_rg_summary.tsv.gz": ("library", "source_index", "RG"),
    "all_libraries_source_yield.tsv": (
        "library", "source_order_index", "source_id"
    ),
    "all_libraries_barcode_summary.tsv.gz": ("library", "CB"),
}


class EvidenceError(RuntimeError):
    """A user-facing evidence, validation, or orchestration error."""


class ReconciliationFailure(EvidenceError):
    def __init__(self, details: dict[str, object]):
        super().__init__(
            "semantic reconciliation failed: "
            + ", ".join(details["failed_equalities"])
        )
        self.details = details


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def slurm_execution_identity() -> dict[str, object]:
    """Return the exact Slurm allocation/task identity for this invocation."""
    job_id = os.environ.get("SLURM_JOB_ID", "")
    array_job_id = os.environ.get("SLURM_ARRAY_JOB_ID", "")
    array_task_id = os.environ.get("SLURM_ARRAY_TASK_ID", "")
    if not (job_id or array_job_id or array_task_id):
        return {"kind": "not_slurm"}
    if not job_id.isdecimal():
        raise EvidenceError("SLURM_JOB_ID is missing or non-numeric")
    if bool(array_job_id) != bool(array_task_id):
        raise EvidenceError(
            "SLURM_ARRAY_JOB_ID and SLURM_ARRAY_TASK_ID must be present together"
        )
    if not array_job_id:
        return {
            "kind": "standalone",
            "SLURM_JOB_ID": job_id,
            "accounting_job_task_id": job_id,
        }
    if not array_job_id.isdecimal() or not array_task_id.isdecimal():
        raise EvidenceError("Slurm array job/task IDs must be numeric")
    return {
        "kind": "array_task",
        "SLURM_JOB_ID": job_id,
        "SLURM_ARRAY_JOB_ID": array_job_id,
        "SLURM_ARRAY_TASK_ID": int(array_task_id),
        "accounting_job_task_id": f"{array_job_id}_{array_task_id}",
    }


def atomic_text(path: Path, text: str, mode: int | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    temporary.write_text(text, encoding="utf-8")
    if mode is not None:
        temporary.chmod(mode)
    os.replace(temporary, path)


def atomic_json(path: Path, payload: object) -> None:
    atomic_text(path, json.dumps(payload, indent=2, sort_keys=True) + "\n")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


class DigestingRawReader(io.RawIOBase):
    """Hash bytes as a sequential binary consumer reads them."""

    def __init__(self, raw: io.RawIOBase, digest: object):
        super().__init__()
        self.raw = raw
        self.digest = digest

    def readable(self) -> bool:
        return True

    def readinto(self, buffer: object) -> int:
        count = self.raw.readinto(buffer)
        if count:
            self.digest.update(memoryview(buffer)[:count])
        return count

    def close(self) -> None:
        if not self.closed:
            self.raw.close()
        super().close()


@contextmanager
def open_hashed_text(path: Path, digest: object) -> Iterator[TextIO]:
    """Open text while hashing the exact on-disk bytes in the same pass."""
    raw = path.open("rb", buffering=0)
    digesting = DigestingRawReader(raw, digest)
    buffered = io.BufferedReader(digesting, buffer_size=1024 * 1024)
    binary: io.BufferedIOBase | gzip.GzipFile
    if path.suffix == ".gz":
        binary = gzip.GzipFile(fileobj=buffered, mode="rb")
    else:
        binary = buffered
    handle = io.TextIOWrapper(binary, encoding="utf-8", newline="")
    completed = False
    try:
        yield handle
        completed = True
    finally:
        handle.close()
        # GzipFile intentionally leaves a caller-owned file object open. Drain
        # any bytes it buffered past the gzip member so the digest covers the
        # complete file, including any invalid trailing bytes.
        if completed and path.suffix == ".gz" and not buffered.closed:
            for _block in iter(lambda: buffered.read(1024 * 1024), b""):
                pass
        if not buffered.closed:
            buffered.close()


@contextmanager
def open_text(path: Path) -> Iterator[TextIO]:
    if path.suffix == ".gz":
        with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
            yield handle
    else:
        with path.open("r", encoding="utf-8", newline="") as handle:
            yield handle


def nonempty_file(path: Path, label: str) -> Path:
    if not path.is_file() or path.stat().st_size == 0:
        raise EvidenceError(f"missing or empty {label}: {path}")
    return path


def read_manifest(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    nonempty_file(path, "evidence manifest")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        missing = sorted(REQUIRED_MANIFEST_COLUMNS - set(fields))
        if missing:
            raise EvidenceError(
                f"manifest lacks required columns: {', '.join(missing)}"
            )
        rows = list(reader)
    if not rows:
        raise EvidenceError(f"manifest has no library rows: {path}")
    libraries = [row["library"] for row in rows]
    if any(not value for value in libraries) or len(libraries) != len(set(libraries)):
        raise EvidenceError(f"manifest has blank or duplicate libraries: {path}")
    return fields, rows


def select_manifest_row(path: Path, row_index: int) -> dict[str, str]:
    _fields, rows = read_manifest(path)
    if row_index < 0 or row_index >= len(rows):
        raise EvidenceError(
            f"manifest row index {row_index} is outside 0..{len(rows) - 1}"
        )
    return rows[row_index]


def canonical_payload_hash(payload: dict[str, object]) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()


def load_project_scientific_config(
    manifest: Path,
    fields: list[str],
    rows: list[dict[str, str]],
) -> dict[str, object]:
    config_paths = {row.get("scientific_config", "") for row in rows}
    if len(config_paths) != 1 or "" in config_paths:
        raise EvidenceError(
            "manifest rows do not reference one frozen scientific configuration"
        )
    config_path = nonempty_file(
        Path(next(iter(config_paths))).expanduser().resolve(),
        "frozen scientific configuration",
    )
    try:
        config = json.loads(config_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise EvidenceError(
            f"could not read frozen scientific configuration: {config_path}"
        ) from exc
    if not isinstance(config, dict):
        raise EvidenceError("frozen scientific configuration is not a JSON object")
    stored_hash = config.get("configuration_hash")
    unhashed = dict(config)
    unhashed.pop("configuration_hash", None)
    if stored_hash != canonical_payload_hash(unhashed):
        raise EvidenceError("frozen scientific configuration hash is invalid")
    if config.get("manifest_fields") != fields or config.get("manifest_rows") != rows:
        raise EvidenceError(
            "current project manifest differs from the pilot scientific configuration"
        )
    source = config.get("source_order")
    if not isinstance(source, dict):
        raise EvidenceError("scientific configuration lacks source-order provenance")
    source_path = nonempty_file(
        Path(str(source.get("path", ""))), "frozen source-order file"
    )
    source_bytes = source_path.read_bytes()
    if (
        hashlib.sha256(source_bytes).hexdigest() != source.get("sha256")
        or source_bytes.decode("utf-8") != source.get("contents")
    ):
        raise EvidenceError("source-order file differs from the pilot configuration")
    classification = config.get("biological_classification")
    if not isinstance(classification, dict):
        raise EvidenceError("scientific configuration lacks classification status")
    if classification.get("status") == "manifest":
        class_path = nonempty_file(
            Path(str(classification.get("path", ""))),
            "frozen biological-classification manifest",
        )
        class_bytes = class_path.read_bytes()
        try:
            class_contents = class_bytes.decode("utf-8")
        except UnicodeDecodeError as exc:
            raise EvidenceError("classification manifest is not UTF-8 text") from exc
        if (
            hashlib.sha256(class_bytes).hexdigest() != classification.get("sha256")
            or class_contents != classification.get("contents")
        ):
            raise EvidenceError(
                "classification manifest differs from the pilot configuration"
            )
    elif classification != {
        "status": "explicit_unavailable",
        "path": "",
        "sha256": "",
        "contents": "",
    }:
        raise EvidenceError("invalid classification-unavailable configuration")
    frozen_profiler_path(config)
    return config


def scientific_parameters(
    row: dict[str, str], config: dict[str, object]
) -> dict[str, object]:
    hash_config = config.get("hash")
    starsolo = config.get("starsolo")
    if not isinstance(hash_config, dict) or not isinstance(starsolo, dict):
        raise EvidenceError("scientific configuration lacks hash/STARsolo semantics")
    expected = {
        "hash_algorithm": HASH_ALGORITHM,
        "hash_seed": hash_config.get("seed"),
        "hash_bins": hash_config.get("bins"),
        "starsolo_feature": row["starsolo_feature"],
        "starsolo_umi_filtering": row["starsolo_umi_filtering"],
        "starsolo_umi_dedup": row["starsolo_umi_dedup"],
        "starsolo_multimappers": row["starsolo_multimappers"],
        "ordinary_countedU_read_definition":
            ORDINARY_COUNTEDU_READ_DEFINITION,
        "ordinary_molecule_definition": ORDINARY_MOLECULE_DEFINITION,
        "multimapper_definition": MULTIMAPPER_DEFINITION,
        "summary_unique_read_metric": SUMMARY_UNIQUE_READ_METRIC,
        "nh_gt1_unique_gene_countedU_definition":
            NH_GT1_UNIQUE_GENE_COUNTEDU_DEFINITION,
        "starsolo_EM_evidence_availability":
            STARSOLO_EM_EVIDENCE_AVAILABILITY,
    }
    config_expected = {
        "hash_algorithm": hash_config.get("algorithm"),
        "hash_seed": hash_config.get("seed"),
        "hash_bins": hash_config.get("bins"),
        "starsolo_feature": starsolo.get("feature"),
        "starsolo_umi_filtering": starsolo.get("umi_filtering"),
        "starsolo_umi_dedup": starsolo.get("umi_dedup"),
        "starsolo_multimappers": starsolo.get("multimappers"),
        "ordinary_countedU_read_definition":
            starsolo.get("ordinary_countedU_read_definition"),
        "ordinary_molecule_definition": starsolo.get("ordinary_molecule_definition"),
        "multimapper_definition": starsolo.get("multimapper_definition"),
        "summary_unique_read_metric": starsolo.get("summary_unique_read_metric"),
        "nh_gt1_unique_gene_countedU_definition":
            starsolo.get("nh_gt1_unique_gene_countedU_definition"),
        "starsolo_EM_evidence_availability":
            starsolo.get("starsolo_EM_evidence_availability"),
    }
    if expected != config_expected:
        raise EvidenceError(
            "manifest/runtime STARsolo or hash semantics differ from the pilot configuration"
        )
    return expected


def roster(path: Path) -> list[str]:
    """Load a filtered-cell roster. Large raw rosters use roster_count()."""
    nonempty_file(path, "barcode roster")
    values: list[str] = []
    seen: set[str] = set()
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            value = line.rstrip("\r\n").split("\t", 1)[0]
            if not value:
                continue
            if value in seen:
                raise EvidenceError(
                    f"duplicate barcode {value!r} in {path}:{line_number}"
                )
            seen.add(value)
            values.append(value)
    if not values:
        raise EvidenceError(f"barcode roster is empty: {path}")
    return values


def roster_count(path: Path) -> int:
    nonempty_file(path, "barcode roster")
    count = 0
    with open_text(path) as handle:
        for line in handle:
            if line.rstrip("\r\n").split("\t", 1)[0]:
                count += 1
    if not count:
        raise EvidenceError(f"barcode roster is empty: {path}")
    return count


def source_order(path: Path) -> list[str]:
    nonempty_file(path, "source-order file")
    values: list[str] = []
    seen: set[str] = set()
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#"):
            continue
        value = line.split("\t", 1)[0]
        if value in {"source_id", "bp_id"}:
            continue
        if value in seen:
            raise EvidenceError(f"duplicate source {value!r} in {path}")
        seen.add(value)
        values.append(value)
    if not values or len(values) > 64:
        raise EvidenceError(f"source order must contain 1..64 sources: {path}")
    return values


def rg_manifest_rows(path: Path, library: str) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"library", "bp_id", "rg_id"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise EvidenceError(
                f"RG metadata lacks columns: {', '.join(sorted(missing))}"
            )
        rows = [dict(row) for row in reader if row["library"] == library]
    ids = [row["rg_id"] for row in rows]
    if not rows or any(not value for value in ids) or len(ids) != len(set(ids)):
        raise EvidenceError(
            f"RG metadata rows for {library} are empty, blank, or duplicated"
        )
    return rows


def matrix_shape_and_counts(
    matrix_path: Path,
    expected_barcodes: int,
    *,
    keep_per_barcode: bool,
) -> tuple[dict[str, int], list[int] | None, list[int] | None]:
    """Stream MatrixMarket entries; retain only fixed-size integer vectors."""
    nonempty_file(matrix_path, "STARsolo matrix")
    with open_text(matrix_path) as handle:
        banner = handle.readline().strip().lower()
        if not banner.startswith("%%matrixmarket matrix coordinate"):
            raise EvidenceError(f"unsupported MatrixMarket banner in {matrix_path}")
        line = handle.readline()
        while line and line.startswith("%"):
            line = handle.readline()
        try:
            n_features, n_barcodes, declared_nnz = map(int, line.split())
        except (AttributeError, TypeError, ValueError) as exc:
            raise EvidenceError(
                f"malformed MatrixMarket dimensions in {matrix_path}"
            ) from exc
        if n_barcodes != expected_barcodes:
            raise EvidenceError(
                f"matrix/barcode dimension mismatch for {matrix_path}: "
                f"{n_barcodes} columns versus {expected_barcodes} barcodes"
            )
        umi = [0] * n_barcodes if keep_per_barcode else None
        genes = [0] * n_barcodes if keep_per_barcode else None
        observed_nnz = 0
        molecule_total = 0
        for line_number, line in enumerate(handle, start=3):
            if not line.strip() or line.startswith("%"):
                continue
            fields = line.split()
            if len(fields) != 3:
                raise EvidenceError(
                    f"malformed matrix entry at {matrix_path}:{line_number}"
                )
            try:
                feature_index, barcode_index, value = map(int, fields)
            except ValueError as exc:
                raise EvidenceError(
                    f"non-integer matrix entry at {matrix_path}:{line_number}"
                ) from exc
            if (
                not 1 <= feature_index <= n_features
                or not 1 <= barcode_index <= n_barcodes
                or value <= 0
            ):
                raise EvidenceError(
                    f"invalid matrix entry at {matrix_path}:{line_number}"
                )
            observed_nnz += 1
            molecule_total += value
            if umi is not None and genes is not None:
                umi[barcode_index - 1] += value
                genes[barcode_index - 1] += 1
        if observed_nnz != declared_nnz:
            raise EvidenceError(
                f"MatrixMarket nnz mismatch in {matrix_path}: "
                f"{observed_nnz} entries versus declared {declared_nnz}"
            )
    return (
        {
            "features": n_features,
            "barcodes": n_barcodes,
            "nnz": declared_nnz,
            "molecules": molecule_total,
        },
        umi,
        genes,
    )


def summary_values(path: Path) -> dict[str, float]:
    nonempty_file(path, "STARsolo Summary.csv")
    result: dict[str, float] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.reader(handle):
            if len(row) < 2:
                continue
            try:
                result[row[0]] = float(row[1])
            except ValueError:
                continue
    return result


def iter_tsv(
    path: Path, required: Iterable[str] = ()
) -> Iterator[dict[str, str]]:
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise EvidenceError(f"TSV has no header: {path}")
        missing = set(required) - set(reader.fieldnames)
        if missing:
            raise EvidenceError(
                f"{path} lacks columns: {', '.join(sorted(missing))}"
            )
        for row in reader:
            yield row


def product_sort_key(name: str, row: dict[str, str]) -> tuple[object, ...]:
    if name == "molecule_source_hash_bins.tsv.gz":
        return (
            row["CB"],
            int(row["source_mask_hex"], 16),
            int(row["nested_min_hash_bin"]),
        )
    if name == "rg_summary.tsv":
        index = int(row["source_index"])
        return (
            row["library"],
            index if index >= 0 else 1 << 30,
            row["RG"],
        )
    if name == "all_libraries_rg_summary.tsv.gz":
        index = int(row["source_index"])
        return (
            row["library"],
            index if index >= 0 else 1 << 30,
            row["RG"],
        )
    if name == "barcode_rg_metrics.tsv.gz":
        return (row["CB"], int(row["source_index"]), row["RG"])
    if name == "rg_contig_class_summary.tsv.gz":
        return (
            row["library"],
            int(row["source_index"]),
            row["RG"],
            int(row["contig_index"]),
        )
    return tuple(row[field] for field in UNIQUE_KEYS[name])


def validate_product(path: Path, expected_name: str) -> dict[str, object]:
    """Validate, count, and hash in one bounded-memory sequential pass."""
    nonempty_file(path, expected_name)
    digest = hashlib.sha256()
    with open_hashed_text(path, digest) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        missing = sorted(SCHEMAS[expected_name] - set(fields))
        if missing:
            raise EvidenceError(
                f"{path} lacks required columns: {', '.join(missing)}"
            )
        previous: tuple[object, ...] | None = None
        rows = 0
        for line_number, row in enumerate(reader, start=2):
            raw_key = tuple(row.get(field, "") for field in UNIQUE_KEYS[expected_name])
            if not all(raw_key):
                raise EvidenceError(
                    f"blank uniqueness key in {path} line {line_number}"
                )
            try:
                key = product_sort_key(expected_name, row)
            except ValueError as exc:
                raise EvidenceError(
                    f"malformed sort key in {path} line {line_number}"
                ) from exc
            if previous is not None and key <= previous:
                kind = "duplicate" if key == previous else "out-of-order"
                raise EvidenceError(
                    f"{kind} key {raw_key!r} in {path} line {line_number}"
                )
            previous = key
            rows += 1
    return {
        "rows": rows,
        "bytes": path.stat().st_size,
        "sha256": digest.hexdigest(),
        "schema": fields,
    }


def validate_recorded_product_metadata(
    path: Path,
    recorded: object,
    required_schema: set[str],
    label: str,
    *,
    require_rows: bool = False,
) -> dict[str, object]:
    """Fast resume check for a product validated when it was published."""
    nonempty_file(path, label)
    if not isinstance(recorded, dict):
        raise EvidenceError(f"{label} has no recorded product audit")
    schema = recorded.get("schema")
    if not isinstance(schema, list) or not required_schema.issubset(schema):
        raise EvidenceError(f"{label} has an invalid recorded schema")
    rows = recorded.get("rows")
    byte_count = recorded.get("bytes")
    digest = recorded.get("sha256")
    if (
        not isinstance(rows, int)
        or rows < (1 if require_rows else 0)
        or not isinstance(byte_count, int)
        or byte_count <= 0
        or not isinstance(digest, str)
        or len(digest) != 64
    ):
        raise EvidenceError(f"{label} has malformed recorded product statistics")
    if path.stat().st_size != byte_count:
        raise EvidenceError(f"{label} no longer matches recorded audit field bytes")
    return dict(recorded)


def parse_profiler_summary(path: Path) -> dict[str, str]:
    nonempty_file(path, "profiler summary")
    result: dict[str, str] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        if next(reader, None) != ["metric", "value"]:
            raise EvidenceError(f"unexpected profiler summary header: {path}")
        for row in reader:
            if len(row) < 2 or not row[0] or row[0] in result:
                raise EvidenceError(f"malformed or duplicate profiler metric in {path}")
            result[row[0]] = row[1]
    return result


def profiler_version(executable: str) -> str:
    try:
        result = subprocess.run(
            [executable, "--version"], capture_output=True, text=True, check=False
        )
    except OSError as exc:
        raise EvidenceError(
            f"could not execute profiler {executable}: {exc}"
        ) from exc
    if result.returncode != 0 or not result.stdout.strip():
        raise EvidenceError(
            f"could not query profiler version from {executable}: "
            f"{result.stderr.strip()}"
        )
    return result.stdout.strip()


def file_signature(path: Path) -> dict[str, int | str]:
    stat = path.stat()
    return {
        "path": str(path.resolve()),
        "bytes": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
    }


def frozen_profiler_path(scientific_config: dict[str, object]) -> Path:
    frozen = scientific_config.get("profiler_executable")
    if not isinstance(frozen, dict):
        raise EvidenceError(
            "scientific configuration lacks the frozen profiler executable"
        )
    raw_path = frozen.get("path")
    if not isinstance(raw_path, str) or not raw_path:
        raise EvidenceError("frozen profiler executable path is invalid")
    path = nonempty_file(Path(raw_path).expanduser().resolve(), "frozen profiler")
    actual = {
        "path": str(path),
        "sha256": sha256(path),
        "bytes": path.stat().st_size,
    }
    if frozen != actual:
        raise EvidenceError(
            "profiler executable differs from the pilot scientific configuration"
        )
    return path


def build_signature(
    row: dict[str, str],
    executable_version: str,
    executable: Path,
    args: argparse.Namespace,
    scientific_config: dict[str, object],
) -> dict[str, object]:
    return {
        "release": RELEASE,
        "library": row["library"],
        "manifest_row": dict(row),
        "inputs": {
            column: file_signature(Path(row[column]))
            for column in SIGNATURE_INPUT_COLUMNS
            if row.get(column)
        },
        "profiler": executable_version,
        "profiler_executable": file_signature(executable),
        "scientific_configuration_hash": scientific_config["configuration_hash"],
        "parameters": scientific_parameters(row, scientific_config),
    }


def preflight_row(row: dict[str, str]) -> dict[str, Path]:
    paths: dict[str, Path] = {}
    required = (
        "library_dir", "bam", "bam_index", "summary", "raw_barcodes",
        "filtered_barcodes", "raw_features", "raw_matrix", "filtered_matrix",
        "rg_metadata", "source_order", "output_dir",
        "scientific_config",
    )
    for column in required:
        if not row.get(column):
            raise EvidenceError(
                f"manifest row for {row.get('library')} has blank {column}"
            )
        paths[column] = Path(row[column]).expanduser().resolve(strict=False)
    for column in required:
        if column not in {"library_dir", "output_dir"}:
            nonempty_file(paths[column], column.replace("_", " "))
    expected_semantics = {
        "starsolo_feature": STARSOLO_FEATURE,
        "starsolo_umi_filtering": STARSOLO_UMI_FILTERING,
        "starsolo_umi_dedup": STARSOLO_UMI_DEDUP,
        "starsolo_multimappers": STARSOLO_MULTIMAPPERS,
    }
    for column, expected in expected_semantics.items():
        if row.get(column) != expected:
            raise EvidenceError(
                f"manifest {column} must be {expected!r}; got "
                f"{row.get(column)!r}"
            )
    classification_intent = row.get("biological_classification_intent")
    if classification_intent == "manifest":
        if not row.get("class_manifest"):
            raise EvidenceError(
                "classification intent is manifest but class_manifest is blank"
            )
    elif classification_intent == "explicit_unavailable":
        if row.get("class_manifest"):
            raise EvidenceError(
                "explicit_unavailable classification intent conflicts with "
                "class_manifest"
            )
    else:
        raise EvidenceError(
            "biological_classification_intent must be 'manifest' or "
            "'explicit_unavailable'"
        )
    for column in (
        "old_raw_barcodes", "old_filtered_barcodes", "class_manifest",
        "native_cell_reads",
    ):
        if row.get(column):
            paths[column] = Path(row[column]).expanduser().resolve(strict=False)
            nonempty_file(paths[column], column.replace("_", " "))
    library_dir = paths["library_dir"]
    if not library_dir.is_dir():
        raise EvidenceError(f"library directory does not exist: {library_dir}")
    expected_output = library_dir / "bam_evidence"
    if paths["output_dir"] != expected_output:
        raise EvidenceError(
            f"output ownership violation: expected {expected_output}, "
            f"manifest requested {paths['output_dir']}"
        )
    for label in (
        "bam", "bam_index", "summary", "raw_barcodes", "filtered_barcodes",
        "raw_features", "raw_matrix", "filtered_matrix",
    ):
        if paths[label] != library_dir and library_dir not in paths[label].parents:
            raise EvidenceError(
                f"input ownership violation for {label}: {paths[label]}"
            )
    for label, path in paths.items():
        if label in {"library_dir", "output_dir"}:
            continue
        if path == expected_output or expected_output in path.parents:
            raise EvidenceError(
                f"input may not be inside evidence output directory: {path}"
            )
    return paths


def run_checked(command: list[str], label: str) -> subprocess.CompletedProcess[str]:
    try:
        result = subprocess.run(
            command, capture_output=True, text=True, check=False
        )
    except OSError as exc:
        raise EvidenceError(f"could not start {label}: {exc}") from exc
    if result.returncode != 0:
        details = (result.stderr or result.stdout).strip()
        raise EvidenceError(f"{label} failed ({result.returncode}): {details}")
    return result


def bam_preflight(paths: dict[str, Path]) -> dict[str, object]:
    quickcheck = run_checked(
        ["samtools", "quickcheck", str(paths["bam"])], "samtools quickcheck"
    )
    header = run_checked(
        ["samtools", "view", "-H", str(paths["bam"])], "BAM header preflight"
    )
    header_rgs: list[str] = []
    header_contigs: list[str] = []
    for line in header.stdout.splitlines():
        fields = dict(
            field.split(":", 1)
            for field in line.split("\t")[1:]
            if ":" in field
        )
        if line.startswith("@RG\t") and fields.get("ID"):
            header_rgs.append(fields["ID"])
        elif line.startswith("@SQ\t") and fields.get("SN"):
            header_contigs.append(fields["SN"])
    if not header_rgs or len(header_rgs) != len(set(header_rgs)):
        raise EvidenceError("BAM header has no RG IDs or has duplicate RG IDs")
    if not header_contigs or len(header_contigs) != len(set(header_contigs)):
        raise EvidenceError("BAM header has no SQ contigs or has duplicate SQ contigs")

    # samtools 1.20 idxstats has no -X option. Probe a real region with
    # samtools view -X instead so the exact manifest-declared index is opened,
    # rather than silently scanning an unindexed BAM or discovering another
    # index beside it.
    index_probe_region = f"{header_contigs[0]}:1-1"
    index_probe = run_checked(
        [
            "samtools", "view", "-c", "-X", str(paths["bam"]),
            str(paths["bam_index"]), index_probe_region,
        ],
        "BAM index preflight",
    )
    try:
        index_probe_records = int(index_probe.stdout.strip())
    except ValueError as exc:
        raise EvidenceError(
            "BAM index preflight returned a non-integer record count"
        ) from exc
    if index_probe_records < 0:
        raise EvidenceError("BAM index preflight returned a negative record count")
    return {
        "quickcheck_stderr": quickcheck.stderr.strip(),
        "header_sha256": hashlib.sha256(
            header.stdout.encode("utf-8")
        ).hexdigest(),
        "header_rg_ids": sorted(header_rgs),
        "indexed_contigs": len(header_contigs),
        "index_probe_region": index_probe_region,
        "index_probe_records": index_probe_records,
    }


def plain_product_stats(staging: Path) -> dict[str, dict[str, object]]:
    return {
        final: validate_product(staging / plain, final)
        for plain, final in {**COMPRESSED_PRODUCTS, **PLAIN_PRODUCTS}.items()
    }


def validate_profiler_contract(
    profiler_summary: dict[str, str],
    args: argparse.Namespace,
    bam_check: dict[str, object],
    manifest_rg_rows: list[dict[str, str]],
    sources: list[str],
    class_manifest_supplied: bool,
) -> None:
    expected = {
        "hash_algorithm": HASH_ALGORITHM,
        "hash_seed": str(args.hash_seed),
        "hash_bins": str(args.hash_bins),
        "header_rg_count": str(len(bam_check["header_rg_ids"])),
        "manifest_rg_count": str(len(manifest_rg_rows)),
        "source_order_count": str(len(sources)),
        "molecule_slot_bytes": "40",
        "starsolo_feature": STARSOLO_FEATURE,
        "starsolo_umi_filtering": STARSOLO_UMI_FILTERING,
        "starsolo_umi_dedup": STARSOLO_UMI_DEDUP,
        "starsolo_multimappers": STARSOLO_MULTIMAPPERS,
        "ordinary_countedU_read_definition":
            ORDINARY_COUNTEDU_READ_DEFINITION,
        "ordinary_matrix_molecule_definition": ORDINARY_MOLECULE_DEFINITION,
        "summary_unique_read_metric": SUMMARY_UNIQUE_READ_METRIC,
        "nh_gt1_unique_gene_countedU_definition":
            NH_GT1_UNIQUE_GENE_COUNTEDU_DEFINITION,
        "starsolo_EM_evidence_availability":
            STARSOLO_EM_EVIDENCE_AVAILABILITY,
        "within_rg_conflict_definition":
            "same_CR_same_RG_has_multiple_corrected_CB_values",
        "cross_rg_conflict_definition":
            "two_or_more_exact_RG_level_mappings_for_same_CR_disagree",
        "global_conflict_definition":
            "same_CR_has_multiple_corrected_CB_values_anywhere",
        "record_contribution_method":
            "single_tag_decode_and_single_cigar_walk_reused_across_accumulators",
        "ordinary_raw_matrix_exact_coordinate_reconciliation": "PASS",
        "ordinary_filtered_matrix_exact_coordinate_reconciliation": "PASS",
        "biological_classification_manifest_supplied":
            "1" if class_manifest_supplied else "0",
    }
    for metric, value in expected.items():
        if profiler_summary.get(metric) != value:
            raise EvidenceError(
                f"profiler summary mismatch for {metric}: "
                f"{profiler_summary.get(metric)!r} != {value!r}"
            )
    if profiler_summary.get("htslib_version") != "1.20":
        raise EvidenceError(
            "rna_bam_evidence must run against htslib 1.20; reported "
            f"{profiler_summary.get('htslib_version')!r}"
        )
    required_memory_metrics = (
        "configured_total_process_memory_bytes",
        "effective_RLIMIT_AS_bytes",
        "conservative_memory_admission_limit_bytes",
        "accounted_peak_bytes",
        "process_MaxRSS_kib",
        "ordinary_molecule_table_allocated_bytes",
        "barcode_rg_table_allocated_bytes",
        "correction_table_allocated_bytes",
        "runtime_allocator_htslib_reserve_bytes",
        "hash_bin_aggregation_additional_heap_bytes",
    )
    for metric in required_memory_metrics:
        try:
            value = int(profiler_summary[metric])
        except (KeyError, ValueError) as exc:
            raise EvidenceError(
                f"profiler omitted numeric memory metric {metric}"
            ) from exc
        if value < 0:
            raise EvidenceError(f"profiler memory metric {metric} is negative")
    configured = int(profiler_summary["configured_total_process_memory_bytes"])
    expected_configured = int(args.max_memory_gb * 1024**3)
    if configured != expected_configured:
        raise EvidenceError(
            f"profiler process memory cap {configured} != requested {expected_configured}"
        )
    effective = int(profiler_summary["effective_RLIMIT_AS_bytes"])
    admission = int(
        profiler_summary["conservative_memory_admission_limit_bytes"]
    )
    runtime_reserve = int(
        profiler_summary["runtime_allocator_htslib_reserve_bytes"]
    )
    if admission + runtime_reserve != effective or effective > configured:
        raise EvidenceError("profiler memory admission and runtime reserve are inconsistent")
    if int(profiler_summary["process_MaxRSS_kib"]) * 1024 > effective:
        raise EvidenceError("profiler MaxRSS exceeded its declared process memory cap")
    if (
        profiler_summary.get("hash_bin_aggregation_method")
        != "in_place_sorted_ordinary_molecule_slots"
        or profiler_summary.get("output_sorting_method")
        != "in_place_high_cardinality_tables_plus_numeric_barcode_order_vectors"
    ):
        raise EvidenceError("profiler did not report bounded in-place aggregation/sorting")
    status = profiler_summary.get("biological_classification_status")
    if not class_manifest_supplied and status != "unavailable":
        raise EvidenceError(
            "missing class manifest must produce unavailable biological metrics"
        )


def integer_field(row: dict[str, str], field: str) -> int:
    value = row.get(field, "")
    if value == "" and field in OPTIONAL_CLASS_COUNT_FIELDS:
        return 0
    try:
        return int(value)
    except ValueError as exc:
        raise EvidenceError(f"non-integer {field} value {value!r}") from exc


def aggregate_barcode_rg(
    path: Path,
) -> Iterator[tuple[str, dict[str, int]]]:
    current_cb: str | None = None
    totals = {field: 0 for field in CORE_ADDITIVE_FIELDS}
    for row in iter_tsv(path, ("CB", "RG", *CORE_ADDITIVE_FIELDS)):
        cb = row["CB"]
        if current_cb is not None and cb != current_cb:
            yield current_cb, totals
            totals = {field: 0 for field in CORE_ADDITIVE_FIELDS}
        current_cb = cb
        for field in CORE_ADDITIVE_FIELDS:
            totals[field] += integer_field(row, field)
    if current_cb is not None:
        yield current_cb, totals


def validate_correction_conflicts(path: Path) -> None:
    """Validate CR groups in constant space from CR,RG,CB sorted input."""
    current_cr: str | None = None
    current_rg: str | None = None
    global_first_cb = ""
    rg_first_cb = ""
    first_exact_rg_cb = ""
    have_exact_rg = False
    global_conflict = False
    cross_conflict = False
    rg_distinct = 0
    last_cb = ""
    within_flags: set[str] = set()
    global_flags: set[str] = set()
    cross_flags: set[str] = set()

    def finish_rg() -> None:
        nonlocal first_exact_rg_cb, have_exact_rg, cross_conflict
        if current_rg is None:
            return
        expected = "1" if rg_distinct > 1 else "0"
        if within_flags != {expected}:
            raise EvidenceError(
                f"inconsistent within-RG correction conflict for "
                f"{current_cr}/{current_rg}"
            )
        if rg_distinct == 1:
            if not have_exact_rg:
                first_exact_rg_cb = rg_first_cb
                have_exact_rg = True
            elif rg_first_cb != first_exact_rg_cb:
                cross_conflict = True

    def finish_cr() -> None:
        if current_cr is None:
            return
        finish_rg()
        expected_global = "1" if global_conflict else "0"
        expected_cross = "1" if cross_conflict else "0"
        if global_flags != {expected_global}:
            raise EvidenceError(
                f"inconsistent global correction conflict for {current_cr}"
            )
        if cross_flags != {expected_cross}:
            raise EvidenceError(
                f"inconsistent cross-RG correction conflict for {current_cr}"
            )

    for row in iter_tsv(
        path,
        (
            "CR", "RG", "CB", "read_count", "within_rg_conflict",
            "cross_rg_conflict", "global_conflict",
        ),
    ):
        cr, rg, cb = row["CR"], row["RG"], row["CB"]
        if int(row["read_count"]) < 1:
            raise EvidenceError(f"nonpositive correction count for {cr}/{rg}/{cb}")
        if cr != current_cr:
            finish_cr()
            current_cr, current_rg = cr, rg
            global_first_cb = cb
            rg_first_cb = cb
            first_exact_rg_cb = ""
            have_exact_rg = False
            global_conflict = False
            cross_conflict = False
            rg_distinct = 1
            last_cb = cb
            within_flags = {row["within_rg_conflict"]}
            global_flags = {row["global_conflict"]}
            cross_flags = {row["cross_rg_conflict"]}
            continue
        global_flags.add(row["global_conflict"])
        cross_flags.add(row["cross_rg_conflict"])
        if cb != global_first_cb:
            global_conflict = True
        if rg != current_rg:
            finish_rg()
            current_rg = rg
            rg_first_cb = cb
            rg_distinct = 1
            last_cb = cb
            within_flags = {row["within_rg_conflict"]}
        else:
            within_flags.add(row["within_rg_conflict"])
            if cb != last_cb:
                rg_distinct += 1
                last_cb = cb
    finish_cr()


def validate_internal_reconciliation(
    staging: Path,
    profiler_summary: dict[str, str],
    expected_hash_bins: int,
) -> dict[str, object]:
    """Reconcile sorted profiler products with bounded streaming state."""
    rg_rows = list(
        iter_tsv(staging / "rg_summary.tsv", ("RG", "source_id", "source_index"))
    )
    source_by_index: dict[int, str] = {}
    for row in rg_rows:
        index = int(row["source_index"])
        if index == -1 and row["RG"] == "__MISSING_OR_UNDECLARED__":
            continue
        if not 0 <= index < 64:
            raise EvidenceError(f"invalid source index {index}")
        prior = source_by_index.setdefault(index, row["source_id"])
        if prior != row["source_id"]:
            raise EvidenceError(
                f"source index {index} maps to both {prior} and {row['source_id']}"
            )
    valid_mask = sum(1 << index for index in source_by_index)
    molecule_bin_total = 0
    for row in iter_tsv(staging / "molecule_source_hash_bins.tsv"):
        mask = int(row["source_mask_hex"], 16)
        bin_index = int(row["nested_min_hash_bin"])
        bins = int(row["n_hash_bins"])
        count = int(row["candidate_matrix_molecules"])
        if mask <= 0 or mask & ~valid_mask:
            raise EvidenceError(
                f"invalid molecule source mask {row['source_mask_hex']}"
            )
        if bins != expected_hash_bins or not 0 <= bin_index < bins or count < 1:
            raise EvidenceError("invalid molecule nested hash-bin metadata")
        molecule_bin_total += count
    expected_molecules = int(profiler_summary["candidate_matrix_molecules"])
    if molecule_bin_total != expected_molecules:
        raise EvidenceError(
            f"molecule-bin total {molecule_bin_total} != {expected_molecules}"
        )

    barcode_molecules = 0
    compare_barcode_rg = (
        int(profiler_summary["records_rg_not_manifest"]) == 0
        and int(profiler_summary["missing_RG"]) == 0
        and int(profiler_summary["malformed_RG"]) == 0
    )
    rg_aggregate = iter(aggregate_barcode_rg(staging / "barcode_rg_metrics.tsv"))
    next_rg = next(rg_aggregate, None)
    for barcode in iter_tsv(
        staging / "barcode_read_metrics.tsv",
        ("CB", "candidate_matrix_molecules", *CORE_ADDITIVE_FIELDS),
    ):
        cb = barcode["CB"]
        barcode_molecules += int(barcode["candidate_matrix_molecules"])
        if not compare_barcode_rg:
            continue
        while next_rg is not None and next_rg[0] < cb:
            raise EvidenceError(f"barcode-RG row has no barcode row: {next_rg[0]}")
        totals = (
            next_rg[1]
            if next_rg is not None and next_rg[0] == cb
            else {field: 0 for field in CORE_ADDITIVE_FIELDS}
        )
        for field in CORE_ADDITIVE_FIELDS:
            observed = totals[field]
            expected = integer_field(barcode, field)
            if observed != expected:
                raise EvidenceError(
                    f"barcode/RG mismatch for {cb} {field}: "
                    f"{observed} != {expected}"
                )
        if next_rg is not None and next_rg[0] == cb:
            next_rg = next(rg_aggregate, None)
    if compare_barcode_rg and next_rg is not None:
        raise EvidenceError(f"trailing barcode-RG row for {next_rg[0]}")
    if barcode_molecules != expected_molecules:
        raise EvidenceError(
            f"barcode molecule total {barcode_molecules} != {expected_molecules}"
        )

    additive = {
        field: sum(integer_field(row, field) for row in rg_rows)
        for field in CORE_ADDITIVE_FIELDS
    }
    for field, summary_key in (
        ("all_records", "total_records"),
        ("primary_mapped_reads", "primary_mapped_reads"),
        ("candidate_countedU_reads", "candidate_countedU_reads"),
        (
            "nh_gt1_unique_gene_countedU_reads",
            "nh_gt1_unique_gene_countedU_reads",
        ),
    ):
        expected = int(profiler_summary[summary_key])
        if additive[field] != expected:
            raise EvidenceError(
                f"RG additive {field} total {additive[field]} != {expected}"
            )
    validate_correction_conflicts(
        staging / "raw_to_corrected_barcode_counts.tsv"
    )
    return {
        "streaming_adjacent_key_validation": True,
        "molecule_bins_equal_library_total": True,
        "barcode_molecules_equal_library_total": True,
        "additive_rg_read_totals_equal_library_totals": True,
        "additive_barcode_rg_read_totals_equal_barcode_totals":
            compare_barcode_rg,
        "within_and_cross_rg_correction_flags_reconciled": True,
    }


def read_native_counted_u(path: Path, filtered: set[str]) -> dict[str, int]:
    result: dict[str, int] = {}
    with open_text(path) as handle:
        header = handle.readline().split()
        if not {"CB", "countedU"}.issubset(header):
            raise EvidenceError(
                f"native STARsolo cell-read table lacks CB/countedU: {path}"
            )
        cb_index = header.index("CB")
        counted_index = header.index("countedU")
        for line_number, line in enumerate(handle, start=2):
            fields = line.split()
            if len(fields) <= max(cb_index, counted_index):
                raise EvidenceError(
                    f"malformed native cell-read row at {path}:{line_number}"
                )
            cb = fields[cb_index]
            if cb not in filtered:
                continue
            if cb in result:
                raise EvidenceError(
                    f"duplicate native cell barcode in {path}: {cb}"
                )
            try:
                result[cb] = int(fields[counted_index])
            except ValueError as exc:
                raise EvidenceError(
                    f"non-integer native countedU at {path}:{line_number}"
                ) from exc
    return result


def reconcile(
    row: dict[str, str],
    paths: dict[str, Path],
    staging: Path,
    profiler_summary: dict[str, str],
) -> tuple[dict[str, object], list[tuple[str, int]]]:
    """Exact ordinary STARsolo gate using only fixed filtered-cell vectors."""
    del row
    filtered = roster(paths["filtered_barcodes"])
    filtered_index = {cb: index for index, cb in enumerate(filtered)}
    raw_count = roster_count(paths["raw_barcodes"])
    matrix_stats, matrix_umis, matrix_genes = matrix_shape_and_counts(
        paths["filtered_matrix"], len(filtered), keep_per_barcode=True
    )
    assert matrix_umis is not None and matrix_genes is not None
    candidate_reads: list[int | None] = [None] * len(filtered)
    candidate_molecules: list[int | None] = [None] * len(filtered)
    for metric in iter_tsv(
        staging / "barcode_read_metrics.tsv",
        ("CB", "candidate_countedU_reads", "candidate_matrix_molecules"),
    ):
        index = filtered_index.get(metric["CB"])
        if index is None:
            continue
        candidate_reads[index] = int(metric["candidate_countedU_reads"])
        candidate_molecules[index] = int(metric["candidate_matrix_molecules"])
    missing = [
        filtered[index]
        for index, value in enumerate(candidate_reads)
        if value is None
    ]
    if missing:
        raise EvidenceError(
            f"profiler omitted {len(missing)} filtered barcodes; examples: "
            + ", ".join(missing[:3])
        )
    reads = [int(value) for value in candidate_reads if value is not None]
    molecules = [
        int(value) for value in candidate_molecules if value is not None
    ]
    mismatch_count = 0
    mismatch_examples: list[dict[str, object]] = []
    for cb, observed, expected in zip(filtered, molecules, matrix_umis):
        if observed == expected:
            continue
        mismatch_count += 1
        if len(mismatch_examples) < 20:
            mismatch_examples.append(
                {"CB": cb, "observed": observed, "expected": expected}
            )

    summary = summary_values(paths["summary"])
    if SUMMARY_UNIQUE_READ_METRIC not in summary:
        raise EvidenceError(
            f"STARsolo summary lacks exact metric {SUMMARY_UNIQUE_READ_METRIC!r}: "
            f"{paths['summary']}"
        )
    for key in ("Estimated Number of Cells", "Median Reads per Cell"):
        if key not in summary:
            raise EvidenceError(f"STARsolo summary lacks metric {key!r}")
    observed_median = sorted(reads)[len(reads) // 2]
    expected_cells = int(round(summary["Estimated Number of Cells"]))
    expected_reads = int(round(summary[SUMMARY_UNIQUE_READ_METRIC]))
    expected_median = int(round(summary["Median Reads per Cell"]))
    equality: dict[str, dict[str, object]] = {
        "compiled_filtered_sparse_coordinate_reconciliation": {
            "observed_molecules": int(
                profiler_summary["ordinary_filtered_matrix_molecules"]
            ),
            "expected_molecules": matrix_stats["molecules"],
            "observed_nnz": int(
                profiler_summary[
                    "ordinary_filtered_matrix_nonzero_coordinates"
                ]
            ),
            "expected_nnz": matrix_stats["nnz"],
            "equal": (
                int(profiler_summary["ordinary_filtered_matrix_molecules"])
                == matrix_stats["molecules"]
                and int(
                    profiler_summary[
                        "ordinary_filtered_matrix_nonzero_coordinates"
                    ]
                ) == matrix_stats["nnz"]
            ),
        },
        "filtered_barcode_count": {
            "observed": len(filtered),
            "expected": expected_cells,
            "equal": len(filtered) == expected_cells,
        },
        "filtered_ordinary_countedU_read_total": {
            "observed": sum(reads),
            "expected": expected_reads,
            "equal": sum(reads) == expected_reads,
        },
        "filtered_ordinary_countedU_read_star_median": {
            "observed": observed_median,
            "expected": expected_median,
            "equal": observed_median == expected_median,
        },
        "filtered_ordinary_molecule_total": {
            "observed": sum(molecules),
            "expected": matrix_stats["molecules"],
            "equal": sum(molecules) == matrix_stats["molecules"],
        },
        "filtered_ordinary_molecule_per_barcode": {
            "mismatches": mismatch_count,
            "equal": mismatch_count == 0,
            "examples": mismatch_examples,
        },
    }
    declaration_failures = {
        "records_rg_not_header":
            int(profiler_summary.get("records_rg_not_header", "-1")),
        "records_rg_not_manifest":
            int(profiler_summary.get("records_rg_not_manifest", "-1")),
        "records_missing_rg": int(profiler_summary.get("missing_RG", "-1")),
        "records_malformed_rg":
            int(profiler_summary.get("malformed_RG", "-1")),
    }
    equality["all_record_rg_declarations"] = {
        **declaration_failures,
        "equal": all(value == 0 for value in declaration_failures.values()),
    }
    rg_rows = list(iter_tsv(staging / "rg_summary.tsv", ("RG", "RG_in_BAM_header")))
    missing_header = [
        item["RG"]
        for item in rg_rows
        if item["RG"] != "__MISSING_OR_UNDECLARED__"
        and item["RG_in_BAM_header"] != "1"
    ]
    equality["manifest_rgs_present_in_header"] = {
        "missing": missing_header,
        "equal": not missing_header,
    }

    native_result: dict[str, object] = {"supplied": False, "equal": None}
    if "native_cell_reads" in paths:
        native = read_native_counted_u(
            paths["native_cell_reads"], set(filtered)
        )
        difference_count = 0
        examples: list[dict[str, object]] = []
        for index, cb in enumerate(filtered):
            if native.get(cb) == reads[index]:
                continue
            difference_count += 1
            if len(examples) < 20:
                examples.append(
                    {"CB": cb, "profiler": reads[index], "native": native.get(cb)}
                )
        native_result = {
            "supplied": True,
            "equal": difference_count == 0 and len(native) == len(filtered),
            "mismatches": difference_count,
            "examples": examples,
        }
    equality["native_countedU_when_supplied"] = native_result
    failed = [
        name for name, result in equality.items()
        if result.get("equal") is False
    ]
    reconciliation = {
        "status": "PASS" if not failed else "FAIL",
        "failed_equalities": failed,
        "semantic_scope": {
            "ordinary_matrix": {
                "molecule_definition": ORDINARY_MOLECULE_DEFINITION,
                "producing_starsolo_options": {
                    "soloFeatures": STARSOLO_FEATURE,
                    "soloUMIfiltering": STARSOLO_UMI_FILTERING,
                    "soloUMIdedup": STARSOLO_UMI_DEDUP,
                    "soloMultiMappers": STARSOLO_MULTIMAPPERS,
                },
                "ub_contract": (
                    "A matrix molecule requires a valid STARsolo corrected UB after "
                    "1MM_CR and MultiGeneUMI_CR; '-' means unavailable/rejected and "
                    "cannot form a matrix molecule"
                ),
                "targets": [
                    str(paths["raw_matrix"]), str(paths["filtered_matrix"]),
                ],
                "raw_sparse_coordinate_gate": "compiled_exact_merge",
                "filtered_per_barcode_gate": "wrapper_streamed_matrix_vector",
                "exact_gate": True,
            },
            "ordinary_countedU_reads": {
                "record_definition": ORDINARY_COUNTEDU_READ_DEFINITION,
                "summary_metric": SUMMARY_UNIQUE_READ_METRIC,
                "target": str(paths["summary"]),
                "ub_contract": (
                    "STARsolo countedU and Summary unique reads are recorded before "
                    "UMI collapse/filtering; a valid UB is not required and UB '-' "
                    "does not exclude an otherwise counted unique-gene read"
                ),
            },
            "nh_gt1_unique_gene_countedU": {
                "definition": NH_GT1_UNIQUE_GENE_COUNTEDU_DEFINITION,
                "reported_metric": "nh_gt1_unique_gene_countedU_reads",
                "label": "ordinary_unique_gene_countedU_NH_gt1_subset",
                "included_in_ordinary_countedU": True,
                "NH_is_EM_membership": False,
            },
            "starsolo_multi_gene_EM": {
                "definition": MULTIMAPPER_DEFINITION,
                "availability": STARSOLO_EM_EVIDENCE_AVAILABILITY,
                "reported_metric": None,
                "ordinary_matrix_promotion": False,
                "em_matrix_reconciliation_performed": False,
                "complete_STAR_EM_evidence_claimed": False,
            },
        },
        "summary_unique_read_metric": SUMMARY_UNIQUE_READ_METRIC,
        "matrix": matrix_stats,
        "raw_barcode_count": raw_count,
        "equalities": equality,
        "candidate_molecule_label":
            "validated_ordinary_matrix_molecules" if not failed
            else "candidate_ordinary_matrix_molecules",
        "per_cell_gene_counts": {
            "claimed_from_bam": False,
            "authoritative_source": str(paths["filtered_matrix"]),
            "matrix_nonzero_gene_counts_loaded": sum(matrix_genes),
        },
    }
    if failed:
        raise ReconciliationFailure(reconciliation)
    return reconciliation, sorted(zip(filtered, reads))


def deterministic_gzip(source: Path, destination: Path) -> None:
    temporary = destination.with_name(f".{destination.name}.tmp.{os.getpid()}")
    with source.open("rb") as input_handle, temporary.open("wb") as raw_output:
        with gzip.GzipFile(
            filename="", mode="wb", fileobj=raw_output, mtime=0
        ) as output:
            shutil.copyfileobj(input_handle, output, length=1024 * 1024)
    os.replace(temporary, destination)


def publish_plain(source: Path, destination: Path) -> None:
    temporary = destination.with_name(f".{destination.name}.tmp.{os.getpid()}")
    shutil.copyfile(source, temporary)
    os.replace(temporary, destination)


def publish_profiler_products(staging: Path, output: Path) -> None:
    output.mkdir(parents=True, exist_ok=True)
    for source_name, destination_name in COMPRESSED_PRODUCTS.items():
        deterministic_gzip(staging / source_name, output / destination_name)
    for source_name, destination_name in PLAIN_PRODUCTS.items():
        publish_plain(staging / source_name, output / destination_name)


def write_compatibility(path: Path, rows: Iterable[tuple[str, int]]) -> None:
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    with temporary.open("wb") as raw_output:
        with gzip.GzipFile(
            filename="", mode="wb", fileobj=raw_output, mtime=0
        ) as compressed:
            with io.TextIOWrapper(
                compressed, encoding="utf-8", newline=""
            ) as text:
                writer = csv.writer(text, delimiter="\t", lineterminator="\n")
                writer.writerow(["CB", "countedU"])
                writer.writerows(rows)
    os.replace(temporary, path)


def validate_published(output: Path) -> dict[str, dict[str, object]]:
    return {
        name: validate_product(output / name, name)
        for name in REQUIRED_PRODUCTS
    }


def validate_compatibility(path: Path) -> dict[str, object]:
    nonempty_file(path, "validated compatibility countedU product")
    digest = hashlib.sha256()
    previous: str | None = None
    rows = 0
    with open_hashed_text(path, digest) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        if not {"CB", "countedU"}.issubset(fields):
            raise EvidenceError(f"compatibility table lacks CB/countedU: {path}")
        for line_number, row in enumerate(reader, start=2):
            cb = row["CB"]
            if not cb or (previous is not None and cb <= previous):
                raise EvidenceError(
                    f"blank, duplicate, or out-of-order CB in {path}:{line_number}"
                )
            previous = cb
            try:
                count = int(row["countedU"])
            except ValueError as exc:
                raise EvidenceError(
                    f"non-integer countedU in {path}:{line_number}"
                ) from exc
            if count < 0:
                raise EvidenceError(f"negative countedU in {path}:{line_number}")
            rows += 1
    if not rows:
        raise EvidenceError(f"compatibility table has no barcode rows: {path}")
    return {
        "rows": rows,
        "bytes": path.stat().st_size,
        "sha256": digest.hexdigest(),
        "schema": fields,
    }


def ensure_stats_match(
    actual: dict[str, object],
    expected: dict[str, object],
    label: str,
) -> None:
    for field in ("rows", "bytes", "sha256", "schema"):
        if actual.get(field) != expected.get(field):
            raise EvidenceError(
                f"{label} no longer matches recorded audit field {field}"
            )


def existing_is_valid(
    output: Path,
    signature: dict[str, object],
    *,
    allow_stale_replacement: bool,
) -> bool:
    marker = output / "BAM_EVIDENCE_COMPLETE.ok"
    audit_path = output / "audit.json"
    if not marker.is_file() or not audit_path.is_file():
        return False
    nonempty_file(marker, "evidence completion marker")
    sha256(marker)  # also prove the marker is readable before accepting resume
    try:
        audit = json.loads(audit_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise EvidenceError(f"could not read existing audit {audit_path}: {exc}") from exc
    if audit.get("resume_signature") != signature:
        if allow_stale_replacement:
            return False
        raise EvidenceError(
            f"completed evidence output is stale for {signature['library']}; "
            "remove it deliberately or use --replace-stale"
        )
    if (
        audit.get("library") != signature["library"]
        or audit.get("reconciliation", {}).get("status") != "PASS"
    ):
        raise EvidenceError(
            f"completed evidence audit is not a passing match: {audit_path}"
        )
    expected_products = audit.get("outputs")
    if not isinstance(expected_products, dict):
        raise EvidenceError(f"audit lacks recorded output inventory: {audit_path}")
    actual_products = validate_published(output)
    if set(actual_products) != set(expected_products):
        raise EvidenceError("published product inventory differs from audit")
    for name, actual in actual_products.items():
        ensure_stats_match(actual, expected_products[name], f"{output / name}")
    actual_compatibility = validate_compatibility(
        output / "CellReads.countedU.from_bam.tsv.gz"
    )
    expected_compatibility = audit.get("compatibility_product", {})
    ensure_stats_match(
        actual_compatibility, expected_compatibility,
        str(output / "CellReads.countedU.from_bam.tsv.gz"),
    )
    return True


@contextmanager
def exclusive_output_lock(output: Path) -> Iterator[None]:
    """Prevent two allocations from publishing one library concurrently."""
    output.mkdir(parents=True, exist_ok=True)
    lock_path = output / ".rna_bam_evidence.lock"
    with lock_path.open("a+", encoding="utf-8") as handle:
        try:
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise EvidenceError(
                "another RNA BAM evidence process holds the library output "
                f"lock: {lock_path}"
            ) from exc
        try:
            yield
        finally:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


def run_one(args: argparse.Namespace) -> None:
    """Resolve output ownership, lock it, and run one manifest row."""
    manifest = Path(args.manifest).expanduser().resolve()
    _manifest_fields, manifest_rows = read_manifest(manifest)
    if args.row_index < 0 or args.row_index >= len(manifest_rows):
        raise EvidenceError(
            f"manifest row index {args.row_index} is outside "
            f"0..{len(manifest_rows) - 1}"
        )
    # preflight_row validates that output_dir is exactly this library's
    # bam_evidence directory before the lock file is created.
    output = preflight_row(manifest_rows[args.row_index])["output_dir"]
    with exclusive_output_lock(output):
        _run_one_locked(args)


def _run_one_locked(args: argparse.Namespace) -> None:
    manifest = Path(args.manifest).expanduser().resolve()
    manifest_fields, manifest_rows = read_manifest(manifest)
    scientific_config = load_project_scientific_config(
        manifest, manifest_fields, manifest_rows
    )
    if args.row_index < 0 or args.row_index >= len(manifest_rows):
        raise EvidenceError(
            f"manifest row index {args.row_index} is outside 0..{len(manifest_rows) - 1}"
        )
    row = manifest_rows[args.row_index]
    paths = preflight_row(row)
    executable = (
        shutil.which(args.profiler)
        if os.sep not in args.profiler
        else args.profiler
    )
    if not executable:
        raise EvidenceError(f"rna_bam_evidence executable not found: {args.profiler}")
    executable_path = Path(executable).resolve()
    if executable_path != frozen_profiler_path(scientific_config):
        raise EvidenceError(
            "requested profiler path differs from the frozen scientific "
            "configuration; refusing to scan the BAM"
        )
    version = profiler_version(str(executable_path))
    expected_parameters = scientific_parameters(row, scientific_config)
    if (
        args.hash_seed != expected_parameters["hash_seed"]
        or args.hash_bins != expected_parameters["hash_bins"]
    ):
        raise EvidenceError(
            "requested hash seed/bins differ from the frozen scientific configuration"
        )
    signature = build_signature(
        row, version, executable_path, args, scientific_config
    )
    execution_identity = slurm_execution_identity()
    output = paths["output_dir"]
    existing_valid = output.exists() and existing_is_valid(
        output, signature, allow_stale_replacement=args.replace_stale
    )
    if existing_valid:
        print(f"Validated existing RNA BAM evidence: {output}")
        return

    output.mkdir(parents=True, exist_ok=True)
    (output / "BAM_EVIDENCE_COMPLETE.ok").unlink(missing_ok=True)
    (output / "audit.json").unlink(missing_ok=True)
    (output / "BAM_EVIDENCE_FAILED.json").unlink(missing_ok=True)
    # This compatibility file is the only compact product intended for legacy
    # direct consumers.  Never leave a previously validated copy beside a new
    # failed/incomplete scan where it could be mistaken for current evidence.
    (output / "CellReads.countedU.from_bam.tsv.gz").unlink(missing_ok=True)
    bam_check = bam_preflight(paths)
    sources = source_order(paths["source_order"])
    manifest_rg_rows = rg_manifest_rows(paths["rg_metadata"], row["library"])
    rg_ids = {item["rg_id"] for item in manifest_rg_rows}
    if rg_ids != set(bam_check["header_rg_ids"]):
        raise EvidenceError(
            f"BAM header RG set does not exactly match manifest for {row['library']}"
        )
    source_set = set(sources)
    if {item["bp_id"] for item in manifest_rg_rows} - source_set:
        raise EvidenceError("RG metadata contains a source absent from source order")

    raw_barcode_count = roster_count(paths["raw_barcodes"])
    raw_feature_count = roster_count(paths["raw_features"])
    raw_matrix_stats, _raw_umis, _raw_genes = matrix_shape_and_counts(
        paths["raw_matrix"], raw_barcode_count, keep_per_barcode=False
    )
    if raw_matrix_stats["features"] != raw_feature_count:
        raise EvidenceError(
            "raw matrix feature dimension does not match ordinary feature roster"
        )

    staging = Path(tempfile.mkdtemp(prefix=".rna_bam_evidence.", dir=output))
    memory_bytes = int(args.max_memory_gb * 1024**3)
    command = [
        str(executable_path),
        "--bam", str(paths["bam"]),
        "--output-dir", str(staging),
        "--library", row["library"],
        "--raw-barcodes", str(paths["raw_barcodes"]),
        "--raw-matrix", str(paths["raw_matrix"]),
        "--filtered-barcodes", str(paths["filtered_barcodes"]),
        "--filtered-matrix", str(paths["filtered_matrix"]),
        "--features", str(paths["raw_features"]),
        "--rg-metadata", str(paths["rg_metadata"]),
        "--source-order", str(paths["source_order"]),
        "--starsolo-feature", row["starsolo_feature"],
        "--starsolo-umi-filtering", row["starsolo_umi_filtering"],
        "--starsolo-umi-dedup", row["starsolo_umi_dedup"],
        "--starsolo-multimappers", row["starsolo_multimappers"],
        "--threads", str(args.threads),
        "--hash-bins", str(args.hash_bins),
        "--hash-seed", str(args.hash_seed),
        "--expected-molecules", str(raw_matrix_stats["molecules"]),
        "--max-memory-bytes", str(memory_bytes),
    ]
    for manifest_key, option in (
        ("old_raw_barcodes", "--old-raw-barcodes"),
        ("old_filtered_barcodes", "--old-filtered-barcodes"),
        ("class_manifest", "--class-manifest"),
    ):
        if manifest_key in paths:
            command.extend([option, str(paths[manifest_key])])

    started = time.monotonic()
    before_usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    try:
        try:
            result = subprocess.run(
                command, capture_output=True, text=True, check=False
            )
        except OSError as exc:
            raise EvidenceError(f"could not start compiled profiler: {exc}") from exc
        elapsed = time.monotonic() - started
        if result.returncode != 0:
            raise EvidenceError(
                f"compiled profiler failed for {row['library']} "
                f"({result.returncode}): {result.stderr.strip()}"
            )
        plain_product_stats(staging)
        profiler_summary = parse_profiler_summary(staging / "profiler_summary.tsv")
        if (
            int(profiler_summary.get("ordinary_raw_matrix_molecules", "-1"))
            != raw_matrix_stats["molecules"]
            or int(
                profiler_summary.get(
                    "ordinary_raw_matrix_nonzero_coordinates", "-1"
                )
            ) != raw_matrix_stats["nnz"]
        ):
            raise EvidenceError(
                "compiled raw-matrix reconciliation metrics differ from "
                "the wrapper's streamed MatrixMarket dimensions"
            )
        validate_profiler_contract(
            profiler_summary,
            args,
            bam_check,
            manifest_rg_rows,
            sources,
            "class_manifest" in paths,
        )
        internal = validate_internal_reconciliation(
            staging, profiler_summary, args.hash_bins
        )
        try:
            reconciliation, compatibility_rows = reconcile(
                row, paths, staging, profiler_summary
            )
        except ReconciliationFailure as exc:
            publish_profiler_products(staging, output)
            atomic_json(
                output / "BAM_EVIDENCE_FAILED.json",
                {
                    "release": RELEASE,
                    "failed_utc": utc_now(),
                    "library": row["library"],
                    "reason": str(exc),
                    "reconciliation": exc.details,
                    "profiler_summary": profiler_summary,
                    "command": command,
                    "limitations": [
                        "Ordinary candidate labels were not promoted.",
                        "Multimapper evidence was never treated as ordinary matrix evidence.",
                        "No completion marker or compatibility table was written.",
                    ],
                },
            )
            raise

        publish_profiler_products(staging, output)
        write_compatibility(
            output / "CellReads.countedU.from_bam.tsv.gz",
            compatibility_rows,
        )
        published_stats = validate_published(output)
        compatibility_stats = validate_compatibility(
            output / "CellReads.countedU.from_bam.tsv.gz"
        )
        after_usage = resource.getrusage(resource.RUSAGE_CHILDREN)
        wrapper_child_peak = max(
            before_usage.ru_maxrss, after_usage.ru_maxrss
        )
        audit = {
            "release": RELEASE,
            "created_utc": utc_now(),
            "library": row["library"],
            "resume_signature": signature,
            "manifest": str(manifest),
            "manifest_row_index": args.row_index,
            "slurm_execution": execution_identity,
            "input_paths": {key: str(value) for key, value in paths.items()},
            "bam": file_signature(paths["bam"]),
            "bam_preflight": bam_check,
            "rg_manifest": str(paths["rg_metadata"]),
            "rg_manifest_rows": manifest_rg_rows,
            "source_order": sources,
            "source_order_provenance":
                row.get("source_order_provenance", "manifest"),
            "profiler_version": version,
            "profiler_summary": profiler_summary,
            "internal_reconciliation": internal,
            "command": command,
            "parameters": signature["parameters"],
            "execution_parameters": {
                "threads": args.threads,
                "total_process_memory_gb": args.max_memory_gb,
            },
            "scientific_configuration": {
                "path": str(paths["scientific_config"]),
                "configuration_hash": scientific_config["configuration_hash"],
            },
            "semantics": reconciliation["semantic_scope"],
            "hash": {
                "algorithm": HASH_ALGORITHM,
                "seed": args.hash_seed,
                "bins": args.hash_bins,
                "molecule_value":
                    "minimum QNAME hash across ordinary contributing reads; "
                    "genomic NH is unrestricted",
            },
            "raw_matrix_presizing": raw_matrix_stats,
            "raw_matrix_reconciliation": {
                "status": "PASS",
                "exact_sparse_coordinates_and_counts": True,
                "matrix": str(paths["raw_matrix"]),
                "nonzero_coordinates": raw_matrix_stats["nnz"],
                "molecules": raw_matrix_stats["molecules"],
            },
            "reconciliation": reconciliation,
            "biological_classification": {
                "manifest_supplied": "class_manifest" in paths,
                "intent": row["biological_classification_intent"],
                "status":
                    profiler_summary["biological_classification_status"],
                "unavailable_values_are_blank": True,
            },
            "outputs": published_stats,
            "compatibility_product": {
                "path": str(output / "CellReads.countedU.from_bam.tsv.gz"),
                "label": "ordinary_STARsolo_countedU_from_bam",
                "accepted_after_exact_ordinary_reconciliation": True,
                **compatibility_stats,
            },
            "performance": {
                "elapsed_seconds": elapsed,
                "bam_bytes_per_second":
                    paths["bam"].stat().st_size / max(elapsed, 1e-9),
                "profiler_process_MaxRSS_kib":
                    int(profiler_summary["process_MaxRSS_kib"]),
                "wrapper_RUSAGE_CHILDREN_MaxRSS_kib": wrapper_child_peak,
                "htslib_bgzf_threads": args.threads,
                "logical_bam_record_loops": 1,
            },
            "limitations": [
                "BAM evidence does not provide all-input valid-barcode rates or raw pool representation.",
                "The BAM duplicate flag is not an RNA UMI or PCR duplicate rate.",
                "Random thinning is conditional on the final BAM and is not an exact historical mapping.",
                "Per-cell gene counts remain authoritative to the STARsolo matrix.",
                "NH>1 singleton-gene reads are ordinary countedU evidence; NH is not STARsolo EM membership.",
                "Exact multi-gene STARsolo EM evidence is unavailable from standard uppercase GX/UB BAM tags.",
            ],
        }
        atomic_json(output / "audit.json", audit)
        (output / "BAM_EVIDENCE_FAILED.json").unlink(missing_ok=True)
        atomic_text(output / "BAM_EVIDENCE_COMPLETE.ok", utc_now() + "\n")
        print(
            f"Published validated ordinary RNA BAM evidence for {row['library']}: "
            f"{profiler_summary.get('total_records', '?')} records in "
            f"{elapsed:.1f}s"
        )
    finally:
        shutil.rmtree(staging, ignore_errors=True)


def write_tsv(
    path: Path,
    fields: list[str],
    rows: Iterable[dict[str, object]],
    *,
    gzip_output: bool = False,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    if gzip_output:
        binary = temporary.open("wb")
        compressed = gzip.GzipFile(
            filename="", mode="wb", fileobj=binary, mtime=0
        )
        handle: TextIO = io.TextIOWrapper(
            compressed, encoding="utf-8", newline=""
        )
    else:
        binary = None
        compressed = None
        handle = temporary.open("w", encoding="utf-8", newline="")
    try:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=fields,
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    finally:
        handle.close()
        if compressed is not None and not compressed.closed:
            compressed.close()
        if binary is not None and not binary.closed:
            binary.close()
    os.replace(temporary, path)


class AtomicTsvWriter:
    def __init__(self, path: Path, fields: list[str], compressed: bool = False):
        self.path = path
        self.fields = fields
        self.compressed_output = compressed
        self.temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
        self.handle: TextIO | None = None
        self.binary = None
        self.compressed = None
        self.writer: csv.DictWriter[str] | None = None
        self.rows = 0

    def __enter__(self) -> "AtomicTsvWriter":
        self.path.parent.mkdir(parents=True, exist_ok=True)
        if self.compressed_output:
            self.binary = self.temporary.open("wb")
            self.compressed = gzip.GzipFile(
                filename="", mode="wb", fileobj=self.binary, mtime=0,
                compresslevel=GATHER_GZIP_LEVEL,
            )
            self.handle = io.TextIOWrapper(
                self.compressed, encoding="utf-8", newline=""
            )
        else:
            self.handle = self.temporary.open(
                "w", encoding="utf-8", newline=""
            )
        self.writer = csv.DictWriter(
            self.handle,
            delimiter="\t",
            fieldnames=self.fields,
            extrasaction="ignore",
            lineterminator="\n",
        )
        self.writer.writeheader()
        return self

    def writerow(self, row: dict[str, object]) -> None:
        assert self.writer is not None
        self.writer.writerow(row)
        self.rows += 1

    def __exit__(self, exc_type: object, exc: object, traceback: object) -> None:
        assert self.handle is not None
        self.handle.close()
        if self.compressed is not None and not self.compressed.closed:
            self.compressed.close()
        if self.binary is not None and not self.binary.closed:
            self.binary.close()
        if exc_type is None:
            os.replace(self.temporary, self.path)
        else:
            self.temporary.unlink(missing_ok=True)


def validate_library_audit(
    row: dict[str, str],
    scientific_config: dict[str, object],
    *,
    expected_manifest: Path | None = None,
    expected_row_index: int | None = None,
    deep_validate: bool = True,
    validated_profiler_path: Path | None = None,
) -> tuple[Path, dict[str, object], dict[str, dict[str, object]], dict[str, object]]:
    library = row["library"]
    evidence = Path(row["output_dir"])
    nonempty_file(
        evidence / "BAM_EVIDENCE_COMPLETE.ok", f"{library} completion marker"
    )
    audit_path = nonempty_file(evidence / "audit.json", f"{library} audit")
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    if (
        audit.get("release") != RELEASE
        or audit.get("library") != library
        or audit.get("reconciliation", {}).get("status") != "PASS"
    ):
        raise EvidenceError(f"library audit is not a passing match for {library}")
    signature = audit.get("resume_signature")
    expected_parameters = scientific_parameters(row, scientific_config)
    expected_inputs = {
        column: file_signature(Path(row[column]))
        for column in SIGNATURE_INPUT_COLUMNS
        if row.get(column)
    }
    profiler_path = validated_profiler_path or frozen_profiler_path(
        scientific_config
    )
    if (
        not isinstance(signature, dict)
        or signature.get("release") != RELEASE
        or signature.get("manifest_row") != row
        or signature.get("inputs") != expected_inputs
        or signature.get("profiler_executable") != file_signature(profiler_path)
        or signature.get("profiler") != audit.get("profiler_version")
        or signature.get("parameters") != expected_parameters
        or signature.get("scientific_configuration_hash")
        != scientific_config.get("configuration_hash")
        or audit.get("parameters") != expected_parameters
    ):
        raise EvidenceError(
            f"library audit scientific configuration differs from the current "
            f"manifest for {library}"
        )
    if expected_manifest is not None:
        manifest_value = audit.get("manifest")
        if (
            not isinstance(manifest_value, str)
            or Path(manifest_value).expanduser().resolve()
            != expected_manifest.expanduser().resolve()
        ):
            raise EvidenceError(
                f"library audit manifest differs from the current project for {library}"
            )
    if (
        expected_row_index is not None
        and audit.get("manifest_row_index") != expected_row_index
    ):
        raise EvidenceError(
            f"library audit row index differs from the current manifest for {library}"
        )
    audit_config = audit.get("scientific_configuration")
    expected_config_path = Path(row["scientific_config"]).expanduser().resolve()
    if (
        not isinstance(audit_config, dict)
        or Path(str(audit_config.get("path", ""))).expanduser().resolve()
        != expected_config_path
        or audit_config.get("configuration_hash")
        != scientific_config.get("configuration_hash")
    ):
        raise EvidenceError(
            f"library audit has stale scientific-configuration provenance for {library}"
        )
    recorded = audit.get("outputs")
    if not isinstance(recorded, dict) or set(recorded) != set(REQUIRED_PRODUCTS):
        raise EvidenceError(f"library audit output inventory mismatch for {library}")
    if deep_validate:
        actual = {}
        for product_index, name in enumerate(REQUIRED_PRODUCTS, start=1):
            print(
                f"    deep verifying {library} product "
                f"{product_index}/{len(REQUIRED_PRODUCTS) + 1}: {name}",
                flush=True,
            )
            actual[name] = validate_product(evidence / name, name)
        for name, stats in actual.items():
            ensure_stats_match(stats, recorded[name], f"{library}/{name}")
        print(
            f"    deep verifying {library} product "
            f"{len(REQUIRED_PRODUCTS) + 1}/{len(REQUIRED_PRODUCTS) + 1}: "
            "CellReads.countedU.from_bam.tsv.gz",
            flush=True,
        )
        compatibility = validate_compatibility(
            evidence / "CellReads.countedU.from_bam.tsv.gz"
        )
        ensure_stats_match(
            compatibility,
            audit.get("compatibility_product", {}),
            f"{library}/CellReads.countedU.from_bam.tsv.gz",
        )
    else:
        actual = {
            name: validate_recorded_product_metadata(
                evidence / name,
                recorded[name],
                SCHEMAS[name],
                f"{library}/{name}",
            )
            for name in REQUIRED_PRODUCTS
        }
        compatibility = validate_recorded_product_metadata(
            evidence / "CellReads.countedU.from_bam.tsv.gz",
            audit.get("compatibility_product"),
            {"CB", "countedU"},
            f"{library}/CellReads.countedU.from_bam.tsv.gz",
            require_rows=True,
        )
    return evidence, audit, actual, compatibility


def validate_gather_output(
    path: Path,
    name: str,
) -> dict[str, object]:
    nonempty_file(path, name)
    digest = hashlib.sha256()
    previous: tuple[object, ...] | None = None
    rows = 0
    with open_hashed_text(path, digest) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        missing = GATHER_SCHEMAS[name] - set(fields)
        if missing:
            raise EvidenceError(
                f"{path} lacks gathered columns: {', '.join(sorted(missing))}"
            )
        for line_number, row in enumerate(reader, start=2):
            raw = [row[field] for field in GATHER_KEYS[name]]
            if not all(raw):
                raise EvidenceError(f"blank gather key in {path}:{line_number}")
            key: tuple[object, ...]
            if name in {
                "all_libraries_rg_summary.tsv.gz",
                "all_libraries_source_yield.tsv",
            }:
                index = int(raw[1])
                key = (
                    raw[0], index if index >= 0 else 1 << 30, *raw[2:]
                )
            else:
                key = tuple(raw)
            if previous is not None and key <= previous:
                raise EvidenceError(
                    f"duplicate or out-of-order gather key in {path}:{line_number}"
                )
            previous = key
            rows += 1
    return {
        "rows": rows,
        "bytes": path.stat().st_size,
        "sha256": digest.hexdigest(),
        "schema": fields,
    }


def gather_resume_is_valid(
    output: Path,
    signature: dict[str, object],
    *,
    deep_verify: bool = False,
) -> bool:
    marker = output / "BAM_EVIDENCE_GATHER_COMPLETE.ok"
    audit_path = output / "project_audit.json"
    if not marker.is_file() or not audit_path.is_file():
        return False
    try:
        audit = json.loads(audit_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise EvidenceError(f"could not read gather audit: {exc}") from exc
    if audit.get("resume_signature") != signature:
        return False
    recorded = audit.get("aggregate_outputs")
    if not isinstance(recorded, dict) or set(recorded) != set(GATHER_SCHEMAS):
        raise EvidenceError("gather audit has an incomplete aggregate inventory")
    for name in GATHER_SCHEMAS:
        if deep_verify:
            actual = validate_gather_output(output / name, name)
            ensure_stats_match(actual, recorded[name], f"gather/{name}")
        else:
            validate_recorded_product_metadata(
                output / name,
                recorded[name],
                GATHER_SCHEMAS[name],
                f"gather/{name}",
            )
    return True


def gather(args: argparse.Namespace) -> None:
    """BAM-free, streaming project gather with fast audited resume."""
    started = time.monotonic()
    deep_verify = bool(getattr(args, "deep_verify", False))
    manifest = Path(args.manifest).expanduser().resolve()
    fields, rows = read_manifest(manifest)
    scientific_config = load_project_scientific_config(manifest, fields, rows)
    profiler_path = frozen_profiler_path(scientific_config)
    output = Path(args.output_dir).expanduser().resolve(strict=False)
    mode = "deep content verification" if deep_verify else "fast audited verification"
    print(
        f"RNA BAM evidence gather: {len(rows)} libraries; {mode}",
        flush=True,
    )
    signature = {
        "release": RELEASE,
        "manifest": str(manifest),
        "manifest_sha256": sha256(manifest),
        "scientific_configuration_hash": scientific_config["configuration_hash"],
        "libraries": [row["library"] for row in rows],
        "library_audit_sha256": {
            row["library"]: (
                sha256(Path(row["output_dir"]) / "audit.json")
                if (Path(row["output_dir"]) / "audit.json").is_file()
                else "MISSING"
            )
            for row in rows
        },
    }
    marker = output / "BAM_EVIDENCE_GATHER_COMPLETE.ok"
    audit_path = output / "project_audit.json"
    validated_libraries: dict[
        str,
        tuple[Path, dict[str, object], dict[str, dict[str, object]], dict[str, object]],
    ] = {}
    # Normal resume checks immutable scientific provenance plus the complete
    # recorded inventory and current file sizes. --deep-verify additionally
    # streams every source and aggregate byte through schema/order/hash checks.
    if marker.is_file() and audit_path.is_file():
        for row_index, manifest_row in enumerate(rows):
            library = manifest_row["library"]
            validated_libraries[library] = validate_library_audit(
                manifest_row,
                scientific_config,
                expected_manifest=manifest,
                expected_row_index=row_index,
                deep_validate=deep_verify,
                validated_profiler_path=profiler_path,
            )
            print(
                f"  library audit {row_index + 1}/{len(rows)}: {library} PASS",
                flush=True,
            )
    try:
        resumed = gather_resume_is_valid(
            output, signature, deep_verify=deep_verify
        )
    except EvidenceError:
        if not args.replace_stale:
            raise
        resumed = False
    if resumed:
        print(
            f"Revalidated existing RNA BAM evidence gather in "
            f"{time.monotonic() - started:.1f}s: {output}",
            flush=True,
        )
        return
    if marker.exists() and not args.replace_stale:
        raise EvidenceError(
            f"stale or corrupt completed gather at {output}; "
            "use --replace-stale only after review"
        )
    output.mkdir(parents=True, exist_ok=True)
    marker.unlink(missing_ok=True)
    audit_path.unlink(missing_ok=True)

    barcode_fields = [
        "library", "CB", "current_raw_member", "current_filtered_member",
        "historical_raw_member", "historical_filtered_member",
        "barcode_category", "primary_mapped_reads",
        "candidate_countedU_reads", "nh_gt1_unique_gene_countedU_reads",
        "candidate_matrix_molecules", "biological_classified_reads",
        "mitochondrial_reads", "rrna_reads",
        "biological_classification_status",
    ]
    rg_fields: list[str] | None = None
    # Read the first small RG header to establish the project output schema.
    first_evidence = Path(rows[0]["output_dir"])
    with open_text(first_evidence / "rg_summary.tsv") as handle:
        rg_fields = list(csv.DictReader(handle, delimiter="\t").fieldnames or [])
    if not rg_fields:
        raise EvidenceError("first RG summary has no schema")
    source_fields = [
        "library", "source_id", "source_order_index", "RGs", "rg_count",
        *CORE_ADDITIVE_FIELDS, "biological_classification_status",
        "candidate_matrix_molecules",
        "marginal_candidate_molecules_first_observed",
    ]
    inventory_fields = [
        "library", "product", "path", "rows", "bytes", "sha256"
    ]
    library_audits: dict[str, object] = {}

    with (
        AtomicTsvWriter(
            output / "all_libraries_barcode_summary.tsv.gz",
            barcode_fields,
            compressed=True,
        ) as barcode_writer,
        AtomicTsvWriter(
            output / "all_libraries_rg_summary.tsv.gz",
            rg_fields,
            compressed=True,
        ) as rg_writer,
        AtomicTsvWriter(
            output / "all_libraries_source_yield.tsv", source_fields
        ) as source_writer,
        AtomicTsvWriter(
            output / "bam_evidence_inventory.tsv", inventory_fields
        ) as inventory_writer,
    ):
        row_indexes = {
            manifest_row["library"]: row_index
            for row_index, manifest_row in enumerate(rows)
        }
        sorted_rows = sorted(rows, key=lambda item: item["library"])
        for gather_index, manifest_row in enumerate(sorted_rows, start=1):
            library = manifest_row["library"]
            print(
                f"  aggregating library {gather_index}/{len(sorted_rows)}: "
                f"{library}",
                flush=True,
            )
            if library in validated_libraries:
                evidence, audit, products, compatibility = (
                    validated_libraries[library]
                )
            else:
                evidence, audit, products, compatibility = validate_library_audit(
                    manifest_row,
                    scientific_config,
                    expected_manifest=manifest,
                    expected_row_index=row_indexes[library],
                    deep_validate=deep_verify,
                    validated_profiler_path=profiler_path,
                )
            library_audits[library] = {
                "audit_sha256": sha256(evidence / "audit.json"),
                "resume_signature": audit.get("resume_signature"),
            }
            inventory_products: dict[str, dict[str, object]] = {
                **products,
                "CellReads.countedU.from_bam.tsv.gz": compatibility,
            }
            for name in sorted(inventory_products):
                inventory_writer.writerow(
                    {
                        "library": library,
                        "product": name,
                        "path": evidence / name,
                        **inventory_products[name],
                    }
                )

            with open_text(evidence / "rg_summary.tsv") as handle:
                rg_reader = csv.DictReader(handle, delimiter="\t")
                library_rg_fields = list(rg_reader.fieldnames or [])
                small_rg_rows = list(rg_reader)
            if library_rg_fields != rg_fields:
                raise EvidenceError(
                    f"RG summary schema differs across libraries at {library}"
                )
            additive = {
                field: sum(integer_field(item, field) for item in small_rg_rows)
                for field in CORE_ADDITIVE_FIELDS
            }
            profiler = audit.get("profiler_summary", {})
            for field, profiler_field in (
                ("all_records", "total_records"),
                ("primary_mapped_reads", "primary_mapped_reads"),
                ("candidate_countedU_reads", "candidate_countedU_reads"),
                ("nh_gt1_unique_gene_countedU_reads",
                 "nh_gt1_unique_gene_countedU_reads"),
            ):
                if additive[field] != int(profiler.get(profiler_field, -1)):
                    raise EvidenceError(
                        f"RG {field} does not reconcile for {library}"
                    )
            for item in small_rg_rows:
                rg_writer.writerow(item)

            source_by_index = {
                int(item["source_index"]): item["source_id"]
                for item in small_rg_rows
                if int(item["source_index"]) >= 0
            }
            source_candidate = {value: 0 for value in source_by_index.values()}
            source_marginal = {value: 0 for value in source_by_index.values()}
            mask_sources: dict[int, tuple[str, ...]] = {}
            marginal_total = 0
            for molecule in iter_tsv(
                evidence / "molecule_source_hash_bins.tsv.gz"
            ):
                mask = int(molecule["source_mask_hex"], 16)
                count = int(molecule["candidate_matrix_molecules"])
                first_index = (mask & -mask).bit_length() - 1
                if first_index not in source_by_index:
                    raise EvidenceError(
                        f"molecule source mask uses unknown bit for {library}"
                    )
                source_marginal[source_by_index[first_index]] += count
                marginal_total += count
                contributing_sources = mask_sources.get(mask)
                if contributing_sources is None:
                    contributing_sources = tuple(
                        source
                        for index, source in source_by_index.items()
                        if mask & (1 << index)
                    )
                    mask_sources[mask] = contributing_sources
                for source in contributing_sources:
                    source_candidate[source] += count
            if marginal_total != int(
                profiler.get("candidate_matrix_molecules", -1)
            ):
                raise EvidenceError(
                    f"marginal molecule total does not reconcile for {library}"
                )

            by_source: dict[str, dict[str, object]] = {}
            for item in small_rg_rows:
                if int(item["source_index"]) < 0:
                    continue
                source = item["source_id"]
                target = by_source.setdefault(
                    source,
                    {
                        "library": library,
                        "source_id": source,
                        "source_order_index": int(item["source_index"]),
                        "RGs": [],
                        **{field: 0 for field in CORE_ADDITIVE_FIELDS},
                    },
                )
                target["RGs"].append(item["RG"])
                for field in CORE_ADDITIVE_FIELDS:
                    target[field] = int(target[field]) + integer_field(item, field)
            for source, target in sorted(
                by_source.items(),
                key=lambda value: (
                    int(value[1]["source_order_index"]), value[0]
                ),
            ):
                rgs = sorted(target["RGs"])
                target["RGs"] = ",".join(rgs)
                target["rg_count"] = len(rgs)
                classified = int(target["biological_classified_reads"])
                primary_mapped = int(target["primary_mapped_reads"])
                if classified == 0:
                    target["mitochondrial_reads"] = ""
                    target["rrna_reads"] = ""
                    target["biological_classification_status"] = "unavailable"
                else:
                    target["biological_classification_status"] = (
                        "complete" if classified == primary_mapped else "partial"
                    )
                target["candidate_matrix_molecules"] = source_candidate[source]
                target["marginal_candidate_molecules_first_observed"] = (
                    source_marginal[source]
                )
                source_writer.writerow(target)

            for barcode in iter_tsv(
                evidence / "barcode_read_metrics.tsv.gz"
            ):
                barcode_writer.writerow(
                    {
                        "library": library,
                        **{
                            field: barcode.get(field, "")
                            for field in barcode_fields
                            if field != "library"
                        },
                    }
                )
            print(
                f"  aggregated library {gather_index}/{len(sorted_rows)}: "
                f"{library} ({time.monotonic() - started:.1f}s elapsed)",
                flush=True,
            )

    aggregate_stats: dict[str, dict[str, object]] = {}
    for output_index, name in enumerate(GATHER_SCHEMAS, start=1):
        print(
            f"  validating aggregate {output_index}/{len(GATHER_SCHEMAS)}: "
            f"{name}",
            flush=True,
        )
        aggregate_stats[name] = validate_gather_output(output / name, name)
        print(
            f"  validated aggregate {output_index}/{len(GATHER_SCHEMAS)}: "
            f"{name}",
            flush=True,
        )
    project_audit = {
        "release": RELEASE,
        "created_utc": utc_now(),
        "resume_signature": signature,
        "library_count": len(rows),
        "libraries": library_audits,
        "validation": {
            "all_completion_markers_present": True,
            "schemas_row_counts_hashes_and_adjacent_keys_validated": deep_verify,
            "source_product_validation": (
                "deep_schema_rows_hashes_and_adjacent_keys"
                if deep_verify
                else "publication_audits_plus_current_inventory_and_sizes"
            ),
            "aggregate_schemas_row_counts_hashes_and_adjacent_keys_validated": True,
            "additive_rg_read_totals_reconciled": True,
            "marginal_first_source_molecules_reconciled": True,
            "resume_revalidates_every_declared_output": True,
            "bam_reopened": False,
        },
        "aggregate_outputs": aggregate_stats,
    }
    atomic_json(audit_path, project_audit)
    atomic_text(marker, utc_now() + "\n")
    print(
        f"Gathered validated RNA BAM evidence for {len(rows)} libraries in "
        f"{time.monotonic() - started:.1f}s: {output}",
        flush=True,
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", action="version", version=f"%(prog)s {RELEASE}")
    subparsers = parser.add_subparsers(dest="command", required=True)

    run = subparsers.add_parser(
        "run", help="profile and validate exactly one manifest row"
    )
    run.add_argument("--manifest", required=True)
    run.add_argument("--row-index", required=True, type=int)
    run.add_argument("--profiler", default="rna_bam_evidence")
    run.add_argument("--threads", type=int, default=4)
    run.add_argument("--hash-bins", type=int, default=100)
    run.add_argument("--hash-seed", type=int, default=DEFAULT_HASH_SEED)
    run.add_argument(
        "--max-memory-gb",
        type=float,
        default=40.0,
        help="hard total-process address-space ceiling for the compiled profiler",
    )
    run.add_argument("--replace-stale", action="store_true")

    gather_parser = subparsers.add_parser(
        "gather",
        help="validate and combine completed libraries without opening BAMs",
    )
    gather_parser.add_argument("--manifest", required=True)
    gather_parser.add_argument("--output-dir", required=True)
    gather_parser.add_argument("--replace-stale", action="store_true")
    gather_parser.add_argument(
        "--deep-verify",
        action="store_true",
        help=(
            "reread, decompress, validate, and checksum every per-library and "
            "aggregate product; ordinary gather uses audited metadata on resume"
        ),
    )

    return parser


def main() -> int:
    args = build_parser().parse_args()
    try:
        if args.command == "run":
            if args.threads < 1 or args.hash_bins < 1 or args.max_memory_gb <= 0:
                raise EvidenceError(
                    "threads, hash bins, and maximum memory must be positive"
                )
            run_one(args)
        else:
            gather(args)
        return 0
    except (EvidenceError, KeyError, ValueError, json.JSONDecodeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
