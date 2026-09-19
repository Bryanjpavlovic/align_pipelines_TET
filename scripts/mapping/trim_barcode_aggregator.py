#!/usr/bin/env python3
"""Join per-read cutadapt records to STARsolo-corrected 10X barcodes.

The master outputs contain every observed barcode assigned by STARsolo or the
conservative exact/unique one-mismatch fallback. STARsolo's filtered cell set
is an annotation and a separate derived view; it is never an inclusion gate
for the master data. An empty-drop view is produced only when an explicit
library/barcode roster is supplied. A barcode not called as a STARsolo cell is
not automatically labelled as an empty drop.

When a profiler ``raw_to_corrected_barcode_counts.tsv.gz`` bridge is present,
its observed ``CR`` -> final ``CB`` mapping is authoritative and no BAM is
opened. The manifest's exact RG resolves differences between source FASTQs;
conflicts within one RG are reported and left unresolved, never guessed. The
historical BAM lookup remains a fallback for runs without the bridge. Reads
whose raw barcode is absent from either evidence source use the conservative
exact/unique one-mismatch whitelist fallback.

Source-level tables retain the BAM read-group ID and exact source FASTQ for
each aggregate. This preserves source provenance without a prohibitively large
per-read output.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import subprocess
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from itertools import zip_longest
from pathlib import Path
from typing import Callable, Iterator, TextIO


RELEASE = "2026-09-05-v9-profiler-barcode-bridge"
FRONT_ADAPTERS = {"TSO", "TSO_5prime_marker", "TSO_5prime_marker_RC"}
MANIFEST_REQUIRED = {
    "library",
    "run_id",
    "read_format",
    "mate",
    "barcode_fastq",
    "info_file",
    "tso_info_file",
}
MANIFEST_OPTIONAL = {
    "bam",
    "barcode_correction_bridge",
    "rg_metadata",
    "source_id",
    "source_fastq",
    "raw_barcodes",
    "filtered_barcodes",
}


class AggregatorError(RuntimeError):
    """A data-integrity or configuration error."""


@dataclass
class InfoRecord:
    read_name: str
    original_length: int
    final_length: int
    adapter_counts: Counter[str] = field(default_factory=Counter)


@dataclass
class BarcodeStats:
    total_reads: int = 0
    barcode_starsolo_cb_reads: int = 0
    barcode_exact_reads: int = 0
    barcode_unique_1mm_reads: int = 0
    reads_with_adapter: int = 0
    reads_no_adapter: int = 0
    reads_too_short: int = 0
    total_bp_trimmed: int = 0
    original_length_sum: int = 0
    final_length_sum: int = 0
    multi_adapter_reads: int = 0
    fixed_tso_trimmed_reads: int = 0
    tso_unrecognized_reads: int = 0
    adapter_counts: Counter[str] = field(default_factory=Counter)
    source_ids: set[str] = field(default_factory=set)
    source_runs: set[str] = field(default_factory=set)
    source_fastqs: set[str] = field(default_factory=set)


@dataclass(frozen=True)
class BarcodeAnnotation:
    is_starsolo_filtered_cell: bool | None
    is_starsolo_raw_barcode: bool | None
    is_explicit_empty_drop: bool | None

    @property
    def population(self) -> str:
        if self.is_explicit_empty_drop is True:
            return "explicit_empty_drop"
        if self.is_starsolo_filtered_cell is True:
            return "filtered_cell"
        if self.is_starsolo_filtered_cell is False:
            return "observed_not_filtered"
        return "not_evaluated"


@dataclass
class StarsoloCBMap:
    bam: Path
    source_kind: str
    rg_filter: str | None
    assignments: dict[str, str]
    conflicts: set[str]
    conflict_read_assignments: dict[str, str]
    conflict_read_assignment_conflicts: set[str]
    alignments_seen: int
    tagged_alignments: int


def open_text(path: Path, mode: str) -> TextIO:
    return gzip.open(path, mode) if path.name.endswith(".gz") else path.open(mode)


def normalize_barcode(value: str) -> str:
    value = value.strip()
    if not value:
        return ""
    value = value.split()[0]
    if value in {"-", "*", "."}:
        return ""
    return value.rsplit("-", 1)[0] if value.endswith("-1") else value


def load_barcodes(path: Path) -> set[str]:
    if not path.is_file():
        raise AggregatorError(f"barcode file does not exist: {path}")
    barcodes: set[str] = set()
    with open_text(path, "rt") as handle:
        for line in handle:
            barcode = normalize_barcode(line)
            if barcode:
                barcodes.add(barcode)
    return barcodes


def load_starsolo_cb_map(path: Path, whitelist: set[str]) -> StarsoloCBMap:
    """Load STARsolo's observed raw-barcode to corrected-barcode decisions.

    STARsolo writes the raw barcode in ``CR`` and its corrected barcode in
    ``CB``. A compact CR->CB map handles raw barcodes with one observed final
    decision. If STARsolo assigned one CR to multiple CBs, a second BAM pass
    retains the exact read-name-level decision instead of guessing.
    """
    if not path.is_file():
        raise AggregatorError(f"STARsolo BAM does not exist: {path}")
    command = ["samtools", "view", "-@", "1", "-F", "2304", str(path)]
    try:
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            encoding="utf-8",
            errors="replace",
        )
    except OSError as exc:
        raise AggregatorError(f"could not start samtools for {path}: {exc}") from exc
    assert process.stdout is not None
    assignments: dict[str, str] = {}
    conflicts: set[str] = set()
    alignments_seen = 0
    tagged_alignments = 0
    for line in process.stdout:
        alignments_seen += 1
        fields = line.rstrip("\r\n").split("\t")
        raw_barcode: str | None = None
        corrected_barcode: str | None = None
        for field in fields[11:]:
            if field.startswith("CR:Z:"):
                raw_barcode = normalize_barcode(field[5:])
            elif field.startswith("CB:Z:"):
                corrected_barcode = normalize_barcode(field[5:])
        if not raw_barcode or not corrected_barcode:
            continue
        tagged_alignments += 1
        if corrected_barcode not in whitelist:
            process.kill()
            process.wait()
            raise AggregatorError(
                f"STARsolo CB is absent from the configured whitelist for {path}: "
                f"{corrected_barcode}"
            )
        if raw_barcode in conflicts:
            continue
        existing = assignments.get(raw_barcode)
        if existing is None:
            assignments[raw_barcode] = corrected_barcode
        elif existing != corrected_barcode:
            assignments.pop(raw_barcode, None)
            conflicts.add(raw_barcode)
    stderr = process.stderr.read() if process.stderr is not None else ""
    return_code = process.wait()
    if return_code != 0:
        raise AggregatorError(
            f"samtools failed while reading {path} (exit {return_code}): "
            f"{stderr.strip() or '[no diagnostic]'}"
        )
    if not assignments and not conflicts:
        raise AggregatorError(
            f"STARsolo BAM contains no usable CR/CB-tagged alignments: {path}"
        )
    conflict_read_assignments: dict[str, str] = {}
    conflict_read_assignment_conflicts: set[str] = set()
    if conflicts:
        try:
            conflict_process = subprocess.Popen(
                command,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                encoding="utf-8",
                errors="replace",
            )
        except OSError as exc:
            raise AggregatorError(
                f"could not restart samtools for conflicting CR assignments in {path}: {exc}"
            ) from exc
        assert conflict_process.stdout is not None
        for line in conflict_process.stdout:
            fields = line.rstrip("\r\n").split("\t")
            raw_barcode: str | None = None
            corrected_barcode: str | None = None
            for field in fields[11:]:
                if field.startswith("CR:Z:"):
                    raw_barcode = normalize_barcode(field[5:])
                elif field.startswith("CB:Z:"):
                    corrected_barcode = normalize_barcode(field[5:])
            if raw_barcode not in conflicts or not corrected_barcode:
                continue
            read_name = fields[0].split()[0]
            if read_name in conflict_read_assignment_conflicts:
                continue
            existing = conflict_read_assignments.get(read_name)
            if existing is None:
                conflict_read_assignments[read_name] = corrected_barcode
            elif existing != corrected_barcode:
                conflict_read_assignments.pop(read_name, None)
                conflict_read_assignment_conflicts.add(read_name)
        conflict_stderr = (
            conflict_process.stderr.read() if conflict_process.stderr is not None else ""
        )
        conflict_return_code = conflict_process.wait()
        if conflict_return_code != 0:
            raise AggregatorError(
                f"samtools failed during conflict resolution for {path} "
                f"(exit {conflict_return_code}): "
                f"{conflict_stderr.strip() or '[no diagnostic]'}"
            )
    return StarsoloCBMap(
        bam=path,
        source_kind="bam_read_name_fallback",
        rg_filter=None,
        assignments=assignments,
        conflicts=conflicts,
        conflict_read_assignments=conflict_read_assignments,
        conflict_read_assignment_conflicts=conflict_read_assignment_conflicts,
        alignments_seen=alignments_seen,
        tagged_alignments=tagged_alignments,
    )


def load_profiler_cb_bridge(
    path: Path,
    whitelist: set[str],
    rg_filter: str | None = None,
) -> StarsoloCBMap:
    """Load the compact profiler bridge without reopening the BAM.

    When an RG is supplied, the bridge is exact for that source FASTQ and a
    conflict in another RG does not make this source ambiguous. Without an RG
    filter, a raw barcode carrying ``global_conflict=1`` is deliberately
    unresolved. With an exact RG filter, only ``within_rg_conflict`` applies;
    the cross-RG status is audit information. The profiler emits no QNAME table.
    """
    if not path.is_file():
        raise AggregatorError(f"barcode-correction bridge does not exist: {path}")
    assignments: dict[str, str] = {}
    conflicts: set[str] = set()
    reads_seen = 0
    with open_text(path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "RG", "CR", "CB", "read_count", "within_rg_conflict",
            "cross_rg_conflict", "global_conflict",
        }
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise AggregatorError(
                f"barcode-correction bridge lacks {', '.join(sorted(missing))}: {path}"
            )
        for line_number, row in enumerate(reader, start=2):
            if rg_filter is not None and row["RG"] != rg_filter:
                continue
            raw = normalize_barcode(row["CR"])
            corrected = normalize_barcode(row["CB"])
            if not raw or not corrected:
                raise AggregatorError(f"blank CR/CB in bridge {path} line {line_number}")
            if corrected not in whitelist:
                raise AggregatorError(
                    f"profiler CB is absent from the configured whitelist: {corrected}"
                )
            try:
                count = int(row["read_count"])
            except ValueError as exc:
                raise AggregatorError(
                    f"non-integer read_count in bridge {path} line {line_number}"
                ) from exc
            if count < 1:
                raise AggregatorError(
                    f"nonpositive read_count in bridge {path} line {line_number}"
                )
            reads_seen += count
            flag_field = "global_conflict" if rg_filter is None else "within_rg_conflict"
            flagged = row[flag_field].strip().lower() in {"1", "true", "yes"}
            existing = assignments.get(raw)
            if flagged or (existing is not None and existing != corrected):
                assignments.pop(raw, None)
                conflicts.add(raw)
            elif raw not in conflicts:
                assignments[raw] = corrected
    if rg_filter is None and not assignments and not conflicts:
        raise AggregatorError(f"barcode-correction bridge has no assignments: {path}")
    return StarsoloCBMap(
        bam=path,
        source_kind="profiler_aggregate_bridge",
        rg_filter=rg_filter,
        assignments=assignments,
        conflicts=conflicts,
        conflict_read_assignments={},
        conflict_read_assignment_conflicts=set(),
        alignments_seen=reads_seen,
        tagged_alignments=reads_seen,
    )


def resolve_barcode(
    read_name: str,
    raw_barcode: str,
    whitelist: set[str],
    starsolo_map: StarsoloCBMap | None,
) -> tuple[str | None, str]:
    """Prefer STARsolo's CB; use a labelled conservative fallback if absent."""
    if starsolo_map is not None:
        if raw_barcode in starsolo_map.conflicts:
            if read_name in starsolo_map.conflict_read_assignment_conflicts:
                return None, "starsolo_conflict_unresolved"
            corrected = starsolo_map.conflict_read_assignments.get(read_name)
            if corrected is not None:
                return corrected, "starsolo_cb_read"
            return None, "starsolo_conflict_unresolved"
        corrected = starsolo_map.assignments.get(raw_barcode)
        if corrected is not None:
            return corrected, "starsolo_cb"
        corrected, status = assign_barcode(raw_barcode, whitelist)
        if status in {"exact", "unique_1mm"}:
            status = f"fallback_{status}"
        return corrected, status
    return assign_barcode(raw_barcode, whitelist)


def load_empty_drop_roster(path: Path | None) -> dict[str, set[str]] | None:
    """Load an explicit two-column library/barcode roster."""
    if path is None:
        return None
    if not path.is_file():
        raise AggregatorError(f"empty-drop roster does not exist: {path}")
    with open_text(path, "rt") as handle:
        header = handle.readline().rstrip("\r\n").split("\t")
        normalized = {name.strip().lower(): index for index, name in enumerate(header)}
        library_index = normalized.get("library")
        barcode_index = next(
            (
                normalized[name]
                for name in ("cell_barcode", "barcode", "cb")
                if name in normalized
            ),
            None,
        )
        if library_index is None or barcode_index is None:
            raise AggregatorError(
                "empty-drop roster must be tab-separated with a library column and "
                "one of: cell_barcode, barcode, CB"
            )
        roster: dict[str, set[str]] = defaultdict(set)
        for number, line in enumerate(handle, start=2):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if max(library_index, barcode_index) >= len(fields):
                raise AggregatorError(f"empty-drop roster row {number} is incomplete")
            library = fields[library_index].strip()
            barcode = normalize_barcode(fields[barcode_index])
            if not library or not barcode:
                raise AggregatorError(f"empty-drop roster row {number} has an empty value")
            roster[library].add(barcode)
    return dict(roster)


def assign_barcode(raw: str, whitelist: set[str]) -> tuple[str | None, str]:
    """Assign exact or unique Hamming-1 matches; never guess ambiguity."""
    if raw in whitelist:
        return raw, "exact"
    if len(raw) != 16 or any(base not in "ACGT" for base in raw):
        return None, "uncorrectable"
    candidates: set[str] = set()
    for position, observed in enumerate(raw):
        for replacement in "ACGT":
            if replacement == observed:
                continue
            candidate = raw[:position] + replacement + raw[position + 1 :]
            if candidate in whitelist:
                candidates.add(candidate)
                if len(candidates) > 1:
                    return None, "ambiguous"
    if len(candidates) == 1:
        return next(iter(candidates)), "unique_1mm"
    return None, "uncorrectable"


def stream_r1(path: Path) -> Iterator[tuple[str, str]]:
    with open_text(path, "rt") as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            sequence = handle.readline().rstrip("\r\n")
            plus = handle.readline()
            quality = handle.readline().rstrip("\r\n")
            if not header.startswith("@") or not plus.startswith("+"):
                raise AggregatorError(f"invalid FASTQ record in {path}")
            if len(sequence) != len(quality):
                raise AggregatorError(f"sequence/quality mismatch in {path}")
            yield header[1:].split()[0], sequence[:16]


def parse_info_group(lines: list[list[str]]) -> InfoRecord:
    first = lines[0]
    read_name = first[0].split()[0]
    if len(first) >= 3 and first[1] == "-1":
        original = len(first[2])
        return InfoRecord(read_name, original, original)
    if len(first) < 8:
        raise AggregatorError(f"malformed cutadapt info row for {read_name}")
    original = len(first[4]) + len(first[5]) + len(first[6])
    previous_observed = original
    previous_retained_options: set[int] | None = None
    final_before = 0
    final_after = 0
    final_adapter = ""
    adapters: Counter[str] = Counter()
    for fields in lines:
        if len(fields) >= 2 and fields[1] == "-1":
            continue
        if len(fields) < 8:
            raise AggregatorError(f"malformed cutadapt info row for {read_name}")
        before = len(fields[4])
        matched = len(fields[5])
        after = len(fields[6])
        observed_current = before + matched + after
        # Cutadapt emits one row per --times round.  The total represented by
        # each row is the actual input length for that round, so validate that
        # progression directly.  Do not first guess which side Cutadapt kept
        # from the adapter name: older drivers produced legacy names such as
        # ``^TSO``, and making validation depend on that spelling caused valid
        # 290 -> 260 -> 92 progressions to be rejected.
        if observed_current > previous_observed:
            raise AggregatorError(
                f"cutadapt info length increased within --times rounds for {read_name}"
            )
        if (
            previous_retained_options is not None
            and observed_current not in previous_retained_options
        ):
            expected = "/".join(str(value) for value in sorted(previous_retained_options))
            raise AggregatorError(
                f"cutadapt info rounds are not contiguous for {read_name}: "
                f"observed {observed_current} bp after a row that could retain "
                f"{expected} bp"
            )
        # Older trimming drivers placed Cutadapt's anchor marker before the
        # adapter name (``^TSO=...``) rather than before the named sequence
        # (``TSO=^...``). Cutadapt consequently recorded ``^TSO`` in column 8.
        # Normalize that legacy spelling so its 5' trimming direction and
        # adapter category remain correct for already-generated info files.
        adapter = fields[7].split(";", 1)[0].strip().lstrip("^")
        adapters[adapter] += 1
        previous_observed = observed_current
        previous_retained_options = {before, after}
        final_before = before
        final_after = after
        final_adapter = adapter
    final_length = final_after if final_adapter in FRONT_ADAPTERS else final_before
    return InfoRecord(read_name, original, final_length, adapters)


def stream_info(path: Path) -> Iterator[InfoRecord]:
    with open_text(path, "rt") as handle:
        current_name: str | None = None
        rows: list[list[str]] = []
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\r\n").split("\t")
            name = fields[0].split()[0]
            if current_name is not None and name != current_name:
                yield parse_info_group(rows)
                rows = []
            current_name = name
            rows.append(fields)
        if rows:
            yield parse_info_group(rows)


def inspect_info_read(path: Path, requested_read_name: str) -> None:
    """Parse and report one named read using the production parser."""
    requested = requested_read_name.split()[0]
    if not path.is_file():
        raise AggregatorError(f"cutadapt info file does not exist: {path}")
    for record in stream_info(path):
        if record.read_name != requested:
            continue
        adapters = ",".join(
            f"{name}:{count}" for name, count in sorted(record.adapter_counts.items())
        )
        print(f"read_name\t{record.read_name}")
        print(f"original_length\t{record.original_length}")
        print(f"final_length\t{record.final_length}")
        print(f"adapters\t{adapters or 'none'}")
        return
    raise AggregatorError(f"read is absent from cutadapt info file {path}: {requested}")


def stream_tso_info(path: Path | None) -> Iterator[tuple[str, int, str]]:
    if path is None:
        return
    with gzip.open(path, "rt") as handle:
        header = handle.readline().rstrip("\r\n").split("\t")
        expected = ["read_name", "raw_length", "tso_start", "tso_length", "errors", "status"]
        if header != expected:
            raise AggregatorError(f"unexpected PE150 TSO audit header in {path}")
        for line in handle:
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 6:
                raise AggregatorError(f"malformed PE150 TSO audit row in {path}")
            yield fields[0], int(fields[3]), fields[5]


def read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8") as handle:
        header = handle.readline().rstrip("\r\n").split("\t")
        if len(header) != len(set(header)):
            raise AggregatorError("manifest contains duplicate column names")
        missing = MANIFEST_REQUIRED - set(header)
        unknown = set(header) - MANIFEST_REQUIRED - MANIFEST_OPTIONAL
        if missing or unknown:
            details: list[str] = []
            if missing:
                details.append("missing " + ", ".join(sorted(missing)))
            if unknown:
                details.append("unknown " + ", ".join(sorted(unknown)))
            raise AggregatorError("invalid manifest columns: " + "; ".join(details))
        rows: list[dict[str, str]] = []
        for number, line in enumerate(handle, start=2):
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != len(header):
                raise AggregatorError(f"manifest row {number} has {len(fields)} fields")
            row = dict(zip(header, fields))
            for name in MANIFEST_OPTIONAL:
                row.setdefault(name, "")
            rows.append(row)
    return rows


def validate_rg_metadata(rows: list[dict[str, str]]) -> list[str]:
    """Verify planned source IDs/FASTQs against mapping-time RG metadata."""
    metadata_paths = sorted(
        {Path(row["rg_metadata"]).resolve() for row in rows if row["rg_metadata"]}
    )
    if not metadata_paths:
        return []
    validated: list[str] = []
    for path in metadata_paths:
        if not path.is_file():
            raise AggregatorError(f"RNA RG metadata does not exist: {path}")
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            required = {"raw_r1", "raw_r2", "rg_id", "bp_id"}
            missing = required - set(reader.fieldnames or [])
            if missing:
                raise AggregatorError(
                    f"RG metadata is missing column(s) {', '.join(sorted(missing))}: {path}"
                )
            records: dict[Path, dict[str, str]] = {}
            for record in reader:
                raw_r1 = Path(record["raw_r1"]).resolve(strict=False)
                if raw_r1 in records:
                    raise AggregatorError(f"duplicate raw R1 in RG metadata: {raw_r1}")
                records[raw_r1] = record
        relevant = [
            row
            for row in rows
            if row["rg_metadata"] and Path(row["rg_metadata"]).resolve() == path
        ]
        for row in relevant:
            raw_r1 = Path(row["barcode_fastq"]).resolve(strict=False)
            record = records.get(raw_r1)
            if record is None:
                raise AggregatorError(
                    f"barcode FASTQ has no mapping-time RG metadata record: {raw_r1}"
                )
            expected = {
                "source_id": record["rg_id"],
                "run_id": record["bp_id"],
            }
            for field, value in expected.items():
                if row[field] != value:
                    raise AggregatorError(
                        f"manifest/RG metadata mismatch for {raw_r1}: "
                        f"{field}={row[field]} but expected {value}"
                    )
            if Path(row["source_fastq"]).resolve(strict=False) != Path(
                record["raw_r2"]
            ).resolve(strict=False):
                raise AggregatorError(
                    f"manifest/RG metadata source R2 mismatch for {raw_r1}: "
                    f"{row['source_fastq']} != {record['raw_r2']}"
                )
        validated.append(str(path))
    return validated


def atomic_gzip_writer(path: Path):
    temp = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    return temp, gzip.open(temp, "wt")


def tri_state(value: bool | None) -> str:
    if value is None:
        return "NA"
    return "1" if value else "0"


def annotation_values(annotation: BarcodeAnnotation) -> list[str]:
    return [
        tri_state(annotation.is_starsolo_filtered_cell),
        tri_state(annotation.is_starsolo_raw_barcode),
        tri_state(annotation.is_explicit_empty_drop),
        annotation.population,
    ]


ANNOTATION_HEADER = [
    "is_starsolo_filtered_cell",
    "is_starsolo_raw_barcode",
    "is_explicit_empty_drop",
    "barcode_population",
]
METRIC_HEADER = [
    "total_reads",
    "barcode_starsolo_cb_reads",
    "barcode_exact_reads",
    "barcode_unique_1mm_reads",
    "reads_with_adapter",
    "reads_no_adapter",
    "reads_too_short",
    "frac_reads_with_adapter",
    "total_bp_trimmed",
    "mean_bp_trimmed",
    "multi_adapter_reads",
    "mean_original_length",
    "mean_final_length",
    "fixed_tso_trimmed_reads",
    "tso_unrecognized_reads",
]


def metric_values(value: BarcodeStats) -> list[str]:
    total = value.total_reads
    return [
        str(total),
        str(value.barcode_starsolo_cb_reads),
        str(value.barcode_exact_reads),
        str(value.barcode_unique_1mm_reads),
        str(value.reads_with_adapter),
        str(value.reads_no_adapter),
        str(value.reads_too_short),
        f"{value.reads_with_adapter / total:.6f}",
        str(value.total_bp_trimmed),
        f"{value.total_bp_trimmed / total:.3f}",
        str(value.multi_adapter_reads),
        f"{value.original_length_sum / total:.3f}",
        f"{value.final_length_sum / total:.3f}",
        str(value.fixed_tso_trimmed_reads),
        str(value.tso_unrecognized_reads),
    ]


def merge_annotation(
    existing: BarcodeAnnotation | None,
    incoming: BarcodeAnnotation,
    key: tuple[str, str, str],
) -> BarcodeAnnotation:
    if existing is None:
        return incoming

    def combine(label: str, left: bool | None, right: bool | None) -> bool | None:
        if left is not None and right is not None and left != right:
            raise AggregatorError(
                f"inconsistent {label} annotation across sources for "
                f"barcode/read-format/mate {key}"
            )
        return left if left is not None else right

    return BarcodeAnnotation(
        combine(
            "STARsolo filtered-cell",
            existing.is_starsolo_filtered_cell,
            incoming.is_starsolo_filtered_cell,
        ),
        combine(
            "STARsolo raw-barcode",
            existing.is_starsolo_raw_barcode,
            incoming.is_starsolo_raw_barcode,
        ),
        combine(
            "explicit empty-drop",
            existing.is_explicit_empty_drop,
            incoming.is_explicit_empty_drop,
        ),
    )


def update_stats(
    value: BarcodeStats,
    info_record: InfoRecord,
    assignment_status: str,
    row: dict[str, str],
    source_id: str,
    source_fastq: str,
    tso_record: tuple[str, int, str] | None,
    min_length: int,
) -> None:
    original = info_record.original_length
    final = info_record.final_length
    adapters = Counter(info_record.adapter_counts)
    if tso_record is not None:
        _, tso_length, tso_status = tso_record
        original += tso_length
        adapters["TSO_5prime_internal_fixed"] += 1
        value.fixed_tso_trimmed_reads += 1
        if tso_status == "unrecognized":
            value.tso_unrecognized_reads += 1
    value.total_reads += 1
    if assignment_status in {"starsolo_cb", "starsolo_cb_read"}:
        value.barcode_starsolo_cb_reads += 1
    elif assignment_status in {"exact", "fallback_exact"}:
        value.barcode_exact_reads += 1
    elif assignment_status in {"unique_1mm", "fallback_unique_1mm"}:
        value.barcode_unique_1mm_reads += 1
    else:
        raise AggregatorError(f"unexpected assigned barcode status: {assignment_status}")
    value.source_ids.add(source_id)
    value.source_runs.add(row["run_id"])
    value.source_fastqs.add(source_fastq)
    value.original_length_sum += original
    value.final_length_sum += final
    value.total_bp_trimmed += max(0, original - final)
    value.adapter_counts.update(adapters)
    n_adapters = sum(adapters.values())
    if n_adapters:
        value.reads_with_adapter += 1
    else:
        value.reads_no_adapter += 1
    if n_adapters > 1:
        value.multi_adapter_reads += 1
    effective_min_length = (
        min_length + 28
        if row["read_format"] == "pe150" and row["mate"] == "R1_cDNA"
        else min_length
    )
    if final < effective_min_length:
        value.reads_too_short += 1


def write_combined_summary(
    path: Path,
    library: str,
    stats: dict[tuple[str, str, str], BarcodeStats],
    annotations: dict[tuple[str, str, str], BarcodeAnnotation],
    include: Callable[[BarcodeAnnotation], bool] | None = None,
) -> int:
    header = [
        "cell_barcode",
        "library",
        "read_format",
        "read_mate",
        "source_count",
        "source_ids",
        "source_runs",
        "source_fastqs",
        *ANNOTATION_HEADER,
        *METRIC_HEADER,
    ]
    written: set[str] = set()
    temp, handle = atomic_gzip_writer(path)
    with handle:
        handle.write("\t".join(header) + "\n")
        for key, value in sorted(stats.items()):
            annotation = annotations[key]
            if include is not None and not include(annotation):
                continue
            barcode, read_format, mate = key
            handle.write(
                "\t".join(
                    [
                        barcode,
                        library,
                        read_format,
                        mate,
                        str(len(value.source_ids)),
                        ";".join(sorted(value.source_ids)),
                        ";".join(sorted(value.source_runs)),
                        ";".join(sorted(value.source_fastqs)),
                        *annotation_values(annotation),
                        *metric_values(value),
                    ]
                )
                + "\n"
            )
            written.add(barcode)
    os.replace(temp, path)
    return len(written)


def write_combined_adapters(
    path: Path,
    library: str,
    stats: dict[tuple[str, str, str], BarcodeStats],
    annotations: dict[tuple[str, str, str], BarcodeAnnotation],
    include: Callable[[BarcodeAnnotation], bool] | None = None,
) -> None:
    header = [
        "cell_barcode",
        "library",
        "read_format",
        "read_mate",
        *ANNOTATION_HEADER,
        "adapter",
        "read_count",
        "fraction_of_reads",
    ]
    temp, handle = atomic_gzip_writer(path)
    with handle:
        handle.write("\t".join(header) + "\n")
        for key, value in sorted(stats.items()):
            annotation = annotations[key]
            if include is not None and not include(annotation):
                continue
            barcode, read_format, mate = key
            for adapter, count in sorted(value.adapter_counts.items()):
                handle.write(
                    "\t".join(
                        [
                            barcode,
                            library,
                            read_format,
                            mate,
                            *annotation_values(annotation),
                            adapter,
                            str(count),
                            f"{count / value.total_reads:.6f}",
                        ]
                    )
                    + "\n"
                )
    os.replace(temp, path)


def process_library(
    manifest: Path,
    library: str,
    whitelist_path: Path,
    output_dir: Path,
    min_length: int,
    empty_drop_roster_path: Path | None,
) -> None:
    rows = [row for row in read_manifest(manifest) if row["library"] == library]
    if not rows:
        raise AggregatorError(f"library {library} is not present in {manifest}")
    validated_rg_metadata = validate_rg_metadata(rows)
    whitelist = load_barcodes(whitelist_path)
    empty_drop_roster = load_empty_drop_roster(empty_drop_roster_path)
    library_empty_drops = (
        empty_drop_roster.get(library, set()) if empty_drop_roster is not None else None
    )
    barcode_cache: dict[Path, set[str]] = {}
    starsolo_cache: dict[Path, StarsoloCBMap] = {}
    bridge_cache: dict[tuple[Path, str], StarsoloCBMap] = {}
    stats: dict[tuple[str, str, str], BarcodeStats] = defaultdict(BarcodeStats)
    annotations: dict[tuple[str, str, str], BarcodeAnnotation] = {}
    qc: Counter[str] = Counter()
    source_qc: list[dict[str, object]] = []
    has_filtered_annotations = any(row["filtered_barcodes"] for row in rows)
    has_raw_annotations = any(row["raw_barcodes"] for row in rows)
    has_starsolo_cb_assignments = any(
        row["barcode_correction_bridge"] or row["bam"] for row in rows
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    source_summary_path = output_dir / f"{library}_trim_by_barcode_by_source.tsv.gz"
    source_adapters_path = output_dir / f"{library}_trim_adapters_by_barcode_by_source.tsv.gz"
    source_summary_temp, source_summary_handle = atomic_gzip_writer(source_summary_path)
    source_adapters_temp, source_adapters_handle = atomic_gzip_writer(source_adapters_path)
    source_summary_header = [
        "cell_barcode",
        "library",
        "read_format",
        "read_mate",
        "source_id",
        "run_id",
        "source_fastq",
        "barcode_fastq",
        *ANNOTATION_HEADER,
        *METRIC_HEADER,
    ]
    source_adapter_header = [
        "cell_barcode",
        "library",
        "read_format",
        "read_mate",
        "source_id",
        "run_id",
        "source_fastq",
        *ANNOTATION_HEADER,
        "adapter",
        "read_count",
        "fraction_of_reads",
    ]

    with source_summary_handle, source_adapters_handle:
        source_summary_handle.write("\t".join(source_summary_header) + "\n")
        source_adapters_handle.write("\t".join(source_adapter_header) + "\n")
        for row in rows:
            barcode_fastq = Path(row["barcode_fastq"])
            info_file = Path(row["info_file"])
            tso_file = Path(row["tso_info_file"]) if row["tso_info_file"] else None
            bam = Path(row["bam"]) if row["bam"] else None
            bridge = (
                Path(row["barcode_correction_bridge"])
                if row["barcode_correction_bridge"]
                else None
            )
            source_fastq = row["source_fastq"] or row["barcode_fastq"]
            source_id = row["source_id"] or f"{row['run_id']}:{Path(source_fastq).name}"
            for path in (barcode_fastq, info_file):
                if not path.is_file():
                    raise AggregatorError(f"required input does not exist: {path}")
            if tso_file is not None and not tso_file.is_file():
                raise AggregatorError(f"required PE150 TSO audit does not exist: {tso_file}")
            starsolo_map: StarsoloCBMap | None = None
            if bridge is not None:
                resolved_bridge = bridge.resolve()
                bridge_key = (resolved_bridge, source_id)
                if bridge_key not in bridge_cache:
                    bridge_cache[bridge_key] = load_profiler_cb_bridge(
                        resolved_bridge, whitelist, rg_filter=source_id
                    )
                starsolo_map = bridge_cache[bridge_key]
            elif bam is not None:
                resolved_bam = bam.resolve()
                if resolved_bam not in starsolo_cache:
                    starsolo_cache[resolved_bam] = load_starsolo_cb_map(
                        resolved_bam, whitelist
                    )
                starsolo_map = starsolo_cache[resolved_bam]

            def optional_set(value: str) -> set[str] | None:
                if not value:
                    return None
                path = Path(value)
                if path not in barcode_cache:
                    barcode_cache[path] = load_barcodes(path)
                return barcode_cache[path]

            filtered = optional_set(row["filtered_barcodes"])
            raw_matrix = optional_set(row["raw_barcodes"])
            if filtered is not None and library_empty_drops is not None:
                overlap = filtered & library_empty_drops
                if overlap:
                    preview = ", ".join(sorted(overlap)[:5])
                    raise AggregatorError(
                        f"explicit empty-drop roster overlaps STARsolo filtered cells for "
                        f"{library}: {len(overlap)} barcode(s), including {preview}"
                    )

            local_stats: dict[str, BarcodeStats] = defaultdict(BarcodeStats)
            local_annotations: dict[str, BarcodeAnnotation] = {}
            local_qc: Counter[str] = Counter()
            tso_stream = stream_tso_info(tso_file) if tso_file else None
            triples = zip_longest(
                stream_r1(barcode_fastq),
                stream_info(info_file),
                tso_stream if tso_stream is not None else [],
                fillvalue=None,
            )
            for fastq_item, info_record, tso_record in triples:
                if fastq_item is None or info_record is None:
                    raise AggregatorError(
                        f"FASTQ and cutadapt info record counts differ: "
                        f"{barcode_fastq}, {info_file}"
                    )
                read_name, raw_barcode = fastq_item
                if read_name != info_record.read_name:
                    raise AggregatorError(
                        f"read-name mismatch: FASTQ {read_name}, info {info_record.read_name}"
                    )
                if tso_file is not None:
                    if tso_record is None or tso_record[0] != read_name:
                        raise AggregatorError(f"TSO audit is not synchronized at read {read_name}")
                elif tso_record is not None:
                    raise AggregatorError("unexpected TSO record for a non-PE150-R1 row")

                qc["input_reads"] += 1
                local_qc["input_reads"] += 1
                corrected, assignment_status = resolve_barcode(
                    read_name, raw_barcode, whitelist, starsolo_map
                )
                qc[f"barcode_{assignment_status}"] += 1
                local_qc[f"barcode_{assignment_status}"] += 1
                if corrected is None:
                    continue
                qc["barcode_assigned"] += 1
                local_qc["barcode_assigned"] += 1
                annotation = BarcodeAnnotation(
                    corrected in filtered if filtered is not None else None,
                    corrected in raw_matrix if raw_matrix is not None else None,
                    (
                        corrected in library_empty_drops
                        if library_empty_drops is not None
                        else None
                    ),
                )
                if annotation.is_starsolo_filtered_cell is True:
                    qc["assigned_reads_filtered_cell"] += 1
                    local_qc["assigned_reads_filtered_cell"] += 1
                elif annotation.is_starsolo_filtered_cell is False:
                    qc["assigned_reads_not_filtered"] += 1
                    local_qc["assigned_reads_not_filtered"] += 1
                else:
                    qc["assigned_reads_filtered_status_unavailable"] += 1
                    local_qc["assigned_reads_filtered_status_unavailable"] += 1
                if annotation.is_explicit_empty_drop is True:
                    qc["assigned_reads_explicit_empty_drop"] += 1
                    local_qc["assigned_reads_explicit_empty_drop"] += 1

                key = (corrected, row["read_format"], row["mate"])
                annotations[key] = merge_annotation(annotations.get(key), annotation, key)
                local_annotations[corrected] = merge_annotation(
                    local_annotations.get(corrected), annotation, key
                )
                update_stats(
                    stats[key],
                    info_record,
                    assignment_status,
                    row,
                    source_id,
                    source_fastq,
                    tso_record,
                    min_length,
                )
                update_stats(
                    local_stats[corrected],
                    info_record,
                    assignment_status,
                    row,
                    source_id,
                    source_fastq,
                    tso_record,
                    min_length,
                )

            for barcode, value in sorted(local_stats.items()):
                annotation = local_annotations[barcode]
                source_summary_handle.write(
                    "\t".join(
                        [
                            barcode,
                            library,
                            row["read_format"],
                            row["mate"],
                            source_id,
                            row["run_id"],
                            source_fastq,
                            row["barcode_fastq"],
                            *annotation_values(annotation),
                            *metric_values(value),
                        ]
                    )
                    + "\n"
                )
                for adapter, count in sorted(value.adapter_counts.items()):
                    source_adapters_handle.write(
                        "\t".join(
                            [
                                barcode,
                                library,
                                row["read_format"],
                                row["mate"],
                                source_id,
                                row["run_id"],
                                source_fastq,
                                *annotation_values(annotation),
                                adapter,
                                str(count),
                                f"{count / value.total_reads:.6f}",
                            ]
                        )
                        + "\n"
                    )
            source_qc.append(
                {
                    "source_id": source_id,
                    "run_id": row["run_id"],
                    "read_format": row["read_format"],
                    "read_mate": row["mate"],
                    "source_fastq": source_fastq,
                    "barcode_fastq": row["barcode_fastq"],
                    "starsolo_bam": str(bam) if bam is not None else None,
                    "barcode_correction_bridge": (
                        str(bridge) if bridge is not None else None
                    ),
                    "barcode_correction_source": (
                        starsolo_map.source_kind if starsolo_map is not None else None
                    ),
                    "counters": dict(sorted(local_qc.items())),
                    "output_barcodes": len(local_stats),
                }
            )
    os.replace(source_summary_temp, source_summary_path)
    os.replace(source_adapters_temp, source_adapters_path)

    master_path = output_dir / f"{library}_trim_by_barcode.tsv.gz"
    adapters_path = output_dir / f"{library}_trim_adapters_by_barcode.tsv.gz"
    output_barcodes = write_combined_summary(master_path, library, stats, annotations)
    write_combined_adapters(adapters_path, library, stats, annotations)

    filtered_output_barcodes: int | None = None
    observed_not_filtered_output_barcodes: int | None = None
    if has_filtered_annotations:
        filtered_output_barcodes = write_combined_summary(
            output_dir / f"{library}_trim_by_barcode_filtered_cells.tsv.gz",
            library,
            stats,
            annotations,
            lambda annotation: annotation.is_starsolo_filtered_cell is True,
        )
        observed_not_filtered_output_barcodes = write_combined_summary(
            output_dir / f"{library}_trim_by_barcode_observed_not_filtered.tsv.gz",
            library,
            stats,
            annotations,
            lambda annotation: annotation.is_starsolo_filtered_cell is False,
        )
        write_combined_adapters(
            output_dir / f"{library}_trim_adapters_by_barcode_observed_not_filtered.tsv.gz",
            library,
            stats,
            annotations,
            lambda annotation: annotation.is_starsolo_filtered_cell is False,
        )
        write_combined_adapters(
            output_dir / f"{library}_trim_adapters_by_barcode_filtered_cells.tsv.gz",
            library,
            stats,
            annotations,
            lambda annotation: annotation.is_starsolo_filtered_cell is True,
        )

    empty_output_barcodes: int | None = None
    if library_empty_drops is not None:
        empty_output_barcodes = write_combined_summary(
            output_dir / f"{library}_trim_by_barcode_empty_drops.tsv.gz",
            library,
            stats,
            annotations,
            lambda annotation: annotation.is_explicit_empty_drop is True,
        )
        write_combined_adapters(
            output_dir / f"{library}_trim_adapters_by_barcode_empty_drops.tsv.gz",
            library,
            stats,
            annotations,
            lambda annotation: annotation.is_explicit_empty_drop is True,
        )

    assigned_statuses = (
        "starsolo_cb",
        "starsolo_cb_read",
        "fallback_exact",
        "fallback_unique_1mm",
        "exact",
        "unique_1mm",
    )
    assigned_total = sum(qc[f"barcode_{name}"] for name in assigned_statuses)
    if assigned_total != qc["barcode_assigned"]:
        raise AggregatorError("assigned barcode accounting did not reconcile")
    statuses = (
        *assigned_statuses,
        "ambiguous",
        "uncorrectable",
        "starsolo_conflict_unresolved",
    )
    if sum(qc[f"barcode_{name}"] for name in statuses) != qc["input_reads"]:
        raise AggregatorError("input barcode accounting did not reconcile")

    observed_barcodes = {key[0] for key in stats}
    qc_payload = {
        "release": RELEASE,
        "library": library,
        "barcode_universe": (
            "all observed barcodes assigned by STARsolo CB or the conservative exact/unique "
            "Hamming-1 fallback; not restricted to STARsolo filtered cells"
        ),
        "barcode_assignment": (
            "The profiler aggregate CR-to-CB bridge is preferred and avoids reopening BAM; "
            "the exact manifest RG resolves cross-source differences, while conflicts "
            "within an RG remain unresolved and are reported. For runs without a bridge, "
            "the historical BAM/read-name resolver remains available. Raw barcodes "
            "absent from evidence use exact or unique Hamming-1 fallback; unresolved, "
            "ambiguous, and uncorrectable reads are counted but not assigned"
        ),
        "starsolo_cb_assignment_available": has_starsolo_cb_assignments,
        "starsolo_cb_maps": [
            {
                "source": str(item.bam),
                "source_kind": item.source_kind,
                "rg_filter": item.rg_filter,
                "raw_barcodes_assigned": len(item.assignments),
                "raw_barcodes_conflicting": len(item.conflicts),
                "conflict_read_assignments": len(item.conflict_read_assignments),
                "conflict_read_assignment_conflicts": len(
                    item.conflict_read_assignment_conflicts
                ),
                "alignments_seen": item.alignments_seen,
                "tagged_alignments": item.tagged_alignments,
            }
            for item in sorted(
                [*bridge_cache.values(), *starsolo_cache.values()],
                key=lambda value: (
                    value.source_kind,
                    str(value.bam),
                    value.rg_filter or "",
                ),
            )
        ],
        "profiler_bridge_preferred_over_bam": True,
        "profiler_bridge_conflicts_are_resolved_only_when_unique_within_rg": True,
        "profiler_bridge_global_conflicts_without_rg_are_unresolved": True,
        "starsolo_filtered_cells_are_annotation_only": True,
        "unfiltered_is_not_empty_drop": True,
        "empty_drop_roster": (
            str(empty_drop_roster_path) if empty_drop_roster_path is not None else None
        ),
        "manifest": str(manifest),
        "validated_rg_metadata": validated_rg_metadata,
        "counters": dict(sorted(qc.items())),
        "source_counters": source_qc,
        "output_barcodes": output_barcodes,
        "filtered_output_barcodes": filtered_output_barcodes,
        "observed_not_filtered_output_barcodes": observed_not_filtered_output_barcodes,
        "explicit_empty_drop_output_barcodes": empty_output_barcodes,
        "explicit_empty_drop_roster_barcodes": (
            len(library_empty_drops) if library_empty_drops is not None else None
        ),
        "explicit_empty_drop_roster_barcodes_observed": (
            len(observed_barcodes & library_empty_drops)
            if library_empty_drops is not None
            else None
        ),
        "starsolo_raw_annotation_available": has_raw_annotations,
        "starsolo_filtered_annotation_available": has_filtered_annotations,
    }
    qc_path = output_dir / f"{library}_trim_by_barcode_qc.json"
    temp_json = qc_path.with_name(f".{qc_path.name}.tmp.{os.getpid()}")
    temp_json.write_text(
        json.dumps(qc_payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    os.replace(temp_json, qc_path)


def merge_gzip_tables(files: list[Path], output: Path) -> None:
    if not files:
        raise AggregatorError(f"no input tables found for {output.name}")
    temp, destination = atomic_gzip_writer(output)
    expected_header: str | None = None
    with destination:
        for path in files:
            with gzip.open(path, "rt") as source:
                header = source.readline()
                if expected_header is None:
                    expected_header = header
                    destination.write(header)
                elif header != expected_header:
                    raise AggregatorError(f"table headers differ while merging {path}")
                for line in source:
                    destination.write(line)
    os.replace(temp, output)


def merge_optional_complete(
    per_library_dir: Path,
    pattern: str,
    output: Path,
    expected_count: int,
) -> None:
    files = sorted(per_library_dir.glob(pattern))
    if not files:
        return
    if len(files) != expected_count:
        raise AggregatorError(
            f"found {len(files)}/{expected_count} expected files matching {pattern}"
        )
    merge_gzip_tables(files, output)


def merge_outputs(per_library_dir: Path, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    summaries = sorted(per_library_dir.glob("*_trim_by_barcode.tsv.gz"))
    adapters = sorted(per_library_dir.glob("*_trim_adapters_by_barcode.tsv.gz"))
    sources = sorted(per_library_dir.glob("*_trim_by_barcode_by_source.tsv.gz"))
    source_adapters = sorted(
        per_library_dir.glob("*_trim_adapters_by_barcode_by_source.tsv.gz")
    )
    if not summaries:
        raise AggregatorError("no per-library trimming-by-barcode summaries were found")
    expected_count = len(summaries)
    for label, files in (
        ("adapter", adapters),
        ("source summary", sources),
        ("source adapter", source_adapters),
    ):
        if len(files) != expected_count:
            raise AggregatorError(
                f"found {len(files)}/{expected_count} per-library {label} files"
            )
    merge_gzip_tables(summaries, output_dir / "all_libraries_trim_by_barcode.tsv.gz")
    merge_gzip_tables(
        adapters, output_dir / "all_libraries_trim_adapters_by_barcode.tsv.gz"
    )
    merge_gzip_tables(
        sources, output_dir / "all_libraries_trim_by_barcode_by_source.tsv.gz"
    )
    merge_gzip_tables(
        source_adapters,
        output_dir / "all_libraries_trim_adapters_by_barcode_by_source.tsv.gz",
    )
    merge_optional_complete(
        per_library_dir,
        "*_trim_by_barcode_filtered_cells.tsv.gz",
        output_dir / "all_libraries_trim_by_barcode_filtered_cells.tsv.gz",
        expected_count,
    )
    merge_optional_complete(
        per_library_dir,
        "*_trim_adapters_by_barcode_filtered_cells.tsv.gz",
        output_dir / "all_libraries_trim_adapters_by_barcode_filtered_cells.tsv.gz",
        expected_count,
    )
    merge_optional_complete(
        per_library_dir,
        "*_trim_by_barcode_observed_not_filtered.tsv.gz",
        output_dir / "all_libraries_trim_by_barcode_observed_not_filtered.tsv.gz",
        expected_count,
    )
    merge_optional_complete(
        per_library_dir,
        "*_trim_adapters_by_barcode_observed_not_filtered.tsv.gz",
        output_dir / "all_libraries_trim_adapters_by_barcode_observed_not_filtered.tsv.gz",
        expected_count,
    )
    merge_optional_complete(
        per_library_dir,
        "*_trim_by_barcode_empty_drops.tsv.gz",
        output_dir / "all_libraries_trim_by_barcode_empty_drops.tsv.gz",
        expected_count,
    )
    merge_optional_complete(
        per_library_dir,
        "*_trim_adapters_by_barcode_empty_drops.tsv.gz",
        output_dir / "all_libraries_trim_adapters_by_barcode_empty_drops.tsv.gz",
        expected_count,
    )
    qc_files = sorted(per_library_dir.glob("*_trim_by_barcode_qc.json"))
    if len(qc_files) != expected_count:
        raise AggregatorError(
            f"found {len(qc_files)}/{expected_count} per-library QC JSON files"
        )
    payload = {
        "release": RELEASE,
        "libraries": [json.loads(path.read_text(encoding="utf-8")) for path in qc_files],
    }
    destination = output_dir / "all_libraries_trim_by_barcode_qc.json"
    temp = destination.with_name(f".{destination.name}.tmp.{os.getpid()}")
    temp.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temp, destination)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Aggregate trimming across the full observed 10X barcode population"
    )
    parser.add_argument("--version", action="version", version=RELEASE)
    subparsers = parser.add_subparsers(dest="command", required=True)
    process = subparsers.add_parser("process")
    process.add_argument("--manifest", required=True)
    process.add_argument("--library", required=True)
    process.add_argument("--whitelist", required=True)
    process.add_argument("--output-dir", required=True)
    process.add_argument("--minimum-length", type=int, default=25)
    process.add_argument(
        "--empty-drop-roster",
        default=None,
        help=(
            "Optional TSV(.gz) with library and barcode/cell_barcode columns. "
            "Empty drops are never inferred from STARsolo filtering."
        ),
    )
    merge = subparsers.add_parser("merge")
    merge.add_argument("--per-library-dir", required=True)
    merge.add_argument("--output-dir", required=True)
    inspect_read = subparsers.add_parser(
        "inspect-read",
        help="validate one cutadapt info record with the production parser",
    )
    inspect_read.add_argument("--info-file", required=True)
    inspect_read.add_argument("--read-name", required=True)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    try:
        if args.command == "process":
            process_library(
                Path(args.manifest),
                args.library,
                Path(args.whitelist),
                Path(args.output_dir),
                args.minimum_length,
                Path(args.empty_drop_roster) if args.empty_drop_roster else None,
            )
        elif args.command == "merge":
            merge_outputs(Path(args.per_library_dir), Path(args.output_dir))
        else:
            inspect_info_read(Path(args.info_file), args.read_name)
        return 0
    except (AggregatorError, OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
