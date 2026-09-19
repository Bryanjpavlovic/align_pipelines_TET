#!/usr/bin/env python3
"""Compare historical/current cells with bounded, evidence-bearing background.

Each library is processed independently through a temporary SQLite shard.
Raw-whitelist strings are streamed and only nonzero matrix columns are retained.
All shared/gained/lost cells are exported. Background detail is a deterministic
bounded sample from barcodes that have nonzero matrix or BAM evidence.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import heapq
import io
import json
import math
import os
import sqlite3
import sys
import tempfile
from collections import defaultdict
from contextlib import ExitStack
from pathlib import Path
from typing import Iterable, Iterator, TextIO

from rna_evidence_common import (
    AnalysisError,
    find_product,
    iter_nonzero_matrix_barcodes,
    iter_tsv,
    library_dirs,
    open_text,
    write_tsv,
)


RELEASE = "2026-09-05-gained-cells-v2"
BACKGROUND_RULE = (
    "neither_current_nor_historical_filtered_and_at_least_one_of:"
    "current_nonzero_raw_matrix,historical_nonzero_raw_matrix,"
    "current_BAM_observed_corrected_barcode,historical_BAM_observed_corrected_barcode"
)
COUNT_BINS = (
    0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000,
    10000, 20000, 50000, 100000, 200000, 500000, 1000000, math.inf,
)
FRACTION_BINS = (
    0.0, 0.001, 0.005, 0.01, 0.025, 0.05, 0.1, 0.2, 0.4, 0.6,
    0.8, 0.9, 0.95, 0.99, 1.0, math.inf,
)
QUANTILES = (0.0, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1.0)
CORE_METRICS = (
    "candidate_countedU_reads",
    "nh_gt1_unique_gene_countedU_reads",
    "primary_mapped_reads",
    "primary_mapping_record_fraction",
    "unique_gene_tagged_reads",
    "unique_gene_tagged_primary_fraction",
    "ordinary_candidate_primary_fraction",
    "current_raw_umis",
    "current_raw_genes",
    "bam_conditional_saturation",
    "mitochondrial_primary_read_fraction",
    "rrna_primary_read_fraction",
    "source_dominance_fraction",
)
TRIM_METRICS = (
    "trim_frac_reads_with_adapter",
    "trim_mean_bp_trimmed",
    "trim_fixed_tso_trimmed_reads",
    "trim_tso_unrecognized_reads",
    "trim_reads_too_short",
    "trim_bam_primary_mapped_per_trim_read_proxy",
    "trim_bam_matrix_conditional_saturation",
    "trim_matrix_umis",
    "trim_matrix_genes",
)


def safe_fraction(numerator: int | None, denominator: int) -> str:
    if numerator is None or denominator == 0:
        return ""
    return f"{numerator / denominator:.8f}"


def conditional_saturation(
    reads: int, molecules: int, authority: str
) -> tuple[str, str]:
    if reads <= 0:
        return "", "unavailable_no_ordinary_candidate_counted_reads"
    if molecules > reads:
        return "", "unavailable_molecule_count_exceeds_read_count"
    return f"{(reads - molecules) / reads:.8f}", authority


def canonical_barcode(value: str) -> str:
    return value[:-2] if value.endswith("-1") else value


def deterministic_hash(*values: str) -> int:
    return int.from_bytes(
        hashlib.sha256("\0".join(values).encode("utf-8")).digest()[:8], "big"
    )


class AtomicWriter:
    def __init__(self, path: Path, fields: list[str], compressed: bool = False):
        self.path = path
        self.fields = fields
        self.compressed_output = compressed
        self.temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
        self.handle: TextIO | None = None
        self.binary = None
        self.compressed = None
        self.writer: csv.DictWriter | None = None
        self.rows = 0

    def __enter__(self) -> "AtomicWriter":
        self.path.parent.mkdir(parents=True, exist_ok=True)
        if self.compressed_output:
            self.binary = self.temporary.open("wb")
            self.compressed = gzip.GzipFile(
                filename="", mode="wb", fileobj=self.binary, mtime=0
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
            fieldnames=self.fields,
            delimiter="\t",
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


class Distribution:
    """Exact fixed-bin histogram plus deterministic bounded quantile sample."""

    def __init__(
        self,
        library: str,
        category: str,
        metric: str,
        sample_size: int,
    ):
        self.library = library
        self.category = category
        self.metric = metric
        self.edges = FRACTION_BINS if (
            "fraction" in metric or "saturation" in metric or metric.startswith("trim_frac")
        ) else COUNT_BINS
        self.bins = [0] * (len(self.edges) - 1)
        self.n = 0
        self.minimum = math.inf
        self.maximum = -math.inf
        self.sample_size = sample_size
        self.heap: list[tuple[int, float]] = []

    def add(self, barcode: str, value: object) -> None:
        try:
            number = float(str(value))
        except (TypeError, ValueError):
            return
        if not math.isfinite(number):
            return
        self.n += 1
        self.minimum = min(self.minimum, number)
        self.maximum = max(self.maximum, number)
        index = len(self.bins) - 1
        for candidate, (left, right) in enumerate(zip(self.edges, self.edges[1:])):
            if left <= number < right or (
                candidate == len(self.bins) - 1 and number == right
            ):
                index = candidate
                break
        self.bins[index] += 1
        score = deterministic_hash(
            self.library, self.category, self.metric, barcode
        )
        item = (-score, number)
        if len(self.heap) < self.sample_size:
            heapq.heappush(self.heap, item)
        elif score < -self.heap[0][0]:
            heapq.heapreplace(self.heap, item)

    def histogram_rows(self) -> Iterator[dict[str, object]]:
        for index, count in enumerate(self.bins):
            right = self.edges[index + 1]
            yield {
                "library": self.library,
                "cell_transition": self.category,
                "metric": self.metric,
                "bin_index": index,
                "bin_left_inclusive": self.edges[index],
                "bin_right_exclusive": "inf" if math.isinf(right) else right,
                "count": count,
                "population_n": self.n,
                "histogram_scope": "all_evidence_bearing_rows_exact",
            }

    def quantile_rows(self) -> Iterator[dict[str, object]]:
        values = sorted(item[1] for item in self.heap)
        for quantile in QUANTILES:
            if values:
                index = min(
                    len(values) - 1,
                    max(0, int(round(quantile * (len(values) - 1)))),
                )
                value: object = values[index]
            else:
                value = ""
            yield {
                "library": self.library,
                "cell_transition": self.category,
                "metric": self.metric,
                "quantile": quantile,
                "value": value,
                "population_n": self.n,
                "sample_n": len(values),
                "method": (
                    "exact" if self.n <= self.sample_size
                    else "deterministic_smallest_sha256_sample"
                ),
            }


def optional_fields(path: Path | None, prefix: str) -> list[str]:
    if path is None:
        return []
    with open_text(path) as handle:
        fields = list(csv.DictReader(handle, delimiter="\t").fieldnames or [])
    return sorted(
        prefix + field
        for field in fields
        if field not in {"library", "CB", "cell_barcode"}
    )


def create_database(path: Path) -> sqlite3.Connection:
    connection = sqlite3.connect(path)
    connection.execute("PRAGMA journal_mode=OFF")
    connection.execute("PRAGMA synchronous=OFF")
    connection.execute("PRAGMA temp_store=FILE")
    connection.executescript(
        """
        CREATE TABLE barcode (
            barcode_key TEXT PRIMARY KEY,
            display_cb TEXT NOT NULL,
            current_filtered INTEGER NOT NULL DEFAULT 0,
            old_filtered INTEGER NOT NULL DEFAULT 0,
            current_umi INTEGER,
            current_genes INTEGER,
            current_rank INTEGER,
            old_umi INTEGER,
            old_genes INTEGER,
            old_rank INTEGER,
            current_bam TEXT,
            old_bam TEXT,
            current_source TEXT,
            old_source TEXT
        ) WITHOUT ROWID;
        CREATE TABLE optional_join (
            kind TEXT NOT NULL,
            barcode_key TEXT NOT NULL,
            payload TEXT NOT NULL,
            PRIMARY KEY (kind, barcode_key)
        ) WITHOUT ROWID;
        """
    )
    return connection


def upsert_filtered(
    connection: sqlite3.Connection, path: Path, column: str
) -> None:
    sql = (
        "INSERT INTO barcode(barcode_key,display_cb,{0}) VALUES(?,?,1) "
        "ON CONFLICT(barcode_key) DO UPDATE SET {0}=1"
    ).format(column)
    batch: list[tuple[str, str]] = []
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            display = line.rstrip("\r\n").split("\t", 1)[0]
            if not display:
                raise AnalysisError(f"blank barcode in {path}:{line_number}")
            batch.append((canonical_barcode(display), display))
            if len(batch) >= 10000:
                connection.executemany(sql, batch)
                batch.clear()
    if batch:
        connection.executemany(sql, batch)


def upsert_matrix(
    connection: sqlite3.Connection,
    directory: Path,
    prefix: str,
) -> None:
    umi_column = f"{prefix}_umi"
    gene_column = f"{prefix}_genes"
    sql = (
        f"INSERT INTO barcode(barcode_key,display_cb,{umi_column},{gene_column}) "
        "VALUES(?,?,?,?) ON CONFLICT(barcode_key) DO UPDATE SET "
        f"{umi_column}=excluded.{umi_column},{gene_column}=excluded.{gene_column}"
    )
    batch: list[tuple[str, str, int, int]] = []
    for display, umi, genes in iter_nonzero_matrix_barcodes(
        find_product(directory, "raw/matrix.mtx.gz"),
        find_product(directory, "raw/barcodes.tsv.gz"),
    ):
        batch.append((canonical_barcode(display), display, umi, genes))
        if len(batch) >= 10000:
            connection.executemany(sql, batch)
            batch.clear()
    if batch:
        connection.executemany(sql, batch)
    rank_column = f"{prefix}_rank"
    updates: list[tuple[int, str]] = []
    cursor = connection.execute(
        f"SELECT barcode_key FROM barcode WHERE {umi_column} IS NOT NULL "
        f"ORDER BY {umi_column} DESC, barcode_key"
    )
    for rank, (barcode_key,) in enumerate(cursor, start=1):
        updates.append((rank, barcode_key))
        if len(updates) >= 10000:
            connection.executemany(
                f"UPDATE barcode SET {rank_column}=? WHERE barcode_key=?",
                updates,
            )
            updates.clear()
    if updates:
        connection.executemany(
            f"UPDATE barcode SET {rank_column}=? WHERE barcode_key=?", updates
        )


def upsert_bam(
    connection: sqlite3.Connection,
    directory: Path,
    prefix: str,
) -> bool:
    path = directory / "bam_evidence" / "barcode_read_metrics.tsv.gz"
    if not path.is_file():
        return False
    column = f"{prefix}_bam"
    sql = (
        f"INSERT INTO barcode(barcode_key,display_cb,{column}) VALUES(?,?,?) "
        f"ON CONFLICT(barcode_key) DO UPDATE SET {column}=excluded.{column}"
    )
    batch: list[tuple[str, str, str]] = []
    for row in iter_tsv(
        path,
        (
            "CB", "primary_mapped_reads", "candidate_countedU_reads",
            "nh_gt1_unique_gene_countedU_reads",
            "candidate_matrix_molecules", "biological_classification_status",
        ),
    ):
        batch.append(
            (
                canonical_barcode(row["CB"]),
                row["CB"],
                json.dumps(row, separators=(",", ":"), sort_keys=True),
            )
        )
        if len(batch) >= 5000:
            connection.executemany(sql, batch)
            batch.clear()
    if batch:
        connection.executemany(sql, batch)
    upsert_source_summary(connection, directory, prefix)
    return True


def iter_source_metrics(
    path: Path,
) -> Iterator[tuple[str, dict[str, int | str]]]:
    current_cb: str | None = None
    by_source: dict[str, int] = {}
    for row in iter_tsv(path, ("CB", "source_id", "candidate_countedU_reads")):
        cb = row["CB"]
        if current_cb is not None and cb != current_cb:
            total = sum(by_source.values())
            dominant, maximum = min(
                by_source.items(), key=lambda item: (-item[1], item[0])
            )
            yield current_cb, {
                "total": total,
                "maximum": maximum,
                "dominant": dominant,
                "sources": sum(value > 0 for value in by_source.values()),
            }
            by_source = {}
        current_cb = cb
        by_source[row["source_id"]] = (
            by_source.get(row["source_id"], 0)
            + int(row["candidate_countedU_reads"])
        )
    if current_cb is not None:
        total = sum(by_source.values())
        dominant, maximum = min(
            by_source.items(), key=lambda item: (-item[1], item[0])
        )
        yield current_cb, {
            "total": total,
            "maximum": maximum,
            "dominant": dominant,
            "sources": sum(value > 0 for value in by_source.values()),
        }


def upsert_source_summary(
    connection: sqlite3.Connection, directory: Path, prefix: str
) -> None:
    path = directory / "bam_evidence" / "barcode_rg_metrics.tsv.gz"
    column = f"{prefix}_source"
    sql = (
        f"INSERT INTO barcode(barcode_key,display_cb,{column}) VALUES(?,?,?) "
        f"ON CONFLICT(barcode_key) DO UPDATE SET {column}=excluded.{column}"
    )
    batch: list[tuple[str, str, str]] = []
    for cb, summary in iter_source_metrics(path):
        batch.append(
            (
                canonical_barcode(cb),
                cb,
                json.dumps(summary, separators=(",", ":"), sort_keys=True),
            )
        )
        if len(batch) >= 5000:
            connection.executemany(sql, batch)
            batch.clear()
    if batch:
        connection.executemany(sql, batch)


def ingest_optional(
    connection: sqlite3.Connection,
    path: Path | None,
    library: str,
    kind: str,
) -> None:
    if path is None:
        return
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or [])
        barcode_field = "CB" if "CB" in fields else "cell_barcode"
        if "library" not in fields or barcode_field not in fields:
            raise AnalysisError(
                f"optional {kind} table lacks library/{barcode_field}: {path}"
            )
        batch: list[tuple[str, str, str]] = []
        for row in reader:
            if row["library"] != library:
                continue
            key = canonical_barcode(row[barcode_field])
            payload = {
                field: value
                for field, value in row.items()
                if field not in {"library", "CB", "cell_barcode"}
            }
            batch.append(
                (kind, key, json.dumps(payload, separators=(",", ":"), sort_keys=True))
            )
            if len(batch) >= 5000:
                try:
                    connection.executemany(
                        "INSERT INTO optional_join(kind,barcode_key,payload) "
                        "VALUES(?,?,?)",
                        batch,
                    )
                except sqlite3.IntegrityError as exc:
                    raise AnalysisError(
                        f"duplicate optional {kind} barcode for {library}"
                    ) from exc
                batch.clear()
        if batch:
            try:
                connection.executemany(
                    "INSERT INTO optional_join(kind,barcode_key,payload) VALUES(?,?,?)",
                    batch,
                )
            except sqlite3.IntegrityError as exc:
                raise AnalysisError(
                    f"duplicate optional {kind} barcode for {library}"
                ) from exc


def optional_payloads(
    connection: sqlite3.Connection, barcode_key: str
) -> dict[str, dict[str, str]]:
    result: dict[str, dict[str, str]] = {}
    for kind, payload in connection.execute(
        "SELECT kind,payload FROM optional_join WHERE barcode_key=?",
        (barcode_key,),
    ):
        result[kind] = json.loads(payload)
    return result


def transition(current: int, old: int) -> str:
    if current and old:
        return "shared"
    if current:
        return "gained"
    if old:
        return "lost"
    return "background"


def optional_int(row: dict[str, str], field: str) -> int | None:
    value = row.get(field, "")
    return None if value == "" else int(value)


def build_row(
    connection: sqlite3.Connection,
    library: str,
    record: tuple[object, ...],
) -> dict[str, object]:
    (
        key, display_cb, current_filtered, old_filtered,
        current_umi, current_genes, current_rank,
        old_umi, old_genes, old_rank,
        current_bam_json, old_bam_json,
        current_source_json, old_source_json,
    ) = record
    current_bam = json.loads(current_bam_json) if current_bam_json else {}
    old_bam = json.loads(old_bam_json) if old_bam_json else {}
    source = (
        json.loads(current_source_json)
        if current_source_json
        else {"total": 0, "maximum": 0, "dominant": "", "sources": 0}
    )
    old_source = (
        json.loads(old_source_json)
        if old_source_json
        else {"total": 0, "maximum": 0, "dominant": "", "sources": 0}
    )
    category = transition(int(current_filtered), int(old_filtered))
    counted = int(current_bam.get("candidate_countedU_reads", 0))
    molecules = int(current_bam.get("candidate_matrix_molecules", 0))
    primary = int(current_bam.get("primary_mapped_reads", 0))
    all_records = int(current_bam.get("all_records", 0))
    unique_gene = int(current_bam.get("unique_gene_tagged_reads", 0))
    old_counted = int(old_bam.get("candidate_countedU_reads", 0))
    old_molecules = int(old_bam.get("candidate_matrix_molecules", 0))
    old_primary = int(old_bam.get("primary_mapped_reads", 0))
    authoritative_molecules = int(current_umi) if current_umi is not None else molecules
    saturation_value, saturation_status = conditional_saturation(
        counted,
        authoritative_molecules,
        (
            "ordinary_bam_reads_with_authoritative_nonzero_raw_matrix_umis"
            if current_umi is not None
            else "ordinary_bam_reads_with_candidate_bam_molecules"
        ),
    )
    if old_bam:
        old_saturation, old_saturation_status = conditional_saturation(
            old_counted,
            int(old_umi) if old_umi is not None else old_molecules,
            (
                "ordinary_bam_reads_with_authoritative_historical_matrix_umis"
                if old_umi is not None
                else "ordinary_bam_reads_with_candidate_historical_molecules"
            ),
        )
    else:
        old_saturation = ""
        old_saturation_status = "unavailable_no_historical_bam_evidence"
    mito = optional_int(current_bam, "mitochondrial_reads")
    rrna = optional_int(current_bam, "rrna_reads")
    old_mito = optional_int(old_bam, "mitochondrial_reads") if old_bam else None
    old_rrna = optional_int(old_bam, "rrna_reads") if old_bam else None
    evidence_sources = [
        name
        for name, present in (
            ("current_nonzero_raw_matrix", current_umi is not None),
            ("historical_nonzero_raw_matrix", old_umi is not None),
            ("current_BAM_observed_corrected_barcode", bool(current_bam)),
            ("historical_BAM_observed_corrected_barcode", bool(old_bam)),
        )
        if present
    ]
    row: dict[str, object] = {
        "library": library,
        "CB": display_cb,
        "cell_transition": category,
        "detail_row_scope": (
            "all_filtered_transition_cells"
            if category != "background"
            else "deterministic_evidence_bearing_background_sample"
        ),
        "background_inclusion_rule": BACKGROUND_RULE,
        "evidence_bearing_sources": ",".join(evidence_sources),
        "historical_filtered_member": int(old_filtered),
        "current_filtered_member": int(current_filtered),
        "current_nonzero_raw_matrix": int(current_umi is not None),
        "historical_nonzero_raw_matrix": int(old_umi is not None),
        "historical_raw_umi_rank": old_rank if old_rank is not None else "",
        "historical_raw_umis": old_umi if old_umi is not None else 0,
        "historical_raw_genes": old_genes if old_genes is not None else 0,
        "current_raw_umi_rank": current_rank if current_rank is not None else "",
        "current_raw_umis": current_umi if current_umi is not None else 0,
        "current_raw_genes": current_genes if current_genes is not None else 0,
        "candidate_countedU_reads": counted,
        "nh_gt1_unique_gene_countedU_reads":
            int(current_bam.get("nh_gt1_unique_gene_countedU_reads", 0)),
        "unique_gene_tagged_reads": unique_gene,
        "unique_gene_tagged_primary_fraction":
            safe_fraction(unique_gene, primary),
        "ordinary_candidate_primary_fraction":
            safe_fraction(counted, primary),
        "validated_counted_reads": counted if current_filtered else "",
        "candidate_matrix_molecules": molecules,
        "primary_mapped_reads": primary,
        "primary_mapping_record_fraction":
            safe_fraction(primary, all_records),
        "bam_conditional_saturation": saturation_value,
        "bam_conditional_saturation_status": saturation_status,
        "historical_bam_conditional_saturation": old_saturation,
        "historical_bam_conditional_saturation_status": old_saturation_status,
        "biological_classification_status":
            current_bam.get("biological_classification_status", "unavailable"),
        "mitochondrial_primary_read_fraction": safe_fraction(mito, primary),
        "rrna_primary_read_fraction": safe_fraction(rrna, primary),
        "historical_mitochondrial_primary_read_fraction":
            safe_fraction(old_mito, old_primary),
        "historical_rrna_primary_read_fraction":
            safe_fraction(old_rrna, old_primary),
        "contributing_sources": source["sources"],
        "dominant_source": source["dominant"],
        "source_dominance_basis":
            "ordinary_candidate_countedU_reads_by_physical_source",
        "source_dominance_fraction":
            safe_fraction(int(source["maximum"]), int(source["total"])),
        "historical_contributing_sources":
            old_source["sources"] if old_bam else "",
        "historical_dominant_source":
            old_source["dominant"] if old_bam else "",
        "historical_source_dominance_fraction": (
            safe_fraction(
                int(old_source["maximum"]), int(old_source["total"])
            )
            if old_bam else ""
        ),
    }
    for kind, payload in optional_payloads(connection, str(key)).items():
        for field, value in payload.items():
            row[f"{kind}_{field}"] = value
    return row


def evidence_records(
    connection: sqlite3.Connection,
) -> Iterator[tuple[object, ...]]:
    yield from connection.execute(
        """
        SELECT barcode_key,display_cb,current_filtered,old_filtered,
               current_umi,current_genes,current_rank,
               old_umi,old_genes,old_rank,
               current_bam,old_bam,current_source,old_source
        FROM barcode
        WHERE current_filtered=1 OR old_filtered=1
           OR current_umi IS NOT NULL OR old_umi IS NOT NULL
           OR current_bam IS NOT NULL OR old_bam IS NOT NULL
        ORDER BY barcode_key
        """
    )


def analyze_library(
    library: str,
    current: Path,
    baseline: Path,
    optional_paths: dict[str, Path | None],
    output_fields: list[str],
    global_writer: AtomicWriter,
    output: Path,
    background_sample_size: int,
    quantile_sample_size: int,
) -> tuple[dict[str, object], list[dict[str, object]], list[dict[str, object]]]:
    with tempfile.TemporaryDirectory(
        prefix=f".gained_{library}_", dir=output
    ) as temporary:
        connection = create_database(Path(temporary) / "barcode.sqlite")
        try:
            upsert_filtered(
                connection,
                find_product(current, "filtered/barcodes.tsv.gz"),
                "current_filtered",
            )
            upsert_filtered(
                connection,
                find_product(baseline, "filtered/barcodes.tsv.gz"),
                "old_filtered",
            )
            upsert_matrix(connection, current, "current")
            upsert_matrix(connection, baseline, "old")
            if not upsert_bam(connection, current, "current"):
                raise AnalysisError(f"current BAM evidence is missing for {library}")
            upsert_bam(connection, baseline, "old")
            for kind, path in optional_paths.items():
                ingest_optional(connection, path, library, kind)
            connection.commit()

            counts = {"shared": 0, "gained": 0, "lost": 0, "background": 0}
            distributions: dict[tuple[str, str], Distribution] = {}
            background_heap: list[tuple[int, str]] = []
            for record in evidence_records(connection):
                row = build_row(connection, library, record)
                category = str(row["cell_transition"])
                counts[category] += 1
                cb = str(row["CB"])
                if category == "background" and background_sample_size:
                    score = deterministic_hash(library, cb, "background_detail")
                    item = (-score, canonical_barcode(cb))
                    if len(background_heap) < background_sample_size:
                        heapq.heappush(background_heap, item)
                    elif score < -background_heap[0][0]:
                        heapq.heapreplace(background_heap, item)
                for metric in (*CORE_METRICS, *TRIM_METRICS):
                    value = row.get(metric, "")
                    if value == "":
                        continue
                    key = (category, metric)
                    target = distributions.setdefault(
                        key,
                        Distribution(
                            library, category, metric, quantile_sample_size
                        ),
                    )
                    target.add(cb, value)

            selected_background = {item[1] for item in background_heap}
            library_path = output / "libraries" / library / "cell_category_metrics.tsv.gz"
            with AtomicWriter(
                library_path, output_fields, compressed=True
            ) as library_writer:
                for record in evidence_records(connection):
                    category = transition(int(record[2]), int(record[3]))
                    if (
                        category == "background"
                        and str(record[0]) not in selected_background
                    ):
                        continue
                    row = build_row(connection, library, record)
                    library_writer.writerow(row)
                    global_writer.writerow(row)

            histogram_rows = [
                row
                for target in distributions.values()
                for row in target.histogram_rows()
            ]
            quantile_rows = [
                row
                for target in distributions.values()
                for row in target.quantile_rows()
            ]
            transition_row = {
                "library": library,
                "comparison_status": "complete",
                **counts,
                "background_detail_rows": len(selected_background),
                "background_detail_cap": background_sample_size,
                "background_inclusion_rule": BACKGROUND_RULE,
                "background_distribution_scope":
                    "all_evidence_bearing_background_barcodes",
            }
            return transition_row, histogram_rows, quantile_rows
        finally:
            connection.close()


def run(args: argparse.Namespace) -> None:
    current_dirs = library_dirs(Path(args.current_root).resolve())
    baseline_dirs = library_dirs(Path(args.baseline_root).resolve())
    requested = set(getattr(args, "library", None) or [])
    selected = sorted(requested or (set(current_dirs) | set(baseline_dirs)))
    unknown = requested - (set(current_dirs) | set(baseline_dirs))
    if unknown:
        raise AnalysisError(
            f"requested libraries are unavailable: {', '.join(sorted(unknown))}"
        )
    optional_paths = {
        "demux": Path(args.demux_tsv).resolve() if args.demux_tsv else None,
        "trim": Path(args.trim_tsv).resolve() if args.trim_tsv else None,
        "ambient": (
            Path(args.ambient_similarity_tsv).resolve()
            if args.ambient_similarity_tsv else None
        ),
    }
    core_fields = [
        "library", "CB", "cell_transition", "detail_row_scope",
        "background_inclusion_rule", "evidence_bearing_sources",
        "historical_filtered_member", "current_filtered_member",
        "current_nonzero_raw_matrix", "historical_nonzero_raw_matrix",
        "historical_raw_umi_rank", "historical_raw_umis",
        "historical_raw_genes", "current_raw_umi_rank", "current_raw_umis",
        "current_raw_genes", "candidate_countedU_reads",
        "nh_gt1_unique_gene_countedU_reads", "validated_counted_reads",
        "unique_gene_tagged_reads", "unique_gene_tagged_primary_fraction",
        "ordinary_candidate_primary_fraction",
        "candidate_matrix_molecules", "primary_mapped_reads",
        "primary_mapping_record_fraction",
        "bam_conditional_saturation", "bam_conditional_saturation_status",
        "historical_bam_conditional_saturation",
        "historical_bam_conditional_saturation_status",
        "biological_classification_status",
        "mitochondrial_primary_read_fraction",
        "rrna_primary_read_fraction",
        "historical_mitochondrial_primary_read_fraction",
        "historical_rrna_primary_read_fraction",
        "contributing_sources", "dominant_source",
        "source_dominance_basis", "source_dominance_fraction",
        "historical_contributing_sources", "historical_dominant_source",
        "historical_source_dominance_fraction",
    ]
    dynamic_fields = sorted(
        {
            field
            for kind, path in optional_paths.items()
            for field in optional_fields(path, kind + "_")
        }
    )
    output_fields = core_fields + dynamic_fields
    output = Path(args.output_dir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    transitions: list[dict[str, object]] = []
    histogram_rows: list[dict[str, object]] = []
    quantile_rows: list[dict[str, object]] = []
    global_path = output / "cell_category_metrics.tsv.gz"
    with AtomicWriter(global_path, output_fields, compressed=True) as global_writer:
        for library in selected:
            if library not in current_dirs or library not in baseline_dirs:
                transitions.append(
                    {
                        "library": library,
                        "comparison_status":
                            "missing_current" if library not in current_dirs
                            else "missing_historical",
                        "shared": "", "gained": "", "lost": "",
                        "background": "", "background_detail_rows": "",
                        "background_detail_cap": getattr(
                            args, "background_sample_size", 10000
                        ),
                        "background_inclusion_rule": BACKGROUND_RULE,
                        "background_distribution_scope": "",
                    }
                )
                continue
            transition_row, hist, quantiles = analyze_library(
                library,
                current_dirs[library],
                baseline_dirs[library],
                optional_paths,
                output_fields,
                global_writer,
                output,
                getattr(args, "background_sample_size", 10000),
                getattr(args, "quantile_sample_size", 4096),
            )
            transitions.append(transition_row)
            # These summaries are bounded: fixed bins and fixed samples per
            # category/metric/library, never one item per raw barcode.
            histogram_rows.extend(hist)
            quantile_rows.extend(quantiles)

    write_tsv(
        output / "library_cell_transitions.tsv",
        [
            "library", "comparison_status", "shared", "gained", "lost",
            "background", "background_detail_rows", "background_detail_cap",
            "background_inclusion_rule", "background_distribution_scope",
        ],
        transitions,
    )
    write_tsv(
        output / "cell_category_histograms.tsv",
        [
            "library", "cell_transition", "metric", "bin_index",
            "bin_left_inclusive", "bin_right_exclusive", "count",
            "population_n", "histogram_scope",
        ],
        sorted(
            histogram_rows,
            key=lambda row: (
                str(row["library"]), str(row["cell_transition"]),
                str(row["metric"]), int(row["bin_index"]),
            ),
        ),
    )
    write_tsv(
        output / "cell_category_quantiles.tsv",
        [
            "library", "cell_transition", "metric", "quantile", "value",
            "population_n", "sample_n", "method",
        ],
        sorted(
            quantile_rows,
            key=lambda row: (
                str(row["library"]), str(row["cell_transition"]),
                str(row["metric"]), float(row["quantile"]),
            ),
        ),
    )
    write_tsv(
        output / "background_contract.tsv",
        [
            "release", "inclusion_rule", "detail_rule",
            "distribution_rule", "raw_whitelist_zero_rows_exported",
        ],
        [{
            "release": RELEASE,
            "inclusion_rule": BACKGROUND_RULE,
            "detail_rule":
                "deterministic_smallest_sha256_sample_per_library_cap_"
                f"{getattr(args, 'background_sample_size', 10000)}",
            "distribution_rule":
                "exact_fixed_histograms_over_all_evidence_bearing_rows",
            "raw_whitelist_zero_rows_exported": 0,
        }],
    )
    print(
        f"Processed {len(selected)} independent library shard(s); "
        f"background detail cap {getattr(args, 'background_sample_size', 10000)} per library"
    )


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--current-root", required=True)
    result.add_argument("--baseline-root", required=True)
    result.add_argument("--output-dir", required=True)
    result.add_argument(
        "--library",
        action="append",
        help="process only this library; repeat for independent shards",
    )
    result.add_argument("--demux-tsv")
    result.add_argument("--trim-tsv")
    result.add_argument("--ambient-similarity-tsv")
    result.add_argument("--background-sample-size", type=int, default=10000)
    result.add_argument("--quantile-sample-size", type=int, default=4096)
    return result


if __name__ == "__main__":
    try:
        arguments = parser().parse_args()
        if arguments.background_sample_size < 0 or arguments.quantile_sample_size < 1:
            raise AnalysisError("sample sizes must be nonnegative/positive")
        run(arguments)
    except AnalysisError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)
