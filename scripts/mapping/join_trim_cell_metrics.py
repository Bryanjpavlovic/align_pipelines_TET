#!/usr/bin/env python3
"""Left-join source/mate trimming rows to BAM evidence with disk-backed state.

``trim_cell_metrics.tsv.gz`` preserves every trimming input row. The companion
``trim_cell_metrics_by_barcode.tsv.gz`` contains exactly one library/barcode
row for ``analyze_gained_cells.py --trim-tsv``. Counts are summed, fractions
are recomputed from summed numerators and denominators, and mean-only measures
are weighted by ``total_reads``. Already-calculated fractions are never summed.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sqlite3
import sys
import tempfile
from pathlib import Path
from typing import Iterator

from rna_evidence_common import (
    AnalysisError,
    find_product,
    iter_tsv,
    library_dirs,
    open_text,
    write_tsv,
)


RELEASE = "2026-09-05-v2-streaming-left-join"
REQUIRED_TRIM_FIELDS = {
    "library", "cell_barcode", "source_id", "read_mate", "total_reads",
}
COUNT_FIELDS = (
    "total_reads", "barcode_starsolo_cb_reads", "barcode_exact_reads",
    "barcode_unique_1mm_reads", "reads_with_adapter", "reads_no_adapter",
    "reads_too_short", "total_bp_trimmed", "multi_adapter_reads",
    "fixed_tso_trimmed_reads", "tso_unrecognized_reads",
)
TRIM_X_METRICS = (
    "frac_reads_with_adapter", "mean_bp_trimmed", "fixed_tso_trimmed_reads",
    "tso_unrecognized_reads", "reads_too_short",
)
STORED_EXTRA_METRICS = (
    "frac_reads_with_adapter", "mean_bp_trimmed",
    "mean_original_length", "mean_final_length",
)
BAM_Y_METRICS = (
    "bam_primary_mapped_per_trim_read_proxy",
    "bam_unique_gene_assigned_per_trim_read_proxy",
    "bam_mitochondrial_primary_fraction",
    "bam_rrna_primary_fraction",
    "bam_matrix_conditional_saturation", "matrix_umis", "matrix_genes",
)
BAM_INTEGER_FIELDS = (
    "primary_mapped_reads", "unique_gene_tagged_reads",
    "candidate_countedU_reads", "candidate_matrix_molecules",
)
DETAIL_ADDITIONS = [
    "bam_match_status", "RG", "barcode_correction_conflict_observed",
    "bam_primary_mapped_reads", "bam_unique_gene_tagged_reads",
    "bam_candidate_countedU_reads", "bam_candidate_matrix_molecules",
    "matrix_umis", "matrix_genes", "matrix_value_authority",
    "bam_primary_mapped_per_trim_read_proxy",
    "bam_unique_gene_assigned_per_trim_read_proxy",
    "bam_mitochondrial_primary_fraction", "bam_rrna_primary_fraction",
    "bam_biological_classification_status",
    "bam_matrix_conditional_saturation",
    "bam_matrix_conditional_saturation_status", "causal_interpretation",
]


def canonical_barcode(value: str) -> str:
    return value[:-2] if value.endswith("-1") else value


def fraction(numerator: int | float | None, denominator: int) -> str:
    if numerator is None or denominator <= 0:
        return ""
    return f"{float(numerator) / denominator:.8f}"


def optional_int(value: str | None, label: str) -> int | None:
    if value in {None, ""}:
        return None
    try:
        result = int(value)
    except ValueError as exc:
        raise AnalysisError(f"non-integer {label}: {value!r}") from exc
    if result < 0:
        raise AnalysisError(f"negative {label}: {value!r}")
    return result


def optional_float(value: str | None, label: str) -> float | None:
    if value in {None, ""}:
        return None
    try:
        result = float(value)
    except ValueError as exc:
        raise AnalysisError(f"non-numeric {label}: {value!r}") from exc
    if not math.isfinite(result):
        raise AnalysisError(f"non-finite {label}: {value!r}")
    return result


def saturation(reads: int, molecules: int | None) -> tuple[str, str]:
    if molecules is None:
        return "", "unavailable_no_matrix_molecule_measurement"
    if reads <= 0:
        return "", "unavailable_no_candidate_counted_reads"
    if molecules > reads:
        return "", "unavailable_molecule_count_exceeds_read_count"
    return (
        f"{(reads - molecules) / reads:.8f}",
        "bam_reads_with_authoritative_matrix_umis",
    )


def create_database(path: Path) -> sqlite3.Connection:
    connection = sqlite3.connect(path)
    connection.execute("PRAGMA journal_mode=OFF")
    connection.execute("PRAGMA synchronous=OFF")
    connection.execute("PRAGMA temp_store=FILE")
    numeric_columns = ",\n".join(
        f"{field} REAL" for field in (*COUNT_FIELDS, *STORED_EXTRA_METRICS)
    )
    connection.executescript(
        f"""
        CREATE TABLE trim_row (
            row_id INTEGER PRIMARY KEY,
            library TEXT NOT NULL,
            barcode TEXT NOT NULL,
            source_id TEXT NOT NULL,
            read_mate TEXT NOT NULL,
            barcode_population TEXT NOT NULL,
            payload TEXT NOT NULL,
            {numeric_columns}
        );
        CREATE INDEX trim_library_key
            ON trim_row(library,barcode,source_id,row_id);
        CREATE TABLE bam_rg (
            library TEXT NOT NULL,
            barcode TEXT NOT NULL,
            rg TEXT NOT NULL,
            primary_mapped_reads INTEGER NOT NULL,
            unique_gene_tagged_reads INTEGER NOT NULL,
            candidate_countedU_reads INTEGER NOT NULL,
            candidate_matrix_molecules INTEGER NOT NULL,
            mitochondrial_reads INTEGER,
            rrna_reads INTEGER,
            biological_classification_status TEXT NOT NULL,
            PRIMARY KEY(library,barcode,rg)
        ) WITHOUT ROWID;
        CREATE TABLE correction_conflict (
            library TEXT NOT NULL,
            barcode TEXT NOT NULL,
            rg TEXT NOT NULL,
            PRIMARY KEY(library,barcode,rg)
        ) WITHOUT ROWID;
        CREATE TABLE bam_barcode_total (
            library TEXT NOT NULL,
            barcode TEXT NOT NULL,
            primary_mapped_reads INTEGER NOT NULL,
            unique_gene_tagged_reads INTEGER NOT NULL,
            candidate_countedU_reads INTEGER NOT NULL,
            candidate_matrix_molecules INTEGER NOT NULL,
            mitochondrial_reads INTEGER,
            rrna_reads INTEGER,
            biological_classification_status TEXT NOT NULL,
            PRIMARY KEY(library,barcode)
        ) WITHOUT ROWID;
        CREATE TABLE matrix_value (
            library TEXT NOT NULL,
            barcode TEXT NOT NULL,
            umis INTEGER NOT NULL,
            genes INTEGER NOT NULL,
            PRIMARY KEY(library,barcode)
        ) WITHOUT ROWID;
        CREATE TABLE joined_metric (
            row_id INTEGER PRIMARY KEY,
            library TEXT NOT NULL,
            barcode TEXT NOT NULL,
            barcode_population TEXT NOT NULL,
            bam_match_status TEXT NOT NULL,
            frac_reads_with_adapter REAL,
            mean_bp_trimmed REAL,
            fixed_tso_trimmed_reads REAL,
            tso_unrecognized_reads REAL,
            reads_too_short REAL,
            bam_primary_mapped_per_trim_read_proxy REAL,
            bam_unique_gene_assigned_per_trim_read_proxy REAL,
            bam_mitochondrial_primary_fraction REAL,
            bam_rrna_primary_fraction REAL,
            bam_matrix_conditional_saturation REAL,
            matrix_umis REAL,
            matrix_genes REAL
        );
        CREATE INDEX joined_group
            ON joined_metric(library,barcode_population,bam_match_status);
        """
    )
    return connection


def read_trim_into_database(
    connection: sqlite3.Connection, path: Path
) -> tuple[list[str], int]:
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        missing = REQUIRED_TRIM_FIELDS - set(fields)
        if missing:
            raise AnalysisError(f"{path} lacks columns: {', '.join(sorted(missing))}")
        if len(fields) != len(set(fields)):
            raise AnalysisError(f"duplicate trim-table columns: {path}")
        columns = [
            "library", "barcode", "source_id", "read_mate",
            "barcode_population", "payload", *COUNT_FIELDS,
            *STORED_EXTRA_METRICS,
        ]
        placeholders = ",".join("?" for _ in columns)
        sql = f"INSERT INTO trim_row({','.join(columns)}) VALUES({placeholders})"
        batch: list[tuple[object, ...]] = []
        total = 0
        for line_number, row in enumerate(reader, start=2):
            library = row["library"]
            barcode = canonical_barcode(row["cell_barcode"])
            source_id = row["source_id"]
            read_mate = row["read_mate"]
            if not library or not barcode or not source_id or not read_mate:
                raise AnalysisError(f"blank trim join key in {path}:{line_number}")
            counts = {
                field: optional_int(row.get(field), f"{field} at {path}:{line_number}")
                for field in COUNT_FIELDS
            }
            if counts["total_reads"] is None:
                raise AnalysisError(f"blank total_reads in {path}:{line_number}")
            metrics = {
                field: optional_float(row.get(field), f"{field} at {path}:{line_number}")
                for field in STORED_EXTRA_METRICS
            }
            batch.append(
                (
                    library, barcode, source_id, read_mate,
                    row.get("barcode_population", "all") or "all",
                    json.dumps(row, sort_keys=True, separators=(",", ":")),
                    *(counts[field] for field in COUNT_FIELDS),
                    *(metrics[field] for field in STORED_EXTRA_METRICS),
                )
            )
            total += 1
            if len(batch) >= 5000:
                connection.executemany(sql, batch)
                batch.clear()
        if batch:
            connection.executemany(sql, batch)
    if total == 0:
        raise AnalysisError(f"trim table has no rows: {path}")
    connection.commit()
    return fields, total


def matrix_subset(directory: Path, targets: set[str]) -> Iterator[tuple[str, int, int]]:
    """Stream one raw matrix and retain only the current library's join targets."""
    roster_path = find_product(directory, "raw/barcodes.tsv.gz")
    matrix_path = find_product(directory, "raw/matrix.mtx.gz")
    columns: dict[int, str] = {}
    seen_targets: set[str] = set()
    roster_count = 0
    with open_text(roster_path) as handle:
        for roster_count, line in enumerate(handle, start=1):
            barcode = canonical_barcode(line.rstrip("\r\n").split("\t", 1)[0])
            if barcode in targets:
                if barcode in seen_targets:
                    raise AnalysisError(f"duplicate canonical barcode in {roster_path}: {barcode}")
                columns[roster_count] = barcode
                seen_targets.add(barcode)
    values = {barcode: [0, 0] for barcode in seen_targets}
    with open_text(matrix_path) as handle:
        banner = handle.readline().lower()
        if not banner.startswith("%%matrixmarket matrix coordinate"):
            raise AnalysisError(f"unsupported MatrixMarket file: {matrix_path}")
        line = handle.readline()
        while line.startswith("%"):
            line = handle.readline()
        try:
            feature_count, barcode_count, declared = map(int, line.split())
        except ValueError as exc:
            raise AnalysisError(f"malformed MatrixMarket dimensions: {matrix_path}") from exc
        if barcode_count != roster_count:
            raise AnalysisError(f"matrix/barcode dimension mismatch: {matrix_path}")
        observed = 0
        for line_number, line in enumerate(handle, start=3):
            if not line.strip() or line.startswith("%"):
                continue
            try:
                feature, column, count = map(int, line.split())
            except ValueError as exc:
                raise AnalysisError(f"malformed matrix row {line_number}: {matrix_path}") from exc
            if not 1 <= feature <= feature_count or not 1 <= column <= barcode_count or count <= 0:
                raise AnalysisError(f"invalid matrix row {line_number}: {matrix_path}")
            observed += 1
            barcode = columns.get(column)
            if barcode is not None:
                values[barcode][0] += count
                values[barcode][1] += 1
        if observed != declared:
            raise AnalysisError(f"MatrixMarket nnz mismatch: {matrix_path}")
    for barcode in sorted(values):
        yield barcode, values[barcode][0], values[barcode][1]


def ingest_library_evidence(
    connection: sqlite3.Connection, library: str, directory: Path
) -> None:
    targets = {
        row[0]
        for row in connection.execute(
            "SELECT DISTINCT barcode FROM trim_row WHERE library=?", (library,)
        )
    }
    evidence = directory / "bam_evidence"
    for row in iter_tsv(
        evidence / "barcode_rg_metrics.tsv.gz",
        (
            "CB", "RG", "primary_mapped_reads", "unique_gene_tagged_reads",
            "candidate_countedU_reads", "candidate_matrix_molecules",
        ),
    ):
        barcode = canonical_barcode(row["CB"])
        if barcode not in targets:
            continue
        values = [optional_int(row.get(field), field) for field in BAM_INTEGER_FIELDS]
        if any(value is None for value in values):
            raise AnalysisError(f"blank required BAM metric for {library}/{barcode}/{row['RG']}")
        try:
            connection.execute(
                "INSERT INTO bam_rg VALUES(?,?,?,?,?,?,?,?,?,?)",
                (
                    library, barcode, row["RG"], *values,
                    optional_int(row.get("mitochondrial_reads"), "mitochondrial_reads"),
                    optional_int(row.get("rrna_reads"), "rrna_reads"),
                    row.get("biological_classification_status", "unavailable"),
                ),
            )
        except sqlite3.IntegrityError as exc:
            raise AnalysisError(
                f"duplicate barcode/RG evidence key: {(library, barcode, row['RG'])!r}"
            ) from exc
    for row in iter_tsv(
        evidence / "barcode_read_metrics.tsv.gz",
        (
            "CB", "primary_mapped_reads", "unique_gene_tagged_reads",
            "candidate_countedU_reads", "candidate_matrix_molecules",
        ),
    ):
        barcode = canonical_barcode(row["CB"])
        if barcode not in targets:
            continue
        values = [optional_int(row.get(field), field) for field in BAM_INTEGER_FIELDS]
        if any(value is None for value in values):
            raise AnalysisError(f"blank required barcode BAM metric for {library}/{barcode}")
        try:
            connection.execute(
                "INSERT INTO bam_barcode_total VALUES(?,?,?,?,?,?,?,?,?)",
                (
                    library, barcode, *values,
                    optional_int(row.get("mitochondrial_reads"), "mitochondrial_reads"),
                    optional_int(row.get("rrna_reads"), "rrna_reads"),
                    row.get("biological_classification_status", "unavailable"),
                ),
            )
        except sqlite3.IntegrityError as exc:
            raise AnalysisError(
                f"duplicate barcode BAM evidence key: {(library, barcode)!r}"
            ) from exc
    for row in iter_tsv(
        evidence / "raw_to_corrected_barcode_counts.tsv.gz",
        ("RG", "CB", "within_rg_conflict", "cross_rg_conflict"),
    ):
        barcode = canonical_barcode(row["CB"])
        if barcode in targets and row["within_rg_conflict"] == "1":
            connection.execute(
                "INSERT OR IGNORE INTO correction_conflict VALUES(?,?,?)",
                (library, barcode, row["RG"]),
            )
    connection.executemany(
        "INSERT INTO matrix_value VALUES(?,?,?,?)",
        (
            (library, barcode, umis, genes)
            for barcode, umis, genes in matrix_subset(directory, targets)
        ),
    )
    connection.commit()


def detail_rows(
    connection: sqlite3.Connection,
    counters: dict[str, int],
) -> Iterator[dict[str, object]]:
    query = """
        SELECT t.row_id,t.payload,t.library,t.barcode,t.source_id,t.total_reads,
               t.barcode_population,
               b.rg,b.primary_mapped_reads,b.unique_gene_tagged_reads,
               b.candidate_countedU_reads,b.candidate_matrix_molecules,
               b.mitochondrial_reads,b.rrna_reads,
               b.biological_classification_status,
               CASE WHEN c.rg IS NULL THEN 0 ELSE 1 END,
               m.umis,m.genes,
               t.frac_reads_with_adapter,t.mean_bp_trimmed,
               t.fixed_tso_trimmed_reads,t.tso_unrecognized_reads,t.reads_too_short
        FROM trim_row AS t
        LEFT JOIN bam_rg AS b
          ON b.library=t.library AND b.barcode=t.barcode AND b.rg=t.source_id
        LEFT JOIN correction_conflict AS c
          ON c.library=t.library AND c.barcode=t.barcode AND c.rg=t.source_id
        LEFT JOIN matrix_value AS m
          ON m.library=t.library AND m.barcode=t.barcode
        ORDER BY t.library,t.row_id
    """
    insert = "INSERT INTO joined_metric VALUES(" + ",".join("?" for _ in range(17)) + ")"
    batch: list[tuple[object, ...]] = []
    for values in connection.execute(query):
        (
            row_id, payload, library, barcode, source_id, total_reads,
            population, rg, primary, unique_gene, candidate, molecules,
            mito, rrna, class_status, conflict, matrix_umis, matrix_genes,
            frac_adapter, mean_trimmed, fixed_tso, tso_unrecognized, too_short,
        ) = values
        row: dict[str, object] = json.loads(payload)
        matched = rg is not None
        status = "matched" if matched else "no_bam_evidence"
        counters[status] = counters.get(status, 0) + 1
        counters["conflict"] = counters.get("conflict", 0) + (conflict if matched else 0)
        if matched:
            sat, sat_status = saturation(
                int(candidate),
                int(matrix_umis) if matrix_umis is not None else int(molecules),
            )
            updates: dict[str, object] = {
                "bam_match_status": status, "RG": rg,
                "barcode_correction_conflict_observed": conflict,
                "bam_primary_mapped_reads": primary,
                "bam_unique_gene_tagged_reads": unique_gene,
                "bam_candidate_countedU_reads": candidate,
                "bam_candidate_matrix_molecules": molecules,
                "matrix_umis": "" if matrix_umis is None else matrix_umis,
                "matrix_genes": "" if matrix_genes is None else matrix_genes,
                "matrix_value_authority": (
                    "STARsolo_raw_matrix"
                    if matrix_umis is not None else "unavailable_no_matrix_value"
                ),
                "bam_primary_mapped_per_trim_read_proxy": fraction(primary, int(total_reads)),
                "bam_unique_gene_assigned_per_trim_read_proxy": fraction(unique_gene, int(total_reads)),
                "bam_mitochondrial_primary_fraction": fraction(mito, int(primary)),
                "bam_rrna_primary_fraction": fraction(rrna, int(primary)),
                "bam_biological_classification_status": class_status,
                "bam_matrix_conditional_saturation": sat,
                "bam_matrix_conditional_saturation_status": sat_status,
                "causal_interpretation": "correlation_only",
            }
        else:
            updates = {field: "" for field in DETAIL_ADDITIONS}
            updates.update({
                "bam_match_status": status,
                "bam_biological_classification_status": "unavailable_no_bam_evidence",
                "bam_matrix_conditional_saturation_status": "unavailable_no_bam_evidence",
                "matrix_value_authority": "unavailable_no_bam_evidence",
                "causal_interpretation": "trimming_only_no_bam_evidence",
            })
        row.update(updates)
        y_values = [
            optional_float(str(row.get(field, "")), field)
            for field in BAM_Y_METRICS
        ]
        batch.append((
            row_id, library, barcode, population, status,
            frac_adapter, mean_trimmed, fixed_tso, tso_unrecognized, too_short,
            *y_values,
        ))
        if len(batch) >= 5000:
            connection.executemany(insert, batch)
            batch.clear()
        yield row
    if batch:
        connection.executemany(insert, batch)
    connection.commit()


def ordered_distinct_lists(connection: sqlite3.Connection) -> None:
    for table, column in (
        ("aggregate_sources", "source_id"),
        ("aggregate_mates", "read_mate"),
    ):
        connection.execute(
            f"""
            CREATE TABLE {table} AS
            SELECT library,barcode,GROUP_CONCAT({column}, ',') AS values_text,
                   COUNT(*) AS value_count
            FROM (
                SELECT DISTINCT library,barcode,{column}
                FROM trim_row ORDER BY library,barcode,{column}
            ) GROUP BY library,barcode
            """
        )
        connection.execute(
            f"CREATE UNIQUE INDEX {table}_key ON {table}(library,barcode)"
        )


def aggregate_rows(connection: sqlite3.Connection) -> Iterator[dict[str, object]]:
    ordered_distinct_lists(connection)
    connection.execute(
        """
        CREATE TEMP TABLE represented_bam AS
        SELECT b.* FROM bam_rg AS b
        JOIN (
            SELECT DISTINCT library,barcode,source_id FROM trim_row
        ) AS t
          ON t.library=b.library AND t.barcode=b.barcode AND t.source_id=b.rg
        """
    )
    connection.execute(
        """
        CREATE TEMP TABLE aggregate_matched_rgs AS
        SELECT library,barcode,GROUP_CONCAT(rg, ',') AS values_text
        FROM (
            SELECT DISTINCT library,barcode,rg FROM represented_bam
            ORDER BY library,barcode,rg
        ) GROUP BY library,barcode
        """
    )
    sums = ",".join(
        (
            "SUM(t.total_reads) AS total_reads"
            if field == "total_reads" else
            f"CASE WHEN COUNT(t.{field})=COUNT(*) THEN SUM(t.{field}) END AS {field}"
        )
        for field in COUNT_FIELDS
    )
    connection.execute(
        f"""
        CREATE TEMP TABLE trim_barcode AS
        SELECT t.library,t.barcode,COUNT(*) AS trim_row_count,
               SUM(CASE WHEN j.bam_match_status='matched' THEN 1 ELSE 0 END) AS matched_rows,
               SUM(CASE WHEN j.bam_match_status='no_bam_evidence' THEN 1 ELSE 0 END) AS unmatched_rows,
               CASE WHEN MIN(t.barcode_population)=MAX(t.barcode_population)
                    THEN MIN(t.barcode_population) ELSE 'mixed' END AS barcode_population,
               {sums},
               CASE WHEN SUM(t.total_reads)>0 AND COUNT(t.frac_reads_with_adapter)>0
                    THEN SUM(t.frac_reads_with_adapter*t.total_reads)/
                         SUM(CASE WHEN t.frac_reads_with_adapter IS NOT NULL THEN t.total_reads END) END AS weighted_adapter_fraction,
               CASE WHEN SUM(t.total_reads)>0 AND COUNT(t.mean_bp_trimmed)>0
                    THEN SUM(t.mean_bp_trimmed*t.total_reads)/
                         SUM(CASE WHEN t.mean_bp_trimmed IS NOT NULL THEN t.total_reads END) END AS weighted_mean_bp_trimmed,
               CASE WHEN SUM(t.total_reads)>0 AND COUNT(t.mean_original_length)>0
                    THEN SUM(t.mean_original_length*t.total_reads)/
                         SUM(CASE WHEN t.mean_original_length IS NOT NULL THEN t.total_reads END) END AS mean_original_length,
               CASE WHEN SUM(t.total_reads)>0 AND COUNT(t.mean_final_length)>0
                    THEN SUM(t.mean_final_length*t.total_reads)/
                         SUM(CASE WHEN t.mean_final_length IS NOT NULL THEN t.total_reads END) END AS mean_final_length
        FROM trim_row AS t JOIN joined_metric AS j ON j.row_id=t.row_id
        GROUP BY t.library,t.barcode
        """
    )
    query = """
        SELECT t.*,s.values_text,s.value_count,m.values_text,m.value_count,
               b.primary_mapped_reads,b.unique_gene_tagged_reads,
               b.candidate_countedU_reads,b.candidate_matrix_molecules,
               b.mitochondrial_reads,b.rrna_reads,b.biological_classification_status,
               v.umis,v.genes,r.values_text,
               EXISTS(
                   SELECT 1 FROM correction_conflict c
                   JOIN represented_bam r
                     ON r.library=c.library AND r.barcode=c.barcode AND r.rg=c.rg
                   WHERE c.library=t.library AND c.barcode=t.barcode
               )
        FROM trim_barcode t
        JOIN aggregate_sources s USING(library,barcode)
        JOIN aggregate_mates m USING(library,barcode)
        LEFT JOIN bam_barcode_total b USING(library,barcode)
        LEFT JOIN matrix_value v USING(library,barcode)
        LEFT JOIN aggregate_matched_rgs r USING(library,barcode)
        ORDER BY t.library,t.barcode
    """
    for values in connection.execute(query):
        (
            library, barcode, trim_count, matched, unmatched, population,
            total_reads, starsolo_cb, exact, one_mm, reads_adapter, reads_no_adapter,
            reads_short, total_bp_trimmed, multi_adapter, fixed_tso, tso_unrecognized,
            weighted_adapter, weighted_trimmed, mean_original, mean_final,
            sources, source_count, mates, mate_count,
            primary, unique_gene, candidate, molecules, mito, rrna, class_status,
            matrix_umis, matrix_genes, matched_rgs, conflict,
        ) = values
        if matched == 0:
            match_status = "no_bam_evidence"
        elif unmatched == 0:
            match_status = "matched"
        else:
            match_status = "partial"
        if matched and primary is None:
            raise AnalysisError(
                f"barcode-level BAM evidence is missing for {library}/{barcode}"
            )
        authoritative_molecules = matrix_umis if matrix_umis is not None else molecules
        sat, sat_status = (
            saturation(int(candidate or 0), authoritative_molecules)
            if matched else ("", "unavailable_no_bam_evidence")
        )
        yield {
            "library": library, "cell_barcode": barcode,
            "barcode_population": population, "bam_match_status": match_status,
            "trim_row_count": trim_count, "bam_matched_trim_rows": matched,
            "bam_unmatched_trim_rows": unmatched, "source_ids": sources,
            "source_count": source_count, "read_mates": mates,
            "read_mate_count": mate_count, "total_reads": int(total_reads),
            "barcode_starsolo_cb_reads": "" if starsolo_cb is None else int(starsolo_cb),
            "barcode_exact_reads": "" if exact is None else int(exact),
            "barcode_unique_1mm_reads": "" if one_mm is None else int(one_mm),
            "reads_with_adapter": "" if reads_adapter is None else int(reads_adapter),
            "reads_no_adapter": "" if reads_no_adapter is None else int(reads_no_adapter),
            "reads_too_short": "" if reads_short is None else int(reads_short),
            "frac_reads_with_adapter": (
                fraction(reads_adapter, int(total_reads))
                if reads_adapter is not None
                else ("" if weighted_adapter is None else f"{weighted_adapter:.8f}")
            ),
            "total_bp_trimmed": "" if total_bp_trimmed is None else int(total_bp_trimmed),
            "mean_bp_trimmed": (
                fraction(total_bp_trimmed, int(total_reads))
                if total_bp_trimmed is not None
                else ("" if weighted_trimmed is None else f"{weighted_trimmed:.8f}")
            ),
            "multi_adapter_reads": "" if multi_adapter is None else int(multi_adapter),
            "mean_original_length": "" if mean_original is None else f"{mean_original:.8f}",
            "mean_final_length": "" if mean_final is None else f"{mean_final:.8f}",
            "fixed_tso_trimmed_reads": "" if fixed_tso is None else int(fixed_tso),
            "tso_unrecognized_reads": "" if tso_unrecognized is None else int(tso_unrecognized),
            "RG": matched_rgs if matched else "",
            "barcode_correction_conflict_observed": "" if not matched else int(conflict),
            "bam_primary_mapped_reads": "" if not matched else int(primary),
            "bam_unique_gene_tagged_reads": "" if not matched else int(unique_gene),
            "bam_candidate_countedU_reads": "" if not matched else int(candidate),
            "bam_candidate_matrix_molecules": "" if not matched else int(molecules),
            "matrix_umis": "" if not matched or matrix_umis is None else int(matrix_umis),
            "matrix_genes": "" if not matched or matrix_genes is None else int(matrix_genes),
            "matrix_value_authority": (
                "STARsolo_raw_matrix" if matched and matrix_umis is not None
                else (
                    "unavailable_no_matrix_value"
                    if matched else "unavailable_no_bam_evidence"
                )
            ),
            "bam_primary_mapped_per_trim_read_proxy": fraction(primary, int(total_reads)) if matched else "",
            "bam_unique_gene_assigned_per_trim_read_proxy": fraction(unique_gene, int(total_reads)) if matched else "",
            "bam_mitochondrial_primary_fraction": fraction(mito, int(primary or 0)) if matched else "",
            "bam_rrna_primary_fraction": fraction(rrna, int(primary or 0)) if matched else "",
            "bam_biological_classification_status": class_status if matched else "unavailable_no_bam_evidence",
            "bam_matrix_conditional_saturation": sat,
            "bam_matrix_conditional_saturation_status": sat_status,
            "causal_interpretation": "correlation_only" if matched else "trimming_only_no_bam_evidence",
        }


def streaming_pearson(rows: Iterator[tuple[float, float]]) -> tuple[int, float | None]:
    n = 0
    mean_x = mean_y = sum_xx = sum_yy = covariance = 0.0
    for x, y in rows:
        n += 1
        delta_x = x - mean_x
        mean_x += delta_x / n
        delta_y = y - mean_y
        mean_y += delta_y / n
        sum_xx += delta_x * (x - mean_x)
        sum_yy += delta_y * (y - mean_y)
        covariance += delta_x * (y - mean_y)
    if n < 3 or sum_xx == 0 or sum_yy == 0:
        return n, None
    return n, covariance / math.sqrt(sum_xx * sum_yy)


def correlation_rows(connection: sqlite3.Connection) -> Iterator[dict[str, object]]:
    groups = list(connection.execute(
        "SELECT DISTINCT library,barcode_population FROM joined_metric ORDER BY 1,2"
    ))
    for library, population in groups:
        for x_name in TRIM_X_METRICS:
            for y_name in BAM_Y_METRICS:
                query = f"""
                    WITH pairs AS (
                        SELECT {x_name} AS x,{y_name} AS y
                        FROM joined_metric
                        WHERE library=? AND barcode_population=?
                          AND bam_match_status='matched'
                          AND {x_name} IS NOT NULL AND {y_name} IS NOT NULL
                    ), ranked AS (
                        SELECT
                          RANK() OVER (ORDER BY x)+(COUNT(*) OVER (PARTITION BY x)-1)/2.0 AS rx,
                          RANK() OVER (ORDER BY y)+(COUNT(*) OVER (PARTITION BY y)-1)/2.0 AS ry
                        FROM pairs
                    ) SELECT rx,ry FROM ranked
                """
                n, rho = streaming_pearson(
                    (float(x), float(y))
                    for x, y in connection.execute(query, (library, population))
                )
                yield {
                    "library": library, "barcode_population": population,
                    "trim_metric": x_name, "bam_metric": y_name, "n": n,
                    "spearman_rho": "" if rho is None else f"{rho:.8f}",
                    "interpretation": "association_only_not_causal",
                }


def trim_distribution_rows(connection: sqlite3.Connection) -> Iterator[dict[str, object]]:
    groups = connection.execute(
        "SELECT DISTINCT library,bam_match_status FROM joined_metric ORDER BY 1,2"
    ).fetchall()
    for library, status in groups:
        for metric in TRIM_X_METRICS:
            count, minimum, maximum, mean = connection.execute(
                f"SELECT COUNT({metric}),MIN({metric}),MAX({metric}),AVG({metric}) "
                "FROM joined_metric WHERE library=? AND bam_match_status=?",
                (library, status),
            ).fetchone()
            if not count:
                continue
            quantiles: list[float] = []
            for q in (0.0, 0.25, 0.5, 0.75, 1.0):
                offset = int(math.floor(q * (count - 1)))
                value = connection.execute(
                    f"SELECT {metric} FROM joined_metric "
                    f"WHERE library=? AND bam_match_status=? AND {metric} IS NOT NULL "
                    f"ORDER BY {metric} LIMIT 1 OFFSET ?",
                    (library, status, offset),
                ).fetchone()[0]
                quantiles.append(float(value))
            yield {
                "library": library, "bam_match_status": status, "metric": metric,
                "n": count, "minimum": minimum, "q25": quantiles[1],
                "median": quantiles[2], "q75": quantiles[3],
                "maximum": maximum, "mean": mean,
                "scope": "all_trimming_rows_including_no_bam_evidence",
            }


AGGREGATE_FIELDS = [
    "library", "cell_barcode", "barcode_population", "bam_match_status",
    "trim_row_count", "bam_matched_trim_rows", "bam_unmatched_trim_rows",
    "source_ids", "source_count", "read_mates", "read_mate_count",
    *COUNT_FIELDS, "frac_reads_with_adapter", "mean_bp_trimmed",
    "mean_original_length", "mean_final_length", *DETAIL_ADDITIONS[1:],
]


def run(args: argparse.Namespace) -> None:
    current = library_dirs(Path(args.current_root).resolve())
    output = Path(args.output_dir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".trim_bam_join_", dir=output) as temporary:
        connection = create_database(Path(temporary) / "join.sqlite")
        try:
            trim_fields, trim_count = read_trim_into_database(
                connection, Path(args.trim_by_source_tsv).resolve()
            )
            libraries = [
                row[0] for row in connection.execute(
                    "SELECT DISTINCT library FROM trim_row ORDER BY library"
                )
            ]
            for library in libraries:
                if library in current:
                    ingest_library_evidence(connection, library, current[library])

            counters: dict[str, int] = {}
            detail_fields = trim_fields + [
                field for field in DETAIL_ADDITIONS if field not in trim_fields
            ]
            write_tsv(
                output / "trim_cell_metrics.tsv.gz",
                detail_fields,
                detail_rows(connection, counters),
                compressed=True,
            )
            aggregate_path = output / "trim_cell_metrics_by_barcode.tsv.gz"
            write_tsv(
                aggregate_path, AGGREGATE_FIELDS, aggregate_rows(connection),
                compressed=True,
            )
            write_tsv(
                output / "trim_bam_correlations.tsv",
                [
                    "library", "barcode_population", "trim_metric", "bam_metric",
                    "n", "spearman_rho", "interpretation",
                ],
                correlation_rows(connection),
            )
            write_tsv(
                output / "trim_only_distributions.tsv",
                [
                    "library", "bam_match_status", "metric", "n", "minimum",
                    "q25", "median", "q75", "maximum", "mean", "scope",
                ],
                trim_distribution_rows(connection),
            )
            aggregate_count = connection.execute(
                "SELECT COUNT(*) FROM (SELECT 1 FROM trim_row GROUP BY library,barcode)"
            ).fetchone()[0]
            audit_rows = []
            for library in libraries:
                total, matched = connection.execute(
                    """
                    SELECT COUNT(*),SUM(CASE WHEN bam_match_status='matched' THEN 1 ELSE 0 END)
                    FROM joined_metric WHERE library=?
                    """, (library,),
                ).fetchone()
                barcodes = connection.execute(
                    "SELECT COUNT(*) FROM (SELECT 1 FROM trim_row WHERE library=? GROUP BY barcode)",
                    (library,),
                ).fetchone()[0]
                conflicts = connection.execute(
                    """
                    SELECT COUNT(*) FROM joined_metric j
                    JOIN trim_row t USING(row_id)
                    JOIN correction_conflict c
                      ON c.library=t.library AND c.barcode=t.barcode AND c.rg=t.source_id
                    WHERE j.library=? AND j.bam_match_status='matched'
                    """, (library,),
                ).fetchone()[0]
                audit_rows.append({
                    "release": RELEASE, "library": library, "trim_rows": total,
                    "matched_rows": matched or 0, "unmatched_rows": total-(matched or 0),
                    "barcode_aggregate_rows": barcodes,
                    "conflict_annotated_rows": conflicts,
                    "join_type": "left_join_all_trimming_rows",
                    "aggregate_contract": (
                        "counts_summed;fractions_recomputed_from_summed_counts;"
                        "mean_only_fields_total_reads_weighted;sources_and_mates_sorted_unique;"
                        "bam_barcode_totals_used_once;matched_rgs_sorted_unique"
                    ),
                })
            write_tsv(
                output / "trim_join_audit.tsv",
                [
                    "release", "library", "trim_rows", "matched_rows",
                    "unmatched_rows", "barcode_aggregate_rows",
                    "conflict_annotated_rows", "join_type", "aggregate_contract",
                ],
                audit_rows,
            )
            if sum(int(row["trim_rows"]) for row in audit_rows) != trim_count:
                raise AnalysisError("trim join audit row count did not reconcile")
            if sum(int(row["barcode_aggregate_rows"]) for row in audit_rows) != aggregate_count:
                raise AnalysisError("trim aggregate row count did not reconcile")
        finally:
            connection.close()
    print(
        f"Left-joined {trim_count} trimming rows; "
        f"{counters.get('matched', 0)} matched and "
        f"{counters.get('no_bam_evidence', 0)} retained without BAM evidence"
    )


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--current-root", required=True)
    result.add_argument("--trim-by-source-tsv", required=True)
    result.add_argument("--output-dir", required=True)
    return result


if __name__ == "__main__":
    try:
        run(parser().parse_args())
    except AnalysisError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)
