#!/usr/bin/env python3
"""Add BAM-derived conversion evidence to an authoritative raw FASTQ inventory."""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path

from rna_evidence_common import AnalysisError, read_tsv, write_tsv


SUM_FIELDS = (
    "all_records",
    "primary_mapped_reads",
    "nh1_reads",
    "unique_gene_tagged_reads",
    "current_cell_associated_reads",
    "candidate_countedU_reads",
    "candidate_matrix_molecules",
)


def ratio(numerator: int, denominator: float | int) -> str:
    return "" if denominator <= 0 else f"{numerator / denominator:.8f}"


def run(args: argparse.Namespace) -> None:
    inventory_path = Path(args.fastq_inventory)
    inventory = read_tsv(inventory_path)
    if not inventory:
        raise AnalysisError(f"raw FASTQ inventory is empty: {inventory_path}")
    for column in (args.library_column, args.source_column, args.raw_reads_column):
        if column not in inventory[0]:
            raise AnalysisError(f"raw FASTQ inventory lacks configured column {column!r}")
    source_rows = read_tsv(
        Path(args.source_yield_tsv),
        ("library", "source_id", "primary_mapped_reads", "candidate_matrix_molecules"),
    )
    totals: dict[tuple[str, str], dict[str, int]] = defaultdict(lambda: defaultdict(int))
    for row in source_rows:
        item = totals[(row["library"], row["source_id"])]
        for field in SUM_FIELDS:
            if row.get(field, ""):
                item[field] += int(row[field])
    output_rows: list[dict[str, object]] = []
    matched = 0
    for row in inventory:
        key = (row[args.library_column], row[args.source_column])
        bam = totals.get(key)
        raw_reads = float(row[args.raw_reads_column])
        out: dict[str, object] = dict(row)
        out["raw_fastq_representation_authority"] = "raw_fastq_inventory"
        out["bam_evidence_match"] = int(bam is not None)
        if bam is not None:
            matched += 1
            for field in SUM_FIELDS:
                out[f"bam_{field}"] = bam[field]
            out["bam_primary_mapped_per_raw_read"] = ratio(bam["primary_mapped_reads"], raw_reads)
            out["bam_unique_gene_assigned_per_raw_read"] = ratio(bam["unique_gene_tagged_reads"], raw_reads)
            out["bam_cell_associated_per_raw_read"] = ratio(bam["current_cell_associated_reads"], raw_reads)
            out["bam_candidate_molecules_per_raw_read"] = ratio(bam["candidate_matrix_molecules"], raw_reads)
        else:
            for field in SUM_FIELDS:
                out[f"bam_{field}"] = ""
            for field in (
                "bam_primary_mapped_per_raw_read", "bam_unique_gene_assigned_per_raw_read",
                "bam_cell_associated_per_raw_read", "bam_candidate_molecules_per_raw_read",
            ):
                out[field] = ""
        output_rows.append(out)
    fields = list(output_rows[0])
    output = Path(args.output_dir).resolve()
    write_tsv(output / "repooling_evidence.tsv", fields, output_rows)
    write_tsv(
        output / "repooling_evidence_audit.tsv",
        ["inventory_rows", "matched_rows", "unmatched_rows", "authority"],
        [{
            "inventory_rows": len(inventory), "matched_rows": matched,
            "unmatched_rows": len(inventory) - matched,
            "authority": "raw FASTQ counts remain authoritative for pool representation",
        }],
    )
    print(f"Added BAM conversion evidence to {matched}/{len(inventory)} raw inventory rows")


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--fastq-inventory", required=True)
    result.add_argument("--source-yield-tsv", required=True)
    result.add_argument("--output-dir", required=True)
    result.add_argument("--library-column", default="library")
    result.add_argument("--source-column", default="bp_id")
    result.add_argument("--raw-reads-column", default="raw_reads")
    return result


if __name__ == "__main__":
    try:
        run(parser().parse_args())
    except (AnalysisError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)
