#!/usr/bin/env python3
"""Export observed, RG-prefix, and final-BAM thinning depth trajectories.

Each library is an independent shard. Per-step accumulators are fixed-width
vectors over the union of current and historical filtered-cell rosters; raw
barcode dictionaries are never built.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from array import array
from collections import defaultdict
from pathlib import Path
from typing import Iterable

from rna_evidence_common import (
    AnalysisError,
    barcodes,
    find_product,
    iter_nonzero_matrix_barcodes,
    iter_tsv,
    library_dirs,
    write_tsv,
)


def upper_median(values: Iterable[int]) -> int:
    ordered = sorted(values)
    return ordered[len(ordered) // 2] if ordered else 0


def fastq_inventory(path: Path | None) -> dict[tuple[str, str], int]:
    if path is None:
        return {}
    result: dict[tuple[str, str], int] = defaultdict(int)
    for row in iter_tsv(path):
        library = row.get("library") or row.get("Library") or row.get("sample")
        source = row.get("bp_id") or row.get("run_id") or row.get("Run")
        raw = next(
            (
                row.get(name)
                for name in ("raw_reads", "ReadsInRun", "reads", "read_pairs")
                if row.get(name)
            ),
            None,
        )
        if library and source and raw is not None:
            result[(library, source)] += int(float(raw))
    return dict(result)


def summary_read_depth(directory: Path) -> int | None:
    path = directory / "Summary.csv"
    if not path.is_file():
        return None
    with path.open("r", encoding="utf-8", newline="") as handle:
        values = {row[0]: row[1] for row in csv.reader(handle) if len(row) >= 2}
    for key in (
        "Number of Reads",
        "Number of Reads With Valid Barcodes",
        "Reads Mapped to Genome: Unique+Multiple",
    ):
        if key in values:
            try:
                return int(round(float(values[key])))
            except ValueError:
                continue
    return None


def category_for(current: bool, historical: bool) -> str:
    if current and historical:
        return "shared"
    if current:
        return "gained"
    return "lost"


def threshold_fields(
    counts: array,
    categories: list[str],
    thresholds: list[int],
) -> dict[str, int]:
    result: dict[str, int] = {}
    for threshold in thresholds:
        result[f"historical_barcodes_ge_{threshold}"] = sum(
            value >= threshold and category in {"shared", "lost"}
            for value, category in zip(counts, categories)
        )
        result[f"current_barcodes_ge_{threshold}"] = sum(
            value >= threshold and category in {"shared", "gained"}
            for value, category in zip(counts, categories)
        )
        result[f"shared_barcodes_ge_{threshold}"] = sum(
            value >= threshold and category == "shared"
            for value, category in zip(counts, categories)
        )
        result[f"gained_barcodes_ge_{threshold}"] = sum(
            value >= threshold and category == "gained"
            for value, category in zip(counts, categories)
        )
        result[f"lost_barcodes_ge_{threshold}"] = sum(
            value >= threshold and category == "lost"
            for value, category in zip(counts, categories)
        )
    return result


def matrix_vectors(
    directory: Path,
    roster_index: dict[str, int],
) -> tuple[array, array]:
    umi = array("Q", [0]) * len(roster_index)
    genes = array("I", [0]) * len(roster_index)
    for cb, count, gene_count in iter_nonzero_matrix_barcodes(
        find_product(directory, "raw/matrix.mtx.gz"),
        find_product(directory, "raw/barcodes.tsv.gz"),
    ):
        index = roster_index.get(cb)
        if index is not None:
            umi[index] = count
            genes[index] = gene_count
    return umi, genes


def bam_read_vector(
    directory: Path,
    roster_index: dict[str, int],
) -> tuple[array, bool]:
    result = array("Q", [0]) * len(roster_index)
    path = directory / "bam_evidence" / "barcode_read_metrics.tsv.gz"
    if not path.is_file():
        return result, False
    for row in iter_tsv(path, ("CB", "candidate_countedU_reads")):
        index = roster_index.get(row["CB"])
        if index is not None:
            result[index] = int(row["candidate_countedU_reads"])
    return result, True


def source_contract(
    evidence: Path,
) -> tuple[list[str], list[dict[str, object]]]:
    audit_path = evidence / "audit.json"
    try:
        audit = json.loads(audit_path.read_text(encoding="utf-8"))
        order = list(audit["source_order"])
    except (OSError, ValueError, KeyError, TypeError) as exc:
        raise AnalysisError(
            f"could not load explicit source order from {audit_path}"
        ) from exc
    if not order or len(order) > 64 or len(order) != len(set(order)):
        raise AnalysisError(f"invalid explicit source order in {audit_path}")
    source_index = {source: index for index, source in enumerate(order)}
    source_rows: dict[str, dict[str, object]] = {
        source: {
            "source_id": source,
            "source_index": index,
            "RGs": [],
            "primary_mapped_reads": 0,
        }
        for index, source in enumerate(order)
    }
    for row in iter_tsv(
        evidence / "rg_summary.tsv",
        ("RG", "source_id", "source_index", "primary_mapped_reads"),
    ):
        if row["RG"] == "__MISSING_OR_UNDECLARED__":
            continue
        source = row["source_id"]
        if source not in source_index:
            raise AnalysisError(f"unknown source {source!r} in RG summary")
        if int(row["source_index"]) != source_index[source]:
            raise AnalysisError(f"source index mismatch for {source}")
        target = source_rows[source]
        target["RGs"].append(row["RG"])
        target["primary_mapped_reads"] = (
            int(target["primary_mapped_reads"])
            + int(row["primary_mapped_reads"])
        )
    return order, [source_rows[source] for source in order]


def analyze_library(
    library: str,
    directory: Path,
    baseline: Path | None,
    inventory: dict[tuple[str, str], int],
    historical_inventory: dict[tuple[str, str], int],
    thresholds: list[int],
    random_steps: int,
) -> list[dict[str, object]]:
    evidence = directory / "bam_evidence"
    current_filtered = set(
        barcodes(find_product(directory, "filtered/barcodes.tsv.gz"))
    )
    historical_filtered = (
        set(barcodes(find_product(baseline, "filtered/barcodes.tsv.gz")))
        if baseline is not None else set()
    )
    fixed_roster = sorted(current_filtered | historical_filtered)
    roster_index = {cb: index for index, cb in enumerate(fixed_roster)}
    categories = [
        category_for(cb in current_filtered, cb in historical_filtered)
        for cb in fixed_roster
    ]
    current_indices = [
        roster_index[cb] for cb in fixed_roster if cb in current_filtered
    ]
    source_order, physical_sources = source_contract(evidence)
    source_count = len(source_order)
    roster_size = len(fixed_roster)

    # Fixed-width, fixed-roster vectors replace dict-per-prefix and
    # dict-per-thinning-fraction accumulation.
    prefix_umis = [array("Q", [0]) * roster_size for _ in range(source_count)]
    random_umis = [array("Q", [0]) * roster_size for _ in range(random_steps)]
    prefix_totals = array("Q", [0]) * source_count
    random_totals = array("Q", [0]) * random_steps
    for row in iter_tsv(
        evidence / "molecule_source_hash_bins.tsv.gz",
        (
            "CB", "source_mask_hex", "nested_min_hash_bin",
            "n_hash_bins", "candidate_matrix_molecules",
        ),
    ):
        mask = int(row["source_mask_hex"], 16)
        bin_index = int(row["nested_min_hash_bin"])
        bins = int(row["n_hash_bins"])
        count = int(row["candidate_matrix_molecules"])
        if mask <= 0 or mask & ~((1 << source_count) - 1):
            raise AnalysisError(
                f"invalid physical-source mask for {library}: "
                f"{row['source_mask_hex']}"
            )
        index = roster_index.get(row["CB"])
        for prefix in range(source_count):
            if mask & ((1 << (prefix + 1)) - 1):
                prefix_totals[prefix] += count
                if index is not None:
                    prefix_umis[prefix][index] += count
        for step in range(1, random_steps + 1):
            if bin_index < (bins * step) // random_steps:
                random_totals[step - 1] += count
                if index is not None:
                    random_umis[step - 1][index] += count

    source_reads = [array("Q", [0]) * roster_size for _ in range(source_count)]
    for row in iter_tsv(
        evidence / "barcode_rg_metrics.tsv.gz",
        ("CB", "source_index", "candidate_countedU_reads"),
    ):
        index = roster_index.get(row["CB"])
        source_index = int(row["source_index"])
        if not 0 <= source_index < source_count:
            raise AnalysisError(f"invalid source index for {library}")
        if index is not None:
            source_reads[source_index][index] += int(
                row["candidate_countedU_reads"]
            )
    prefix_reads = [array("Q", values) for values in source_reads]
    for source_index in range(1, source_count):
        previous = prefix_reads[source_index - 1]
        current = prefix_reads[source_index]
        for index in range(roster_size):
            current[index] += previous[index]

    final_reads, final_reads_available = bam_read_vector(
        directory, roster_index
    )
    current_umis, current_genes = matrix_vectors(directory, roster_index)
    if baseline is not None:
        old_umis, old_genes = matrix_vectors(baseline, roster_index)
        old_reads, old_reads_available = bam_read_vector(
            baseline, roster_index
        )
    else:
        old_umis = array("Q", [0]) * roster_size
        old_genes = array("I", [0]) * roster_size
        old_reads = array("Q", [0]) * roster_size
        old_reads_available = False

    rows: list[dict[str, object]] = []

    def add_point(
        trajectory_type: str,
        depth_label: str,
        umi_counts: array,
        read_counts: array | None,
        total_molecules: int,
        primary_reads: int | str,
        raw_reads: int | str,
        caller_cells: int | str,
        cell_value_type: str,
        read_value_type: str,
        source_prefix: str = "",
        random_fraction: str = "",
        gene_counts: array | None = None,
    ) -> None:
        row: dict[str, object] = {
            "library": library,
            "trajectory_type": trajectory_type,
            "depth_label": depth_label,
            "source_prefix": source_prefix,
            "random_fraction": random_fraction,
            "raw_input_reads": raw_reads,
            "primary_mapped_reads": primary_reads,
            "candidate_or_validated_molecules": total_molecules,
            "fixed_current_roster_size": len(current_indices),
            "fixed_current_roster_median_umis": upper_median(
                umi_counts[index] for index in current_indices
            ),
            "fixed_current_roster_median_candidate_reads": (
                upper_median(read_counts[index] for index in current_indices)
                if read_counts is not None else ""
            ),
            "fixed_current_roster_median_genes": (
                upper_median(gene_counts[index] for index in current_indices)
                if gene_counts is not None else ""
            ),
            "fixed_roster_read_value_type": read_value_type,
            "caller_cell_count": caller_cells,
            "cell_count_value_type": cell_value_type,
        }
        row.update(threshold_fields(umi_counts, categories, thresholds))
        rows.append(row)

    cumulative_primary = 0
    cumulative_raw = 0
    prefix_names: list[str] = []
    for index, source in enumerate(physical_sources):
        cumulative_primary += int(source["primary_mapped_reads"])
        cumulative_raw += inventory.get(
            (library, str(source["source_id"])), 0
        )
        prefix_names.append(str(source["source_id"]))
        add_point(
            "rg_prefix_threshold_proxy",
            f"prefix_{index + 1}",
            prefix_umis[index],
            prefix_reads[index],
            int(prefix_totals[index]),
            cumulative_primary,
            cumulative_raw or "",
            "",
            "posthoc_threshold_proxy_not_native_cell_call",
            "exact_ordinary_candidate_reads_in_rg_prefix",
            source_prefix="|".join(prefix_names),
        )

    final_primary = sum(
        int(source["primary_mapped_reads"]) for source in physical_sources
    )
    final_raw = sum(
        inventory.get((library, source), 0) for source in source_order
    )
    for step, counts in enumerate(random_umis, start=1):
        fraction = step / random_steps
        random_read_proxy = array(
            "Q", (int(round(value * fraction)) for value in final_reads)
        )
        add_point(
            "final_bam_random_thinning_threshold_proxy",
            f"final_bam_hash_{fraction:.3f}",
            counts,
            random_read_proxy if final_reads_available else None,
            int(random_totals[step - 1]),
            int(round(final_primary * fraction)),
            int(round(final_raw * fraction)) if final_raw else "",
            "",
            "posthoc_threshold_proxy_not_native_cell_call",
            (
                "linear_scaled_final_ordinary_candidate_read_proxy"
                if final_reads_available else "unavailable"
            ),
            random_fraction=f"{fraction:.6f}",
        )

    add_point(
        "observed_historical_cell_count",
        "current_native_mapping",
        current_umis,
        final_reads if final_reads_available else None,
        int(sum(current_umis)),
        final_primary,
        final_raw or summary_read_depth(directory) or "",
        len(current_filtered),
        "native_observed_STARsolo_cell_count",
        (
            "exact_final_ordinary_candidate_reads"
            if final_reads_available else "unavailable"
        ),
        gene_counts=current_genes,
    )
    if baseline is not None:
        historical_raw = sum(
            value
            for (item_library, _source), value in historical_inventory.items()
            if item_library == library
        )
        add_point(
            "observed_historical_cell_count",
            "historical_native_mapping",
            old_umis,
            old_reads if old_reads_available else None,
            int(sum(old_umis)),
            "",
            historical_raw or summary_read_depth(baseline) or "",
            len(historical_filtered),
            "native_observed_STARsolo_cell_count",
            (
                "exact_historical_ordinary_candidate_reads"
                if old_reads_available else "unavailable"
            ),
            gene_counts=old_genes,
        )
    return rows


def output_fields(thresholds: list[int]) -> list[str]:
    threshold_columns = [
        f"{category}_barcodes_ge_{threshold}"
        for threshold in thresholds
        for category in ("historical", "current", "shared", "gained", "lost")
    ]
    return [
        "library", "trajectory_type", "depth_label", "source_prefix",
        "random_fraction", "raw_input_reads", "primary_mapped_reads",
        "candidate_or_validated_molecules", "fixed_current_roster_size",
        "fixed_current_roster_median_umis",
        "fixed_current_roster_median_candidate_reads",
        "fixed_current_roster_median_genes",
        "fixed_roster_read_value_type", "caller_cell_count",
        "cell_count_value_type", *threshold_columns,
    ]


def run(args: argparse.Namespace) -> None:
    current = library_dirs(Path(args.current_root).resolve())
    baseline_dirs = (
        library_dirs(Path(args.baseline_root).resolve())
        if args.baseline_root else {}
    )
    inventory = fastq_inventory(
        Path(args.fastq_inventory) if args.fastq_inventory else None
    )
    historical_inventory = fastq_inventory(
        Path(args.historical_fastq_inventory)
        if getattr(args, "historical_fastq_inventory", None) else None
    )
    thresholds = sorted(set(args.thresholds))
    requested = set(getattr(args, "library", None) or [])
    selected = sorted(requested or current)
    unknown = requested - set(current)
    if unknown:
        raise AnalysisError(
            f"requested libraries are unavailable: {', '.join(sorted(unknown))}"
        )
    fields = output_fields(thresholds)
    output = Path(args.output_dir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    destination = output / "depth_trajectories.tsv"
    temporary = destination.with_name(f".{destination.name}.tmp.{os.getpid()}")
    row_count = 0
    try:
        with temporary.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=fields, delimiter="\t",
                extrasaction="ignore", lineterminator="\n",
            )
            writer.writeheader()
            for library in selected:
                rows = analyze_library(
                    library,
                    current[library],
                    baseline_dirs.get(library),
                    inventory,
                    historical_inventory,
                    thresholds,
                    args.random_steps,
                )
                library_output = (
                    output / "libraries" / library / "depth_trajectories.tsv"
                )
                write_tsv(library_output, fields, rows)
                for row in rows:
                    writer.writerow(row)
                    row_count += 1
        os.replace(temporary, destination)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise
    write_tsv(
        output / "depth_trajectory_definitions.tsv",
        ["trajectory_type", "interpretation"],
        [
            {
                "trajectory_type": "observed_historical_cell_count",
                "interpretation":
                    "actual retained STARsolo endpoint, observed depth, and native cell count",
            },
            {
                "trajectory_type": "rg_prefix_threshold_proxy",
                "interpretation":
                    "chronological physical-source evidence threshold curve conditional on final BAM; not a native lower-depth call",
            },
            {
                "trajectory_type":
                    "final_bam_random_thinning_threshold_proxy",
                "interpretation":
                    "nested minimum-QNAME-hash molecule threshold curve conditional on final BAM; not a native lower-depth call",
            },
        ],
    )
    print(
        f"Wrote {row_count} explicitly separated depth points from "
        f"{len(selected)} library shard(s) to {output}"
    )


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--current-root", required=True)
    result.add_argument("--baseline-root")
    result.add_argument("--fastq-inventory")
    result.add_argument("--historical-fastq-inventory")
    result.add_argument("--output-dir", required=True)
    result.add_argument(
        "--library", action="append",
        help="process only this library; repeat for independent shards",
    )
    result.add_argument(
        "--thresholds", type=int, nargs="+", default=[100, 500, 1000]
    )
    result.add_argument("--random-steps", type=int, default=10)
    return result


if __name__ == "__main__":
    try:
        arguments = parser().parse_args()
        if (
            arguments.random_steps < 1
            or any(value < 0 for value in arguments.thresholds)
        ):
            raise AnalysisError(
                "random steps must be positive and thresholds nonnegative"
            )
        run(arguments)
    except AnalysisError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)
