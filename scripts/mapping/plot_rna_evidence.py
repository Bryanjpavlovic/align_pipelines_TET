#!/usr/bin/env python3
"""Plot RNA evidence distributions and explicitly separated depth curves."""

from __future__ import annotations

import argparse
import math
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from analyze_gained_cells import (
    CORE_METRICS,
    COUNT_BINS,
    FRACTION_BINS,
    QUANTILES,
    TRIM_METRICS,
    Distribution,
)
from rna_evidence_common import AnalysisError, iter_tsv, read_tsv, write_tsv


def numeric(value: object) -> float | None:
    try:
        result = float(str(value))
        return result if math.isfinite(result) else None
    except (TypeError, ValueError):
        return None


def derive_distribution_tables(
    detail_path: Path,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Fallback for bounded detail fixtures; production uses all-row summaries."""
    accumulators: dict[tuple[str, str, str], Distribution] = {}
    for row in iter_tsv(detail_path, ("library", "CB", "cell_transition")):
        for metric in (*CORE_METRICS, *TRIM_METRICS):
            if numeric(row.get(metric)) is None:
                continue
            key = (row["library"], row["cell_transition"], metric)
            target = accumulators.setdefault(
                key,
                Distribution(
                    row["library"], row["cell_transition"], metric, 4096
                ),
            )
            target.add(row["CB"], row[metric])
    histograms = [
        row
        for target in accumulators.values()
        for row in target.histogram_rows()
    ]
    quantiles = [
        row
        for target in accumulators.values()
        for row in target.quantile_rows()
    ]
    return histograms, quantiles


def plot_category_distributions(
    detail_path: Path,
    histogram_path: Path | None,
    quantile_path: Path | None,
    output: Path,
) -> None:
    if histogram_path is not None and histogram_path.is_file():
        histograms = read_tsv(histogram_path)
        quantiles = (
            read_tsv(quantile_path)
            if quantile_path is not None and quantile_path.is_file()
            else []
        )
        scope = "all evidence-bearing rows"
    else:
        histograms, quantiles = derive_distribution_tables(detail_path)
        scope = "bounded detail rows"
    if not histograms:
        raise AnalysisError("cell-category histogram table contains no rows")
    write_tsv(
        output / "plot_values_cell_category_histograms.tsv",
        [
            "library", "cell_transition", "metric", "bin_index",
            "bin_left_inclusive", "bin_right_exclusive", "count",
            "population_n", "histogram_scope",
        ],
        histograms,
    )
    write_tsv(
        output / "plot_values_cell_category_quantiles.tsv",
        [
            "library", "cell_transition", "metric", "quantile", "value",
            "population_n", "sample_n", "method",
        ],
        quantiles,
    )

    aggregate: dict[tuple[str, str, int, str, str], int] = defaultdict(int)
    for row in histograms:
        key = (
            row["metric"],
            row["cell_transition"],
            int(row["bin_index"]),
            row["bin_left_inclusive"],
            row["bin_right_exclusive"],
        )
        aggregate[key] += int(row["count"])
    metrics = [
        metric
        for metric in (*CORE_METRICS, *TRIM_METRICS)
        if any(key[0] == metric for key in aggregate)
    ]
    categories = ("shared", "gained", "lost", "background")
    colors = {
        "shared": "#4C78A8",
        "gained": "#59A14F",
        "lost": "#E15759",
        "background": "#9C755F",
    }
    columns = 4
    rows = max(1, math.ceil(len(metrics) / columns))
    figure, axes = plt.subplots(
        rows, columns, figsize=(4.3 * columns, 3.5 * rows), squeeze=False
    )
    for axis, metric in zip(axes.flat, metrics):
        for category in categories:
            points = sorted(
                (
                    (key[2], key[3], key[4], count)
                    for key, count in aggregate.items()
                    if key[0] == metric and key[1] == category
                ),
                key=lambda value: value[0],
            )
            total = sum(point[3] for point in points)
            if not total:
                continue
            x = list(range(len(points)))
            y = [point[3] / total for point in points]
            axis.step(
                x, y, where="mid", label=category,
                color=colors[category], linewidth=1.6,
            )
        axis.set_title(metric.replace("_", " "), fontsize=9)
        axis.set_xlabel("fixed histogram bin")
        axis.set_ylabel("fraction")
        axis.set_yscale("log")
        axis.grid(alpha=0.2)
    for axis in list(axes.flat)[len(metrics):]:
        axis.axis("off")
    handles, labels = axes.flat[0].get_legend_handles_labels()
    if handles:
        figure.legend(handles, labels, loc="upper center", ncol=4)
    figure.suptitle(
        f"RNA evidence distributions by cell category ({scope})", y=1.002
    )
    figure.tight_layout()
    figure.savefig(output / "cell_category_distributions.png", dpi=180)
    plt.close(figure)


def depth_x(row: dict[str, str], fallback: int) -> float:
    for field in ("raw_input_reads", "primary_mapped_reads"):
        value = numeric(row.get(field))
        if value is not None and value > 0:
            return value
    fraction = numeric(row.get("random_fraction"))
    return fraction if fraction is not None else float(fallback)


def threshold_columns(rows: list[dict[str, str]]) -> list[str]:
    return sorted(
        {
            field
            for row in rows
            for field in row
            if "_barcodes_ge_" in field
        },
        key=lambda value: (
            value.split("_barcodes_ge_", 1)[0],
            int(value.rsplit("_", 1)[1]),
        ),
    )


def plot_depth(rows: list[dict[str, str]], output: Path) -> None:
    if not rows:
        raise AnalysisError("depth trajectory table contains no rows")
    fields = list(rows[0])
    write_tsv(output / "plot_values_depth_trajectories.tsv", fields, rows)
    observed_name = "observed_historical_cell_count"
    prefix_name = "rg_prefix_threshold_proxy"
    random_name = "final_bam_random_thinning_threshold_proxy"
    trajectories = (observed_name, prefix_name, random_name)
    titles = (
        "Observed historical cell count vs depth",
        "RG-prefix threshold/proxy cell counts",
        "Final-BAM random-thinning threshold/proxy cell counts",
    )
    thresholds = threshold_columns(rows)
    figure, axes = plt.subplots(1, 3, figsize=(18, 5.5))
    for axis, trajectory, title in zip(axes, trajectories, titles):
        by_library: dict[str, list[dict[str, str]]] = defaultdict(list)
        for row in rows:
            if row["trajectory_type"] == trajectory:
                by_library[row["library"]].append(row)
        for library, points in sorted(by_library.items()):
            if trajectory == observed_name:
                plotted = [
                    (depth_x(point, index + 1), numeric(point["caller_cell_count"]))
                    for index, point in enumerate(points)
                ]
                plotted = sorted(
                    item for item in plotted if item[1] is not None
                )
                if plotted:
                    axis.plot(
                        [item[0] for item in plotted],
                        [item[1] for item in plotted],
                        marker="o", label=library, linewidth=1.2,
                    )
            else:
                for field in thresholds:
                    if not field.startswith("current_barcodes_ge_"):
                        continue
                    plotted = [
                        (depth_x(point, index + 1), numeric(point.get(field)))
                        for index, point in enumerate(points)
                    ]
                    plotted = sorted(
                        item for item in plotted if item[1] is not None
                    )
                    if plotted:
                        threshold = field.rsplit("_", 1)[1]
                        axis.plot(
                            [item[0] for item in plotted],
                            [item[1] for item in plotted],
                            alpha=0.55, linewidth=1,
                            label=f"{library} ≥{threshold}",
                        )
        axis.set_title(title)
        axis.set_xlabel("raw input reads (primary reads/fraction fallback)")
        axis.set_ylabel(
            "native caller cells" if trajectory == observed_name
            else "fixed-roster barcodes above UMI threshold"
        )
        axis.set_xscale("symlog")
        axis.set_yscale("symlog")
        axis.grid(alpha=0.2)
    figure.suptitle(
        "Observed endpoints, source-prefix proxies, and thinning proxies are separate"
    )
    figure.tight_layout()
    figure.savefig(output / "depth_cell_count_curves.png", dpi=180)
    plt.close(figure)

    fixed_fields = [
        "library", "trajectory_type", "depth_label", "source_prefix",
        "random_fraction", "raw_input_reads", "primary_mapped_reads",
        "fixed_current_roster_size", "fixed_current_roster_median_umis",
        "fixed_current_roster_median_candidate_reads",
        "fixed_roster_read_value_type",
    ]
    fixed_rows = [
        {field: row.get(field, "") for field in fixed_fields} for row in rows
    ]
    write_tsv(
        output / "plot_values_fixed_roster_trajectories.tsv",
        fixed_fields,
        fixed_rows,
    )
    figure, axes = plt.subplots(1, 2, figsize=(13, 5))
    for axis, metric, label in (
        (axes[0], "fixed_current_roster_median_umis", "median UMIs"),
        (
            axes[1],
            "fixed_current_roster_median_candidate_reads",
            "median candidate reads",
        ),
    ):
        by_series: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
        for row in rows:
            by_series[(row["library"], row["trajectory_type"])].append(row)
        for (library, trajectory), points in sorted(by_series.items()):
            plotted = [
                (depth_x(point, index + 1), numeric(point.get(metric)))
                for index, point in enumerate(points)
            ]
            plotted = sorted(item for item in plotted if item[1] is not None)
            if plotted:
                axis.plot(
                    [item[0] for item in plotted],
                    [item[1] for item in plotted],
                    marker=".", linewidth=1,
                    label=f"{library} {trajectory}",
                )
        axis.set_xlabel("raw input reads (primary reads/fraction fallback)")
        axis.set_ylabel(label)
        axis.set_xscale("symlog")
        axis.set_yscale("symlog")
        axis.grid(alpha=0.2)
    figure.suptitle("Fixed current-roster UMI and read trajectories")
    figure.tight_layout()
    figure.savefig(output / "fixed_roster_trajectories.png", dpi=180)
    plt.close(figure)


def plot_sources(rows: list[dict[str, str]], output: Path) -> None:
    fields = [
        name for name in (
            "library", "RGs", "source_id", "source_order_index",
            "primary_mapped_reads", "unique_gene_tagged_reads",
            "nh_gt1_unique_gene_countedU_reads",
            "current_cell_associated_reads", "candidate_matrix_molecules",
            "marginal_candidate_molecules_first_observed",
        ) if rows and name in rows[0]
    ]
    write_tsv(output / "plot_values_source_yield.tsv", fields, rows)
    source_totals: dict[str, dict[str, int]] = defaultdict(
        lambda: defaultdict(int)
    )
    for row in rows:
        source = row["source_id"]
        for metric in (
            "primary_mapped_reads", "unique_gene_tagged_reads",
            "candidate_matrix_molecules",
            "marginal_candidate_molecules_first_observed",
        ):
            source_totals[source][metric] += int(row.get(metric, 0) or 0)
    sources = sorted(source_totals)
    figure, axes = plt.subplots(1, 2, figsize=(13, 5))
    axes[0].bar(
        sources,
        [source_totals[source]["unique_gene_tagged_reads"] for source in sources],
        color="#59A14F",
    )
    axes[0].set_title("Source-specific ordinary unique-gene reads")
    axes[1].bar(
        sources,
        [
            source_totals[source][
                "marginal_candidate_molecules_first_observed"
            ]
            for source in sources
        ],
        color="#F28E2B",
    )
    axes[1].set_title("Ordinary molecules first observed at source")
    for axis in axes:
        axis.tick_params(axis="x", rotation=40)
        axis.set_yscale("symlog")
    figure.tight_layout()
    figure.savefig(output / "source_yield_and_molecules.png", dpi=180)
    plt.close(figure)


LIBRARY_PANEL_ALIASES = (
    (
        "mapping",
        (
            "Reads Mapped to Genome: Unique+Multiple",
            "Reads Mapped to Genome: Unique + Multiple",
            "Mapping Rate",
        ),
    ),
    (
        "unique mapping",
        ("Reads Mapped to Genome: Unique", "Unique Mapping Rate"),
    ),
    (
        "gene assignment",
        (
            "Reads Mapped to GeneFull_Ex50pAS: Unique GeneFull_Ex50pAS",
            "Gene Assignment Rate",
        ),
    ),
    (
        "valid barcode",
        ("Reads With Valid Barcodes", "Valid Barcode Fraction"),
    ),
    (
        "cell associated",
        ("Fraction of Unique Reads in Cells", "Cell Associated Fraction"),
    ),
    ("sequencing saturation", ("Sequencing Saturation",)),
)


def plot_library_context(path: Path, output: Path) -> None:
    rows = read_tsv(path)
    if not rows:
        raise AnalysisError(f"processed statistics are empty: {path}")
    library_field = next(
        (
            name for name in ("Library", "Shortname", "Library_Number")
            if name in rows[0]
        ),
        None,
    )
    if library_field is None:
        raise AnalysisError(f"processed statistics lack a library field: {path}")
    resolved: list[tuple[str, str]] = []
    for label, aliases in LIBRARY_PANEL_ALIASES:
        field = next(
            (
                alias for alias in aliases
                if alias in rows[0]
                and any(numeric(row.get(alias)) is not None for row in rows)
            ),
            None,
        )
        if field is None:
            raise AnalysisError(
                f"processed statistics lack required {label} panel metric"
            )
        resolved.append((label, field))
    values = [
        {
            "library": row[library_field],
            "panel": label,
            "source_field": field,
            "value": row[field],
        }
        for row in rows
        for label, field in resolved
    ]
    write_tsv(
        output / "plot_values_library_mapping_context.tsv",
        ["library", "panel", "source_field", "value"],
        values,
    )
    libraries = [row[library_field] for row in rows]
    figure, axes = plt.subplots(2, 3, figsize=(16, 9), squeeze=False)
    for axis, (label, field) in zip(axes.flat, resolved):
        axis.bar(
            range(len(rows)),
            [numeric(row.get(field)) or 0 for row in rows],
            color="#76B7B2",
        )
        axis.set_title(label)
        axis.set_xticks(range(len(rows)), libraries, rotation=90, fontsize=6)
    figure.suptitle("STARsolo library-level mapping and cell-association context")
    figure.tight_layout()
    figure.savefig(output / "library_mapping_context.png", dpi=180)
    plt.close(figure)


def run(args: argparse.Namespace) -> None:
    output = Path(args.output_dir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    created = 0
    if args.cell_category_tsv:
        detail = Path(args.cell_category_tsv).resolve()
        histogram_arg = getattr(args, "cell_category_histograms_tsv", None)
        quantile_arg = getattr(args, "cell_category_quantiles_tsv", None)
        histogram = (
            Path(histogram_arg).resolve() if histogram_arg
            else detail.parent / "cell_category_histograms.tsv"
        )
        quantiles = (
            Path(quantile_arg).resolve() if quantile_arg
            else detail.parent / "cell_category_quantiles.tsv"
        )
        plot_category_distributions(
            detail,
            histogram if histogram.is_file() else None,
            quantiles if quantiles.is_file() else None,
            output,
        )
        created += 1
    if args.depth_trajectories_tsv:
        plot_depth(
            read_tsv(
                Path(args.depth_trajectories_tsv),
                ("library", "trajectory_type"),
            ),
            output,
        )
        created += 1
    if args.source_yield_tsv:
        plot_sources(
            read_tsv(
                Path(args.source_yield_tsv), ("library", "source_id")
            ),
            output,
        )
        created += 1
    if args.processedstats_tsv:
        plot_library_context(Path(args.processedstats_tsv), output)
        created += 1
    if not created:
        raise AnalysisError("select at least one input table to plot")
    print(
        f"Created {created} RNA evidence figure families and value tables in {output}"
    )


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--cell-category-tsv")
    result.add_argument("--cell-category-histograms-tsv")
    result.add_argument("--cell-category-quantiles-tsv")
    result.add_argument("--depth-trajectories-tsv")
    result.add_argument("--source-yield-tsv")
    result.add_argument("--processedstats-tsv")
    result.add_argument("--output-dir", required=True)
    return result


if __name__ == "__main__":
    try:
        run(parser().parse_args())
    except AnalysisError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)
