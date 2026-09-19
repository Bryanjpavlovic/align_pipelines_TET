#!/usr/bin/env python3
"""
RNA mapping statistics visualization tool.
Generates focused plots from Summary.csv files across library folders and can
compare the current mapping statistics with a prior processedstats.tsv.
V10 - Longitudinal deltas plus profiler-backed cell-level read distributions

Run with: python plot_mapping_stats_V5.py --base-path ./
"""

import gzip
import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, Circle, Wedge
from matplotlib.gridspec import GridSpec
import matplotlib.patheffects as path_effects
import seaborn as sns
from pathlib import Path
import warnings
from scipy import stats
from scipy.interpolate import make_interp_spline
import argparse


RELEASE = "2026-09-05-v10-bam-evidence-layout"

# Try to import adjustText, set flag if not available
try:
    from adjustText import adjust_text
    HAS_ADJUST_TEXT = True
except ImportError:
    HAS_ADJUST_TEXT = False
    print("Note: adjustText not installed. Install with 'pip install adjustText' for better label placement.")

warnings.filterwarnings('ignore')

# Clean aesthetic setup
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'Liberation Sans']
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 11
plt.rcParams['axes.titlesize'] = 13
plt.rcParams['xtick.labelsize'] = 9
plt.rcParams['ytick.labelsize'] = 9
plt.rcParams['legend.fontsize'] = 9
plt.rcParams['figure.titlesize'] = 14
plt.rcParams['axes.spines.top'] = False
plt.rcParams['axes.spines.right'] = False
plt.rcParams['axes.grid'] = False
plt.rcParams['axes.edgecolor'] = '#666666'
plt.rcParams['axes.linewidth'] = 1.0


def get_text_color_for_background(value, vmin=0, vmax=1, cmap_name='RdYlGn'):
    """
    Calculate optimal text color (black or white) based on background luminance.
    Uses relative luminance formula for perceptual accuracy.
    """
    cmap = plt.cm.get_cmap(cmap_name)
    # Normalize value to 0-1 range
    norm_val = (value - vmin) / (vmax - vmin) if vmax != vmin else 0.5
    norm_val = np.clip(norm_val, 0, 1)
    
    # Get RGB from colormap
    rgba = cmap(norm_val)
    r, g, b = rgba[0], rgba[1], rgba[2]
    
    # Calculate relative luminance (WCAG formula)
    def to_linear(c):
        return c / 12.92 if c <= 0.03928 else ((c + 0.055) / 1.055) ** 2.4
    
    luminance = 0.2126 * to_linear(r) + 0.7152 * to_linear(g) + 0.0722 * to_linear(b)
    
    # Return white for dark backgrounds, black for light
    return 'white' if luminance < 0.4 else 'black'


def add_text_with_outline(ax, x, y, text, fontsize=9, fontweight='bold', 
                          text_color='white', outline_color='black', outline_width=2):
    """Add text with outline for universal readability on any background."""
    txt = ax.text(x, y, text, ha='center', va='center', 
                  fontsize=fontsize, fontweight=fontweight, color=text_color)
    txt.set_path_effects([
        path_effects.Stroke(linewidth=outline_width, foreground=outline_color),
        path_effects.Normal()
    ])
    return txt


class MappingStatsVisualizer:
    def __init__(
        self,
        base_path,
        output_dir="plots",
        interactive=False,
        baseline_stats=None,
        baseline_label="Baseline",
        current_label="Current",
        cell_reads_work_root=None,
        baseline_base_path=None,
        baseline_cell_reads_work_root=None,
    ):
        """Initialize visualizer with data path and output settings."""
        self.base_path = Path(base_path)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.interactive = interactive
        self.data = None
        self.baseline_stats = Path(baseline_stats) if baseline_stats else None
        self.baseline_data = None
        self.baseline_label = baseline_label
        self.current_label = current_label
        self.cell_reads_work_root = (
            Path(cell_reads_work_root) if cell_reads_work_root else None
        )
        self.baseline_base_path = (
            Path(baseline_base_path) if baseline_base_path else None
        )
        self.baseline_cell_reads_work_root = (
            Path(baseline_cell_reads_work_root)
            if baseline_cell_reads_work_root
            else None
        )
        self.cell_reads_data = None
        self.library_pattern = re.compile(r'Tet_2025_Multiome-RNA_(\d+)')
        
        # Define metric categories
        self.quality_metrics = [
            'Q30 Bases in CB+UMI',
            'Q30 Bases in RNA read',
            'Reads With Valid Barcodes'
        ]
        
        self.efficiency_metrics = [
            'Sequencing Saturation',
            'Fraction of Unique Reads in Cells',
            'Reads Mapped to Genome: Unique',
            'Reads Mapped to GeneFull_Ex50pAS: Unique GeneFull_Ex50pAS'
        ]
        
        self.cell_metrics = [
            'Estimated Number of Cells',
            'Median Reads per Cell',
            'Median UMI per Cell',
            'Median GeneFull_Ex50pAS per Cell'
        ]
        
        self.count_metrics = [
            'Number of Reads',
            'UMIs in Cells',
            'Total GeneFull_Ex50pAS Detected'
        ]
        
        # Color palettes
        self.colors = {
            'quality': ['#2E7D32', '#FDD835', '#E65100'],
            'good': '#00C853',
            'warning': '#FFB300', 
            'bad': '#D50000',
            'neutral': '#607D8B',
            'gradient': plt.cm.viridis,
            'diverging': plt.cm.RdBu_r,
            'categorical': plt.cm.tab10
        }
    
    def extract_library_number(self, lib_name):
        """Extract numeric ID from library name."""
        match = self.library_pattern.search(lib_name)
        return int(match.group(1)) if match else float('inf')
    
    def load_data(self):
        """Load and categorize the data."""
        data_list = []
        library_dirs = [d for d in self.base_path.iterdir() 
                       if d.is_dir() and d.name.startswith('Tet_')]
        
        if not library_dirs:
            print(f"No library directories found in {self.base_path}")
            return None
        
        print(f"Found {len(library_dirs)} library directories")
        
        for lib_dir in library_dirs:
            summary_file = lib_dir / "Summary.csv"
            if summary_file.exists():
                try:
                    df = pd.read_csv(summary_file, header=None, names=['Metric', 'Value'])
                    data_dict = dict(zip(df['Metric'], df['Value']))
                    data_dict['Library'] = lib_dir.name
                    data_dict['Library_Number'] = self.extract_library_number(lib_dir.name)
                    data_list.append(data_dict)
                    print(f"  Loaded: {lib_dir.name}")
                except Exception as e:
                    print(f"  Error with {lib_dir.name}: {e}")
        
        if not data_list:
            return None
        
        self.data = pd.DataFrame(data_list).sort_values('Library_Number')
        
        # Convert numeric columns
        for col in self.data.columns.difference(['Library', 'Library_Number']):
            self.data[col] = pd.to_numeric(self.data[col], errors='coerce')
        
        print(f"\nSuccessfully loaded {len(self.data)} libraries")
        
        # Get the range of expected library numbers
        self.min_lib_num = self.data['Library_Number'].min()
        self.max_lib_num = max(40, self.data['Library_Number'].max())
        
        return self.data

    def _load_processed_stats_table(self, path):
        """Load a processedstats.tsv into the same shape as Summary.csv data."""
        path = Path(path)
        table = pd.read_csv(path, sep='\t')
        if 'Library' not in table.columns:
            raise ValueError(f"Processed statistics file lacks Library column: {path}")

        if 'Shortname' in table.columns:
            library_numbers = pd.to_numeric(table['Shortname'], errors='coerce')
            table = table.drop(columns=['Shortname'])
        else:
            library_numbers = table['Library'].map(self.extract_library_number)

        if library_numbers.isna().any() or np.isinf(library_numbers).any():
            raise ValueError(f"Could not determine every library number in: {path}")

        table['Library_Number'] = library_numbers.astype(int)
        if table['Library_Number'].duplicated().any():
            duplicates = sorted(
                table.loc[table['Library_Number'].duplicated(False), 'Library_Number']
                .astype(int)
                .unique()
            )
            raise ValueError(
                f"Duplicate library number(s) in {path}: "
                + ', '.join(str(value) for value in duplicates)
            )

        for col in table.columns.difference(['Library', 'Library_Number']):
            table[col] = pd.to_numeric(table[col], errors='coerce')
        return table.sort_values('Library_Number').reset_index(drop=True)

    def load_processed_stats(self, path):
        """Use an existing processedstats.tsv as the current dataset."""
        self.data = self._load_processed_stats_table(path)
        self.min_lib_num = self.data['Library_Number'].min()
        self.max_lib_num = max(40, self.data['Library_Number'].max())
        print(f"Loaded current statistics for {len(self.data)} libraries from {path}")
        return self.data

    def load_baseline(self):
        """Load and validate the optional prior processedstats.tsv."""
        if self.baseline_stats is None:
            return None
        self.baseline_data = self._load_processed_stats_table(self.baseline_stats)
        current_libraries = set(self.data['Library_Number'])
        baseline_libraries = set(self.baseline_data['Library_Number'])
        common = sorted(current_libraries & baseline_libraries)
        if not common:
            raise ValueError(
                "The baseline and current statistics contain no common library numbers"
            )
        only_current = sorted(current_libraries - baseline_libraries)
        only_baseline = sorted(baseline_libraries - current_libraries)
        if only_current:
            print(
                "  Comparison excludes current-only libraries: "
                + ', '.join(str(value) for value in only_current)
            )
        if only_baseline:
            print(
                "  Comparison excludes baseline-only libraries: "
                + ', '.join(str(value) for value in only_baseline)
            )
        print(
            f"Loaded baseline statistics for {len(self.baseline_data)} libraries "
            f"from {self.baseline_stats}; {len(common)} libraries overlap"
        )
        return self.baseline_data

    @staticmethod
    def _open_text(path):
        """Open plain-text or gzip-compressed STARsolo output."""
        path = Path(path)
        if path.suffix == '.gz':
            return gzip.open(path, 'rt', encoding='utf-8')
        return path.open('r', encoding='utf-8')

    @staticmethod
    def _summary_values(path):
        """Return numeric values from a STARsolo Summary.csv."""
        table = pd.read_csv(path, header=None, names=['Metric', 'Value'])
        values = {}
        for metric, value in zip(table['Metric'], table['Value']):
            try:
                values[str(metric)] = float(value)
            except (TypeError, ValueError):
                continue
        return values

    def _index_cached_cell_reads(self, work_root):
        """Index CellReads.stats files retained in a Nextflow work directory."""
        index = {}
        if work_root is None:
            return index
        work_root = Path(work_root)
        if not work_root.is_dir():
            raise FileNotFoundError(
                f"Nextflow work directory does not exist: {work_root}"
            )
        for name in ('CellReads.stats', 'CellReads.stats.gz'):
            for path in work_root.rglob(name):
                solo_dir = next(
                    (parent for parent in path.parents if parent.name.endswith('Solo.out')),
                    None,
                )
                if solo_dir is None:
                    continue
                library = solo_dir.name[:-len('Solo.out')]
                index.setdefault(library, []).append(path)
        for paths in index.values():
            paths.sort(key=lambda path: path.stat().st_mtime, reverse=True)
        return index

    def _read_filtered_barcodes(self, library_dir):
        """Read STARsolo's called-cell barcode set for one library."""
        candidates = [
            library_dir / 'filtered' / 'barcodes.tsv.gz',
            library_dir / 'filtered' / 'barcodes.tsv',
        ]
        barcode_path = next((path for path in candidates if path.is_file()), None)
        if barcode_path is None:
            raise FileNotFoundError(
                f"Missing filtered/barcodes.tsv[.gz] below {library_dir}"
            )
        with self._open_text(barcode_path) as handle:
            barcodes = [line.rstrip('\r\n').split('\t', 1)[0] for line in handle]
        barcodes = [barcode for barcode in barcodes if barcode]
        if not barcodes:
            raise ValueError(f"No filtered barcodes found in {barcode_path}")
        if len(barcodes) != len(set(barcodes)):
            raise ValueError(f"Duplicate filtered barcodes found in {barcode_path}")
        return barcodes

    def _read_and_validate_cell_counts(self, stats_path, library_dir):
        """Load countedU for called cells and verify it against Summary.csv."""
        barcodes = self._read_filtered_barcodes(library_dir)
        counts = {barcode: None for barcode in barcodes}
        with self._open_text(stats_path) as handle:
            header = handle.readline().split()
            if not header:
                raise ValueError(f"Empty STARsolo cell-read file: {stats_path}")
            try:
                barcode_index = header.index('CB')
                counted_unique_index = header.index('countedU')
            except ValueError as exc:
                raise ValueError(
                    f"{stats_path} lacks required CB and countedU columns"
                ) from exc
            required_fields = max(barcode_index, counted_unique_index) + 1
            for line_number, line in enumerate(handle, start=2):
                fields = line.split()
                if len(fields) < required_fields:
                    raise ValueError(
                        f"Malformed {stats_path} line {line_number}: expected at least "
                        f"{required_fields} fields"
                    )
                barcode = fields[barcode_index]
                if barcode in counts:
                    if counts[barcode] is not None:
                        raise ValueError(
                            f"Duplicate called-cell barcode {barcode} in {stats_path}"
                        )
                    try:
                        counts[barcode] = int(fields[counted_unique_index])
                    except ValueError as exc:
                        raise ValueError(
                            f"Non-integer countedU in {stats_path} line {line_number}"
                        ) from exc

        missing = [barcode for barcode, value in counts.items() if value is None]
        if missing:
            example = ', '.join(missing[:3])
            raise ValueError(
                f"{stats_path} lacks {len(missing)} filtered barcode(s); examples: {example}"
            )

        values = [counts[barcode] for barcode in barcodes]
        summary_path = library_dir / 'Summary.csv'
        if not summary_path.is_file():
            raise FileNotFoundError(f"Missing STARsolo summary: {summary_path}")
        summary = self._summary_values(summary_path)
        unique_keys = [
            key for key in summary if key.startswith('Unique Reads in Cells Mapped to')
        ]
        required_summary = [
            'Estimated Number of Cells',
            'Median Reads per Cell',
        ]
        missing_summary = [key for key in required_summary if key not in summary]
        if len(unique_keys) != 1 or missing_summary:
            details = []
            if len(unique_keys) != 1:
                details.append(
                    f"expected one Unique Reads in Cells metric, found {len(unique_keys)}"
                )
            if missing_summary:
                details.append('missing ' + ', '.join(missing_summary))
            raise ValueError(f"Cannot validate {summary_path}: " + '; '.join(details))

        expected_cells = int(round(summary['Estimated Number of Cells']))
        expected_total = int(round(summary[unique_keys[0]]))
        expected_median = int(round(summary['Median Reads per Cell']))
        observed_total = int(sum(values))
        observed_median = int(sorted(values)[len(values) // 2])
        failures = []
        if len(values) != expected_cells:
            failures.append(f"cells {len(values):,} != {expected_cells:,}")
        if observed_total != expected_total:
            failures.append(f"read total {observed_total:,} != {expected_total:,}")
        if observed_median != expected_median:
            failures.append(f"median {observed_median:,} != {expected_median:,}")
        if failures:
            raise ValueError(
                f"{stats_path} does not reproduce {summary_path}: " + '; '.join(failures)
            )
        return barcodes, values

    def _load_cell_reads_snapshot(self, base_path, work_root, snapshot):
        """Load validated per-called-cell unique gene-assigned read counts."""
        base_path = Path(base_path)
        library_dirs = sorted(
            (
                path for path in base_path.iterdir()
                if path.is_dir() and self.library_pattern.search(path.name)
            ),
            key=lambda path: self.extract_library_number(path.name),
        )
        if not library_dirs:
            raise ValueError(f"No RNA library directories found in {base_path}")
        cached = self._index_cached_cell_reads(work_root)
        rows = []
        for library_dir in library_dirs:
            published = [
                library_dir / 'CellReads.stats.gz',
                library_dir / 'CellReads.stats',
                library_dir / 'bam_evidence' / 'CellReads.countedU.from_bam.tsv.gz',
                library_dir / 'CellReads.countedU.from_bam.tsv.gz',
            ]
            candidates = [path for path in published if path.is_file()]
            candidates.extend(cached.get(library_dir.name, []))
            if not candidates:
                raise FileNotFoundError(
                    f"No native or BAM-derived cell-read table found for "
                    f"{library_dir.name}; checked the "
                    f"published library directory and optional legacy Nextflow work "
                    f"cache. Historical mappings run without "
                    f"--soloCellReadStats Standard require the orchestrator's "
                    f"validated BAM-backfill mode."
                )
            errors = []
            for stats_path in candidates:
                try:
                    barcodes, values = self._read_and_validate_cell_counts(
                        stats_path, library_dir
                    )
                    break
                except Exception as exc:
                    errors.append(f"{stats_path}: {exc}")
            else:
                raise ValueError(
                    f"No valid cell-read table for {library_dir.name}:\n  "
                    + '\n  '.join(errors)
                )
            library_number = self.extract_library_number(library_dir.name)
            rows.extend(
                {
                    'Snapshot': snapshot,
                    'Library': library_dir.name,
                    'Shortname': int(library_number),
                    'Cell Barcode': barcode,
                    'Unique Gene-Assigned Reads': value,
                }
                for barcode, value in zip(barcodes, values)
            )
            print(
                f"  Validated cell reads: {library_dir.name} "
                f"({len(values):,} cells; median {sorted(values)[len(values) // 2]:,})"
            )
        return pd.DataFrame(rows)

    def load_cell_reads(self):
        """Load the current and optional baseline cell-level read distributions."""
        frames = [
            self._load_cell_reads_snapshot(
                self.base_path,
                self.cell_reads_work_root,
                self.current_label,
            )
        ]
        if self.baseline_base_path is not None:
            frames.insert(
                0,
                self._load_cell_reads_snapshot(
                    self.baseline_base_path,
                    self.baseline_cell_reads_work_root,
                    self.baseline_label,
                ),
            )
        self.cell_reads_data = pd.concat(frames, ignore_index=True)
        return self.cell_reads_data

    def export_cell_reads(self):
        """Export the validated observations plotted in the cell distribution."""
        if self.cell_reads_data is None:
            return
        output_file = self.output_dir / 'reads_per_cell.tsv.gz'
        self.cell_reads_data.to_csv(
            output_file,
            sep='\t',
            index=False,
            compression='gzip',
        )
        print("  Saved: reads_per_cell.tsv.gz")

    @staticmethod
    def _kde_sample(values, maximum=100000):
        """Return a deterministic, bounded sample for rendering a smooth KDE."""
        values = np.asarray(values, dtype=float)
        if len(values) <= maximum:
            return values
        indices = np.linspace(0, len(values) - 1, maximum, dtype=int)
        return np.sort(values)[indices]

    @staticmethod
    def _starsolo_median(values):
        """Match STARsolo's upper-middle median for even-sized cell sets."""
        ordered = np.sort(np.asarray(values, dtype=float))
        if len(ordered) == 0:
            return np.nan
        return float(ordered[len(ordered) // 2])

    def plot_reads_per_cell_distribution(self):
        """Plot called-cell read distributions and optional longitudinal change."""
        if self.cell_reads_data is None:
            raise ValueError('Cell-level read counts have not been loaded')
        value_column = 'Unique Gene-Assigned Reads'
        snapshots = list(dict.fromkeys(self.cell_reads_data['Snapshot']))
        current = self.cell_reads_data[
            self.cell_reads_data['Snapshot'] == self.current_label
        ]
        all_values = self.cell_reads_data[value_column].to_numpy(dtype=float)
        x_limit = max(52500, float(np.nanpercentile(all_values, 99.5)))
        bin_edges = np.linspace(0, x_limit, 90)
        colors = {
            self.baseline_label: '#78909C',
            self.current_label: '#00897B',
        }

        fig = plt.figure(figsize=(16, 10), facecolor='white')
        gs = GridSpec(2, 1, figure=fig, height_ratios=[1.15, 1], hspace=0.34)
        fig.suptitle(
            'Unique Gene-Assigned Reads per Called Cell',
            fontsize=18,
            fontweight='700',
            y=0.98,
        )

        ax = fig.add_subplot(gs[0])
        for snapshot in snapshots:
            values = self.cell_reads_data.loc[
                self.cell_reads_data['Snapshot'] == snapshot,
                value_column,
            ].to_numpy(dtype=float)
            color = colors.get(snapshot, '#5E35B1')
            ax.hist(
                values,
                bins=bin_edges,
                density=True,
                histtype='stepfilled',
                alpha=0.18,
                color=color,
                edgecolor=color,
                linewidth=1.1,
                label=f'{snapshot} histogram (n={len(values):,})',
            )
            kde_values = self._kde_sample(values)
            if len(np.unique(kde_values)) > 1:
                sns.kdeplot(
                    x=kde_values,
                    ax=ax,
                    color=color,
                    linewidth=2.3,
                    clip=(0, x_limit),
                    cut=0,
                    label=f'{snapshot} KDE',
                )
            median = self._starsolo_median(values)
            ax.axvline(median, color=color, linestyle='--', linewidth=1.8)
            ax.text(
                median,
                ax.get_ylim()[1] * (0.91 if snapshot == self.current_label else 0.80),
                f'{snapshot} median\n{median / 1000:.1f}k',
                color=color,
                ha='center',
                va='top',
                fontsize=9,
                fontweight='600',
            )
        ax.axvline(50000, color='#C62828', linestyle=':', linewidth=2)
        ax.text(
            50000,
            ax.get_ylim()[1] * 0.98,
            '50k target',
            color='#C62828',
            ha='center',
            va='top',
            fontsize=9,
            fontweight='600',
        )
        ax.set_xlim(0, x_limit)
        ax.set_xlabel('Unique GeneFull_Ex50pAS-assigned reads per called cell')
        ax.set_ylabel('Density')
        ax.set_title('Histogram and KDE across called cells', fontweight='600')
        ax.legend(loc='upper right', frameon=True, ncol=min(2, len(snapshots)))
        ax.grid(True, axis='y', alpha=0.22, linestyle=':')

        ax = fig.add_subplot(gs[1])
        if len(snapshots) == 2:
            baseline = self.cell_reads_data[
                self.cell_reads_data['Snapshot'] == self.baseline_label
            ]
            keys = ['Shortname', 'Cell Barcode']
            paired = baseline[keys + [value_column]].merge(
                current[keys + [value_column]],
                on=keys,
                how='inner',
                suffixes=('_baseline', '_current'),
                validate='one_to_one',
            )
            if paired.empty:
                raise ValueError(
                    'Baseline and current cell-read tables have no shared library/barcode pairs'
                )
            delta = (
                paired[f'{value_column}_current']
                - paired[f'{value_column}_baseline']
            ).to_numpy(dtype=float)
            limit = float(np.nanpercentile(np.abs(delta), 99.5))
            limit = max(limit, 1.0)
            ax.hist(
                delta,
                bins=np.linspace(-limit, limit, 100),
                density=True,
                color='#5E35B1',
                alpha=0.25,
                edgecolor='#4527A0',
                linewidth=0.8,
                label=f'Histogram (n={len(delta):,} shared cells)',
            )
            kde_values = self._kde_sample(delta)
            if len(np.unique(kde_values)) > 1:
                sns.kdeplot(
                    x=kde_values,
                    ax=ax,
                    color='#4527A0',
                    linewidth=2.3,
                    clip=(-limit, limit),
                    cut=0,
                    label='KDE',
                )
            median_delta = self._starsolo_median(delta)
            ax.axvline(0, color='#424242', linewidth=1.2)
            ax.axvline(median_delta, color='#4527A0', linestyle='--', linewidth=2)
            ax.set_xlim(-limit, limit)
            ax.set_xlabel(
                f'{self.current_label} − {self.baseline_label} reads for the same called cell'
            )
            ax.set_ylabel('Density')
            ax.set_title(
                f'Same-cell gain distribution · median {median_delta / 1000:+.1f}k reads',
                fontweight='600',
            )
            ax.legend(loc='upper right', frameon=True)
            ax.grid(True, axis='y', alpha=0.22, linestyle=':')
        else:
            summaries = []
            for library_number, group in current.groupby('Shortname', sort=True):
                values = group[value_column].to_numpy(dtype=float)
                percentiles = np.percentile(values, [10, 25, 75, 90])
                summaries.append(
                    (
                        int(library_number),
                        percentiles[0],
                        percentiles[1],
                        self._starsolo_median(values),
                        percentiles[2],
                        percentiles[3],
                    )
                )
            summaries = np.asarray(summaries, dtype=float)
            positions = np.arange(len(summaries))
            ax.vlines(
                positions,
                summaries[:, 1],
                summaries[:, 5],
                color='#80CBC4',
                linewidth=2,
                label='10th–90th percentile',
            )
            ax.vlines(
                positions,
                summaries[:, 2],
                summaries[:, 4],
                color='#00897B',
                linewidth=7,
                label='Interquartile range',
            )
            ax.scatter(
                positions,
                summaries[:, 3],
                color='#004D40',
                s=24,
                zorder=3,
                label='Median',
            )
            ax.axhline(50000, color='#C62828', linestyle=':', linewidth=2, label='50k target')
            ax.set_xticks(positions)
            ax.set_xticklabels(summaries[:, 0].astype(int), fontsize=8)
            ax.set_xlabel('Library number')
            ax.set_ylabel('Unique gene-assigned reads per called cell')
            ax.set_title('Cell-level distribution by library', fontweight='600')
            ax.legend(loc='upper right', frameon=True, ncol=2)
            ax.grid(True, axis='y', alpha=0.22, linestyle=':')

        fig.text(
            0.5,
            0.015,
            'Histograms use all called cells; KDE curves use a deterministic maximum of 100,000 cells per snapshot. '
            'Axes stop at the 99.5th percentile for readability.',
            ha='center',
            va='bottom',
            fontsize=9,
            color='#616161',
        )
        plt.tight_layout(rect=[0, 0.035, 1, 0.955])
        self._save_and_show('reads_per_cell_distribution')
    
    def plot_quality_control_dashboard(self):
        """QC dashboard with expanded saturation plot (3 rows)."""
        fig = plt.figure(figsize=(16, 18), facecolor='white')
        
        # 4-row layout: Performance Matrix, Complexity, Saturation, Table
        # Increased hspace from 0.3 to 0.45 for more breathing room between plots
        gs = GridSpec(4, 1, figure=fig, height_ratios=[1, 1.2, 1, 0.5], hspace=0.45)
        
        # Title - lowered from y=0.98 to y=0.965 to reduce top gap
        fig.suptitle('RNA Mapping Quality Control Dashboard', 
                    fontsize=18, fontweight='700', y=0.965, color='#1a1a1a')
        
        # Row 1: Performance Matrix (full width) - NOW WITH ABSOLUTE VALUES
        ax1 = fig.add_subplot(gs[0])
        self._create_performance_heatmap_absolute(ax1, self.data)
        
        # Row 2: Library Complexity Analysis (full width)
        if 'Mean UMI per Cell' in self.data.columns and 'Mean GeneFull_Ex50pAS per Cell' in self.data.columns:
            ax2 = fig.add_subplot(gs[1])
            self._create_complexity_scatter(ax2, self.data)
        
        # Row 3: Saturation Distribution
        if 'Sequencing Saturation' in self.data.columns:
            ax3 = fig.add_subplot(gs[2])
            ax_table = fig.add_subplot(gs[3])
            self._create_saturation_plot_expanded(ax3, ax_table, self.data)
        
        # Adjusted rect: top margin reduced from 0.96 to 0.94 to use more vertical space
        plt.tight_layout(rect=[0, 0.02, 1, 0.94])
        self._save_and_show('quality_control_dashboard')
    
    def plot_cell_quality_matrix(self):
        """Matrix plot showing cell-related metrics relationships."""
        cell_metrics = ['Median Reads per Cell', 'Median UMI per Cell',
                       'Median GeneFull_Ex50pAS per Cell', 'Estimated Number of Cells']
        available = [m for m in cell_metrics if m in self.data.columns]
        
        if len(available) < 2:
            return
        
        fig, axes = plt.subplots(1, 3, figsize=(16, 5), facecolor='white')
        
        # 1. Correlation heatmap - improved contrast logic
        ax = axes[0]
        corr_matrix = self.data[available].corr()
        
        im = ax.imshow(corr_matrix, cmap='RdBu_r', vmin=-1, vmax=1, aspect='auto')
        
        # Add correlation values with proper contrast
        for i in range(len(available)):
            for j in range(len(available)):
                val = corr_matrix.iloc[i, j]
                # Use luminance-based contrast
                color = get_text_color_for_background(val, vmin=-1, vmax=1, cmap_name='RdBu_r')
                ax.text(j, i, f'{val:.2f}', ha='center', va='center',
                       color=color, fontweight='600', fontsize=11)
        
        ax.set_xticks(range(len(available)))
        ax.set_yticks(range(len(available)))
        
        x_labels = [m.replace(' per Cell', '').replace('Mean ', '').replace('Estimated Number of', 'N') 
                   for m in available]
        ax.set_xticklabels(x_labels, rotation=45, ha='right', fontsize=9)
        
        y_labels = [m.replace(' per Cell', '\nper Cell').replace('Estimated Number of', 'Number of') 
                   for m in available]
        ax.set_yticklabels(y_labels, fontsize=9)
        
        ax.set_title('Metric Correlations', fontweight='600', pad=10)
        cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label('Correlation', rotation=270, labelpad=15)
        
        # 2. UMI vs Reads efficiency plot - IMPROVED: green=good, red=bad duplication
        if 'Median Reads per Cell' in available and 'Median UMI per Cell' in available:
            ax = axes[1]
            x = self.data['Median Reads per Cell']
            y = self.data['Median UMI per Cell']
            
            # Calculate duplication rate (lower is better)
            dup_rate = x / y
            
            # Use RdYlGn_r: green for LOW duplication (good), red for HIGH (bad)
            scatter = ax.scatter(x, y, c=dup_rate, s=120, cmap='RdYlGn_r', 
                               edgecolors='#333333', linewidth=1.5, alpha=0.9,
                               vmin=dup_rate.min(), vmax=dup_rate.max())
            
            # Add trend line
            z = np.polyfit(x, y, 1)
            p = np.poly1d(z)
            x_trend = np.linspace(x.min(), x.max(), 100)
            ax.plot(x_trend, p(x_trend), "--", color='#555555', alpha=0.7, linewidth=2)
            
            ax.set_xlabel('Median Reads per Cell', fontweight='500')
            ax.set_ylabel('Median UMI per Cell', fontweight='500')
            ax.set_title('Sequencing Depth vs Unique Molecules', fontweight='600', pad=10)
            # Removed grid for cleaner look
            
            cbar = plt.colorbar(scatter, ax=ax, fraction=0.046, pad=0.04)
            cbar.set_label('Duplication Rate\n(Reads/UMI)\n← Better | Worse →', rotation=270, labelpad=40)
        
        # 3. Cell count vs complexity - IMPROVED: single color with outlined text
        if 'Estimated Number of Cells' in available and 'Median GeneFull_Ex50pAS per Cell' in available:
            ax = axes[2]
            x = self.data['Estimated Number of Cells']
            y = self.data['Median GeneFull_Ex50pAS per Cell']
            
            # Single color scheme - teal/cyan for all points
            scatter = ax.scatter(x, y, s=200, c='#00838F', 
                               edgecolors='white', linewidth=2, alpha=0.9)
            
            # Annotate with outlined text for universal readability
            for i, txt in enumerate(self.data['Library_Number']):
                add_text_with_outline(ax, x.iloc[i], y.iloc[i], str(txt),
                                     fontsize=9, fontweight='bold',
                                     text_color='white', outline_color='#004D40', outline_width=2)
            
            ax.set_xlabel('Estimated Number of Cells', fontweight='500')
            ax.set_ylabel('Median Genes per Cell', fontweight='500')
            ax.set_title('Cell Recovery vs Complexity', fontweight='600', pad=10)
            # Removed grid for cleaner look
        
        fig.suptitle('Cell Quality Analysis Matrix', fontsize=15, fontweight='600', y=1.02)
        plt.tight_layout()
        self._save_and_show('cell_quality_matrix')
    
    def plot_outlier_detection_report(self):
        """Identify and visualize outlier libraries across multiple metrics."""
        metrics = self.quality_metrics + self.efficiency_metrics + ['Median UMI per Cell']
        available = [m for m in metrics if m in self.data.columns]
        
        if len(available) < 2:
            return
        
        # Calculate outlier scores
        outlier_scores = pd.DataFrame()
        outlier_details = {}
        
        for metric in available:
            values = self.data[metric].values
            z_scores = np.abs(stats.zscore(values))
            outlier_scores[metric] = z_scores
            
            # Find outliers (z-score > 2)
            outlier_mask = z_scores > 2
            if outlier_mask.any():
                outlier_libs = self.data.loc[outlier_mask, 'Library_Number'].values
                outlier_vals = values[outlier_mask]
                outlier_details[metric] = list(zip(outlier_libs, outlier_vals, z_scores[outlier_mask]))
        
        # Create visualization
        fig = plt.figure(figsize=(16, 10), facecolor='white')
        gs = GridSpec(3, 3, figure=fig, hspace=0.3, wspace=0.3)
        
        # 1. Outlier heatmap - IMPROVED: better colormap with more contrast at low end
        ax = fig.add_subplot(gs[:2, :2])
        
        matrix = outlier_scores.T
        lib_nums = self.data['Library_Number'].values
        
        # Custom colormap: light gray -> yellow -> orange -> red
        # This gives better contrast for low values
        from matplotlib.colors import LinearSegmentedColormap
        colors_list = ['#F5F5F5', '#FFF9C4', '#FFE082', '#FFB74D', '#FF8A65', '#E53935']
        cmap_outlier = LinearSegmentedColormap.from_list('outlier', colors_list)
        
        im = ax.imshow(matrix, cmap=cmap_outlier, aspect='auto', vmin=0, vmax=3)
        
        ax.set_xticks(np.arange(len(lib_nums)))
        ax.set_yticks(np.arange(len(available)))
        ax.set_xticklabels([str(n) for n in lib_nums], rotation=0, ha='center', fontsize=8)
        ax.set_xlabel('Library Number', fontsize=10, fontweight='500')
        
        # Shortened y-axis labels
        y_labels = []
        for m in available:
            if 'Q30 Bases' in m:
                y_labels.append(m.replace('Q30 Bases in ', 'Q30 '))
            elif 'Reads Mapped to' in m:
                y_labels.append(m.replace('Reads Mapped to GeneFull_Ex50pAS: Unique GeneFull_Ex50pAS', 'Mapped to Genes')
                              .replace('Reads Mapped to Genome: Unique', 'Mapped to Genome'))
            elif 'Mean' in m:
                y_labels.append(m.replace('Mean ', '').replace(' per Cell', '/Cell'))
            else:
                y_labels.append(m.replace('Fraction of ', ''))
        
        ax.set_yticklabels(y_labels, fontsize=9)
        
        # Show ALL z-scores, but bold significant ones
        for i in range(len(available)):
            for j in range(len(lib_nums)):
                val = matrix.iloc[i, j]
                if val > 2:
                    # Significant outlier - bold white text
                    ax.text(j, i, f'{val:.1f}', ha='center', va='center',
                           color='white', fontsize=8, fontweight='bold')
                elif val > 1:
                    # Mild concern - show value in dark text
                    ax.text(j, i, f'{val:.1f}', ha='center', va='center',
                           color='#333333', fontsize=7, fontweight='normal')
        
        ax.set_title('Outlier Detection Heatmap\n(Z-scores, * = |z| > 2)', 
                    fontsize=12, fontweight='600', pad=10)
        
        cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label('|Z-score|', rotation=270, labelpad=15)
        
        # 2. Outlier summary (top right)
        ax = fig.add_subplot(gs[0, 2])
        
        outliers_per_lib = (outlier_scores > 2).sum(axis=1)
        problem_libs = outliers_per_lib[outliers_per_lib > 0].sort_values(ascending=False)
        
        if len(problem_libs) > 0:
            top_problems = problem_libs.head(5)
            y_pos = np.arange(len(top_problems))
            
            bars = ax.barh(y_pos, top_problems.values, 
                         color=plt.cm.Reds(0.3 + 0.7 * top_problems.values / max(top_problems.max(), 1)),
                         edgecolor='#8B0000', linewidth=1)
            
            ax.set_yticks(y_pos)
            ax.set_yticklabels([f'Lib {self.data.iloc[i]["Library_Number"]}' 
                               for i in top_problems.index], fontsize=9)
            ax.set_xlabel('Number of Outlier Metrics', fontsize=10)
            ax.set_title('Most Problematic Libraries', fontsize=11, fontweight='600')
            ax.grid(axis='x', alpha=0.3, linestyle=':')
            
            for bar, val in zip(bars, top_problems.values):
                ax.text(val + 0.1, bar.get_y() + bar.get_height()/2, 
                       f'{int(val)}', va='center', fontsize=9, fontweight='bold')
        else:
            ax.text(0.5, 0.5, 'No significant\noutliers detected', 
                   ha='center', va='center', fontsize=12,
                   color='#2E7D32', fontweight='bold', transform=ax.transAxes)
            ax.set_title('Outlier Summary', fontsize=11, fontweight='600')
        
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        # 3. Outlier details table (middle right)
        ax = fig.add_subplot(gs[1, 2])
        ax.axis('tight')
        ax.axis('off')
        
        if outlier_details:
            text_lines = ['Key Outliers Detected:\n' + '='*25]
            for metric, details in list(outlier_details.items())[:3]:
                # Shorten metric name
                short_name = metric
                if 'Q30 Bases' in metric:
                    short_name = metric.replace('Q30 Bases in ', 'Q30 ')
                text_lines.append(f'\n{short_name}:')
                for lib, val, z in details[:2]:
                    direction = 'high' if z > 0 else 'low'
                    text_lines.append(f'  Lib {lib}: {val:.3f} ({direction}, z={z:.1f})')
            
            ax.text(0.05, 0.95, '\n'.join(text_lines), transform=ax.transAxes,
                   fontsize=9, verticalalignment='top', fontfamily='monospace',
                   bbox=dict(boxstyle='round,pad=0.5', facecolor='#FFF3E0', edgecolor='#E65100', linewidth=1.5))
        else:
            ax.text(0.5, 0.5, 'No outliers to report', 
                   ha='center', va='center', transform=ax.transAxes,
                   fontsize=11, color='#2E7D32', fontweight='bold')
        
        # 4. Global distribution (bottom panel) - IMPROVED: cleaner, no gray placeholders
        ax = fig.add_subplot(gs[2, :])
        
        overall_scores = outlier_scores.mean(axis=1)
        lib_nums_actual = self.data['Library_Number'].values
        
        # Only plot actual libraries, no gray placeholders
        colors = ['#4CAF50' if s < 1 else '#FF9800' if s < 2 else '#D32F2F' 
                 for s in overall_scores]
        
        bars = ax.bar(lib_nums_actual, overall_scores, color=colors, edgecolor='white', 
                      linewidth=1, alpha=0.9, width=0.8)
        
        ax.axhline(1, color='#FF9800', linestyle='--', linewidth=2, alpha=0.8, 
                  label='Mild concern (z=1)')
        ax.axhline(2, color='#D32F2F', linestyle='--', linewidth=2, alpha=0.8, 
                  label='Significant outlier (z=2)')
        
        ax.set_xticks(lib_nums_actual)
        ax.set_xticklabels([str(n) for n in lib_nums_actual], rotation=0, fontsize=8)
        
        ax.set_xlabel('Library Number', fontsize=11, fontweight='500')
        ax.set_ylabel('Mean |Z-score| Across Metrics', fontsize=11, fontweight='500')
        ax.set_title('Overall Outlier Score by Library', fontsize=11, fontweight='600')
        ax.legend(loc='upper right', frameon=True, fontsize=9, framealpha=0.95)
        ax.grid(True, axis='y', alpha=0.3, linestyle=':')
        ax.set_xlim(lib_nums_actual.min() - 0.5, lib_nums_actual.max() + 0.5)
        
        fig.suptitle('Outlier Detection & Quality Assurance Report', 
                    fontsize=15, fontweight='600', y=0.98)
        
        plt.tight_layout()
        self._save_and_show('outlier_detection_report')
    
    def _create_saturation_plot_expanded(self, ax, ax_table, data):
        """Create expanded saturation bar plot with median reads indicator and below-target table."""
        saturation = data['Sequencing Saturation'].values
        lib_nums = data['Library_Number'].values
        reads_per_cell = data['Median Reads per Cell'].values if 'Median Reads per Cell' in data.columns else None
        n_cells = data['Estimated Number of Cells'].values if 'Estimated Number of Cells' in data.columns else None
        
        # Sort by saturation for visualization
        sorted_idx = np.argsort(saturation)
        saturation_sorted = saturation[sorted_idx]
        lib_nums_sorted = lib_nums[sorted_idx]
        reads_sorted = reads_per_cell[sorted_idx] if reads_per_cell is not None else None
        cells_sorted = n_cells[sorted_idx] if n_cells is not None else None
        
        # Create color gradient based on saturation level
        colors = ['#C62828' if s < 0.5 else '#F57C00' if s < 0.7 else '#388E3C' 
                 for s in saturation_sorted]
        
        bar_width = 0.7
        x_positions = np.arange(len(saturation_sorted))
        
        # Bar plot
        bars = ax.bar(x_positions, saturation_sorted, width=bar_width, color=colors, 
                      edgecolor='#333333', linewidth=1, alpha=0.9)
        
        # Add threshold line
        ax.axhline(0.7, color='#333333', linestyle='--', linewidth=2.5, alpha=0.8)
        
        # Fixed y-position for all read labels
        READS_LABEL_Y = 0.35
        
        for i, (bar, lib) in enumerate(zip(bars, lib_nums_sorted)):
            # Bold library number above bar
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
                   str(lib), ha='center', va='bottom', fontsize=9, fontweight='bold',
                   color='#1a1a1a')
            
            # Median reads indicator inside bar at FIXED position
            if reads_sorted is not None:
                reads_k = reads_sorted[i] / 1000
                reads_label = f'{reads_k:.0f}k' if reads_k >= 1 else f'{reads_k:.1f}k'
                
                ax.text(bar.get_x() + bar.get_width()/2, READS_LABEL_Y,
                       reads_label, ha='center', va='center', fontsize=7,
                       color='black', fontweight='500',
                       bbox=dict(boxstyle='round,pad=0.15', facecolor='white', 
                                alpha=0.9, edgecolor='none'))
        
        ax.set_xticks([])
        ax.set_xlabel('Libraries Sorted by Saturation', fontsize=12, fontweight='500')
        ax.set_ylabel('Sequencing Saturation', fontsize=12, fontweight='500')
        ax.set_title('Sequencing Saturation Distribution\n(values in bars = median reads/cell)', 
                    fontsize=14, fontweight='600')
        ax.set_ylim(0, 1.05)
        
        ax.legend(['Typical target: 0.7'], loc='lower right', frameon=True, fontsize=10,
                 fancybox=True, framealpha=0.95, edgecolor='#333333')
        ax.grid(True, axis='y', alpha=0.3, linestyle=':')
        
        # --- Create below-target table ---
        ax_table.axis('off')
        
        # Find libraries below 0.7 target
        below_target_mask = saturation < 0.7
        
        if below_target_mask.any():
            below_libs = lib_nums[below_target_mask]
            below_sat = saturation[below_target_mask]
            below_cells = n_cells[below_target_mask] if n_cells is not None else None
            below_reads = reads_per_cell[below_target_mask] if reads_per_cell is not None else None
            
            # Sort by library number
            sort_by_num = np.argsort(below_libs)
            below_libs = below_libs[sort_by_num]
            below_sat = below_sat[sort_by_num]
            if below_cells is not None:
                below_cells = below_cells[sort_by_num]
            if below_reads is not None:
                below_reads = below_reads[sort_by_num]
            
            col_labels = ['Library', 'Saturation', 'Total Cells', 'Median Reads/Cell']
            
            table_data = []
            for i in range(len(below_libs)):
                row = [
                    str(int(below_libs[i])),
                    f'{below_sat[i]:.3f}',
                    f'{int(below_cells[i]):,}' if below_cells is not None else 'N/A',
                    f'{below_reads[i]/1000:.1f}k' if below_reads is not None else 'N/A'
                ]
                table_data.append(row)
            
            table = ax_table.table(
                cellText=table_data,
                colLabels=col_labels,
                loc='center',
                cellLoc='center',
                colWidths=[0.15, 0.15, 0.15, 0.18]
            )
            
            table.auto_set_font_size(False)
            table.set_fontsize(10)
            table.scale(1.2, 1.6)
            
            for key, cell in table.get_celld().items():
                cell.set_edgecolor('#666666')
                cell.set_linewidth(1)
                row, col = key
                if row == 0:
                    cell.set_text_props(fontweight='bold', color='white')
                    cell.set_facecolor('#C62828')
                else:
                    cell.set_facecolor('#FFF8F8' if row % 2 == 0 else 'white')
            
            ax_table.set_title('Libraries Below Target Saturation (< 0.7)', 
                              fontsize=12, fontweight='600', color='#C62828', pad=30)
        else:
            ax_table.text(0.5, 0.5, '✓ All libraries meet target saturation (≥ 0.7)', 
                         ha='center', va='center', fontsize=12,
                         color='#388E3C', fontweight='bold', transform=ax_table.transAxes)
    
    def _create_complexity_scatter(self, ax, data):
        """Create UMI vs Genes scatter showing library complexity - full width."""
        x = data['Median UMI per Cell'].values
        y = data['Median GeneFull_Ex50pAS per Cell'].values
        lib_nums = data['Library_Number'].values
        
        # Calculate complexity score (genes per 1000 UMIs)
        complexity = (y / x) * 1000
        
        # Use RdYlGn colormap
        scatter = ax.scatter(x, y, c=complexity, s=220, cmap='RdYlGn',
                           edgecolors='#333333', linewidth=2, alpha=0.9,
                           vmin=np.percentile(complexity, 5),
                           vmax=np.percentile(complexity, 95))
        
        # Create text labels with outline for readability
        texts = []
        for i, lib in enumerate(lib_nums):
            txt = ax.text(x[i], y[i], str(lib), fontsize=8, ha='center', va='center',
                         fontweight='bold', color='#1a1a1a')
            txt.set_path_effects([
                path_effects.Stroke(linewidth=2, foreground='white'),
                path_effects.Normal()
            ])
            texts.append(txt)
        
        # Use adjustText to move overlapping labels
        if HAS_ADJUST_TEXT:
            adjust_text(texts, x=x, y=y, ax=ax,
                       arrowprops=dict(arrowstyle='-', color='#666666', lw=0.5),
                       expand_points=(1.5, 1.5),
                       force_points=(0.5, 0.5),
                       force_text=(0.3, 0.3))
        
        # Add ideal complexity lines - IMPROVED: thicker, more visible
        umi_range = np.array([x.min() * 0.9, x.max() * 1.1])
        for ratio, label, color in [(0.3, 'Low complexity', '#E65100'),
                                     (0.4, 'Good complexity', '#2E7D32'),
                                     (0.5, 'High complexity', '#1565C0')]:
            ax.plot(umi_range, umi_range * ratio, '--', 
                   color=color, alpha=0.8, linewidth=2.5, label=label)
        
        ax.set_xlabel('Median UMI per Cell', fontsize=12, fontweight='500')
        ax.set_ylabel('Median Genes per Cell', fontsize=12, fontweight='500')
        ax.set_title('Library Complexity Analysis', fontsize=14, fontweight='600')
        
        ax.legend(loc='lower right', frameon=True, fontsize=10, framealpha=0.95)
        ax.grid(True, alpha=0.2, linestyle=':')
        
        # Larger colorbar
        cbar = plt.colorbar(scatter, ax=ax, fraction=0.025, pad=0.015)
        cbar.set_label('Complexity (Genes/1000 UMIs)', rotation=270, labelpad=18, fontsize=10)
    
    def _create_performance_heatmap_absolute(self, ax, data):
        """
        Create performance heatmap with GLOBAL normalization for coloring.
        Uses global min/max across all values for consistent color scaling.
        """
        # Select key metrics
        metrics = ['Reads With Valid Barcodes', 'Sequencing Saturation',
                  'Q30 Bases in CB+UMI', 'Q30 Bases in RNA read',
                  'Fraction of Unique Reads in Cells']
        available = [m for m in metrics if m in data.columns]
        
        # Reset index and sort by Library_Number
        data_sorted = data.sort_values('Library_Number').reset_index(drop=True)
        lib_nums = data_sorted['Library_Number'].values
        
        # Create matrix with metrics as rows, libraries as columns
        matrix = data_sorted[available].T
        
        # Global normalization: find min/max across ALL values in the matrix
        global_min = matrix.values.min()
        global_max = matrix.values.max()
        
        # Normalize all values to 0-1 range using global min/max
        matrix_norm = (matrix.values - global_min) / (global_max - global_min) if global_max != global_min else np.full_like(matrix.values, 0.5)
        
        im = ax.imshow(matrix_norm, cmap='RdYlGn', aspect='auto', vmin=0, vmax=1)
        
        ax.grid(False)
        
        # Add actual values as text with proper contrast (using globally normalized values)
        for i in range(len(available)):
            for j in range(len(lib_nums)):
                val = matrix.iloc[i, j]
                norm_val = matrix_norm[i, j]
                color = get_text_color_for_background(norm_val, vmin=0, vmax=1, cmap_name='RdYlGn')
                ax.text(j, i, f'{val:.2f}', ha='center', va='center',
                       color=color, fontsize=8, fontweight='500')
        
        # Set ticks
        ax.set_xticks(np.arange(len(lib_nums)))
        ax.set_yticks(np.arange(len(available)))
        ax.set_xticklabels([str(n) for n in lib_nums], rotation=0, fontsize=9)
        ax.set_xlabel('Library Number', fontsize=11, fontweight='500')
        
        # Shortened y-axis labels
        y_labels = []
        for m in available:
            if 'Q30 Bases' in m:
                y_labels.append(m.replace('Q30 Bases in ', 'Q30 '))
            else:
                y_labels.append(m.replace('Fraction of ', ''))
        ax.set_yticklabels(y_labels, fontsize=9)
        
        # Title explains global coloring
        ax.set_title('Performance Matrix (Global normalization: red=low, green=high)', 
                    fontsize=12, fontweight='600', pad=10)
        
        # Add colorbar
        cbar = plt.colorbar(im, ax=ax, fraction=0.02, pad=0.02)
        cbar.set_label('Quality Score', rotation=270, labelpad=15, fontsize=9)
        cbar.set_ticks([0, 0.5, 1])
        cbar.set_ticklabels(['Poor', 'OK', 'Good'])
    
    def _save_and_show(self, name):
        """Save figure and optionally show it."""
        # IMPROVED: Higher DPI for publication quality
        plt.savefig(self.output_dir / f'{name}.png', dpi=200, bbox_inches='tight',
                   facecolor='white', edgecolor='none')
        print(f"  Saved: {name}.png")
        
        if self.interactive:
            plt.show()
        else:
            plt.close()
    
    def export_tsv(self):
        """Export all processed data to a TSV file."""
        if self.data is None:
            return
        
        output_df = pd.DataFrame()
        output_df['Library'] = self.data['Library']
        output_df['Shortname'] = self.data['Library_Number']
        
        data_cols = [col for col in self.data.columns if col not in ['Library', 'Library_Number']]
        
        for col in data_cols:
            output_df[col] = self.data[col]
        
        output_df = output_df.sort_values('Shortname')
        
        output_file = self.output_dir / 'processedstats.tsv'
        output_df.to_csv(output_file, sep='\t', index=False)
        print(f"  Saved: processedstats.tsv")

    def _comparison_frames(self):
        """Return current and baseline numeric tables on their shared libraries."""
        if self.data is None or self.baseline_data is None:
            raise ValueError("Delta output requires both current and baseline statistics")

        current = self.data.set_index('Library_Number').sort_index()
        baseline = self.baseline_data.set_index('Library_Number').sort_index()
        common_libraries = current.index.intersection(baseline.index).sort_values()
        current = current.loc[common_libraries]
        baseline = baseline.loc[common_libraries]
        common_metrics = [
            column
            for column in current.columns
            if column != 'Library'
            and column in baseline.columns
            and pd.api.types.is_numeric_dtype(current[column])
            and pd.api.types.is_numeric_dtype(baseline[column])
        ]
        if not common_metrics:
            raise ValueError("Baseline and current tables have no shared numeric metrics")
        return current, baseline, common_metrics

    def export_delta_tsv(self):
        """Export an exact long-form baseline/current/delta table."""
        current, baseline, metrics = self._comparison_frames()
        rows = []
        for library_number in current.index:
            library_name = current.loc[library_number, 'Library']
            for metric in metrics:
                baseline_value = baseline.loc[library_number, metric]
                current_value = current.loc[library_number, metric]
                if pd.isna(baseline_value) or pd.isna(current_value):
                    delta = np.nan
                    percent_delta = np.nan
                else:
                    delta = current_value - baseline_value
                    percent_delta = (
                        (delta / baseline_value) * 100
                        if baseline_value != 0
                        else np.nan
                    )
                rows.append(
                    {
                        'Library': library_name,
                        'Shortname': int(library_number),
                        'Metric': metric,
                        'Baseline': baseline_value,
                        'Current': current_value,
                        'Delta': delta,
                        'Percent Delta': percent_delta,
                    }
                )

        output_file = self.output_dir / 'mapping_stats_deltas.tsv'
        pd.DataFrame(rows).to_csv(output_file, sep='\t', index=False)
        print("  Saved: mapping_stats_deltas.tsv")

    def plot_mapping_delta_dashboard(self):
        """Plot longitudinal changes in depth, cell recovery, and complexity."""
        current, baseline, _metrics = self._comparison_frames()
        required = [
            'Median Reads per Cell',
            'Estimated Number of Cells',
            'Sequencing Saturation',
            'Median UMI per Cell',
            'Median GeneFull_Ex50pAS per Cell',
        ]
        missing = [
            metric
            for metric in required
            if metric not in current.columns or metric not in baseline.columns
        ]
        if missing:
            raise ValueError(
                "Delta dashboard is missing required metric(s): " + ', '.join(missing)
            )

        libraries = current.index.to_numpy(dtype=int)
        x = np.arange(len(libraries))

        baseline_reads = baseline['Median Reads per Cell'].to_numpy(dtype=float) / 1000
        current_reads = current['Median Reads per Cell'].to_numpy(dtype=float) / 1000
        added_reads = current_reads - baseline_reads
        baseline_cells = baseline['Estimated Number of Cells'].to_numpy(dtype=float)
        current_cells = current['Estimated Number of Cells'].to_numpy(dtype=float)
        cell_delta = current_cells - baseline_cells
        baseline_saturation = baseline['Sequencing Saturation'].to_numpy(dtype=float)
        current_saturation = current['Sequencing Saturation'].to_numpy(dtype=float)

        baseline_umi = baseline['Median UMI per Cell'].to_numpy(dtype=float)
        current_umi = current['Median UMI per Cell'].to_numpy(dtype=float)
        baseline_genes = baseline['Median GeneFull_Ex50pAS per Cell'].to_numpy(dtype=float)
        current_genes = current['Median GeneFull_Ex50pAS per Cell'].to_numpy(dtype=float)
        umi_percent = np.divide(
            current_umi - baseline_umi,
            baseline_umi,
            out=np.full_like(current_umi, np.nan),
            where=baseline_umi != 0,
        ) * 100
        gene_percent = np.divide(
            current_genes - baseline_genes,
            baseline_genes,
            out=np.full_like(current_genes, np.nan),
            where=baseline_genes != 0,
        ) * 100

        reads_percent = np.divide(
            added_reads,
            baseline_reads,
            out=np.full_like(added_reads, np.nan),
            where=baseline_reads != 0,
        ) * 100
        cells_percent = np.divide(
            cell_delta,
            baseline_cells,
            out=np.full_like(cell_delta, np.nan),
            where=baseline_cells != 0,
        ) * 100
        saturation_delta = current_saturation - baseline_saturation

        fig = plt.figure(figsize=(18, 16), facecolor='white')
        gs = GridSpec(4, 1, figure=fig, hspace=0.43)
        fig.suptitle(
            'RNA Mapping Longitudinal Delta Dashboard',
            fontsize=19,
            fontweight='700',
            y=0.985,
            color='#1a1a1a',
        )
        summary = (
            f"{len(libraries)} matched libraries · median change: "
            f"reads/cell {np.nanmedian(reads_percent):+.1f}% · "
            f"estimated cells {np.nanmedian(cells_percent):+.1f}% · "
            f"saturation {np.nanmedian(saturation_delta):+.3f} · "
            f"UMIs/cell {np.nanmedian(umi_percent):+.1f}% · "
            f"genes/cell {np.nanmedian(gene_percent):+.1f}%"
        )
        fig.text(0.5, 0.962, summary, ha='center', va='center', fontsize=11, color='#424242')

        # 1. Stacked depth: the baseline plus the increment observed in the current run.
        ax = fig.add_subplot(gs[0])
        ax.bar(
            x,
            baseline_reads,
            width=0.78,
            color='#B0BEC5',
            edgecolor='#455A64',
            linewidth=0.8,
            label=self.baseline_label,
        )
        positive = np.clip(added_reads, 0, None)
        negative = np.clip(-added_reads, 0, None)
        ax.bar(
            x,
            positive,
            bottom=baseline_reads,
            width=0.78,
            color='#00897B',
            edgecolor='#00695C',
            linewidth=0.8,
            label=f'Added since {self.baseline_label}',
        )
        if np.any(negative > 0):
            ax.bar(
                x,
                negative,
                bottom=current_reads,
                width=0.78,
                color='#EF5350',
                edgecolor='#B71C1C',
                linewidth=0.8,
                label='Decrease',
            )
        ax.scatter(x, current_reads, marker='_', s=110, linewidth=1.8, color='#102027', zorder=4)
        ax.set_title('Median Reads per Cell: Baseline plus Added Sequencing', fontweight='600')
        ax.set_ylabel('Median reads per cell (thousands)')
        ax.legend(loc='upper left', ncol=3, frameon=True)
        ax.grid(True, axis='y', alpha=0.25, linestyle=':')

        # 2. Cell recovery is deliberately shown as a delta, not an additive read component.
        ax = fig.add_subplot(gs[1])
        cell_colors = np.where(cell_delta >= 0, '#2E7D32', '#C62828')
        ax.bar(x, cell_delta, width=0.78, color=cell_colors, edgecolor='white', linewidth=0.5)
        ax.axhline(0, color='#424242', linewidth=1)
        ax.set_title('Change in Estimated Number of Cells', fontweight='600')
        ax.set_ylabel(f'{self.current_label} − {self.baseline_label} (cells)')
        ax.grid(True, axis='y', alpha=0.25, linestyle=':')

        # 3. Saturation is a bounded rate and therefore is compared, never stacked.
        ax = fig.add_subplot(gs[2])
        ax.plot(
            x,
            baseline_saturation,
            color='#78909C',
            marker='o',
            markersize=4,
            linewidth=1.5,
            label=self.baseline_label,
        )
        ax.plot(
            x,
            current_saturation,
            color='#1565C0',
            marker='o',
            markersize=4,
            linewidth=1.8,
            label=self.current_label,
        )
        ax.axhline(0.7, color='#616161', linestyle='--', linewidth=1.5, label='Target: 0.7')
        ax.set_ylim(0, 1.02)
        ax.set_title('Sequencing Saturation Before and After Additional Sequencing', fontweight='600')
        ax.set_ylabel('Sequencing saturation')
        ax.legend(loc='lower right', ncol=3, frameon=True)
        ax.grid(True, axis='y', alpha=0.25, linestyle=':')

        # 4. Molecule and gene gains show whether added reads improved usable complexity.
        ax = fig.add_subplot(gs[3])
        width = 0.38
        ax.bar(x - width / 2, umi_percent, width, color='#5E35B1', label='Median UMIs/cell')
        ax.bar(x + width / 2, gene_percent, width, color='#F9A825', label='Median genes/cell')
        ax.axhline(0, color='#424242', linewidth=1)
        ax.set_title('Per-Cell Complexity Gain', fontweight='600')
        ax.set_ylabel('Change from baseline (%)')
        ax.legend(loc='upper right', ncol=2, frameon=True)
        ax.grid(True, axis='y', alpha=0.25, linestyle=':')

        for ax in fig.axes:
            ax.set_xlim(-0.7, len(libraries) - 0.3)
            ax.set_xticks(x)
            ax.set_xticklabels([str(value) for value in libraries], fontsize=8)
            ax.set_xlabel('Library number')

        plt.tight_layout(rect=[0, 0.01, 1, 0.95])
        self._save_and_show('mapping_delta_dashboard')
    
    def generate_all_plots(self):
        """Generate all plots."""
        print("\nGenerating visualizations...")
        print("=" * 60)
        
        plot_methods = [
            ('Quality Control Dashboard', self.plot_quality_control_dashboard),
            ('Cell Quality Matrix', self.plot_cell_quality_matrix),
            ('Outlier Detection Report', self.plot_outlier_detection_report),
        ]
        if self.baseline_data is not None:
            plot_methods.append(
                ('Longitudinal Delta Dashboard', self.plot_mapping_delta_dashboard)
            )
        if self.cell_reads_data is not None:
            plot_methods.append(
                ('Reads per Called Cell Distribution', self.plot_reads_per_cell_distribution)
            )

        failures = []
        for name, method in plot_methods:
            try:
                print(f"\nCreating: {name}")
                method()
            except Exception as e:
                print(f"  Error with {name}: {e}")
                import traceback
                traceback.print_exc()
                failures.append(name)

        if failures:
            raise RuntimeError(
                "Plot generation failed for: " + ', '.join(failures)
            )
        
        print("\n" + "=" * 60)
        print(f"All visualizations saved to: {self.output_dir.absolute()}")
        print("=" * 60)


def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description='Create visualizations from RNA mapping statistics',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s --base-path /path/to/libraries
  %(prog)s --base-path /path/to/libraries --interactive
  %(prog)s --base-path /path/to/libraries --output-dir results/plots
  %(prog)s --current-stats new/processedstats.tsv --baseline-stats old/processedstats.tsv --plot-type deltas
  %(prog)s --base-path /path/to/mapping_output --cell-read-distribution --cell-reads-work-root /path/to/nextflow/work
        """
    )
    
    parser.add_argument('--base-path', '-b', 
                       default='/mnt/beegfs/tetmultiome_rna_mapped/mapping_output/',
                       help='Base path containing library folders')
    parser.add_argument('--current-stats',
                       default=None,
                       help='Use this processedstats.tsv as the current dataset instead of reading Summary.csv files')
    parser.add_argument('--baseline-stats',
                       default=None,
                       help='Prior processedstats.tsv used to produce longitudinal delta output')
    parser.add_argument('--baseline-label',
                       default='Baseline',
                       help='Display label for the prior statistics')
    parser.add_argument('--current-label',
                       default='Current',
                       help='Display label for the current statistics')
    parser.add_argument('--cell-read-distribution',
                       action='store_true',
                       help='Create a validated distribution of unique gene-assigned reads per called cell')
    parser.add_argument('--cell-reads-work-root',
                       default=None,
                       help='Optional legacy Nextflow work directory; usable only if STAR generated CellReads.stats')
    parser.add_argument('--baseline-base-path',
                       default=None,
                       help='Optional prior mapping_output root for an old-to-new cell-read distribution')
    parser.add_argument('--baseline-cell-reads-work-root',
                       default=None,
                       help='Optional prior work directory; usable only if STAR generated CellReads.stats')
    parser.add_argument('--output-dir', '-o', 
                       default='mapping_plots',
                       help='Output directory for plots')
    parser.add_argument('--interactive', '-i', action='store_true',
                       help='Show plots interactively (matplotlib window)')
    parser.add_argument('--plot-type', '-p',
                       choices=['all', 'qc', 'matrix', 'outliers', 'deltas', 'cell-reads'],
                       default='all',
                       help='Type of plot to generate')
    
    args = parser.parse_args()
    
    if args.plot_type == 'deltas' and not args.baseline_stats:
        parser.error('--plot-type deltas requires --baseline-stats')
    if args.plot_type == 'cell-reads':
        args.cell_read_distribution = True
    if args.baseline_cell_reads_work_root and not args.baseline_base_path:
        parser.error('--baseline-cell-reads-work-root requires --baseline-base-path')

    print(f"RNA Mapping Visualization Tool {RELEASE}")
    if args.current_stats:
        print(f"Current statistics: {args.current_stats}")
    else:
        print(f"Base path: {args.base_path}")
    
    visualizer = MappingStatsVisualizer(
        args.base_path,
        args.output_dir,
        args.interactive,
        baseline_stats=args.baseline_stats,
        baseline_label=args.baseline_label,
        current_label=args.current_label,
        cell_reads_work_root=args.cell_reads_work_root,
        baseline_base_path=args.baseline_base_path,
        baseline_cell_reads_work_root=args.baseline_cell_reads_work_root,
    )
    
    data = (
        visualizer.load_processed_stats(args.current_stats)
        if args.current_stats
        else visualizer.load_data()
    )
    
    if data is None:
        print("Failed to load data. Exiting.")
        return 1

    if args.baseline_stats:
        visualizer.load_baseline()
    if args.cell_read_distribution:
        visualizer.load_cell_reads()
    
    print("\nExporting data to TSV...")
    visualizer.export_tsv()
    if visualizer.baseline_data is not None:
        visualizer.export_delta_tsv()
    if visualizer.cell_reads_data is not None:
        visualizer.export_cell_reads()
    
    if args.plot_type == 'all':
        visualizer.generate_all_plots()
    else:
        plot_map = {
            'qc': visualizer.plot_quality_control_dashboard,
            'matrix': visualizer.plot_cell_quality_matrix,
            'outliers': visualizer.plot_outlier_detection_report,
            'deltas': visualizer.plot_mapping_delta_dashboard,
            'cell-reads': visualizer.plot_reads_per_cell_distribution,
        }
        
        if args.plot_type in plot_map:
            print(f"\nGenerating {args.plot_type} visualization...")
            plot_map[args.plot_type]()
            print(f"Saved to: {visualizer.output_dir.absolute()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
