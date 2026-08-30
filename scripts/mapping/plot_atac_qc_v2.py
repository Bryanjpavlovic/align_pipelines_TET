#!/usr/bin/env python3
"""
ATAC-seq QC Visualization V2 - Boxplots & Aggregated Fragment Distribution

Produces:
  1. atac_qc_boxplots.png - QC metrics as clean boxplots
  2. atac_fragment_distribution.png - Median fragment size with SD ribbon

Usage:
  python plot_atac_qc_v2.py --stats-dir /mnt/beegfs/tetmultiome_atac/mapping_output/atac_qc_stats --output-dir /mnt/beegfs/tetmultiome_atac/mapping_output/PlotsV2
"""

import argparse
import json
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from pathlib import Path
import warnings

warnings.filterwarnings('ignore')

# Styling
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 0.8
plt.rcParams['axes.edgecolor'] = '#333333'

# Color palette - semantic by data type
COLOR_COUNT = '#f0c93a'      # Yellow - counts/cells
COLOR_PERCENT = '#6ab04c'    # Green - percentages  
COLOR_DEPTH = '#22a6b3'      # Cyan - sequencing depth per cell
COLOR_RATE = '#eb8c34'       # Orange - rates (mapping, duplication)

COLORS = {
    'num_cells': COLOR_COUNT,
    'frac_reads_in_cells': COLOR_PERCENT,
    'median_reads_per_cell': COLOR_DEPTH,
    'median_frags_per_cell': COLOR_DEPTH,
    'mapping_rate': COLOR_RATE,
    'duplicate_rate': COLOR_RATE,
    'mito_fraction': COLOR_PERCENT,
    'frac_frags_in_cells': COLOR_PERCENT,
}

BLUE = '#4a90d9'


def extract_library_number(lib_name):
    """Extract numeric ID from library name."""
    match = re.search(r'(\d+)$', str(lib_name))
    return int(match.group(1)) if match else float('inf')


def load_all_stats(stats_dir):
    """Load all JSON stats files."""
    stats_dir = Path(stats_dir)
    data_list = []
    fragment_hists = {}
    
    json_files = list(stats_dir.glob("*_stats.json"))
    print(f"Found {len(json_files)} stats files")
    
    for json_file in sorted(json_files):
        try:
            with open(json_file) as f:
                stats = json.load(f)
            
            lib_name = stats.get('Library', stats.get('atac_library', json_file.stem))
            lib_num = stats.get('Library_Number', stats.get('library_number', extract_library_number(lib_name)))
            
            # Extract fragment histogram
            hist_key = 'valid_bc_no_mito_frag_size_hist'
            if hist_key in stats and stats[hist_key]:
                fragment_hists[lib_num] = np.array(stats[hist_key])
            
            # Build row with consistent naming
            row = {
                'Library': lib_name,
                'Library_Number': lib_num,
                'num_cells': stats.get('valid_bc_no_mito_num_cells', stats.get('num_cells', 0)),
                'frac_reads_in_cells': stats.get('valid_bc_no_mito_frac_reads', stats.get('frac_reads_in_cells', 0)),
                'median_reads_per_cell': stats.get('valid_bc_no_mito_median_reads_per_cell', stats.get('median_reads_per_cell', 0)),
                'median_frags_per_cell': stats.get('valid_bc_no_mito_median_frags_per_cell', stats.get('median_frags_per_cell', 0)),
                'mapping_rate': stats.get('mapping_rate', 0),
                'duplicate_rate': stats.get('duplicate_rate', 0),
                'mito_fraction': stats.get('mito_fraction', 0),
                'frac_frags_in_cells': stats.get('valid_bc_no_mito_frac_frags', stats.get('frac_frags_in_cells', 0)),
                'nfr_fraction': stats.get('valid_bc_no_mito_nfr_fraction', stats.get('nfr_fraction', 0)),
                'mono_nuc_fraction': stats.get('valid_bc_no_mito_mono_nuc_fraction', stats.get('mono_nuc_fraction', 0)),
                'di_nuc_fraction': stats.get('valid_bc_no_mito_di_nuc_fraction', stats.get('di_nuc_fraction', 0)),
                'total_reads': stats.get('total_reads', 0),
                'total_frags': stats.get('total_fragments', stats.get('total_frags', 0)),
                'reads_in_cells': stats.get('valid_bc_no_mito_total_reads', stats.get('reads_in_cells', 0)),
                'frags_in_cells': stats.get('valid_bc_no_mito_total_frags', stats.get('frags_in_cells', 0)),
            }
            data_list.append(row)
            
        except Exception as e:
            print(f"  Error loading {json_file.name}: {e}")
    
    if not data_list:
        return None, None
    
    df = pd.DataFrame(data_list)
    df = df.sort_values('Library_Number').reset_index(drop=True)
    
    return df, fragment_hists


def plot_qc_boxplots(df, output_dir):
    """Create clean boxplot dashboard for QC metrics."""
    
    metrics = [
        ('num_cells', 'Number of Cells', None, 1),
        ('frac_reads_in_cells', 'Fraction Reads in Cells (%)', '%', 1),
        ('median_reads_per_cell', 'Median Reads per Cell', None, 1000),
        ('median_frags_per_cell', 'Median Fragments per Cell', None, 1000),
        ('mapping_rate', 'Mapping Rate (%)', '%', 1),
        ('duplicate_rate', 'Duplicate Rate (%)', '%', 1),
        ('mito_fraction', 'Mitochondrial Fraction (%)', '%', 1),
        ('frac_frags_in_cells', 'Fraction Fragments in Cells (%)', '%', 1),
    ]
    
    fig, axes = plt.subplots(2, 4, figsize=(11, 6))
    axes = axes.flatten()
    
    for idx, (col, title, suffix, divisor) in enumerate(metrics):
        ax = axes[idx]
        color = COLORS.get(col, BLUE)
        
        # Gray header bar for title (mimicking facet style)
        ax.set_title(title, fontsize=9, fontweight='600', 
                    bbox=dict(boxstyle='square,pad=0.3', facecolor='#e0e0e0', edgecolor='none'),
                    loc='center', y=1.0)
        
        if col not in df.columns:
            ax.text(0.5, 0.5, f'No data', ha='center', va='center', transform=ax.transAxes)
            ax.set_xticks([])
            continue
        
        values = df[col].dropna() / divisor
        
        if len(values) == 0:
            ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
            ax.set_xticks([])
            continue
        
        # Clean boxplot - narrow
        bp = ax.boxplot(values, vert=True, widths=0.4, patch_artist=True,
                       boxprops=dict(facecolor=color, edgecolor='#333333', linewidth=1, alpha=0.9),
                       whiskerprops=dict(color='#333333', linewidth=1),
                       capprops=dict(color='#333333', linewidth=1),
                       medianprops=dict(color='#1a1a1a', linewidth=1.5),
                       flierprops=dict(marker='o', markerfacecolor='#555555', markersize=4, 
                                      markeredgecolor='none', alpha=0.6))
        
        # Clean up axes
        ax.set_xticks([])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        
        # Y-axis label
        if divisor == 1000:
            ax.set_ylabel('Thousands', fontsize=8)
        elif suffix == '%':
            ax.set_ylabel('Percent', fontsize=8)
        else:
            ax.set_ylabel('Count', fontsize=8)
        
        # Subtle grid
        ax.yaxis.grid(True, linestyle='-', alpha=0.15, color='#666666')
        ax.set_axisbelow(True)
        
        # Tighten x limits to make box appear narrower
        ax.set_xlim(0.4, 1.6)
    
    plt.tight_layout(h_pad=1.5, w_pad=0.8)
    fig.suptitle('ATAC-seq QC Dashboard: RNA-Filtered Cells (excluding chrM)', 
                 fontsize=11, fontweight='bold', y=1.02)
    
    output_path = output_dir / 'atac_qc_boxplots.png'
    plt.savefig(output_path, dpi=200, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"  Saved: {output_path.name}")


def plot_fragment_distribution_aggregated(fragment_hists, output_dir):
    """Create aggregated fragment size distribution - smooth curves, dark blue styling."""
    
    if not fragment_hists:
        print("  No fragment histogram data available")
        return
    
    from scipy.ndimage import gaussian_filter1d
    
    # Bin centers (10bp bins from 0-1000)
    bin_centers = np.arange(5, 1000, 10)
    
    # Normalize each histogram and collect
    all_densities = []
    for lib_num, hist in sorted(fragment_hists.items()):
        hist = np.array(hist)
        if len(hist) == 0 or hist.sum() == 0:
            continue
        # Normalize to density
        density = hist / hist.sum()
        # Pad or truncate to match bin_centers length
        if len(density) < len(bin_centers):
            density = np.pad(density, (0, len(bin_centers) - len(density)), mode='constant')
        else:
            density = density[:len(bin_centers)]
        all_densities.append(density)
    
    if len(all_densities) < 1:
        print(f"  No libraries with histogram data")
        return
    
    all_densities = np.array(all_densities)
    n_libs = len(all_densities)
    
    # Calculate statistics
    mean_density = np.mean(all_densities, axis=0)
    std_density = np.std(all_densities, axis=0)
    
    # Smooth the curves
    sigma = 2  # Smoothing factor
    mean_smooth = gaussian_filter1d(mean_density, sigma=sigma)
    std_smooth = gaussian_filter1d(std_density, sigma=sigma)
    
    # Dark blue color scheme
    DARK_BLUE = '#1e3a5f'
    LIGHT_BLUE = '#6fa8dc'
    
    # Plot
    fig, ax = plt.subplots(figsize=(12, 4.5))
    
    # SD ribbon (±1 SD) - smooth
    lower = mean_smooth - std_smooth
    upper = mean_smooth + std_smooth
    lower = np.maximum(lower, 0)  # No negative density
    
    ax.fill_between(bin_centers, lower, upper, 
                   alpha=0.35, color=LIGHT_BLUE, label='±1 SD', linewidth=0)
    
    # Mean line - smooth, dark blue
    ax.plot(bin_centers, mean_smooth, color=DARK_BLUE, linewidth=2.5, label='Mean')
    
    # Nucleosome markers - gray dashed/dotted lines
    ax.axvline(x=147, color='#555555', linestyle='--', linewidth=1.2, alpha=0.8)
    ax.axvline(x=294, color='#555555', linestyle=':', linewidth=1.2, alpha=0.8)
    ax.axvline(x=441, color='#555555', linestyle=':', linewidth=1.2, alpha=0.6)
    
    # Get y limits for annotation placement
    ymax = upper.max() * 1.12
    ax.set_ylim(0, ymax)
    
    # Annotations at top
    ax.text(80, ymax * 0.94, 'NFR', fontsize=12, fontweight='bold', ha='center', va='top')
    ax.text(147, ymax * 0.94, 'Mono-nuc', fontsize=10, ha='center', va='top', color='#444444')
    ax.text(294, ymax * 0.94, 'Di-nuc', fontsize=10, ha='center', va='top', color='#444444')
    ax.text(441, ymax * 0.94, 'Tri-nuc', fontsize=10, ha='center', va='top', color='#555555')
    
    ax.set_xlabel('Fragment Size (bp)', fontsize=12)
    ax.set_ylabel('Density', fontsize=12)
    ax.set_title(f'Fragment Size Distribution (n={n_libs} libraries, RNA-filtered cells, no chrM)', 
                fontsize=13, fontweight='bold')
    ax.set_xlim(0, 800)
    
    # Clean axes
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Legend
    ax.legend(loc='upper right', framealpha=0.95, edgecolor='none', fontsize=10)
    
    plt.tight_layout()
    
    output_path = output_dir / 'atac_fragment_distribution.png'
    plt.savefig(output_path, dpi=200, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"  Saved: {output_path.name}")


def main():
    parser = argparse.ArgumentParser(
        description='ATAC-seq QC Visualization V2 - Boxplots & Aggregated Fragment Distribution'
    )
    parser.add_argument('--stats-dir', '-s', required=True,
                       help='Directory containing *_stats.json files')
    parser.add_argument('--output-dir', '-o', 
                       help='Output directory (default: stats_dir/../atac_plots_v2)')
    args = parser.parse_args()
    
    stats_dir = Path(args.stats_dir)
    output_dir = Path(args.output_dir) if args.output_dir else stats_dir.parent / 'atac_plots_v2'
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"Stats directory: {stats_dir}")
    print(f"Output directory: {output_dir}")
    
    # Load data
    print("\nLoading data...")
    df, fragment_hists = load_all_stats(stats_dir)
    
    if df is None or len(df) == 0:
        print("ERROR: No data loaded")
        return 1
    
    print(f"Loaded {len(df)} libraries")
    print(f"Fragment histograms available for {len(fragment_hists)} libraries")
    
    # Generate plots
    print("\nGenerating plots...")
    plot_qc_boxplots(df, output_dir)
    plot_fragment_distribution_aggregated(fragment_hists, output_dir)
    
    print(f"\nDone! Outputs in: {output_dir}")
    return 0


if __name__ == "__main__":
    exit(main())
