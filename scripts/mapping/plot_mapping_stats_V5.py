#!/usr/bin/env python3
"""
RNA mapping statistics visualization tool.
Generates focused plots from Summary.csv files across library folders.
V5 - Enhanced visual clarity, better contrast, absolute value heatmaps, cleaner aesthetics

Run with: python plot_mapping_stats_V5.py --base-path ./
"""

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
    def __init__(self, base_path, output_dir="plots", interactive=False):
        """Initialize visualizer with data path and output settings."""
        self.base_path = Path(base_path)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        self.interactive = interactive
        self.data = None
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
            'Mean Reads per Cell',
            'Mean UMI per Cell', 
            'Mean GeneFull_Ex50pAS per Cell'
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
        cell_metrics = ['Mean Reads per Cell', 'Mean UMI per Cell', 
                       'Mean GeneFull_Ex50pAS per Cell', 'Estimated Number of Cells']
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
        if 'Mean Reads per Cell' in available and 'Mean UMI per Cell' in available:
            ax = axes[1]
            x = self.data['Mean Reads per Cell']
            y = self.data['Mean UMI per Cell']
            
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
        if 'Estimated Number of Cells' in available and 'Mean GeneFull_Ex50pAS per Cell' in available:
            ax = axes[2]
            x = self.data['Estimated Number of Cells']
            y = self.data['Mean GeneFull_Ex50pAS per Cell']
            
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
        metrics = self.quality_metrics + self.efficiency_metrics + ['Mean UMI per Cell']
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
        reads_per_cell = data['Mean Reads per Cell'].values if 'Mean Reads per Cell' in data.columns else None
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
        x = data['Mean UMI per Cell'].values
        y = data['Mean GeneFull_Ex50pAS per Cell'].values
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
    
    def generate_all_plots(self):
        """Generate all plots."""
        print("\nGenerating visualizations...")
        print("=" * 60)
        
        plot_methods = [
            ('Quality Control Dashboard', self.plot_quality_control_dashboard),
            ('Cell Quality Matrix', self.plot_cell_quality_matrix),
            ('Outlier Detection Report', self.plot_outlier_detection_report),
        ]
        
        for name, method in plot_methods:
            try:
                print(f"\nCreating: {name}")
                method()
            except Exception as e:
                print(f"  Error with {name}: {e}")
                import traceback
                traceback.print_exc()
        
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
        """
    )
    
    parser.add_argument('--base-path', '-b', 
                       default='/mnt/beegfs/tetmultiome_rna_mapped/mapping_output/',
                       help='Base path containing library folders')
    parser.add_argument('--output-dir', '-o', 
                       default='mapping_plots',
                       help='Output directory for plots')
    parser.add_argument('--interactive', '-i', action='store_true',
                       help='Show plots interactively (matplotlib window)')
    parser.add_argument('--plot-type', '-p',
                       choices=['all', 'qc', 'matrix', 'outliers'],
                       default='all',
                       help='Type of plot to generate')
    
    args = parser.parse_args()
    
    print(f"RNA Mapping Visualization Tool v5.0")
    print(f"Base path: {args.base_path}")
    
    visualizer = MappingStatsVisualizer(args.base_path, args.output_dir, args.interactive)
    
    data = visualizer.load_data()
    
    if data is None:
        print("Failed to load data. Exiting.")
        return
    
    print("\nExporting data to TSV...")
    visualizer.export_tsv()
    
    if args.plot_type == 'all':
        visualizer.generate_all_plots()
    else:
        plot_map = {
            'qc': visualizer.plot_quality_control_dashboard,
            'matrix': visualizer.plot_cell_quality_matrix,
            'outliers': visualizer.plot_outlier_detection_report
        }
        
        if args.plot_type in plot_map:
            print(f"\nGenerating {args.plot_type} visualization...")
            plot_map[args.plot_type]()
            print(f"Saved to: {visualizer.output_dir.absolute()}")


if __name__ == "__main__":
    main()
