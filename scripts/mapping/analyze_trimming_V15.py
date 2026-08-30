#!/usr/bin/env python3
"""
Comprehensive visualization suite for cutadapt JSON reports.
Combines global and library-level analysis with aggregation across sequencing runs.
Features intelligent library name truncation for cleaner visualizations.

Version 13 Updates:
- Handles both long R2 format (290bp R2) and paired-end format (151bp × 151bp)
- For long R2 format: Aggregates R2 files by library name
- For paired-end format: Keeps R1 and R2 as SEPARATE libraries (no aggregation)
- Automatically detects format based on filename patterns

Usage:
    python analyze_trimming_V13.py                    # Default: looks in logs/
    python analyze_trimming_V13.py --path /my/data/   # Custom path
    python analyze_trimming_V13.py -p .               # Current directory
    python analyze_trimming_V13.py --full-names       # Use full library names

Output:
    - cutadapt_figures/ directory with all visualizations
      - global/ - Cross-library analysis
      - by_library/ - Per-library aggregated analysis
    - adapter_detailed.csv, adapter_by_sample.csv, adapter_summary.csv

File Naming Examples:
    Long R2 format:
      Sample_39_S19_L001_R2_001_cutadapt.json → Library: Sample_39
      Sample_39_S19_L002_R2_001_cutadapt.json → Library: Sample_39 (aggregated)
    
    Paired-end format:
      Sample_5p_S01_L001_R1_001_cutadapt.json → Library: Sample_5p_R1
      Sample_5p_S01_L001_R2_001_cutadapt.json → Library: Sample_5p_R2 (separate)
"""

import json
import glob
import re
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from collections import defaultdict, Counter
from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist, squareform

# Set style with fallback
try:
    plt.style.use('seaborn-v0_8-darkgrid')
except:
    plt.style.use('default')
    sns.set_style("darkgrid")

sns.set_palette("husl")

# Define consistent adapter colors to avoid reuse across plots
ADAPTER_COLORS = {
    # TSO-related (blues/cyans)
    'TSO': '#1f77b4',
    '^TSO': '#0d5a8f',  # Darker blue for anchored TSO
    'TSO_RC': '#aec7e8',
    'TSO_5prime_marker': '#17becf',
    'TSO_5prime_marker_RC': '#9edae5',
    
    # TruSeq adapters (oranges/reds)
    'TruSeq_Read2_Universal': '#ff7f0e',
    'TruSeq_Read1_RC': '#ffbb78',
    'TruSeq_Read2_RC_Primer': '#d62728',
    'TruSeq Read2 Primer': '#ff9896',
    'TruSeq_Read2_Primer': '#ff9896',
    
    # P7/P5 adapters (purples/pinks)
    'P7_adapter': '#9467bd',
    'P7': '#9467bd',
    'P7_RC': '#c5b0d5',
    'P5_adapter_RC': '#e377c2',
    'P5 RC': '#f7b6d2',
    
    # cDNA/Library primers (greens/browns)
    'cDNA_Reverse_V3_V4': '#2ca02c',
    'cDNA_Reverse_5prime': '#98df8a',
    'Library_PCR_P1_Full': '#8c564b',
    'Library_PCR_P2_Core': '#c49c94',
    
    # RT primers (yellows/browns)
    'PolyT_RT_Primer_5prime': '#bcbd22',
    
    # Poly-G (gray)
    'G{20}': '#7f7f7f',
    'Poly-G': '#7f7f7f',
    '1': '#7f7f7f',  # Poly-G is reported as "1" in cutadapt output
    
    # Other/fallback
    'other': '#c7c7c7',
}

def get_adapter_color(adapter_name, idx=0):
    """
    Get consistent color for an adapter.
    Falls back to a distinct color from a large palette if not in dictionary.
    """
    if adapter_name in ADAPTER_COLORS:
        return ADAPTER_COLORS[adapter_name]
    
    # Fallback: use tab20 colormap for consistent distinct colors
    import matplotlib.pyplot as plt
    try:
        tab20 = plt.get_cmap('tab20')
    except AttributeError:
        # Fallback for older matplotlib versions
        from matplotlib import cm
        tab20 = cm.get_cmap('tab20')
    return tab20(idx % 20)


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Analyze cutadapt JSON reports and generate visualizations.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s --path /path/to/logs/                      # Analyze JSON files in directory
  %(prog)s -p /mnt/beegfs/trimming/BP12952/long_r2/logs/
  %(prog)s -p ./logs/                                 # Current directory's logs
  %(prog)s --path /data/trimming/ --full-names        # Don't truncate library names
  %(prog)s -p /data/logs/ --no-clustering             # Disable hierarchical clustering
        """
    )
    parser.add_argument(
        '-p', '--path',
        required=True,
        help='Path to directory containing cutadapt JSON files'
    )
    parser.add_argument(
        '--full-names',
        action='store_true',
        help='Use full library names without truncation'
    )
    parser.add_argument(
        '--no-clustering',
        action='store_true',
        help='Disable hierarchical clustering in heatmaps'
    )
    return parser.parse_args()


def extract_library_and_run(sample_name):
    """
    Extract library name and run identifier from sample filename.
    
    For long R2 format (only R2 files exist, no R1):
        - Example: Sample_39_S19_L001_R2_001 → ('Sample_39', 'S19_L001')
    
    For paired-end format (both R1 and R2 exist):
        - Keeps R1 and R2 as SEPARATE libraries
        - Example: Sample_5p_S01_L001_R1_001 → ('Sample_5p_R1', 'S01_L001')
        - Example: Sample_5p_S01_L001_R2_001 → ('Sample_5p_R2', 'S01_L001')
    """
    name = sample_name.replace('_cutadapt', '')
    
    # Detect if this is R1 or R2
    is_r1 = '_R1_' in name or name.endswith('_R1')
    is_r2 = '_R2_' in name or name.endswith('_R2')
    
    # Remove R1/R2 suffix from name
    name_clean = name.replace('_R2_001', '').replace('_R1_001', '')
    name_clean = name_clean.replace('_R2', '').replace('_R1', '')
    
    # Extract library base and run info
    match = re.match(r'(.+?)_(S\d+_L\d+)', name_clean)
    if match:
        library_base = match.group(1)
        run_info = match.group(2)
    else:
        library_base = name_clean
        run_info = 'Run1'
    
    # For paired-end format, append _R1 or _R2 to keep them separate
    # For long R2 format, R1 won't exist so just use base name
    if is_r1:
        library_name = f"{library_base}_R1"
    elif is_r2:
        library_name = f"{library_base}_R2"
    else:
        library_name = library_base
    
    return library_name, run_info


def find_common_prefix(strings):
    """Find the longest common prefix among a list of strings."""
    if not strings:
        return ""
    
    # Find the shortest string length
    min_length = min(len(s) for s in strings)
    
    # Find common prefix
    for i in range(min_length):
        char = strings[0][i]
        if not all(s[i] == char for s in strings):
            # Check if we're breaking at a word boundary
            prefix = strings[0][:i]
            # If prefix ends with underscore or dash, include it
            if prefix and prefix[-1] in ['_', '-']:
                return prefix
            # Otherwise, try to break at last underscore or dash
            last_sep = max(prefix.rfind('_'), prefix.rfind('-'))
            if last_sep > 0:
                return prefix[:last_sep + 1]
            return prefix
    
    # All strings start the same up to min_length
    return strings[0][:min_length]


def find_common_suffix(strings):
    """Find the longest common suffix among a list of strings."""
    if not strings:
        return ""
    
    # Reverse strings to find suffix
    reversed_strings = [s[::-1] for s in strings]
    common_prefix = find_common_prefix(reversed_strings)
    return common_prefix[::-1]


def create_display_names(library_names, use_full_names=False):
    """Create shortened display names for libraries."""
    if use_full_names or len(library_names) <= 1:
        return {name: name for name in library_names}
    
    # Find common prefix and suffix
    common_prefix = find_common_prefix(library_names)
    common_suffix = find_common_suffix(library_names)
    
    display_names = {}
    
    for name in library_names:
        # Remove common prefix and suffix
        display = name
        
        # Remove prefix if it's substantial (more than just a delimiter)
        if len(common_prefix) > 1:
            display = name[len(common_prefix):]
        
        # Remove suffix if it exists and is substantial
        if common_suffix and len(common_suffix) > 1:
            # Make sure we don't remove everything
            suffix_start = len(display) - len(common_suffix)
            if suffix_start > 0:
                display = display[:suffix_start]
        
        # If display name is empty or just numbers, keep some context
        if not display or display.isdigit():
            # Try to keep the last meaningful part after the prefix
            parts = name[len(common_prefix):].split('_')
            display = parts[0] if parts else name
        
        # Clean up leading/trailing underscores or dashes
        display = display.strip('_-')
        
        # If still empty, use the full name
        if not display:
            display = name
            
        display_names[name] = display
    
    # Check for duplicates and add minimal context if needed
    display_counts = Counter(display_names.values())
    for name, display in display_names.items():
        if display_counts[display] > 1:
            # Add more context from the original name
            # Try to add the next segment after the common prefix
            remaining = name[len(common_prefix):] if common_prefix else name
            parts = remaining.split('_')
            if len(parts) > 1:
                display_names[name] = '_'.join(parts[:2]).strip('_-')
            else:
                display_names[name] = remaining.strip('_-')
    
    return display_names


def natural_sort_key(text):
    """Create a key for natural sorting that handles numbers correctly."""
    import re
    def atoi(text):
        return int(text) if text.isdigit() else text
    return [atoi(c) for c in re.split(r'(\d+)', text)]


def group_samples_by_library(json_files):
    """
    Group JSON files by library name.
    
    For long R2 format: Groups all R2 files from the same library together (no R1 exists)
    For paired-end format: Keeps R1 and R2 as separate libraries
    """
    library_groups = defaultdict(list)
    
    # Track what we're detecting for diagnostics
    r1_files = []
    r2_only_files = []  # R2 files that will be treated as library base (no R1 counterpart)
    r1_r2_pairs = []    # R2 files that have R1 counterparts (paired-end)
    
    # First pass: detect all R1 and R2
    all_r1_bases = set()
    all_r2_bases = set()
    
    for json_file in json_files:
        sample = Path(json_file).stem.replace('_cutadapt', '')
        
        # Get the base name without R1/R2
        base_clean = sample.replace('_R2_001', '').replace('_R1_001', '')
        base_clean = base_clean.replace('_R2', '').replace('_R1', '')
        
        if '_R1_' in sample or sample.endswith('_R1'):
            all_r1_bases.add(base_clean)
        if '_R2_' in sample or sample.endswith('_R2'):
            all_r2_bases.add(base_clean)
    
    # Now process files and determine library names
    for json_file in json_files:
        sample = Path(json_file).stem.replace('_cutadapt', '')
        
        # Get the base name to check if it's paired
        base_clean = sample.replace('_R2_001', '').replace('_R1_001', '')
        base_clean = base_clean.replace('_R2', '').replace('_R1', '')
        
        is_r1 = '_R1_' in sample or sample.endswith('_R1')
        is_r2 = '_R2_' in sample or sample.endswith('_R2')
        
        # Determine if this base has both R1 and R2 (paired-end) or just R2 (long R2 format)
        has_paired_r1 = base_clean in all_r1_bases
        
        if is_r1:
            r1_files.append(sample)
        elif is_r2 and has_paired_r1:
            r1_r2_pairs.append(sample)
        elif is_r2 and not has_paired_r1:
            r2_only_files.append(sample)
        
        library, run = extract_library_and_run(sample)
        library_groups[library].append((json_file, run))
    
    # Print diagnostic information
    if r1_files and r1_r2_pairs:
        print(f"\n✓ Detected PAIRED-END format:")
        print(f"  - Found {len(r1_files)} R1 files")
        print(f"  - Found {len(r1_r2_pairs)} R2 files (with R1 counterparts)")
        print(f"  → R1 and R2 will be treated as SEPARATE libraries")
    elif r2_only_files and not r1_files:
        print(f"\n✓ Detected LONG R2 format:")
        print(f"  - Found {len(r2_only_files)} R2 files")
        print(f"  - No R1 files found")
        print(f"  → R2 files will be aggregated by library name")
    elif r2_only_files and r1_files:
        print(f"\n✓ Detected MIXED format:")
        print(f"  - Found {len(r1_files)} R1 files (paired-end)")
        print(f"  - Found {len(r1_r2_pairs)} R2 files (paired-end)")
        print(f"  - Found {len(r2_only_files)} R2 files (long R2, no R1)")
        print(f"  → Will handle both formats appropriately")
    else:
        print(f"\n⚠ WARNING: Unusual file pattern detected")
        print(f"  - R1 files: {len(r1_files)}")
        print(f"  - R2 files (paired): {len(r1_r2_pairs)}")
        print(f"  - R2 files (long): {len(r2_only_files)}")
    
    for library in library_groups:
        library_groups[library].sort(key=lambda x: x[1])
    
    # Sort libraries naturally (handles numbers correctly)
    sorted_libraries = sorted(library_groups.keys(), key=natural_sort_key)
    return {lib: library_groups[lib] for lib in sorted_libraries}


def parse_cutadapt_json(json_file):
    """Extract key adapter statistics from cutadapt JSON."""
    with open(json_file) as f:
        data = json.load(f)
    
    sample_name = Path(json_file).stem.replace('_cutadapt', '')
    total_reads = data['read_counts']['input']
    
    # Check both read1 and read2 adapters (handles both single-end and paired-end)
    adapters = data.get('adapters_read1', []) or data.get('adapters_read2', [])
    
    adapter_stats = []
    for adapter in adapters:
        adapter_stats.append({
            'Sample': sample_name,
            'Adapter': adapter['name'],
            'Matches': adapter['total_matches'],
            'Percent': (adapter['total_matches'] / total_reads) * 100 if total_reads > 0 else 0
        })
    
    return pd.DataFrame(adapter_stats)


def aggregate_library_metrics(json_files_with_runs):
    """Aggregate metrics across all runs for a library."""
    total_reads = 0
    adapter_totals = defaultdict(int)
    reads_with_adapters = 0
    reads_too_short = 0
    reads_passing = 0
    reads_untrimmed = 0
    
    for json_file, _ in json_files_with_runs:
        with open(json_file) as f:
            data = json.load(f)
        
        total_reads += data['read_counts']['input']
        
        # Handle both single-end and paired-end formats
        r1_adapt = data['read_counts'].get('read1_with_adapter')
        r2_adapt = data['read_counts'].get('read2_with_adapter')
        if r1_adapt is not None:
            reads_with_adapters += r1_adapt
        elif r2_adapt is not None:
            reads_with_adapters += r2_adapt
        
        reads_passing += data['read_counts'].get('output', 0)
        
        if 'filtered' in data['read_counts']:
            for key, value in data['read_counts']['filtered'].items():
                if value is not None:
                    if key == 'too_short':
                        reads_too_short += value
                    elif key == 'discard_untrimmed':
                        reads_untrimmed += value
        
        # Check both read1 and read2 adapters
        adapters = data.get('adapters_read1', []) or data.get('adapters_read2', [])
        for adapter in adapters:
            adapter_totals[adapter['name']] += adapter['total_matches']
    
    metrics = {
        'total_reads': total_reads,
        'reads_with_adapters_pct': (reads_with_adapters / total_reads) * 100 if total_reads > 0 else 0,
        'reads_too_short_pct': (reads_too_short / total_reads) * 100 if total_reads > 0 else 0,
        'reads_passing_pct': (reads_passing / total_reads) * 100 if total_reads > 0 else 0,
        'reads_untrimmed_pct': (reads_untrimmed / total_reads) * 100 if total_reads > 0 else 0,
        'adapter_percentages': {
            name: (count / total_reads) * 100 if total_reads > 0 else 0
            for name, count in adapter_totals.items()
        }
    }
    
    return metrics


def create_adapter_summary(json_files):
    """Create summary table across all samples."""
    all_data = []
    for json_file in json_files:
        df = parse_cutadapt_json(json_file)
        if not df.empty:
            all_data.append(df)
    
    if not all_data:
        return pd.DataFrame(), pd.DataFrame()
    
    combined = pd.concat(all_data, ignore_index=True)
    pivot = combined.pivot(index='Sample', columns='Adapter', values='Percent')
    
    return combined, pivot


def calculate_expected_values(json_files):
    """Calculate median values from actual data across all samples."""
    all_metrics = {
        'TSO_5p': [], 'TSO_RC': [], 'Poly_G': [], 'TruSeq': [],
        'With_Adapters': [], 'Too_Short': [], 'Passing': []
    }
    
    for json_file in json_files:
        with open(json_file) as f:
            data = json.load(f)
        
        total = data['read_counts']['input']
        if total == 0:
            continue
        
        # Check both read1 and read2 adapters
        adapters_list = data.get('adapters_read1', []) or data.get('adapters_read2', [])
        adapters = {a['name']: a['total_matches'] for a in adapters_list}    
        
        all_metrics['TSO_5p'].append((adapters.get('^TSO', 0) / total) * 100)
        all_metrics['TSO_RC'].append((adapters.get('TSO_RC', 0) / total) * 100)
        all_metrics['Poly_G'].append((adapters.get('1', 0) / total) * 100)
        
        truseq = adapters.get('TruSeq_Read1_RC', 0) + adapters.get('TruSeq_Read2_Universal', 0)
        all_metrics['TruSeq'].append((truseq / total) * 100)
        
        # Handle both single-end and paired-end
        r1_adapt = data['read_counts'].get('read1_with_adapter')
        r2_adapt = data['read_counts'].get('read2_with_adapter')
        if r1_adapt is not None:
            all_metrics['With_Adapters'].append((r1_adapt / total) * 100)
        elif r2_adapt is not None:
            all_metrics['With_Adapters'].append((r2_adapt / total) * 100)
        
        too_short = 0
        if 'filtered' in data['read_counts'] and 'too_short' in data['read_counts']['filtered']:
            too_short = data['read_counts']['filtered']['too_short'] or 0
        all_metrics['Too_Short'].append((too_short / total) * 100)
        all_metrics['Passing'].append((data['read_counts'].get('output', 0) / total) * 100)
    
    return {k: np.median(v) if v else 0 for k, v in all_metrics.items()}


def plot_adapter_heatmap(pivot_df, output_file):
    """Create heatmap showing adapter frequencies across samples."""
    if pivot_df.empty:
        print("  ⚠ Skipping heatmap - no data")
        return
        
    adapter_order = ['^TSO', '1', 'TSO', 'TSO_RC',
                    'TSO_5prime_marker', 'TSO_5prime_marker_RC',
                    'TruSeq_Read2_Universal', 'TruSeq_Read1_RC', 
                    'TruSeq_Read2_Primer', 'TruSeq_Read2_RC_Primer',
                    'P7_adapter', 'P7_RC', 'P5_adapter_RC',
                    'cDNA_Reverse_V3_V4', 'cDNA_Reverse_5prime',
                    'Library_PCR_P1_Full', 'Library_PCR_P2_Core',
                    'PolyT_RT_Primer_5prime']
    
    # Only include adapters that exist in data AND have non-zero values
    ordered_cols = [col for col in adapter_order 
                    if col in pivot_df.columns and pivot_df[col].sum() > 0]
    if not ordered_cols:
        print("  ⚠ Skipping heatmap - no matching adapters with non-zero values")
        return
        
    pivot_ordered = pivot_df[ordered_cols]
    
    rename_dict = {
        '^TSO': 'TSO (5\' Anchored)',
        'TSO': 'TSO (Unanchored)',
        'TSO_RC': 'TSO_RC (3\')',
        '1': 'Poly-G (TSO artifact)',
        'TruSeq_Read1_RC': 'TruSeq Read1 RC',
        'TruSeq_Read2_Universal': 'TruSeq Read2',
        'TruSeq_Read2_Primer': 'TruSeq Read2 Primer',
        'TruSeq_Read2_RC_Primer': 'TruSeq R2 RC Primer',
        'P7_adapter': 'P7',
        'P7_RC': 'P7 RC',
        'P5_adapter_RC': 'P5 RC',
        'TSO_5prime_marker': 'TSO 5prime marker',
        'TSO_5prime_marker_RC': 'TSO 5prime marker RC',
        'cDNA_Reverse_V3_V4': 'cDNA Reverse V3_V4',
        'cDNA_Reverse_5prime': 'cDNA Reverse 5prime',
        'Library_PCR_P1_Full': 'Library PCR P1 Full',
        'Library_PCR_P2_Core': 'Library PCR P2 Core',
        'PolyT_RT_Primer_5prime': 'PolyT RT Primer 5prime'
    }
    pivot_ordered = pivot_ordered.rename(columns=rename_dict)
    
    plt.figure(figsize=(12, max(8, len(pivot_ordered) * 0.4)))
    sns.heatmap(pivot_ordered, annot=True, fmt='.1f', cmap='YlOrRd',
                cbar_kws={'label': '% of Reads'})
    plt.title('Adapter Detection Frequencies Across All Samples and Runs\n(% of total reads)', 
              fontsize=14, fontweight='bold')
    plt.xlabel('Adapter Type', fontsize=12)
    plt.ylabel('Sample', fontsize=12)
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  ✓ Heatmap saved to {output_file}")


def plot_library_aggregated_heatmap(library_groups, display_names, output_file, use_clustering=True):
    """Create heatmap with libraries aggregated, with optional hierarchical clustering."""
    library_data = []
    library_names_list = []
    
    # Sort libraries naturally (if not already sorted)
    sorted_libraries = sorted(library_groups.keys(), key=natural_sort_key)
    
    for library in sorted_libraries:
        json_files_with_runs = library_groups[library]
        total_reads = 0
        adapter_totals = defaultdict(int)
        
        for json_file, _ in json_files_with_runs:
            with open(json_file) as f:
                data = json.load(f)
            total_reads += data['read_counts']['input']
            # Check both read1 and read2 adapters
            adapters = data.get('adapters_read1', []) or data.get('adapters_read2', [])
            for adapter in adapters:
                adapter_totals[adapter['name']] += adapter['total_matches']
        
        if total_reads == 0:
            continue
            
        library_names_list.append(display_names[library])
        library_row = {}
        for adapter_name, count in adapter_totals.items():
            library_row[adapter_name] = (count / total_reads) * 100
        library_data.append(library_row)
    
    if not library_data:
        print("  ⚠ Skipping library heatmap - no data")
        return
        
    df = pd.DataFrame(library_data, index=library_names_list)
    
    adapter_order = ['^TSO', '1', 'TSO', 'TSO_RC',
                    'TSO_5prime_marker', 'TSO_5prime_marker_RC',
                    'TruSeq_Read2_Universal', 'TruSeq_Read1_RC', 
                    'TruSeq_Read2_Primer', 'TruSeq_Read2_RC_Primer',
                    'P7_adapter', 'P7_RC', 'P5_adapter_RC',
                    'cDNA_Reverse_V3_V4', 'cDNA_Reverse_5prime',
                    'Library_PCR_P1_Full', 'Library_PCR_P2_Core',
                    'PolyT_RT_Primer_5prime']
    
    # Only include adapters that exist in data AND have non-zero values
    ordered_cols = [col for col in adapter_order 
                    if col in df.columns and df[col].sum() > 0]
    
    if not ordered_cols:
        print("  ⚠ Skipping library heatmap - no matching adapters with non-zero values")
        return
        
    df_ordered = df[ordered_cols]                
    
    rename_dict = {
        '^TSO': 'TSO (5\' Anchored)',
        'TSO': 'TSO (Unanchored)',
        'TSO_RC': 'TSO_RC (3\')',
        '1': 'Poly-G (TSO artifact)',
        'TruSeq_Read1_RC': 'TruSeq Read1 RC',
        'TruSeq_Read2_Universal': 'TruSeq Read2',
        'TruSeq_Read2_Primer': 'TruSeq Read2 Primer',
        'TruSeq_Read2_RC_Primer': 'TruSeq R2 RC Primer',
        'P7_adapter': 'P7',
        'P7_RC': 'P7 RC',
        'P5_adapter_RC': 'P5 RC',
        'TSO_5prime_marker': 'TSO 5prime marker',
        'TSO_5prime_marker_RC': 'TSO 5prime marker RC',
        'cDNA_Reverse_V3_V4': 'cDNA Reverse V3_V4',
        'cDNA_Reverse_5prime': 'cDNA Reverse 5prime',
        'Library_PCR_P1_Full': 'Library PCR P1 Full',
        'Library_PCR_P2_Core': 'Library PCR P2 Core',
        'PolyT_RT_Primer_5prime': 'PolyT RT Primer 5prime'
    }
    df_ordered = df_ordered.rename(columns=rename_dict)
    
    # Perform hierarchical clustering if requested and possible
    if use_clustering and len(df_ordered) > 1:  # Need at least 2 samples to cluster
        try:
            # Calculate linkage for rows (libraries)
            row_linkage = hierarchy.linkage(pdist(df_ordered, metric='euclidean'), method='ward')
            
            # Create figure with dendrogram
            fig = plt.figure(figsize=(14, max(8, len(df_ordered) * 0.5)))
            
            # Create grid for dendrogram + heatmap
            gs = fig.add_gridspec(1, 2, width_ratios=[1, 4], wspace=0.01)
            
            # Plot dendrogram
            ax_dendro = fig.add_subplot(gs[0])
            dendro = hierarchy.dendrogram(row_linkage, orientation='left', 
                                         labels=df_ordered.index.tolist(),
                                         ax=ax_dendro, color_threshold=0,
                                         above_threshold_color='gray')
            ax_dendro.set_xticks([])
            ax_dendro.set_ylabel('Library', fontsize=12)
            
            # Reorder dataframe based on dendrogram
            reordered_idx = dendro['leaves']
            df_clustered = df_ordered.iloc[reordered_idx]
            
            # Plot heatmap
            ax_heatmap = fig.add_subplot(gs[1])
            sns.heatmap(df_clustered, annot=True, fmt='.1f', cmap='YlOrRd',
                        cbar_kws={'label': '% of Reads'}, ax=ax_heatmap,
                        yticklabels=True)
            ax_heatmap.set_xlabel('Adapter Type', fontsize=12)
            ax_heatmap.set_ylabel('')  # Y-label already on dendrogram
            ax_heatmap.tick_params(axis='y', left=False)  # Remove y-axis ticks
            
            # Rotate x-labels
            plt.setp(ax_heatmap.get_xticklabels(), rotation=45, ha='right')
            
            plt.suptitle('Adapter Detection Frequencies by Library (Hierarchically Clustered)\n(Aggregated across all sequencing runs)', 
                        fontsize=14, fontweight='bold', y=0.98)
            
            plt.subplots_adjust(left=0.05, right=0.98, top=0.94, bottom=0.1)
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            plt.close()
            print(f"  ✓ Library-aggregated clustered heatmap saved to {output_file}")
            
        except Exception as e:
            print(f"  ⚠ Clustering failed, falling back to sorted heatmap: {str(e)}")
            use_clustering = False
    
    # If not clustering (or clustering failed), create regular sorted heatmap
    if not use_clustering or len(df_ordered) <= 1:
        plt.figure(figsize=(12, max(8, len(df_ordered) * 0.5)))
        sns.heatmap(df_ordered, annot=True, fmt='.1f', cmap='YlOrRd',
                    cbar_kws={'label': '% of Reads'})
        title_suffix = '' if use_clustering else ' (Sorted by Library Number)'
        plt.title(f'Adapter Detection Frequencies by Library{title_suffix}\n(Aggregated across all sequencing runs)', 
                  fontsize=14, fontweight='bold')
        plt.xlabel('Adapter Type', fontsize=12)
        plt.ylabel('Library', fontsize=12)
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  ✓ Library-aggregated heatmap saved to {output_file}")


def plot_library_other_adapters_breakdown(library_groups, display_names, output_file):
    """Create detailed breakdown of adapters - 2x2 grid with normalized and non-normalized versions."""
    
    libraries = []
    other_adapters = {}
    all_adapters = {}  # For the plots with all categories
    
    # Define the main adapters we already show
    main_adapters = ['^TSO', 'TSO_RC', '1', 'TruSeq_Read1_RC', 'TruSeq_Read2_Universal']
    
    # Sort libraries naturally
    sorted_libraries = sorted(library_groups.keys(), key=natural_sort_key)
    
    # Aggregate by library
    for library in sorted_libraries:
        json_files_with_runs = library_groups[library]
        display_name = display_names[library]
        libraries.append(display_name)
        
        total_reads = 0
        adapter_totals = defaultdict(int)
        all_adapter_totals = defaultdict(int)  # For all adapters including main ones
        
        # Aggregate across all runs for this library
        for json_file, _ in json_files_with_runs:
            with open(json_file) as f:
                data = json.load(f)
            
            total_reads += data['read_counts']['input']
            # Check both read1 and read2 adapters
            adapters = data.get('adapters_read1', []) or data.get('adapters_read2', [])
            for adapter in adapters:
                name = adapter['name']
                all_adapter_totals[name] += adapter['total_matches']  # Collect all adapters
                if name not in main_adapters:
                    adapter_totals[name] += adapter['total_matches']
              
        # Calculate percentages for minor adapters
        for name, count in adapter_totals.items():
            if name not in other_adapters:
                other_adapters[name] = []
            pct = (count / total_reads) * 100 if total_reads > 0 else 0
            other_adapters[name].append(pct)
        
        # Calculate percentages for all adapters
        # Group TSO 5', TSO_RC 3', Poly-G, and TruSeq adapters
        tso_5p_count = all_adapter_totals.get('^TSO', 0)
        tso_rc_count = all_adapter_totals.get('TSO_RC', 0)
        polyg_count = all_adapter_totals.get('1', 0)
        truseq_count = all_adapter_totals.get('TruSeq_Read1_RC', 0) + all_adapter_totals.get('TruSeq_Read2_Universal', 0)
        
        # Add main categories to all_adapters dictionary
        if 'TSO 5\'' not in all_adapters:
            all_adapters['TSO 5\''] = []
        all_adapters['TSO 5\''].append((tso_5p_count / total_reads) * 100 if total_reads > 0 else 0)
        
        if 'TSO_RC 3\'' not in all_adapters:
            all_adapters['TSO_RC 3\''] = []
        all_adapters['TSO_RC 3\''].append((tso_rc_count / total_reads) * 100 if total_reads > 0 else 0)
        
        if 'Poly-G' not in all_adapters:
            all_adapters['Poly-G'] = []
        all_adapters['Poly-G'].append((polyg_count / total_reads) * 100 if total_reads > 0 else 0)
        
        if 'TruSeq' not in all_adapters:
            all_adapters['TruSeq'] = []
        all_adapters['TruSeq'].append((truseq_count / total_reads) * 100 if total_reads > 0 else 0)
        
        # Add other adapters (minor ones) as a group
        other_total = sum(count for name, count in adapter_totals.items())
        if 'Other Adapters' not in all_adapters:
            all_adapters['Other Adapters'] = []
        all_adapters['Other Adapters'].append((other_total / total_reads) * 100 if total_reads > 0 else 0)
    
    if not other_adapters:
        print("  ⚠ No minor adapters found")
        # Still create the plot with just main adapters
    
    # Adjust figure size based on number of libraries
    num_libs = len(libraries)
    fig_width = max(20, num_libs * 0.6)  # Scale width with number of libraries
    fig_height = 12  # Fixed height for 2x2 grid
    
    # Create figure with 2x2 subplots
    fig, axes = plt.subplots(2, 2, figsize=(fig_width, fig_height))
    
    # Prepare data for all adapters
    all_adapter_names_full = ['TSO 5\'', 'TSO_RC 3\'', 'Poly-G', 'TruSeq', 'Other Adapters']
    all_adapter_data_full = np.array([all_adapters.get(name, [0]*len(libraries)) for name in all_adapter_names_full])
    
    # Filter out adapters with all zeros
    non_zero_mask = [all_adapter_data_full[i].sum() > 0 for i in range(len(all_adapter_names_full))]
    all_adapter_names = [name for name, keep in zip(all_adapter_names_full, non_zero_mask) if keep]
    all_adapter_data = all_adapter_data_full[non_zero_mask]
    
    if len(all_adapter_names) == 0:
        print("  ⚠ No adapters detected in any library")
        return
    
    # Normalize data for normalized plots
    normalized_all_data = []
    for i in range(len(all_adapter_names)):
        normalized_row = []
        for j in range(len(libraries)):
            total_adapter_pct = sum(all_adapter_data[k][j] for k in range(len(all_adapter_names)))
            if total_adapter_pct > 0:
                normalized_value = (all_adapter_data[i][j] / total_adapter_pct) * 100
            else:
                normalized_value = 0
            normalized_row.append(normalized_value)
        normalized_all_data.append(normalized_row)
    normalized_all_data = np.array(normalized_all_data)
    
    # Colors for all adapter categories - use consistent mapping
    colors_all = [get_adapter_color(name, idx) for idx, name in enumerate(all_adapter_names)]
    
    # TOP LEFT: Adapter composition normalized to 100%
    ax1 = axes[0, 0]
    bottom = np.zeros(len(libraries))
    for i, (name, color) in enumerate(zip(all_adapter_names, colors_all)):
        ax1.bar(range(len(libraries)), normalized_all_data[i], bottom=bottom, 
               label=name, alpha=0.85, color=color, edgecolor='black', linewidth=0.5)
        bottom += normalized_all_data[i]
    
    ax1.set_ylabel('% of Adapter-Containing Reads', fontsize=11)
    ax1.set_title('Adapter Composition (Normalized to 100%)', fontsize=12, fontweight='bold')
    ax1.set_xticks(range(len(libraries)))
    ax1.set_ylim([0, 100])
    
    # Adjust labels
    if num_libs > 20:
        ax1.set_xticklabels(libraries, rotation=90, ha='center', fontsize=max(6, 10 - num_libs//10))
    else:
        ax1.set_xticklabels(libraries, rotation=45, ha='right', fontsize=9)
    
    # Legend outside the plot
    ax1.legend(loc='upper left', bbox_to_anchor=(1.02, 1), fontsize=9)
    ax1.grid(True, alpha=0.3, axis='y')
    
    # TOP RIGHT: Adapter composition by % of total reads (not normalized)
    ax2 = axes[0, 1]
    bottom = np.zeros(len(libraries))
    for i, (name, color) in enumerate(zip(all_adapter_names, colors_all)):
        ax2.bar(range(len(libraries)), all_adapter_data[i], bottom=bottom, 
               label=name, alpha=0.85, color=color, edgecolor='black', linewidth=0.5)
        bottom += all_adapter_data[i]
    
    ax2.set_ylabel('% of Total Reads', fontsize=11)
    ax2.set_title('Adapter Composition (% of Total Reads)', fontsize=12, fontweight='bold')
    ax2.set_xticks(range(len(libraries)))
    
    # Adjust labels
    if num_libs > 20:
        ax2.set_xticklabels(libraries, rotation=90, ha='center', fontsize=max(6, 10 - num_libs//10))
    else:
        ax2.set_xticklabels(libraries, rotation=45, ha='right', fontsize=9)
    
    # Legend outside the plot
    ax2.legend(loc='upper left', bbox_to_anchor=(1.02, 1), fontsize=9)
    ax2.grid(True, alpha=0.3, axis='y')
    
    # Process minor adapters data
    if other_adapters:
        # Filter out adapters with all zeros
        adapter_names = []
        adapter_data_list = []
        for name, data in other_adapters.items():
            if sum(data) > 0:  # Only include if there's at least some non-zero value
                adapter_names.append(name)
                adapter_data_list.append(data)
        
        if not adapter_names:
            # All minor adapters are zero
            other_adapters = None
        else:
            adapter_data = np.array(adapter_data_list)
        
    if other_adapters:
        # Rename for clarity
        rename_dict = {
            '^TSO': 'TSO (5\' Anchored)',
            'TSO': 'TSO (Unanchored)',
            'TSO_RC': 'TSO_RC (3\')',
            '1': 'Poly-G (TSO artifact)',
            'P7_adapter': 'P7',
            'P7_RC': 'P7 RC',
            'P5_adapter_RC': 'P5 RC',
            'TruSeq_Read2_RC_Primer': 'TruSeq R2 RC Primer',
            'TruSeq_Read2_Primer': 'TruSeq Read2 Primer',
            'TruSeq_Read2_Universal': 'TruSeq Read2',
            'TruSeq_Read1_RC': 'TruSeq Read1 RC',
            'TSO_5prime_marker': 'TSO 5prime marker',
            'TSO_5prime_marker_RC': 'TSO 5prime marker RC',
            'cDNA_Reverse_V3_V4': 'cDNA Reverse V3_V4',
            'cDNA_Reverse_5prime': 'cDNA Reverse 5prime',
            'Library_PCR_P1_Full': 'Library PCR P1 Full',
            'Library_PCR_P2_Core': 'Library PCR P2 Core',
            'PolyT_RT_Primer_5prime': 'PolyT RT Primer 5prime'
        }
        display_adapter_names = [rename_dict.get(name, name) for name in adapter_names]
        
        # Colors for minor adapters - use consistent mapping
        colors_minor = [get_adapter_color(name, idx) for idx, name in enumerate(adapter_names)]
        
        # Normalize minor adapter data
        normalized_minor_data = []
        for i in range(len(adapter_names)):
            normalized_row = []
            for j in range(len(libraries)):
                total_minor_pct = sum(adapter_data[k][j] for k in range(len(adapter_names)))
                if total_minor_pct > 0:
                    normalized_value = (adapter_data[i][j] / total_minor_pct) * 100
                else:
                    normalized_value = 0
                normalized_row.append(normalized_value)
            normalized_minor_data.append(normalized_row)
        normalized_minor_data = np.array(normalized_minor_data)
        
        # BOTTOM LEFT: Minor adapter breakdown normalized to 100%
        ax3 = axes[1, 0]
        bottom = np.zeros(len(libraries))
        for i, (name, color) in enumerate(zip(display_adapter_names, colors_minor)):
            ax3.bar(range(len(libraries)), normalized_minor_data[i], bottom=bottom, 
                   label=name, alpha=0.85, color=color, edgecolor='black', linewidth=0.5)
            bottom += normalized_minor_data[i]
        
        ax3.set_ylabel('% of Minor Adapter-Containing Reads', fontsize=11)
        ax3.set_title('Minor Adapter Breakdown (Normalized to 100%)', fontsize=12, fontweight='bold')
        ax3.set_ylim([0, 100])
        ax3.legend(loc='upper left', bbox_to_anchor=(1.02, 1), fontsize=9)
        
        # BOTTOM RIGHT: Minor adapter breakdown by % of total reads
        ax4 = axes[1, 1]
        bottom = np.zeros(len(libraries))
        for i, (name, color) in enumerate(zip(display_adapter_names, colors_minor)):
            ax4.bar(range(len(libraries)), adapter_data[i], bottom=bottom, 
                   label=name, alpha=0.85, color=color, edgecolor='black', linewidth=0.5)
            bottom += adapter_data[i]
        
        ax4.set_ylabel('% of Total Reads', fontsize=11)
        ax4.set_title('Minor Adapter Breakdown (% of Total Reads)', fontsize=12, fontweight='bold')
        ax4.legend(loc='upper left', bbox_to_anchor=(1.02, 1), fontsize=9)
    else:
        # If no minor adapters, show messages
        ax3 = axes[1, 0]
        ax3.text(0.5, 0.5, 'No minor adapters detected', 
                ha='center', va='center', fontsize=14, transform=ax3.transAxes)
        ax3.set_title('Minor Adapter Breakdown (Normalized to 100%)', fontsize=12, fontweight='bold')
        
        ax4 = axes[1, 1]
        ax4.text(0.5, 0.5, 'No minor adapters detected', 
                ha='center', va='center', fontsize=14, transform=ax4.transAxes)
        ax4.set_title('Minor Adapter Breakdown (% of Total Reads)', fontsize=12, fontweight='bold')
    
    # Set x-axis for bottom plots
    for ax in [ax3, ax4]:
        ax.set_xticks(range(len(libraries)))
        if num_libs > 20:
            ax.set_xticklabels(libraries, rotation=90, ha='center', fontsize=max(6, 10 - num_libs//10))
        else:
            ax.set_xticklabels(libraries, rotation=45, ha='right', fontsize=9)
        ax.grid(True, alpha=0.3, axis='y')
    
    plt.suptitle('Comprehensive Adapter Breakdown by Library', 
                 fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  ✓ Library-aggregated other adapters breakdown saved to {output_file}")


def plot_library_multi_comparison(library_groups, display_names, output_file):
    """Create comprehensive comparison across all libraries - LIBRARY AGGREGATED VERSION."""
    
    libraries = []
    metrics = {
        'Total Reads': [],
        'Reads with Adapters (%)': [],
        'Reads Too Short (%)': [],
        'Reads Passing (%)': [],
        'TSO 5\' (%)': [],
        'TSO_RC 3\' (%)': [],
        'Poly-G (%)': [],
        'TruSeq Adapters (%)': []
    }
    
    # Sort libraries naturally 
    sorted_libraries = sorted(library_groups.keys(), key=natural_sort_key)
    
    # Aggregate metrics by library
    for library in sorted_libraries:
        json_files_with_runs = library_groups[library]
        display_name = display_names[library]
        libraries.append(display_name)
        
        agg_metrics = aggregate_library_metrics(json_files_with_runs)
        adapter_pcts = agg_metrics['adapter_percentages']
        
        metrics['Total Reads'].append(agg_metrics['total_reads'] / 1e6)
        metrics['Reads with Adapters (%)'].append(agg_metrics['reads_with_adapters_pct'])
        metrics['Reads Too Short (%)'].append(agg_metrics['reads_too_short_pct'])
        metrics['Reads Passing (%)'].append(agg_metrics['reads_passing_pct'])
        
        metrics['TSO 5\' (%)'].append(adapter_pcts.get('^TSO', 0))
        metrics['TSO_RC 3\' (%)'].append(adapter_pcts.get('TSO_RC', 0))
        metrics['Poly-G (%)'].append(adapter_pcts.get('1', 0))
        
        truseq = adapter_pcts.get('TruSeq_Read1_RC', 0) + adapter_pcts.get('TruSeq_Read2_Universal', 0)
        metrics['TruSeq Adapters (%)'].append(truseq)
    
    if not libraries:
        print("  ⚠ No libraries to compare")
        return
    
    num_libs = len(libraries)
    # Adjust font size based on number of libraries
    if num_libs > 30:
        label_fontsize = 7
        rotation = 90
        ha = 'center'
    elif num_libs > 20:
        label_fontsize = 8
        rotation = 90
        ha = 'center'
    else:
        label_fontsize = 9
        rotation = 45
        ha = 'right'
    
    # Calculate expected values (medians)
    expected = {k: np.median(v) for k, v in metrics.items() if k not in ['Total Reads'] and v}
    
    # Create figure with subplots
    fig = plt.figure(figsize=(20, 14))
    gs = fig.add_gridspec(3, 3, hspace=0.5, wspace=0.4)
    
    # Plot 1: Total reads processed
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.bar(range(len(libraries)), metrics['Total Reads'], alpha=0.7, color='#4682b4', edgecolor='black')  # steelblue for non-adapter metric
    ax1.set_ylabel('Millions of Reads', fontsize=11)
    ax1.set_title('Total Reads Processed by Library', fontsize=12, fontweight='bold')
    ax1.set_xticks(range(len(libraries)))
    ax1.set_xticklabels(libraries, rotation=rotation, ha=ha, fontsize=label_fontsize)
    ax1.grid(True, alpha=0.3, axis='y')
    
    # Plot 2: Reads with adapters
    ax2 = fig.add_subplot(gs[0, 1])
    adapters_key = 'Reads with Adapters (%)'
    ax2.bar(range(len(libraries)), metrics[adapters_key], alpha=0.7, color='#ff7f50', edgecolor='black')  # coral for general adapter metric
    if adapters_key in expected:
        ax2.axhline(expected[adapters_key], color='red', linestyle='--', alpha=0.5, 
                    label=f'Median: {expected[adapters_key]:.1f}%')
    ax2.set_ylabel('% of Reads', fontsize=11)
    ax2.set_title('Reads with Adapters Detected', fontsize=12, fontweight='bold')
    ax2.set_xticks(range(len(libraries)))
    ax2.set_xticklabels(libraries, rotation=rotation, ha=ha, fontsize=label_fontsize)
    ax2.legend()
    ax2.grid(True, alpha=0.3, axis='y')
    
    # Plot 3: Reads passing filters
    ax3 = fig.add_subplot(gs[0, 2])
    passing_key = 'Reads Passing (%)'
    colors = ['green' if x > 80 else 'orange' if x > 70 else 'red' 
              for x in metrics[passing_key]]
    ax3.bar(range(len(libraries)), metrics[passing_key], alpha=0.7, color=colors, edgecolor='black')
    if passing_key in expected:
        ax3.axhline(expected[passing_key], color='green', linestyle='--', alpha=0.5, 
                    label=f'Median: {expected[passing_key]:.1f}%')
    ax3.set_ylabel('% of Reads', fontsize=11)
    ax3.set_title('Reads Passing QC Filters', fontsize=12, fontweight='bold')
    ax3.set_xticks(range(len(libraries)))
    ax3.set_xticklabels(libraries, rotation=rotation, ha=ha, fontsize=label_fontsize)
    ax3.legend()
    ax3.grid(True, alpha=0.3, axis='y')
    ax3.set_ylim([0, 100])
    
    # Plot 4: TSO contamination (apoptotic marker)
    ax4 = fig.add_subplot(gs[1, 0])
    tso_5p_key = 'TSO 5\' (%)'
    if tso_5p_key in expected and expected[tso_5p_key] > 0:
        colors_tso = ['red' if x > expected[tso_5p_key] * 1.5 else 
                      'orange' if x > expected[tso_5p_key] * 1.2 else 'green' 
                      for x in metrics[tso_5p_key]]
    else:
        colors_tso = 'green'
    ax4.bar(range(len(libraries)), metrics[tso_5p_key], alpha=0.7, color=colors_tso, edgecolor='black')
    if tso_5p_key in expected:
        ax4.axhline(expected[tso_5p_key], color='blue', linestyle='--', alpha=0.5, 
                    label=f'Median: {expected[tso_5p_key]:.1f}%')
    ax4.set_ylabel('% of Reads', fontsize=11)
    ax4.set_title('TSO at 5\' End (Apoptotic Cells)', fontsize=12, fontweight='bold')
    ax4.set_xticks(range(len(libraries)))
    ax4.set_xticklabels(libraries, rotation=rotation, ha=ha, fontsize=label_fontsize)
    ax4.legend()
    ax4.grid(True, alpha=0.3, axis='y')
    
    # Plot 5: TSO_RC at 3' end
    ax5 = fig.add_subplot(gs[1, 1])
    tso_rc_key = 'TSO_RC 3\' (%)'
    ax5.bar(range(len(libraries)), metrics[tso_rc_key], alpha=0.7, color=get_adapter_color('TSO_RC'), edgecolor='black')
    if tso_rc_key in expected:
        ax5.axhline(expected[tso_rc_key], color='blue', linestyle='--', alpha=0.5, 
                    label=f'Median: {expected[tso_rc_key]:.3f}%')
    ax5.set_ylabel('% of Reads', fontsize=11)
    ax5.set_title('TSO_RC at 3\' End', fontsize=12, fontweight='bold')
    ax5.set_xticks(range(len(libraries)))
    ax5.set_xticklabels(libraries, rotation=rotation, ha=ha, fontsize=label_fontsize)
    ax5.legend()
    ax5.grid(True, alpha=0.3, axis='y')
    
    # Plot 6: Poly-G contamination
    ax6 = fig.add_subplot(gs[1, 2])
    polyg_key = 'Poly-G (%)'
    if polyg_key in expected and expected[polyg_key] > 0:
        colors_polyg = ['red' if x > expected[polyg_key] * 1.5 else 
                        'orange' if x > expected[polyg_key] * 1.2 else 'green' 
                        for x in metrics[polyg_key]]
    else:
        colors_polyg = 'green'
    ax6.bar(range(len(libraries)), metrics[polyg_key], alpha=0.7, color=colors_polyg, edgecolor='black')
    if polyg_key in expected:
        ax6.axhline(expected[polyg_key], color='blue', linestyle='--', alpha=0.5, 
                    label=f'Median: {expected[polyg_key]:.1f}%')
    ax6.set_ylabel('% of Reads', fontsize=11)
    ax6.set_title('Poly-G (TSO Artifact)', fontsize=12, fontweight='bold')
    ax6.set_xticks(range(len(libraries)))
    ax6.set_xticklabels(libraries, rotation=rotation, ha=ha, fontsize=label_fontsize)
    ax6.legend()
    ax6.grid(True, alpha=0.3, axis='y')
    
    # Plot 7: TruSeq adapters (read-through)
    ax7 = fig.add_subplot(gs[2, 0])
    truseq_key = 'TruSeq Adapters (%)'
    ax7.bar(range(len(libraries)), metrics[truseq_key], alpha=0.7, color=get_adapter_color('TruSeq_Read2_Universal'), edgecolor='black')
    if truseq_key in expected:
        ax7.axhline(expected[truseq_key], color='blue', linestyle='--', alpha=0.5, 
                    label=f'Median: {expected[truseq_key]:.1f}%')
    ax7.set_ylabel('% of Reads', fontsize=11)
    ax7.set_title('TruSeq Adapters (Read-through)', fontsize=12, fontweight='bold')
    ax7.set_xticks(range(len(libraries)))
    ax7.set_xticklabels(libraries, rotation=rotation, ha=ha, fontsize=label_fontsize)
    ax7.legend()
    ax7.grid(True, alpha=0.3, axis='y')
    
    # Plot 8: Stacked bar chart of adapter types
    ax8 = fig.add_subplot(gs[2, 1:])
    
    tso_5p_key = 'TSO 5\' (%)'
    tso_rc_key = 'TSO_RC 3\' (%)'
    polyg_key = 'Poly-G (%)'
    truseq_key = 'TruSeq Adapters (%)'
    
    adapter_types = ['TSO 5\'', 'TSO_RC 3\'', 'Poly-G', 'TruSeq', 'No Major Adapters']
    adapter_data = np.array([
        metrics[tso_5p_key],
        metrics[tso_rc_key],
        metrics[polyg_key],
        metrics[truseq_key],
        [max(0, 100 - sum(x)) for x in zip(metrics[tso_5p_key], metrics[tso_rc_key], 
                                   metrics[polyg_key], metrics[truseq_key])]
    ])
    
    # Use consistent color mapping
    colors_stack = [
        get_adapter_color('TSO', 0),
        get_adapter_color('TSO_RC', 1),
        get_adapter_color('Poly-G', 2),
        get_adapter_color('TruSeq_Read2_Universal', 3),
        '#bdbdbd'  # Gray for "No Major Adapters"
    ]
    bottom = np.zeros(len(libraries))
    
    for i, (adapter_type, color) in enumerate(zip(adapter_types, colors_stack)):
        ax8.bar(range(len(libraries)), adapter_data[i], bottom=bottom, 
               label=adapter_type, alpha=0.8, color=color, edgecolor='black', linewidth=0.5)
        bottom += adapter_data[i]
    
    ax8.set_ylabel('% of Total Reads', fontsize=11)
    ax8.set_title('Adapter Composition Breakdown by Library', fontsize=12, fontweight='bold')
    ax8.set_xticks(range(len(libraries)))
    ax8.set_xticklabels(libraries, rotation=rotation, ha=ha, fontsize=label_fontsize)
    ax8.legend(loc='upper left', bbox_to_anchor=(1, 1))
    ax8.set_ylim([0, 100])
    ax8.grid(True, alpha=0.3, axis='y')
    
    plt.suptitle('Multi-Library Adapter Trimming Comparison\n(Aggregated across all sequencing runs)', 
                 fontsize=16, fontweight='bold', y=0.995)
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  ✓ Library-aggregated multi-comparison saved to {output_file}")


def create_summary_table(pivot_df, output_file):
    """Create summary statistics table."""
    if pivot_df.empty:
        print(f"  ⚠ Skipping summary table - no data")
        return pd.DataFrame()
        
    summary = pivot_df.describe().T
    summary['median_%'] = pivot_df.median()
    summary['mean_%'] = summary['mean']
    summary['std_%'] = summary['std']
    summary['min_%'] = summary['min']
    summary['max_%'] = summary['max']
    summary = summary[['median_%', 'mean_%', 'std_%', 'min_%', 'max_%']].round(3)
    summary.to_csv(output_file)
    print(f"  ✓ Summary statistics saved to {output_file}")
    return summary


def plot_library_aggregated_summary(library_name, display_name, agg_metrics, json_files_with_runs, output_dir):
    """Create summary visualization for aggregated library metrics."""
    fig, axes = plt.subplots(2, 2, figsize=(16, 10))
    
    adapter_pcts = agg_metrics['adapter_percentages']
    total_reads = agg_metrics['total_reads']
    
    reads_passing = int(total_reads * agg_metrics['reads_passing_pct'] / 100)
    reads_too_short = int(total_reads * agg_metrics['reads_too_short_pct'] / 100)
    reads_untrimmed = int(total_reads * agg_metrics['reads_untrimmed_pct'] / 100)
    
    # Plot 1: Pie chart
    ax1 = axes[0, 0]
    sizes, labels, colors = [], [], []
    if reads_passing > 0:
        sizes.append((reads_passing / total_reads) * 100)
        labels.append(f'Trimmed & Passed\n{(reads_passing / total_reads) * 100:.1f}%')
        colors.append('#2ca02c')
    if reads_untrimmed > 0:
        sizes.append((reads_untrimmed / total_reads) * 100)
        labels.append(f'Untrimmed\n{(reads_untrimmed / total_reads) * 100:.1f}%')
        colors.append('#1f77b4')
    if reads_too_short > 0:
        sizes.append((reads_too_short / total_reads) * 100)
        labels.append(f'Too Short\n{(reads_too_short / total_reads) * 100:.1f}%')
        colors.append('#d62728')
    
    if sizes:
        ax1.pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%', startangle=90)
    ax1.set_title(f'Read Fate Distribution\n(Total: {total_reads/1e6:.1f}M reads)', 
                  fontsize=12, fontweight='bold')
    
    # Plot 2: Adapter rates
    ax2 = axes[0, 1]
    main_adapters = {
        'TSO 5\'': adapter_pcts.get('^TSO', 0),
        'TSO_RC 3\'': adapter_pcts.get('TSO_RC', 0),
        'Poly-G': adapter_pcts.get('1', 0),
        'TruSeq R1': adapter_pcts.get('TruSeq_Read1_RC', 0),
        'TruSeq R2': adapter_pcts.get('TruSeq_Read2_Universal', 0)
    }
    # Use consistent colors for each adapter type
    adapter_color_map = {
        'TSO 5\'': get_adapter_color('^TSO'),
        'TSO_RC 3\'': get_adapter_color('TSO_RC'),
        'Poly-G': get_adapter_color('1'),
        'TruSeq R1': get_adapter_color('TruSeq_Read1_RC'),
        'TruSeq R2': get_adapter_color('TruSeq_Read2_Universal')
    }
    colors_for_bars = [adapter_color_map[name] for name in main_adapters.keys()]
    ax2.barh(list(main_adapters.keys()), list(main_adapters.values()), 
             alpha=0.7, color=colors_for_bars, edgecolor='black')
    ax2.set_xlabel('% of Total Reads', fontsize=11)
    ax2.set_title('Main Adapter Detection Rates', fontsize=12, fontweight='bold')
    ax2.grid(True, alpha=0.3, axis='x')
    for i, (name, val) in enumerate(main_adapters.items()):
        ax2.text(val + 0.1, i, f'{val:.2f}%', va='center', fontsize=10)
    
    # Plot 3: Summary table
    ax3 = axes[1, 0]
    ax3.axis('tight')
    ax3.axis('off')
    
    table_data = [
        ['Metric', 'Count', 'Percentage'],
        ['Total Reads', f'{int(total_reads):,}', '100.0%'],
        ['Trimmed & Passed QC', f'{reads_passing:,}', f'{(reads_passing/total_reads)*100:.2f}%'],
        ['Untrimmed', f'{reads_untrimmed:,}', f'{(reads_untrimmed/total_reads)*100:.2f}%'],
        ['Too Short', f'{reads_too_short:,}', f'{(reads_too_short/total_reads)*100:.2f}%'],
        ['─' * 20, '─' * 12, '─' * 12],
        ['Sequencing Runs', f'{len(json_files_with_runs)}', '─']
    ]
    
    table = ax3.table(cellText=table_data, cellLoc='left', loc='center',
                     colWidths=[0.5, 0.25, 0.25])
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 1.8)
    
    for i in range(3):
        table[(0, i)].set_facecolor('#4CAF50')
        table[(0, i)].set_text_props(weight='bold', color='white')
        table[(1, i)].set_facecolor('#E8F5E9')
        table[(1, i)].set_text_props(weight='bold')
        table[(5, i)].set_facecolor('#F5F5F5')
    
    ax3.set_title('Read Processing Summary', fontsize=12, fontweight='bold', pad=20)
    
    # Plot 4: Adapter composition
    ax4 = axes[1, 1]
    categories = {
        'TSO (5\' & 3\')': adapter_pcts.get('^TSO', 0) + adapter_pcts.get('TSO_RC', 0),
        'Poly-G': adapter_pcts.get('1', 0),
        'TruSeq': adapter_pcts.get('TruSeq_Read1_RC', 0) + adapter_pcts.get('TruSeq_Read2_Universal', 0),
        'No Major Adapters': sum([v for k, v in adapter_pcts.items() 
                     if k not in ['^TSO', 'TSO_RC', '1', 'TruSeq_Read1_RC', 'TruSeq_Read2_Universal']])
    }
    
    colors_cat = ['#d62728', '#2ca02c', '#9467bd', '#8c564b']
    wedges, texts, autotexts = ax4.pie(list(categories.values()), labels=list(categories.keys()), 
                                         colors=colors_cat, autopct='%1.2f%%', startangle=90)
    ax4.set_title('Adapter Composition', fontsize=12, fontweight='bold')
    for autotext in autotexts:
        autotext.set_color('white')
        autotext.set_fontweight('bold')
        autotext.set_fontsize(10)
    
    plt.suptitle(f'{display_name}\nAggregated Summary (All Runs Combined)', 
                 fontsize=16, fontweight='bold', y=0.995)
    plt.tight_layout()
    
    # Use safe filename (replace any problematic characters)
    safe_filename = library_name.replace('/', '_').replace('\\', '_')
    output_file = output_dir / f'{safe_filename}_aggregated_summary.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()


def plot_library_read_length_distribution(library_name, display_name, json_files_with_runs, output_dir):
    """Plot aggregated read length distribution for a library."""
    all_lengths = {}
    
    for json_file, _ in json_files_with_runs:
        with open(json_file) as f:
            data = json.load(f)
        # Check both read1 and read2 adapters
        adapters = data.get('adapters_read1', []) or data.get('adapters_read2', [])
        for adapter in adapters:
            if adapter.get('three_prime_end'):
                trimmed = adapter['three_prime_end'].get('trimmed_lengths', [])
                for item in trimmed:
                    length = item['len']
                    count = sum(item['counts'])
                    if length not in all_lengths:
                        all_lengths[length] = 0
                    all_lengths[length] += count
    
    if not all_lengths:
        print(f"  ⚠ No trimmed length data for {display_name}")
        return
    
    lengths = sorted(all_lengths.keys())
    counts = [all_lengths[l] for l in lengths]
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    
    # Histogram
    ax1.bar(lengths, counts, alpha=0.7, color='steelblue', edgecolor='black', linewidth=0.5)
    ax1.axvline(25, color='red', linestyle='--', linewidth=2, label='Min length cutoff (25bp)')
    ax1.set_xlabel('Trimmed Length (bp)', fontsize=12)
    ax1.set_ylabel('Number of Reads', fontsize=12)
    ax1.set_title('Distribution of Adapter Trim Lengths', fontsize=14, fontweight='bold')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Cumulative
    cumulative = np.cumsum(counts)
    cumulative_pct = (cumulative / cumulative[-1]) * 100
    ax2.plot(lengths, cumulative_pct, linewidth=2, color='darkgreen')
    ax2.axhline(50, color='red', linestyle='--', alpha=0.5, label='50th percentile')
    ax2.axhline(90, color='orange', linestyle='--', alpha=0.5, label='90th percentile')
    ax2.axvline(25, color='red', linestyle='--', linewidth=2, label='Min length cutoff')
    ax2.set_xlabel('Trimmed Length (bp)', fontsize=12)
    ax2.set_ylabel('Cumulative % of Reads', fontsize=12)
    ax2.set_title('Cumulative Distribution', fontsize=14, fontweight='bold')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim([0, 100])
    
    plt.suptitle(f'{display_name} - Read Length Distribution (Aggregated)', 
                 fontsize=16, fontweight='bold')
    plt.tight_layout()
    
    safe_filename = library_name.replace('/', '_').replace('\\', '_')
    plt.savefig(output_dir / f'{safe_filename}_read_lengths.png', dpi=300, bbox_inches='tight')
    plt.close()


def plot_library_tso_analysis(library_name, display_name, json_files_with_runs, output_dir):
    """Plot aggregated TSO position distribution for a library."""
    tso_5p_lengths, tso_5p_counts = [], []
    tso_rc_lengths, tso_rc_counts = [], []
    tso_rc_adj_bases = defaultdict(int)
    
    for json_file, _ in json_files_with_runs:
        with open(json_file) as f:
            data = json.load(f)
        # Check both read1 and read2 adapters
        adapters = data.get('adapters_read1', []) or data.get('adapters_read2', [])
        for adapter in adapters:
            if adapter['name'] == '^TSO' and adapter.get('five_prime_end'):
                trimmed = adapter['five_prime_end'].get('trimmed_lengths', [])
                for item in trimmed:
                    tso_5p_lengths.append(item['len'])
                    tso_5p_counts.append(sum(item['counts']))
            elif adapter['name'] == 'TSO_RC' and adapter.get('three_prime_end'):
                trimmed = adapter['three_prime_end'].get('trimmed_lengths', [])
                for item in trimmed:
                    tso_rc_lengths.append(item['len'])
                    tso_rc_counts.append(sum(item['counts']))
                adj_bases = adapter['three_prime_end'].get('adjacent_bases', {})
                for base, count in adj_bases.items():
                    tso_rc_adj_bases[base] += count
    
    if not tso_5p_lengths and not tso_rc_lengths:
        print(f"  ⚠ No TSO data for {display_name}")
        return
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    # 5' TSO
    if tso_5p_lengths:
        axes[0, 0].bar(tso_5p_lengths, tso_5p_counts, alpha=0.7, color=get_adapter_color('^TSO'))
        axes[0, 0].set_xlabel('Length of TSO Match (bp)', fontsize=12)
        axes[0, 0].set_ylabel('Number of Reads', fontsize=12)
        axes[0, 0].set_title('5\' TSO (Apoptotic Cell Marker)', fontsize=14, fontweight='bold')
        axes[0, 0].grid(True, alpha=0.3)
        
        axes[0, 1].bar(tso_5p_lengths, tso_5p_counts, alpha=0.7, color=get_adapter_color('^TSO'))
        axes[0, 1].set_yscale('log')
        axes[0, 1].set_xlabel('Length of TSO Match (bp)', fontsize=12)
        axes[0, 1].set_ylabel('Number of Reads (log scale)', fontsize=12)
        axes[0, 1].set_title('5\' TSO Distribution (Log Scale)', fontsize=14, fontweight='bold')
        axes[0, 1].grid(True, alpha=0.3)
    
    # 3' TSO_RC
    if tso_rc_lengths:
        axes[1, 0].bar(tso_rc_lengths, tso_rc_counts, alpha=0.7, color=get_adapter_color('TSO_RC'))
        axes[1, 0].set_xlabel('Length of TSO_RC Match (bp)', fontsize=12)
        axes[1, 0].set_ylabel('Number of Reads', fontsize=12)
        axes[1, 0].set_title('3\' TSO_RC (After Biological Sequence)', fontsize=14, fontweight='bold')
        axes[1, 0].grid(True, alpha=0.3)
        
        bases = ['A', 'C', 'G', 'T']
        base_counts = [tso_rc_adj_bases.get(b, 0) for b in bases]
        colors_bases = ['green', 'blue', 'orange', 'red']
        axes[1, 1].bar(bases, base_counts, alpha=0.7, color=colors_bases)
        axes[1, 1].set_xlabel('Base Before TSO_RC', fontsize=12)
        axes[1, 1].set_ylabel('Number of Reads', fontsize=12)
        axes[1, 1].set_title('Base Composition Before TSO_RC', fontsize=14, fontweight='bold')
        axes[1, 1].grid(True, alpha=0.3)
        
        total = sum(base_counts)
        if total > 0:
            for i, (base, count) in enumerate(zip(bases, base_counts)):
                pct = (count / total) * 100
                axes[1, 1].text(i, count, f'{pct:.1f}%', ha='center', va='bottom', fontsize=10)
    
    plt.suptitle(f'{display_name} - TSO Adapter Detection Patterns (Aggregated)', 
                 fontsize=16, fontweight='bold')
    plt.tight_layout()
    
    safe_filename = library_name.replace('/', '_').replace('\\', '_')
    plt.savefig(output_dir / f'{safe_filename}_tso_positions.png', dpi=300, bbox_inches='tight')
    plt.close()


def plot_library_polyg_analysis(library_name, display_name, json_files_with_runs, output_dir):
    """Plot aggregated poly-G distribution for a library."""
    all_lengths = defaultdict(int)
    adj_bases = defaultdict(int)
    
    for json_file, _ in json_files_with_runs:
        with open(json_file) as f:
            data = json.load(f)
        # Check both read1 and read2 adapters
        adapters = data.get('adapters_read1', []) or data.get('adapters_read2', [])
        for adapter in adapters:
            if adapter['name'] == '1' and adapter.get('three_prime_end'):
        
                trimmed = adapter['three_prime_end'].get('trimmed_lengths', [])
                for item in trimmed:
                    all_lengths[item['len']] += sum(item['counts'])
                
                bases_dict = adapter['three_prime_end'].get('adjacent_bases', {})
                for base, count in bases_dict.items():
                    adj_bases[base] += count
    
    if not all_lengths:
        print(f"  ⚠ No poly-G data for {display_name}")
        return
    
    lengths = sorted(all_lengths.keys())
    counts = [all_lengths[l] for l in lengths]
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    # Full distribution
    axes[0, 0].bar(lengths, counts, alpha=0.7, color=get_adapter_color('Poly-G'), edgecolor='black', linewidth=0.3)
    axes[0, 0].axvline(20, color='red', linestyle='--', linewidth=2, label='Trim threshold (20bp)')
    axes[0, 0].set_xlabel('Poly-G Stretch Length (bp)', fontsize=12)
    axes[0, 0].set_ylabel('Number of Reads', fontsize=12)
    axes[0, 0].set_title('Poly-G Stretch Length Distribution', fontsize=14, fontweight='bold')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)
    
    # Long poly-G
    long_polyg = [(l, c) for l, c in zip(lengths, counts) if l >= 100]
    if long_polyg:
        long_lengths, long_counts = zip(*long_polyg)
        axes[0, 1].bar(long_lengths, long_counts, alpha=0.7, color='#5f5f5f', edgecolor='black', linewidth=0.5)  # Darker gray for long poly-G
        axes[0, 1].set_xlabel('Poly-G Stretch Length (bp)', fontsize=12)
        axes[0, 1].set_ylabel('Number of Reads', fontsize=12)
        axes[0, 1].set_title('Long Poly-G Stretches (>=100bp)', fontsize=14, fontweight='bold')
        axes[0, 1].grid(True, alpha=0.3)
    else:
        axes[0, 1].text(0.5, 0.5, 'No long poly-G stretches detected\n(Good!)', 
                       ha='center', va='center', fontsize=14, transform=axes[0, 1].transAxes)
        axes[0, 1].set_title('Long Poly-G Stretches (>=100bp)', fontsize=14, fontweight='bold')
    
    # Cumulative
    cumulative = np.cumsum(counts)
    cumulative_pct = (cumulative / cumulative[-1]) * 100
    axes[1, 0].plot(lengths, cumulative_pct, linewidth=3, color=get_adapter_color('Poly-G'))
    axes[1, 0].axvline(20, color='red', linestyle='--', linewidth=2, label='Trim threshold')
    axes[1, 0].axhline(90, color='orange', linestyle='--', alpha=0.5)
    axes[1, 0].set_xlabel('Poly-G Stretch Length (bp)', fontsize=12)
    axes[1, 0].set_ylabel('Cumulative % of Reads', fontsize=12)
    axes[1, 0].set_title('Cumulative Distribution', fontsize=14, fontweight='bold')
    axes[1, 0].legend()
    axes[1, 0].grid(True, alpha=0.3)
    axes[1, 0].set_ylim([0, 100])
    
    # Base composition
    bases = ['A', 'C', 'G', 'T']
    base_counts = [adj_bases.get(b, 0) for b in bases]
    colors_bases = ['green', 'blue', 'orange', 'red']
    axes[1, 1].bar(bases, base_counts, alpha=0.7, color=colors_bases, edgecolor='black', linewidth=1)
    axes[1, 1].set_xlabel('Base Before Poly-G', fontsize=12)
    axes[1, 1].set_ylabel('Number of Reads', fontsize=12)
    axes[1, 1].set_title('Base Composition Before Poly-G', fontsize=14, fontweight='bold')
    axes[1, 1].grid(True, alpha=0.3)
    
    total = sum(base_counts)
    if total > 0:
        for i, (base, count) in enumerate(zip(bases, base_counts)):
            pct = (count / total) * 100
            axes[1, 1].text(i, count, f'{pct:.1f}%', ha='center', va='bottom', fontsize=11, fontweight='bold')
    
    plt.suptitle(f'{display_name} - Poly-G Contamination Analysis (Aggregated)', 
                 fontsize=16, fontweight='bold')
    plt.tight_layout()
    
    safe_filename = library_name.replace('/', '_').replace('\\', '_')
    plt.savefig(output_dir / f'{safe_filename}_polyg.png', dpi=300, bbox_inches='tight')
    plt.close()


def separate_r1_r2_files(json_files):
    """
    Separate JSON files into R1 and R2 groups.
    
    Returns:
        r1_files: List of R1 JSON files
        r2_files: List of R2 JSON files (includes both paired-end R2 and long R2)
        has_r1: Boolean indicating if any R1 files exist
    """
    r1_files = []
    r2_files = []
    
    for json_file in json_files:
        sample = Path(json_file).stem.replace('_cutadapt', '')
        
        if '_R1_' in sample or sample.endswith('_R1'):
            r1_files.append(json_file)
        elif '_R2_' in sample or sample.endswith('_R2'):
            r2_files.append(json_file)
        else:
            # If no R1 or R2 in name, assume it's R2 (shouldn't happen with new pipeline)
            r2_files.append(json_file)
    
    return r1_files, r2_files, len(r1_files) > 0


def plot_global_adapter_positions(json_files, output_file):
    """
    Create global aggregate visualization of where adapters are found across all files.
    Shows 5' vs 3' position distribution and typical match lengths.
    """
    print("\nAnalyzing adapter positions across all files...")
    
    # Aggregate data across all files
    adapter_data = defaultdict(lambda: {
        'five_prime_matches': 0,
        'three_prime_matches': 0,
        'five_prime_lengths': defaultdict(int),
        'three_prime_lengths': defaultdict(int),
        'total_matches': 0
    })
    
    total_reads = 0
    total_bp = 0
    
    for json_file in json_files:
        with open(json_file) as f:
            data = json.load(f)
        
        # Accumulate total reads
        total_reads += data['read_counts']['input']
        total_bp += data['basepair_counts']['input_read1']
        
        # Process each adapter
        for adapter in data.get('adapters_read1', []):
            name = adapter['name']
            
            # 5' end matches
            five_prime = adapter.get('five_prime_end', {})
            five_matches = five_prime.get('matches', 0)
            adapter_data[name]['five_prime_matches'] += five_matches
            
            # Aggregate 5' trimmed lengths
            for item in five_prime.get('trimmed_lengths', []):
                length = item['len']
                counts = item['counts']
                total_count = sum(counts) if isinstance(counts, list) else counts
                adapter_data[name]['five_prime_lengths'][length] += total_count
            
            # 3' end matches
            three_prime = adapter.get('three_prime_end', {})
            three_matches = three_prime.get('matches', 0)
            adapter_data[name]['three_prime_matches'] += three_matches
            
            # Aggregate 3' trimmed lengths
            for item in three_prime.get('trimmed_lengths', []):
                length = item['len']
                counts = item['counts']
                total_count = sum(counts) if isinstance(counts, list) else counts
                adapter_data[name]['three_prime_lengths'][length] += total_count
            
            adapter_data[name]['total_matches'] += adapter.get('total_matches', 0)
    
    # Calculate average read length
    avg_read_length = total_bp / total_reads if total_reads > 0 else 0
    
    # Filter to significant adapters (>0.01% of reads) and sort
    min_threshold = total_reads * 0.0001  # 0.01% threshold
    significant = [(name, data) for name, data in adapter_data.items() 
                   if data['total_matches'] > min_threshold]
    significant.sort(key=lambda x: x[1]['total_matches'], reverse=True)
    
    if not significant:
        print("  ⚠ No adapters with >0.01% of reads found")
        return
    
    # Take top 15 for visualization
    top_adapters = significant[:15]
    
    # Create figure with 3 subplots
    fig = plt.figure(figsize=(16, 14))
    gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 1.3], hspace=0.35, wspace=0.3)
    
    # Top: 5' vs 3' distribution
    ax1 = fig.add_subplot(gs[0, :])
    
    adapter_names = [name for name, _ in top_adapters]
    five_counts = [data['five_prime_matches'] for _, data in top_adapters]
    three_counts = [data['three_prime_matches'] for _, data in top_adapters]
    
    x = np.arange(len(adapter_names))
    width = 0.35
    
    colors = [get_adapter_color(name, i) for i, name in enumerate(adapter_names)]
    
    bars1 = ax1.barh(x - width/2, five_counts, width, label="5' end", 
                     color=colors, alpha=0.7, edgecolor='black')
    bars2 = ax1.barh(x + width/2, three_counts, width, label="3' end",
                     color=colors, alpha=0.4, edgecolor='black', hatch='//')
    
    ax1.set_xlabel('Total Matches Across All Files', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Adapter', fontsize=12, fontweight='bold')
    ax1.set_title("Global Adapter Position Distribution (5' vs 3' end of read)", 
                  fontsize=14, fontweight='bold')
    ax1.set_yticks(x)
    ax1.set_yticklabels(adapter_names, fontsize=10)
    ax1.legend(fontsize=11, loc='upper left', bbox_to_anchor=(1.01, 1))
    ax1.grid(True, alpha=0.3, axis='x')
    
    # Add percentage labels - only if >=0.01%
    for i, (five, three) in enumerate(zip(five_counts, three_counts)):
        five_pct = (five / total_reads) * 100
        three_pct = (three / total_reads) * 100
        if five_pct >= 0.01:  # Only show if >=0.01%
            ax1.text(five, i - width/2, f'{five_pct:.1f}%', 
                    va='center', ha='left', fontsize=8)
        if three_pct >= 0.01:  # Only show if >=0.01%
            ax1.text(three, i + width/2, f'{three_pct:.1f}%',
                    va='center', ha='left', fontsize=8)
    
    # Middle left: Top 5' adapters match length distribution
    ax2 = fig.add_subplot(gs[1, 0])
    
    # Find top 5 adapters by 5' matches
    top_five_prime = sorted(top_adapters, key=lambda x: x[1]['five_prime_matches'], reverse=True)[:5]
    
    plotted_any = False
    for name, data in top_five_prime:
        if data['five_prime_matches'] < 100:
            continue
        lengths = sorted(data['five_prime_lengths'].keys())
        counts = [data['five_prime_lengths'][l] for l in lengths]
        if lengths:
            color = get_adapter_color(name)
            ax2.plot(lengths, counts, marker='o', linewidth=2, label=name, 
                    color=color, alpha=0.7, markersize=4)
            plotted_any = True
    
    if plotted_any:
        ax2.set_xlabel('Adapter Match Length (bp)', fontsize=11)
        ax2.set_ylabel('Number of Reads', fontsize=11)
        ax2.set_title("5' End Match Length Distribution (Top 5 Adapters)", 
                      fontsize=12, fontweight='bold')
        ax2.legend(fontsize=9, loc='upper left', bbox_to_anchor=(1.01, 1))
        ax2.grid(True, alpha=0.3)
        ax2.set_yscale('log')
    else:
        ax2.text(0.5, 0.5, 'No significant 5\' matches', ha='center', va='center',
                transform=ax2.transAxes, fontsize=12)
        ax2.set_title("5' End Match Length Distribution", fontsize=12, fontweight='bold')
    
    # Middle right: Top 3' adapters match length distribution
    ax3 = fig.add_subplot(gs[1, 1])
    
    # Find top 5 adapters by 3' matches
    top_three_prime = sorted(top_adapters, key=lambda x: x[1]['three_prime_matches'], reverse=True)[:5]
    
    plotted_any = False
    for name, data in top_three_prime:
        if data['three_prime_matches'] < 100:
            continue
        lengths = sorted(data['three_prime_lengths'].keys())
        counts = [data['three_prime_lengths'][l] for l in lengths]
        if lengths:
            color = get_adapter_color(name)
            ax3.plot(lengths, counts, marker='o', linewidth=2, label=name,
                    color=color, alpha=0.7, markersize=4)
            plotted_any = True
    
    if plotted_any:
        ax3.set_xlabel('Adapter Match Length (bp)', fontsize=11)
        ax3.set_ylabel('Number of Reads', fontsize=11)
        ax3.set_title("3' End Match Length Distribution (Top 5 Adapters)", 
                      fontsize=12, fontweight='bold')
        ax3.legend(fontsize=9, loc='upper left', bbox_to_anchor=(1.01, 1))
        ax3.grid(True, alpha=0.3)
        ax3.set_yscale('log')
    else:
        ax3.text(0.5, 0.5, 'No significant 3\' matches', ha='center', va='center',
                transform=ax3.transAxes, fontsize=12)
        ax3.set_title("3' End Match Length Distribution", fontsize=12, fontweight='bold')
    
    # Bottom: Position schematic
    ax4 = fig.add_subplot(gs[2, :])
    
    read_length = int(avg_read_length)
    
    # Draw read as a line
    ax4.plot([0, read_length], [0, 0], 'k-', linewidth=10, alpha=0.3, label='Read', zorder=1)
    
    # Filter adapters for schematic - only show those with significant matches
    schematic_adapters = [(name, data) for name, data in top_adapters 
                          if data['five_prime_matches'] >= 100 or data['three_prime_matches'] >= 100]
    schematic_adapters = schematic_adapters[:10]  # Top 10 only
    
    y_offset = 0
    y_spacing = 1.2  # Increased spacing to reduce overlap
    
    for i, (name, data) in enumerate(schematic_adapters):
        five_matches = data['five_prime_matches']
        three_matches = data['three_prime_matches']
        
        color = get_adapter_color(name, i)
        y_pos = y_offset
        
        # Calculate average match lengths for positioning
        if five_matches > 0 and data['five_prime_lengths']:
            five_lengths = data['five_prime_lengths']
            avg_five_len = sum(l * c for l, c in five_lengths.items()) / sum(five_lengths.values())
        else:
            avg_five_len = 0
            
        if three_matches > 0 and data['three_prime_lengths']:
            three_lengths = data['three_prime_lengths']
            avg_three_len = sum(l * c for l, c in three_lengths.items()) / sum(three_lengths.values())
        else:
            avg_three_len = 0
        
        # Plot position indicators
        if five_matches > three_matches * 2:  # Predominantly 5'
            # Draw at 5' end with size proportional to avg match length
            match_len = min(avg_five_len, read_length * 0.3)  # Cap at 30% of read
            ax4.scatter([match_len/2], [y_pos], s=600, c=[color], 
                       alpha=0.7, marker='>', edgecolors='black', linewidths=2, zorder=3)
            ax4.plot([0, match_len], [y_pos, y_pos], color=color, 
                    linewidth=5, alpha=0.3, zorder=2)
        elif three_matches > five_matches * 2:  # Predominantly 3'
            # Draw at 3' end
            match_len = min(avg_three_len, read_length * 0.3)  # Cap at 30% of read
            end_pos = read_length - match_len/2
            ax4.scatter([end_pos], [y_pos], s=600, c=[color],
                       alpha=0.7, marker='<', edgecolors='black', linewidths=2, zorder=3)
            ax4.plot([read_length - match_len, read_length], [y_pos, y_pos], 
                    color=color, linewidth=5, alpha=0.3, zorder=2)
        else:  # Both ends
            match_len_5 = min(avg_five_len, read_length * 0.15)
            match_len_3 = min(avg_three_len, read_length * 0.15)
            ax4.scatter([match_len_5/2, read_length - match_len_3/2], [y_pos, y_pos], 
                       s=600, c=[color, color], alpha=0.7, 
                       marker='o', edgecolors='black', linewidths=2, zorder=3)
        
        # Add label - with background box to avoid overlap
        ax4.text(-read_length * 0.1, y_pos, name, va='center', ha='right', 
                fontsize=9, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', 
                         edgecolor='none', alpha=0.8))
        
        y_offset += y_spacing
    
    ax4.set_xlim(-read_length * 0.15, read_length * 1.05)
    ax4.set_ylim(-1, y_offset + 0.5)
    ax4.set_xlabel(f'Position in Read (avg length: {read_length:.0f}bp)', 
                   fontsize=12, fontweight='bold')
    ax4.set_title('Typical Adapter Positions Along Read\n(Symbol position and line length indicate average match length)', 
                  fontsize=12, fontweight='bold', pad=15)
    ax4.set_yticks([])
    ax4.spines['left'].set_visible(False)
    ax4.spines['right'].set_visible(False)
    ax4.spines['top'].set_visible(False)
    ax4.grid(True, alpha=0.2, axis='x')
    
    # Add markers for read position
    ax4.axvline(0, color='green', linestyle='--', alpha=0.5, linewidth=2, 
                label='5\' end', zorder=1)
    ax4.axvline(read_length, color='red', linestyle='--', alpha=0.5, linewidth=2, 
                label='3\' end', zorder=1)
    ax4.legend(loc='upper right', fontsize=10, framealpha=0.9)
    
    plt.suptitle(f'Global Adapter Position Analysis\n({len(json_files)} files, {total_reads:,} total reads)', 
                 fontsize=16, fontweight='bold', y=0.995)
    
    plt.tight_layout(rect=[0, 0, 0.98, 0.99])  # Leave room for legends on right
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  ✓ Global adapter position analysis saved to {output_file}")


def main():
    """Generate all visualizations."""
    
    # Parse command line arguments
    args = parse_arguments()
    
    # Create output directories
    output_dir = Path('cutadapt_figures')
    
    print("\n" + "="*70)
    print("CUTADAPT COMPREHENSIVE ANALYSIS SUITE")
    print("="*70)
    
    # Find JSON files
    print(f"\nSearching for JSON files in: {args.path}/")
    search_path = Path(args.path)
    
    if not search_path.exists():
        print(f"\nERROR: Directory not found: {search_path}")
        sys.exit(1)
    
    json_files = sorted(glob.glob(str(search_path / '*_cutadapt.json')))
    
    if not json_files:
        print(f"\nERROR: No *_cutadapt.json files found in {search_path}/")
        print("\nPlease ensure your cutadapt JSON files are in the specified directory.")
        sys.exit(1)
    
    print(f"Found {len(json_files)} cutadapt JSON reports")
    
    # Separate R1 and R2 files
    r1_files, r2_files, has_r1 = separate_r1_r2_files(json_files)
    
    print(f"\nFile breakdown:")
    if r1_files:
        print(f"  - R1 files: {len(r1_files)}")
    if r2_files:
        print(f"  - R2 files: {len(r2_files)}")
    
    # Process R1 files if they exist
    if r1_files:
        print("\n" + "="*70)
        print("PROCESSING R1 FILES")
        print("="*70)
        
        global_dir_r1 = output_dir / 'global_R1'
        library_dir_r1 = output_dir / 'by_library_R1'
        
        for dir_path in [output_dir, global_dir_r1, library_dir_r1]:
            dir_path.mkdir(exist_ok=True)
        
        process_file_group(r1_files, global_dir_r1, library_dir_r1, "R1", args)
    
    # Process R2 files if they exist
    if r2_files:
        print("\n" + "="*70)
        print("PROCESSING R2 FILES")
        print("="*70)
        
        global_dir_r2 = output_dir / 'global_R2'
        library_dir_r2 = output_dir / 'by_library_R2'
        
        for dir_path in [output_dir, global_dir_r2, library_dir_r2]:
            dir_path.mkdir(exist_ok=True)
        
        process_file_group(r2_files, global_dir_r2, library_dir_r2, "R2", args)
    
    # Final summary
    print("\n" + "="*70)
    print("ANALYSIS COMPLETE!")
    print("="*70)
    print(f"\nAll results saved to {output_dir}/")
    
    if r1_files:
        print(f"\n📁 R1 Analysis:")
        print(f"   - Global plots: {output_dir / 'global_R1'}/")
        print(f"   - Library plots: {output_dir / 'by_library_R1'}/")
    
    if r2_files:
        print(f"\n📁 R2 Analysis:")
        print(f"   - Global plots: {output_dir / 'global_R2'}/")
        print(f"   - Library plots: {output_dir / 'by_library_R2'}/")
    
    print("\n" + "="*70 + "\n")


def process_file_group(json_files, global_dir, library_dir, read_type, args):
    """Process a group of files (either R1 or R2)."""
    
    # Group samples by library
    library_groups = group_samples_by_library(json_files)
    
    # Create display names for libraries
    library_names = list(library_groups.keys())
    display_names = create_display_names(library_names, args.full_names)
    
    print(f"\nIdentified {len(library_groups)} unique {read_type} libraries:")
    for library, runs in library_groups.items():
        display = display_names[library]
        if display != library:
            print(f"  • {display} ({library}): {len(runs)} sequencing runs")
        else:
            print(f"  • {library}: {len(runs)} sequencing runs")
    
    if not args.full_names and len(set(display_names.values())) < len(display_names):
        print("\nNote: Some library names may not be unique after truncation.")
        print("      Use --full-names flag to use full library names.")
    
    # === GLOBAL ANALYSIS ===
    print(f"\n{'='*70}")
    print(f"PHASE 1: {read_type} GLOBAL ANALYSIS")
    print("="*70)
    
    # Calculate expected values
    print("\nCalculating global statistics...")
    expected = calculate_expected_values(json_files)
    print(f"  • Median TSO (5'): {expected['TSO_5p']:.2f}%")
    print(f"  • Median TSO_RC (3'): {expected['TSO_RC']:.3f}%")
    print(f"  • Median Poly-G: {expected['Poly_G']:.2f}%")
    print(f"  • Median reads passing: {expected['Passing']:.2f}%")
    
    # Create summary tables
    print("\nCreating summary tables...")
    combined_df, pivot_df = create_adapter_summary(json_files)
    
    if not combined_df.empty:
        combined_df.to_csv(f'adapter_detailed_{read_type}.csv', index=False)
        print(f"  ✓ adapter_detailed_{read_type}.csv")
    
    if not pivot_df.empty:
        pivot_df.to_csv(f'adapter_by_sample_{read_type}.csv')
        print(f"  ✓ adapter_by_sample_{read_type}.csv")
        
        summary = create_summary_table(pivot_df, f'adapter_summary_{read_type}.csv')
    
    # Generate global heatmaps
    print("\nGenerating global heatmaps...")
    plot_adapter_heatmap(pivot_df, global_dir / f'adapter_heatmap_all_runs_{read_type}.png')
    plot_library_aggregated_heatmap(library_groups, display_names, 
                                   global_dir / f'adapter_heatmap_by_library_{read_type}.png',
                                   use_clustering=(not args.no_clustering))
    
    # Generate library-aggregated comparison plots
    print("\nGenerating library-aggregated comparison plots...")
    plot_library_other_adapters_breakdown(library_groups, display_names, 
                                          global_dir / f'library_other_adapters_breakdown_{read_type}.png')
    plot_library_multi_comparison(library_groups, display_names, 
                                  global_dir / f'library_multi_comparison_{read_type}.png')
    
    # Generate global adapter position analysis
    print("\nGenerating global adapter position analysis...")
    plot_global_adapter_positions(json_files, 
                                  global_dir / f'adapter_positions_global_{read_type}.png')
    
    # === LIBRARY-LEVEL ANALYSIS ===
    print(f"\n{'='*70}")
    print(f"PHASE 2: {read_type} LIBRARY-LEVEL ANALYSIS (Aggregated)")
    print("="*70)
    
    for i, (library, json_files_with_runs) in enumerate(library_groups.items(), 1):
        display_name = display_names[library]
        print(f"\n[{i}/{len(library_groups)}] Processing library: {display_name}")
        print(f"  Runs: {len(json_files_with_runs)}")
        
        # Aggregate metrics
        agg_metrics = aggregate_library_metrics(json_files_with_runs)
        
        # Create all plots for this library
        print("  → Creating aggregated summary panel...", end=" ")
        plot_library_aggregated_summary(library, display_name, agg_metrics, json_files_with_runs, library_dir)
        print("✓")
        
        print("  → Creating read length distribution...", end=" ")
        plot_library_read_length_distribution(library, display_name, json_files_with_runs, library_dir)
        print("✓")
        
        print("  → Creating TSO position analysis...", end=" ")
        plot_library_tso_analysis(library, display_name, json_files_with_runs, library_dir)
        print("✓")
        
        print("  → Creating poly-G analysis...", end=" ")
        plot_library_polyg_analysis(library, display_name, json_files_with_runs, library_dir)
        print("✓")


if __name__ == '__main__':
    main()
