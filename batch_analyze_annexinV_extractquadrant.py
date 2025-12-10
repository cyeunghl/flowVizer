#!/usr/bin/env python3
"""
Batch Flow Cytometry Quadrant Extraction for AnnexinV Data

Automatically generates Singlets scatter plots (B1-A vs R2-A) for multiple time points without
manual intervention.

Usage:
    python batch_analyze_annexinV_extractquadrant.py

Output:
    Generates 4 HTML files in /Users/clarenceyeung/Downloads/annexinV/:
    - 20251209_annexinV_v2_scatter_Singlets_B1-A_R2-A_Time_point_hr_48.html
    - 20251209_annexinV_v2_scatter_Singlets_B1-A_R2-A_Time_point_hr_72.html
    - 20251209_annexinV_v2_scatter_Singlets_B1-A_R2-A_Time_point_hr_96.html
    - 20251209_annexinV_v2_scatter_Singlets_B1-A_R2-A_Time_point_hr_144.html
"""

import sys
import os
from pathlib import Path

# Add current directory to path for importing analyze_flow
sys.path.insert(0, str(Path(__file__).parent))

from analyze_flow import FlowAnalyzer

# Gate/scatter configuration derived from interactive prompt
SCATTER_PARAMETERS = ['B1-A', 'R2-A']
GATE_NAME = 'Singlets'
GATE_PATH = ('root', 'Cells_withDebris', 'Singlets')
WELL_ID_KEYWORD = '$WELLID'


def build_output_filename(wsp_path, selections):
    """Construct an output filename matching analyze_flow interactive defaults."""
    base_name = os.path.splitext(os.path.basename(wsp_path))[0]

    plot_configs = selections.get('plot_configs', [])
    if not plot_configs:
        plot_configs = [{
            'plot_type': selections.get('plot_type', 'histogram'),
            'gate_name': selections.get('gate_name', 'Ungated'),
            'parameters': selections.get('parameters', [])
        }]

    if len(plot_configs) == 1:
        plot_type = plot_configs[0].get('plot_type', 'plot')
        gate_name = plot_configs[0].get('gate_name', 'gate').replace(' ', '_').replace('/', '_')
        params = plot_configs[0].get('parameters', [])
        params_str = '_'.join(params).replace(' ', '_').replace('/', '_')
    else:
        plot_type = 'multi'
        gate_name = 'multi'
        params_str = f"{len(plot_configs)}plots"

    filter_str = ""
    keyword_filter = selections.get('keyword_filter')
    if keyword_filter:
        key_clean = keyword_filter['key'].replace(' ', '_').replace('/', '_').replace('(', '').replace(')', '')
        value_clean = str(keyword_filter['value']).replace(' ', '_').replace('/', '_')
        filter_str = f"_{key_clean}_{value_clean}"

    return f"{base_name}_{plot_type}_{gate_name}_{params_str}{filter_str}.html"

def main():
    """Main batch processing function."""

    # Configuration
    wsp_path = "/Users/clarenceyeung/Downloads/annexinV/20251209_annexinV_v2.wsp"
    fcs_dir = "/Users/clarenceyeung/Downloads/annexinV/forAnalysis"
    output_dir = "/Users/clarenceyeung/Downloads/annexinV"
    os.makedirs(output_dir, exist_ok=True)

    # Time points to iterate over
    time_points = ['48', '72', '96', '144']

    print(f"\n{'='*70}")
    print("Batch Flow Cytometry Analysis - AnnexinV Time Series")
    print(f"{'='*70}\n")
    print(f"Workspace: {wsp_path}")
    print(f"FCS Directory: {fcs_dir}")
    print(f"Output Directory: {output_dir}")
    print(f"Time Points: {', '.join(time_points)} hr\n")

    # Initialize FlowAnalyzer
    print("Loading workspace...")
    try:
        analyzer = FlowAnalyzer(wsp_path, fcs_dir)
        print("✓ Workspace loaded successfully\n")
    except Exception as e:
        print(f"✗ Failed to load workspace: {e}")
        return 1

    # Get all sample IDs for "All Samples" group
    try:
        all_sample_ids = analyzer.workspace.get_sample_ids(loaded_only=True)
        print(f"Found {len(all_sample_ids)} total samples in workspace\n")
    except Exception as e:
        print(f"✗ Failed to get sample IDs: {e}")
        return 1

    if not all_sample_ids:
        print("✗ No samples available in the workspace")
        return 1

    reference_sample_id = all_sample_ids[0]
    gates_to_visualize = []
    try:
        gate_polygons = analyzer._extract_gate_polygons(reference_sample_id, SCATTER_PARAMETERS[0], SCATTER_PARAMETERS[1])
        gates_to_visualize = [gate['name'] for gate in gate_polygons]
        if gates_to_visualize:
            print(f"Gates available for visualization on {SCATTER_PARAMETERS[0]} vs {SCATTER_PARAMETERS[1]}: {', '.join(gates_to_visualize)}\n")
        else:
            print(f"No gates available for visualization on {SCATTER_PARAMETERS[0]} vs {SCATTER_PARAMETERS[1]}\n")
    except Exception as e:
        print(f"Warning: Unable to extract gates for visualization ({e})\n")
        gates_to_visualize = []

    # Process each time point
    successful_outputs = []
    failed_outputs = []

    for time_point in time_points:
        print(f"\n{'─'*70}")
        print(f"Processing Time Point: {time_point} hr")
        print(f"{'─'*70}\n")

        try:
            # Filter samples by time point
            # Note: The keyword filter will be applied during plot generation
            # We just need to construct the selections dictionary

            # Construct selections dictionary (matches structure from interactive_plot_prompt)
            selections = {
                # Sample selection - use all samples, filtering happens via keyword_filter
                'selected_groups': ['All Samples'],

                # Display options - no additional information
                'show_keywords': False,
                'show_statistics': False,
                'selected_keywords_to_show': [],

                # Keyword filter - THIS is what filters samples by time point
                'keyword_filter': {
                    'key': 'Time point (hr)',
                    'value': time_point
                },

                # Well ID extraction - use $WELLID keyword
                'well_id_source': 'keyword',
                'well_id_keyword': WELL_ID_KEYWORD,

                # Plot configuration
                'num_plots': 1,
                'plot_configs': [
                    {
                        'gate_path': GATE_PATH,
                        'gate_name': GATE_NAME,
                        'plot_type': 'scatter',
                        'parameters': SCATTER_PARAMETERS,
                        'show_gates': bool(gates_to_visualize),
                        'gates_to_visualize': gates_to_visualize
                    }
                ]
            }

            # Generate deterministic output path matching analyze_flow defaults
            output_filename = build_output_filename(wsp_path, selections)
            output_path = os.path.join(output_dir, output_filename)

            # Generate plots
            print(f"Generating scatter plots for time point {time_point} hr...")
            analyzer.generate_interactive_plots(selections, output_path)

            # The actual filename is auto-generated, so we reconstruct it for tracking
            expected_filename = output_filename
            expected_path = output_path

            if os.path.exists(expected_path):
                successful_outputs.append(expected_filename)
                print(f"✓ Generated: {expected_filename}\n")
            else:
                # File might have been generated with slightly different name
                # Look for any file matching the pattern
                import glob
                pattern = os.path.join(output_dir, f"*{time_point}.html")
                matches = glob.glob(pattern)
                if matches:
                    actual_filename = os.path.basename(matches[-1])  # Get most recent
                    successful_outputs.append(actual_filename)
                    print(f"✓ Generated: {actual_filename}\n")
                else:
                    failed_outputs.append((time_point, "Output file not found"))
                    print(f"⚠ Output file not found for time point {time_point}\n")

        except Exception as e:
            failed_outputs.append((time_point, str(e)))
            print(f"✗ Failed to process time point {time_point}: {e}\n")
            import traceback
            traceback.print_exc()

    # Summary
    print(f"\n{'='*70}")
    print("Batch Processing Complete")
    print(f"{'='*70}\n")
    print(f"Successful: {len(successful_outputs)}/{len(time_points)}")

    if successful_outputs:
        print("\nGenerated files:")
        for filename in successful_outputs:
            print(f"  ✓ {filename}")

    if failed_outputs:
        print(f"\nFailed: {len(failed_outputs)}/{len(time_points)}")
        for time_point, error in failed_outputs:
            print(f"  ✗ Time point {time_point}: {error}")
        return 1

    print("\n✓ All time points processed successfully!\n")
    return 0

if __name__ == "__main__":
    sys.exit(main())
