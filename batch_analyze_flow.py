#!/usr/bin/env python3
"""
Batch Flow Cytometry Analysis Script

Automatically generates flow cytometry plots for multiple keyword filter values
without manual intervention. This script wraps the FlowAnalyzer class to enable
batch processing with custom configurations.
Repository: https://github.com/cyeunghl/flowViz

Author: Generated with Claude Code
License: MIT
Repository: https://github.com/yourusername/flowViz

Usage:
    python batch_analyze_flow.py

Configuration:
    Edit the CONFIGURATION section below to customize for your experiment:
    - wsp_path: Path to your FlowJo workspace (.wsp) file
    - fcs_dir: Directory containing your FCS files
    - output_dir: Directory where HTML output files will be saved
    - filter_values: List of keyword filter values to iterate over
    - gate_path: Gate path tuple (e.g., ('root',) for root gates)
    - gate_name: Name of the gate to visualize
    - channels: X and Y axis channel names for scatter plots
    - keyword_filter_key: Keyword to filter samples by
    - well_id_keyword: Keyword containing well ID information

Output:
    Generates one HTML file per filter value in the output directory.
    Each file contains a 96-well plate layout with plots and gate overlays.

Example Configuration:
    For time series analysis with time points [48, 72, 96, 144]:
    - filter_values = ['48', '72', '96', '144']
    - keyword_filter_key = 'Time point (hr)'

    This generates 4 HTML files:
    - workspace_scatter_GateName_FSC-A_SSC-A_Time_point_hr_48.html
    - workspace_scatter_GateName_FSC-A_SSC-A_Time_point_hr_72.html
    - workspace_scatter_GateName_FSC-A_SSC-A_Time_point_hr_96.html
    - workspace_scatter_GateName_FSC-A_SSC-A_Time_point_hr_144.html
"""

import sys
import os
from pathlib import Path

# Add current directory to path for importing analyze_flow
sys.path.insert(0, str(Path(__file__).parent))

from analyze_flow import FlowAnalyzer

def main():
    """Main batch processing function."""

    # ========================================================================
    # CONFIGURATION - Edit these values for your experiment
    # ========================================================================

    # File paths - UPDATE THESE for your system
    wsp_path = "path/to/your/workspace.wsp"  # e.g., "/data/experiment/workspace.wsp"
    fcs_dir = "path/to/your/fcs_files"       # e.g., "/data/experiment/fcs"
    output_dir = "path/to/output"            # e.g., "/data/experiment/output"

    # Analysis parameters - UPDATE THESE for your experiment
    filter_values = ['48', '72', '96', '144']  # Values to iterate over (e.g., time points)
    gate_path = ('root',)                      # Gate path as tuple, e.g., ('root',) or ('root', 'Cells', 'Singlets')
    gate_name = 'Cells_withDebris'             # Gate to filter data and visualize
    x_channel = 'FSC-A'                        # X-axis channel for scatter plot
    y_channel = 'SSC-A'                        # Y-axis channel for scatter plot
    keyword_filter_key = 'Time point (hr)'     # Keyword to filter by (e.g., 'Time point (hr)', 'Treatment', 'Plate')
    well_id_keyword = '$WELLID'                # Keyword containing well ID (e.g., '$WELLID', 'Well ID', 'WellID')

    # ========================================================================
    # END CONFIGURATION
    # ========================================================================

    print(f"\n{'='*70}")
    print("Batch Flow Cytometry Analysis")
    print(f"{'='*70}\n")
    print(f"Workspace: {wsp_path}")
    print(f"FCS Directory: {fcs_dir}")
    print(f"Output Directory: {output_dir}")
    print(f"Filter values: {', '.join(filter_values)}\n")

    # Validate paths
    if not os.path.exists(wsp_path):
        print(f"✗ Error: Workspace file not found: {wsp_path}")
        print("Please update the 'wsp_path' variable in the CONFIGURATION section.")
        return 1

    if not os.path.exists(fcs_dir):
        print(f"✗ Error: FCS directory not found: {fcs_dir}")
        print("Please update the 'fcs_dir' variable in the CONFIGURATION section.")
        return 1

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Initialize FlowAnalyzer
    print("Loading workspace...")
    try:
        analyzer = FlowAnalyzer(wsp_path, fcs_dir)
        print("✓ Workspace loaded successfully\n")
    except Exception as e:
        print(f"✗ Failed to load workspace: {e}")
        return 1

    # Get all sample IDs
    try:
        all_sample_ids = analyzer.workspace.get_sample_ids(loaded_only=True)
        print(f"Found {len(all_sample_ids)} total samples in workspace\n")
    except Exception as e:
        print(f"✗ Failed to get sample IDs: {e}")
        return 1

    # Process each filter value
    successful_outputs = []
    failed_outputs = []

    for filter_value in filter_values:
        print(f"\n{'─'*70}")
        print(f"Processing {keyword_filter_key}: {filter_value}")
        print(f"{'─'*70}\n")

        try:
            # Get workspace basename for output filename
            base_name = os.path.splitext(os.path.basename(wsp_path))[0]

            # Build output filename following analyze_flow.py naming convention
            gate_name_clean = gate_name.replace(' ', '_').replace('/', '_')
            params_str = f"{x_channel}_{y_channel}".replace(' ', '_').replace('/', '_')
            key_clean = keyword_filter_key.replace(' ', '_').replace('/', '_').replace('(', '').replace(')', '')
            value_clean = str(filter_value).replace(' ', '_').replace('/', '_')

            output_filename = f"{base_name}_scatter_{gate_name_clean}_{params_str}_{key_clean}_{value_clean}.html"
            output_path = os.path.join(output_dir, output_filename)

            # Construct selections dictionary (matches structure from interactive_plot_prompt)
            selections = {
                # Sample selection - use all samples, filtering happens via keyword_filter
                'selected_groups': ['All Samples'],

                # Display options - no additional information
                'show_keywords': False,
                'show_statistics': False,
                'selected_keywords_to_show': [],

                # Keyword filter - THIS is what filters samples by the specified value
                'keyword_filter': {
                    'key': keyword_filter_key,
                    'value': filter_value
                },

                # Well ID extraction - use specified keyword
                'well_id_source': 'keyword',
                'well_id_keyword': well_id_keyword,

                # Plot configuration
                'num_plots': 1,
                'plot_configs': [
                    {
                        'gate_path': gate_path,  # Gate path as tuple
                        'gate_name': gate_name,
                        'plot_type': 'scatter',
                        'parameters': [x_channel, y_channel],  # X and Y channels
                        'show_gates': True,
                        'gates_to_visualize': [gate_name]
                    }
                ]
            }

            # Generate plots
            print(f"Generating scatter plots for {keyword_filter_key} = {filter_value}...")
            analyzer.generate_interactive_plots(selections, output_path)

            if os.path.exists(output_path):
                successful_outputs.append(output_filename)
                print(f"✓ Generated: {output_filename}\n")
            else:
                failed_outputs.append((filter_value, "Output file not found"))
                print(f"⚠ Output file not found: {output_filename}\n")

        except Exception as e:
            failed_outputs.append((filter_value, str(e)))
            print(f"✗ Failed to process {keyword_filter_key} = {filter_value}: {e}\n")
            import traceback
            traceback.print_exc()

    # Summary
    print(f"\n{'='*70}")
    print("Batch Processing Complete")
    print(f"{'='*70}\n")
    print(f"Successful: {len(successful_outputs)}/{len(filter_values)}")

    if successful_outputs:
        print("\nGenerated files:")
        for filename in successful_outputs:
            print(f"  ✓ {filename}")

    if failed_outputs:
        print(f"\nFailed: {len(failed_outputs)}/{len(filter_values)}")
        for filter_val, error in failed_outputs:
            print(f"  ✗ {keyword_filter_key} = {filter_val}: {error}")
        return 1

    print("\n✓ All filter values processed successfully!\n")
    return 0

if __name__ == "__main__":
    sys.exit(main())
