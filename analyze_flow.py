import argparse
import os
import sys
import flowkit as fk
import pandas as pd
import numpy as np
from scipy.stats import gaussian_kde
from bokeh.plotting import figure, save, output_file
from bokeh.layouts import column, row, gridplot
from bokeh.models import ColumnDataSource, HoverTool, Div
from bokeh.io import curdoc
try:
    from matplotlib.path import Path as MPLPath
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MPLPath = None
    MATPLOTLIB_AVAILABLE = False

"""
FlowJo Analysis Tool - Interactive Plotting Mode
================================================

This script analyzes FlowJo workspaces (.wsp) and associated FCS files to generate
interactive HTML plots. It leverages `flowkit` for parsing flow cytometry data and
`bokeh` for creating interactive visualizations.

Features:
---------
- Batch generation of histograms or scatter plots across all samples
- Optional keyword filtering to subset samples (e.g., filter by time point)
- Flexible well ID extraction (from keywords, filename, or auto-detect)
- Histogram options: linear/log scale, median/mean statistic display
- Organizes plots in a 96-well plate grid layout
- Styled with minimalist Tailwind CSS
- All visualizations use raw (untransformed) data values

Usage:
------
    # Interactive plotting mode (required)
    python analyze_flow.py --wsp <path_to_workspace> --fcs_dir <path_to_fcs_files> --interactive [--output <output_html>]

Interactive Mode Workflow:
--------------------------
The interactive mode guides you through the following steps:

1. Keyword Filtering (Optional):
   - Option to filter samples by a keyword (e.g., "Time point (hr)")
   - Lists all available keywords and their values
   - Select keyword and value to filter samples before plotting
   - Useful when workspace contains multiple files per well (e.g., different time points)

2. Well ID Extraction:
   - Choose where to extract well ID information:
     * Auto: Try keywords first, then filename
     * From keyword: Select a specific keyword containing well IDs
     * From filename: Extract from filename only
   - Validates extraction by listing all matched well IDs
   - Stops if no well IDs can be extracted

3. Gate Selection:
   - Select gate path from available gating hierarchies
   - Select specific gate name within that path

4. Plot Configuration:
   - Choose plot type: Histogram or Scatter
   - For Histogram:
     * Enter channel keyword (e.g., "RL1-A", "pHrodo")
     * Select statistic to display: Median or Mean
     * Select x-axis scale: Linear or Log
   - For Scatter:
     * Enter X-axis channel keyword (e.g., "FSC-A")
     * Enter Y-axis channel keyword (e.g., "SSC-A")

5. Output:
   - Generates plots for all samples matching the filter criteria
   - Organizes plots in a 96-well plate grid (8 rows x 12 columns)
   - Saves to HTML file with filename based on: input_name + plot_type + gate_name + parameters + filter
   - Includes modal popup with plot configuration and processing summary

Examples:
---------
    # Interactive mode - histogram of RL1-A for Singlets gate
    python analyze_flow.py --wsp workspace.wsp --fcs_dir ./fcs_files --interactive
    
    # Interactive mode with custom output
    python analyze_flow.py --wsp workspace.wsp --fcs_dir ./fcs_files --interactive --output my_plots.html

Output Files:
-------------
- Single HTML file with grid of plots organized by well position
- Filename format: {workspace_name}_{plot_type}_{gate_name}_{parameters}_{filter}.html
- Example: phrodo011_histogram_Singlets_RL1-A_Time_point_hr_0.html
- If --output is specified and not "report.html", uses that filename instead

Dependencies:
-------------
    - flowkit: FlowJo workspace parsing and FCS file handling
    - bokeh: Interactive HTML visualizations
    - pandas: Data manipulation and statistics
    - numpy: Numerical operations
    - matplotlib: Optional, for gate filtering with raw values

Notes:
------
- All data extraction uses raw (untransformed) values for accuracy
- Well IDs are parsed using fuzzy matching (handles "WellID", "Well ID", etc.)
- Well ID formats "A1" and "A01" are treated as equivalent
- Log scale histograms have fixed x-axis range: 0.1 to 10^5
- Linear scale histograms use dynamic x-axis based on data with outlier filtering
"""

def parse_args():
    parser = argparse.ArgumentParser(description="Analyze FlowJo Workspace and generate interactive HTML plots or inspect gate hierarchy.")
    parser.add_argument("--wsp", required=True, help="Path to the FlowJo .wsp file")
    parser.add_argument("--fcs_dir", help="Directory containing .fcs files (if different from wsp directory)")
    parser.add_argument("--output", default="report.html", help="Output HTML file path (auto-generated if not specified)")
    parser.add_argument("--interactive", action="store_true", help="Run interactive plotting mode")
    parser.add_argument("--inspect", action="store_true", help="Run inspect mode to view gate hierarchy")
    parser.add_argument("--sample", help="Optional sample ID to inspect (inspect mode only)")
    return parser.parse_args()

class FlowAnalyzer:
    """
    Main class for analyzing FlowJo workspaces and generating interactive HTML plots.
    
    This class provides batch plot generation capabilities for FlowJo workspaces. All data
    processing uses raw (untransformed) values for accuracy.
    
    Attributes:
        wsp_path (str): Path to the FlowJo workspace (.wsp) file.
        fcs_dir (str): Directory containing FCS files referenced by the workspace.
        workspace (flowkit.Workspace): Loaded FlowKit workspace object providing access to
            samples, gates, and events.
        sample_groups (list): List of sample group names found in the workspace
            (e.g., ['All Samples', 'Group1']).
    
    Main Methods:
        - interactive_plot_prompt(): Interactive CLI for batch plot configuration
        - generate_interactive_plots(): Generate batch plots and save to HTML
        - parse_well_id(): Extract well position from sample metadata
        - _get_all_gates_info(): Extract gate information for interactive selection
        - _plot_histogram(): Generate histogram plot for a single channel
        - _plot_scatter(): Generate scatter plot for two channels
    
    Example:
        >>> analyzer = FlowAnalyzer("workspace.wsp", "./fcs_files")
        >>> selections = analyzer.interactive_plot_prompt()
        >>> analyzer.generate_interactive_plots(selections, "plots.html")
    """
    def __init__(self, wsp_path, fcs_dir=None):
        """
        Initialize the FlowAnalyzer.

        Args:
            wsp_path (str): Path to the FlowJo workspace file.
            fcs_dir (str, optional): Directory containing FCS files. Defaults to the workspace directory.
        """
        self.wsp_path = wsp_path
        self.fcs_dir = fcs_dir if fcs_dir else os.path.dirname(wsp_path)
        
        print(f"Loading workspace: {self.wsp_path}")
        try:
            self.workspace = fk.Workspace(self.wsp_path, fcs_samples=self.fcs_dir)
        except Exception as e:
            print(f"Error loading workspace: {e}")
            # Fallback: try loading without specifying fcs_samples immediately if that fails, 
            # though usually it's better to provide it.
            # If fcs_dir is wrong, FlowKit might complain.
            raise e

        self.sample_groups = self.workspace.get_sample_groups()
        print(f"Found sample groups: {self.sample_groups}")

    # Standard report functions removed - only interactive mode is supported

    def parse_well_id(self, sample, source='auto', keyword_name=None, return_method=False, sample_id=None):
        """
        Parse the Well ID (e.g., 'A01', 'H12') from sample metadata.
        
        Extracts well position information from FlowJo sample metadata using flexible matching.
        Supports multiple extraction methods and fuzzy keyword matching. Handles various well ID
        formats including "A1", "A01", "A_1", etc., treating them as equivalent.
        
        Args:
            sample (flowkit.Sample): The sample object to parse. Must have either a 'keywords'
                attribute or be accessible via workspace.get_keywords().
            source (str, optional): Where to extract well ID from. Options:
                - 'auto': Try keywords first (fuzzy match for 'wellid'), then filename (default)
                - 'keyword': Extract from a specific keyword only
                - 'filename': Extract from filename only
                Defaults to 'auto'.
            keyword_name (str, optional): If source='keyword', the name of the keyword to use.
                Supports fuzzy matching (e.g., "WellID" matches "Well ID", "wellid", etc.).
                If None and source='keyword', will search for 'wellid' using fuzzy match.
                Defaults to None.
            return_method (bool, optional): If True, returns (row, col, method_used) instead of
                (row, col). The method string describes how the well ID was extracted.
                Defaults to False.
            sample_id (str, optional): Sample ID to use for workspace.get_keywords() if
                sample.keywords is not available. If None, attempts to extract from sample
                attributes. Defaults to None.

        Returns:
            tuple: Well position information:
                - If return_method=False: (Row, Col) where:
                    * Row is a string (e.g., 'A', 'B', 'H')
                    * Col is an int (e.g., 1, 12)
                    * Returns (None, None) if parsing fails
                - If return_method=True: (Row, Col, method) where:
                    * Row and Col are as above
                    * method is a string like "keyword 'WellID'", "filename", etc.
                    * Returns (None, None, None) if parsing fails
        
        Raises:
            AttributeError: If sample doesn't have required attributes and sample_id is not provided.
        
        Examples:
            >>> # Auto mode - tries keywords then filename
            >>> row, col = analyzer.parse_well_id(sample)
            >>> print(f"Well: {row}{col:02d}")  # Well: A01
            
            >>> # Extract from specific keyword with method tracking
            >>> row, col, method = analyzer.parse_well_id(
            ...     sample, source='keyword', keyword_name='WellID', return_method=True
            ... )
            >>> print(f"Well: {row}{col:02d} extracted from {method}")
            >>> # Well: A01 extracted from keyword 'WellID'
            
            >>> # Extract from filename only
            >>> row, col = analyzer.parse_well_id(sample, source='filename')
            >>> print(f"Well: {row}{col:02d}")  # Well: B12
        
        Notes:
            - Well ID formats "A1" and "A01" are treated as equivalent (both become row='A', col=1)
            - Keyword matching is case-insensitive and handles spaces/symbols
              (e.g., "WellID", "Well ID", "well_id" all match)
            - Filename parsing uses regex pattern: letter A-H followed by 1-2 digits
        """
        import re
        
        # Helper function to parse well ID from a string
        def parse_well_string(s):
            """Extract (row, col) from a string like 'A01', 'A1', 'A_01', etc.
            Handles both 'A1' and 'A01' formats - they are treated as the same.
            """
            s_str = str(s).strip()
            # Pattern: letter A-H (case insensitive) followed by optional separator and 1-2 digits
            # This matches: A1, A01, A_1, A-01, etc.
            match = re.search(r'([A-H])\W*(\d{1,2})', s_str, re.IGNORECASE)
            if match:
                row = match.group(1).upper()
                col = int(match.group(2))  # Convert to int (A1 and A01 both become 1)
                return row, col
            return None, None
        
        # Normalize keys: remove symbols, spaces, lowercase - for fuzzy matching
        def normalize(s):
            """Normalize string for fuzzy matching: remove all non-alphanumeric, lowercase"""
            return re.sub(r'[^a-zA-Z0-9]', '', s).lower()
        
        # Check if a keyword name matches "wellid" (fuzzy match)
        def matches_wellid(key):
            """Check if keyword matches 'wellid' with fuzzy matching"""
            normalized = normalize(key)
            # Check for exact match or contains "well" and "id"
            return normalized == "wellid" or (normalized.startswith("well") and "id" in normalized)
        
        # 1. Extract from keyword
        if source in ['auto', 'keyword']:
            # Get keywords - try multiple methods
            keywords = None
            # First try: sample.keywords attribute
            if hasattr(sample, 'keywords') and sample.keywords:
                keywords = sample.keywords
            # Second try: workspace.get_keywords() method (most reliable for FlowKit)
            elif hasattr(self, 'workspace'):
                try:
                    # Use provided sample_id or try to get it from sample
                    sid = sample_id
                    if not sid:
                        if hasattr(sample, 'id'):
                            sid = sample.id
                        elif hasattr(sample, 'sample_id'):
                            sid = sample.sample_id
                    if sid:
                        keywords = self.workspace.get_keywords(sid)
                except Exception:
                    pass
            # Third try: sample.get_keywords() method
            if not keywords and hasattr(sample, 'get_keywords'):
                try:
                    keywords = sample.get_keywords()
                except:
                    pass
            
            if keywords:
                # If specific keyword name provided, use it (with fuzzy matching)
                if source == 'keyword' and keyword_name:
                    # Try exact match first
                    if keyword_name in keywords:
                        value = keywords[keyword_name]
                        row, col = parse_well_string(value)
                        if row and col:
                            method = f"keyword '{keyword_name}'"
                            return (row, col, method) if return_method else (row, col)
                    else:
                        # Try fuzzy match - find keyword that matches the provided name
                        # First try normalized exact match
                        keyword_normalized = normalize(keyword_name)
                        for key, value in keywords.items():
                            if normalize(key) == keyword_normalized:
                                row, col = parse_well_string(value)
                                if row and col:
                                    method = f"keyword '{key}' (matched '{keyword_name}')"
                                    return (row, col, method) if return_method else (row, col)
                        
                        # Then try case-insensitive partial match
                        for key, value in keywords.items():
                            if (keyword_name.lower() in key.lower() or 
                                key.lower() in keyword_name.lower() or
                                normalize(keyword_name) in normalize(key) or
                                normalize(key) in normalize(keyword_name)):
                                row, col = parse_well_string(value)
                                if row and col:
                                    method = f"keyword '{key}' (matched '{keyword_name}')"
                                    return (row, col, method) if return_method else (row, col)
                else:
                    # Search for 'wellid' (fuzzy match) - original behavior
                    for key, value in keywords.items():
                        if matches_wellid(key):
                            row, col = parse_well_string(value)
                            if row and col:
                                method = f"keyword '{key}' (auto-detected)"
                                return (row, col, method) if return_method else (row, col)
                            break  # Found the key, even if parsing failed
        
        # 2. Extract from filename (if auto mode or explicitly requested)
        if source in ['auto', 'filename']:
            filename = sample.original_filename
            row, col = parse_well_string(filename)
            if row and col:
                method = "filename"
                return (row, col, method) if return_method else (row, col)
        
        return (None, None, None) if return_method else (None, None)

    # Standard report function removed: generate_interactive_heatmap
    # Standard report function removed: get_gate_polygons
    # Standard report function removed: generate_multi_sample_plot
    def _get_all_gates_info(self, sample_id):
        """
        Extract all gate information from a sample for interactive selection.
        
        Traverses the gating hierarchy for a sample and groups gates by their hierarchical
        path. Includes "Ungated" as a special root-level option. Used internally by
        interactive_plot_prompt() to display available gates to the user.
        
        Args:
            sample_id (str): The sample ID to extract gates from. Must be a valid sample
                ID in the workspace.
        
        Returns:
            dict: Dictionary mapping gate path strings to lists of (gate_name, gate_path) tuples.
                - Keys are string representations of gate paths (e.g., "root", "root → Cells")
                - Values are lists of tuples: (gate_name, gate_path_tuple)
                - Example: {
                    "root": [("Ungated", ())],
                    "root → Cells": [("Cells", ("root", "Cells"))],
                    "root → Cells → Singlets": [("Singlets", ("root", "Cells", "Singlets"))]
                  }
        
        Raises:
            ValueError: If sample_id is not found in the workspace.
        
        Example:
            >>> gates_info = analyzer._get_all_gates_info("A1.fcs")
            >>> for path, gates in gates_info.items():
            ...     print(f"{path}: {[g[0] for g in gates]}")
            >>> # root: ['Ungated']
            >>> # root → Cells: ['Cells']
            >>> # root → Cells → Singlets: ['Singlets']
        
        Notes:
            - "Ungated" is always included as the root-level option
            - Gate paths are represented as tuples for programmatic use
            - Path strings use " → " separator for display purposes
        """
        gate_ids = self.workspace.get_gate_ids(sample_id)
        gates_by_path = {}
        
        # Add "Ungated" option
        gates_by_path["root"] = [("Ungated", ())]
        
        for gate_name, gate_path in gate_ids:
            path_str = " → ".join(gate_path) if gate_path else "root"
            if path_str not in gates_by_path:
                gates_by_path[path_str] = []
            gates_by_path[path_str].append((gate_name, gate_path))
        
        return gates_by_path

    def _extract_quadrant_gate(self, sample_id, gate_name, gate_path, x_channel, y_channel):
        """
        Extract quadrant gate divider information for visualization on scatter plots.
        
        Args:
            sample_id (str): Sample ID to extract gate from
            gate_name (str): Name of the quadrant gate
            gate_path (tuple): Path tuple for the gate
            x_channel (str): X-axis channel name (e.g., "FSC-A")
            y_channel (str): Y-axis channel name (e.g., "SSC-A")
        
        Returns:
            dict or None: Dictionary containing quadrant gate data, or None if not a quadrant gate:
                - 'name': Gate name
                - 'type': 'quadrant'
                - 'x_dim': X-axis channel name
                - 'y_dim': Y-axis channel name
                - 'dividers': List of divider dictionaries with:
                    - 'dimension': Channel name
                    - 'value': Divider value (in raw data space)
                    - 'orientation': 'vertical' or 'horizontal'
        """
        try:
            # Handle "Ungated" case
            if gate_name == "Ungated":
                return None

            # Get the gate using the parent path
            get_gate_path = gate_path[:-1] if len(gate_path) > 1 else ()
            gate = self.workspace.get_gate(sample_id, gate_name, gate_path=get_gate_path)
            
            # Check if this is a quadrant gate
            # QuadrantGate has 'dividers' and/or 'quadrants' attributes
            has_dividers = hasattr(gate, 'dividers') and len(gate.dividers) > 0
            has_quadrants = hasattr(gate, 'quadrants') and len(gate.quadrants) > 0
            
            # If this gate doesn't have dividers but looks like a quadrant region (Q1-Q4),
            # try to find a parent QuadrantGate or sibling gate with dividers
            if not (has_dividers or has_quadrants):
                # Check if this is a quadrant region gate (Q1, Q2, Q3, Q4)
                is_quadrant_region = False
                try:
                    import re
                    if re.match(r'^Q\d+:', gate_name):
                        is_quadrant_region = True
                except:
                    pass
                
                if is_quadrant_region:
                    print(f"    DEBUG: Gate '{gate_name}' is a quadrant region, searching for QuadrantGate with dividers...")
                    # Search for gates that might have dividers
                    # Check: 1) Same path level (siblings), 2) Parent path level
                    all_gate_ids = self.workspace.get_gate_ids(sample_id)
                    parent_quadrant_gate = None
                    parent_path = gate_path[:-1] if len(gate_path) > 0 else ()
                    
                    # First, try to find a gate at the same path level with dividers
                    for other_gate_name, other_gate_path in all_gate_ids:
                        # Skip the current gate
                        if other_gate_name == gate_name and other_gate_path == gate_path:
                            continue
                            
                        # Check gates at same path level OR parent path level
                        if other_gate_path == gate_path or other_gate_path == parent_path:
                            try:
                                other_get_gate_path = other_gate_path[:-1] if len(other_gate_path) > 1 else ()
                                other_gate = self.workspace.get_gate(sample_id, other_gate_name, gate_path=other_get_gate_path)
                                
                                # Check if this gate has dividers and matches our channels
                                if hasattr(other_gate, 'dividers') and len(other_gate.dividers) > 0:
                                    if len(other_gate.dimensions) == 2:
                                        other_dims = [d.id for d in other_gate.dimensions]
                                        if ((other_dims[0] == x_channel and other_dims[1] == y_channel) or
                                            (other_dims[0] == y_channel and other_dims[1] == x_channel)):
                                            parent_quadrant_gate = other_gate
                                            print(f"    DEBUG: Found QuadrantGate '{other_gate_name}' at path {other_gate_path} with {len(other_gate.dividers)} divider(s)")
                                            break
                            except Exception as e:
                                print(f"    DEBUG: Error checking gate '{other_gate_name}': {e}")
                                continue
                    
                    if parent_quadrant_gate:
                        gate = parent_quadrant_gate
                        has_dividers = hasattr(gate, 'dividers') and len(gate.dividers) > 0
                        has_quadrants = hasattr(gate, 'quadrants') and len(gate.quadrants) > 0
                    else:
                        print(f"    DEBUG: Could not find QuadrantGate with dividers for '{gate_name}'")
                        print(f"           Attempting to infer dividers from Q1-Q4 quadrant regions...")
                        # Fallback: Try to infer dividers from Q1-Q4 gates themselves
                        inferred_dividers = self._infer_quadrant_dividers_from_regions(
                            sample_id, gate_path, x_channel, y_channel
                        )
                        if inferred_dividers:
                            print(f"    DEBUG: Successfully inferred {len(inferred_dividers)} divider(s) from quadrant regions")
                            # Get dimensions from the current gate
                            if len(gate.dimensions) == 2:
                                dims = [d.id for d in gate.dimensions]
                                return {
                                    'name': f"Quadrants ({gate_name})",
                                    'type': 'quadrant',
                                    'x_dim': dims[0],
                                    'y_dim': dims[1],
                                    'dividers': inferred_dividers
                                }
                        print(f"    DEBUG: Could not infer dividers from quadrant regions")
                        return None
                else:
                    print(f"    DEBUG: Gate '{gate_name}' does not have dividers or quadrants attributes")
                    return None
            
            divider_count = len(gate.dividers) if has_dividers else 0
            quadrant_count = len(gate.quadrants) if has_quadrants else 0
            print(f"    DEBUG: Detected quadrant gate '{gate_name}' with {divider_count} divider(s) and {quadrant_count} quadrant(s)")
            
            # Check dimensions
            if len(gate.dimensions) != 2:
                return None
            
            dims = [d.id for d in gate.dimensions]  # Channel names
            
            # Check if gate matches our axes (in either order)
            if not ((dims[0] == x_channel and dims[1] == y_channel) or
                   (dims[0] == y_channel and dims[1] == x_channel)):
                return None
            
            # Extract dividers
            dividers = []
            divider_list = gate.dividers if hasattr(gate, 'dividers') else []
            for divider in divider_list:
                try:
                    # Get divider dimension and value
                    # Try multiple possible attribute names for divider properties
                    divider_dim = None
                    divider_value = None
                    
                    # Try different possible attribute names for dimension reference
                    for attr_name in ['dimension_ref', 'dimension', 'dim_ref', 'dim']:
                        if hasattr(divider, attr_name):
                            attr_value = getattr(divider, attr_name)
                            # Could be a string or a Dimension object
                            if isinstance(attr_value, str):
                                divider_dim = attr_value
                            elif hasattr(attr_value, 'id'):
                                divider_dim = attr_value.id
                            elif hasattr(attr_value, '__str__'):
                                divider_dim = str(attr_value)
                            if divider_dim:
                                break
                    
                    # Try different possible attribute names for divider value
                    for attr_name in ['value', 'position', 'divider_value', 'threshold']:
                        if hasattr(divider, attr_name):
                            divider_value = getattr(divider, attr_name)
                            if divider_value is not None:
                                break
                    
                    if divider_dim and divider_value is not None:
                        # Determine orientation
                        if divider_dim == x_channel:
                            orientation = 'vertical'
                        elif divider_dim == y_channel:
                            orientation = 'horizontal'
                        else:
                            # Divider doesn't match our axes
                            continue
                        
                        # Convert divider value from display space to raw data space if needed
                        # Similar to polygon gate conversion
                        try:
                            # Divider values might be in display space (0-1) or already in raw space
                            # Try to convert if it's in display space (values between 0 and 1)
                            if 0 <= divider_value <= 1:
                                # Convert from display space to raw values
                                display_min_log = 0  # log10(1)
                                display_max_log = 5  # log10(100000)
                                divider_value = 10 ** (display_min_log + divider_value * (display_max_log - display_min_log))
                        except:
                            # If conversion fails, use value as-is (might already be in raw space)
                            pass
                        
                        dividers.append({
                            'dimension': divider_dim,
                            'value': divider_value,
                            'orientation': orientation
                        })
                except Exception as e:
                    print(f"    DEBUG: Could not extract divider: {e}")
                    continue
            
            if dividers:
                return {
                    'name': gate_name,
                    'type': 'quadrant',
                    'x_dim': dims[0],
                    'y_dim': dims[1],
                    'dividers': dividers
                }
        except Exception as e:
            print(f"    DEBUG: Could not extract quadrant gate '{gate_name}': {e}")
        
        return None
    
    def _infer_quadrant_dividers_from_regions(self, sample_id, gate_path, x_channel, y_channel):
        """
        Infer quadrant divider positions from Q1-Q4 region gates.
        
        Args:
            sample_id (str): Sample ID
            gate_path (tuple): Path tuple for the quadrant gates
            x_channel (str): X-axis channel name
            y_channel (str): Y-axis channel name
        
        Returns:
            list or None: List of divider dictionaries, or None if inference fails
        """
        try:
            import re
            import numpy as np

            print(f"\n    DEBUG _infer_quadrant_dividers_from_regions:")
            print(f"           sample_id: {sample_id}")
            print(f"           gate_path: {gate_path}")
            print(f"           x_channel: {x_channel}")
            print(f"           y_channel: {y_channel}")

            # Find all Q1-Q4 gates at this path
            all_gate_ids = self.workspace.get_gate_ids(sample_id)
            quadrant_gates = {}
            get_gate_path = gate_path[:-1] if len(gate_path) > 1 else ()

            print(f"           Checking {len(all_gate_ids)} total gates...")
            gates_at_target_path = [(name, path) for name, path in all_gate_ids if path == gate_path]
            print(f"           Found {len(gates_at_target_path)} gates at target path {gate_path}")
            
            for gate_name, other_gate_path in all_gate_ids:
                if other_gate_path == gate_path:
                    # Check if it's a Q1-Q4 gate
                    if re.match(r'^Q\d+:', gate_name):
                        try:
                            gate = self.workspace.get_gate(sample_id, gate_name, gate_path=get_gate_path)
                            if len(gate.dimensions) == 2:
                                dims = [d.id for d in gate.dimensions]
                                # Check if dimensions match
                                if ((dims[0] == x_channel and dims[1] == y_channel) or
                                    (dims[0] == y_channel and dims[1] == x_channel)):
                                    # Extract quadrant number
                                    match = re.match(r'^Q(\d+):', gate_name)
                                    if match:
                                        q_num = int(match.group(1))
                                        quadrant_gates[q_num] = {
                                            'name': gate_name,
                                            'gate': gate,
                                            'dims': dims
                                        }
                        except Exception as e:
                            print(f"    DEBUG: Error getting gate '{gate_name}': {e}")
                            continue
            
            if len(quadrant_gates) < 2:
                print(f"    DEBUG: Found only {len(quadrant_gates)} quadrant gates, need at least 2")
                return None
            
            # Debug: Print what attributes the Q gates have
            print(f"    DEBUG: Found {len(quadrant_gates)} quadrant gates: {list(quadrant_gates.keys())}")

            # NEW APPROACH: Parse XML directly to extract min/max values
            # Import the helper function
            from quadrant_xml_parser import infer_quadrant_dividers_from_xml

            result = infer_quadrant_dividers_from_xml(
                self.wsp_path,
                sample_id,
                gate_path,
                x_channel,
                y_channel,
                quadrant_gates
            )

            return result

            # FALLBACK: If XML parsing fails, try other methods
            gate_boundaries = {}
            for q_num, q_info in quadrant_gates.items():
                gate = q_info['gate']
                print(f"    DEBUG: Q{q_num} ('{q_info['name']}'):")

                # Try to get min/max from the gate
                if hasattr(gate, 'vertices') and gate.vertices is not None and len(gate.vertices) > 0:
                    # Polygon-style gate - extract bounding box
                    verts = np.array(gate.vertices)
                    gate_boundaries[q_num] = {
                        'x_min': float(verts[:, 0].min()),
                        'x_max': float(verts[:, 0].max()),
                        'y_min': float(verts[:, 1].min()),
                        'y_max': float(verts[:, 1].max()),
                        'source': 'vertices'
                    }
                    print(f"           Got boundaries from vertices: X=[{gate_boundaries[q_num]['x_min']:.4f}, {gate_boundaries[q_num]['x_max']:.4f}], Y=[{gate_boundaries[q_num]['y_min']:.4f}, {gate_boundaries[q_num]['y_max']:.4f}]")
                elif hasattr(gate, 'min') and hasattr(gate, 'max'):
                    # Rectangle gate with min/max attributes
                    dims = q_info['dims']
                    try:
                        min_vals = gate.min
                        max_vals = gate.max

                        # Extract min/max for our dimensions
                        if isinstance(min_vals, dict) and isinstance(max_vals, dict):
                            x_min = min_vals.get(x_channel, None)
                            x_max = max_vals.get(x_channel, None)
                            y_min = min_vals.get(y_channel, None)
                            y_max = max_vals.get(y_channel, None)
                        else:
                            # Assume it's array-like [x, y]
                            x_idx = 0 if dims[0] == x_channel else 1
                            y_idx = 1 - x_idx
                            x_min = min_vals[x_idx] if hasattr(min_vals, '__getitem__') else None
                            x_max = max_vals[x_idx] if hasattr(max_vals, '__getitem__') else None
                            y_min = min_vals[y_idx] if hasattr(min_vals, '__getitem__') else None
                            y_max = max_vals[y_idx] if hasattr(max_vals, '__getitem__') else None

                        if all(v is not None for v in [x_min, x_max, y_min, y_max]):
                            gate_boundaries[q_num] = {
                                'x_min': float(x_min),
                                'x_max': float(x_max),
                                'y_min': float(y_min),
                                'y_max': float(y_max),
                                'source': 'min/max'
                            }
                            print(f"           Got boundaries from min/max: X=[{x_min:.4f}, {x_max:.4f}], Y=[{y_min:.4f}, {y_max:.4f}]")
                    except Exception as e:
                        print(f"           Failed to extract min/max: {e}")

            if len(gate_boundaries) >= 2:
                print(f"    DEBUG: Successfully extracted {len(gate_boundaries)} gate boundaries")
                print(f"    DEBUG: Converting from display space (0-1) to raw data space...")

                # Transform from display space to raw values (same as polygon gates)
                # FlowJo typically displays log scale from 10^0 to 10^5 (1 to 100,000)
                display_min_log = 0  # log10(1)
                display_max_log = 5  # log10(100000)

                transformed_boundaries = {}
                for q_num, bounds in gate_boundaries.items():
                    transformed_boundaries[q_num] = {
                        'x_min': 10 ** (display_min_log + bounds['x_min'] * (display_max_log - display_min_log)),
                        'x_max': 10 ** (display_min_log + bounds['x_max'] * (display_max_log - display_min_log)),
                        'y_min': 10 ** (display_min_log + bounds['y_min'] * (display_max_log - display_min_log)),
                        'y_max': 10 ** (display_min_log + bounds['y_max'] * (display_max_log - display_min_log)),
                    }
                    print(f"    DEBUG: Q{q_num} transformed: X=[{transformed_boundaries[q_num]['x_min']:.2f}, {transformed_boundaries[q_num]['x_max']:.2f}], Y=[{transformed_boundaries[q_num]['y_min']:.2f}, {transformed_boundaries[q_num]['y_max']:.2f}]")

                # Now infer dividers from the transformed boundaries
                quadrant_data = transformed_boundaries
            else:
                print(f"    DEBUG: Could not extract enough gate boundaries ({len(gate_boundaries)})")
                print(f"    DEBUG: Direct boundary extraction failed - cannot infer quadrant dividers")
                return None

            # quadrant_data now contains the transformed boundaries
            # Continue with divider inference below
            if len(quadrant_data) < 2:
                print(f"    DEBUG: Could not get data from enough quadrants ({len(quadrant_data)})")
                return None

            # Infer divider positions from quadrant boundaries
            # Q1: FITC-A- (X-), APC-Vio770-A+ (Y+) -> X < divider_x, Y > divider_y
            # Q2: FITC-A+ (X+), APC-Vio770-A+ (Y+) -> X > divider_x, Y > divider_y
            # Q3: FITC-A+ (X+), APC-Vio770-A- (Y-) -> X > divider_x, Y < divider_y
            # Q4: FITC-A- (X-), APC-Vio770-A- (Y-) -> X < divider_x, Y < divider_y
            
            dividers = []
            
            # X divider: boundary between Q1/Q4 (X-) and Q2/Q3 (X+)
            # Should be between max X of Q1/Q4 and min X of Q2/Q3
            x_negative_max = []  # Q1, Q4
            x_positive_min = []  # Q2, Q3
            
            for q_num, data in quadrant_data.items():
                if q_num in [1, 4]:  # Negative X
                    x_negative_max.append(data['x_max'])
                elif q_num in [2, 3]:  # Positive X
                    x_positive_min.append(data['x_min'])
            
            if x_negative_max and x_positive_min:
                x_divider = (max(x_negative_max) + min(x_positive_min)) / 2.0
                dividers.append({
                    'dimension': x_channel,
                    'value': float(x_divider),
                    'orientation': 'vertical'
                })
                print(f"    DEBUG: Inferred X divider at {x_divider:.2f} for {x_channel}")
            
            # Y divider: boundary between Q1/Q2 (Y+) and Q3/Q4 (Y-)
            # Should be between max Y of Q3/Q4 and min Y of Q1/Q2
            y_negative_max = []  # Q3, Q4
            y_positive_min = []  # Q1, Q2
            
            for q_num, data in quadrant_data.items():
                if q_num in [3, 4]:  # Negative Y
                    y_negative_max.append(data['y_max'])
                elif q_num in [1, 2]:  # Positive Y
                    y_positive_min.append(data['y_min'])
            
            if y_negative_max and y_positive_min:
                y_divider = (max(y_negative_max) + min(y_positive_min)) / 2.0
                dividers.append({
                    'dimension': y_channel,
                    'value': float(y_divider),
                    'orientation': 'horizontal'
                })
                print(f"    DEBUG: Inferred Y divider at {y_divider:.2f} for {y_channel}")
            
            if len(dividers) > 0:
                return dividers
            else:
                return None
        except Exception as e:
            print(f"    DEBUG: Error in _infer_quadrant_dividers_from_regions: {e}")
            import traceback
            traceback.print_exc()
            return None
                    
    
    def _extract_quadrant_thresholds_from_parent(self, sample_id, gate_path, x_channel, y_channel):
        """
        Extract X and Y thresholds from a parent QuadrantGate for quadrant mode visualization.
        
        Searches for a QuadrantGate at the same path level or parent path level that matches
        the specified channels. Extracts divider values and converts them from display space
        to raw data space.
        
        Also checks if any gate at the path has quadrant-related attributes that might contain
        divider information.
        
        Args:
            sample_id (str): Sample ID to extract gate from
            gate_path (tuple): Path tuple where Q1-Q4 gates are located (e.g., ('root', 'Cells', 'Singlets'))
            x_channel (str): X-axis channel name (e.g., "B1-A")
            y_channel (str): Y-axis channel name (e.g., "R2-A")
        
        Returns:
            dict or None: Dictionary with 'x_threshold' and 'y_threshold' in raw data space, or None if not found
        """
        try:
            # Search for a QuadrantGate with dividers at the same path or parent path
            all_gate_ids = self.workspace.get_gate_ids(sample_id)
            parent_path = gate_path[:-1] if len(gate_path) > 0 else ()
            
            # Debug: Print what gates we're checking
            print(f"    DEBUG: Searching for QuadrantGate with dividers...")
            print(f"           Looking at path: {gate_path}")
            print(f"           Parent path: {parent_path}")
            print(f"           Channels: {x_channel} × {y_channel}")
            
            # Search at same path level and parent path level
            checked_gates = []
            for other_gate_name, other_gate_path in all_gate_ids:
                # Check gates at same path level OR parent path level
                if other_gate_path == gate_path or other_gate_path == parent_path:
                    try:
                        other_get_gate_path = other_gate_path[:-1] if len(other_gate_path) > 1 else ()
                        other_gate = self.workspace.get_gate(sample_id, other_gate_name, gate_path=other_get_gate_path)
                        checked_gates.append((other_gate_name, other_gate_path, hasattr(other_gate, 'dividers')))
                        
                        # Check if this gate has dividers and matches our channels
                        has_dividers = hasattr(other_gate, 'dividers') and len(other_gate.dividers) > 0
                        if has_dividers:
                            print(f"    DEBUG: Found gate '{other_gate_name}' with dividers at path {other_gate_path}")
                            if len(other_gate.dimensions) == 2:
                                other_dims = [d.id for d in other_gate.dimensions]
                                if ((other_dims[0] == x_channel and other_dims[1] == y_channel) or
                                    (other_dims[0] == y_channel and other_dims[1] == x_channel)):
                                    
                                    # Extract dividers
                                    x_threshold = None
                                    y_threshold = None
                                    
                                    divider_list = other_gate.dividers
                                    for divider in divider_list:
                                        try:
                                            # Get divider dimension and value
                                            divider_dim = None
                                            divider_value = None
                                            
                                            # Try different possible attribute names for dimension reference
                                            for attr_name in ['dimension_ref', 'dimension', 'dim_ref', 'dim']:
                                                if hasattr(divider, attr_name):
                                                    attr_value = getattr(divider, attr_name)
                                                    if isinstance(attr_value, str):
                                                        divider_dim = attr_value
                                                    elif hasattr(attr_value, 'id'):
                                                        divider_dim = attr_value.id
                                                    elif hasattr(attr_value, '__str__'):
                                                        divider_dim = str(attr_value)
                                                    if divider_dim:
                                                        break
                                            
                                            # Try different possible attribute names for divider value
                                            for attr_name in ['value', 'position', 'divider_value', 'threshold']:
                                                if hasattr(divider, attr_name):
                                                    divider_value = getattr(divider, attr_name)
                                                    if divider_value is not None:
                                                        break
                                            
                                            if divider_dim and divider_value is not None:
                                                # Convert divider value from display space to raw data space if needed
                                                try:
                                                    # Divider values might be in display space (0-1) or already in raw space
                                                    if 0 <= divider_value <= 1:
                                                        # Convert from display space to raw values
                                                        display_min_log = 0  # log10(1)
                                                        display_max_log = 5  # log10(100000)
                                                        divider_value = 10 ** (display_min_log + divider_value * (display_max_log - display_min_log))
                                                except:
                                                    # If conversion fails, use value as-is (might already be in raw space)
                                                    pass
                                                
                                                # Assign to x or y threshold based on dimension
                                                if divider_dim == x_channel:
                                                    x_threshold = divider_value
                                                elif divider_dim == y_channel:
                                                    y_threshold = divider_value
                                        except Exception as e:
                                            continue
                                    
                                    # Return thresholds if we found both
                                    if x_threshold is not None and y_threshold is not None:
                                        return {
                                            'x_threshold': float(x_threshold),
                                            'y_threshold': float(y_threshold)
                                        }
                    except Exception as e:
                        continue
            
            # Debug output
            if checked_gates:
                print(f"    DEBUG: Checked {len(checked_gates)} gates:")
                for gate_name, gate_path, has_div in checked_gates[:5]:  # Show first 5
                    print(f"           - '{gate_name}' at {gate_path}: has_dividers={has_div}")
            else:
                print(f"    DEBUG: No gates found at path {gate_path} or parent path {parent_path}")
            
            return None
        except Exception as e:
            return None
    
    def _extract_gate_polygons(self, sample_id, x_channel, y_channel):
        """
        Extract gate polygon coordinates for visualization on scatter plots.
        
        Uses get_gate_ids() to get all gates, then extracts polygon coordinates
        for gates that match the specified x and y channels. Also extracts quadrant gates.
        
        Args:
            sample_id (str): Sample ID to extract gates from
            x_channel (str): X-axis channel name (e.g., "FSC-A")
            y_channel (str): Y-axis channel name (e.g., "SSC-A")
        
        Returns:
            list: List of dictionaries, each containing:
                - 'name': Gate name
                - 'type': Gate type ('polygon', 'rectangle', or 'quadrant')
                - 'x_dim': X-axis channel name
                - 'y_dim': Y-axis channel name
                - 'xs': List of X coordinates (for polygon/rectangle gates)
                - 'ys': List of Y coordinates (for polygon/rectangle gates)
                - 'dividers': List of dividers (for quadrant gates)
        """
        gates_data = []
        gate_ids = self.workspace.get_gate_ids(sample_id)
        
        for gate_name, gate_path in gate_ids:
            try:
                gate = self.workspace.get_gate(sample_id, gate_name, gate_path=gate_path[:-1] if len(gate_path) > 1 else ())
                # Check dimensions
                if len(gate.dimensions) == 2:
                    dims = [d.id for d in gate.dimensions]  # Channel names
                    
                    # Only include gates that match our axes (in either order)
                    if (dims[0] == x_channel and dims[1] == y_channel) or \
                       (dims[0] == y_channel and dims[1] == x_channel):
                        
                        # Check if this is a quadrant gate first
                        if hasattr(gate, 'dividers') and hasattr(gate, 'quadrants'):
                            quadrant_gate = self._extract_quadrant_gate(sample_id, gate_name, gate_path, x_channel, y_channel)
                            if quadrant_gate:
                                gates_data.append(quadrant_gate)
                                continue
                        
                        # Get vertices for polygon/rectangle gates
                        xs = []
                        ys = []
                        
                        if hasattr(gate, 'vertices'):
                            # Polygon gate
                            verts = np.array(gate.vertices)

                            # Gate vertices are in transformed space, convert to raw data coordinates
                            try:
                                # Get the transforms from the workspace
                                x_transform = self.workspace.get_transform(sample_id, dims[0])
                                y_transform = self.workspace.get_transform(sample_id, dims[1])

                                print(f"    DEBUG (_extract_gate_polygons): Gate '{gate_name}'")
                                print(f"           Original verts range: X=[{verts[:, 0].min():.4f}, {verts[:, 0].max():.4f}], Y=[{verts[:, 1].min():.4f}, {verts[:, 1].max():.4f}]")
                                print(f"           X transform: {x_transform}")
                                print(f"           Y transform: {y_transform}")

                                # Gate vertices are in DISPLAY SPACE (0-1 normalized over visible range)
                                # NOT in transform space! We need to map to the display range.
                                # FlowJo typically displays log scale from 10^0 to 10^5 (1 to 100,000)

                                # Map from display space (0-1) to raw values
                                display_min_log = 0  # log10(1)
                                display_max_log = 5  # log10(100000)

                                # Convert normalized display coordinates to log space, then to linear
                                xs = (10 ** (display_min_log + verts[:, 0] * (display_max_log - display_min_log))).tolist()
                                ys = (10 ** (display_min_log + verts[:, 1] * (display_max_log - display_min_log))).tolist()

                                print(f"           Mapped from display space (0-1) to raw values")
                                print(f"           X: [{min(xs):.2f}, {max(xs):.2f}]")
                                print(f"           Y: [{min(ys):.2f}, {max(ys):.2f}]")
                            except Exception as e:
                                print(f"    DEBUG: Transform failed: {e}")
                                import traceback
                                traceback.print_exc()
                                # Fall back to original vertices
                                xs = verts[:, 0].tolist()
                                ys = verts[:, 1].tolist()

                            # Close the loop
                            xs.append(xs[0])
                            ys.append(ys[0])
                        elif hasattr(gate, 'min') and hasattr(gate, 'max'):
                            # Rectangle gate - create polygon from min/max
                            try:
                                min_x = None
                                min_y = None
                                max_x = None
                                max_y = None
                                
                                if isinstance(gate.min, dict):
                                    min_x = gate.min.get(dims[0], None)
                                    min_y = gate.min.get(dims[1], None)
                                elif hasattr(gate.min, '__getitem__'):
                                    try:
                                        min_x = gate.min[0]
                                        min_y = gate.min[1] if len(gate.min) > 1 else gate.min[0]
                                    except:
                                        pass
                                
                                if isinstance(gate.max, dict):
                                    max_x = gate.max.get(dims[0], None)
                                    max_y = gate.max.get(dims[1], None)
                                elif hasattr(gate.max, '__getitem__'):
                                    try:
                                        max_x = gate.max[0]
                                        max_y = gate.max[1] if len(gate.max) > 1 else gate.max[0]
                                    except:
                                        pass
                                
                                if min_x is not None and max_x is not None and min_y is not None and max_y is not None:
                                    # Create rectangle polygon (closed loop)
                                    xs = [float(min_x), float(max_x), float(max_x), float(min_x), float(min_x)]
                                    ys = [float(min_y), float(min_y), float(max_y), float(max_y), float(min_y)]
                            except Exception as e:
                                # If rectangle extraction fails, skip this gate
                                pass
                        
                        if xs and ys:
                            gates_data.append({
                                'name': gate_name,
                                'type': 'polygon' if hasattr(gate, 'vertices') else 'rectangle',
                                'x_dim': dims[0],
                                'y_dim': dims[1],
                                'xs': xs,
                                'ys': ys
                            })
            except Exception as e:
                # Log error but continue processing other gates
                print(f"    Warning: Could not extract gate '{gate_name}' for channels {x_channel}, {y_channel}: {e}")
        return gates_data

    def _extract_selected_gate(self, sample_id, gate_name, gate_path, x_channel, y_channel):
        """
        Extract polygon or quadrant gate coordinates for a specific gate for visualization on scatter plots.
        
        Args:
            sample_id (str): Sample ID to extract gate from
            gate_name (str): Name of the gate to extract
            gate_path (tuple): Path tuple for the gate (e.g., ('root', 'Cells'))
            x_channel (str): X-axis channel name (e.g., "FSC-A")
            y_channel (str): Y-axis channel name (e.g., "SSC-A")
        
        Returns:
            dict or None: Dictionary containing gate data, or None if gate not found:
                For polygon/rectangle gates:
                - 'name': Gate name
                - 'type': 'polygon' or 'rectangle'
                - 'x_dim': X-axis channel name
                - 'y_dim': Y-axis channel name
                - 'xs': List of X coordinates
                - 'ys': List of Y coordinates
                For quadrant gates:
                - 'name': Gate name
                - 'type': 'quadrant'
                - 'x_dim': X-axis channel name
                - 'y_dim': Y-axis channel name
                - 'dividers': List of divider dictionaries
        """
        try:
            # Handle "Ungated" case
            if gate_name == "Ungated":
                return None

            # Get the gate using the parent path
            get_gate_path = gate_path[:-1] if len(gate_path) > 1 else ()
            print(f"    DEBUG: Getting gate '{gate_name}' with gate_path={get_gate_path}")
            gate = self.workspace.get_gate(sample_id, gate_name, gate_path=get_gate_path)

            # Check if gate is 2D and matches our channels
            if len(gate.dimensions) != 2:
                print(f"    DEBUG: Gate is {len(gate.dimensions)}D, not 2D")
                return None

            dims = [d.id for d in gate.dimensions]  # Channel names
            print(f"    DEBUG: Gate dimensions: {dims}")
            print(f"    DEBUG: Expected channels: [{x_channel}, {y_channel}]")

            # Check if gate matches our axes (in either order)
            if not ((dims[0] == x_channel and dims[1] == y_channel) or
                   (dims[0] == y_channel and dims[1] == x_channel)):
                print(f"    DEBUG: Channel mismatch!")
                return None
            
            # Check if this is a quadrant gate first (comprehensive detection)
            is_quadrant = False
            
            # Method 1: Check for dividers/quadrants attributes
            if hasattr(gate, 'dividers') and hasattr(gate, 'quadrants'):
                if len(gate.dividers) > 0 or len(gate.quadrants) > 0:
                    is_quadrant = True
            elif hasattr(gate, 'dividers') and len(gate.dividers) > 0:
                is_quadrant = True
            elif hasattr(gate, 'quadrants') and len(gate.quadrants) > 0:
                is_quadrant = True
            
            # Method 2: Check gate class name
            if not is_quadrant:
                try:
                    gate_class_name = gate.__class__.__name__
                    if 'Quadrant' in gate_class_name:
                        is_quadrant = True
                except:
                    pass
            
            # Method 2b: Try isinstance check with flowkit gates
            if not is_quadrant:
                try:
                    import flowkit.gates as fk_gates
                    if isinstance(gate, (fk_gates.QuadrantGate, fk_gates.Quadrant)):
                        is_quadrant = True
                except (ImportError, AttributeError):
                    pass
            
            # Method 3: Heuristic based on gate name pattern (Q1:, Q2:, etc.)
            if not is_quadrant:
                try:
                    if gate_name.startswith('Q') and ':' in gate_name:
                        # Check if it's Q followed by a number
                        import re
                        if re.match(r'^Q\d+:', gate_name):
                            # Also check that it doesn't have vertices (polygon) or min/max (rectangle)
                            has_vertices = hasattr(gate, 'vertices') and len(gate.vertices) > 0
                            has_min_max = (hasattr(gate, 'min') and hasattr(gate, 'max'))
                            if not has_vertices and not has_min_max:
                                is_quadrant = True
                except:
                    pass
            
            if is_quadrant:
                quadrant_gate = self._extract_quadrant_gate(sample_id, gate_name, gate_path, x_channel, y_channel)
                if quadrant_gate:
                    print(f"    DEBUG: Extracted quadrant gate '{gate_name}' with {len(quadrant_gate['dividers'])} divider(s)")
                    return quadrant_gate
                else:
                    print(f"    DEBUG: Detected as quadrant but extraction failed")
            
            # Get vertices for polygon/rectangle gates
            xs = []
            ys = []
            
            if hasattr(gate, 'vertices'):
                # Polygon gate
                verts = np.array(gate.vertices)

                # Gate vertices are in transformed space, but we need raw data coordinates
                # Get the sample to access transforms
                sample = self.workspace.get_sample(sample_id)

                # Get transforms for the channels
                try:
                    # Apply inverse transform to convert from display coordinates to raw values
                    # The vertices are in transformed/display space (0-1 or similar)
                    # We need to convert them to raw data space

                    # Get channel indices
                    x_channel_idx = sample.pnn_labels.index(dims[0])
                    y_channel_idx = sample.pnn_labels.index(dims[1])

                    # Gate vertices are in DISPLAY SPACE (0-1 normalized over visible range)
                    # NOT in transform space! We need to map to the display range.
                    # FlowJo typically displays log scale from 10^0 to 10^5 (1 to 100,000)

                    # Map from display space (0-1) to raw values
                    display_min_log = 0  # log10(1)
                    display_max_log = 5  # log10(100000)

                    # Convert normalized display coordinates to log space, then to linear
                    xs = (10 ** (display_min_log + verts[:, 0] * (display_max_log - display_min_log))).tolist()
                    ys = (10 ** (display_min_log + verts[:, 1] * (display_max_log - display_min_log))).tolist()

                    print(f"    DEBUG: Mapped gate from display space to raw values")
                    print(f"           X: [{min(xs):.2f}, {max(xs):.2f}]")
                    print(f"           Y: [{min(ys):.2f}, {max(ys):.2f}]")
                except Exception as e:
                    print(f"    DEBUG: Could not apply inverse transform: {e}")
                    import traceback
                    traceback.print_exc()
                    # Fall back to original vertices
                    xs = verts[:, 0].tolist()
                    ys = verts[:, 1].tolist()

                # Close the loop
                xs.append(xs[0])
                ys.append(ys[0])
            elif hasattr(gate, 'min') and hasattr(gate, 'max'):
                # Rectangle gate - create polygon from min/max
                try:
                    min_x = None
                    min_y = None
                    max_x = None
                    max_y = None
                    
                    if isinstance(gate.min, dict):
                        min_x = gate.min.get(dims[0], None)
                        min_y = gate.min.get(dims[1], None)
                    elif hasattr(gate.min, '__getitem__'):
                        try:
                            min_x = gate.min[0]
                            min_y = gate.min[1] if len(gate.min) > 1 else gate.min[0]
                        except:
                            pass
                    
                    if isinstance(gate.max, dict):
                        max_x = gate.max.get(dims[0], None)
                        max_y = gate.max.get(dims[1], None)
                    elif hasattr(gate.max, '__getitem__'):
                        try:
                            max_x = gate.max[0]
                            max_y = gate.max[1] if len(gate.max) > 1 else gate.max[0]
                        except:
                            pass
                    
                    if min_x is not None and max_x is not None and min_y is not None and max_y is not None:
                        # Create rectangle polygon (closed loop)
                        xs = [float(min_x), float(max_x), float(max_x), float(min_x), float(min_x)]
                        ys = [float(min_y), float(min_y), float(max_y), float(max_y), float(min_y)]
                except Exception:
                    pass
            
            if xs and ys:
                # Validate coordinates before returning
                if len(xs) != len(ys):
                    print(f"    ⚠️  Warning: Gate '{gate_name}' has mismatched coordinates (xs={len(xs)}, ys={len(ys)})")
                    return None
                if len(xs) < 3:
                    print(f"    ⚠️  Warning: Gate '{gate_name}' has insufficient points ({len(xs)} < 3)")
                    return None

                # Show coordinate ranges for debugging
                print(f"    DEBUG: Gate '{gate_name}' coordinates:")
                print(f"           X range: [{min(xs):.2f}, {max(xs):.2f}]")
                print(f"           Y range: [{min(ys):.2f}, {max(ys):.2f}]")

                return {
                    'name': gate_name,
                    'type': 'polygon' if hasattr(gate, 'vertices') else 'rectangle',
                    'x_dim': dims[0],
                    'y_dim': dims[1],
                    'xs': xs,
                    'ys': ys
                }
        except Exception as e:
            print(f"    ⚠️  Warning: Could not extract gate '{gate_name}': {e}")

        return None

    def interactive_plot_prompt(self):
        """
        Interactive command-line prompt to gather user selections for plotting.
        
        Guides the user through selecting:
        1. Optional keyword filtering (yes/no, then select keyword and value)
        2. Gate path (hierarchical path in gating tree)
        3. Gate name (specific gate within the path)
        4. Plot type (histogram or scatter)
        5. Channel parameters (1 for histogram, 2 for scatter)
        
        Example interaction:
            Filter by keyword? (yes/no): yes
            Available keywords:
              1. Time point (hr): 0
              2. ...
            Select keyword: Time point (hr)
            Enter keyword value to filter: 0
            
            Available Gate Paths:
            1. root (gates: Ungated)
            2. root → Cells (gates: Cells)
            ...
        
        Returns:
            dict: Dictionary containing user selections:
                - 'gate_path': tuple of gate path (e.g., ("root", "Cells", "Singlets"))
                - 'gate_name': str gate name (e.g., "Singlets")
                - 'plot_type': "histogram" or "scatter"
                - 'parameters': list of channel keywords
                - 'keyword_filter': dict with 'key' and 'value' if filtering was applied
                - 'filtered_sample_ids': list of sample IDs after filtering
            Returns None if no samples found or user cancels.
        """
        if not self.sample_groups:
            print("No sample groups found in workspace.")
            return None
        
        # Step 0: Sample group selection
        print("\n" + "="*60)
        print("Sample Group Selection")
        print("="*60)
        print("Available sample groups:")
        for i, group_name in enumerate(self.sample_groups, 1):
            # Count samples in this group
            try:
                group_samples = self.workspace.get_sample_ids(group_name, loaded_only=True)
                print(f"  {i}. {group_name} ({len(group_samples)} samples)")
            except:
                print(f"  {i}. {group_name}")
        
        print(f"  {len(self.sample_groups) + 1}. All groups")
        
        selected_groups = []
        while True:
            try:
                group_choice = input(f"\nSelect sample group(s) (1-{len(self.sample_groups) + 1}, comma-separated for multiple, or 'all'): ").strip()
                if group_choice.lower() == 'all':
                    selected_groups = self.sample_groups
                    break
                else:
                    # Parse comma-separated numbers
                    choices = [c.strip() for c in group_choice.split(',')]
                    for choice in choices:
                        idx = int(choice) - 1
                        if 0 <= idx < len(self.sample_groups):
                            selected_groups.append(self.sample_groups[idx])
                        elif idx == len(self.sample_groups):  # "All groups" option
                            selected_groups = self.sample_groups
                            break
                    if selected_groups:
                        break
                    else:
                        print(f"Please enter numbers between 1 and {len(self.sample_groups) + 1}")
            except ValueError:
                print("Please enter valid numbers separated by commas, or 'all'")
        
        # Collect all sample IDs from selected groups
        all_sample_ids = []
        for group in selected_groups:
            try:
                group_samples = self.workspace.get_sample_ids(group, loaded_only=True)
                all_sample_ids.extend(group_samples)
            except Exception as e:
                print(f"Warning: Could not get samples from group '{group}': {e}")
        
        # Remove duplicates while preserving order
        seen = set()
        all_sample_ids = [sid for sid in all_sample_ids if not (sid in seen or seen.add(sid))]
        
        if not all_sample_ids:
            print("No samples found in selected groups.")
            return None
        
        print(f"\n✓ Selected {len(selected_groups)} group(s) with {len(all_sample_ids)} total samples")
        
        # Step 0.5: Gate mode selection (quadrant vs polygon)
        print("\n" + "="*60)
        print("Gate Mode Selection")
        print("="*60)
        print("Choose the type of gates to visualize:")
        print("  1. Polygon gates (default)")
        print("  2. Quadrant gates")
        
        gate_mode = 'polygon'  # Default
        while True:
            try:
                mode_choice = input(f"\nSelect gate mode (1-2, default=1): ").strip()
                if not mode_choice:
                    mode_choice = "1"
                if mode_choice == "1":
                    gate_mode = 'polygon'
                    print("✓ Selected: Polygon gates")
                    break
                elif mode_choice == "2":
                    gate_mode = 'quadrant'
                    print("✓ Selected: Quadrant gates")
                    break
                else:
                    print("Please enter 1 for Polygon gates or 2 for Quadrant gates")
            except KeyboardInterrupt:
                return None
        
        # Step 0.6: Ask how many plots to generate
        print("\n" + "="*60)
        print("Multiple Plots Configuration")
        print("="*60)
        while True:
            try:
                num_plots_input = input("\nHow many plots do you want to generate? (default=1): ").strip()
                if not num_plots_input:
                    num_plots = 1
                else:
                    num_plots = int(num_plots_input)
                if num_plots > 0:
                    break
                else:
                    print("Please enter a positive number")
            except ValueError:
                print("Please enter a valid number")
        
        # Step 0.6: Ask about keyword/statistic display
        print("\n" + "="*60)
        print("Additional Information Display")
        print("="*60)
        print("You can display additional information (statistics or keywords) in the bottom right of each plot.")
        show_keywords = False
        selected_keywords_to_show = None
        show_statistics = False
        
        while True:
            info_choice = input("\nDisplay additional information? (1=Statistics, 2=Keywords, 3=None, default=3): ").strip()
            if not info_choice:
                info_choice = "3"
            
            if info_choice == "1":
                # Statistics option
                show_keywords = True  # Reuse this flag for displaying info
                show_statistics = True
                print("\nAvailable statistics:")
                print("  1. Event count only")
                print("  2. Median")
                print("  3. Mean")
                print("  4. Median + Mean")
                print("  5. Event count + Median")
                print("  6. Event count + Mean")
                print("  7. Event count + Median + Mean")
                
                while True:
                    stat_display_choice = input("\nSelect statistics to display (1-7): ").strip()
                    if stat_display_choice == "1":
                        selected_keywords_to_show = ["count"]
                        break
                    elif stat_display_choice == "2":
                        selected_keywords_to_show = ["median"]
                        break
                    elif stat_display_choice == "3":
                        selected_keywords_to_show = ["mean"]
                        break
                    elif stat_display_choice == "4":
                        selected_keywords_to_show = ["median", "mean"]
                        break
                    elif stat_display_choice == "5":
                        selected_keywords_to_show = ["count", "median"]
                        break
                    elif stat_display_choice == "6":
                        selected_keywords_to_show = ["count", "mean"]
                        break
                    elif stat_display_choice == "7":
                        selected_keywords_to_show = ["count", "median", "mean"]
                        break
                    else:
                        print("Please enter a number between 1 and 7")
                break
            elif info_choice == "2":
                # Keywords option
                show_keywords = True
                show_statistics = False
                # Prompt user to select which keywords to display
                print("\nSelect keywords to display on plots:")
                # Get keywords from first sample
                first_sample_id = all_sample_ids[0] if all_sample_ids else None
                if first_sample_id:
                    try:
                        sample_keywords = self.workspace.get_keywords(first_sample_id)
                        if sample_keywords:
                            keyword_list = list(sample_keywords.keys())
                            print("\nAvailable keywords:")
                            for i, key in enumerate(keyword_list, 1):
                                print(f"  {i}. {key}")
                            
                            selected_keywords_to_show = []
                            print("\nEnter keyword numbers to display (comma-separated, e.g., '1,2,3' or 'all' for all):")
                            while True:
                                try:
                                    keyword_selection = input("Selection: ").strip()
                                    if keyword_selection.lower() == 'all':
                                        selected_keywords_to_show = keyword_list
                                        break
                                    
                                    # Parse comma-separated numbers
                                    indices = [int(x.strip()) - 1 for x in keyword_selection.split(',')]
                                    valid_indices = [i for i in indices if 0 <= i < len(keyword_list)]
                                    if valid_indices:
                                        selected_keywords_to_show = [keyword_list[i] for i in valid_indices]
                                        print(f"Selected keywords: {', '.join(selected_keywords_to_show)}")
                                        break
                                    else:
                                        print(f"Please enter valid numbers between 1 and {len(keyword_list)}")
                                except ValueError:
                                    print("Please enter comma-separated numbers (e.g., '1,2,3') or 'all'")
                                except KeyboardInterrupt:
                                    return None
                        else:
                            print("No keywords found in samples. Keywords will not be displayed.")
                            selected_keywords_to_show = None
                            show_keywords = False
                    except Exception as e:
                        print(f"Warning: Could not get keywords: {e}. Keywords will not be displayed.")
                        selected_keywords_to_show = None
                        show_keywords = False
                else:
                    selected_keywords_to_show = None
                    show_keywords = False
                break
            elif info_choice == "3":
                # None option
                show_keywords = False
                selected_keywords_to_show = None
                show_statistics = False
                break
            else:
                print("Please enter 1 for Statistics, 2 for Keywords, or 3 for None")
        
        # Step 1: Optional keyword filtering
        keyword_filter = None
        filtered_sample_ids = all_sample_ids
        
        print("\n" + "="*60)
        print("Keyword Filtering (Optional)")
        print("="*60)
        print("Some workspaces have multiple files per well (e.g., different time points).")
        print("You can filter samples by a keyword before plotting.")
        
        while True:
            filter_choice = input("\nFilter by keyword? (yes/no): ").strip().lower()
            if filter_choice in ['yes', 'y']:
                # Get keywords from first sample (should be universal)
                first_sample_id = all_sample_ids[0]
                try:
                    keywords = self.workspace.get_keywords(first_sample_id)
                    
                    if not keywords:
                        print("No keywords found in samples. Skipping keyword filter.")
                        break
                    
                    # Display available keywords
                    print("\nAvailable keywords:")
                    keyword_list = list(keywords.items())
                    for i, (key, value) in enumerate(keyword_list, 1):
                        print(f"  {i}. {key}: {value}")
                    
                    # Prompt for keyword selection
                    while True:
                        try:
                            key_choice = input(f"\nSelect keyword (1-{len(keyword_list)}) or type name: ").strip()
                            
                            # Try as number first
                            try:
                                key_idx = int(key_choice) - 1
                                if 0 <= key_idx < len(keyword_list):
                                    selected_key = keyword_list[key_idx][0]
                                    break
                                else:
                                    print(f"Please enter a number between 1 and {len(keyword_list)}")
                            except ValueError:
                                # Try as keyword name (case-insensitive partial match)
                                selected_key = None
                                for key, _ in keyword_list:
                                    if key_choice.lower() in key.lower() or key.lower() in key_choice.lower():
                                        selected_key = key
                                        break
                                
                                if selected_key:
                                    break
                                else:
                                    print(f"Keyword '{key_choice}' not found. Please try again.")
                        except KeyboardInterrupt:
                            return None
                    
                    # Get unique values for this keyword across all samples
                    # Store both original and normalized (stripped) versions
                    unique_values_dict = {}  # normalized -> original
                    for sid in all_sample_ids:
                        try:
                            kw = self.workspace.get_keywords(sid)
                            if selected_key in kw:
                                orig_value = str(kw[selected_key])
                                normalized = orig_value.strip()
                                unique_values_dict[normalized] = orig_value
                        except:
                            pass
                    
                    unique_values_normalized = sorted(list(unique_values_dict.keys()))
                    # Display with quotes around values
                    values_display = [f"'{v}'" for v in unique_values_normalized]
                    print(f"\nAvailable values for '{selected_key}': {', '.join(values_display)}")
                    
                    # Prompt for value
                    while True:
                        filter_value_input = input(f"\nEnter keyword value to filter (or 'all' for no filter): ").strip()
                        if filter_value_input.lower() == 'all':
                            break
                        
                        # Normalize input (strip whitespace) for comparison
                        filter_value_normalized = filter_value_input.strip()
                        
                        # Check if normalized value matches any unique value
                        if filter_value_normalized in unique_values_normalized:
                            # Use the normalized value for filtering
                            keyword_filter = {'key': selected_key, 'value': filter_value_normalized}
                            
                            # Filter sample IDs (compare with normalized values)
                            filtered_sample_ids = []
                            for sid in all_sample_ids:
                                try:
                                    kw = self.workspace.get_keywords(sid)
                                    if selected_key in kw:
                                        sample_value_normalized = str(kw[selected_key]).strip()
                                        if sample_value_normalized == filter_value_normalized:
                                            filtered_sample_ids.append(sid)
                                except:
                                    pass
                            
                            if not filtered_sample_ids:
                                print(f"Warning: No samples found with {selected_key} = '{filter_value_normalized}'")
                                print("Try a different value or 'all' to skip filtering.")
                                continue
                            
                            print(f"Filtered to {len(filtered_sample_ids)} samples with {selected_key} = '{filter_value_normalized}'")
                            break
                        else:
                            print(f"Value '{filter_value_input}' not found. Available: {', '.join(values_display)}")
                            print("Or type 'all' to skip filtering.")
                    
                    break
                except Exception as e:
                    print(f"Error getting keywords: {e}")
                    print("Skipping keyword filter.")
                    break
            elif filter_choice in ['no', 'n']:
                break
            else:
                print("Please enter 'yes' or 'no'")
        
        # Use filtered sample IDs for rest of the process
        sample_ids = filtered_sample_ids if filtered_sample_ids else all_sample_ids
        
        if not sample_ids:
            print("No samples available after filtering.")
            return None
        
        # Step 2: Interactive well ID source selection
        well_id_source = 'auto'
        well_id_keyword = None
        
        print("\n" + "="*60)
        print("Well ID Extraction (Optional)")
        print("="*60)
        print("Choose where to extract well ID information from:")
        print("1. Auto (try keywords first, then filename)")
        print("2. From keyword")
        print("3. From filename only")
        
        while True:
            source_choice = input("\nSelect well ID source (1-3, default=1): ").strip()
            if not source_choice:
                source_choice = "1"
            
            if source_choice == "1":
                well_id_source = 'auto'
                break
            elif source_choice == "2":
                well_id_source = 'keyword'
                # Get keywords from first sample
                first_sample_id = sample_ids[0]
                try:
                    keywords = self.workspace.get_keywords(first_sample_id)
                    
                    if not keywords:
                        print("No keywords found. Falling back to auto mode.")
                        well_id_source = 'auto'
                        break
                    
                    # Display available keywords
                    print("\nAvailable keywords:")
                    keyword_list = list(keywords.items())
                    for i, (key, value) in enumerate(keyword_list, 1):
                        print(f"  {i}. {key}: {value}")
                    
                    # Prompt for keyword selection
                    while True:
                        try:
                            key_choice = input(f"\nSelect keyword for well ID (1-{len(keyword_list)}) or type name: ").strip()
                            
                            # Try as number first
                            try:
                                key_idx = int(key_choice) - 1
                                if 0 <= key_idx < len(keyword_list):
                                    well_id_keyword = keyword_list[key_idx][0]
                                    break
                                else:
                                    print(f"Please enter a number between 1 and {len(keyword_list)}")
                            except ValueError:
                                # Try as keyword name (case-insensitive partial match)
                                selected_key = None
                                for key, _ in keyword_list:
                                    if key_choice.lower() in key.lower() or key.lower() in key_choice.lower():
                                        selected_key = key
                                        break
                                
                                if selected_key:
                                    well_id_keyword = selected_key
                                    break
                                else:
                                    print(f"Keyword '{key_choice}' not found. Please try again.")
                        except KeyboardInterrupt:
                            return None
                    
                    # Validate: Test the selected keyword on all samples and list matched well IDs
                    print(f"\nValidating well ID extraction using keyword '{well_id_keyword}'...")
                    matched_well_ids = []
                    samples_with_well_ids = []
                    
                    for sid in sample_ids:
                        try:
                            sample = self.workspace.get_sample(sid)
                            # Pass sample_id to parse_well_id so it can use workspace.get_keywords()
                            r, c, method = self.parse_well_id(sample, source='keyword', keyword_name=well_id_keyword, return_method=True, sample_id=sid)
                            if r and c:
                                well_id = f"{r}{c:02d}"
                                if well_id not in matched_well_ids:
                                    matched_well_ids.append(well_id)
                                samples_with_well_ids.append((sid, well_id))
                        except Exception as e:
                            # Debug: print error for first few failures to help diagnose
                            if len(samples_with_well_ids) < 3:
                                print(f"    Debug: Failed to parse well ID for {sid}: {e}")
                            pass
                    
                    # Display results
                    if matched_well_ids:
                        matched_well_ids.sort()  # Sort for easier reading
                        print(f"\n✓ Found {len(matched_well_ids)} unique well IDs from {len(samples_with_well_ids)} samples:")
                        # Display in rows of 12 (like a 96-well plate)
                        for i in range(0, len(matched_well_ids), 12):
                            row = matched_well_ids[i:i+12]
                            print(f"  {', '.join(row)}")
                        print()
                    else:
                        print(f"\n✗ ERROR: No well IDs could be extracted using keyword '{well_id_keyword}'")
                        print("Please check that:")
                        print("  1. The keyword name is correct")
                        print("  2. The keyword values contain well IDs in format like 'A1', 'A01', 'B12', etc.")
                        print("\nStopping. Please try again with a different keyword or use 'auto' mode.")
                        return None
                    
                    break
                except Exception as e:
                    print(f"Error getting keywords: {e}")
                    print("Falling back to auto mode.")
                    well_id_source = 'auto'
                    break
            elif source_choice == "3":
                well_id_source = 'filename'
                break
            else:
                print("Please enter 1, 2, or 3")
        
        # Collect plot configurations
        plot_configs = []
        
        for plot_num in range(1, num_plots + 1):
            if num_plots > 1:
                print("\n" + "="*60)
                print(f"Plot Configuration {plot_num} of {num_plots}")
                print("="*60)
            
            # Use first sample to get gate structure
            first_sample_id = sample_ids[0]
            gates_by_path = self._get_all_gates_info(first_sample_id)
            
            # Display available gate paths
            print("\nAvailable Gate Paths:")
            path_list = sorted(gates_by_path.keys())
            for i, path in enumerate(path_list, 1):
                gate_names = [g[0] for g in gates_by_path[path]]
                print(f"{i}. {path} (gates: {', '.join(gate_names)})")
            
            # Prompt for gate path selection
            while True:
                try:
                    path_choice = input(f"\nSelect gate path (1-{len(path_list)}): ").strip()
                    path_idx = int(path_choice) - 1
                    if 0 <= path_idx < len(path_list):
                        selected_path_str = path_list[path_idx]
                        break
                    else:
                        print(f"Please enter a number between 1 and {len(path_list)}")
                except ValueError:
                    print("Please enter a valid number")
            
            # Get gate names for selected path
            gate_options = gates_by_path[selected_path_str]
            print(f"\nAvailable gates in path '{selected_path_str}':")
            
            # Check if this path has a parent (not root)
            use_parent_option = False
            parent_gate_name = None
            parent_gate_path = None
            
            if selected_path_str != "root":
                # Find the parent gate by looking at the path structure
                # e.g., "root → Cells → Singlets" has parent "Singlets" at "root → Cells"
                path_parts = selected_path_str.split(" → ")
                if len(path_parts) > 1:
                    # Parent path is everything except the last part
                    parent_path_str = " → ".join(path_parts[:-1])
                    if parent_path_str in gates_by_path:
                        parent_gates = gates_by_path[parent_path_str]
                        # The parent gate is the last gate in the parent path
                        if parent_gates:
                            parent_gate_name, parent_gate_path = parent_gates[-1]
                            use_parent_option = True
            
            # For quadrant mode, skip individual gate selection (just use parent gate for data filtering)
            selected_gate_name = None
            selected_gate_path = None
            use_gate_data_directly = False
            
            if gate_mode == 'quadrant':
                # In quadrant mode, we only need the parent gate for data filtering
                # No need to select individual Q1-Q4 gates
                if use_parent_option:
                    selected_gate_name = parent_gate_name
                    selected_gate_path = parent_gate_path
                    print(f"✓ Using parent gate '{parent_gate_name}' for data filtering (quadrant mode)")
                else:
                    # If no parent, use the first gate in the path
                    if gate_options:
                        selected_gate_name, selected_gate_path = gate_options[0]
                        print(f"✓ Using gate '{selected_gate_name}' for data filtering (quadrant mode)")
                    else:
                        print("✗ No gates available in selected path")
                        continue
            else:
                # Polygon mode: show gate selection as before
                # Display option to use parent gate if available
                if use_parent_option:
                    print(f"0. Use parent gate '{parent_gate_name}' (no further gating)")
                
                for i, (gate_name, gate_path) in enumerate(gate_options, 1):
                    print(f"{i}. {gate_name}")
                
                # Prompt for gate name selection
                while True:
                    try:
                        max_option = len(gate_options)
                        if use_parent_option:
                            prompt_text = f"\nSelect gate name (0 to use parent '{parent_gate_name}', 1-{max_option}): "
                        else:
                            prompt_text = f"\nSelect gate name (1-{max_option}): "
                        
                        gate_choice = input(prompt_text).strip()
                        gate_idx = int(gate_choice)
                        
                        if use_parent_option and gate_idx == 0:
                            # User wants to use parent gate
                            selected_gate_name = parent_gate_name
                            selected_gate_path = parent_gate_path
                            print(f"✓ Selected parent gate '{parent_gate_name}' for data filtering (no further gating)")
                            break
                        elif 1 <= gate_idx <= max_option:
                            selected_gate_name, selected_gate_path = gate_options[gate_idx - 1]
                            break
                        else:
                            if use_parent_option:
                                print(f"Please enter 0 to use parent gate, or a number between 1 and {max_option}")
                            else:
                                print(f"Please enter a number between 1 and {max_option}")
                    except ValueError:
                        print("Please enter a valid number")
            
            # Get available channels from first sample
            first_sample = self.workspace.get_sample(first_sample_id)
            available_channels = first_sample.pnn_labels
            
            print(f"\nAvailable channels:")
            for i, channel in enumerate(available_channels, 1):
                print(f"  {i}. {channel}")
            
            # Prompt for plot type
            print("\nPlot types:")
            print("1. Histogram (requires 1 parameter)")
            print("2. Scatter (requires 2 parameters: x and y)")
            print("3. Contour (density plot, requires 2 parameters: x and y)")

            while True:
                plot_choice = input("\nSelect plot type (1, 2, or 3): ").strip()
                if plot_choice == "1":
                    plot_type = "histogram"
                    break
                elif plot_choice == "2":
                    plot_type = "scatter"
                    break
                elif plot_choice == "3":
                    plot_type = "contour"
                    break
                else:
                    print("Please enter 1 for Histogram, 2 for Scatter, or 3 for Contour")
            
            # Prompt for parameters
            parameters = []
            statistic_choice = None
            scale_choice = "linear"  # Default scale
            plot_config_show_gates = False  # Initialize, will be set for scatter plots
            if plot_type == "histogram":
                while True:
                    try:
                        param_choice = input(f"\nSelect channel for histogram (1-{len(available_channels)}) or type name: ").strip()
                        
                        # Try as number first
                        try:
                            param_idx = int(param_choice) - 1
                            if 0 <= param_idx < len(available_channels):
                                param = available_channels[param_idx]
                                parameters.append(param)
                                break
                            else:
                                print(f"Please enter a number between 1 and {len(available_channels)}")
                        except ValueError:
                            # Try as channel name (case-insensitive partial match)
                            selected_channel = None
                            for channel in available_channels:
                                if param_choice.lower() in channel.lower() or channel.lower() in param_choice.lower():
                                    selected_channel = channel
                                    break
                            
                            if selected_channel:
                                parameters.append(selected_channel)
                                break
                            else:
                                print(f"Channel '{param_choice}' not found. Please try again.")
                    except KeyboardInterrupt:
                        return None
                
                # Prompt for statistic to display
                print("\nStatistic to display:")
                print("1. Median")
                print("2. Mean")
                while True:
                    stat_choice = input("\nSelect statistic (1 or 2): ").strip()
                    if stat_choice == "1":
                        statistic_choice = "median"
                        break
                    elif stat_choice == "2":
                        statistic_choice = "mean"
                        break
                    else:
                        print("Please enter 1 for Median or 2 for Mean")
                
                # Prompt for scale (linear or log)
                print("\nX-axis scale:")
                print("1. Linear")
                print("2. Log")
                while True:
                    scale_choice = input("\nSelect scale (1 or 2, default=1): ").strip()
                    if not scale_choice:
                        scale_choice = "1"
                    if scale_choice == "1":
                        scale_choice = "linear"
                        break
                    elif scale_choice == "2":
                        scale_choice = "log"
                        break
                    else:
                        print("Please enter 1 for Linear or 2 for Log")
            else:  # scatter or contour (both require 2 parameters)
                while True:
                    try:
                        x_choice = input(f"\nSelect channel for X-axis (1-{len(available_channels)}) or type name: ").strip()
                        
                        # Try as number first
                        try:
                            x_idx = int(x_choice) - 1
                            if 0 <= x_idx < len(available_channels):
                                x_param = available_channels[x_idx]
                                parameters.append(x_param)
                                break
                            else:
                                print(f"Please enter a number between 1 and {len(available_channels)}")
                        except ValueError:
                            # Try as channel name (case-insensitive partial match)
                            selected_channel = None
                            for channel in available_channels:
                                if x_choice.lower() in channel.lower() or channel.lower() in x_choice.lower():
                                    selected_channel = channel
                                    break
                            
                            if selected_channel:
                                parameters.append(selected_channel)
                                break
                            else:
                                print(f"Channel '{x_choice}' not found. Please try again.")
                    except KeyboardInterrupt:
                        return None
                
                while True:
                    try:
                        y_choice = input(f"Select channel for Y-axis (1-{len(available_channels)}) or type name: ").strip()
                        
                        # Try as number first
                        try:
                            y_idx = int(y_choice) - 1
                            if 0 <= y_idx < len(available_channels):
                                y_param = available_channels[y_idx]
                                parameters.append(y_param)
                                break
                            else:
                                print(f"Please enter a number between 1 and {len(available_channels)}")
                        except ValueError:
                            # Try as channel name (case-insensitive partial match)
                            selected_channel = None
                            for channel in available_channels:
                                if y_choice.lower() in channel.lower() or channel.lower() in y_choice.lower():
                                    selected_channel = channel
                                    break
                            
                            if selected_channel:
                                parameters.append(selected_channel)
                                break
                            else:
                                print(f"Channel '{y_choice}' not found. Please try again.")
                    except KeyboardInterrupt:
                        return None
                
                # For quadrant mode, extract thresholds from parent QuadrantGate or infer from Q1-Q4
                quadrant_thresholds = None
                if gate_mode == 'quadrant':
                    x_channel = parameters[0]
                    y_channel = parameters[1]
                    # Convert selected_path_str to tuple format
                    path_parts = selected_path_str.split(" → ")
                    gate_path_tuple = tuple([part.strip() for part in path_parts])
                    
                    # First try to extract from parent QuadrantGate
                    quadrant_thresholds = self._extract_quadrant_thresholds_from_parent(
                        first_sample_id, gate_path_tuple, x_channel, y_channel
                    )
                    
                    # If not found, try to infer from Q1-Q4 gates
                    if not quadrant_thresholds:
                        print(f"\n⚠️  No parent QuadrantGate found. Attempting to infer thresholds from Q1-Q4 gates...")
                        inferred_dividers = self._infer_quadrant_dividers_from_regions(
                            first_sample_id, gate_path_tuple, x_channel, y_channel
                        )
                        if inferred_dividers:
                            # Convert dividers to thresholds format
                            x_threshold = None
                            y_threshold = None
                            for divider in inferred_dividers:
                                if divider.get('orientation') == 'vertical':
                                    x_threshold = divider.get('value')
                                elif divider.get('orientation') == 'horizontal':
                                    y_threshold = divider.get('value')
                            
                            if x_threshold is not None and y_threshold is not None:
                                quadrant_thresholds = {
                                    'x_threshold': float(x_threshold),
                                    'y_threshold': float(y_threshold)
                                }
                                print(f"✓ Inferred quadrant thresholds from Q1-Q4 gates:")
                                print(f"   X threshold: {quadrant_thresholds['x_threshold']:.2f}")
                                print(f"   Y threshold: {quadrant_thresholds['y_threshold']:.2f}")
                    
                    if quadrant_thresholds:
                        print(f"\n✓ Using quadrant thresholds:")
                        print(f"   X threshold: {quadrant_thresholds['x_threshold']:.2f}")
                        print(f"   Y threshold: {quadrant_thresholds['y_threshold']:.2f}")
                    else:
                        print(f"\n⚠️  Could not extract or infer quadrant thresholds. Switching to polygon mode.")
                        gate_mode = 'polygon'
                
                # Ask about gate visualization for scatter plots (only for polygon mode)
                # Initialize these variables before the if statement to avoid UnboundLocalError
                available_gates = []
                available_gates_with_info = []

                if gate_mode == 'polygon':
                    print("\n" + "="*60)
                    print("Gate Visualization")
                    print("="*60)
                    print("You can visualize gate boundaries on scatter plots.")

                    # Check what gates are available for the selected channels
                    x_channel = parameters[0]
                    y_channel = parameters[1]
                    try:
                        # Get gates from first sample to see what's available
                        gates_data = self._extract_gate_polygons(first_sample_id, x_channel, y_channel)
                        available_gates = [g['name'] for g in gates_data]
                        available_gates_with_info = gates_data
                        
                        # Also check for child gates and sibling gates of the selected gate
                        if selected_gate_name and selected_gate_name != "Ungated":
                            # Get all gate IDs to find children and siblings
                            all_gate_ids = self.workspace.get_gate_ids(first_sample_id)
                            # The full path to the selected gate (for finding siblings)
                            selected_gate_full_path = selected_gate_path + (selected_gate_name,)
                            # The path prefix for finding children
                            child_gate_path_prefix = selected_gate_full_path
                            
                            print(f"    DEBUG: Looking for gates related to '{selected_gate_name}'")
                            print(f"           Selected gate path: {selected_gate_path}")
                            print(f"           Selected gate full path: {selected_gate_full_path}")
                            
                            for other_gate_name, other_gate_path in all_gate_ids:
                                # Skip if this is the selected gate itself
                                if other_gate_name == selected_gate_name and other_gate_path == selected_gate_path:
                                    continue
                                
                                is_child_or_sibling = False
                                
                                # Check if this is a child of the selected gate
                                # Child path should start with selected gate's full path and be longer
                                if len(other_gate_path) > len(selected_gate_full_path):
                                    if other_gate_path[:len(selected_gate_full_path)] == selected_gate_full_path:
                                        is_child_or_sibling = True
                                        print(f"    DEBUG: Found child gate '{other_gate_name}' at path {other_gate_path}")
                                
                                # Also check if this is a sibling at the same path level as the selected gate
                                # (e.g., Q1-Q4 are siblings of Singlets, all at "root → Cells → Singlets")
                                # Siblings have the same path as the selected gate's full path
                                if other_gate_path == selected_gate_full_path:
                                    is_child_or_sibling = True
                                    print(f"    DEBUG: Found sibling gate '{other_gate_name}' at same path {other_gate_path}")
                                
                                if is_child_or_sibling:
                                    # Try to extract gate if it matches channels
                                    other_gate_data = self._extract_selected_gate(first_sample_id, other_gate_name, other_gate_path, x_channel, y_channel)
                                    if other_gate_data and other_gate_data['name'] not in available_gates:
                                        available_gates.append(other_gate_data['name'])
                                        available_gates_with_info.append(other_gate_data)
                                        print(f"    DEBUG: Added '{other_gate_name}' to available gates (type: {other_gate_data.get('type', 'unknown')})")
                                    elif not other_gate_data:
                                        print(f"    DEBUG: Could not extract '{other_gate_name}' (doesn't match channels or extraction failed)")
                    except Exception as e:
                        pass  # Will show empty list
                
                # Store which gates to visualize
                gates_to_visualize = []

                if available_gates:
                    print(f"\nAvailable gates for channels {x_channel} / {y_channel}:")
                    for i, gate_info in enumerate(available_gates_with_info, 1):
                        gate_name = gate_info['name']
                        gate_type = gate_info.get('type', 'unknown')
                        if gate_type == 'quadrant':
                            num_dividers = len(gate_info.get('dividers', []))
                            print(f"  {i}. {gate_name} (quadrant, {num_dividers} divider(s))")
                        else:
                            num_vertices = len(gate_info.get('xs', []))
                            print(f"  {i}. {gate_name} ({gate_type}, {num_vertices} vertices)")

                    print("\nYou can visualize gate boundaries on this scatter plot.")
                    print("Quadrant gates will be shown as divider lines, polygon/rectangle gates as filled regions.")
                    print("This is useful to see gates overlaid on Ungated populations.")

                    while True:
                        gate_viz_choice = input("\nSelect gates to visualize (comma-separated numbers, 'all', or 'none'): ").strip().lower()

                        if gate_viz_choice in ['none', 'no', 'n', '']:
                            plot_config_show_gates = False
                            break
                        elif gate_viz_choice in ['all', 'yes', 'y']:
                            plot_config_show_gates = True
                            gates_to_visualize = available_gates
                            print(f"Will visualize: {', '.join(gates_to_visualize)}")
                            break
                        else:
                            # Parse comma-separated numbers
                            try:
                                indices = [int(x.strip()) - 1 for x in gate_viz_choice.split(',')]
                                selected_gates = []
                                for idx in indices:
                                    if 0 <= idx < len(available_gates):
                                        selected_gates.append(available_gates[idx])

                                if selected_gates:
                                    plot_config_show_gates = True
                                    gates_to_visualize = selected_gates
                                    print(f"Will visualize: {', '.join(gates_to_visualize)}")
                                    break
                                else:
                                    print(f"Please enter valid numbers between 1 and {len(available_gates)}")
                            except ValueError:
                                print("Please enter comma-separated numbers (e.g., '1,2'), 'all', or 'none'")
                else:
                    print(f"\nNo gates found for channels {x_channel} / {y_channel}.")
                    print("Gate visualization will not be available for this plot configuration.")
                    plot_config_show_gates = False
            
            # Create plot configuration
            # Check if user selected "use parent gate" option
            use_parent_for_data = use_parent_option and (selected_gate_name == parent_gate_name)
            
            plot_config = {
                'gate_path': selected_gate_path,
                'gate_name': selected_gate_name,
                'plot_type': plot_type,
                'parameters': parameters,
                'use_gate_data_directly': use_parent_for_data,  # Flag to use gate's data directly, not parent's
                'gate_mode': gate_mode
            }
            
            # Add quadrant thresholds if in quadrant mode
            if gate_mode == 'quadrant' and quadrant_thresholds:
                plot_config['quadrant_thresholds'] = quadrant_thresholds
            
            # Add statistic and scale choices for histograms
            if plot_type == "histogram":
                if statistic_choice:
                    plot_config['statistic'] = statistic_choice
                if scale_choice:
                    plot_config['scale'] = scale_choice
                plot_config['show_gates'] = False  # Gates not applicable for histograms
            else:  # scatter
                # Use explicit check for scatter plot type instead of unreliable locals() check
                plot_config['show_gates'] = plot_config_show_gates
                # Store which specific gates to visualize
                if 'gates_to_visualize' in locals() and gates_to_visualize:
                    plot_config['gates_to_visualize'] = gates_to_visualize
            
            plot_configs.append(plot_config)
        
        # Combine base config with plot configs
        result = {
            'selected_groups': selected_groups,
            'filtered_sample_ids': sample_ids,
            'well_id_source': well_id_source,
            'well_id_keyword': well_id_keyword,
            'show_keywords': show_keywords,
            'selected_keywords_to_show': selected_keywords_to_show if 'selected_keywords_to_show' in locals() else None,
            'show_statistics': show_statistics if 'show_statistics' in locals() else False,
            'num_plots': num_plots,
            'plot_configs': plot_configs
        }
        
        # Add keyword filter if applied
        if keyword_filter:
            result['keyword_filter'] = keyword_filter
        
        return result

    def _plot_histogram(self, df, channel_keyword, sample_id, gate_name, bins=200, statistic='median', scale='linear', keywords=None, show_keywords=False):
        """
        Generate a histogram plot for a single channel using raw (untransformed) data.
        
        Creates a Bokeh histogram visualization showing the frequency distribution of events
        for a selected flow cytometry parameter. Supports both linear and log scales with
        automatic outlier filtering for better visualization.
        
        Args:
            df (pd.DataFrame): DataFrame containing gate events with raw (untransformed) data.
                Must have columns matching the channel_keyword. Typically obtained from
                workspace.get_gate_events(..., source="raw").
            channel_keyword (str): Keyword substring to identify the channel column.
                Uses substring matching, so "RL1-A" will match "RL1-A" or "RL1-A pHrodo-APC-A".
                Examples: "RL1-A", "pHrodo", "FSC-A", "SSC-A", "B1-A".
            sample_id (str): Sample identifier to include in plot title. Can be sample ID or
                well position (e.g., "A01", "Sample001").
            gate_name (str): Gate name to include in plot title (e.g., "Singlets", "Cells").
            bins (int, optional): Number of bins for histogram calculation. Defaults to 200.
            statistic (str, optional): Statistic to display as a vertical line on the plot.
                Options: 'median' or 'mean'. Defaults to 'median'.
            scale (str, optional): X-axis scale type. Options: 'linear' or 'log'.
                - Linear: Dynamic range based on data with outlier filtering
                - Log: Fixed range from 0.1 to 10^5 (log scale cannot start at 0)
                Defaults to 'linear'.
        
        Returns:
            bokeh.plotting.figure: Bokeh figure object with histogram rendered. The figure
                includes:
                - Histogram bars (skyblue fill, white borders)
                - Vertical dashed line indicating the selected statistic (median/mean)
                - Text label showing statistic name and value (rounded to 1 decimal)
                - Appropriate x-axis scale (linear or log)
                - Standard Bokeh tools (pan, zoom, reset, save)
        
        Raises:
            ValueError: If:
                - No channel column matches the keyword
                - No valid data found for the channel
                - Log scale selected but no positive values exist
        
        Example:
            >>> # Create histogram with linear scale and median
            >>> df = workspace.get_gate_events("A1.fcs", "Singlets", source="raw")
            >>> p = analyzer._plot_histogram(df, "RL1-A", "A01", "Singlets", 
            ...                              statistic='median', scale='linear')
            
            >>> # Create histogram with log scale and mean
            >>> p = analyzer._plot_histogram(df, "B1-A", "A01", "Singlets",
            ...                              statistic='mean', scale='log')
        
        Notes:
            - Uses raw (untransformed) data for accurate representation
            - Applies IQR-based outlier filtering (3×IQR) to improve visualization
            - For log scale, automatically filters out zero and negative values
            - Statistic calculation uses full data (before outlier filtering)
            - Histogram uses filtered data (after outlier removal) for better display
            - Log scale has fixed x-axis range (0.1 to 10^5) for consistency across samples
        """
        from bokeh.plotting import figure
        
        # Find matching channel
        matching_channels = [c for c in df.columns if channel_keyword in c]
        if not matching_channels:
            raise ValueError(f"No channel found matching keyword '{channel_keyword}'")
        channel = matching_channels[0]
        
        # Create histogram
        # Filter out NaN/inf values
        valid_data = df[channel].dropna()
        if len(valid_data) == 0:
            raise ValueError(f"No valid data for channel '{channel}'")
        
        # For log scale, filter out zero and negative values first
        if scale == 'log':
            valid_data = valid_data[valid_data > 0]
            if len(valid_data) == 0:
                raise ValueError(f"No positive values for channel '{channel}' (required for log scale)")
        
        # Filter out small values (<10) to prevent log transform issues and giant bars
        # This is especially important for log scale but also helps with linear scale
        valid_data = valid_data[valid_data >= 10]
        if len(valid_data) == 0:
            raise ValueError(f"No values >= 10 for channel '{channel}' (filtered out small values)")
        
        # Filter outliers using IQR method for better visualization
        # This helps when data has extreme outliers that compress the main distribution
        q1 = valid_data.quantile(0.25)
        q3 = valid_data.quantile(0.75)
        iqr = q3 - q1
        
        # Use a more conservative approach: filter only extreme outliers (beyond 3*IQR)
        # This preserves most of the data while removing extreme outliers
        lower_bound = q1 - 3 * iqr
        upper_bound = q3 + 3 * iqr
        
        # For log scale, ensure lower bound is positive
        if scale == 'log' and lower_bound <= 0:
            lower_bound = valid_data.min()
        
        # Filter data for plotting (but keep original for statistics)
        filtered_data = valid_data[(valid_data >= lower_bound) & (valid_data <= upper_bound)]
        
        # If filtering removed too much (>50%), use original data
        # This handles cases where the distribution itself is very wide
        if len(filtered_data) < len(valid_data) * 0.5:
            filtered_data = valid_data
        
        # Compute KDE (kernel density estimate) on filtered data
        # For log scale, compute KDE on log-transformed data for better density estimation
        if scale == "log":
            # Transform to log space for KDE
            log_data = np.log10(filtered_data)
            kde = gaussian_kde(log_data)
            # Generate x grid in log space (from log10(0.1) = -1 to log10(1e5) = 5)
            x_min_log = -1
            x_max_log = 5
            x_grid_log = np.linspace(x_min_log, x_max_log, 500)
            # Evaluate KDE on log grid
            kde_values = kde(x_grid_log)
            # Transform grid back to linear space for plotting
            x_grid = 10 ** x_grid_log
            # Adjust density for log scale: multiply by derivative of log transform
            # d/dx log10(x) = 1/(x * ln(10)), so we need to multiply by x * ln(10)
            kde_values = kde_values * (x_grid * np.log(10))
        else:
            # For linear scale, compute KDE directly on filtered data
            kde = gaussian_kde(filtered_data)
            # Generate x grid over valid range with padding
            x_min = filtered_data.min()
            x_max = filtered_data.max()
            x_padding = (x_max - x_min) * 0.1  # 10% padding
            x_grid = np.linspace(max(0, x_min - x_padding), x_max + x_padding, 500)
            # Evaluate KDE on grid
            kde_values = kde(x_grid)
        
        # Set x-axis type based on scale parameter
        x_axis_type = "log" if scale == "log" else "linear"
        
        # For log scale, set fixed range from ~0 to 10^5
        # Note: Log scale can't start at exactly 0, so we use 0.1 (10^-1) as minimum
        if scale == "log":
            p = figure(title=f"{sample_id} - {gate_name}",
                       x_axis_label=channel, y_axis_label="Density",
                       x_axis_type=x_axis_type,
                       x_range=(0.1, 1e5),  # Fixed range: ~0 (0.1) to 10^5
                       width=400, height=300, tools="pan,wheel_zoom,box_zoom,reset,save",
                       sizing_mode="fixed")
        else:
            p = figure(title=f"{sample_id} - {gate_name}",
                       x_axis_label=channel, y_axis_label="Density",
                       x_axis_type=x_axis_type,
                       width=400, height=300, tools="pan,wheel_zoom,box_zoom,reset,save",
                       sizing_mode="fixed")
        
        # Initialize max_density for use in keyword display
        max_density = 1
        
        # Plot KDE as smooth silhouette with olive green fill
        if len(kde_values) > 0 and len(x_grid) > 0:
            # Create area under curve (filled silhouette)
            p.patch(np.concatenate([[x_grid[0]], x_grid, [x_grid[-1]]]),
                   np.concatenate([[0], kde_values, [0]]),
                   fill_color="#6B8E23", fill_alpha=0.6, line_color="#556B2F", line_width=2)
            
            # Add line on top for clearer silhouette
            p.line(x_grid, kde_values, line_color="#556B2F", line_width=2.5, alpha=0.9)
            
            # Calculate max density for positioning
            max_density = np.max(kde_values) if len(kde_values) > 0 else 1
            
            # Calculate and draw statistic line (median or mean)
            if statistic == 'mean':
                stat_val = valid_data.mean()
                stat_label = 'Mean'
            else:  # default to median
                stat_val = valid_data.median()
                stat_label = 'Median'
            
            from bokeh.models import Span, Label
            stat_line = Span(location=stat_val, dimension='height', 
                           line_color='gray', line_width=2, line_dash='dashed', line_alpha=0.7)
            p.add_layout(stat_line)
            
            # Add text label with value (rounded to 1 decimal)
            stat_text = f"{stat_label}: {stat_val:.1f}"
            # Position label at top of the line, slightly to the right
            stat_label_obj = Label(x=stat_val, y=max_density * 0.95, text=stat_text,
                                  text_font_size="9pt", text_color="gray",
                                  text_align="center", background_fill_color="white",
                                  background_fill_alpha=0.8, border_line_color="gray",
                                  border_line_width=1)
            p.add_layout(stat_label_obj)
            
            # Set ranges explicitly
            # For log scale, keep the fixed range (0 to 10^5)
            if scale != "log":
                # For linear scale, use the grid range
                p.x_range.start = max(0, x_grid[0])
                p.x_range.end = x_grid[-1]
            # For log scale, x_range is already set to (0.1, 1e5) in figure creation
            
            # For y-axis, add padding at top
            p.y_range.start = 0
            p.y_range.end = max_density * 1.1 if max_density > 0 else 1
        else:
            # If no valid data, add a text message
            p.x_range.start = 0
            p.x_range.end = 1
            p.y_range.start = 0
            p.y_range.end = 1
            p.text(x=[0.5], y=[0.5], text=["No valid data"], 
                   text_font_size="12pt", text_color="red",
                   text_align="center", text_baseline="middle")
        
        # Add keyword display in bottom right if requested
        if show_keywords and keywords:
            from bokeh.models import Label
            # Format keywords as text (one per line)
            # Keywords are already filtered by selected_keywords_to_show in generate_interactive_plots
            keyword_lines = []
            for key, value in keywords.items():
                keyword_lines.append(f"{key}: {value}")
            keyword_text = "\n".join(keyword_lines)
            
            # Get plot dimensions for positioning
            if scale == "log":
                x_min = 1
                x_max = 1e5
            else:
                x_min = x_grid[0] if len(x_grid) > 0 else 0
                x_max = x_grid[-1] if len(x_grid) > 0 else 1
            y_max = max_density * 1.1 if max_density > 0 else 1

            # Position in top left with padding
            keyword_label = Label(x=x_min + (x_max - x_min) * 0.02, y=y_max * 0.98,
                                 text=keyword_text,
                                 text_font_size="8pt", text_color="black",
                                 text_align="left", text_baseline="top",
                                 background_fill_color="white",
                                 background_fill_alpha=0.7, border_line_color="gray",
                                 border_line_width=1, x_units="data", y_units="data")
            p.add_layout(keyword_label)
        
        return p

    def _plot_scatter(self, df, x_keyword, y_keyword, sample_id, gate_name, show_gates=False, gates=None, selected_gate=None, show_keywords=False, keywords=None, gate_mode='polygon', quadrant_thresholds=None):
        """
        Generate a scatter plot for two channels using raw (untransformed) data.
        
        Creates a 2D scatter plot showing the relationship between two parameters.
        Automatically downsamples to 10,000 points if the dataset is larger for performance.
        
        Args:
            df (pd.DataFrame): DataFrame with gate events (raw data, source="raw")
            x_keyword (str): Keyword substring to identify X-axis channel column
                            (e.g., "FSC-A", "SSC-A", "RL1-A")
            y_keyword (str): Keyword substring to identify Y-axis channel column
                            (e.g., "FSC-A", "SSC-A", "VL1-A")
            sample_id (str): Sample ID to include in plot title
            gate_name (str): Gate name to include in plot title
            
        Returns:
            bokeh.plotting.figure: Bokeh figure object with scatter plot rendered
            
        Raises:
            ValueError: If no channel column matches either keyword
            
        Example:
            >>> p = analyzer._plot_scatter(df, "FSC-A", "SSC-A", "A1.fcs", "Cells")
            >>> # Creates FSC-A vs SSC-A scatter plot for Cells gate in sample A1.fcs
        """
        from bokeh.plotting import figure
        
        # Find matching channels
        x_channels = [c for c in df.columns if x_keyword in c]
        y_channels = [c for c in df.columns if y_keyword in c]
        
        if not x_channels:
            raise ValueError(f"No channel found matching keyword '{x_keyword}'")
        if not y_channels:
            raise ValueError(f"No channel found matching keyword '{y_keyword}'")
        
        x_channel = x_channels[0]
        y_channel = y_channels[0]
        
        # Filter out values <= 0 for log scale (log can't handle 0 or negative values)
        df = df[(df[x_channel] > 0) & (df[y_channel] > 0)]

        # Downsample if too many points for performance
        max_points = 10000
        if len(df) > max_points:
            df_plot = df.sample(n=max_points, random_state=42)
        else:
            df_plot = df

        # Convert to ColumnDataSource for Bokeh
        from bokeh.models import ColumnDataSource
        source = ColumnDataSource(df_plot)

        # Use log scale for both axes to match FlowJo display
        # Set fixed ranges: 1 to 10^5 (matching FlowJo's typical range)
        p = figure(title=f"{sample_id} - {gate_name}",
                   x_axis_label=x_channel, y_axis_label=y_channel,
                   x_axis_type="log", y_axis_type="log",
                   x_range=(1, 1e5), y_range=(1, 1e5),
                   width=400, height=300, tools="pan,wheel_zoom,box_zoom,reset,save",
                   sizing_mode="fixed")

        p.scatter(x=x_channel, y=y_channel, source=source, size=2, alpha=0.6, color="navy")

        # Render quadrant dividers if in quadrant mode
        if gate_mode == 'quadrant' and quadrant_thresholds:
            print(f"    → Rendering quadrant dividers on log scale (1 to 10^5)")
            from bokeh.models import Span
            
            # Use a distinct color for quadrant dividers
            divider_color = "#DC143C"  # Crimson
            
            # Render vertical divider (X threshold)
            if 'x_threshold' in quadrant_thresholds:
                x_thresh = quadrant_thresholds['x_threshold']
                span = Span(location=x_thresh, dimension='height',
                           line_color=divider_color, line_width=2.5,
                           line_dash='dashed', line_alpha=0.9)
                p.add_layout(span)
                print(f"      ✓ Rendered vertical divider at x={x_thresh:.2f}")
            
            # Render horizontal divider (Y threshold)
            if 'y_threshold' in quadrant_thresholds:
                y_thresh = quadrant_thresholds['y_threshold']
                span = Span(location=y_thresh, dimension='width',
                           line_color=divider_color, line_width=2.5,
                           line_dash='dashed', line_alpha=0.9)
                p.add_layout(span)
                print(f"      ✓ Rendered horizontal divider at y={y_thresh:.2f}")

        # Render gate boundaries if requested (polygon mode)
        elif show_gates and gates:
            print(f"    → Rendering {len(gates)} gate(s) on log scale (1 to 10^5)")

            # Color palette for gates (distinct colors/shades)
            gate_colors = [
                ("#6B8E23", "#556B2F"),  # Olive green
                ("#4682B4", "#36648B"),  # Steel blue
                ("#CD853F", "#8B5A2B"),  # Peru/tan
                ("#9370DB", "#7B68EE"),  # Medium purple
                ("#DC143C", "#B22222"),  # Crimson
                ("#20B2AA", "#008B8B"),  # Light sea green
                ("#FF8C00", "#CD6600"),  # Dark orange
                ("#8B4789", "#68387E"),  # Dark orchid
            ]

            from bokeh.models import Span

            # gates should be a list of gate dictionaries with 'name', 'x_dim', 'y_dim', 'xs', 'ys' (or 'dividers' for quadrant gates)
            for gate_idx, gate in enumerate(gates):
                # Check if gate matches current axes
                if (gate.get('x_dim') == x_channel and gate.get('y_dim') == y_channel) or \
                   (gate.get('x_dim') == y_channel and gate.get('y_dim') == x_channel):
                    
                    gate_type = gate.get('type', 'polygon')  # Default to polygon for backward compatibility
                    
                    # Handle quadrant gates
                    if gate_type == 'quadrant':
                        dividers = gate.get('dividers', [])
                        if dividers:
                            # Get color for this gate
                            fill_color, line_color = gate_colors[gate_idx % len(gate_colors)]
                            
                            # Render divider lines
                            for divider in dividers:
                                orientation = divider.get('orientation')
                                value = divider.get('value')
                                
                                if orientation == 'vertical' and value is not None:
                                    # Vertical line (x = value)
                                    span = Span(location=value, dimension='height',
                                               line_color=line_color, line_width=2.5,
                                               line_dash='dashed', line_alpha=0.9)
                                    p.add_layout(span)
                                    print(f"      ✓ Rendered vertical divider at x={value:.2f} for '{gate.get('name')}'")
                                elif orientation == 'horizontal' and value is not None:
                                    # Horizontal line (y = value)
                                    span = Span(location=value, dimension='width',
                                               line_color=line_color, line_width=2.5,
                                               line_dash='dashed', line_alpha=0.9)
                                    p.add_layout(span)
                                    print(f"      ✓ Rendered horizontal divider at y={value:.2f} for '{gate.get('name')}'")
                            
                            print(f"      ✓ Rendered quadrant gate '{gate.get('name')}' with {len(dividers)} divider(s), color: {line_color}")
                        else:
                            print(f"      ⚠️  Skipped '{gate.get('name')}': no dividers found")
                    
                    # Handle polygon/rectangle gates
                    else:
                        xs = gate.get('xs', [])
                        ys = gate.get('ys', [])
                        if xs and ys and len(xs) == len(ys) and len(xs) >= 3:
                            # Swap if dimensions are reversed
                            if gate.get('x_dim') == y_channel:
                                xs, ys = ys, xs

                            print(f"      DEBUG: Gate '{gate.get('name')}' coordinates:")
                            print(f"             X: [{min(xs):.2f}, {max(xs):.2f}]")
                            print(f"             Y: [{min(ys):.2f}, {max(ys):.2f}]")

                            # Ensure polygon is closed
                            if xs[0] != xs[-1] or ys[0] != ys[-1]:
                                xs = list(xs) + [xs[0]]
                                ys = list(ys) + [ys[0]]

                            # Get color for this gate
                            fill_color, line_color = gate_colors[gate_idx % len(gate_colors)]

                            # Plot gates with distinct colors (no legend_label - will add global legend)
                            # Gates are already in raw data space from inverse transform
                            # They will be displayed correctly on log-scale axes
                            p.patch(xs, ys, fill_color=fill_color, fill_alpha=0.3, line_color=line_color,
                                   line_width=2.5, line_alpha=0.9)
                            print(f"      ✓ Rendered '{gate.get('name')}' ({len(xs)} vertices, color: {fill_color})")
                        else:
                            print(f"      ⚠️  Skipped '{gate.get('name')}': invalid coordinates")
                else:
                    print(f"      ⚠️  Skipped '{gate.get('name')}': channel mismatch")
        
        # Add keyword display in bottom right if requested
        if show_keywords and keywords:
            from bokeh.models import Label
            # Format keywords as text (one per line)
            # Keywords are already filtered by selected_keywords_to_show in generate_interactive_plots
            keyword_lines = []
            for key, value in keywords.items():
                keyword_lines.append(f"{key}: {value}")
            keyword_text = "\n".join(keyword_lines)
            
            # Get plot dimensions for positioning
            x_data = df_plot[x_channel].dropna()
            y_data = df_plot[y_channel].dropna()
            if len(x_data) > 0 and len(y_data) > 0:
                x_min = x_data.min()
                x_max = x_data.max()
                y_max = y_data.max()

                # Position in top left with padding
                keyword_label = Label(x=x_min + (x_max - x_min) * 0.02, y=y_max * 0.98,
                                     text=keyword_text,
                                     text_font_size="8pt", text_color="black",
                                     text_align="left", text_baseline="top",
                                     background_fill_color="white",
                                     background_fill_alpha=0.7, border_line_color="gray",
                                     border_line_width=1, x_units="data", y_units="data")
                p.add_layout(keyword_label)
        
        return p
    def _plot_contour(self, df, x_keyword, y_keyword, sample_id, gate_name, show_gates=False, gates=None, selected_gate=None, show_keywords=False, keywords=None, gate_mode='polygon', quadrant_thresholds=None):
        """
        Generate a contour/density plot for two channels using raw (untransformed) data.

        Creates a 2D density plot with iso-density contours similar to traditional
        flow cytometry plots (e.g., FlowJo). Uses 2D Kernel Density Estimation (KDE)
        to calculate density and matplotlib for contour extraction, then renders
        contours in Bokeh for interactive visualization.

        The method works in log space to ensure proper density calculation for
        log-scale axes, and renders 6 contour levels at percentiles [10, 25, 50,
        75, 90, 95] to show population density gradients.

        Args:
            df (pd.DataFrame): DataFrame with gate events (raw data, source="raw")
            x_keyword (str): Keyword substring to identify X-axis channel column
                            (e.g., "FSC-A", "SSC-A", "B1-A")
            y_keyword (str): Keyword substring to identify Y-axis channel column
                            (e.g., "FSC-A", "SSC-A", "R2-A")
            sample_id (str): Sample ID to include in plot title
            gate_name (str): Gate name to include in plot title
            show_gates (bool, optional): Whether to display gate boundaries. Defaults to False.
            gates (list, optional): List of gate dictionaries to display. Defaults to None.
            selected_gate (str, optional): Name of the selected gate. Defaults to None.
            show_keywords (bool, optional): Whether to display keywords. Defaults to False.
            keywords (dict, optional): Dictionary of keywords to display. Defaults to None.
            gate_mode (str, optional): 'polygon' or 'quadrant'. Defaults to 'polygon'.
            quadrant_thresholds (dict, optional): Quadrant threshold values
                                                 {'x_threshold': float, 'y_threshold': float}.
                                                 Defaults to None.

        Returns:
            bokeh.plotting.figure: Bokeh figure object with contour plot rendered.
                                  Figure uses log-log axes (1 to 10⁵) matching FlowJo.

        Raises:
            ValueError: If no channel column matches either x_keyword or y_keyword.

        Technical Details:
            - Filters out values ≤ 0 for log scale compatibility
            - Downsamples to 10,000 points maximum for KDE performance
            - Creates 100×100 evaluation grid in log space
            - Uses scipy.stats.gaussian_kde for 2D density estimation
            - Uses matplotlib.pyplot.contour for contour line extraction
            - Renders contour paths in Bokeh (matplotlib only for calculation)
            - Supports gate overlays (quadrant dividers and polygon gates)
            - Annotations positioned in top-left corner for visibility

        Example:
            >>> # Create contour plot for FITC vs APC channels
            >>> df = workspace.get_gate_events("A1.fcs", "Singlets", source="raw")
            >>> p = analyzer._plot_contour(df, "B1-A", "R2-A", "A1.fcs", "Singlets")

            >>> # Create contour plot with quadrant dividers
            >>> thresholds = {'x_threshold': 3112.94, 'y_threshold': 444.11}
            >>> p = analyzer._plot_contour(df, "B1-A", "R2-A", "A1.fcs", "Singlets",
            ...                           gate_mode='quadrant',
            ...                           quadrant_thresholds=thresholds)

        See Also:
            _plot_scatter: For point-based visualization of 2D data
            _plot_histogram: For 1D density visualization
        """
        from bokeh.plotting import figure
        from scipy.stats import gaussian_kde
        import matplotlib.pyplot as plt
        import numpy as np

        # Find matching channels
        x_channels = [c for c in df.columns if x_keyword in c]
        y_channels = [c for c in df.columns if y_keyword in c]

        if not x_channels:
            raise ValueError(f"No channel found matching keyword '{x_keyword}'")
        if not y_channels:
            raise ValueError(f"No channel found matching keyword '{y_keyword}'")

        x_channel = x_channels[0]
        y_channel = y_channels[0]

        # Filter out values <= 0 for log scale
        df_filtered = df[(df[x_channel] > 0) & (df[y_channel] > 0)]

        if len(df_filtered) < 10:
            # Not enough data points for KDE
            print(f"    ⚠️  Warning: Only {len(df_filtered)} data points - cannot generate contours")
            # Return empty plot with message
            p = figure(title=f"{sample_id} - {gate_name} (Insufficient data)",
                      x_axis_label=x_channel, y_axis_label=y_channel,
                      x_axis_type="log", y_axis_type="log",
                      x_range=(1, 1e5), y_range=(1, 1e5),
                      width=400, height=300, tools="pan,wheel_zoom,box_zoom,reset,save",
                      sizing_mode="fixed")
            return p

        # Downsample if too many points (for KDE performance)
        max_points = 10000
        if len(df_filtered) > max_points:
            df_kde = df_filtered.sample(n=max_points, random_state=42)
        else:
            df_kde = df_filtered

        # Extract data
        x_data = df_kde[x_channel].values
        y_data = df_kde[y_channel].values

        # Work in log space (since axes are log scale)
        log_x = np.log10(x_data)
        log_y = np.log10(y_data)

        # Create 2D KDE
        xy = np.vstack([log_x, log_y])
        try:
            kde = gaussian_kde(xy)
        except Exception as e:
            print(f"    ⚠️  Warning: KDE calculation failed: {e}")
            # Return empty plot
            p = figure(title=f"{sample_id} - {gate_name} (KDE failed)",
                      x_axis_label=x_channel, y_axis_label=y_channel,
                      x_axis_type="log", y_axis_type="log",
                      x_range=(1, 1e5), y_range=(1, 1e5),
                      width=400, height=300, tools="pan,wheel_zoom,box_zoom,reset,save",
                      sizing_mode="fixed")
            return p

        # Create evaluation grid
        # Use 100x100 grid for smooth contours
        grid_points = 100
        x_grid_log = np.linspace(0, 5, grid_points)  # log10(1) to log10(100000)
        y_grid_log = np.linspace(0, 5, grid_points)
        X, Y = np.meshgrid(x_grid_log, y_grid_log)
        positions = np.vstack([X.ravel(), Y.ravel()])

        # Evaluate KDE on grid
        Z = kde(positions).reshape(X.shape)

        # Use matplotlib to calculate contours
        # We'll use percentile-based levels for consistent appearance
        sorted_density = np.sort(Z.ravel())
        sorted_density = sorted_density[sorted_density > 0]  # Remove zeros

        if len(sorted_density) == 0:
            print(f"    ⚠️  Warning: Zero density - cannot generate contours")
            p = figure(title=f"{sample_id} - {gate_name} (Zero density)",
                      x_axis_label=x_channel, y_axis_label=y_channel,
                      x_axis_type="log", y_axis_type="log",
                      x_range=(1, 1e5), y_range=(1, 1e5),
                      width=400, height=300, tools="pan,wheel_zoom,box_zoom,reset,save",
                      sizing_mode="fixed")
            return p

        # Define contour levels at specific percentiles
        # Higher percentiles = higher density = inner contours
        percentiles = [10, 25, 50, 75, 90, 95]
        levels = [np.percentile(sorted_density, p) for p in percentiles]

        # Create matplotlib figure for contour calculation (don't display)
        fig, ax = plt.subplots()
        cs = ax.contour(X, Y, Z, levels=levels)
        plt.close(fig)  # Close immediately - we just need the paths

        # Create Bokeh figure
        p = figure(title=f"{sample_id} - {gate_name}",
                  x_axis_label=x_channel, y_axis_label=y_channel,
                  x_axis_type="log", y_axis_type="log",
                  x_range=(1, 1e5), y_range=(1, 1e5),
                  width=400, height=300, tools="pan,wheel_zoom,box_zoom,reset,save",
                  sizing_mode="fixed")

        # Extract and render contour paths
        contour_count = 0
        for collection in cs.collections:
            for path in collection.get_paths():
                vertices = path.vertices
                if len(vertices) > 0:
                    # Convert from log space to linear space
                    contour_x = 10 ** vertices[:, 0]
                    contour_y = 10 ** vertices[:, 1]

                    # Render contour line in Bokeh
                    p.line(contour_x, contour_y, line_width=1.5,
                          color='navy', alpha=0.7)
                    contour_count += 1

        print(f"    ✓ Rendered {contour_count} contour lines")

        # Render quadrant dividers if in quadrant mode
        if gate_mode == 'quadrant' and quadrant_thresholds:
            print(f"    → Rendering quadrant dividers on log scale (1 to 10^5)")
            from bokeh.models import Span

            divider_color = "#DC143C"  # Crimson

            # Render vertical divider (X threshold)
            if 'x_threshold' in quadrant_thresholds:
                x_thresh = quadrant_thresholds['x_threshold']
                span = Span(location=x_thresh, dimension='height',
                           line_color=divider_color, line_width=2.5,
                           line_dash='dashed', line_alpha=0.9)
                p.add_layout(span)
                print(f"      ✓ Rendered vertical divider at x={x_thresh:.2f}")

            # Render horizontal divider (Y threshold)
            if 'y_threshold' in quadrant_thresholds:
                y_thresh = quadrant_thresholds['y_threshold']
                span = Span(location=y_thresh, dimension='width',
                           line_color=divider_color, line_width=2.5,
                           line_dash='dashed', line_alpha=0.9)
                p.add_layout(span)
                print(f"      ✓ Rendered horizontal divider at y={y_thresh:.2f}")

        # Render gate boundaries if requested (polygon mode)
        elif show_gates and gates:
            print(f"    → Rendering {len(gates)} gate(s) on log scale (1 to 10^5)")

            # Color palette for gates
            gate_colors = [
                ("#6B8E23", "#556B2F"),  # Olive green
                ("#4682B4", "#36648B"),  # Steel blue
                ("#CD853F", "#8B5A2B"),  # Peru/tan
                ("#9370DB", "#7B68EE"),  # Medium purple
                ("#DC143C", "#B22222"),  # Crimson
                ("#20B2AA", "#008B8B"),  # Light sea green
                ("#FF8C00", "#CD6600"),  # Dark orange
                ("#8B4789", "#68387E"),  # Dark orchid
            ]

            from bokeh.models import Span

            for gate_idx, gate in enumerate(gates):
                if (gate.get('x_dim') == x_channel and gate.get('y_dim') == y_channel) or \
                   (gate.get('x_dim') == y_channel and gate.get('y_dim') == x_channel):

                    gate_type = gate.get('type', 'polygon')

                    # Handle quadrant gates
                    if gate_type == 'quadrant':
                        dividers = gate.get('dividers', [])
                        if dividers:
                            fill_color, line_color = gate_colors[gate_idx % len(gate_colors)]

                            for divider in dividers:
                                orientation = divider.get('orientation')
                                value = divider.get('value')

                                if orientation == 'vertical' and value is not None:
                                    span = Span(location=value, dimension='height',
                                               line_color=line_color, line_width=2.5,
                                               line_dash='dashed', line_alpha=0.9)
                                    p.add_layout(span)
                                    print(f"      ✓ Rendered vertical divider at x={value:.2f} for '{gate.get('name')}'")
                                elif orientation == 'horizontal' and value is not None:
                                    span = Span(location=value, dimension='width',
                                               line_color=line_color, line_width=2.5,
                                               line_dash='dashed', line_alpha=0.9)
                                    p.add_layout(span)
                                    print(f"      ✓ Rendered horizontal divider at y={value:.2f} for '{gate.get('name')}'")

                            print(f"      ✓ Rendered quadrant gate '{gate.get('name')}' with {len(dividers)} divider(s), color: {line_color}")
                        else:
                            print(f"      ⚠️  Skipped '{gate.get('name')}': no dividers found")

                    # Handle polygon/rectangle gates
                    else:
                        xs = gate.get('xs', [])
                        ys = gate.get('ys', [])
                        if xs and ys and len(xs) == len(ys) and len(xs) >= 3:
                            # Swap if dimensions are reversed
                            if gate.get('x_dim') == y_channel:
                                xs, ys = ys, xs

                            print(f"      DEBUG: Gate '{gate.get('name')}' coordinates:")
                            print(f"             X: [{min(xs):.2f}, {max(xs):.2f}]")
                            print(f"             Y: [{min(ys):.2f}, {max(ys):.2f}]")

                            # Ensure polygon is closed
                            if xs[0] != xs[-1] or ys[0] != ys[-1]:
                                xs = list(xs) + [xs[0]]
                                ys = list(ys) + [ys[0]]

                            fill_color, line_color = gate_colors[gate_idx % len(gate_colors)]

                            p.patch(xs, ys, fill_color=fill_color, fill_alpha=0.3, line_color=line_color,
                                   line_width=2.5, line_alpha=0.9)
                            print(f"      ✓ Rendered '{gate.get('name')}' ({len(xs)} vertices, color: {fill_color})")
                        else:
                            print(f"      ⚠️  Skipped '{gate.get('name')}': invalid coordinates")
                else:
                    print(f"      ⚠️  Skipped '{gate.get('name')}': channel mismatch")

        # Add keyword display in top left if requested
        if show_keywords and keywords:
            from bokeh.models import Label
            keyword_lines = []
            for key, value in keywords.items():
                keyword_lines.append(f"{key}: {value}")
            keyword_text = "\n".join(keyword_lines)

            # Get plot dimensions for positioning
            x_min = x_data.min()
            x_max = x_data.max()
            y_max = y_data.max()

            # Position in top left with padding
            keyword_label = Label(x=x_min + (x_max - x_min) * 0.02, y=y_max * 0.98,
                                 text=keyword_text,
                                 text_font_size="8pt", text_color="black",
                                 text_align="left", text_baseline="top",
                                 background_fill_color="white",
                                 background_fill_alpha=0.7, border_line_color="gray",
                                 border_line_width=1, x_units="data", y_units="data")
            p.add_layout(keyword_label)

        return p

    def generate_interactive_plots(self, selections, output_path):
        """
        Generate plots for all samples based on user selections and save to HTML.
        
        Iterates through all samples in the workspace, extracts gate events using raw data
        (source="raw"), and generates plots according to user selections. All plots are
        combined into a single HTML file with a grid layout (2 columns).
        
        For "Ungated" gates, uses sample.get_events(source='raw') instead of get_gate_events().
        For other gates, uses workspace.get_gate_events() with source="raw".
        
        Args:
            selections (dict): Dictionary from interactive_plot_prompt() containing:
                - 'gate_path': tuple of gate path
                - 'gate_name': str gate name
                - 'plot_type': "histogram" or "scatter"
                - 'parameters': list of channel keywords
            output_path (str): Path where HTML file will be saved (e.g., "plots.html")
            
        Output:
            Creates an HTML file with:
            - Summary section showing plot configuration
            - Grid of plots (2 columns) for all samples
            - Each plot titled with sample ID, gate name, and parameters
            
        Example:
            >>> selections = {
            ...     'gate_path': ('root', 'Cells', 'Singlets'),
            ...     'gate_name': 'Singlets',
            ...     'plot_type': 'histogram',
            ...     'parameters': ['RL1-A']
            ... }
            >>> analyzer.generate_interactive_plots(selections, "output.html")
            >>> # Generates histograms of RL1-A for Singlets gate across all samples
        """
        from bokeh.layouts import gridplot, column
        from bokeh.plotting import output_file, save
        from bokeh.models import Div, TabPanel, Tabs
        
        # Get selected groups and sample IDs
        selected_groups = selections.get('selected_groups', self.sample_groups[:1])
        show_keywords = selections.get('show_keywords', False)
        selected_keywords_to_show = selections.get('selected_keywords_to_show', None)
        show_statistics = selections.get('show_statistics', False)
        plot_configs = selections.get('plot_configs', [])
        
        # If plot_configs is not present, use old format (backward compatibility)
        if not plot_configs:
            # Convert old format to new format
            plot_configs = [{
                'gate_path': selections.get('gate_path'),
                'gate_name': selections.get('gate_name'),
                'plot_type': selections.get('plot_type'),
                'parameters': selections.get('parameters'),
                'statistic': selections.get('statistic', 'median'),
                'scale': selections.get('scale', 'linear')
            }]
        
        # Collect all sample IDs from selected groups
        all_sample_ids = []
        for group in selected_groups:
            try:
                group_samples = self.workspace.get_sample_ids(group, loaded_only=True)
                all_sample_ids.extend(group_samples)
            except Exception as e:
                print(f"Warning: Could not get samples from group '{group}': {e}")
        
        # Remove duplicates while preserving order
        seen = set()
        sample_ids = [sid for sid in all_sample_ids if not (sid in seen or seen.add(sid))]
        
        if not sample_ids:
            print("No samples found in selected groups.")
            return
        
        # Analyze samples if not already done
        print("Analyzing samples...")
        for group in selected_groups:
            try:
                self.workspace.analyze_samples(group_name=group)
            except Exception as e:
                print(f"Warning: Could not analyze group '{group}': {e}")
        
        # Generate plots for each configuration
        tabs = []
        tab_titles = []
        
        for plot_idx, plot_config in enumerate(plot_configs, 1):
            gate_path = plot_config['gate_path']
            gate_name = plot_config['gate_name']
            plot_type = plot_config['plot_type']
            parameters = plot_config['parameters']
            
            print(f"\n{'='*60}")
            print(f"Generating Plot {plot_idx} of {len(plot_configs)}")
            print(f"{'='*60}")
            print(f"Gate: {gate_name}")
            print(f"Type: {plot_type}")
            print(f"Parameters: {', '.join(parameters)}")
            print(f"\nGenerating {plot_type} plots for {len(sample_ids)} samples...")
            print(f"Using well ID source: {selections.get('well_id_source', 'auto')}")
            if selections.get('well_id_keyword'):
                print(f"Using keyword: '{selections.get('well_id_keyword')}'")
            print()
            
            # Organize plots by well position (Row, Col)
            plots_by_well = {}  # (row, col) -> plot
            successful_samples = []
            failed_samples = []
            sample_info = {}  # (row, col) -> (sample_id, filename)
            
            for i, sample_id in enumerate(sample_ids, 1):
                try:
                    sample = self.workspace.get_sample(sample_id)
                    print(f"  [{i}/{len(sample_ids)}] Processing {sample_id}...")
                    
                    # Parse well ID using user-selected source
                    well_id_source = selections.get('well_id_source', 'auto')
                    well_id_keyword = selections.get('well_id_keyword')
                    r, c, method_used = self.parse_well_id(sample, source=well_id_source, keyword_name=well_id_keyword, return_method=True, sample_id=sample_id)
                    
                    if not r or not c:
                        # If no well ID, use sample index as fallback
                        print(f"    ⚠️  Well ID not found, using fallback position")
                        r = f"Row{i//12 + 1}"
                        c = (i % 12) + 1
                        method_used = "fallback (sample index)"
                    else:
                        print(f"    ✓ Found well ID: {r}{c:02d} (extracted from {method_used})")
                    
                    # Get gate events with raw data
                    # Check if we should use the selected gate's data directly (when "use parent gate" option was selected)
                    use_gate_directly = plot_config.get('use_gate_data_directly', False)
                    
                    if use_gate_directly:
                        # Use the selected gate's data directly (for visualizing child gates on top)
                        print(f"    → Loading data from gate '{gate_name}'...", end=" ")
                        if gate_name == "Ungated" or (not gate_path or gate_path == ()):
                            # Ungated case - use raw events
                            events = sample.get_events(source='raw')
                            if events is not None:
                                if not isinstance(events, pd.DataFrame):
                                    df = pd.DataFrame(events, columns=sample.pnn_labels)
                                else:
                                    df = events
                            else:
                                df = None
                        else:
                            # Get the selected gate's events directly
                            try:
                                df = self.workspace.get_gate_events(sample_id, gate_name, gate_path=gate_path, source="raw")
                            except Exception as e:
                                print(f"\n    ⚠️  Could not load gate data: {e}")
                                df = None
                    else:
                        # IMPORTANT: Load PRE-FILTERED data (parent gate), not the selected gate's data
                        # This allows us to visualize the gate boundary overlaid on ungated/parent populations
                        print(f"    → Loading PRE-FILTERED data (before '{gate_name}' gate)...", end=" ")
                        if gate_name == "Ungated" or (not gate_path or gate_path == ()):
                            # Ungated case - use raw events
                            events = sample.get_events(source='raw')
                            if events is not None:
                                # Convert to DataFrame if needed
                                if not isinstance(events, pd.DataFrame):
                                    df = pd.DataFrame(events, columns=sample.pnn_labels)
                                else:
                                    df = events
                            else:
                                df = None
                        elif gate_path == ('root',):
                            # First-level gate (like Cells_withDebris) - parent is raw ungated data
                            events = sample.get_events(source='raw')
                            if events is not None:
                                if not isinstance(events, pd.DataFrame):
                                    df = pd.DataFrame(events, columns=sample.pnn_labels)
                                else:
                                    df = events
                            else:
                                df = None
                        else:
                            # Multi-level gate - get parent gate's data
                            parent_gate_path = gate_path[:-1]
                            # Get the parent gate name from the workspace
                            try:
                                # Navigate to parent gate
                                if parent_gate_path == ('root',):
                                    # Parent is a root-level gate, need to find it by name
                                    # For now, use raw events as fallback
                                    events = sample.get_events(source='raw')
                                    if events is not None:
                                        if not isinstance(events, pd.DataFrame):
                                            df = pd.DataFrame(events, columns=sample.pnn_labels)
                                        else:
                                            df = events
                                    else:
                                        df = None
                                else:
                                    # Use the parent path to get parent gate events
                                    parent_gate_name = parent_gate_path[-1]
                                    df = self.workspace.get_gate_events(sample_id, parent_gate_name, gate_path=parent_gate_path, source="raw")
                            except Exception as e:
                                print(f"\n    ⚠️  Could not load parent gate data: {e}")
                                df = None
                    
                    if df is None or df.empty:
                        print("SKIPPED")
                        print(f"    ✗ DataFrame is empty - no events in gate '{gate_name}'")
                        failed_samples.append((sample_id, f"No events in gate '{gate_name}'"))
                        continue
                    
                    print(f"OK ({len(df)} events)")
                    
                    # Parse well ID for labeling
                    well_id = f"{r}{c:02d}"
                    
                    # Get keywords/statistics if requested (filtered by selection)
                    keywords = None
                    if show_keywords:
                        if show_statistics:
                            # Display statistics instead of keywords
                            # Calculate statistics from the data
                            event_count = len(df)  # Total number of events/cells
                            
                            if plot_type == "histogram":
                                channel = [c for c in df.columns if parameters[0] in c][0]
                                valid_data = df[channel].dropna()
                                valid_data = valid_data[valid_data >= 10]  # Apply same filter as histogram
                                if len(valid_data) > 0:
                                    stats_dict = {}
                                    if 'count' in selected_keywords_to_show:
                                        stats_dict['Events'] = f"{event_count:,}"
                                    if 'median' in selected_keywords_to_show:
                                        stats_dict['Median'] = f"{valid_data.median():.1f}"
                                    if 'mean' in selected_keywords_to_show:
                                        stats_dict['Mean'] = f"{valid_data.mean():.1f}"
                                    keywords = stats_dict
                            else:  # scatter
                                # For scatter, calculate statistics for both channels
                                x_channel = [c for c in df.columns if parameters[0] in c][0]
                                y_channel = [c for c in df.columns if parameters[1] in c][0]
                                x_data = df[x_channel].dropna()
                                y_data = df[y_channel].dropna()
                                stats_dict = {}
                                if 'count' in selected_keywords_to_show:
                                    stats_dict['Events'] = f"{event_count:,}"
                                if 'median' in selected_keywords_to_show:
                                    if len(x_data) > 0:
                                        stats_dict[f'{x_channel} Median'] = f"{x_data.median():.1f}"
                                    if len(y_data) > 0:
                                        stats_dict[f'{y_channel} Median'] = f"{y_data.median():.1f}"
                                if 'mean' in selected_keywords_to_show:
                                    if len(x_data) > 0:
                                        stats_dict[f'{x_channel} Mean'] = f"{x_data.mean():.1f}"
                                    if len(y_data) > 0:
                                        stats_dict[f'{y_channel} Mean'] = f"{y_data.mean():.1f}"
                                keywords = stats_dict
                        else:
                            # Display keywords
                            try:
                                all_keywords = self.workspace.get_keywords(sample_id)
                                # Filter to only selected keywords if specified
                                if selected_keywords_to_show:
                                    keywords = {k: v for k, v in all_keywords.items() if k in selected_keywords_to_show}
                                else:
                                    keywords = all_keywords
                            except:
                                keywords = {}
                    
                    # Extract gates for visualization (for scatter plots)
                    # This can be different from the gate used for filtering
                    gates_for_viz = []
                    if plot_type == "scatter":
                        try:
                            # Find matching channels for gate extraction
                            x_channels = [c for c in df.columns if parameters[0] in c]
                            y_channels = [c for c in df.columns if parameters[1] in c]
                            if x_channels and y_channels:
                                x_channel = x_channels[0]
                                y_channel = y_channels[0]

                                # Get list of gates to visualize from plot config
                                gates_to_viz_list = plot_config.get('gates_to_visualize', [])
                                show_gates = plot_config.get('show_gates', False)

                                print(f"    DEBUG: show_gates={show_gates}, gates_to_visualize={gates_to_viz_list}")

                                if gates_to_viz_list:
                                    # Extract specific gates by name
                                    print(f"    → Extracting gates for visualization: {', '.join(gates_to_viz_list)}")

                                    # Check if the currently selected gate is in the visualization list
                                    if gate_name in gates_to_viz_list and gate_name != "Ungated":
                                        # Extract the selected gate itself
                                        print(f"      → Attempting to extract selected gate '{gate_name}'...")
                                        print(f"         Gate path: {gate_path}")
                                        print(f"         X channel: '{x_channel}' (type: {type(x_channel).__name__})")
                                        print(f"         Y channel: '{y_channel}' (type: {type(y_channel).__name__})")
                                        selected_gate_data = self._extract_selected_gate(sample_id, gate_name, gate_path, x_channel, y_channel)
                                        if selected_gate_data:
                                            gates_for_viz.append(selected_gate_data)
                                            # Handle different gate types in logging
                                            if selected_gate_data.get('type') == 'quadrant':
                                                num_dividers = len(selected_gate_data.get('dividers', []))
                                                print(f"      ✓ Extracted '{gate_name}' (selected gate, quadrant with {num_dividers} divider(s))")
                                            else:
                                                num_vertices = len(selected_gate_data.get('xs', []))
                                                print(f"      ✓ Extracted '{gate_name}' (selected gate, {num_vertices} vertices)")
                                        else:
                                            print(f"      ⚠️  Failed to extract selected gate '{gate_name}' (returned None)")
                                            print(f"         Check warnings above for details")
                                    
                                    # Extract child gates of the selected gate (e.g., Q1-Q4 quadrants under Singlets)
                                    # Child gates have a path that starts with the selected gate's path + gate name
                                    if gate_path:
                                        child_gate_path_prefix = gate_path + (gate_name,)
                                        print(f"      → Looking for child gates under '{gate_name}' (path prefix: {' → '.join(child_gate_path_prefix)})")
                                        all_gate_ids = self.workspace.get_gate_ids(sample_id)
                                        for child_gate_name, child_gate_path in all_gate_ids:
                                            # Check if this is a child of the selected gate
                                            # Child path should be longer and start with parent path + parent name
                                            if (len(child_gate_path) > len(child_gate_path_prefix) and 
                                                child_gate_path[:len(child_gate_path_prefix)] == child_gate_path_prefix):
                                                # This is a child gate - try to extract it if it matches channels
                                                if child_gate_name in gates_to_viz_list or 'all' in gates_to_viz_list:
                                                    child_gate_data = self._extract_selected_gate(sample_id, child_gate_name, child_gate_path, x_channel, y_channel)
                                                    if child_gate_data:
                                                        # Avoid duplicates
                                                        if not any(g['name'] == child_gate_data['name'] for g in gates_for_viz):
                                                            gates_for_viz.append(child_gate_data)
                                                            if child_gate_data.get('type') == 'quadrant':
                                                                num_dividers = len(child_gate_data.get('dividers', []))
                                                                print(f"      ✓ Extracted child gate '{child_gate_name}' (quadrant with {num_dividers} divider(s))")
                                                            else:
                                                                num_vertices = len(child_gate_data.get('xs', []))
                                                                print(f"      ✓ Extracted child gate '{child_gate_name}' ({num_vertices} vertices)")

                                    # Also extract any other gates from the workspace
                                    all_gates = self._extract_gate_polygons(sample_id, x_channel, y_channel)
                                    for gate_data in all_gates:
                                        if gate_data['name'] in gates_to_viz_list:
                                            # Avoid duplicates
                                            if not any(g['name'] == gate_data['name'] for g in gates_for_viz):
                                                gates_for_viz.append(gate_data)
                                                # Handle different gate types in logging
                                                if gate_data.get('type') == 'quadrant':
                                                    num_dividers = len(gate_data.get('dividers', []))
                                                    print(f"      ✓ Extracted '{gate_data['name']}' (quadrant with {num_dividers} divider(s))")
                                                else:
                                                    num_vertices = len(gate_data.get('xs', []))
                                                    print(f"      ✓ Extracted '{gate_data['name']}' ({num_vertices} vertices)")

                                    if not gates_for_viz:
                                        print(f"      ⚠️  Could not extract any of the requested gates")
                        except Exception as e:
                            print(f"    ⚠️  Error extracting gates: {e}")
                            gates_for_viz = []

                    # Generate plot based on type (pass well_id for title)
                    print(f"    → Generating {plot_type} plot...", end=" ")
                    if plot_type == "histogram":
                        # Get statistic and scale choices
                        statistic = plot_config.get('statistic', 'median')
                        scale = plot_config.get('scale', 'linear')
                        p = self._plot_histogram(df, parameters[0], well_id, gate_name,
                                               statistic=statistic, scale=scale,
                                               show_keywords=show_keywords, keywords=keywords)
                        # Shortened title - just well ID and gate name
                        p.title.text = f"{well_id} - {gate_name}"
                    elif plot_type == "scatter":
                        # Get gate_mode and quadrant_thresholds from plot_config
                        plot_gate_mode = plot_config.get('gate_mode', 'polygon')
                        plot_quadrant_thresholds = plot_config.get('quadrant_thresholds', None)

                        # Pass the extracted gates for visualization
                        plot_show_gates = len(gates_for_viz) > 0
                        p = self._plot_scatter(df, parameters[0], parameters[1], well_id, gate_name,
                                             gates=gates_for_viz, show_gates=plot_show_gates,
                                             selected_gate=None,  # Not using selected_gate anymore
                                             keywords=keywords, show_keywords=show_keywords,
                                             gate_mode=plot_gate_mode,
                                             quadrant_thresholds=plot_quadrant_thresholds)
                        # Shortened title - just well ID and gate name
                        p.title.text = f"{well_id} - {gate_name}"
                    else:  # contour
                        # Get gate_mode and quadrant_thresholds from plot_config
                        plot_gate_mode = plot_config.get('gate_mode', 'polygon')
                        plot_quadrant_thresholds = plot_config.get('quadrant_thresholds', None)

                        # Pass the extracted gates for visualization
                        plot_show_gates = len(gates_for_viz) > 0
                        p = self._plot_contour(df, parameters[0], parameters[1], well_id, gate_name,
                                             gates=gates_for_viz, show_gates=plot_show_gates,
                                             selected_gate=None,
                                             keywords=keywords, show_keywords=show_keywords,
                                             gate_mode=plot_gate_mode,
                                             quadrant_thresholds=plot_quadrant_thresholds)
                        # Shortened title - just well ID and gate name
                        p.title.text = f"{well_id} - {gate_name}"
                    
                    # Store plot by well position
                    plots_by_well[(r, c)] = p
                    sample_info[(r, c)] = (sample_id, sample.original_filename, well_id)
                    successful_samples.append(sample_id)
                    print("OK")
                    print(f"    ✓ Plot saved for well {well_id}")
                    print()
                    
                except Exception as e:
                    print(f"FAILED")
                    print(f"    ✗ Error: {e}")
                    failed_samples.append((sample_id, str(e)))
                    print()
                    continue
        
            # Print summary for this plot
            print("="*60)
            print(f"Processing Summary - Plot {plot_idx}")
            print("="*60)
            print(f"Total samples processed: {len(sample_ids)}")
            print(f"Successful plots: {len(successful_samples)}")
            print(f"Failed/Skipped: {len(failed_samples)}")
            
            if not plots_by_well:
                print("\nERROR: No plots were generated successfully for this configuration.")
                if failed_samples:
                    print("\nFailed samples:")
                    for sid, reason in failed_samples:
                        print(f"  - {sid}: {reason}")
                # Create empty tab
                empty_div = Div(text=f"<h2>No plots generated</h2><p>All samples failed for this configuration.</p>")
                tabs.append(TabPanel(child=empty_div, title=f"Plot {plot_idx}"))
                tab_titles.append(f"Plot {plot_idx}")
                continue
            
            if failed_samples:
                print("\nFailed/Skipped samples:")
                for sid, reason in failed_samples:
                    print(f"  - {sid}: {reason}")
            print()
            
            # Create 96-well plate layout (8 rows x 12 columns)
            # Rows: A-H, Columns: 1-12
            rows = list("ABCDEFGH")
            cols = list(range(1, 13))
            
            # Create grid of plots organized by well
            plot_grid = []
            for row in rows:
                row_plots = []
                for col in cols:
                    well_key = (row, col)
                    well_id = f"{row}{col:02d}"
                    if well_key in plots_by_well:
                        # Get plot (title already includes well ID and filename from generation step)
                        p = plots_by_well[well_key]
                        
                        # Ensure title is clear and visible
                        p.title.text_font_size = "11pt"
                        p.title.text_font_style = "bold"
                        p.title.align = "center"
                        
                        row_plots.append(p)
                    else:
                        # Empty well - create placeholder with clear label and a renderer
                        from bokeh.plotting import figure
                        empty_p = figure(width=400, height=300, 
                                       title=f"Well {well_id}: No Data",
                                       tools="", toolbar_location=None,
                                       x_range=(0, 1), y_range=(0, 1))
                        empty_p.axis.visible = False
                        empty_p.grid.visible = False
                        empty_p.title.text_font_size = "11pt"
                        empty_p.title.text_font_style = "bold"
                        empty_p.title.text_color = "gray"
                        empty_p.title.align = "center"
                        # Add a text renderer to avoid missing renderers warning
                        empty_p.text(x=[0.5], y=[0.5], text=["No Data"], 
                                   text_font_size="14pt", text_color="gray",
                                   text_align="center", text_baseline="middle")
                        row_plots.append(empty_p)
                plot_grid.append(row_plots)
            
            # Use fixed sizing for gridplot to prevent overlapping
            # Note: Individual plots already have width=400, height=300 set when created
            plot_layout = gridplot(plot_grid, toolbar_location="right", sizing_mode="fixed")
            
            # Create tab title
            plot_type_display = plot_config['plot_type'].capitalize()
            gate_display = plot_config['gate_name']
            params_display = '_'.join(plot_config['parameters'])
            tab_title = f"{plot_type_display} - {gate_display}"
            if len(tab_title) > 30:
                tab_title = f"Plot {plot_idx}"
            
            tabs.append(TabPanel(child=plot_layout, title=tab_title))
            tab_titles.append(tab_title)
        
        # Create tabs layout if multiple plots, otherwise just use single layout
        if len(tabs) > 1:
            layout = Tabs(tabs=tabs)
        else:
            layout = tabs[0].child if tabs else column()
        
        # Add summary with Tailwind CSS styling
        from bokeh.models import Div
        # Minimalistic Tailwind CSS styling (greptile-inspired)
        tailwind_cdn = """
        <script src="https://cdn.tailwindcss.com"></script>
        <style>
            body { 
                font-family: 'Inter', -apple-system, BlinkMacSystemFont, sans-serif; 
                background-color: #ffffff; 
                margin: 0;
                padding: 0;
            }
            .bk-root {
                margin: 0;
                padding: 0;
            }
            .bk-plot {
                margin: 1px;
                border: 1px solid #e5e7eb;
                border-radius: 2px;
            }
            /* Modal styles */
            .modal {
                display: none;
                position: fixed;
                z-index: 1000;
                left: 0;
                top: 0;
                width: 100%;
                height: 100%;
                background-color: rgba(0,0,0,0.4);
            }
            .modal-content {
                background-color: #ffffff;
                margin: 10% auto;
                padding: 2rem;
                border: 1px solid #e5e7eb;
                border-radius: 0.5rem;
                width: 90%;
                max-width: 600px;
                box-shadow: 0 10px 25px rgba(0,0,0,0.1);
            }
            .close {
                color: #6b7280;
                float: right;
                font-size: 28px;
                font-weight: bold;
                cursor: pointer;
                line-height: 20px;
            }
            .close:hover {
                color: #000;
            }
            .close:focus {
                outline: none;
            }
        </style>
        """
        
        # Create modal HTML for configuration/summary
        keyword_info = ""
        if 'keyword_filter' in selections:
            kf = selections['keyword_filter']
            keyword_info = f'<p class="text-sm"><span class="font-semibold text-gray-700">Keyword Filter:</span> <span class="text-gray-900">{kf["key"]} = {kf["value"]}</span></p>'
        
        modal_html = f"""
        {tailwind_cdn}
        <div id="infoModal" class="modal" style="display: none;">
            <div class="modal-content">
                <span class="close" onclick="closeModal()">&times;</span>
                <h2 class="text-2xl font-bold text-gray-900 mb-4">Plot Configuration</h2>
                <div class="space-y-3 mb-6">
                    {keyword_info}
                    <p class="text-sm"><span class="font-semibold text-gray-700">Gate Path:</span> <span class="text-gray-900">{' → '.join(gate_path) if gate_path else 'root'}</span></p>
                    <p class="text-sm"><span class="font-semibold text-gray-700">Gate Name:</span> <span class="text-gray-900">{gate_name}</span></p>
                    <p class="text-sm"><span class="font-semibold text-gray-700">Plot Type:</span> <span class="capitalize text-gray-900">{plot_type}</span></p>
                    <p class="text-sm"><span class="font-semibold text-gray-700">Parameters:</span> <span class="text-gray-900">{', '.join(parameters)}</span></p>
                    {f'<p class="text-sm"><span class="font-semibold text-gray-700">Statistic:</span> <span class="text-gray-900 capitalize">{selections.get("statistic", "median")}</span></p>' if plot_type == "histogram" else ''}
                </div>
                <h3 class="text-xl font-semibold text-gray-900 mb-3">Processing Summary</h3>
                <div class="space-y-2">
                    <p class="text-sm"><span class="font-semibold text-gray-700">Total Samples:</span> <span class="text-gray-900">{len(sample_ids)}</span></p>
                    <p class="text-sm"><span class="font-semibold text-gray-700">Successful:</span> <span class="text-green-600 font-semibold">{len(successful_samples)}</span> <span class="text-gray-500">samples</span></p>
                    <p class="text-sm"><span class="font-semibold text-gray-700">Failed:</span> <span class="text-red-600 font-semibold">{len(failed_samples)}</span> <span class="text-gray-500">samples</span></p>
                    <p class="text-sm"><span class="font-semibold text-gray-700">Wells with Data:</span> <span class="text-blue-600 font-semibold">{len(plots_by_well)}</span> <span class="text-gray-500">wells</span></p>
                </div>
            </div>
        </div>
        """
        
        # Get base filename and plot type for title
        base_filename = os.path.splitext(os.path.basename(self.wsp_path))[0]
        plot_type_display = plot_type.title() if plot_type else "Analysis"
        
        # Add keyword filter to title if present
        title_suffix = ""
        if 'keyword_filter' in selections:
            kf = selections['keyword_filter']
            title_suffix = f" [{kf['key']}={kf['value']}]"
        
        # Minimalistic header with info button
        header_html = f"""
        {tailwind_cdn}
        <div class="w-full border-b border-gray-200 bg-white sticky top-0 z-40">
            <div class="max-w-7xl mx-auto flex items-center justify-between py-3 px-4">
                <h1 class="text-xl font-bold text-gray-900">{base_filename} Analysis - {plot_type_display}{title_suffix}</h1>
                <button id="infoButton" onclick="openInfoModal()" 
                        class="px-4 py-2 text-sm font-medium text-gray-700 hover:text-gray-900 border border-gray-300 rounded-lg hover:bg-gray-50 transition">
                    ℹ️ Info
                </button>
            </div>
        </div>
        {modal_html}
        <script>
            function openInfoModal() {{
                var modal = document.getElementById('infoModal');
                if (modal) {{
                    modal.style.display = 'block';
                }}
            }}
            function closeModal() {{
                document.getElementById('infoModal').style.display = 'none';
            }}
            // Close modal when clicking outside of it
            window.onclick = function(event) {{
                var modal = document.getElementById('infoModal');
                if (event.target == modal) {{
                    closeModal();
                }}
            }}
            // Close modal with Escape key
            document.addEventListener('keydown', function(event) {{
                if (event.key === 'Escape') {{
                    closeModal();
                }}
            }});
        </script>
        """
        
        header_div = Div(text=header_html, width=1400, sizing_mode="fixed")

        # Create global gate legend (if any gates are visualized across any plot config)
        # Collect all unique gates from all plot configs
        all_gates_to_viz = []
        for pc in plot_configs:
            if pc.get('plot_type') == 'scatter' and pc.get('show_gates', False):
                gates_list = pc.get('gates_to_visualize', [])
                for gate in gates_list:
                    if gate not in all_gates_to_viz:
                        all_gates_to_viz.append(gate)

        # Define the same color palette used in scatter plot rendering
        gate_colors = [
            ("#6B8E23", "#556B2F"),  # Olive green
            ("#4682B4", "#36648B"),  # Steel blue
            ("#CD853F", "#8B5A2B"),  # Peru/tan
            ("#9370DB", "#7B68EE"),  # Medium purple
            ("#DC143C", "#B22222"),  # Crimson
            ("#20B2AA", "#008B8B"),  # Light sea green
            ("#FF8C00", "#CD6600"),  # Dark orange
            ("#8B4789", "#68387E"),  # Dark orchid
        ]

        # Build legend HTML
        legend_items_html = ""
        if all_gates_to_viz:
            for idx, gate_name in enumerate(all_gates_to_viz):
                fill_color, line_color = gate_colors[idx % len(gate_colors)]
                legend_items_html += f"""
                <div class="flex items-center mr-6">
                    <div style="width: 20px; height: 20px; background-color: {fill_color}; border: 2px solid {line_color}; border-radius: 3px; margin-right: 8px;"></div>
                    <span class="text-sm text-gray-700">{gate_name}</span>
                </div>
                """

            legend_html = f"""
            {tailwind_cdn}
            <div class="w-full bg-gray-50 border-b border-gray-200">
                <div class="max-w-7xl mx-auto py-3 px-4">
                    <div class="flex items-center">
                        <span class="text-sm font-semibold text-gray-700 mr-4">Gates:</span>
                        <div class="flex flex-wrap items-center">
                            {legend_items_html}
                        </div>
                    </div>
                </div>
            </div>
            """
            legend_div = Div(text=legend_html, width=1400, sizing_mode="fixed")
        else:
            legend_div = None

        # Use fixed sizing to prevent layout issues and overlapping
        # Header with modal, optional legend, then plots (no summary div - it's in the modal)
        if legend_div:
            final_layout = column(
                header_div,
                legend_div,
                layout,
                sizing_mode="fixed"
            )
        else:
            final_layout = column(
                header_div,
                layout,
                sizing_mode="fixed"
            )
        
        # Save to HTML
        output_file(output_path)
        save(final_layout)
        print(f"\n✓ Plots saved to: {output_path}")
        print(f"  Generated {len(successful_samples)} plots successfully")
        if failed_samples:
            print(f"  {len(failed_samples)} samples failed (see output above)")

    # Standard report function removed: generate_interactive_plot_from_table
def main():
    """
    Main entry point for the FlowJo Analysis Tool.
    
    Supports two modes:
    - Interactive plotting mode: prompts user for gate selection, plot type, and parameters
    - Inspect mode: displays gate hierarchy without generating plots
    
    Examples:
        # Interactive plotting
        python analyze_flow.py --wsp workspace.wsp --fcs_dir ./fcs --interactive
        
        # Inspect gates
        python analyze_flow.py --wsp workspace.wsp --fcs_dir ./fcs --inspect
        python analyze_flow.py --wsp workspace.wsp --fcs_dir ./fcs --inspect --sample sample1.fcs
    """
    args = parse_args()
    
    if not os.path.exists(args.wsp):
        print(f"Error: Workspace file not found at {args.wsp}")
        sys.exit(1)
    
    # Check that exactly one mode is selected
    if args.interactive and args.inspect:
        print("Error: Cannot use both --interactive and --inspect. Choose one mode.")
        sys.exit(1)
    
    if not args.interactive and not args.inspect:
        print("Error: Must specify either --interactive or --inspect mode.")
        print("Usage:")
        print("  Interactive mode: python analyze_flow.py --wsp <workspace.wsp> --fcs_dir <fcs_dir> --interactive")
        print("  Inspect mode:     python analyze_flow.py --wsp <workspace.wsp> --fcs_dir <fcs_dir> --inspect")
        sys.exit(1)
    
    if args.inspect:
        # Inspect mode - display gate hierarchy
        # Import here to avoid circular dependencies
        from inspect_gates import inspect_gates
        inspect_gates(args.wsp, args.fcs_dir, args.sample)
    else:
        # Interactive plotting mode
        analyzer = FlowAnalyzer(args.wsp, args.fcs_dir)
        
        # Prompts user for gate selection, plot type, and parameters
        # Then generates plots for all samples and saves to HTML
        selections = analyzer.interactive_plot_prompt()
        if selections:
            # Generate output filename based on input name + type + parameters + keyword filter
            base_name = os.path.splitext(os.path.basename(args.wsp))[0]
            
            # Get plot configs (new format) or use old format for backward compatibility
            plot_configs = selections.get('plot_configs', [])
            if not plot_configs:
                # Old format - create a single plot config
                plot_configs = [{
                    'plot_type': selections.get('plot_type', 'histogram'),
                    'gate_name': selections.get('gate_name', 'Ungated'),
                    'parameters': selections.get('parameters', [])
                }]
            
            # Use first plot config for filename (or combine all if multiple)
            if len(plot_configs) == 1:
                plot_type = plot_configs[0]['plot_type']
                gate_name = plot_configs[0]['gate_name'].replace(' ', '_').replace('/', '_')
                params_str = '_'.join(plot_configs[0]['parameters']).replace(' ', '_').replace('/', '_')
            else:
                # Multiple plots - use generic name
                plot_type = 'multi'
                gate_name = 'multi'
                params_str = f"{len(plot_configs)}plots"
            
            # Add keyword filter to filename if present
            filter_str = ""
            if 'keyword_filter' in selections:
                kf = selections['keyword_filter']
                key_clean = kf['key'].replace(' ', '_').replace('/', '_').replace('(', '').replace(')', '')
                value_clean = str(kf['value']).replace(' ', '_').replace('/', '_')
                filter_str = f"_{key_clean}_{value_clean}"
            
            if args.output != "report.html":
                output_path = args.output
            else:
                output_path = f"{base_name}_{plot_type}_{gate_name}_{params_str}{filter_str}.html"
            
            analyzer.generate_interactive_plots(selections, output_path)
        else:
            print("No selections made. Exiting.")

if __name__ == "__main__":
    main()
