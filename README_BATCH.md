# Batch Flow Cytometry Analysis Scripts

Automated batch processing scripts for flow cytometry data analysis using FlowJo workspaces and the FlowKit library.

## Overview

These scripts enable automated generation of interactive HTML plots from FlowJo workspace (.wsp) files without manual intervention. Designed for experiments with repeated measurements (e.g., time series, dose-response curves, replicate plates) where you need to generate the same plot configuration across multiple conditions.

## Features

- **Automated batch processing**: Process multiple keyword filter values without manual input
- **Interactive HTML output**: Generate Bokeh-based interactive plots with pan/zoom/save tools
- **96-well plate layout**: Automatic organization by well position
- **Gate visualization**: Overlay gate boundaries on scatter plots
- **Flexible configuration**: Customize via simple configuration variables
- **Global legend**: Single legend showing all visualized gates with distinct colors
- **Error handling**: Continues processing on failure, reports summary at end

## Requirements

```bash
pip install flowkit bokeh pandas numpy scipy matplotlib
```

### Python Packages
- `flowkit` - FlowJo workspace parsing and FCS file handling
- `bokeh` - Interactive HTML visualizations
- `pandas` - Data manipulation
- `numpy` - Numerical operations
- `scipy` - Statistical functions (KDE for histograms)
- `matplotlib` - Gate polygon operations (optional)

## Files

### `analyze_flow.py`
Main analysis script with interactive CLI mode. Handles:
- FlowJo workspace loading
- Gate extraction and visualization
- Plot generation (histograms and scatter plots)
- Well ID parsing from keywords or filenames
- Log-scale transformations
- 96-well plate layout organization

**Usage:**
```bash
python analyze_flow.py --wsp workspace.wsp --fcs_dir ./fcs_files --interactive
```

### `batch_analyze_flow.py`
Batch processing wrapper that automates `analyze_flow.py` for multiple keyword filter values.

**Key Features:**
- Iterates over specified filter values (e.g., time points)
- Generates one HTML file per filter value
- No manual intervention required
- Configurable via simple variables at top of script

## Quick Start

### 1. Configure the Script

Edit the **CONFIGURATION** section in `batch_analyze_flow.py`:

```python
# File paths
wsp_path = "/path/to/your/workspace.wsp"
fcs_dir = "/path/to/your/fcs_files"
output_dir = "/path/to/output"

# Analysis parameters
filter_values = ['48', '72', '96', '144']  # e.g., time points in hours
gate_path = ('root',)                      # Root-level gates
gate_name = 'Cells'                        # Gate to visualize
x_channel = 'FSC-A'                        # X-axis channel
y_channel = 'SSC-A'                        # Y-axis channel
keyword_filter_key = 'Time point (hr)'     # Keyword for filtering
well_id_keyword = '$WELLID'                # Keyword with well positions
```

### 2. Run the Script

```bash
python batch_analyze_flow.py
```

### 3. Check Output

The script will generate HTML files in your output directory:
```
workspace_scatter_Cells_FSC-A_SSC-A_Time_point_hr_48.html
workspace_scatter_Cells_FSC-A_SSC-A_Time_point_hr_72.html
workspace_scatter_Cells_FSC-A_SSC-A_Time_point_hr_96.html
workspace_scatter_Cells_FSC-A_SSC-A_Time_point_hr_144.html
```

Each HTML file contains:
- 96-well plate grid layout (8 rows × 12 columns)
- Interactive scatter plots with gate overlays
- Global legend showing gate colors
- Pan/zoom/save tools
- Only samples matching the specified filter value

## Configuration Guide

### Gate Path Format

Gates are specified as tuples representing the gating hierarchy:

```python
gate_path = ('root',)                        # Root-level gate
gate_path = ('root', 'Cells')               # Cells gate under root
gate_path = ('root', 'Cells', 'Singlets')   # Nested: root → Cells → Singlets
```

### Filter Values

The `filter_values` list determines what keyword values to iterate over:

```python
# Time series
filter_values = ['0', '24', '48', '72']

# Treatment groups
filter_values = ['Control', 'Drug_A', 'Drug_B', 'Combination']

# Dose response
filter_values = ['0.1', '1.0', '10', '100']

# Plate replicates
filter_values = ['Plate_1', 'Plate_2', 'Plate_3']
```

### Keyword Filter Key

This is the FlowJo keyword that contains the values you want to filter by:

```python
keyword_filter_key = 'Time point (hr)'  # For time series
keyword_filter_key = 'Treatment'        # For treatment groups
keyword_filter_key = 'Concentration'    # For dose response
keyword_filter_key = 'Plate Name'       # For replicate plates
```

**Finding your keywords:** Run `analyze_flow.py` in interactive mode once to see all available keywords:
```bash
python analyze_flow.py --wsp workspace.wsp --fcs_dir ./fcs --interactive
```
Then select "Filter by keyword? yes" to see the complete list.

### Well ID Keyword

Specifies which keyword contains well position information (e.g., "A01", "H12"):

```python
well_id_keyword = '$WELLID'  # Common FlowJo keyword
well_id_keyword = 'Well ID'  # Alternative format
well_id_keyword = 'WellID'   # Another variant
```

The script uses fuzzy matching, so `$WELLID`, `Well ID`, and `wellid` will all work.

### Channels

Specify which parameters to plot:

```python
# Forward/side scatter
x_channel = 'FSC-A'
y_channel = 'SSC-A'

# Fluorescence channels
x_channel = 'FITC-A'
y_channel = 'PE-A'

# Custom markers
x_channel = 'CD4-APC'
y_channel = 'CD8-PE-Cy7'
```

**Note:** Channel names should match exactly as they appear in your FCS files.

## Use Cases

### Example 1: Time Series Analysis

**Scenario:** You have samples collected at 4 time points (48, 72, 96, 144 hours) across a 96-well plate.

```python
filter_values = ['48', '72', '96', '144']
keyword_filter_key = 'Time point (hr)'
```

**Output:** 4 HTML files, each showing all wells for one time point.

### Example 2: Treatment Groups

**Scenario:** You have 3 treatment conditions (Control, Drug, Combo) with 32 samples each.

```python
filter_values = ['Control', 'Drug', 'Combo']
keyword_filter_key = 'Treatment'
```

**Output:** 3 HTML files, each showing all samples for one treatment.

### Example 3: Replicate Plates

**Scenario:** You ran 3 replicate plates and want to visualize each separately.

```python
filter_values = ['Plate1', 'Plate2', 'Plate3']
keyword_filter_key = 'Plate'
```

**Output:** 3 HTML files, one per plate.

### Example 4: Multi-Gate Visualization

**Scenario:** You want to visualize multiple gates on the same plot (e.g., show both "Cells" and "Singlets" gates on ungated data).

Modify the `gates_to_visualize` list:

```python
'plot_configs': [
    {
        'gate_path': ('root',),
        'gate_name': 'Ungated',  # Show ungated data
        'plot_type': 'scatter',
        'parameters': ['FSC-A', 'SSC-A'],
        'show_gates': True,
        'gates_to_visualize': ['Cells', 'Singlets']  # Show both gates
    }
]
```

Each gate will be rendered in a different color (8-color palette cycles).

## Output Format

### HTML File Structure

Each generated HTML file contains:

1. **Header**: Experiment name and info button
2. **Global Legend**: Shows all gates with their colors
3. **Plot Grid**: 96-well layout (8 rows × 12 columns)
   - Row labels: A-H
   - Column labels: 1-12
   - Empty wells show "No Data"
4. **Interactive Tools**: Pan, zoom, box zoom, reset, save
5. **Modal**: Click "ℹ️ Info" for configuration details

### Plot Features

- **Log scale axes**: Both X and Y use log scale (10^0 to 10^5)
- **Gate overlays**: Semi-transparent filled polygons with colored borders
- **Downsampling**: Automatically samples 10,000 points if dataset is larger
- **Well ID in title**: Each plot shows well position and sample name

### File Naming Convention

```
{workspace}_{plot_type}_{gate}_{channels}_{filter_key}_{filter_value}.html
```

Example:
```
experiment_scatter_Cells_FSC-A_SSC-A_Time_point_hr_48.html
```

## Troubleshooting

### Error: "Workspace file not found"

**Problem:** The `wsp_path` is incorrect or file doesn't exist.

**Solution:**
- Check the file path is absolute (starts with `/` on Unix or `C:\` on Windows)
- Verify the file exists: `ls /path/to/workspace.wsp` (Unix) or `dir C:\path\to\workspace.wsp` (Windows)

### Error: "No samples found in selected groups"

**Problem:** The workspace doesn't contain a sample group named "All Samples", or samples aren't loaded.

**Solution:**
- Open your workspace in FlowJo and check the group names
- Update the `'selected_groups'` list in the script:
  ```python
  'selected_groups': ['Your Group Name'],
  ```

### Error: "No well IDs could be extracted"

**Problem:** The `well_id_keyword` doesn't exist or doesn't contain valid well positions.

**Solution:**
- Run interactive mode to see available keywords:
  ```bash
  python analyze_flow.py --wsp workspace.wsp --fcs_dir ./fcs --interactive
  ```
- Update `well_id_keyword` to match the correct keyword name
- Verify the keyword values are in format "A01", "B12", etc.

### Error: "stat: path should be string, bytes, os.PathLike or integer, not NoneType"

**Problem:** The `output_path` cannot be `None` (this was fixed in the provided script).

**Solution:** Use the corrected `batch_analyze_flow.py` script that builds proper output filenames.

### Warning: "No values >= 10 for channel"

**Problem:** The channel data contains only very small values (< 10), which are filtered out.

**Solution:**
- Check if you're using the correct channel name
- Verify the data is in raw (untransformed) format
- This is normal for some channels that legitimately have low values

### Gates not rendering

**Problem:** Gates don't appear on scatter plots.

**Possible causes:**
1. **Wrong channels**: Gate is defined on different channels than plot
   - **Solution**: Check gate dimensions match X and Y channels
2. **Wrong gate name**: Typo in `gate_name` or `gates_to_visualize`
   - **Solution**: Run interactive mode to see available gate names
3. **Gate not 2D**: Gate is 1D (range/interval gate) not 2D (polygon gate)
   - **Solution**: Only polygon/rectangle gates can be visualized on scatter plots

## Advanced Customization

### Modifying Plot Appearance

Edit the `_plot_scatter()` method in `analyze_flow.py`:

```python
# Change point color and size
p.scatter(x=x_channel, y=y_channel, source=source,
          size=3, alpha=0.5, color="darkblue")  # Default: size=2, color="navy"

# Change axis ranges
p = figure(..., x_range=(10, 1e4), y_range=(10, 1e4))  # Default: (1, 1e5)
```

### Adding Custom Statistics

Modify the selections dictionary to include statistics:

```python
selections = {
    # ...
    'show_keywords': True,
    'show_statistics': True,
    'selected_keywords_to_show': ['count', 'median', 'mean'],
    # ...
}
```

This will display event count, median, and mean in the bottom-right of each plot.

### Changing Gate Colors

Edit the `gate_colors` palette in `analyze_flow.py` (lines 1774-1783):

```python
gate_colors = [
    ("#FF0000", "#CC0000"),  # Red
    ("#00FF00", "#00CC00"),  # Green
    ("#0000FF", "#0000CC"),  # Blue
    # Add more (fill_color, line_color) tuples
]
```

### Processing Histogram Instead of Scatter

Modify `plot_configs` in the batch script:

```python
'plot_configs': [
    {
        'gate_path': ('root', 'Cells'),
        'gate_name': 'Cells',
        'plot_type': 'histogram',           # Change to histogram
        'parameters': ['FITC-A'],           # Single channel
        'statistic': 'median',              # 'median' or 'mean'
        'scale': 'log',                     # 'linear' or 'log'
        'show_gates': False                 # Gates not applicable for histograms
    }
]
```

## Performance Considerations

### Large Datasets

- **Downsampling**: Scatter plots automatically downsample to 10,000 points
- **Gate extraction**: Extracting gates for many samples can be slow
- **HTML size**: Files with 96 plots can be large (5-20 MB)

**Tips for better performance:**
- Process one filter value at a time for very large datasets
- Use `max_points` parameter to reduce downsampling threshold
- Close other applications to free memory

### Parallel Processing

The current script processes filter values sequentially. For faster processing, you can modify it to use Python's `multiprocessing`:

```python
from multiprocessing import Pool

def process_filter_value(filter_value):
    # ... processing code ...
    return (filter_value, output_filename, success)

with Pool(processes=4) as pool:
    results = pool.map(process_filter_value, filter_values)
```

## License

MIT License - Feel free to modify and distribute

## Contributing

Contributions are welcome! Areas for improvement:
- Additional plot types (contour plots, heatmaps)
- Support for more gate types (ellipse, quadrant)
- Performance optimization for large datasets
- Unit tests and CI/CD
- Documentation improvements

## Citation

If you use these scripts in your research, please cite:
- FlowKit: https://github.com/whitews/FlowKit
- Bokeh: https://bokeh.org

## Support

For issues or questions:
1. Check this README thoroughly
2. Review error messages in the console output
3. Try running interactive mode first: `python analyze_flow.py --wsp workspace.wsp --fcs_dir ./fcs --interactive`
4. Open an issue on GitHub with:
   - Python version (`python --version`)
   - Package versions (`pip list | grep -E "flowkit|bokeh"`)
   - Error message and full traceback
   - Example workspace file structure (if possible)

## Changelog

### v1.0.0 (2025-12-10)
- Initial release
- Batch processing support
- Interactive scatter plots with gate visualization
- 96-well plate layout
- Global legend
- Multi-color gate rendering
- Time series and treatment group analysis
