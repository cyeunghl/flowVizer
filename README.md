# flowViz

Python-based tools for FlowJo workspace analysis with interactive HTML visualizations and automated batch processing.

## Overview

flowViz analyzes FlowJo workspaces and associated FCS files to generate interactive HTML plots organized in a 96-well plate grid. It leverages `flowkit` for parsing flow cytometry data and `bokeh` for creating interactive visualizations. All visualizations use raw (untransformed) data values for accuracy.

**Three modes of operation:**
1. **Interactive Mode** (`analyze_flow.py --interactive`): CLI-guided plot generation with manual configuration
2. **Inspect Mode** (`analyze_flow.py --inspect`): Display gate hierarchy and gate information without generating plots
3. **Batch Mode** (`batch_analyze_flow.py`): Automated processing for time series, dose-response, and replicate analysis

## Quick Start

**For automated batch processing** (recommended for experiments with multiple conditions):
```bash
# See QUICKSTART.md for 5-minute setup guide
python batch_analyze_flow.py
```

**For manual/exploratory analysis:**
```bash
python analyze_flow.py --wsp workspace.wsp --fcs_dir ./fcs_files --interactive
```

**For inspecting gate structure:**
```bash
python analyze_flow.py --wsp workspace.wsp --fcs_dir ./fcs_files --inspect
```

## Documentation

- **[QUICKSTART.md](QUICKSTART.md)** - Get started in 5 minutes with batch processing
- **[README_BATCH.md](README_BATCH.md)** - Comprehensive batch processing guide
- **README.md** (this file) - Interactive mode documentation

## Features

- **Automated Batch Processing**: Process multiple experimental conditions without manual intervention (NEW!)
- **Time Series Analysis**: Iterate over time points, treatments, or other experimental variables (NEW!)
- **Multi-Gate Visualization**: Display multiple gates with distinct colors and global legend (NEW!)
- **Batch Plot Generation**: Generate histograms or scatter plots across all samples in a workspace
- **96-Well Plate Layout**: Automatically organizes plots in a 8×12 grid matching plate positions
- **Interactive HTML Output**: Pan, zoom, and explore plots in your web browser
- **Flexible Well ID Extraction**: Auto-detect well IDs from keywords, filenames, or manual selection
- **Gate Visualization**: Display gate boundaries on scatter plots with distinct colors
- **Quadrant Gate Support**: Visualize quadrant gates as divider lines
- **Multiple Plot Types**:
  - **Histograms**: Density plots (KDE) with linear or log scale
  - **Scatter Plots**: 2D parameter visualization with optional gate overlays
  - **Contour Plots**: Iso-density contours similar to traditional flow cytometry plots (NEW!)
- **Statistics Display**: Show event counts, median, and/or mean values on plots
- **Smart Annotation Positioning**: Annotations automatically positioned in top-left corner for better visibility (NEW!)
- **Keyword Filtering**: Filter samples by metadata keywords (e.g., time points)
- **Multiple Plot Configurations**: Generate multiple plot types in separate tabs

## Installation

### Dependencies

Install the required Python packages:

```bash
pip install flowkit bokeh pandas numpy scipy
```

**Required for contour plots:**
```bash
pip install matplotlib
```

**Note**: Matplotlib is required for contour plot generation (used for density calculations). It's optional if you only use histogram and scatter plots.

### Requirements

- Python 3.7+
- FlowJo workspace file (.wsp)
- Associated FCS files referenced in the workspace

## Usage

### Basic Command

```bash
python analyze_flow.py --wsp <path_to_workspace.wsp> --fcs_dir <path_to_fcs_files> --interactive
```

### Command-Line Arguments

- `--wsp` (required): Path to the FlowJo workspace (.wsp) file
- `--fcs_dir` (optional): Directory containing FCS files. If not specified, uses the workspace directory
- `--interactive` (optional): Run in interactive plotting mode (mutually exclusive with `--inspect`)
- `--inspect` (optional): Run in inspect mode to view gate hierarchy (mutually exclusive with `--interactive`)
- `--sample` (optional, inspect mode only): Specific sample ID to inspect. If not specified, uses first sample
- `--output` (optional, interactive mode only): Custom output HTML filename. If not specified, auto-generates based on plot configuration

**Note**: Either `--interactive` or `--inspect` must be specified, but not both.

### Examples

**Interactive Mode:**
```bash
# Basic usage with auto-generated output filename
python analyze_flow.py --wsp workspace.wsp --fcs_dir ./fcs_files --interactive

# Custom output filename
python analyze_flow.py --wsp workspace.wsp --fcs_dir ./fcs_files --interactive --output my_analysis.html

# FCS files in same directory as workspace
python analyze_flow.py --wsp /path/to/workspace.wsp --interactive
```

**Inspect Mode:**
```bash
# Inspect all gates in workspace (uses first sample)
python analyze_flow.py --wsp workspace.wsp --fcs_dir ./fcs_files --inspect

# Inspect gates for a specific sample
python analyze_flow.py --wsp workspace.wsp --fcs_dir ./fcs_files --inspect --sample sample1.fcs

# FCS files in same directory as workspace
python analyze_flow.py --wsp /path/to/workspace.wsp --inspect
```

## Inspect Mode

The inspect mode displays the gate hierarchy and gate information without generating plots. This is useful for:
- Understanding the gating structure before plotting
- Identifying available gates and their types
- Finding gate dimensions and channel combinations
- Debugging gate extraction issues

### Inspect Mode Output

The inspect mode displays:

1. **Gate Hierarchy**: All gates organized by hierarchical path
   - Shows gate names, types, dimensions, and vertex counts
   - Displays gates in tree structure (root → parent → child)

2. **Gate Details**:
   - **Gate Type**: polygon, rectangle, quadrant, 1D, or unknown
   - **Dimensions**: Channel names (e.g., "FSC-A × SSC-A")
   - **Vertices**: Number of vertices for polygon gates
   - **Dividers**: For quadrant gates, shows divider information

3. **Summary Statistics**:
   - Total number of gates
   - Number of gate paths
   - Gates grouped by type
   - Gates grouped by dimension count

4. **Available Channels**: Lists all channels available in the sample data

### Example Inspect Output

```
Gate Hierarchy:
================================================================================

Path: root
  Number of gates: 2

  Gate: 'Ungated'
    Type: ungated
    Dimensions: (Ungated - all channels)

  Gate: 'Cells'
    Type: polygon
    Dimensions (2D): FSC-A × SSC-A
    Vertices: 4

Path: root → Cells
  Number of gates: 1

  Gate: 'Singlets'
    Type: polygon
    Dimensions (2D): FSC-A × SSC-A
    Vertices: 4

Path: root → Cells → Singlets
  Number of gates: 4

  Gate: 'Q1: FITC-A- , APC-Vio770-A+'
    Type: quadrant
    Dimensions (2D): B1-A × R2-A

  Gate: 'Q2: FITC-A+ , APC-Vio770-A+'
    Type: quadrant
    Dimensions (2D): B1-A × R2-A

  Gate: 'Q3: FITC-A+ , APC-Vio770-A-'
    Type: quadrant
    Dimensions (2D): B1-A × R2-A

  Gate: 'Q4: FITC-A- , APC-Vio770-A-'
    Type: quadrant
    Dimensions (2D): B1-A × R2-A

Summary Statistics
================================================================================

Total gates: 7
Number of gate paths: 3

Gates by type:
  polygon: 2
  quadrant: 4
  ungated: 1

Gates by dimension count:
  Ungated: 1
  2D: 6
```

## Interactive Workflow

The interactive mode guides you through the following steps:

### 1. Sample Group Selection
- Choose which sample groups to plot (subset or all available groups)
- Useful for large workspaces with multiple experimental groups

### 2. Additional Information Display
Choose what information to display on plots:
- **Statistics**: Event count, median, mean, or combinations
- **Keywords**: Display sample metadata keywords
- **None**: No additional information

### 3. Number of Plots
- Specify how many different plot configurations to generate
- Each configuration appears in a separate tab in the output HTML

### 4. Keyword Filtering (Optional)
- Filter samples by keyword values (e.g., "Time point (hr)" = "0")
- Useful when workspace contains multiple files per well
- Lists all available keywords and their unique values

### 5. Well ID Extraction
Choose how to extract well IDs:
- **Auto**: Try keywords first, then filename
- **From keyword**: Select a specific keyword containing well IDs
- **From filename**: Extract from filename only

The script validates extraction and lists all matched well IDs. If no well IDs are found, the process stops with an error message.

**Note**: Well IDs are parsed with fuzzy matching (handles "WellID", "Well ID", etc.) and formats "A1" and "A01" are treated as equivalent.

### 6. Gate Selection
- Select gate path from available gating hierarchies
- Select specific gate name within that path
- **Use Parent Gate Option**: When selecting a gate path with children, option "0" allows you to use the parent gate for data filtering without further gating
  - Example: Select path "root → Cells → Singlets", then choose "0" to use "Singlets" gate data
  - This allows visualizing child gates (e.g., Q1-Q4 quadrants) on top of the parent population
- "Ungated" option available for all events

### 7. Plot Configuration

For each plot configuration:

#### Histogram Plots
- Select channel (e.g., "RL1-A", "B1-A")
- Choose statistic to display: Median or Mean
- Select x-axis scale: Linear or Log
  - **Log scale**: Fixed range from 0.1 to 10⁵
  - **Linear scale**: Dynamic range with outlier filtering
- Filters out values < 10 to prevent log transform issues
- Displays as smooth density curve (KDE) with olive green fill

#### Scatter Plots
- Select X-axis channel (e.g., "FSC-A")
- Select Y-axis channel (e.g., "SSC-A")
- Displays individual data points (downsampled to 10,000 points for performance)
- Log-log scale axes (1 to 10⁵) matching FlowJo defaults

#### Contour Plots (NEW!)
- Select X-axis channel (e.g., "B1-A" for FITC-A)
- Select Y-axis channel (e.g., "R2-A" for APC-Vio770-A)
- Displays iso-density contours similar to traditional flow cytometry plots
- Uses 2D Kernel Density Estimation (KDE) for smooth density calculation
- 6 contour levels at percentiles [10, 25, 50, 75, 90, 95]
- Log-log scale axes (1 to 10⁵) matching FlowJo defaults
- Perfect for visualizing cell population density distributions
- **Gate Visualization** (optional):
  - Shows available gates matching the selected channels
  - Automatically includes child gates of the selected parent gate (e.g., Q1-Q4 quadrants)
  - Choose to display gate boundaries
  - **Gate Types**:
    - **Polygon/Rectangle gates**: Displayed as filled polygons with semi-transparent fill
    - **Quadrant gates**: Displayed as dashed divider lines (vertical/horizontal)
  - Selected gate (used for filtering) is automatically plotted
  - Additional gates can be displayed if enabled
  - Can select "all" to visualize all matching gates, or specific gates by number

### 8. Output Generation
- Processes all samples matching filter criteria
- Organizes plots in 96-well plate grid (8 rows × 12 columns)
- Generates HTML file with interactive plots
- Includes summary information and processing statistics

## Output Format

### File Naming

Output filenames are auto-generated based on:
- Workspace name
- Plot type (histogram/scatter)
- Gate name
- Parameters (channels)
- Keyword filter (if applied)

**Format**: `{workspace_name}_{plot_type}_{gate_name}_{parameters}_{filter}.html`

**Example**: `phrodo011_histogram_Singlets_RL1-A_Time_point_hr_0.html`

If `--output` is specified, that filename is used instead.

### HTML Structure

- **Multiple Tabs**: Each plot configuration appears in a separate tab
- **96-Well Grid**: Plots organized by well position (A01-H12)
- **Interactive Tools**: Pan, zoom, box zoom, reset, save
- **Summary Section**: Plot configuration and processing statistics
- **Empty Wells**: Displayed as "No Data" placeholders

### Plot Features

- **Histograms**:
  - Smooth density curves (KDE) with olive green fill
  - Median/mean statistic lines (if enabled)
  - Event count display (if enabled)
  - Fixed log scale range (0.1 to 10⁵) or dynamic linear range

- **Scatter Plots**:
  - Downsampled to 10,000 points for performance
  - Gate boundaries as olive green polygons (if enabled)
  - Selected gate highlighted with thicker border
  - Statistics for both X and Y channels (if enabled)

## Advanced Features

### Gate Visualization

- **Selected Gate**: Automatically plotted on scatter plots
- **Additional Gates**: Can be enabled to show all gates matching the plot channels
- **Child Gates**: When a parent gate is selected, child gates (e.g., Q1-Q4 quadrants) are automatically discovered and available for visualization
- **Gate Matching**: Only gates matching the X/Y channel dimensions are displayed
- **Gate Types**:
  - **Polygon/Rectangle Gates**: Displayed as filled polygons with semi-transparent fill and colored borders
  - **Quadrant Gates**: Displayed as dashed divider lines spanning the full plot extent
    - Vertical dividers: Lines at specific X values
    - Horizontal dividers: Lines at specific Y values
    - No fill or clipping regions
- **Color Palette**: Distinct colors for each gate type (olive green, steel blue, peru, purple, crimson, etc.)

### Statistics Display

Choose from 7 options:
1. Event count only
2. Median only
3. Mean only
4. Median + Mean
5. Event count + Median
6. Event count + Mean
7. Event count + Median + Mean

Statistics are calculated from the actual plot data with the same filtering applied.

### Multiple Plot Configurations

Generate multiple plot types in a single session:
- Each configuration creates a separate tab
- Useful for comparing different gates, channels, or plot types
- All plots share the same sample filtering and well ID extraction settings

## Technical Details

### Data Processing

- **Raw Data**: All plots use raw (untransformed) event data
- **Outlier Filtering**: IQR-based filtering (3×IQR) for better visualization
- **Value Filtering**: Values ≤ 0 filtered out for log scale plots (prevents log transform issues)
- **Downsampling**:
  - Scatter plots: automatically downsample to 10,000 points if dataset is larger
  - Contour plots: downsample to 10,000 points for KDE calculation if dataset is larger

### Well ID Parsing

- Fuzzy matching for keyword names (e.g., "WellID", "Well ID")
- Normalizes well formats (A1 = A01)
- Validates extraction across all samples before proceeding

### Gate Extraction

- Supports multiple gate types:
  - **Polygon gates**: Extracted from vertex coordinates
  - **Rectangle gates**: Created from min/max bounds
  - **Quadrant gates**: Extracted from divider positions (X and/or Y values)
- Automatically handles dimension swapping
- Closes polygon loops for proper rendering
- Converts gate coordinates from FlowJo display space (0-1 normalized) to raw data space
- Quadrant gates return divider metadata (not vertices) for line rendering

## Troubleshooting

### Common Issues

**No well IDs found**:
- Check that keyword names match (try "auto" mode)
- Verify FCS filenames contain well IDs
- Ensure workspace is properly loaded

**Gates not displaying**:
- Verify gate dimensions match the selected X/Y channels
- Check that gate is 2D (not 1D)
- Ensure "Show gates" option is enabled for scatter plots
- For quadrant gates: Verify dividers are extracted correctly (check debug output)
- Child gates: Ensure parent gate is selected and child gates match plot channels

**Empty plots**:
- Check that gate contains events
- Verify channel names match workspace channels
- Check data filtering (values < 10 for histograms)

**Import errors**:
- Ensure all dependencies are installed: `pip install flowkit bokeh pandas numpy scipy`
- Check Python version (3.7+)

## License

[Add your license information here]

## Contributing

[Add contribution guidelines here]

## Contact

[Add contact information here]
