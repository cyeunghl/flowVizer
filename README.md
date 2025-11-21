# flowVizer

Python-based script to plot FlowJo workspace (.wsp) files in a 96-well plate layout with interactive HTML visualizations.

## Overview

flowVizer analyzes FlowJo workspaces and associated FCS files to generate interactive HTML plots organized in a 96-well plate grid. It leverages `flowkit` for parsing flow cytometry data and `bokeh` for creating interactive visualizations. All visualizations use raw (untransformed) data values for accuracy.

## Features

- **Batch Plot Generation**: Generate histograms or scatter plots across all samples in a workspace
- **96-Well Plate Layout**: Automatically organizes plots in a 8×12 grid matching plate positions
- **Interactive HTML Output**: Pan, zoom, and explore plots in your web browser
- **Flexible Well ID Extraction**: Auto-detect well IDs from keywords, filenames, or manual selection
- **Gate Visualization**: Display gate boundaries on scatter plots with olive green polygons
- **Multiple Plot Types**:
  - **Histograms**: Density plots (KDE) with linear or log scale
  - **Scatter Plots**: 2D parameter visualization with optional gate overlays
- **Statistics Display**: Show event counts, median, and/or mean values on plots
- **Keyword Filtering**: Filter samples by metadata keywords (e.g., time points)
- **Multiple Plot Configurations**: Generate multiple plot types in separate tabs

## Installation

### Dependencies

Install the required Python packages:

```bash
pip install flowkit bokeh pandas numpy scipy
```

Optional (for gate filtering):
```bash
pip install matplotlib
```

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
- `--interactive` (required): Run in interactive plotting mode
- `--output` (optional): Custom output HTML filename. If not specified, auto-generates based on plot configuration

### Examples

```bash
# Basic usage with auto-generated output filename
python analyze_flow.py --wsp workspace.wsp --fcs_dir ./fcs_files --interactive

# Custom output filename
python analyze_flow.py --wsp workspace.wsp --fcs_dir ./fcs_files --interactive --output my_analysis.html

# FCS files in same directory as workspace
python analyze_flow.py --wsp /path/to/workspace.wsp --interactive
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
- **Gate Visualization** (optional):
  - Shows available gates matching the selected channels
  - Choose to display gate boundaries
  - Selected gate (used for filtering) is automatically plotted in olive green
  - Additional gates can be displayed if enabled

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

- **Selected Gate**: Automatically plotted on scatter plots in olive green
- **Additional Gates**: Can be enabled to show all gates matching the plot channels
- **Gate Matching**: Only gates matching the X/Y channel dimensions are displayed
- **Polygon Rendering**: Gates displayed as filled polygons with semi-transparent fill

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
- **Value Filtering**: Values < 10 filtered out for histograms (prevents log transform issues)
- **Downsampling**: Scatter plots automatically downsample to 10,000 points if dataset is larger

### Well ID Parsing

- Fuzzy matching for keyword names (e.g., "WellID", "Well ID")
- Normalizes well formats (A1 = A01)
- Validates extraction across all samples before proceeding

### Gate Extraction

- Supports polygon and rectangle gates
- Automatically handles dimension swapping
- Closes polygon loops for proper rendering

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

**Empty plots**:
- Check that gate contains events
- Verify channel names match workspace channels
- Check data filtering (values < 10 for histograms)

**Import errors**:
- Ensure all dependencies are installed: `pip install flowkit bokeh pandas numpy scipy`
- Check Python version (3.7+)

## License

[Add your license information here]

## Contact

If you have any questions or suggestions for improvement, feel free to reach out at clarencehlyeung at gmail dot com.
