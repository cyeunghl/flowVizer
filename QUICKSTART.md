# Quick Start Guide - Batch Flow Cytometry Analysis

Get up and running in 5 minutes!

## Prerequisites

```bash
pip install flowkit bokeh pandas numpy scipy matplotlib
```

## Step 1: Prepare Your Data

You need:
1. FlowJo workspace file (`.wsp`)
2. FCS files directory
3. Keywords in your workspace (e.g., "Time point", "Treatment", "Well ID")

## Step 2: Configure the Script

Open `batch_analyze_flow.py` and edit the CONFIGURATION section:

```python
# File paths
wsp_path = "/path/to/your/workspace.wsp"
fcs_dir = "/path/to/your/fcs_files"
output_dir = "/path/to/output"

# What to iterate over
filter_values = ['48', '72', '96', '144']  # Your experimental values
keyword_filter_key = 'Time point (hr)'     # Your keyword name
```

### How to Find Your Keywords

Run the interactive tool once:

```bash
python analyze_flow.py --wsp your_workspace.wsp --fcs_dir ./fcs --interactive
```

Follow the prompts and select "Filter by keyword? yes" to see all available keywords and their values.

## Step 3: Run the Batch Script

```bash
python batch_analyze_flow.py
```

Expected output:
```
======================================================================
Batch Flow Cytometry Analysis
======================================================================

Workspace: /path/to/workspace.wsp
FCS Directory: /path/to/fcs_files
Output Directory: /path/to/output
Filter values: 48, 72, 96, 144

Loading workspace...
✓ Workspace loaded successfully

Found 96 total samples in workspace

──────────────────────────────────────────────────────────────────────
Processing Time point (hr): 48
──────────────────────────────────────────────────────────────────────

Generating scatter plots for Time point (hr) = 48...
  [1/24] Processing A01.fcs...
    ✓ Found well ID: A01 (extracted from keyword '$WELLID')
    → Loading PRE-FILTERED data (before 'Cells' gate)... ✓
    ✓ Created scatter plot
  [2/24] Processing A02.fcs...
  ...

✓ Generated: workspace_scatter_Cells_FSC-A_SSC-A_Time_point_hr_48.html

...

======================================================================
Batch Processing Complete
======================================================================

Successful: 4/4

Generated files:
  ✓ workspace_scatter_Cells_FSC-A_SSC-A_Time_point_hr_48.html
  ✓ workspace_scatter_Cells_FSC-A_SSC-A_Time_point_hr_72.html
  ✓ workspace_scatter_Cells_FSC-A_SSC-A_Time_point_hr_96.html
  ✓ workspace_scatter_Cells_FSC-A_SSC-A_Time_point_hr_144.html

✓ All filter values processed successfully!
```

## Step 4: View Your Results

Open any generated HTML file in a web browser:

```bash
open output/workspace_scatter_Cells_FSC-A_SSC-A_Time_point_hr_48.html
```

You'll see:
- 96-well plate grid layout
- Interactive scatter plots with gate overlays
- Pan/zoom/save tools
- Legend showing gate colors
- "ℹ️ Info" button for configuration details

## Common Configurations

### Time Series

```python
filter_values = ['0', '24', '48', '72']
keyword_filter_key = 'Time point (hr)'
```

### Treatment Groups

```python
filter_values = ['Control', 'Drug_A', 'Drug_B']
keyword_filter_key = 'Treatment'
```

### Dose Response

```python
filter_values = ['0.1', '1', '10', '100']
keyword_filter_key = 'Concentration (uM)'
```

## Troubleshooting

### "Workspace file not found"
- Check your file path is correct and absolute
- Use forward slashes `/` even on Windows

### "No samples found"
- Check your workspace has a group named "All Samples"
- Or change `'selected_groups': ['Your Group Name']`

### "No well IDs could be extracted"
- Run interactive mode to find the correct keyword name
- Update `well_id_keyword = 'Your Well ID Keyword'`

### Gates not rendering
- Verify gate name: `gate_name = 'YourGateName'`
- Check channels match gate dimensions
- Make sure gate is 2D (polygon/rectangle)

## Next Steps

- Read the full [README_BATCH.md](README_BATCH.md) for advanced features
- Customize plot appearance in `analyze_flow.py`
- Add statistics to plots
- Process histograms instead of scatter plots

## Getting Help

1. Check error messages in terminal
2. Review [README_BATCH.md](README_BATCH.md) troubleshooting section
3. Try interactive mode first to validate your workspace
4. Open a GitHub issue with error details

Happy analyzing! 🔬
