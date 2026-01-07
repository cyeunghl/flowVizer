# Changelog

All notable changes to flowViz will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- **Contour Plot Visualization**: New plot type option for creating iso-density contours similar to traditional flow cytometry plots
  - Uses 2D Kernel Density Estimation (KDE) with scipy.stats.gaussian_kde
  - Renders 6 contour levels at percentiles [10, 25, 50, 75, 90, 95]
  - Supports quadrant gate dividers and polygon gate overlays
  - Works in log space for proper log-scale visualization
  - Accessible as option 3 in interactive plot type selection
- **Smart Annotation Positioning**: Annotations now automatically positioned in top-left corner for better visibility
  - Applies to all plot types: histograms, scatter plots, and contour plots
  - No longer obscures data in dense regions
  - Consistent positioning across all visualization types
- **Quadrant Gate XML Parsing**: Direct extraction of quadrant gate dividers from FlowJo workspace XML
  - Parses .wsp files to extract min/max boundaries for quadrant gates
  - Correctly handles raw data space values (no transformation needed)
  - Supports complex quadrant gate hierarchies
  - Helper module `quadrant_xml_parser.py` for XML processing

### Changed
- Annotation positioning moved from bottom-right/center-right to top-left for improved visibility
- Histogram annotation positioning updated to use consistent calculation (2% from left, 98% from top)
- Scatter plot annotation positioning updated to use consistent calculation (2% from left, 98% from top)
- Plot type selection now offers 3 options instead of 2 (Histogram, Scatter, Contour)

### Fixed
- Quadrant gate dividers now render correctly at precise boundary positions
- XML namespace handling improved for FlowJo workspace files
- Gate boundary extraction works with RectangleGate objects lacking public min/max attributes

### Technical Details
- New method `_plot_contour()` in `FlowAnalyzer` class (283 lines)
- Uses matplotlib for contour calculation (extraction only, no display)
- Maintains Bokeh for all rendering (no matplotlib rendering)
- Compatible with existing gate visualization infrastructure
- Approximately 235 lines of code added/modified across 5 locations

## Previous Versions

See git history for details on earlier releases.

---

## Migration Guide

### For Users
- No breaking changes - all existing functionality preserved
- Contour plots are a new option, existing plots work exactly as before
- Annotations will appear in top-left instead of bottom-right (improved visibility)

### For Developers
- New dependency: `matplotlib` (required for contour plot generation)
- New method: `_plot_contour()` with same signature as `_plot_scatter()`
- Plot type now accepts `'contour'` in addition to `'histogram'` and `'scatter'`
- Annotation positioning logic updated in both `_plot_histogram()` and `_plot_scatter()`
