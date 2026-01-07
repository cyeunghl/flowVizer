#!/usr/bin/env python3
"""
Quadrant Gate XML Parser for FlowJo Workspace Files

This module provides helper functions to parse quadrant gate boundaries directly
from FlowJo workspace (.wsp) XML files. This is necessary because FlowKit's
RectangleGate objects do not expose min/max boundary attributes through the
public API.

The parser extracts gate definitions from the GatingML structure within .wsp files
and returns divider positions in raw data space, which can be directly used for
visualization on log-scale axes.

Author: flowViz Contributors
License: MIT
"""

import xml.etree.ElementTree as ET
import re


def infer_quadrant_dividers_from_xml(workspace_file, sample_id, gate_path, x_channel, y_channel, quadrant_gates):
    """
    Parse the .wsp XML file to extract quadrant gate boundaries.

    Args:
        workspace_file (str): Path to the .wsp workspace file
        sample_id (str): Sample ID
        gate_path (tuple): Path tuple for the quadrant gates
        x_channel (str): X-axis channel name
        y_channel (str): Y-axis channel name
        quadrant_gates (dict): Dict of quadrant gates {q_num: {'name': gate_name, ...}}

    Returns:
        list: List of divider dicts with 'dimension', 'value', 'orientation' keys
        None: If extraction fails
    """
    print(f"    DEBUG: Parsing .wsp XML file to extract quadrant gate boundaries...")
    print(f"    DEBUG: Workspace file: {workspace_file}")

    # Parse the XML
    tree = ET.parse(workspace_file)
    root = tree.getroot()

    # Define XML namespaces
    ns = {
        'gating': 'http://www.isac-net.org/std/Gating-ML/v2.0/gating',
        'data-type': 'http://www.isac-net.org/std/Gating-ML/v2.0/datatypes',
        'transforms': 'http://www.isac-net.org/std/Gating-ML/v2.0/transformations'
    }

    # Extract gate boundaries from XML
    gate_boundaries = {}

    for q_num, q_info in quadrant_gates.items():
        gate_name = q_info['name']
        print(f"    DEBUG: Looking for Q{q_num} ('{gate_name}') in XML...")

        # Search for Population elements with matching names
        for population in root.findall('.//Population[@name]'):
            pop_name = population.get('name')

            # Check if this population matches our gate name
            if pop_name == gate_name:
                print(f"    DEBUG: Found Population: {pop_name}")

                # Find the RectangleGate within this Population
                rect_gate = population.find('.//gating:RectangleGate', ns)
                if rect_gate is not None:
                    gate_id = rect_gate.get('{http://www.isac-net.org/std/Gating-ML/v2.0/gating}id')
                    print(f"    DEBUG: Found RectangleGate with id: {gate_id}")

                    # Extract dimensions
                    dimensions = rect_gate.findall('gating:dimension', ns)

                    bounds = {}
                    for dim in dimensions:
                        dim_min = dim.get('{http://www.isac-net.org/std/Gating-ML/v2.0/gating}min')
                        dim_max = dim.get('{http://www.isac-net.org/std/Gating-ML/v2.0/gating}max')

                        # Get the dimension name
                        fcs_dim = dim.find('data-type:fcs-dimension', ns)
                        if fcs_dim is not None:
                            dim_name = fcs_dim.get('{http://www.isac-net.org/std/Gating-ML/v2.0/datatypes}name')

                            if dim_name:
                                print(f"           Dimension: {dim_name}, min={dim_min}, max={dim_max}")

                                if dim_name == x_channel:
                                    bounds['x_min'] = float(dim_min) if dim_min else None
                                    bounds['x_max'] = float(dim_max) if dim_max else None
                                elif dim_name == y_channel:
                                    bounds['y_min'] = float(dim_min) if dim_min else None
                                    bounds['y_max'] = float(dim_max) if dim_max else None

                    # Check if we got all boundaries (some might be None if unbounded)
                    if any(v is not None for v in bounds.values()):
                        gate_boundaries[q_num] = bounds
                        print(f"           Got boundaries: {bounds}")
                    break

    if len(gate_boundaries) < 2:
        print(f"    DEBUG: Could not extract enough gate boundaries from XML ({len(gate_boundaries)})")
        return None

    print(f"    DEBUG: Successfully extracted {len(gate_boundaries)} gate boundaries from XML")
    print(f"    DEBUG: Raw boundaries (from XML): {gate_boundaries}")

    # These values are in RAW data space, not display space
    # We can use them directly for plotting on raw-scale axes
    quadrant_data = gate_boundaries

    # Infer divider positions from quadrant boundaries
    # Q1: FITC-A- (X-), APC-Vio770-A+ (Y+) -> X < divider_x, Y > divider_y
    # Q2: FITC-A+ (X+), APC-Vio770-A+ (Y+) -> X > divider_x, Y > divider_y
    # Q3: FITC-A+ (X+), APC-Vio770-A- (Y-) -> X > divider_x, Y < divider_y
    # Q4: FITC-A- (X-), APC-Vio770-A- (Y-) -> X < divider_x, Y < divider_y

    dividers = []

    # X divider: boundary between Q1/Q4 (X-) and Q2/Q3 (X+)
    # Should be at the max X of Q1/Q4 or min X of Q2/Q3
    x_negative_max = []  # Q1, Q4
    x_positive_min = []  # Q2, Q3

    for q_num, data in quadrant_data.items():
        if q_num in [1, 4]:  # Negative X
            if data.get('x_max') is not None:
                x_negative_max.append(data['x_max'])
        elif q_num in [2, 3]:  # Positive X
            if data.get('x_min') is not None:
                x_positive_min.append(data['x_min'])

    # For quadrant gates, the divider should be exactly at the boundary
    # Use the max of Q1/Q4 or min of Q2/Q3 (they should be the same)
    if x_negative_max or x_positive_min:
        if x_negative_max and x_positive_min:
            # They should be equal or very close - take the average
            x_divider = (max(x_negative_max) + min(x_positive_min)) / 2.0
        elif x_negative_max:
            x_divider = max(x_negative_max)
        else:
            x_divider = min(x_positive_min)

        dividers.append({
            'dimension': x_channel,
            'value': float(x_divider),
            'orientation': 'vertical'
        })
        print(f"    DEBUG: Inferred X divider at {x_divider:.2f} for {x_channel}")

    # Y divider: boundary between Q1/Q2 (Y+) and Q3/Q4 (Y-)
    # Should be at the max Y of Q3/Q4 or min Y of Q1/Q2
    y_negative_max = []  # Q3, Q4
    y_positive_min = []  # Q1, Q2

    for q_num, data in quadrant_data.items():
        if q_num in [3, 4]:  # Negative Y
            if data.get('y_max') is not None:
                y_negative_max.append(data['y_max'])
        elif q_num in [1, 2]:  # Positive Y
            if data.get('y_min') is not None:
                y_positive_min.append(data['y_min'])

    if y_negative_max or y_positive_min:
        if y_negative_max and y_positive_min:
            # They should be equal or very close - take the average
            y_divider = (max(y_negative_max) + min(y_positive_min)) / 2.0
        elif y_negative_max:
            y_divider = max(y_negative_max)
        else:
            y_divider = min(y_positive_min)

        dividers.append({
            'dimension': y_channel,
            'value': float(y_divider),
            'orientation': 'horizontal'
        })
        print(f"    DEBUG: Inferred Y divider at {y_divider:.2f} for {y_channel}")

    if len(dividers) > 0:
        print(f"    DEBUG: Final dividers (raw space): {dividers}")
        return dividers
    else:
        return None
