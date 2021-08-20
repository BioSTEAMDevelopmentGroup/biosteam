# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from warnings import warn
from biosteam.exceptions import GraphicsWarning
from .utils import colors

__all__ = ('UnitGraphics',
           'box_graphics',
           'mixer_graphics',
           'splitter_graphics',
           'vertical_column_graphics',
           'vertical_vessel_graphics',
           'utility_heat_exchanger_graphics',
           'process_heat_exchanger_graphics',
           'process_specification_graphics',
           'system_unit',
           'stream_unit',
           'junction_graphics')

# %% Base class for unit graphics
    
class UnitGraphics:
    """Create a UnitGraphics object that contains specifications for 
    Graphviz node and edge styles."""
    __slots__ = ('node', 'edge_in', 'edge_out', 'tailor_node_to_unit')
    
    def __init__(self, edge_in, edge_out, node, tailor_node_to_unit=None):
        # [dict] Input stream edge settings
        self.edge_in = edge_in
        
        # [dict] Output stream edge settings
        self.edge_out = edge_out
        
        #: [dict] Node settings
        self.node = node
        
        # [function(node, unit)] tailor node to unit.
        self.tailor_node_to_unit = tailor_node_to_unit
    
    def get_inlet_options(self, sink, sink_index):
        edge_in = self.edge_in
        try:
            options = edge_in[sink_index]
        except IndexError:
            if sink._ins_size_is_fixed: 
                N_inlets = len(edge_in)
                warn(f'inlet #{sink_index} at {repr(sink)} missing graphics options; '
                     f'expected at most {N_inlets} inlet' + ('' if N_inlets == 1 else 's'),
                      GraphicsWarning)
            options = {'headport': 'c'}
        return options
    
    def get_outlet_options(self, source, source_index):
        edge_out = self.edge_out
        try:
            options = edge_out[source_index]
        except IndexError:
            if source._outs_size_is_fixed: 
                N_outlets = len(edge_out)
                warn(f'outlet #{source_index} at {repr(source)} missing graphics options; '
                     f'expected at most {N_outlets} outlet' + ('' if N_outlets == 1 else 's'),
                      GraphicsWarning)
            options = {'tailport': 'c'}
        return options
    
    @classmethod
    def box(cls, N_ins, N_outs):
        edge_in = [{'headport': 'c'} for i in range(N_ins)]
        edge_out = [{'tailport': 'c'} for i in range(N_outs)]
        return cls(edge_in, edge_out, box_node)
    
    def get_minimal_node(self, unit):
        """Return minmal node (a single dot)."""
        return dict(
            name = unit.ID,
            width = '0.1',
            shape = 'oval',
            fillcolor = "#d2d2d2:white",
            style = 'filled',
        )
    
    def get_node_tailored_to_unit(self, unit): # pragma: no coverage
        """Return node tailored to unit specifications"""
        node = self.node
        node['name'] = unit.ID + '\n' + unit.line
        tailor_node_to_unit = self.tailor_node_to_unit
        if tailor_node_to_unit:
            tailor_node_to_unit(node, unit)
        return node
        
    def __repr__(self): # pragma: no coverage
        return f'{type(self).__name__}(node={self.node}, edge_in={self.edge_in}, edge_out={self.edge_out})'


# %% UnitGraphics components

single_edge_in = ({'headport': 'c'},)
single_edge_out = ({'tailport': 'c'},)
multi_edge_in = 20 * single_edge_in
multi_edge_out = 20 * single_edge_out
right_edge_out = ({'tailport': 'e'},)
left_edge_in = ({'headport': 'w'},)
top_bottom_edge_out = ({'tailport': 'n'}, {'tailport': 's'})

box_node = {'shape': 'box',
            'fillcolor': "white:#CDCDCD",
            'style': 'filled',
            'gradientangle': '0',
            'width': '0.6',
            'height': '0.6',
            'orientation': '0.0',
            'color': 'black',
            'peripheries': '1',
            'margin': 'default'}

box_graphics = UnitGraphics(single_edge_in, single_edge_out, box_node)


# %% All graphics objects used in BioSTEAM

# Create mixer graphics
node = box_node.copy()
node['shape'] = 'triangle'
node['orientation'] = '270'
mixer_graphics = UnitGraphics(multi_edge_in, right_edge_out, node)

# Create splitter graphics
node = box_node.copy()
node['shape'] = 'triangle'
node['orientation'] = '90'
node['fillcolor'] = "#bfbfbf:white"
splitter_graphics = UnitGraphics(left_edge_in, 6 * single_edge_out, node)

# Create distillation column graphics
node = box_node.copy()
node['width'] = '1'
node['height'] = '1.2'
vertical_column_graphics = UnitGraphics(single_edge_in, top_bottom_edge_out, node)

# Create flash column graphics
node = box_node.copy()
node['height'] = '1.1'
vertical_vessel_graphics = UnitGraphics(single_edge_in, top_bottom_edge_out, node)

# Mixer-Settler graphics
node = box_node.copy()
node['width'] = '1.2'
mixer_settler_graphics = UnitGraphics(multi_edge_in, top_bottom_edge_out, node)

# Single stream heat exchanger node
node = box_node.copy()
node['shape'] = 'circle'
node['color'] = 'none'
node['margin'] = '0'
def tailor_utility_heat_exchanger_node(node, unit): # pragma: no coverage
    try:
        si = unit.ins[0]
        so = unit.outs[0]
        H_in = si.H
        H_out = so.H
        if H_in > H_out:
            node['fillcolor'] = '#cfecf0'
            node['gradientangle'] = '0'
            line = 'Cooling'
        elif H_in < H_out:
            node['gradientangle'] = '0'
            node['fillcolor'] = '#fad6d8'
            line = 'Heating'
        else:
            node['gradientangle'] = '90'
            node['fillcolor'] = '#cfecf0:#fad6d8'
            line = 'Heat exchanger'
    except:
        line = 'Heat exchanger'
    node['name'] = unit.ID + "\n" + line

utility_heat_exchanger_graphics = UnitGraphics(single_edge_in, single_edge_out, node,
                                               tailor_utility_heat_exchanger_node)

# Process heat exchanger network
node = node.copy()
node['shape'] = 'circle'
node['color'] = 'none'
node['margin'] = '0'
node['gradientangle'] = '90'
node['fillcolor'] = '#cfecf0:#fad6d8'
def tailor_process_heat_exchanger_node(node, unit): # pragma: no coverage
    node['name'] = unit.ID + "\n Heat exchanger"

process_heat_exchanger_graphics = UnitGraphics(2 * single_edge_in, 2 *single_edge_out, node,
                                               tailor_process_heat_exchanger_node)

# Process specification graphics
orange = colors.orange_tint.tint(50)
orange_tint = orange.tint(75)
node = box_node.copy()
node['fillcolor'] = orange_tint.HEX + ':' + orange.HEX
node['shape'] = 'note'
node['margin'] = '0.2'
def tailor_process_specification_node(node, unit): # pragma: no coverage
    node['name'] = (f"{unit.ID} - {unit.description}\n"
                    f"{unit.line}")

process_specification_graphics = UnitGraphics(single_edge_in, single_edge_out, node,
                                              tailor_process_specification_node)

# System unit for creating diagrams
node = box_node.copy()
node['peripheries'] = '2'    
system_unit = UnitGraphics(multi_edge_in, multi_edge_out, node)

node = box_node.copy()
node['fillcolor'] = 'white:#79dae8'
stream_unit = UnitGraphics(multi_edge_in, multi_edge_out, node)


node = box_node.copy()
def tailor_junction_node(node, unit): # pragma: no coverage
    if not any(unit._ins + unit._outs):
        node['fontsize'] = '18'
        node['shape'] = 'plaintext'
        node['fillcolor'] = 'none'
    else:
        node['width'] = '0.1'
        node['shape'] = 'point'
        node['color'] = node['fillcolor'] = 'black'

junction_graphics = UnitGraphics(single_edge_in, single_edge_out, node,
                                 tailor_junction_node)