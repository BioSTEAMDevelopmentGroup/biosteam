# -*- coding: utf-8 -*-
"""
Created on Sat Aug 18 14:36:58 2018

@author: yoelr
"""
from biosteam.utils import colors

__all__ = ('Graphics',)

# %% Base class for unit graphics
    
class Graphics:
    """Create a Graphics object that contains specifications for 
    Graphviz node and edge styles."""
    __slots__ = ('node', 'edge_in', 'edge_out', 'taylor_node_to_unit')
    
    def __init__(self, edge_in, edge_out, node, taylor_node_to_unit=None):
        # [dict] Input stream edge settings
        self.edge_in = edge_in
        
        # [dict] Output stream edge settings
        self.edge_out = edge_out
        
        #: [dict] Node settings
        self.node = node
        
        # [function(node, unit)] Taylor node to unit.
        self.taylor_node_to_unit = taylor_node_to_unit
    
    @classmethod
    def box(cls, N_ins, N_outs):
        edge_in = [{'headport': 'c'} for i in range(N_ins)]
        edge_out = [{'tailport': 'c'} for i in range(N_outs)]
        return cls(edge_in, edge_out, box_node)
    
    def get_node_taylored_to_unit(self, unit):
        """Return node taylored to unit specifications"""
        node = self.node
        node['name'] = unit.ID + '\n' + unit.line
        taylor_node_to_unit = self.taylor_node_to_unit
        if taylor_node_to_unit:
            taylor_node_to_unit(node, unit)
        return node
        
    def __repr__(self):
        return f'{type(self).__name__}(node={self.node}, edge_in={self.edge_in}, edge_out={self.edge_out})'


# %% Graphics components

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

box_graphics = Graphics(single_edge_in, single_edge_out, box_node)


# %% All graphics objects used in BioSTEAM

# Create mixer graphics
node = box_node.copy()
node['shape'] = 'triangle'
node['orientation'] = '270'
mixer_graphics = Graphics(6 * single_edge_in, right_edge_out, node)

# Create splitter graphics
node = box_node.copy()
node['shape'] = 'triangle'
node['orientation'] = '90'
node['fillcolor'] = "#bfbfbf:white"
splitter_graphics = Graphics(left_edge_in, 6 * single_edge_out, node)

# Create distillation column graphics
node = box_node.copy()
node['width'] = '1'
node['height'] = '1.2'
vertical_column_graphics = Graphics(single_edge_in, top_bottom_edge_out, node)

# Create flash column graphics
node = node.copy()
node['height'] = '1.1'
vertical_vessel_graphics = Graphics(single_edge_in, top_bottom_edge_out, node)

# Single stream heat exchanger node
node = box_node.copy()
node['shape'] = 'circle'
node['color'] = 'none'
node['margin'] = '0'
def taylor_utility_heat_exchanger_node(node, unit):
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

utility_heat_exchanger_graphics = Graphics(single_edge_in, single_edge_out, node,
                                           taylor_utility_heat_exchanger_node)

# Process heat exchanger network
node = node.copy()
node['shape'] = 'circle'
node['color'] = 'none'
node['margin'] = '0'
node['gradientangle'] = '90'
node['fillcolor'] = '#cfecf0:#fad6d8'
def taylor_process_heat_exchanger_node(node, unit):
    node['name'] = unit.ID + "\n Heat exchanger"

process_heat_exchanger_graphics = Graphics(2 * single_edge_in, 2 *single_edge_out, node,
                                           taylor_process_heat_exchanger_node)

# Process specification graphics
orange = colors.orange_tint.tint(50)
orange_tint = orange.tint(75)
node = box_node.copy()
node['fillcolor'] = orange_tint.HEX + ':' + orange.HEX
node['shape'] = 'note'
node['margin'] = '0.2'
def taylor_process_specification_node(node, unit):
    node['name'] = (f"{unit.ID} - {unit.description}\n"
                    f"{unit.line}")

process_specification_graphics = Graphics(single_edge_in, single_edge_out, node,
                                          taylor_process_specification_node)

# System unit for creating diagrams
node = box_node.copy()
node['peripheries'] = '2'    
system_unit = Graphics(multi_edge_in, multi_edge_out, node)

node = box_node.copy()
node['fillcolor'] = 'white:#79dae8'
stream_unit = Graphics(multi_edge_in, multi_edge_out, node)


node = box_node.copy()
def taylor_junction_node(node, unit):
    if not any(unit._get_streams()):
        node['fontsize'] = '18'
        node['shape'] = 'plaintext'
        node['fillcolor'] = 'none'
    else:
        node['width'] = '0.1'
        node['shape'] = 'point'
        node['color'] = node['fillcolor'] = 'black'

junction_graphics = Graphics(single_edge_in, single_edge_out, node,
                             taylor_junction_node)