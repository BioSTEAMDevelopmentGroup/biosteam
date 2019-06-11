# -*- coding: utf-8 -*-
"""
Created on Sat Aug 18 14:36:58 2018

@author: yoelr
"""
import copy

__all__ = ('Graphics',)


# %%

_node_function = lambda unit: None
_edge_in = [{'headport': 'c'} for i in range(3)]
_edge_out = [{'tailport': 'c'} for i in range(2)]
_node = {'shape': 'box',
         'fillcolor': "white:#CDCDCD",
         'style': 'filled',
         'gradientangle': '0',
         'width': '0.6',
         'height': '0.6',
         'orientation': '0.0',
         'color': 'black',
         'peripheries': '1',
         'margin': 'default'}

# %%
    
class Graphics:
    """Create a Graphics object that contains specifications for Graphviz node and edge styles"""
    def __init__(self, edge_in=None,
                 edge_out=None,
                 node=None,
                 node_function=None,
                 in_system=True):
        # [dict] Input stream edge settings
        self.edge_in = edge_in or copy.deepcopy(_edge_in)
        
        # [dict] Output stream edge settings
        self.edge_out = edge_out or copy.deepcopy(_edge_out)
        
        #: [dict] Node settings
        self.node = node or copy.copy(_node)
        
        #: [function] Is called to update node settings
        self.node_function = node_function or _node_function
        
        #: [bool] True for Unit object to appear within a system diagram
        self.in_system = in_system

    def __repr__(self):
        return 'f<{type(self).__name__}>'

default_graphics = Graphics()
