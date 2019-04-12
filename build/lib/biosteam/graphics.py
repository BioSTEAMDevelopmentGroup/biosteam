# -*- coding: utf-8 -*-
"""
Created on Sat Aug 18 14:36:58 2018

@author: yoelr
"""
import copy

__all__ = ('Graphics',)


# %%

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
    in_system = True #: [bool] True for Unit object to appear within a system diagram
    def __init__(self, edge_in=_edge_in,
                 edge_out=_edge_out,
                 node=_node,
                 node_function=lambda unit: None,
                 name=None):
        self.edge_in = copy.deepcopy(edge_in)
        self.edge_out = copy.deepcopy(edge_out)
        self.node = copy.deepcopy(node)
        self.node_function = node_function
        self.name = None

default_graphics = Graphics()
