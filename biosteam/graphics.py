# -*- coding: utf-8 -*-
"""
Created on Sat Aug 18 14:36:58 2018

@author: yoelr
"""
import copy
import weakref
head_port = {'headport': 'c'}
tail_port = {'tailport': 'c'}
edge_in = [head_port for i in range(20)]
edge_out = [tail_port for i in range(20)]
node = {'shape': 'box',
        'fillcolor': "white:#CDCDCD",
        'style': 'filled',
        'gradientangle': '0',
        'width': '0.6',
        'height': '0.6',
        'orientation': '0.0',
        'color': 'black',
        'peripheries': '1',
        'margin': 'default'}


class Graphics:
    """Create a Graphics object that contains specifications for Graphviz node and edge styles"""
    _all = []
    in_system = True #: [bool] True for Unit object to appear within a system diagram
    def __init__(self, edge_in=edge_in,
                 edge_out=edge_out,
                 node=node,
                 node_function=lambda unit: None,
                 name=None):
        self.edge_in = copy.deepcopy(edge_in)
        self.edge_out = copy.deepcopy(edge_out)
        self.node = copy.deepcopy(node)
        self.node_function = node_function
        self.name = None
        Graphics._all.append(weakref.ref(self))

    def default(self):
        """Set all attribute to default attributes"""
        self.edge_in = copy.deepcopy(edge_in)
        self.edge_out = copy.deepcopy(edge_out)
        self.node = copy.deepcopy(node)
        self.node_function = lambda unit: None

    @staticmethod
    def default_all():
        """Set attributes of all objects in this class to default"""
        for ref in Graphics._all:
            if ref() is None:
                del ref
            else:
                ref().default()


default_graphics = Graphics()
