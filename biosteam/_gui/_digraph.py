# -*- coding: utf-8 -*-
"""
Created on Sun Jan  5 09:59:29 2020

@author: yoelr
"""
from graphviz import Digraph
from IPython import display
from collections import namedtuple

__all__ = ('new_digraph',
           'make_digraph',
           'update_digraph',
           'save_digraph',
           'get_connections')

Connection = namedtuple('Connection', 
                        ('source', 'source_index', 'stream', 'sink_index', 'sink'),
                        module=__name__)

def new_digraph(format='svg', **graph_attrs):
    # Create a digraph and set direction left to right
    f = Digraph(format=format)
    f.attr(rankdir='LR', **graph_attrs)
    return f

def make_digraph(units, streams, **graph_attrs):
    """Return digraph of units and streams."""
    f = new_digraph()
    connections = get_connections(streams)
    update_digraph(f, units, connections)
    return f

def update_digraph(f: Digraph, units, connections):
    """Update digraph of units and streams."""
    # Set up unit nodes
    unit_names = {}  # Contains full description (ID and line) by unit
    for u in units:
        try: u._load_stream_links()
        except: pass
        graphics = u._graphics
        
        # Initialize graphics and make Unit node with attributes
        node = graphics.get_node_taylored_to_unit(u)
        f.node(**node)
        unit_names[u] = node['name']

    # Set attributes for graph and streams
    f.attr('node', shape='rarrow', fillcolor='#79dae8',
           style='filled', orientation='0', width='0.6',
           height='0.6', color='black', peripheries='1')
    f.attr('graph', splines='normal', overlap='orthoyx',
           outputorder='edgesfirst', nodesep='0.15', maxiter='1000000')
    f.attr('edge', dir='foward')
    for connection in connections:
        add_connection(f, unit_names, connection)

def get_stream_connection(stream):
    source = stream._source
    source_index = source._outs.index(stream) if source else None
    sink = stream._sink
    sink_index = sink._ins.index(stream) if sink else None
    return Connection(source, source_index, stream, sink_index, sink)

def get_connections(streams):
    return {get_stream_connection(s) for s in streams if s and (s._source or s._sink)}

def add_connection(f: Digraph, unit_names, connection):
    source, source_index, stream, sink_index, sink = connection
    has_source = source in unit_names
    has_sink = sink in unit_names
    # Make stream nodes / unit-stream edges / unit-unit edges
    if has_sink and not has_source:
        # Feed stream case
        f.node(stream.ID)
        edge_in = sink._graphics.edge_in
        f.attr('edge', arrowtail='none', arrowhead='none',
               tailport='e', **edge_in[sink_index])
        f.edge(stream.ID, unit_names[sink])
    elif has_source and not has_sink:
        # Product stream case
        f.node(stream.ID)
        edge_out = source._graphics.edge_out
        f.attr('edge', arrowtail='none', arrowhead='none',
               headport='w', **edge_out[source_index])
        f.edge(unit_names[source], stream.ID)
    elif has_sink and has_source:
        # Process stream case
        edge_in = sink._graphics.edge_in
        edge_out = source._graphics.edge_out
        f.attr('edge', arrowtail='none', arrowhead='normal',
               **edge_in[sink_index], **edge_out[source_index])
        f.edge(unit_names[source], unit_names[sink], label=stream.ID)

def save_digraph(digraph, file, format):
    if not file:
        if format == 'svg':
            x = display.SVG(digraph.pipe(format=format))
        else:
            x = display.Image(digraph.pipe(format='png'))
        display.display(x)
    else:
        if '.' not in file:
            file += '.' + format
        img = digraph.pipe(format=format)
        f = open(file, 'wb')
        f.write(img)
        f.close()
    