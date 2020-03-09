# -*- coding: utf-8 -*-
"""
Created on Sun Jan  5 09:59:29 2020

@author: yoelr
"""

from graphviz import Digraph
from IPython import display

def make_digraph(units, streams, format='svg',
                 maxiter='500000', **graph_attrs):
    """Return digraph of units and streams."""
    # Create a digraph and set direction left to right
    f = Digraph(format=format)
    f.attr(rankdir='LR', **graph_attrs)
    # Set up unit nodes
    UD = {}  # Contains full description (ID and line) by ID
    for u in units:
        u._load_stream_links()
        graphics = u._graphics
        if not graphics.in_system:
            continue  # Ignore Unit
        
        # Initialize graphics and make Unit node with attributes
        Type = graphics.node_function(u) or u.line
        name = u.ID + '\n' + Type
        f.attr('node', **u._graphics.node)
        f.node(name)
        UD[u] = name
        
    keys = UD.keys()

    # Set attributes for graph and streams
    f.attr('node', shape='rarrow', fillcolor='#79dae8',
           style='filled', orientation='0', width='0.6',
           height='0.6', color='black', peripheries='1')
    f.attr('graph', splines='normal', overlap='orthoyx',
           outputorder='edgesfirst', nodesep='0.15', maxiter='1000000')
    f.attr('edge', dir='foward')
    for s in streams:
        if not s: continue  # Ignore stream

        oU = s._source
        if oU:
            oi = oU._outs.index(s) 
        
        dU = s._sink
        if dU:
            di = dU._ins.index(s)
        
        # Make stream nodes / unit-stream edges / unit-unit edges
        if oU not in keys and dU not in keys: pass
            # Stream is not attached to anything
        elif oU not in keys:
            # Feed stream case
            f.node(s.ID)
            edge_in = dU._graphics.edge_in
            f.attr('edge', arrowtail='none', arrowhead='none',
                   tailport='e', **edge_in[di])
            f.edge(s.ID, UD[dU])
        elif dU not in keys:
            # Product stream case
            f.node(s.ID)
            edge_out = oU._graphics.edge_out
            f.attr('edge', arrowtail='none', arrowhead='none',
                   headport='w', **edge_out[oi])
            f.edge(UD[oU], s.ID)
        else:
            # Process stream case
            edge_in = dU._graphics.edge_in
            edge_out = oU._graphics.edge_out
            f.attr('edge', arrowtail='none', arrowhead='normal',
                   **edge_in[di], **edge_out[oi])
            f.edge(UD[oU], UD[dU], label=s.ID)
    return f

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