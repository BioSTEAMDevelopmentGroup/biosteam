# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import biosteam as bst
from graphviz import Digraph
from IPython import display
from collections import namedtuple
from thermosteam import Stream

__all__ = ('digraph_from_units',
           'digraph_from_units_and_streams',
           'digraph_from_units_and_connections',
           'update_digraph_from_units_and_connections',
           'blank_digraph',
           'minimal_digraph',
           'surface_digraph',
           'finalize_digraph',
           'update_surface_units',
           'display_digraph',
           'save_digraph',
           'get_all_connections',
           'get_stream_connection')

Connection = namedtuple('Connection', 
                        ('source', 'source_index', 'stream', 'sink_index', 'sink'),
                        module=__name__)

def has_path(obj): # pragma: no coverage
    return hasattr(obj, 'path')

def blank_digraph(format='svg', maxiter='10000000', 
                  Damping='0.995', K='0.5', **graph_attrs): # pragma: no coverage
    # Create a digraph and set direction left to right
    f = Digraph(format=format)
    f.attr(rankdir='LR', maxiter=maxiter, Damping=Damping, K=K, **graph_attrs)
    return f

def get_section_inlets_and_outlets(units, streams): # pragma: no coverage
    outs = []
    ins = []
    units = tuple(units)
    for s in streams:
        source = s._source
        sink = s._sink
        if source in units and sink not in units:
            outs.append(s)
        elif sink in units and source not in units:
            ins.append(s)
    return ins, outs

def minimal_digraph(ID, units, streams, **graph_attrs): # pragma: no coverage
    ins, outs = get_section_inlets_and_outlets(units, streams)
    product = Stream(None)
    product._ID = ''
    feed = Stream(None)
    feed._ID = ''
    feed_box = bst.units.DiagramOnlyStreamUnit('\n'.join([i.ID for i in ins]),
                                               None, feed)
    product_box = bst.units.DiagramOnlyStreamUnit('\n'.join([i.ID for i in outs]),
                                                  product, None)
    system_box = bst.units.DiagramOnlySystemUnit(ID, feed, product)
    return digraph_from_units((feed_box, system_box, product_box),
                              **graph_attrs)

def surface_digraph(path): # pragma: no coverage
    surface_units = set()  
    old_unit_connections = set()
    isa = isinstance
    Unit = bst.Unit
    for i in path:
        if isa(i, Unit):
            surface_units.add(i)
        elif has_path(i):
            update_surface_units(i.ID, i.streams, i.units, 
                                 surface_units, old_unit_connections)
    f = digraph_from_units(surface_units)
    for u, ins, outs in old_unit_connections:
        u._ins[:] = ins
        u._outs[:] = outs    
    return f

def update_surface_units(ID, streams, units, surface_units, old_unit_connections): # pragma: no coverage
    outs = []
    ins = []
    feeds = []
    products = []
    StreamUnit = bst.units.DiagramOnlyStreamUnit
    SystemUnit = bst.units.DiagramOnlySystemUnit
    
    for s in streams:
        source = s._source
        sink = s._sink
        if source in units and sink not in units:
            if sink: outs.append(s)
            else: products.append(s)
            u_io = (source, tuple(source.ins), tuple(source.outs))
            old_unit_connections.add(u_io)
        elif sink in units and source not in units:
            if source: ins.append(s)
            else: feeds.append(s)
            u_io = (sink, tuple(sink.ins), tuple(sink.outs))
            old_unit_connections.add(u_io)
    
    if len(feeds) > 1:
        feed = Stream(None)
        feed._ID = ''
        surface_units.add(StreamUnit('\n'.join([i.ID for i in feeds]),
                                     None, feed))
        ins.append(feed)
    else: ins += feeds
    
    if len(products) > 1:
        product = Stream(None)
        product._ID = ''
        surface_units.add(StreamUnit('\n'.join([i.ID for i in products]),
                                     product, None))
        outs.append(product)
    else: outs += products
    
    subsystem_unit = SystemUnit(ID, ins, outs)
    surface_units.add(subsystem_unit)

def digraph_from_units(units, **graph_attrs): # pragma: no coverage
    streams = bst.utils.streams_from_units(units)
    return digraph_from_units_and_streams(units, streams, **graph_attrs)

def digraph_from_units_and_streams(units, streams, **graph_attrs): # pragma: no coverage
    connections = get_all_connections(streams)
    return digraph_from_units_and_connections(units, connections, **graph_attrs)

def digraph_from_units_and_connections(units, connections, **graph_attrs): # pragma: no coverage
    f = blank_digraph(**graph_attrs)
    update_digraph_from_units_and_connections(f, units, connections)
    return f

def update_digraph_from_units_and_connections(f: Digraph, units, connections): # pragma: no coverage
    # Set up unit nodes
    unit_names = {}  # Contains full description (ID and line) by unit
    for u in units:
        node = u.get_node()
        f.node(**node)
        unit_names[u] = node['name']
    add_connections(f, connections, unit_names)    

def get_stream_connection(stream): # pragma: no coverage
    source = stream._source
    source_index = source._outs.index(stream) if source else None
    sink = stream._sink
    sink_index = sink._ins.index(stream) if sink else None
    return Connection(source, source_index, stream, sink_index, sink)

def get_all_connections(streams): # pragma: no coverage
    return {get_stream_connection(s)
            for s in streams 
            if (s._source or s._sink)}

def add_connection(f: Digraph, connection, unit_names=None): # pragma: no coverage
    source, source_index, stream, sink_index, sink = connection
    if unit_names is None:
        has_source = bool(source)
        has_sink = bool(sink)
    else:
        has_source = source in unit_names
        has_sink = sink in unit_names
    
    if stream and stream.ID:
        # Make stream nodes / unit-stream edges / unit-unit edges
        if has_sink and not has_source:
            # Feed stream case
            f.node(stream.ID)
            edge_in = sink._graphics.edge_in
            try:
                f.attr('edge', arrowtail='none', arrowhead='none',
                   tailport='e', **edge_in[sink_index])
            except:
                print(stream)
                f.attr('edge', arrowtail='none', arrowhead='none',
                   tailport='e', **edge_in[-1])
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
            try:
                f.attr('edge', arrowtail='none', arrowhead='normal',
                       **edge_in[sink_index], **edge_out[source_index])
            except:
                print(stream)
                f.attr('edge', arrowtail='none', arrowhead='normal',
                       **edge_in[-1], **edge_out[-1])
            f.edge(unit_names[source], unit_names[sink], label=stream.ID)
        else:
            f.node(stream.ID)
    elif has_sink and has_source:
        # Missing process stream case
        edge_in = sink._graphics.edge_in
        edge_out = source._graphics.edge_out
        f.attr('edge', arrowtail='none', arrowhead='normal',
               **edge_in[sink_index], **edge_out[source_index])
        f.edge(unit_names[source], unit_names[sink], style='dashed')

def add_connections(f: Digraph, connections, unit_names=None): # pragma: no coverage
    # Set attributes for graph and streams
    f.attr('node', shape='rarrow', fillcolor='#79dae8',
           style='filled', orientation='0', width='0.6',
           height='0.6', color='black', peripheries='1')
    f.attr('graph', overlap='orthoyx',
           outputorder='edgesfirst', nodesep='0.15', maxiter='1000000')
    f.attr('edge', dir='foward')
    for connection in connections:
        add_connection(f, connection, unit_names)

def display_digraph(digraph, format): # pragma: no coverage
    if format == 'svg':
        x = display.SVG(digraph.pipe(format=format))
    else:
        x = display.Image(digraph.pipe(format='png'))
    display.display(x)

def save_digraph(digraph, file, format): # pragma: no coverage
    if '.' not in file:
        file += '.' + format
    img = digraph.pipe(format=format)
    f = open(file, 'wb')
    f.write(img)
    f.close()
    
def finalize_digraph(digraph, file, format): # pragma: no coverage
    if file:
        save_digraph(digraph, file, format)
    else:
        display_digraph(digraph, format)