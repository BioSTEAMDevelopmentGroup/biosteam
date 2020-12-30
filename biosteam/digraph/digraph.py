# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from biosteam.utils.piping import ignore_docking_warnings
from warnings import warn
import biosteam as bst
from graphviz import Digraph
from IPython import display
from collections import namedtuple
from thermosteam import Stream

__all__ = ('digraph_from_system',
           'digraph_from_units',
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

@ignore_docking_warnings
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

@ignore_docking_warnings
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

@ignore_docking_warnings
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

def digraph_from_system(system, **graph_attrs): # pragma: no coverage
    f = blank_digraph(**graph_attrs) 
    unit_names = get_unit_names(f, system._unit_path + [i for i in system.facilities if isinstance(i, bst.Unit)])
    update_digraph_from_path(f, tuple(system.path) + system.facilities, 
                             system.recycle, 0, unit_names, set())
    return f

def update_digraph_from_path(f, path, recycle, depth, unit_names, 
                             excluded_connections): # pragma: no coverage
    streams = set()
    subsystems = []
    isa = isinstance
    System = bst.System
    Unit = bst.Unit
    for i in path:
        if isa(i, Unit):
            streams.update(i._ins + i._outs)
        elif isa(i, System): 
            subsystems.append(i)
    if isa(recycle, bst.Stream): 
        recycles = [recycle] 
    elif recycle:
        recycles = recycle
    else:
        recycles = []
    connections = get_all_connections(recycles)
    excluded_connections.update(connections)
    add_connections(f, connections, unit_names, color='#d71622')
    connections = get_all_connections(streams).difference(excluded_connections)
    excluded_connections.update(connections)
    add_connections(f, connections, unit_names, color='black')    
    depth += 1
    for i in subsystems:
        with f.subgraph(name='cluster_' + i.ID) as c:
            c.attr(label=i.ID + f' [DEPTH {depth}]', fontname="times-bold", 
                   style='dashed', bgcolor='#79bf823f')
            update_digraph_from_path(c, i.path, i.recycle, depth, unit_names, excluded_connections)

def digraph_from_units_and_connections(units, connections, **graph_attrs): # pragma: no coverage
    f = blank_digraph(**graph_attrs)
    update_digraph_from_units_and_connections(f, units, connections)
    return f

def get_unit_names(f: Digraph, units):
    unit_names = {}  # Contains full description (ID and line) by unit
    for u in units:
        node = u.get_node()
        f.node(**node)
        unit_names[u] = node['name']
    return unit_names

def update_digraph_from_units_and_connections(f: Digraph, units, connections): # pragma: no coverage
    add_connections(f, connections, get_unit_names(f, units))    

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

def add_connection(f: Digraph, connection, unit_names, edge_options): # pragma: no coverage
    source, source_index, stream, sink_index, sink = connection
    has_source = source in unit_names
    has_sink = sink in unit_names
    if stream and stream.ID:
        # Make stream nodes / unit-stream edges / unit-unit edges
        if has_sink and not has_source:
            # Feed stream case
            f.node(stream.ID)
            inlet_options = sink._graphics.get_inlet_options(sink, sink_index)
            f.attr('edge', arrowtail='none', arrowhead='none',
                   tailport='e', **inlet_options, **edge_options)
            f.edge(stream.ID, unit_names[sink])
        elif has_source and not has_sink:
            # Product stream case
            f.node(stream.ID)
            outlet_options = source._graphics.get_outlet_options(source, source_index)
            f.attr('edge', arrowtail='none', arrowhead='none',
                   headport='w', **outlet_options, **edge_options)
            f.edge(unit_names[source], stream.ID)
        elif has_sink and has_source:
            # Process stream case
            inlet_options = sink._graphics.get_inlet_options(sink, sink_index)
            outlet_options = source._graphics.get_outlet_options(source, source_index)
            f.attr('edge', arrowtail='none', arrowhead='normal',
                   **inlet_options, **outlet_options, **edge_options)
            label = stream.ID if bst.LABEL_PROCESS_STREAMS_IN_DIAGRAMS else ''
            f.edge(unit_names[source], unit_names[sink], label=label)
        else:
            f.node(stream.ID)
    elif has_sink and has_source:
        # Missing process stream case
        inlet_options = sink._graphics.get_inlet_options(sink, sink_index)
        outlet_options = source._graphics.get_outlet_options(source, source_index)
        f.attr('edge', arrowtail='none', arrowhead='normal',
               **inlet_options, **outlet_options, **edge_options)
        f.edge(unit_names[source], unit_names[sink], style='dashed')

def add_connections(f: Digraph, connections, unit_names, **edge_options): # pragma: no coverage
    # Set attributes for graph and streams
    f.attr('node', shape='rarrow', fillcolor='#79dae8',
           style='filled', orientation='0', width='0.6',
           height='0.6', color='black', peripheries='1')
    f.attr('graph', overlap='orthoyx',
           outputorder='edgesfirst', nodesep='0.15', maxiter='1000000')
    f.attr('edge', dir='foward')
    for connection in connections:
        add_connection(f, connection, unit_names, edge_options)

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