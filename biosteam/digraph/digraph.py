# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import numpy as np
from biosteam.utils.piping import ignore_docking_warnings
from warnings import warn
import biosteam as bst
from graphviz import Digraph
from IPython import display
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
           'extend_surface_units',
           'display_digraph',
           'save_digraph',
           'get_all_connections')

stream_node = dict(
    fillcolor='#55a8b5',
    fontcolor='white', 
    style='filled', 
    orientation='0',
    width='0.6', 
    height='0.6', 
    color='#90918e', 
    margin='0',
    peripheries='1',
    fontname="Arial",
)

class PenWidth:
    __slots__ = ('name', 'percentiles')
    def __init__(self, name, streams):
        self.name = name
        self.percentiles = np.percentile(
            [s.get_property(name) for s in streams if not s.isempty()] or [0], 
            [33, 66, 100]
        )
    
    def __call__(self, stream):
        value = stream.get_property(self.name)
        for width, percentile in enumerate(self.percentiles, 1):
            if value <= percentile: return str(width * 0.6 + 0.4)
        raise Exception(f'{self.name} beyond maximum')
    
    def __repr__(self):
        return f"{type(self).__name__}(name={repr(self.name)}, scale={self.scale:.3g})"


preferences = bst.preferences

def sort_streams(streams):
    return sorted(streams, key=lambda x: -x.F_mass + len(x.ID))

def has_path(obj):
    return hasattr(obj, 'path')

def blank_digraph(format='svg', maxiter='10000000', 
                  Damping='0.995', K='0.5', **graph_attrs):
    # Create a digraph and set direction left to right
    f = Digraph(format=format)
    f.attr(rankdir='LR', maxiter=maxiter, Damping=Damping, K=K,
           penwidth='0', color='none', bgcolor=preferences.background_color,
           fontcolor=preferences.label_color, fontname="Arial",
           labeljust='l', labelloc='t', fontsize='24',
           **graph_attrs)
    return f

def get_section_inlets_and_outlets(units, streams):
    outs = []
    ins = []
    for s in streams:
        source = s._source
        sink = s._sink
        if source in units and sink not in units:
            outs.append(s)
        elif sink in units and source not in units:
            ins.append(s)
    return ins, outs

@ignore_docking_warnings
def minimal_digraph(ID, units, streams, **graph_attrs):
    ins, outs = get_section_inlets_and_outlets(units, streams)
    product = Stream(None)
    product._ID = ''
    feed = Stream(None)
    feed._ID = ''
    ins = sort_streams(ins)
    outs = sort_streams(outs)
    feed_box = bst.units.DiagramOnlyStreamUnit('\n'.join([i.ID for i in ins]) or '-',
                                               None, feed)
    product_box = bst.units.DiagramOnlyStreamUnit('\n'.join([i.ID for i in outs]) or '-',
                                                  product, None)
    system_box = bst.units.DiagramOnlySystemUnit(ID, feed, product)
    return digraph_from_units([feed_box, system_box, product_box],
                              **graph_attrs)

@ignore_docking_warnings
def surface_digraph(path, **graph_attrs):
    surface_units = []
    old_unit_connections = set()
    isa = isinstance
    Unit = bst.Unit
    done = set()
    for i in path:
        if i in done: continue
        done.add(i)
        if isa(i, Unit):
            surface_units.append(i)
        elif has_path(i):
            extend_surface_units(i.ID, i.streams, i.units, 
                                 surface_units, old_unit_connections)
    f = digraph_from_units(surface_units, **graph_attrs)
    for u, ins, outs in old_unit_connections:
        u._ins[:] = ins
        u._outs[:] = outs    
    return f

@ignore_docking_warnings
def extend_surface_units(ID, streams, units, surface_units, old_unit_connections):
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
        feeds = sort_streams(feeds)
        feed_box = StreamUnit('\n'.join([i.ID for i in feeds]) or '-', None, feed)
        ins.append(feed)
    else: 
        feed_box = None
        ins += feeds
    
    if len(products) > 1:
        product = Stream(None)
        product._ID = ''
        products = sort_streams(products)
        product_box = StreamUnit('\n'.join([i.ID for i in products]) or '-', product, None)
        outs.append(product)
    else: 
        product_box = None
        outs += products
    
    subsystem_unit = SystemUnit(ID, ins, outs)
    for i in (feed_box, subsystem_unit, product_box):
        if i: surface_units.append(i)

def digraph_from_units(units, **graph_attrs):
    streams = bst.utils.streams_from_units(units)
    return digraph_from_units_and_streams(units, streams, **graph_attrs)

def digraph_from_units_and_streams(units, streams, **graph_attrs):
    connections = get_all_connections(streams)
    return digraph_from_units_and_connections(units, connections, **graph_attrs)

def digraph_from_system(system, **graph_attrs):
    f = blank_digraph(**graph_attrs) 
    other_streams = set()
    excluded_connections = set()
    unit_names = get_unit_names(f, system.unit_path)
    update_digraph_from_path(f, tuple(system.path) + system.facilities, 
                             system.recycle, 0, unit_names, excluded_connections,
                             other_streams)
    connections = get_all_connections(other_streams).difference(excluded_connections)
    add_connections(f, connections, unit_names)
    return f


def update_digraph_from_path(f, path, recycle, depth, unit_names,
                             excluded_connections,
                             other_streams):
    all_streams = set()
    units = set()
    subsystems = []
    isa = isinstance
    System = bst.System
    Unit = bst.Unit
    for i in path:
        if isa(i, Unit):
            units.add(i)
            all_streams.update(i._ins + i._outs)
        elif isa(i, System): 
            subsystems.append(i)
    if isa(recycle, bst.Stream): 
        recycles = [recycle] 
    elif recycle:
        recycles = recycle
    else:
        recycles = []
    streams = [i for i in all_streams if (not i.sink or i.sink in units) and (not i.source or i.source in units)]
    other_streams.update(all_streams.difference(streams))
    connections = get_all_connections(recycles)
    add_connections(f, connections, unit_names, color='#f98f60', fontcolor='#f98f60')
    excluded_connections.update(connections)
    connections = get_all_connections(streams).difference(excluded_connections)
    add_connections(f, connections, unit_names)
    excluded_connections.update(connections)
    depth += 1
    N_colors = len(preferences.depth_colors)
    color = preferences.depth_colors[(depth - 1) % N_colors]
    if preferences.fill_cluster:
        kwargs = dict(bgcolor=color, penwidth='0.2', color=preferences.stream_color)
    else:
        kwargs = dict(color=color, bgcolor='none', penwidth='5', style='dashed')
    for i in subsystems:
        with f.subgraph(name='cluster_' + i.ID) as c:
            c.attr(label=i.ID + f' [DEPTH {depth}]', fontname="Arial", 
                   labeljust='l', fontcolor=preferences.label_color, **kwargs)
            update_digraph_from_path(c, i.path, i.recycle, depth, unit_names, excluded_connections, other_streams)

def digraph_from_units_and_connections(units, connections, **graph_attrs):
    f = blank_digraph(**graph_attrs)
    update_digraph_from_units_and_connections(f, units, connections)
    return f

def get_unit_names(f: Digraph, units):
    unit_names = {}  # Contains full description (ID and line) by unit
    number = preferences.number_path
    profile = preferences.profile
    TicToc = bst.utils.TicToc
    info_by_unit = {}
    N_junctions = 0
    isa = isinstance
    for i, u in enumerate(units):
        if isa(u, bst.Junction):
            N_junctions += 1
            info_by_unit[u] = [[], None]
            continue
        if u in info_by_unit:
            old_data = info_by_unit[u]
            if number: old_data[0].append(str(i - N_junctions))
        else:
            if profile: # pragma: no cover
                t = TicToc()
                for n in range(10):
                    t.tic(); u.simulate(); t.toc()
                    if n > 1 and sum(t.record) > 0.2: break 
                time = f"{1000 * t.mean:.2g} ms"
            else:
                time = None
            index = [str(i - N_junctions)] if number else []
            info_by_unit[u] = [index, time] 
    for u, (index, time) in info_by_unit.items():
        node = u.get_node()
        name = node['name']
        info = ', '.join(index) 
        if time is not None:
            if info:
                info = f"{info}; {time}"
            else:
                info = time
        if info: name = f"[{info}] {name}"
        unit_names[u] = node['name'] = name
        f.node(**node)
    return unit_names

def update_digraph_from_units_and_connections(f: Digraph, units, connections):
    add_connections(f, connections, get_unit_names(f, units))    

def get_all_connections(streams):
    return {s.get_connection()
            for s in streams 
            if (s._source or s._sink)}

def add_connection(f: Digraph, connection, unit_names, pen_width=None, **edge_options):
    source, source_index, stream, sink_index, sink = connection
    has_source = source in unit_names
    has_sink = sink in unit_names
    style = 'dashed' if (stream.isempty() and not isinstance(stream.source, bst.units.DiagramOnlyStreamUnit)) else 'solid'
    f.attr('edge', label='', taillabel='', headlabel='', labeldistance='2',
           **edge_options)
    if stream:
        lines = []
        line = ''
        for word in stream.ID.split('_'):
            line += ' ' + word
            if len(line) > 10: 
                lines.append(line)
                line = ''
        if line: lines.append(line)
        ID = '\n'.join(lines)
        penwidth = pen_width(stream) if pen_width else '1.0'
        # Make stream nodes / unit-stream edges / unit-unit edges
        if has_sink and not has_source:
            # Feed stream case
            f.node(ID,
                   width='0.15', 
                   height='0.15',
                   shape='diamond',
                   fillcolor='#f98f60',
                   color=preferences.stream_color,
                   label='')
            inlet_options = sink._graphics.get_inlet_options(sink, sink_index)
            f.attr('edge', arrowtail='none', arrowhead='none', label=ID,
                   tailport='e', style=style, penwidth=penwidth, **inlet_options)
            f.edge(ID, unit_names[sink])
        elif has_source and not has_sink:
            # Product stream case
            f.node(ID, 
                   width='0.15', 
                   height='0.2',
                   shape='triangle',
                   orientation='270',
                   fillcolor='#f98f60',
                   color=preferences.stream_color,
                   label='')
            outlet_options = source._graphics.get_outlet_options(source, source_index)
            f.attr('edge', arrowtail='none', arrowhead='none', label=ID,
                   headport='w', style=style, penwidth=penwidth, **outlet_options)
            f.edge(unit_names[source], ID)
        elif has_sink and has_source:
            # Process stream case
            inlet_options = sink._graphics.get_inlet_options(sink, sink_index)
            outlet_options = source._graphics.get_outlet_options(source, source_index)
            f.attr('edge', arrowtail='none', arrowhead='normal', style=style, 
                   **inlet_options, penwidth=penwidth, **outlet_options)
            label = ID if preferences.label_streams else ''
            f.edge(unit_names[source], unit_names[sink], label=label)
        else:
            f.node(ID)
    elif has_sink and has_source:
        # Missing process stream case
        inlet_options = sink._graphics.get_inlet_options(sink, sink_index)
        outlet_options = source._graphics.get_outlet_options(source, source_index)
        f.attr('edge', arrowtail='none', arrowhead='normal',
               **inlet_options, **outlet_options)
        f.edge(unit_names[source], unit_names[sink], style='dashed')

def add_connections(f: Digraph, connections, unit_names, color=None, fontcolor=None, **edge_options):
    stream_width = preferences.stream_width
    if stream_width:
        pen_width = PenWidth(stream_width, [i.stream for i in connections])
    else:
        pen_width = None
    
    # Set attributes for graph and streams
    f.attr('graph', overlap='orthoyx', fontname="Arial",
           outputorder='edgesfirst', nodesep='0.5', ranksep='0.15', maxiter='1000000')
    f.attr('edge', dir='foward', fontname='Arial')
    f.attr('node', **stream_node)
    index = {j: i for i, j in unit_names.items()}
    length = len(index)
    def key(x):
        value = index.get(x.source, 0) + index.get(x.sink, length)
        if x.source_index:
            value += 1e-3 * x.source_index / len(x.source.outs)
        if x.sink_index:
            value += 1e-6 * x.sink_index / len(x.sink.ins)
        return value
    connections = sorted(connections, key=key)
    for connection in connections:
        add_connection(f, connection, unit_names, 
                       color=color or preferences.stream_color,
                       fontcolor=fontcolor or preferences.label_color,
                       pen_width=pen_width,
                       **edge_options)

def display_digraph(digraph, format): # pragma: no coverage
    if format == 'svg':
        x = display.SVG(digraph.pipe(format=format))
    else:
        x = display.Image(digraph.pipe(format='png'))
    display.display(x)

def save_digraph(digraph, file, format): # pragma: no coverage
    if '.' not in file:
        if format is None: format='svg'
        file += '.' + format
    img = digraph.pipe(format=format)
    f = open(file, 'wb')
    f.write(img)
    f.close()
    
def finalize_digraph(digraph, file, format): # pragma: no coverage
    if preferences.raise_exception: 
        if file: save_digraph(digraph, file, format)
        else: display_digraph(digraph, format)
    else:
        try:
            if file: save_digraph(digraph, file, format)
            else: display_digraph(digraph, format)
        except (OSError, TypeError) as exp:
            raise exp from None
        except Exception as exp: 
            warn(
                f"a '{type(exp).__name__}' was raised when generating "
                "graphviz diagram, possibly due to graphviz installation issues, "
                "make sure Graphviz executables are on your systems' PATH",
                RuntimeWarning
            )