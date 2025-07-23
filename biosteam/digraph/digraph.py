# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
# Copyright (C) 2022, Yoel Cortes-Pena <yoelcortes@gmail.com> and Ben Portner <https://github.com/BenPortner>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import numpy as np
from biosteam.utils.piping import (
    ignore_docking_warnings, Connection, 
)
from warnings import warn
import biosteam as bst
from graphviz import Digraph
from IPython import display
from thermosteam import AbstractStream, AbstractUnit
from xml.etree import ElementTree
from typing import Optional
import urllib
import re

__all__ = ('digraph_from_system',
           'digraph_from_units',
           'digraph_from_units_and_connections',
           'update_digraph_from_units_and_connections',
           'blank_digraph',
           'minimal_digraph',
           'surface_digraph',
           'finalize_digraph',
           'extend_surface_units',
           'display_digraph',
           'save_digraph')

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
        values = [s.get_property(name) for s in streams if not s.isempty()] or [0] 
        try:
            self.percentiles = np.percentile(
                values, 
                [33, 66, 100]
            )
        except:
            self.percentiles = 3 * [max(values)]
    
    def __call__(self, stream):
        value = stream.get_property(self.name)
        for width, percentile in enumerate(self.percentiles, 1):
            if value <= percentile: return str(width * 0.6 + 0.4)
        return str(width * 0.6 + 0.4)
    
    def __repr__(self):
        return f"{type(self).__name__}(name={repr(self.name)}, scale={self.scale:.3g})"


preferences = bst.preferences

def sort_streams(streams):
    return sorted(streams, key=lambda x: -x.F_mass + len(x.ID))

def has_path(obj):
    return hasattr(obj, 'path')

def blank_digraph(format='svg', maxiter='10000000000000000000', 
                  Damping='0.995', K='0.5', **graph_attrs):
    # Create a digraph and set direction left to right
    f = Digraph(format=format)
    f.attr(rankdir='LR', maxiter=maxiter, Damping=Damping, K=K, 
           penwidth='0', color='none', bgcolor=preferences.background_color,
           fontcolor=preferences.label_color, fontname="Arial",
           labeljust='l', labelloc='t', fontsize='24', constraint='false',
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
def minimal_digraph(ID, units, streams, auxiliaries=None, **graph_attrs):
    ins, outs = get_section_inlets_and_outlets(units, streams)
    thermo_water = bst.HeatUtility.thermo_water
    product = bst.Stream('.', thermo=thermo_water, Water=sum([i.F_mass for i in outs]), units='kg/hr')
    feed = bst.Stream('.', thermo=thermo_water, Water=sum([i.F_mass for i in ins]), units='kg/hr')
    ins = sort_streams(ins)
    outs = sort_streams(outs)
    system_box = bst.units.DiagramOnlySystemUnit(ID, feed, product)
    boxes = [system_box]
    if preferences.show_all_streams:
        name = '\n'.join([i.ID for i in ins])
    else:
        name = '\n'.join([i.ID for i in ins if not i.isempty()])
    if name:
        feed_box = bst.units.DiagramOnlyStreamUnit(name, None, feed)
        boxes.append(feed_box)
    if preferences.show_all_streams:
        name = '\n'.join([i.ID for i in outs])
    else:
        name = '\n'.join([i.ID for i in outs if not i.isempty()])
    if name:
        product_box = bst.units.DiagramOnlyStreamUnit(name, product, None)
        boxes.append(product_box)
    return digraph_from_units(boxes, auxiliaries=False,
                              **graph_attrs)

@ignore_docking_warnings
def surface_digraph(path, auxiliaries=None, **graph_attrs):
    surface_units = []
    old_unit_connections = set()
    isa = isinstance
    Unit = AbstractUnit
    done = set()
    for i in path:
        if i in done: continue
        done.add(i)
        if isa(i, Unit):
            surface_units.append(i)
        elif has_path(i):
            extend_surface_units(i.ID, i.streams, i.units, 
                                 surface_units, old_unit_connections)
    f = digraph_from_units(surface_units, auxiliaries=False, **graph_attrs)
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
        thermo_water = bst.HeatUtility.thermo_water
        feed = bst.Stream('.', thermo=thermo_water, Water=sum([i.F_mass for i in feeds]), units='kg/hr')
        feeds = sort_streams(feeds)
        if preferences.show_all_streams:
            name = '\n'.join([i.ID for i in feeds])
        else:
            name = '\n'.join([i.ID for i in feeds if not i.isempty()])
        if name:
            feed_box = bst.units.DiagramOnlyStreamUnit(name.replace('_', ' '), None, feed)
            ins.append(feed)
        else:
            feed_box = None
    else: 
        feed_box = None
        ins += feeds
    
    if len(products) > 1:
        thermo_water = bst.HeatUtility.thermo_water
        product = bst.Stream('.', thermo=thermo_water, Water=sum([i.F_mass for i in products]), units='kg/hr')
        products = sort_streams(products)
        if preferences.show_all_streams:
            name = '\n'.join([i.ID for i in products])
        else:
            name = '\n'.join([i.ID for i in products if not i.isempty()])
        if name:
            product_box = bst.units.DiagramOnlyStreamUnit(name.replace('_', ' '), product, None)
            outs.append(product)
        else:
            product_box = None
    else: 
        product_box = None
        outs += products
    
    subsystem_unit = SystemUnit(ID, ins, outs)
    for i in (feed_box, subsystem_unit, product_box):
        if i: surface_units.append(i)

def digraph_from_units(units, streams=None, auxiliaries=None, **graph_attrs):
    if auxiliaries:
        all_units = []
        streams = []
        stream_set = set()
        for unit in tuple(units):
            all_units.append(unit)
            for s in unit.ins + unit.outs:
                if s in stream_set: continue
                streams.append(s)
                stream_set.add(s)
            for name, auxunit in unit._diagram_nested_auxiliary_units_with_names(depth=auxiliaries):
                auxunit.owner = unit # In case units not created correctly through Unit.auxiliary method
                if isinstance(auxunit, AbstractUnit): 
                    auxunit._ID = name
                    all_units.append(auxunit)
                    for s in auxunit.ins + auxunit.outs:
                        if s in stream_set: continue
                        streams.append(s)
                        stream_set.add(s)
        units = all_units
    else:
        if streams is None:
            streams = []
            stream_set = set()
            for u in units:
                for s in u._ins + u._outs:
                    if s in stream_set: continue
                    streams.append(s)
                    stream_set.add(s)
        else:
            stream_set = set(streams)
    if auxiliaries:
        for stream in tuple(streams):
            if hasattr(stream, 'port'):
                other = stream.port.get_stream()
                if other in stream_set:
                    streams.remove(other)
                    stream_set.remove(other)
    digraph = digraph_from_units_and_connections(units, get_all_connections(streams), auxiliaries, **graph_attrs)
    return digraph

def digraph_from_system(system, path=None, auxiliaries=None, **graph_attrs):
    f = blank_digraph(**graph_attrs) 
    other_streams = set()
    excluded_connections = set()
    if not path: path = system.path
    unit_names = get_unit_names(f, path, with_auxiliaries=False)
    update_digraph_from_path(f, path, 
                             system.recycle, 0, unit_names, excluded_connections,
                             other_streams)
    connections = get_all_connections(other_streams, excluded_connections)
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
    Unit = AbstractUnit
    for i in path:
        if isa(i, Unit):
            units.add(i)
            all_streams.update(i._ins + i._outs)
        elif isa(i, System): 
            subsystems.append(i)
    if isa(recycle, AbstractStream): 
        recycles = [recycle] 
    elif recycle:
        recycles = recycle
    else:
        recycles = []
    streams = [i for i in all_streams if (not i.sink or i.sink in units) and (not i.source or i.source in units)]
    other_streams.update(all_streams.difference(streams))
    connections = get_all_connections(recycles, excluded_connections)
    add_connections(f, connections, unit_names, color='#f1777f', fontcolor='#f1777f')
    connections = get_all_connections(streams, excluded_connections)
    add_connections(f, connections, unit_names)
    depth += 1
    N_colors = len(preferences.depth_colors)
    color = preferences.depth_colors[(depth - 1) % N_colors]
    if preferences.fill_cluster:
        kwargs = dict(bgcolor=color, penwidth='0.2', color=preferences.stream_color)
    else:
        kwargs = dict(color=color, bgcolor='none', penwidth='0.75', style='solid')
    for i in subsystems:
        with f.subgraph(name='cluster_' + i.ID) as c:
            c.attr(label=i.ID, fontname="Arial", 
                   labeljust='l', fontcolor=preferences.label_color, **kwargs)
            update_digraph_from_path(c, i.path, i.recycle, depth, unit_names, excluded_connections, other_streams)

def digraph_from_units_and_connections(units, connections, with_auxiliaries, **graph_attrs):
    f = blank_digraph(**graph_attrs)
    update_digraph_from_units_and_connections(f, units, connections, with_auxiliaries)
    return f

def fill_info_from_path(path, indices, info_by_unit):
    isa = isinstance
    number = preferences.number_path
    profile = preferences.profile
    Timer = bst.utils.Timer
    for u in path:
        if isa(u, bst.Junction):
            continue # Do no include junctions
        if isa(u, bst.System):
            fill_info_from_path(u.path, [*indices, 0], info_by_unit)
            indices[-1] += 1
            continue
        index = '.'.join([str(i) for i in indices])
        indices[-1] += 1
        if u in info_by_unit:
            old_data = info_by_unit[u]
            if number: old_data[0].append(index)
        else:
            if profile: # pragma: no cover
                t = Timer()
                for n in range(1):
                    t.start(); u.simulate(); t.elapsed_time()
                    if n > 1 and sum(t.record) > 0.2: break 
                time = f"{1000 * t.mean:.2g} ms"
            else:
                time = None
            info_by_unit[u] = [[index] if number else [], time] 

def get_unit_names(f: Digraph, path, with_auxiliaries):
    unit_names = {}  # Contains full description (ID and line) by unit
    info_by_unit = {}
    fill_info_from_path(path, [0], info_by_unit)
    for u, (index, time) in info_by_unit.items():
        node = u.get_node()
        if 'xlabel' in node:
            key = 'xlabel'
        else:
            key = 'label'
        label = node[key]
        info = ', '.join(index) 
        if time is not None:
            if info:
                info = f"{info}; {time}"
            else:
                info = time
        if info: label = f"[{info}] {label}"
        node[key] = label
        unit_names[u] = node['name']
        if (with_auxiliaries and all([i in info_by_unit for i in u.auxiliary_units])
            and u._assembled_from_auxiliary_units()): continue
        f.node(**node)
    return unit_names

def update_digraph_from_units_and_connections(f: Digraph, units, connections, with_auxiliaries):
    add_connections(f, connections, get_unit_names(f, units, with_auxiliaries))    

def get_all_connections(streams, added_connections=None):
    if added_connections is None: added_connections = set()
    connections = []
    originals = {}
    superinlet = lambda s: s.__class__.__name__ == 'SuperpositionInlet'
    superoutlet = lambda s: s.__class__.__name__ == 'SuperpositionOutlet'
    for s in streams:
        original = s
        while hasattr(original, 'port'):
            original = original.port.get_stream()
        if (s._source or s._sink): 
            connection = s.get_connection(junction=False)
            if original in originals:
                source0, source_index0, stream0, sink_index0, sink0 = old_connection = originals[original]
                source1, source_index1, stream1, sink_index1, sink1 = connection
                if superoutlet(stream1) and superinlet(stream0):
                    connections.remove(old_connection)
                    added_connections.remove(old_connection)
                    stream = stream1.copy(ID='.' + s.ID)
                    connection = Connection(source1, source_index1, stream, sink_index0, sink0)
                elif superoutlet(stream0) and superinlet(stream1):
                    connections.remove(old_connection)
                    added_connections.remove(old_connection)
                    stream = stream1.copy(ID='.' + s.ID)
                    connection = Connection(source0, source_index0, stream, sink_index1, sink1)
            originals[original] = connection
            if connection and connection not in added_connections:
                connections.append(connection)
                added_connections.add(connection)
    return connections

def add_connection(f: Digraph, connection, unit_names, pen_width=None, **edge_options):
    source, source_index, stream, sink_index, sink = connection
    has_source = source in unit_names
    has_sink = sink in unit_names
    style = 'dashed' if (
        stream.isempty()
        and not isinstance(stream.source, bst.units.DiagramOnlyStreamUnit)
        and not isinstance(stream.sink, bst.units.DiagramOnlyStreamUnit)
    ) else 'solid'
    f.attr('edge', label='', taillabel='', headlabel='', labeldistance='2',
           **edge_options)
    tooltip = stream._get_tooltip_string(bst.preferences.graphviz_format, bst.preferences.tooltips_full_results)
    if stream:
        lines = []
        line = ''
        ID = stream.ID
        if ID:
            for word in ID.split('_'):
                line += ' ' + word
                if len(line) > 10: 
                    lines.append(line)
                    line = ''
            if line: lines.append(line)
            ID = '\n'.join(lines)
        else:
            ID = ' '
        ref = str(hash(stream))
        penwidth = pen_width(stream) if pen_width else '1.0'
        # Make stream nodes / unit-stream edges / unit-unit edges
        if has_sink and not has_source:
            # Feed stream case
            if not preferences.show_all_streams and stream.isempty(): return
            f.node(ref,
                   width='0.15', 
                   height='0.15',
                   shape='diamond',
                   fillcolor='#f98f60',
                   color=preferences.stream_color,
                   label='')
            inlet_options = sink._graphics.get_inlet_options(sink, sink_index)
            f.attr('edge', arrowtail='none', arrowhead='none', label=ID,
                   tailport='e', style=style, penwidth=penwidth, **inlet_options)
            f.edge(ref, unit_names[sink], labeltooltip=tooltip)
        elif has_source and not has_sink:
            # Product stream case
            if not preferences.show_all_streams and stream.isempty(): return
            f.node(ref, 
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
            f.edge(unit_names[source], ref, labeltooltip=tooltip)
        elif has_sink and has_source:
            # Process stream case
            inlet_options = sink._graphics.get_inlet_options(sink, sink_index)
            outlet_options = source._graphics.get_outlet_options(source, source_index)
            f.attr('edge', arrowtail='none', arrowhead='normal', style=style, 
                   **inlet_options, penwidth=penwidth, **outlet_options)
            label = ID if preferences.label_streams else ' '
            f.edge(unit_names[source], unit_names[sink], label=label, labeltooltip=tooltip)
    elif has_sink and has_source:
        # Missing process stream case
        inlet_options = sink._graphics.get_inlet_options(sink, sink_index)
        outlet_options = source._graphics.get_outlet_options(source, source_index)
        f.attr('edge', arrowtail='none', arrowhead='normal', label=' ',
               **inlet_options, **outlet_options)
        f.edge(unit_names[source], unit_names[sink], style='dashed', labeltooltip=tooltip)

def add_connections(f: Digraph, connections, unit_names, color=None, fontcolor=None, **edge_options):
    stream_width = preferences.stream_width
    if stream_width:
        pen_width = PenWidth(stream_width, [i.stream for i in connections])
    else:
        pen_width = None
    
    # Set attributes for graph and streams
    f.attr('graph', overlap='orthoyx', fontname="Arial", outputorder='edgesfirst', 
           nodesep='0.5', ranksep='0.15', maxiter='1000000')
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

def fix_valve_symbol_in_svg_output(
        img:bytes,
        unit_color: Optional[str] = None,
        unit_periphery_color: Optional[str] = None,
        label_color: Optional[str] = None,
):
    """Fix valve symbols because images cannot be loaded when choosing `format='svg'`"""
    if unit_color is None: unit_color = bst.preferences.unit_color
    if unit_periphery_color is None: unit_periphery_color = bst.preferences.unit_periphery_color
    if label_color is None: label_color = bst.preferences.unit_label_color
    # get all image tags
    tree = ElementTree.fromstring(img)
    images = [e for e in tree.iter() if 'image' in e.tag]
    # make polygon
    parent_map = {c: p for p in tree.iter() for c in p}
    getchildren = lambda pm: pm.getchildren() if hasattr(pm, 'getchildren') else list(pm)
    polygons = [c for i in images for c in getchildren(parent_map[i]) if 'polygon' in c.tag]
    for p in polygons:
        points = p.attrib["points"].split(' ')
        # turn rect into valve symbol
        buffer = points[1]
        points[1] = points[2]
        points[2] = buffer
        p.attrib["points"] = ' '.join(points)
        # fix color
        p.attrib["fill"] = unit_color.split(':')[-1] # In case of gradiant color (white:#CDCDCD)
        p.attrib["stroke"] = unit_periphery_color
    # fix label text color and position
    label_image = [(c, i) for i in images for c in getchildren(parent_map[parent_map[parent_map[i]]]) if 'text' in c.tag]
    for l,i in label_image:
        l.attrib["fill"] = label_color
        width = float(i.attrib['width'].replace('px',''))
        x = float(i.attrib["x"])
        l.attrib["x"] = f"{x+width/2}"
    # delete image tags
    for i in images:
        parent_map[i].remove(i)
    return ElementTree.tostring(tree)

def inject_javascript(img:bytes):
    html = ElementTree.Element('html')
    head = ElementTree.SubElement(html, 'head')
    # insert css
    links = [
        "https://unpkg.com/tippy.js@6.3.7/themes/translucent.css",
        "https://rawcdn.githack.com/BioSTEAMDevelopmentGroup/biosteam/e065aca079c216d72b75949bbcbb74a3bbddb75d/biosteam/digraph/digraph.css",
    ]
    for href in links:
        link = ElementTree.SubElement(head, 'link')
        link.set("rel", "stylesheet")
        link.set("href", href)
    # insert javascript
    srcs = [
        "https://unpkg.com/@popperjs/core@2",
        "https://unpkg.com/tippy.js@6",
        "https://rawcdn.githack.com/BioSTEAMDevelopmentGroup/biosteam/e065aca079c216d72b75949bbcbb74a3bbddb75d/biosteam/digraph/digraph.js",
    ]
    for src in srcs:
        script = ElementTree.SubElement(head, 'script')
        script.set("src", src)
    # body
    body = ElementTree.SubElement(html, 'body')
    svg = ElementTree.fromstring(img)
    
    # getiterator is deprecated in Python 3.9
    getiter = lambda etree: (getattr(etree, 'getiterator', None) or getattr(etree, 'iter'))()
    
    # remove namespaces from tags and attributes
    for e in getiter(svg):
        e.tag = re.sub("{.*?}","",e.tag)
        for key, value in e.attrib.copy().items():
            if "{" in key:
                clean = re.sub("{.*?}","",key)
                e.attrib[clean] = value
                del e.attrib[key]
    # make tippy tooltips
    for e in getiter(svg):
        for key, value in e.attrib.copy().items():
            if key == "class" and value in ["node", "edge"]:
                title = e.find("./title")
                if title is not None:
                    default_tooltip = title.text
                else:
                    default_tooltip = e.text
                custom_tooltip = (e.find("./g/a") or e).attrib.get("title", None)
                tooltip = custom_tooltip or default_tooltip
                if tooltip is not None:
                    e.attrib["data-tippy-content"] = tooltip.strip()
    # remove default tooltips
    for e in getiter(svg):
        t = e.find("./title")
        if t is not None:
            e.remove(t)
    body.append(svg)
    # add docstring declaration
    s = ElementTree.tostring(html, encoding='utf8', method='html')
    s = b"<!DOCTYPE html>"+s
    return s

def display_digraph(digraph, format, height=None): # pragma: no coverage
    if format is None: format = preferences.graphviz_format
    if height is None: height = '400px'
    if format == 'svg':
        img = digraph.pipe(format='svg')
        img = fix_valve_symbol_in_svg_output(img)
        x = display.SVG(img)
        display.display(x)
    # TODO: Consult about this complicated Javascript injection
    elif format == 'html':
        img = digraph.pipe(format='svg')
        img = fix_valve_symbol_in_svg_output(img)
        img = inject_javascript(img)
        data_uri = 'data:text/html;charset=utf-8,' + urllib.parse.quote(img)
        x = display.IFrame(src=data_uri, width='100%', height=height,
                            extras=['allowtransparency="true"'])
        display.display(x)
    else:
        x = display.Image(digraph.pipe(format='png'))
        display.display(x)

def save_digraph(digraph, file, format): # pragma: no coverage
    if '.' not in file:
        if format is None: format = preferences.graphviz_format
        file += '.' + format
    elif format is None:
        format = file.split()[-1]
    else:
        raise ValueError(
            "cannot specify format extension; file already has format "
           f"extension '{file.split()[-1]}'"
        )
    if format == 'html':
        try:
            img = digraph.pipe(format='svg')
        except Exception as e:
            try:
                from signal import signal, SIGPIPE, SIG_DFL
                signal(SIGPIPE, SIG_DFL)
            except ImportError:
                raise e
        img = fix_valve_symbol_in_svg_output(img)
        img = inject_javascript(img)
    else:
        img = digraph.pipe(format=format)
        if format == 'svg': img = fix_valve_symbol_in_svg_output(img)
    f = open(file, 'wb')
    f.write(img)
    f.close()
    
def finalize_digraph(digraph, file, format, height=None): # pragma: no coverage
    if preferences.raise_exception: 
        if file: save_digraph(digraph, file, format)
        else: display_digraph(digraph, format, height)
    else:
        try:
            if file: save_digraph(digraph, file, format)
            else: display_digraph(digraph, format, height)
        except (OSError, TypeError) as exp:
            raise exp from None
        except Exception as exp: 
            warn(
                f"a '{type(exp).__name__}' was raised when generating "
                "graphviz diagram, possibly due to graphviz installation issues, "
                "make sure Graphviz executables are on your systems' PATH",
                RuntimeWarning
            )