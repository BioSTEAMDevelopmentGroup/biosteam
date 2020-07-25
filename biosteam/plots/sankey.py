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

__all__ = ('Handle', 
           'MassStreamHandle', 'MolarStreamHandle', 'VolumetricStreamHandle', 
           'CapitalNodeHandle', 'CarbonColorHandle', 'CarbonHandle', 
           'StreamNode', 'FeedNode', 
           'ProductNode', 'ProcessNode', 
           'sankey_figure', 'sankey_data',
           'node_dict', 'link_dict')

cycle_colors = ('aliceblue', 'chocolate', 'peachpuff', 
                'palevioletred', 'orange', 'springgreen', 'slategrey', 'goldenrod', 
                'mintcream', 
                'azure', 'darkseagreen', 'chartreuse', 'darkgray', 'khaki', 
                'plum', 'floralwhite', 'darkturquoise', 'steelblue', 'lightcyan', 
                'mistyrose', 'lightgray', 'lime', 'forestgreen',
                'darkblue', 'lightslategrey', 'burlywood', 'lightskyblue',
                'lightgrey', 'palegreen', 'indianred', 'royalblue', 
                'darkslateblue', 'blanchedalmond', 'seashell', 'dimgrey',
                'orchid', 'lightsalmon', 'beige', 'seagreen', 'teal',
                'lightblue', 'yellow', 'cornflowerblue', 'rosybrown', 'orangered',
                'peru', 'mediumorchid', 'mediumslateblue', 'darkgrey', 
                'slateblue', 'purple', 'lightslategray', 'limegreen', 'olivedrab', 
                'lightpink', 'olive', 'lightgoldenrodyellow', 'dimgray', 'mediumblue',
                'greenyellow', 'cyan', 'skyblue', 'green', 'ghostwhite',
                'hotpink', 'mediumvioletred', 
                'lawngreen', 'turquoise', 'powderblue', 'navy', 'gray', 
                'deeppink', 'lightcoral', 'mediumseagreen', 'maroon', 'honeydew', 
                'lavender', 'mediumturquoise', 'darkcyan', 'darksalmon',
                'paleturquoise', 'tomato', 'darkkhaki', 'darkgreen',
                'firebrick', 'ivory', 'mediumspringgreen', 'oldlace',
                'papayawhip', 'lemonchiffon', 'lightyellow', 'aquamarine', 
                'mediumpurple', 'aqua', 'blue', 'salmon', 'blueviolet',
                'lightseagreen', 'whitesmoke', 'linen', 'mediumaquamarine',
                'rebeccapurple', 'deepskyblue', 'sienna', 'violet', 'black',
                'darkgoldenrod', 'darkorchid', 'yellowgreen', 'darkviolet', 
                'pink', 'slategray', 'magenta', 'gainsboro', 'wheat', 'dodgerblue', 
                'fuchsia', 'cornsilk', 'palegoldenrod', 'saddlebrown', 
                'darkslategrey', 'indigo', 'snow', 'darkslategray', 
                'red', 'gold', 'coral', 'bisque', 'midnightblue', 
                'navajowhite', 'tan', 'moccasin', 'silver', 'brown', 
                'darkorange', 'darkred', 'antiquewhite', 'grey', 'crimson', 
                'white', 'sandybrown', 'darkmagenta', 'lavenderblush', 
                'lightgreen', 'thistle', 'darkolivegreen', 'lightsteelblue',
                'cadetblue')

# %% Main class used to handle the creating of Sankey plots

class Handle:
    """
    Create a Handle object for creating Sankey diagrams. You can subclass and
    implement the following to customize streams and processes.
    
    **Abstract methods**
    
    stream_width(stream) -> float
        Return link width.
    stream_color(stream) -> tuple[float, float, float]
        Return link color in RGB. Defaults to gray.
    process_color(unit, index) -> tuple[float, float, float]
        Return node color in RGB. By default, the color of proces nodes cycle 
        through colors defined in the variable, `biosteam.plots.sankey.color_cycle`.
    
    """
    __slots__ = ('size', 'nodes_index')
    
    def __init__(self):
        self.nodes_index = {} #: Dict[Object, ProcessNode] Object - node pairs.
        self.size = 0 #: [int] Number of process nodes created.
    
    def filter_streams(self, streams):
        has_width = self.stream_width
        return [i for i in streams if has_width(i)]
    
    stream_width = bst.utils.NotImplementedMethod
    stream_color = bst.utils.NotImplementedMethod
    process_color = bst.utils.NotImplementedMethod
    
    def next_index(self):
        index = self.size
        self.size = index + 1
        return index
    
    def feed_node(self, feed, name=None):
        feed_node = FeedNode(self, feed, self.next_index(), name)
        self.nodes_index[feed] = feed_node
        return feed_node
    
    def product_node(self, product, name=None):
        product_node = ProductNode(self, product, self.next_index(), name)
        self.nodes_index[product] = product_node
        return product_node
    
    def process_node(self, name, units):
        process_node = ProcessNode(self, name, self.next_index(), units)
        for i in units: self.nodes_index[i] = process_node 
        return process_node
    
    def nodes(self, *unit_groups, **units):
        nodes = []
        streams = set()
        streams_from_units = bst.utils.streams_from_units
        for group in unit_groups:
            node = self.process_node(group.name, group.units)
            nodes.append(node)
            streams.update(streams_from_units(group.units))
        for name, units in units.items():
            node = self.process_node(name.replace('_', ' '), units)
            nodes.append(node)
            streams.update(streams_from_units(units))
        streams = self.filter_streams(streams)
        nodes_index = self.nodes_index
        feeds = [i for i in streams if i._source not in nodes_index and i._sink]
        products = [i for i in streams if i._sink not in nodes_index and i._source]
        for feed in feeds:
            node = self.feed_node(feed)
            nodes.append(node)
        for product in products:
            node = self.product_node(product)
            nodes.append(node)
        return nodes

    def sankey_data(self, nodes, node_kwargs=None, link_kwargs=None):
        return sankey_data(self, nodes, node_kwargs, link_kwargs)

    def sankey_figure(self, nodes, node_kwargs=None, link_kwargs=None, sankey_data_kwargs=None, **kwargs):
        return sankey_figure(self, nodes, node_kwargs, link_kwargs, sankey_data_kwargs, **kwargs)

# %% Custum hanldes

class CapitalNodeHandle(Handle):
    """
    Create a CapitalHandle object that represents node colors by installed 
    equipment cost of the process.
    
    """
    __slots__ = ('max_installed_cost', 'process_color_source')
    def __init__(self, max_installed_cost=None, process_color_source=None):
        super().__init__()
        self.max_installed_cost = max_installed_cost
        self.process_color_source = process_color_source or bst.utils.colors.CABBI_orange

    def installed_cost_color(self, installed_cost):
        scale = 75 * installed_cost / self.max_installed_cost
        return self.process_color_source.shade(scale)

    def process_color(self, units, index):
        if self.max_installed_cost:
            installed_cost = sum([i.installed_cost for i in units])
            RGB = self.installed_cost_color(installed_cost).RGB
            return "rgba(%d, %d, %d, 1.0)" %tuple(RGB)
        else:
            return cycle_colors[index % len(cycle_colors)]
        
    def process_colorbar(self, N_levels=25, orientation='vertical'):
        colors = [self.installed_cost_color(0.).RGBn,
                  self.installed_cost_color(self.max_installed_cost).RGBn]
        return bst.plots.color_bar(colors, label='Installed equipment cost [million USD]',
                                   vmin=0, vmax=self.max_installed_cost / 1e6,
                                   N_levels=25, orientation=orientation)

class MassStreamHandle(Handle):
    """
    Create a MassHandle object that represents stream widths by mass flow.
    
    """
    __slots__ = ()
    def stream_width(self, stream): return stream.F_mass


class MolarStreamHandle(Handle):
    """
    Create a MolarHandle object that represents stream widths by molar flow.
    
    """
    __slots__ = ()
    def stream_width(self, stream): return stream.F_mol


class VolumetricStreamHandle(Handle):
    """
    Create a VolumetricHandle object that represents stream widths by volumetric 
    flow.
    
    """
    __slots__ = ()
    def stream_width(self, stream): return stream.F_vol


class CarbonStreamHandle(Handle):
    """
    Create a CarbonStreamHandle object that represents stream widths by
    carbon flow by weight.
    
    """
    __slots__ = ()
    def stream_width(self, stream): return stream.get_atomic_flow('C')


class CarbonColorHandle(Handle):
    """
    Create a CarbonColorHandle object that represents stream color by
    carbon content by weight.
    
    """
    __slots__ = ()
    
    def stream_carbon_content(self, stream):
        return stream.get_atomic_flow('C') * 12.01 / stream.F_mass
       
    def carbon_content_color(self, carbon_content):
        scale = 85. * (1 - carbon_content)
        return bst.utils.colors.grey.tint(scale)
    
    def stream_color(self, stream):
        carbon_content = self.stream_carbon_content(stream)
        RGB = self.carbon_content_color(carbon_content).RGB
        return "rgba(%d, %d, %d, 1.0)" %tuple(RGB)

    def stream_colorbar(self, N_levels=25, orientation='vertical'):
        colors = [self.carbon_content_color(0.).RGBn,
                  self.carbon_content_color(1.).RGBn]
        return bst.plots.color_bar(colors, label='Carbon content [wt. %]',
                                   vmin=0, vmax=100, N_levels=25,
                                   orientation=orientation)

class CarbonHandle(CarbonColorHandle, MassStreamHandle, CapitalNodeHandle):
    __slots__ = ()


# %% Main classes used by Handle objects to create Sankey plots

class StreamNode:
    """
    Abstract class for stream nodes.
    
    """
    __slots__ = ('handle', 'name', 'stream', 'index')
    
    def __init__(self, handle, stream, index, name=None):
        self.handle = handle
        self.name = name or stream.ID.replace('_', ' ')
        self.index = index
        self.stream = stream

    def color(self):
        return self.handle.stream_color(self.stream)
    
    def __repr__(self):
        return f"<{type(self).__name__}: {self.name}>"


class FeedNode(StreamNode):
    """
    Create a FeedNode object that represents a feed in a Sankey diagram.
    
    """
    __slots__ = ()
    
    def links(self):
        stream = self.stream
        sanky_group = self.handle.nodes_index[stream.sink]
        return [Link(self.handle, self, sanky_group, [stream])]
    
    
class ProductNode(StreamNode):
    """
    Create a ProductNode object that represents a product in a Sankey diagram.
    
    """
    __slots__ = ()
    
    def links(self):
        stream = self.stream
        sanky_group = self.handle.nodes_index[stream.source]
        return [Link(self.handle, sanky_group, self, [stream])]
    

class ProcessNode:
    """
    Create a ProcessNode object that represents a process in a Sankey diagram.
    
    """
    __slots__ = ('handle', 'name', 'index', 'units')
    
    def __init__(self, handle, name, index, units):
        self.handle = handle
        self.name = name
        self.index = index
        self.units = units
        
    def links(self):
        handle = self.handle
        nodes_index = handle.nodes_index
        units = frozenset(self.units)
        streams = bst.utils.outlets(units)
        streams = handle.filter_streams(streams)
        process_streams = [i for i in streams if i.sink not in units and i.sink in nodes_index]
        all_process_nodes = {nodes_index[i.sink] for i in process_streams}
        streams_by_process_node = {i:[] for i in all_process_nodes}
        for i in process_streams:
            process_node = nodes_index[i.sink]
            streams_by_process_node[process_node].append(i)
        return [Link(handle, self, sink, streams) for sink, streams in streams_by_process_node.items()]
  
    def color(self):
        return self.handle.process_color(self.units, self.index)

    __repr__ = StreamNode.__repr__


class Link:
    """
    Create a Link object that represents a link in a Sankey diagram.
    
    """
    __slots__ = ('handle', 'source', 'sink', 'streams')
    
    def __init__(self, handle, source, sink, streams):
        self.handle = handle #: [Handle]
        self.source = source #: [ProcessNode or StreamNode]
        self.sink = sink #: [ProcessNode or StreamNode]
        self.streams = streams #: [Stream]

    @property
    def stream(self):
        streams = self.streams
        if len(streams) == 1: return streams[0]
        return bst.Stream.sum(streams)

    def value(self):
        return self.handle.stream_width(self.stream)
    
    def color(self):
        return self.handle.stream_color(self.stream)
        
    def __repr__(self):
        return f"<{type(self).__name__}: {self.source.name} - {self.sink.name}>"

# %% Functions that use plotly to create diagrams

def sankey_figure(handle, nodes, node_kwargs=None, link_kwargs=None, sankey_data_kwargs=None, **kwargs):
    import plotly.graph_objects as go
    sankey_data_kwargs = sankey_data_kwargs or {}
    return go.Figure(data=sankey_data(handle, nodes, node_kwargs, link_kwargs, **sankey_data_kwargs), **kwargs)

def sankey_data(handle, nodes, arrangement = 'snap', node_kwargs=None, link_kwargs=None):
    import plotly.graph_objects as go
    node_kwargs = node_kwargs or {}
    link_kwargs = link_kwargs or {}
    links = sum([i.links() for i in nodes], [])
    return go.Sankey(
        arrangement = arrangement,
        node = node_dict(handle, nodes, **node_kwargs),
        link = link_dict(handle, links, **link_kwargs)
    )

def node_dict(handle, nodes, **kwargs):
    nodes = sorted(nodes, key=lambda x: x.index)
    dct = {'label': [i.name for i in nodes],
           **kwargs}
    if handle.process_color: dct['color'] = [i.color() for i in nodes]
    return dct
                                             
def link_dict(handle, links, **kwargs):
    dct = {
        'source': [i.source.index for i in links],
        'target': [i.sink.index for i in links],
        'value': [i.value() for i in links],
        **kwargs
    }
    if handle.stream_color: dct['color'] = [i.color() for i in links]
    return dct
