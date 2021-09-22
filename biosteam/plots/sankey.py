# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import biosteam as bst
from colorpalette import Color

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

# %% Helpful functions

def stream_cost(stream):
    return abs(stream.cost)

def abbreviate_word(word):
    if len(word) < 5: return word 
    original_word = word
    word = word.lower()
    vocals = 'aeiouy'
    for index, letter in enumerate(word):
        if letter in vocals: break
    for index, letter in enumerate(word[index:], index + 1):
        if letter not in vocals: break
    for index, letter in enumerate(word[index:], index + 1):
        if letter in vocals: 
            index -= 1
            break
    if index < len(word):
        return original_word[:index] + '.'
    else:
        return original_word
            
def abbreviate_name(name):
    words = name.split(' ')
    return ' '.join([abbreviate_word(i) for i in words])

def reduced_feed(feeds, unit_group):
    feed = bst.Stream.sum(feeds)
    feed._ID = ''
    feed._sink = unit_group.units[0]
    return feed

def reduced_product(products, unit_group):
    product = bst.Stream.sum(products)
    product._ID = ''
    product._source = unit_group.units[0]
    return product

def get_unit_group(stream, unit_groups):
    for group in unit_groups:        
        if stream in bst.utils.streams_from_units(group.units): return group
        
def reduced_streams(streams, unit_groups):
    same_IDs = {}
    for s in streams:
        ID = s.ID
        if ID in same_IDs:
            originals = same_IDs[ID]
            append = True
            unit_group = get_unit_group(s, unit_groups)
            for original in originals:
                if get_unit_group(original, unit_groups) is unit_group:
                    original.mix_from([original, s])
                    append = False
                    break
            if append: originals.append(s)
        else:
            same_IDs[ID] = [s]
    return sum(same_IDs.values(), [])


# %% Main class used to handle the creating of Sankey plots

class Handle: # pragma: no coverage
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
    __slots__ = (
        'size', 'nodes_index',
        'main_feeds', 'max_feeds',
        'main_products', 'max_products',
        'ignore',
    )
    
    def __init__(self, main_feeds=None, main_products=None, 
                 max_feeds=3, max_products=3, ignore=None, 
                 **kwargs):
        self.nodes_index = {} #: Dict[Object, ProcessNode] Object - node pairs.
        self.size = 0 #: [int] Number of process nodes created.
        self.main_feeds = main_feeds
        self.main_products = main_products
        self.max_feeds = max_feeds
        self.max_products = max_products
        self.ignore = ignore or ()
        self._init(**kwargs)
    
    def _init(self): pass
    
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
    
    def process_node(self, unit_group):
        units = unit_group.units
        process_node = ProcessNode(self, unit_group.name, self.next_index(), units)
        for i in units: self.nodes_index[i] = process_node 
        return process_node
    
    def nodes(self, unit_groups):
        if isinstance(unit_groups, dict):
            unit_groups = [bst.UnitGroup(i, j) for i,j in unit_groups.items()]
        nodes = []
        streams = set()
        streams_from_units = bst.utils.streams_from_units
        for group in unit_groups:
            nodes.append(self.process_node(group))
            streams.update(streams_from_units(group.units))
        streams = self.filter_streams(streams)
        nodes_index = self.nodes_index
        feeds = [i for i in streams if i._source not in nodes_index and i._sink and i not in self.ignore and not i.isempty()]
        main_feeds = self.main_feeds or sorted(feeds, key=stream_cost, reverse=True)[:self.max_feeds]
        feeds = reduced_streams(main_feeds, unit_groups)
        for feed in feeds:
            node = self.feed_node(feed)
            nodes.append(node)
        products = [i for i in streams if i._sink not in nodes_index and i._source and i not in self.ignore and not i.isempty()]
        main_products = self.main_products or sorted(products, key=stream_cost, reverse=True)[:self.max_products]
        products = reduced_streams(main_products, unit_groups)
        for product in products:
            node = self.product_node(product)
            nodes.append(node)
        return nodes

    def sankey_data(self, nodes, node_kwargs=None, link_kwargs=None):
        return sankey_data(self, nodes, node_kwargs, link_kwargs)

    def sankey_figure(self, nodes, node_kwargs=None, link_kwargs=None, sankey_data_kwargs=None, **kwargs):
        return sankey_figure(self, nodes, node_kwargs, link_kwargs, sankey_data_kwargs, **kwargs)

# %% Custum hanldes

class CapitalNodeHandle(Handle): # pragma: no coverage
    """
    Create a CapitalNodeHandle object that represents node colors by installed 
    equipment cost of the process.
    
    """
    __slots__ = ('max_installed_cost', 'process_color_source', 'shade')
    
    def _init(self, max_installed_cost=None, process_color_source=None, shade=0.75):
        self.max_installed_cost = max_installed_cost
        self.process_color_source = process_color_source or bst.utils.colors.CABBI_orange
        self.shade = shade

    def installed_cost_color(self, installed_cost):
        scale = 100 * self.shade * installed_cost / self.max_installed_cost
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

class MassStreamHandle(Handle): # pragma: no coverage
    """
    Create a MassHandle object that represents stream widths by mass flow.
    
    """
    __slots__ = ()
    def stream_width(self, stream): return stream.F_mass


class MolarStreamHandle(Handle): # pragma: no coverage
    """
    Create a MolarHandle object that represents stream widths by molar flow.
    
    """
    __slots__ = ()
    def stream_width(self, stream): return stream.F_mol


class VolumetricStreamHandle(Handle): # pragma: no coverage
    """
    Create a VolumetricHandle object that represents stream widths by volumetric 
    flow.
    
    """
    __slots__ = ()
    def stream_width(self, stream): return stream.F_vol


class CarbonStreamHandle(Handle): # pragma: no coverage
    """
    Create a CarbonStreamHandle object that represents stream widths by
    carbon flow by weight.
    
    """
    __slots__ = ()
    def stream_width(self, stream): return stream.get_atomic_flow('C')


class CarbonColorHandle(Handle): # pragma: no coverage
    """
    Create a CarbonColorHandle object that represents stream color by
    carbon content by weight.
    
    """
    __slots__ = ()
    
    def stream_carbon_content(self, stream):
        return stream.get_atomic_flow('C') * 12.01 / stream.F_mass
       
    def carbon_content_color(self, carbon_content):
        scale = 95. * (1 - carbon_content)
        return bst.utils.colors.grey.shade(10).tint(scale)
    
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


class CarbonHandle(CarbonColorHandle, MassStreamHandle, CapitalNodeHandle): # pragma: no coverage
    __slots__ = ()


# %% Main classes used by Handle objects to create Sankey plots

class StreamNode: # pragma: no coverage
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


class FeedNode(StreamNode): # pragma: no coverage
    """
    Create a FeedNode object that represents a feed in a Sankey diagram.
    
    """
    __slots__ = ()
    
    def links(self):
        stream = self.stream
        sanky_group = self.handle.nodes_index[stream.sink]
        return [Link(self.handle, self, sanky_group, [stream])]
    
    
class ProductNode(StreamNode): # pragma: no coverage
    """
    Create a ProductNode object that represents a product in a Sankey diagram.
    
    """
    __slots__ = ()
    
    def links(self):
        stream = self.stream
        sanky_group = self.handle.nodes_index[stream.source]
        return [Link(self.handle, sanky_group, self, [stream])]
    

class ProcessNode: # pragma: no coverage
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


class Link: # pragma: no coverage
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

def sankey_figure(handle, nodes, node_kwargs=None, link_kwargs=None, sankey_data_kwargs=None, **kwargs): # pragma: no coverage
    import plotly.graph_objects as go
    sankey_data_kwargs = sankey_data_kwargs or {}
    return go.Figure(data=sankey_data(handle, nodes, node_kwargs, link_kwargs, **sankey_data_kwargs), **kwargs)

def sankey_data(handle, nodes, arrangement = 'snap', node_kwargs=None, link_kwargs=None): # pragma: no coverage
    import plotly.graph_objects as go
    node_kwargs = node_kwargs or {}
    link_kwargs = link_kwargs or {}
    links = sum([i.links() for i in nodes], [])
    return go.Sankey(
        arrangement = arrangement,
        node = node_dict(handle, nodes, **node_kwargs),
        link = link_dict(handle, links, **link_kwargs)
    )

def node_dict(handle, nodes, **kwargs): # pragma: no coverage
    nodes = sorted(nodes, key=lambda x: x.index)
    dct = {'label': [i.name for i in nodes],
           **kwargs}
    if handle.process_color: dct['color'] = [i.color() for i in nodes]
    return dct
                                             
def link_dict(handle, links, **kwargs): # pragma: no coverage
    dct = {
        'source': [i.source.index for i in links],
        'target': [i.sink.index for i in links],
        'value': [i.value() for i in links],
        **kwargs
    }
    if handle.stream_color: dct['color'] = [i.color() for i in links]
    return dct
