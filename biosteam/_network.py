# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from ._unit import Unit
from ._facility import Facility
from .digraph import digraph_from_units_and_streams, finalize_digraph
from thermosteam import Stream

# %% Path tools

def find_linear_and_cyclic_paths_with_recycle(feed, ends):
    paths_with_recycle, linear_paths = find_paths_with_and_without_recycle(
        feed, ends)
    cyclic_paths_with_recycle = set()
    for path_with_recycle in paths_with_recycle:
        cyclic_path_with_recycle = path_with_recycle_to_cyclic_path_with_recycle(path_with_recycle)
        cyclic_paths_with_recycle.add(cyclic_path_with_recycle)
    cyclic_paths_with_recycle = sorted(cyclic_paths_with_recycle, key=lambda x: -len(x[0]))
    return simplify_linear_paths(linear_paths), cyclic_paths_with_recycle

def find_paths_with_and_without_recycle(feed, ends):
    path = []
    paths_without_recycle  = set()
    paths_with_recycle = set()
    fill_path(feed, path, paths_with_recycle, paths_without_recycle, ends)
    return paths_with_recycle, paths_without_recycle

def fill_path(feed, path, paths_with_recycle,
              paths_without_recycle,
              ends):
    unit = feed.sink
    has_recycle = None
    if feed in ends:
        has_recycle = False
        if unit in path:
            for other_path, recycle in paths_with_recycle:
                has_recycle = recycle.sink is unit
                if has_recycle: break
    if not unit or isinstance(unit, Facility) or has_recycle is False:
        paths_without_recycle.add(tuple(path))
    elif has_recycle or unit in path: 
        path_with_recycle = tuple(path), feed
        paths_with_recycle.add(path_with_recycle)
        ends.add(feed)
    else:
        path.append(unit)
        first_outlet, *other_outlets = unit.outs
        for outlet in other_outlets:
            new_path = path.copy()
            fill_path(outlet, new_path,
                      paths_with_recycle,
                      paths_without_recycle,
                      ends)
        fill_path(first_outlet, path,
                  paths_with_recycle,
                  paths_without_recycle,
                  ends)

def path_with_recycle_to_cyclic_path_with_recycle(path_with_recycle):
    path, recycle = path_with_recycle
    unit = recycle.sink
    recycle_index = path.index(unit)
    return (path[recycle_index:], recycle)

def simplify_linear_paths(linear_paths):
    simplified_linear_paths = []
    linear_paths = sorted(linear_paths, key=len)
    while linear_paths:
        smaller_path, *linear_paths = linear_paths
        simplified_path = simplify_linear_path(smaller_path, linear_paths)
        if simplified_path:
            simplified_linear_paths.append(simplified_path)
    return simplified_linear_paths
    
def simplify_linear_path(path, other_paths):
    simplified_path = list(path)
    for unit in path:
        for other_path in other_paths:
            if unit in other_path:
                simplified_path.remove(unit)
                break
    return simplified_path

def load_network_components(path, units, streams, feeds,
                            products, subnetworks):
    isa = isinstance
    for i in path:
        if isa(i, Unit): 
            units.add(i)
            streams.update(i._ins + i._outs)
            feeds.update([i for i in i._ins if not i._source])
            products.update([i for i in i._outs if not i._sink])
        elif isa(i, Network):
            feeds.update(i.feeds)
            streams.update(i.streams)
            products.update(i.products)
            units.update(i.units)
            subnetworks.append(i)
        else:
            raise ValueError("path elements must be either Unit or Network "
                            f"objects not '{type(i).__name__}' objects")


# %% Network

class Network:
    """
    Create a Network object that defines a network of unit operations.
    
    Parameters
    ----------
    path : Iterable[:class:`~biosteam.Unit` or :class:`~biosteam.Network`]
        A path of unit operations and subnetworks.
    recycle : :class:`~thermosteam.Stream`
        A recycle stream, if any.
    
    Examples
    --------
    Create a network representing two nested recycle loops:
        
    >>> from biosteam import (
    ...     main_flowsheet as f,
    ...     Pump, Mixer, Splitter,
    ...     Stream, settings, Network,
    ... )
    >>> f.set_flowsheet('two_nested_recycle_loops')
    >>> settings.set_thermo(['Water'], cache=True)
    >>> feedstock = Stream('feedstock', Water=1000)
    >>> water = Stream('water', Water=10)
    >>> recycle = Stream('recycle')
    >>> inner_recycle = Stream('inner_recycle')
    >>> product = Stream('product')
    >>> inner_water = Stream('inner_water', Water=10)
    >>> P1 = Pump('P1', feedstock)
    >>> P2 = Pump('P2', water)
    >>> P3 = Pump('P3', inner_water)
    >>> M1 = Mixer('M1', [P1-0, P2-0, recycle])
    >>> M2 = Mixer('M2', [M1-0, P3-0, inner_recycle])
    >>> S2 = Splitter('S2', M2-0, ['', inner_recycle], split=0.5)
    >>> S1 = Splitter('S1', S2-0, [product, recycle], split=0.5)
    >>> network = Network(
    ... [P1,
    ...  P2,
    ...  P3,
    ...  Network(
    ...     [M1,
    ...      Network(
    ...          [M2,
    ...           S2],
    ...          recycle=inner_recycle),
    ...      S1],
    ...     recycle=recycle)])
    >>> network.show()
    Network(
        [P1,
         P2,
         P3,
         Network(
            [M1,
             Network(
                [M2,
                 S2],
                recycle=S2-1),
             S1],
            recycle=S1-1)])
    
    """
    
    __slots__ = ('path', 'recycle', 'units', 'subnetworks',
                 'feeds', 'products', 'streams')
    
    def __init__(self, path, recycle=None):
        self.path = list(path)
        self.recycle = recycle
        self.units = units = set()
        self.subnetworks = subnetworks = []
        self.streams = streams = set()
        self.feeds = feeds = set()
        self.products = products = set()
        load_network_components(path, units, streams, feeds,
                                products, subnetworks)
     
    def __eq__(self, other):
        return isinstance(other, Network) and self.path == other.path and self.recycle is other.recycle
        
    @classmethod
    def from_feedstock(cls, feedstock, feeds=(), ends=None):
        """
        Create a Network object from a feedstock.
        
        Parameters
        ----------
        feedstock : :class:`~thermosteam.Stream`
            Main feedstock of the process.
        feeds : Iterable[:class:`~thermosteam.Stream`]
            Additional feeds to the process.
        facilities : Iterable[:class:`~biosteam.Facility`]
            Offsite facilities that are simulated only after 
            completing the path simulation.
        ends : Iterable[:class:`~thermosteam.Stream`]
            Streams that not products, but are ultimately specified through
            process requirements and not by its unit source.
            
        """
        ends = set(ends) or set()
        linear_paths, cyclic_paths_with_recycle = find_linear_and_cyclic_paths_with_recycle(
            feedstock, ends)
        network = Network(sum(reversed(linear_paths), []))
        recycle_networks = [Network(path, recycle) for path, recycle
                            in cyclic_paths_with_recycle]
        for recycle_network in recycle_networks:
            network.join_network(recycle_network)
        isa = isinstance
        for feed in feeds:
            streams = network.streams
            if feed in streams or isa(feed.sink, Facility):
                continue
            if ends:
                new_ends = streams.union(ends)
            else:
                new_ends = streams
            upstream_network = cls.from_feedstock(feed, (), new_ends)
            connections = upstream_network.streams.intersection(streams)
            connecting_units = {stream._sink for stream in connections
                                if stream._source and stream._sink}
            N_connections = len(connecting_units)
            if N_connections == 1:
                connecting_unit, = connecting_units
                network.join_network_at_unit(upstream_network,
                                             connecting_unit)
            elif N_connections == 0:
                network._append_network(upstream_network)
            else:
                network.join_network(upstream_network)
        return network

    def __contains__(self, other):
        if isinstance(other, Unit):
            return other in self.units
        elif isinstance(other, Network):
            return other in self.subnetworks
        else:
            return False
    
    def isdisjoint(self, network):
        return self.units.isdisjoint(network.units)
    
    def join_network(self, network):
        if self.isdisjoint(network):
            # Always join downstream
            self._append_network(network)
        else:
            self._add_subnetwork(network)
    
    def join_network_at_unit(self, network, unit):
        isa = isinstance
        for index, item in enumerate(self.path):
            if isa(item, Network) and unit in item.units:
                self._remove_overlap(network)
                if network.recycle:
                    item.join_network_at_unit(network, unit)
                    self._update_from_newly_added_network(network)
                else:
                    self._insert_network(index, network, False)
                return
            elif unit == item:
                self._remove_overlap(network)
                self._insert_network(index, network, False)
                return
        raise RuntimeError('unit not in path')
    
    def _update_from_newly_added_network(self, network):
        self.subnetworks.append(network)
        self.units.update(network.units)
        self.streams.update(network.streams)
        self.feeds.update(network.feeds)
        self.products.update(network.products)
    
    def _append_network(self, network):
        if network.recycle:
            self.path.append(network)
        else:
            self.path.extend(network.path)
        self._update_from_newly_added_network(network)
    
    def _insert_network(self, index, network, has_overlap=True):
        path = self.path
        if has_overlap:
            self._remove_overlap(network)
        if network.recycle:
            path.insert(index, network)
        else:
            for item in reversed(network.path):
                path.insert(index, item)
        self._update_from_newly_added_network(network)
    
    def _add_subnetwork(self, subnetwork):
        path = self.path
        isa = isinstance
        done = False
        subnetworks = self.subnetworks
        has_overlap = True
        path_tuple = tuple(path)
        recycle = self.recycle
        if recycle and subnetwork.recycle and recycle.sink is subnetwork.recycle.sink:
            subnetwork.recycle = None # Feed forward scenario
            self._add_subnetwork(subnetwork)
            return
        subunits = subnetwork.units
        for index, i in enumerate(path_tuple):
            if isa(i, Network) and not subnetwork.isdisjoint(i):
                i._add_subnetwork(subnetwork)
                self._update_from_newly_added_network(subnetwork)
                done = True
                break
        if not done:
            for index, item in enumerate(path_tuple):
                if isa(item, Unit) and item in subunits:
                    self._insert_network(index, subnetwork)
                    has_overlap = False
                    done = True
                    break
        if has_overlap:
            self._remove_overlap(subnetwork)
        if not done:
            self._append_network(subnetwork)
        if len(path) == 1:
            subnetwork = path[0]
            if isa(subnetwork, Network):
                self.path = subnetwork.path
                self.recycle = subnetwork.recycle

    def _remove_overlap(self, subnetwork):
        path = self.path
        for item in tuple(path):
            if item in subnetwork: 
                path.remove(item)

    def diagram(self, file=None, format='png'):
        units = self.units
        f = digraph_from_units_and_streams(units, self.streams)
        finalize_digraph(f, file, format)

    def __repr__(self):
        return f"{type(self).__name__}(path={self.path}, recycle={self.recycle})"
    
    def _info(self, spaces):
        info = f"{type(self).__name__}("
        spaces += 4 * " "
        end = ',\n' + spaces
        path_info = []
        path = self.path
        isa = isinstance
        info += '\n' + spaces
        for i in path:
            if isa(i, Unit):
                path_info.append(str(i))
            else:
                path_info.append(i._info(spaces))
        info += '[' + (end + " ").join(path_info) + ']'
        recycle = self.recycle
        if recycle:
            if isinstance(recycle, Stream):
                recycle = recycle._source_info()
            else:
                recycle = ", ".join([i._source_info() for i in recycle])
                recycle = '[' + recycle + ']'
            info += end + f"recycle={recycle})"
        else:
            info += ')'
        return info
    
    def _ipython_display_(self):
        self.show()
    
    def show(self):
        print(self._info(spaces=""))



    
        