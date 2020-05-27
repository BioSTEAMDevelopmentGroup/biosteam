# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 09:01:47 2019

@author: yoelr
"""
from ._unit import Unit
from ._facility import Facility
from ._digraph import digraph_from_units_and_streams, finalize_digraph

# %% Path checking

def feed_forward_recycle(path, other_path, recycle):
    try: path = path[path.index(recycle.sink):]
    except: return False
    other_path = other_path[other_path.index(recycle.sink):]
    return path != other_path


# %% Path tools

def find_linear_and_cyclic_paths_with_recycle(feed, ends):
    paths_with_recycle, linear_paths = find_paths_with_and_without_recycle(
        feed, ends)
    cyclic_paths_with_recycle = set()
    for path_with_recycle in paths_with_recycle:
        cyclic_path_with_recycle = path_with_recycle_to_cyclic_path_with_recycle(path_with_recycle)
        cyclic_paths_with_recycle.add(cyclic_path_with_recycle)
    return simplify_linear_paths(linear_paths), cyclic_paths_with_recycle

def find_paths_with_and_without_recycle(feed, ends):
    path = []
    paths_without_recycle  = set()
    paths_with_recycle = set()
    fill_path(feed, path, paths_with_recycle, paths_without_recycle,
              ends)
    return paths_with_recycle, paths_without_recycle

def fill_path(feed, path, paths_with_recycle,
              paths_without_recycle,
              ends):
    unit = feed.sink
    has_recycle = None
    if feed in ends:
        has_recycle = False
        for other_path, recycle in paths_with_recycle:
            if recycle is feed:
                has_recycle = feed_forward_recycle(path, other_path, recycle)
                if has_recycle: break
    if not unit or isinstance(unit, Facility) or has_recycle is False:
        paths_without_recycle.add(tuple(path))
    elif unit in path or has_recycle: 
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
            subnetworks.add(i)
        else:
            raise ValueError("path elements must be either Unit or Network "
                            f"objects not '{type(i).__name__}' objects")


# %% Network

class Network:
    """
    Create a Network object that defines a network of unit operations.
    
    Parameters
    ----------
    path : Iterable[Unit or Network]
        A path of unit operations and subnetworks.
    recycle : :class:`~thermosteam.Stream`
        A recycle stream, if any.
    
    """
    
    __slots__ = ('path', 'recycle', 'units', 'subnetworks',
                 'feeds', 'products', 'streams')
    
    def __init__(self, path, recycle=None):
        self.path = list(path)
        self.recycle = recycle
        self.units = units = set()
        self.subnetworks = subnetworks = set()
        self.streams = streams = set()
        self.feeds = feeds = set()
        self.products = products = set()
        load_network_components(path, units, streams, feeds,
                                products, subnetworks)
        
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
        facilities : Iterable[Facility]
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

    def copy_like(self, other):
        self.path = other.path
        self.recycle = other.recycle
        self.units = other.units
        self.subnetworks = other.subnetworks
        self.streams = other.streams
        self.feeds = other.feeds
        self.products = other.products

    def __contains__(self, other):
        if isinstance(other, Unit):
            return other in self.units
        elif isinstance(other, Network):
            return other in self.subnetworks
        else:
            return False
    
    def issubset(self, network):
        return self.units.issubset(network.units)
    
    def isdisjoint(self, network):
        return self.units.isdisjoint(network.units)
    
    def join_network(self, network, downstream=True):
        if self.isdisjoint(network):
            if downstream:
                self._append_network(network)
            else:
                self._appendleft_network(network)
        else:
            self._add_subnetwork(network)
    
    def join_network_at_unit(self, network, unit):
        isa = isinstance
        self._remove_overlap(network)
        has_overlap = False
        for index, item in enumerate(self.path):
            if isa(item, Network) and unit in item.units:
                item.join_network_at_unit(network, unit)
                self._update_from_newly_added_network(network)
                break
            elif unit == item:
                self._insert_network(index, network, has_overlap)
                break
    
    def _append_unit(self, unit):
        self.units.add(unit)
        self.products.update([i for i in unit._outs if i and not i._sink])
        self.feeds.update([i for i in unit._ins if i and not i._source])
        self.streams.update(unit._ins + unit._outs)
    
    def _update_from_newly_added_network(self, network):
        self.subnetworks.add(network)
        self.units.update(network.units)
        self.streams.update(network.streams)
        self.feeds.update(network.feeds)
        self.products.update(network.products)
    
    def _appendleft_network(self, network):
        if network.recycle:
            self.path.insert(0, network)
        else:
            for i in reversed(network.path): self.path.insert(0, i)
        self._update_from_newly_added_network(network)
    
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
        index_found = done = False
        subnetworks = self.subnetworks
        has_overlap = True
        path_tuple = tuple(path)
        recycle = self.recycle
        if recycle is subnetwork.recycle:
            # Feed forward scenario
            subpath = subnetwork.path
            for i, item in enumerate(subpath):
                if item not in path_tuple:
                    subpath = subpath[i:]
                    break
            if not subpath: return # Its the exact same network
            for i, item in enumerate(reversed(subpath)):
                if item not in path_tuple:
                    subpath = subpath[:-i]
                    break
            unit = path_tuple[-i]
            self.join_network_at_unit(Network(subpath), unit)
            return
        subunits = subnetwork.units
        for i in path_tuple:
            if isa(i, Unit):
                continue
            elif i.issubset(subnetwork):
                subnetwork._add_subnetwork(i)
                index = path.index(i)
                path.remove(i)
                subnetworks.remove(i)
                index_found = True
            elif subnetwork.issubset(i):
                i._add_subnetwork(subnetwork)
                done = True
        if index_found:
            self._insert_network(index, subnetwork)
        elif not done:
            for index, item in enumerate(path_tuple):
                if isa(item, Unit):
                    if item not in subunits: continue
                    self._insert_network(index, subnetwork)
                    has_overlap = False
                    done = True
                    break
                elif isa(item, Network):
                    if item.isdisjoint(subnetwork):
                        continue
                    else:
                        item._add_subnetwork(subnetwork)
                        done = True
                        break
        if has_overlap:
            self._remove_overlap(subnetwork)
        if not done:
            self._append_network(subnetwork)
        if len(path) == 1 and isa(path[0], Network):
            self.copy_like(path[0])

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
        if self.recycle:
            info += end + f"recycle={self.recycle})"
        else:
            info += ')'
        return info
    
    def _ipython_display_(self):
        self.show()
    
    def show(self):
        print(self._info(spaces=""))



    
        