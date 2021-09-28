# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from ._unit import Unit
from ._facility import Facility
from .utils import streams_from_units
from warnings import warn
from thermosteam import Stream
import biosteam as bst

# %% Other tools

def get_recycle_sink(recycle):
    if isinstance(recycle, Stream):
        return recycle._sink
    elif isinstance(recycle, set):
        for i in recycle: return i._sink
    else: # pragma: no cover
        raise ValueError('recycle must be either a stream or a set; not a '
                        f"'{type(recycle).__name__}' object")

# %% Path tools

class PathSource:
    
    __slots__ = ('source', 'units')
    
    def __init__(self, source, ends=None):
        self.source = source
        if isinstance(source, bst.Unit):
            self.units = units = source.get_downstream_units(ends=ends, facilities=False)
        elif isinstance(source, Network):
            self.units = units = set()
            for i in source.units: units.update(i.get_downstream_units(ends=ends, facilities=False))
        else:
            raise TypeError('source must be a unit or a network; '
                           f'not a {type(source).__name__} object')
        assert not self.downstream_from(self), (
            'recycle system encountered in network sorting algorithm; '
            'please report problem to the biosteam bugtracker'
        )
            
    def downstream_from(self, other):
        source = self.source
        if isinstance(source, bst.Unit):
            return source in other.units
        else:
            return self is not other and any([i in other.units for i in source.units])
        
    def __repr__(self):
        return f"{type(self).__name__}({str(self.source)})"
        

def find_linear_and_cyclic_paths_with_recycle(feed, ends):
    paths_with_recycle, linear_paths = find_paths_with_and_without_recycle(
        feed, ends)
    cyclic_paths_with_recycle = []
    for path_with_recycle in paths_with_recycle:
        cyclic_path_with_recycle = path_with_recycle_to_cyclic_path_with_recycle(path_with_recycle)
        cyclic_paths_with_recycle.append(cyclic_path_with_recycle)
    cyclic_paths_with_recycle.sort(key=lambda x: -len(x[0]))
    return simplified_linear_paths(linear_paths), cyclic_paths_with_recycle

def find_paths_with_and_without_recycle(feed, ends):
    paths_without_recycle  = []
    paths_with_recycle = []
    fill_path(feed, [], paths_with_recycle, paths_without_recycle, ends)
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
        paths_without_recycle.append(path)
    elif has_recycle or unit in path: 
        path_with_recycle = path, feed
        paths_with_recycle.append(path_with_recycle)
        ends.add(feed)
    else:
        path.append(unit)
        first_outlet, *other_outlets = unit._outs
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

def simplified_linear_paths(linear_paths):
    if not linear_paths: return linear_paths
    linear_paths.sort(key=len)
    units, *unit_sets = [set(i) for i in linear_paths]
    for i in unit_sets: units.update(i)
    simplified_paths = []
    for i, path in enumerate(linear_paths):
        simplify_linear_path(path, unit_sets[i:])
        if path:
            add_back_ends(path, units)
            simplified_paths.append(path)
    simplified_paths.reverse()
    return simplified_paths
    
def simplify_linear_path(path, unit_sets):
    if path and unit_sets:
        for unit in path.copy():
            for unit_set in unit_sets:
                if unit in unit_set:
                    path.remove(unit)
                    break

def add_back_ends(path, units):
    for outlet in path[-1]._outs:
        sink = outlet._sink 
        if sink in units: 
            path.append(sink)

def nested_network_units(path):
    units = set()
    isa = isinstance
    for i in path:
        if isa(i, Unit): 
            units.add(i)
        elif isa(i, Network):
            units.update(i.units)
        else: # pragma: no cover
            raise ValueError("path elements must be either Unit or Network "
                            f"objects not '{type(i).__name__}' objects")
    return units


# %% Network

class Network:
    """
    Create a Network object that defines a network of unit operations.
    
    Parameters
    ----------
    path : Iterable[:class:`~biosteam.Unit` or :class:`~biosteam.Network`]
        A path of unit operations and subnetworks.
    recycle : :class:`~thermosteam.Stream` or set[:class:`~thermosteam.Stream`]
        A recycle stream(s), if any.
    
    Examples
    --------
    Create a network representing two nested recycle loops:
        
    >>> from biosteam import (
    ...     main_flowsheet as f,
    ...     Pump, Mixer, Splitter,
    ...     Stream, settings
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
    
    __slots__ = ('path', 'units', 'recycle', 'recycle_sink')
    
    def __init__(self, path, recycle=None):
        self.path = path
        self.recycle = recycle
        self.recycle_sink = get_recycle_sink(recycle) if recycle else None
        try: self.units = set(path)
        except: self.units = nested_network_units(path)
    
    def __eq__(self, other):
        return isinstance(other, Network) and self.path == other.path
    
    def get_all_recycles(self, all_recycles=None):
        if all_recycles is None:
            all_recycles = set()
        recycle = self.recycle
        if recycle:
            if isinstance(recycle, Stream):
                all_recycles.add(recycle)
            else:
                all_recycles.update(recycle)
        for i in self.path:
            if isinstance(i, Network): i.get_all_recycles(all_recycles)
        return all_recycles
    
    def sort(self, ends):
        path_sources = [PathSource(i, ends) for i in self.path]
        N = len(path_sources)
        if not N: return
        for _ in range(N * N):
            stop = True
            for i in range(N - 1):
                upstream_source = path_sources[i]
                for j in range(i + 1, N):
                    downstream_source = path_sources[j]
                    if upstream_source.downstream_from(downstream_source): 
                        path_sources.remove(downstream_source)
                        path_sources.insert(i, downstream_source)
                        upstream_source = downstream_source
                        stop = False
            if stop: 
                self.path = [i.source for i in path_sources]
                return
        warn(RuntimeWarning('network path could not be determined'))
    
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
        ends : Iterable[:class:`~thermosteam.Stream`]
            Streams that not products, but are ultimately specified through
            process requirements and not by its unit source.
            
        """
        ends = set(ends) if ends else set()
        recycle_ends = ends.copy()
        linear_paths, cyclic_paths_with_recycle = find_linear_and_cyclic_paths_with_recycle(
            feedstock, ends)
        network, *linear_networks = [Network(i) for i in linear_paths]
        for linear_network in linear_networks:
            network.join_linear_network(linear_network) 
        recycle_networks = [Network(*i) for i in cyclic_paths_with_recycle]
        for recycle_network in recycle_networks:
            network.join_recycle_network(recycle_network)
        isa = isinstance
        ends.update(network.streams)
        for feed in feeds:
            if feed in ends or isa(feed.sink, Facility): continue
            downstream_network = cls.from_feedstock(feed, (), ends)
            new_streams = downstream_network.streams
            connections = ends.intersection(new_streams)
            connecting_units = {stream._sink for stream in connections
                                if stream._source and stream._sink}
            ends.update(new_streams)
            N_connections = len(connecting_units)
            if N_connections == 0:
                network._append_network(downstream_network)
            elif N_connections == 1:
                connecting_unit, = connecting_units
                network.join_network_at_unit(downstream_network,
                                             connecting_unit)
            else:
                connecting_unit = network.first_unit(connecting_units)
                network.join_network_at_unit(downstream_network,
                                             connecting_unit)
        recycle_ends.update(network.get_all_recycles())
        recycle_ends.update(bst.utils.products_from_units(network.units))
        network.sort(recycle_ends)
        network.add_process_heat_exchangers()
        return network
    
    def add_process_heat_exchangers(self, excluded=None):
        isa = isinstance
        path = self.path
        if excluded is None: excluded = set()
        for i, u in enumerate(path):
            if isa(u, Unit):
                if isa(u, bst.HXprocess): excluded.add(u)
                for s in u.outs:
                    sink = s.sink
                    if u in excluded: continue
                    if isa(sink, bst.HXprocess) and sink in path[:i] and not any([(sink is i if isa(i, Unit) else sink in i.units) for i in path[i+1:]]):
                        excluded.add(sink)
                        path.insert(i+1, sink)
            else:
                u.add_process_heat_exchangers(excluded)
        if len(path) > 1 and path[-1] is path[0]: path.pop()
    
    @property
    def streams(self):
        return streams_from_units(self.units)
    
    def first_unit(self, units):
        isa = isinstance
        for i in self.path:
            if isa(i, Unit) and i in units:
                return i
            elif isa(i, Network) and not i.units.isdisjoint(units):
                return i.first_unit(units)
        raise ValueError('network does not contain any of the given units') # pragma: no cover
    
    def isdisjoint(self, network):
        return self.units.isdisjoint(network.units)
        
    def join_network_at_unit(self, network, unit):
        isa = isinstance
        path_tuple = tuple(self.path)
        self._remove_overlap(network, path_tuple)
        for index, item in enumerate(self.path):
            if isa(item, Network) and unit in item.units:
                if network.recycle:
                    item.join_network_at_unit(network, unit)
                    self.units.update(network.units)
                else:
                    self._insert_linear_network(index, network)
                return
            elif unit == item:
                self._insert_linear_network(index, network)
                return
        raise RuntimeError(f'{repr(unit)} not in path') # pragma: no cover
    
    def join_linear_network(self, linear_network):
        path = self.path
        path_tuple = tuple(path)
        units = linear_network.units
        self._remove_overlap(linear_network, path_tuple)
        for index, item in enumerate(path_tuple):
            if item in units:
                self._insert_linear_network(index, linear_network)
                return
        self._append_linear_network(linear_network)
    
    def join_recycle_network(self, network):
        if self.recycle_sink is network.recycle_sink:
            # Feed forward scenario
            self.add_recycle(network.recycle)
            network.recycle_sink = network.recycle = None 
            self._add_linear_network(network)
            return
        path = self.path
        isa = isinstance
        path_tuple = tuple(path)
        self._remove_overlap(network, path_tuple)
        subunits = network.units
        for index, i in enumerate(path_tuple):
            if isa(i, Network) and not network.isdisjoint(i):
                i.join_recycle_network(network)
                self.units.update(subunits)
                return
        for index, item in enumerate(path_tuple):
            if isa(item, Unit) and item in subunits:
                self._insert_recycle_network(index, network)
                return
        raise ValueError('networks must have units in common to join') # pragma: no cover
    
    def add_recycle(self, stream):
        recycle = self.recycle
        if recycle is stream: return 
        isa = isinstance
        if isa(recycle, Stream):
            if isa(stream, Stream):
                self.recycle = {self.recycle, stream}
            elif isa(stream, set):
                self.recycle = {self.recycle, *stream}
            else: # pragma: no cover
                raise ValueError(f'recycles must be stream objects; not {type(stream).__name__}')
        elif isa(recycle, set):
            if isa(stream, Stream):
                recycle.add(stream)
            elif isa(stream, set):
                recycle.update(stream)
            else: # pragma: no cover
                raise ValueError(f'recycles must be stream objects; not {type(stream).__name__}')
        else: # pragma: no cover
            raise RuntimeError(f"invalid recycle of type '{type(recycle).__name__}' encountered")
     
    def _remove_overlap(self, network, path_tuple):
        path = self.path
        units = network.units
        isa = isinstance
        for i in path_tuple:
            if (isa(i, Unit) and i in units): path.remove(i)
    
    def _append_linear_network(self, network):
        self.path.extend(network.path)
        self.units.update(network.units)
    
    def _append_recycle_network(self, network):
        self.path.append(network)
        self.units.update(network.units)
    
    def _append_network(self, network):
        if self.recycle:
            cls = type(self)
            new = cls.__new__(cls)
            new.path = self.path; new.units = self.units
            new.recycle = self.recycle; new.recycle_sink = self.recycle_sink
            self.recycle = self.recycle_sink = None
            self.path = [new, network] if network.recycle else [new, *network.path]
            self.units = self.units.union(network.units)
        elif network.recycle:
            self._append_recycle_network(network)
        else:
            self._append_linear_network(network)
    
    def _insert_linear_network(self, index, network):
        path = self.path
        self.path = [*path[:index], *network.path, *path[index:]]
        self.units.update(network.units)
    
    def _insert_recycle_network(self, index, network):
        path = self.path
        path.insert(index, network)
        self.units.update(network.units)
        if len(path) == 1:
            network = path[0]
            if isinstance(network, Network):
                self.path = network.path
                self.recycle = network.recycle
                self.recycle_sink = network.recycle_sink

    def _add_linear_network(self, network):
        path = self.path
        isa = isinstance
        path_tuple = tuple(path)
        self._remove_overlap(network, path_tuple)
        subunits = network.units
        for index, i in enumerate(path_tuple):
            if isa(i, Network) and not network.isdisjoint(i):
                i._add_linear_network(network)
                self.units.update(subunits)
                return
        for index, item in enumerate(path_tuple):
            if isa(item, Unit) and item in subunits:
                self._insert_linear_network(index, network)
                return
        self._append_linear_network(network)
    
    def __repr__(self): # pragma: no cover
        recycle = self.recycle
        if recycle:
            return f"{type(self).__name__}(path={self.path}, recycle={self.recycle})"
        else:
            return f"{type(self).__name__}(path={self.path})"
    
    def _info(self, spaces):
        info = f"{type(self).__name__}("
        spaces += 4 * " "
        end = ',\n' + spaces
        path_info = []
        path = self.path
        isa = isinstance
        info += '\n' + spaces
        for i in path:
            path_info.append(i._info(spaces) if isa(i, Network) else str(i))
        info += '[' + (end + " ").join(path_info) + ']'
        recycle = self.recycle
        if recycle:
            if isinstance(recycle, Stream):
                recycle = recycle._source_info()
            else:
                recycle = ", ".join([i._source_info() for i in recycle])
                recycle = '{' + recycle + '}'
            info += end + f"recycle={recycle})"
        else:
            info += ')'
        return info
    
    def _ipython_display_(self):
        self.show()
    
    def show(self):
        print(self._info(spaces=""))