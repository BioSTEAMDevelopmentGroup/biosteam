# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 09:01:47 2019

@author: yoelr
"""
import biosteam as bst
from .misc import strtuple

__all__ = ('_add_upstream_neighbors', '_add_downstream_neighbors',
           'Thread', 'build_network')

def _add_upstream_neighbors(unit, set):
    """Add upsteam neighboring units to set."""
    for s in unit._ins:
        u_source = s._source
        if u_source: set.add(u_source)

def _add_downstream_neighbors(unit, set):
    """Add downstream neighboring units to set."""
    for s in unit._outs:
        u_sink = s._sink
        if u_sink: set.add(u_sink)

def build_network(sources, sinks=None):
    thread = Thread.from_units(sources, sinks)
    return thread.build_network()

def traverse(iterable):
    L = list
    inst = isinstance
    path = []
    for iter in iterable:
        if inst(iter, L):
            path += traverse(iter)
        else:
            path.append(iter)
    return path

def is_ascending(nested_index, other):
    for depth, jk in enumerate(zip(nested_index, other)): 
        j, k = jk
        if j != k: break
    return depth, j < k

def get_nested_index(iterable, item):
    L = list
    inst = isinstance
    index = []
    for i, iter in enumerate(iterable):
        if inst(iter, L):
            index += get_nested_index(iter, item)
        else:
            if iter is item:
                index.append(i)
                break
    return index

def get_nested_item(iterable, nested_index):
    for i in nested_index:
        iterable = iterable[i]
    return iterable

class Thread:
    __slots__ = ('path', 'source', 'sink')
    
    def __init__(self, path, source, sink):
        #: list[Unit or ...] Ordered units in the thread.
        self.path = path
        
        #: [Unit or None] Origin of thread.
        self.source = source
        
        #: [Unit or None] End of thread
        self.sink = sink
        
    @classmethod
    def single_thread(cls, start, source, sinks, thread_args, past_sources):
        path = [start]
        next_unit = start
        sinks = set(sinks)
        neighbors = set()
        while next_unit:
            next_units = [i.sink for i in next_unit.outs if i.sink]
            if next_units:
                last_unit = next_unit
                next_unit, *other_units = next_units
                while True:
                    new_neighbors = [i for i in other_units if i not in neighbors]
                    if not new_neighbors: break
                    neighbors.update(new_neighbors)
                    if next_unit not in past_sources:
                        past_sources.add(next_unit)
                        thread_args.append((next_unit, last_unit))
                    next_unit, *other_units = new_neighbors
                if next_unit in sinks:
                    break
                elif next_unit in path:
                    index = path.index(next_unit)
                    recycle_units = path[index:]
                    path = path[:index]
                    path.append(recycle_units)
                    break
                path.append(next_unit)
            else:
                break
        self = cls(path, source, sink=next_unit)
        sinks.update(traverse(path))
        return self
        
    @classmethod
    def web(cls, start, source, sinks, thread_args, threads, past_sources):
        more = []
        threads.extend([cls.single_thread(*i, sinks, more, past_sources) for i in thread_args])
        if more:
            cls.web(start, source, sinks, more, threads, past_sources)
        return threads
        
    @classmethod
    def from_unit(cls, start, source, sinks):
        thread_args = []
        past_sources = set()
        self = cls.single_thread(start, source, sinks, thread_args, past_sources)
        threads = []
        cls.web(start, source, sinks, thread_args, threads, past_sources)
        for thread in threads: self.knit(thread)
        return self

    @classmethod
    def from_units(cls, sources, sinks):
        self, *threads = [cls.from_unit(i, None, sinks) for i in sources]
        for thread in threads: self.knit(thread)
    
    def knit(self, other):
        source = other.source
        sink = other.sink
        path = self.path
        source_index = get_nested_index(path, source)
        sink_index = get_nested_index(path, sink)
        self.show()
        other.show()
        if not sink_index:
            self.path.extend(other.path)
            self.show()
            return
        depth, check = is_ascending(source_index, sink_index)
        if check:
            *nest_index, item_index = sink_index
            sink_nest = get_nested_item(path, nest_index)
            sink_nest[:] = sink_nest[:item_index] + other.path + sink_nest[item_index:]
        else:
            *nest_index, sink_item_index = sink_index[:depth+1]
            source_item_index = source_index[depth]
            nest = get_nested_item(path, nest_index)
            stop = source_item_index+1
            segment = slice(sink_item_index, stop)
            nest[segment] = [[nest[segment] + [other.path]]] + nest[stop:]
        self.show()
    
    def build_network(self):
        units = self.path + [self.sink]
        network = []
        for segment in units:
            if isinstance(segment, list):
                source = segment[0]
                for recycle in segment[-1].outs:
                    if recycle.sink is source: break
                system = bst.System(segment, recycle=recycle)
                network.append(system)
            else:
                network.append(segment)
        return network

    def show(self):
        path = strtuple(self.path)
        i = 1; last_i = 0
        while True:
            i += 2
            i = path.find(', ', i)
            i_next = path.find(', ', i+2)
            if (i_next-last_i) > 35:
                path = (path[:i] + '%' + path[i:])
                last_i = i
            elif i == -1: break
        path = path.replace('%, ', ',\n'+' '*8)
        print(f"{type(self).__name__}:\n"
             + f" source: {self.source}\n"
             + f" sink: {self.sink}\n"
             + f" path: {path}")

# Example path
# >>> threads = [
# ...     Thread(units=[U0, U1, [U2, U3, U4, U5]],   source=None,   sink=None),
# ...     Thread(units=[U6, U7],                     source=U3,     sink=U5),
# ...     Thread(units=[U8, U9, U0],                 source=U6,     sink=U0),
# ... ]
# >>> path = Thread.path(threads)
# >>> path
# [System(path=[U0, U1, System(path=[U2, U3, U4, U6, U7, U5], recycle=stream), U8, U9])]

# First knit threads toghether

# step 1:
# [Thread(units=[U0, U1, [U2, U3, U4, U6, U7, U5]], source=None, sink=None),
#  Thread(units=[U8, U9, U0],                       source=U3,   sink=U0]
# step 2:
# [Thread(units=[[U0, U1, [U2, U3, U4, U6, U7, U5], U8, U9]], source=None, sink=None)]

# Create systems and subsystems
# [System(path=[U0, U1, System(path=[U2, U3, U4, U6, U7, U5], recycle=stream), U8, U9])]

