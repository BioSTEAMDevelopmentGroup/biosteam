# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 09:01:47 2019

@author: yoelr
"""
# import biosteam as bst
# from .misc import strtuple

# __all__ = ('Thread', 'build_network')

class RecycleNetwork:
    __slots__ = ('recycle', 'network')
    
    def __int__(self, recycle, network):
        self.recycle = recycle
        self.network = network


class Network:
    __slots__ = ('start', 'middle', 'ends', 'units')
    
    def __init__(self, start=None, middle=None, ends=None):
        raise NotImplementedError("Network is avalilable yet")
        self.start = start or []
        self.middle = middle = middle or []
        self.ends = ends or []
        self.units = set(middle)
        
    def insert_cyclic_path_in_middle(self, cyclic_path, recycle):
        middle = self.middle
        unit = cyclic_path[0]
        start_index = middle.index(unit)
        segment = middle[start_index:]
        end_index = start_index
        for unit_cyclic, unit_linear in zip(cyclic_path, segment):
            if unit_cyclic is unit_linear:
                end_index += 1
            else: 
                break
        self.middle = (middle[:start_index],
                       RecycleNetwork(cyclic_path, recycle),
                       middle[end_index:])        

def form_network(linear_paths):
    *ends, linear_path = linear_paths
    return Network(linear_path, ends=ends)

# def form_networks(linear_paths, cyclic_paths_with_recycle):
#     networks = []
#     for linear_path in linear_paths:
#         for cyclic_path_with_recycle in cyclic_paths_with_recycle:
#             cyclic_path, recycle = cyclic_path_with_recycle
#             recycle_network = RecycleNetwork(cyclic_path, recycle)
#             unit = cyclic_path[0]
#             try:
#                 start_index = linear_path.index(unit)
#             except:
                
#                 continue
#             cyclic_lenght = len(cyclic_path)
#             end_index = cyclic_lenght + start_index
#             if len(linear_path) >= end_index and linear_path[start_index:end_index] == cyclic_path:
#                 network = Network(linear_path,
#                                   linear_path[:start_index],
#                                   recycle_network,
#                                   linear_path[end_index:])
#             networks.append(network)
#     return networks

def find_linear_and_cyclic_paths_with_recycle(feed):
    paths_with_recycle, linear_paths = find_paths_with_and_without_recycle(feed)
    cyclic_paths_with_recycle = set()
    for path_with_recycle in paths_with_recycle:
        cyclic_path_with_recycle = path_with_recycle_to_cyclic_path_with_recycle(path_with_recycle)
        cyclic_paths_with_recycle.add(cyclic_path_with_recycle)
    return simplify_linear_paths(linear_paths), cyclic_paths_with_recycle

def find_paths_with_and_without_recycle(feed):
    path = []
    paths_without_recycle  = set()
    paths_with_recycle = set()
    has_recycle = fill_path(feed, path, paths_with_recycle, paths_without_recycle)
    path = tuple(path)
    if not has_recycle:
        paths_without_recycle.add(path)
    return paths_with_recycle, paths_without_recycle

def fill_path(feed, path, paths_with_recycle, paths_without_recycle):
    unit = feed.sink
    if not unit: return 
    if unit in path: 
        path_with_recycle = tuple(path), feed
        paths_with_recycle.add(path_with_recycle)
        has_recycle = True
        return 
    path.append(unit)
    outlet, *other_outlets = unit.outs
    fill_path(outlet, path, paths_with_recycle, paths_without_recycle)
    for outlet in other_outlets:
        new_path = path.copy()
        if not fill_path(outlet, new_path, paths_with_recycle, paths_without_recycle):
            paths_without_recycle.add(tuple(new_path))
    has_recycle = False
    return has_recycle

def path_with_recycle_to_cyclic_path_with_recycle(path_with_recycle):
    path, recycle = path_with_recycle
    unit = recycle.sink
    recycle_index = path.index(unit)
    return (path[recycle_index:], recycle)

def simplify_linear_paths(linear_paths):
    simplified_linear_paths = []
    linear_paths = sorted(linear_paths, key=len)
    while linear_paths:
        bigger_path, *linear_paths = linear_paths
        simplified_path = simplify_linear_path(bigger_path, linear_paths)
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

# def simplify_paths(paths):
#     simplified_paths = set()
#     while paths:
#         path = paths.pop()
#         if is_unique_path(path, paths):
#             simplified_paths.add(path)
#     return simplified_paths
    
# def is_unique_path(path, other_paths):
#     for other_path in other_paths:
#         if path in other_path: return False
#     return True
    
# def build_network(sources, sinks=None):
#     thread = Thread.from_units(sources, sinks)
#     return thread.build_network()

# def traverse(iterable):
#     L = list
#     inst = isinstance
#     path = []
#     for iter in iterable:
#         if inst(iter, L):
#             path += traverse(iter)
#         else:
#             path.append(iter)
#     return path

# def is_ascending(nested_index, other):
#     for depth, jk in enumerate(zip(nested_index, other)): 
#         j, k = jk
#         if j != k: break
#     return depth, j < k

# def get_nested_index(iterable, item):
#     L = list
#     inst = isinstance
#     index = []
#     for i, iter in enumerate(iterable):
#         if inst(iter, L):
#             index += get_nested_index(iter, item)
#         else:
#             if iter is item:
#                 index.append(i)
#                 break
#     return index

# def get_nested_item(iterable, nested_index):
#     for i in nested_index:
#         iterable = iterable[i]
#     return iterable

# class Thread:
#     __slots__ = ('path', 'source', 'sink')
    
#     def __init__(self, path, source, sink):
#         #: list[Unit or ...] Ordered units in the thread.
#         self.path = path
        
#         #: [Unit or None] Origin of thread.
#         self.source = source
        
#         #: [Unit or None] End of thread
#         self.sink = sink
        
#     @classmethod
#     def single_thread(cls, start, source, sinks, thread_args, past_sources):
#         path = [start]
#         next_unit = start
#         sinks = set(sinks)
#         neighbors = set()
#         while next_unit:
#             next_units = [i.sink for i in next_unit.outs if i.sink]
#             if next_units:
#                 last_unit = next_unit
#                 next_unit, *other_units = next_units
#                 while True:
#                     new_neighbors = [i for i in other_units if i not in neighbors]
#                     if not new_neighbors: break
#                     neighbors.update(new_neighbors)
#                     if next_unit not in past_sources:
#                         past_sources.add(next_unit)
#                         thread_args.append((next_unit, last_unit))
#                     next_unit, *other_units = new_neighbors
#                 if next_unit in sinks:
#                     break
#                 elif next_unit in path:
#                     index = path.index(next_unit)
#                     recycle_units = path[index:]
#                     path = path[:index]
#                     path.append(recycle_units)
#                     break
#                 path.append(next_unit)
#             else:
#                 break
#         self = cls(path, source, sink=next_unit)
#         sinks.update(traverse(path))
#         return self
        
#     @classmethod
#     def web(cls, start, source, sinks, thread_args, threads, past_sources):
#         more = []
#         threads.extend([cls.single_thread(*i, sinks, more, past_sources) for i in thread_args])
#         if more:
#             cls.web(start, source, sinks, more, threads, past_sources)
#         return threads
        
#     @classmethod
#     def from_unit(cls, start, source, sinks):
#         thread_args = []
#         past_sources = set()
#         self = cls.single_thread(start, source, sinks, thread_args, past_sources)
#         threads = []
#         cls.web(start, source, sinks, thread_args, threads, past_sources)
#         for thread in threads: self.knit(thread)
#         return self

#     @classmethod
#     def from_units(cls, sources, sinks):
#         self, *threads = [cls.from_unit(i, None, sinks) for i in sources]
#         for thread in threads: self.knit(thread)
    
#     def knit(self, other):
#         source = other.source
#         sink = other.sink
#         path = self.path
#         source_index = get_nested_index(path, source)
#         sink_index = get_nested_index(path, sink)
#         self.show()
#         other.show()
#         if not sink_index:
#             self.path.extend(other.path)
#             self.show()
#             return
#         depth, check = is_ascending(source_index, sink_index)
#         if check:
#             *nest_index, item_index = sink_index
#             sink_nest = get_nested_item(path, nest_index)
#             sink_nest[:] = sink_nest[:item_index] + other.path + sink_nest[item_index:]
#         else:
#             *nest_index, sink_item_index = sink_index[:depth+1]
#             source_item_index = source_index[depth]
#             nest = get_nested_item(path, nest_index)
#             stop = source_item_index+1
#             segment = slice(sink_item_index, stop)
#             nest[segment] = [[nest[segment] + [other.path]]] + nest[stop:]
#         self.show()
    
#     def build_network(self):
#         units = self.path + [self.sink]
#         network = []
#         for segment in units:
#             if isinstance(segment, list):
#                 source = segment[0]
#                 for recycle in segment[-1].outs:
#                     if recycle.sink is source: break
#                 system = bst.System(segment, recycle=recycle)
#                 network.append(system)
#             else:
#                 network.append(segment)
#         return network

#     def show(self):
#         path = strtuple(self.path)
#         i = 1; last_i = 0
#         while True:
#             i += 2
#             i = path.find(', ', i)
#             i_next = path.find(', ', i+2)
#             if (i_next-last_i) > 35:
#                 path = (path[:i] + '%' + path[i:])
#                 last_i = i
#             elif i == -1: break
#         path = path.replace('%, ', ',\n'+' '*8)
#         print(f"{type(self).__name__}:\n"
#              + f" source: {self.source}\n"
#              + f" sink: {self.sink}\n"
#              + f" path: {path}")

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

