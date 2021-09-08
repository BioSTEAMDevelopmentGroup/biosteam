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

__all__ = ('streams_from_units',
           'process_streams',
           'streams_from_path', 
           'feeds',
           'products',
           'inlets',
           'outlets',
           'filter_out_missing_streams',
           'sort_feeds_big_to_small',
           'feeds_from_units',
           'products_from_units',
           'get_inlet_origin')

feed_priorities = {}

def inlets(units):
    return set(sum([i._ins for i in units], []))

def outlets(units):
    return set(sum([i._outs for i in units], []))

def streams_from_units(units):
    return set(sum([i._ins + i._outs for i in units], []))

def streams_from_path(path):
    isa = isinstance
    streams = set()
    System = bst.System
    Unit = bst.Unit
    for i in path:
        if isa(i, System):
            streams.update(i.streams)
        elif isa(i, Unit):
            streams.update(i._ins + i._outs)
    return streams

def get_inlet_origin(inlet):
    source = inlet.source
    while source:
        if len(source.ins) == len(source.outs) == 1 and 'processing' not in source.line.lower():
            inlet = source.ins[0]
        elif isinstance(source, bst.HXprocess):
            index = source.outs.index(inlet)
            inlet = source.ins[index]
        else:
            break
        source = inlet.source
    return inlet
        

def process_streams(streams):
    return {i for i in streams if i._source and i._sink}

def feeds(streams):
    return [s for s in streams if not s._source]

def products(streams):
    return [s for s in streams if not s._sink]

def filter_out_missing_streams(streams):
    streams.intersection_update([i for i in streams if i])

def sort_feeds_big_to_small(feeds):
    if feeds:
        def feed_priority(feed):
            if feed in feed_priorities:
                return feed_priorities[feed]
            elif feed:
                return 1. - feed.F_mass / F_mass_max if F_mass_max else 1.
            else:
                return 2.
        F_mass_max = max([i.F_mass for i in feeds])
        feeds.sort(key=feed_priority)

def feeds_from_units(units):
    unit_set = set(units)
    return sum([[i for i in u._ins if i._source not in unit_set]
                 for u in units], [])

def products_from_units(units):
    unit_set = set(units)
    return sum([[i for i in u._outs if i._sink not in unit_set]
                 for u in units], [])

def get_feed_priority(stream):
    if stream.isfeed():
        return feed_priorities.get(stream)
    else:
        raise RuntimeError(f"stream '{stream}' is not a feed")

def set_feed_priority(stream, value):
    if stream.isfeed():
        feed_priorities[stream] = value
    else:
        raise RuntimeError(f"stream '{stream}' is not a feed")
    
bst.Stream.get_feed_priority = get_feed_priority
bst.Stream.set_feed_priority = set_feed_priority