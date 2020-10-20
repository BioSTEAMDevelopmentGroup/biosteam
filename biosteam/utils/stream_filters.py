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

__all__ = ('streams_from_units',
           'process_streams',
           'streams_from_path', 
           'feeds',
           'products',
           'inlets',
           'outlets',
           'filter_out_missing_streams',)

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
            streams.add(i.streams)
        elif isa(i, Unit):
            streams.update(i._ins + i._outs)
    return streams

def process_streams(streams):
    return {i for i in streams if i._source and i._sink}

def feeds(streams):
    return {s for s in streams if not s._source}

def products(streams):
    return {s for s in streams if not s._sink}

def filter_out_missing_streams(streams):
    streams.intersection_update([i for i in streams if i])

