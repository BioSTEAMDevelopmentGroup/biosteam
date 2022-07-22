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

__all__ = (
    'streams_from_units',
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
    'get_inlet_origin',
    'get_streams_from_context_level',
    'get_fresh_process_water_streams',
    'get_unaccounted_waste_streams',
    'get_unaccounted_combustible_waste_streams',
    'get_streams_from_context_level',
)

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

def get_streams_from_context_level(level=None):
    """
    Return streams created in given context level.
    
    Parameters
    ----------
    level : int, optional
        If given, only filter through streams created in the given context 
        level. For example, use:
        * 0: to account for all streams created in the current context level.
        * 1: to account for streams created in the next outer context level.
    
    Examples
    --------
    Get all waste streams created (as registerd in the flowsheet):
        
    >>> import biosteam as bst
    >>> bst.main_flowsheet.clear()
    >>> bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    >>> waste = bst.Stream('waste')
    >>> mixer = bst.Mixer('mixer', outs=waste)
    >>> streams = bst.get_streams_from_context_level()
    >>> assert len(streams) == 1 and streams[0] == waste
    
    Get waste streams created within a context:
        
    >>> with bst.System('sys') as sys:
    ...     new_waste = bst.Stream('new_waste')
    ...     tank = bst.MixTank('tank', outs=new_waste)
    ...     streams = bst.get_streams_from_context_level(level=0)
    ...     assert len(streams) == 1 and streams[0] == new_waste
    
    """
    if level is None:
        streams = list(bst.main_flowsheet.stream)
    else:
        context_levels = bst.main_flowsheet.unit.context_levels
        N_levels = len(context_levels)
        if level >= N_levels:
            streams = bst.main_flowsheet.stream
        else:
            index = N_levels - level - 1
            units = context_levels[index]
            streams = [i for i in streams_from_units(units) if i]
    return streams

def get_unaccounted_combustible_waste_streams(streams=None):
    """
    Return all product streams that do not have a price and have a lower 
    heating value (LHV) greater than 1 kJ/g. These are assumed to be waste 
    streams that can potentially be used to generate heat through combustion 
    and have not been accounted for.
    
    Parameters
    ----------
    streams : Iterable[Stream], optional
        Stream to filter through.
    
    Examples
    --------
    Get all waste streams created (as registerd in the flowsheet):
        
    >>> import biosteam as bst
    >>> bst.main_flowsheet.clear()
    >>> bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    >>> waste = bst.Stream('waste', Ethanol=1)
    >>> mixer = bst.Mixer('mixer', outs=waste)
    >>> streams = bst.get_unaccounted_waste_streams()
    >>> assert len(streams) == 1 and streams[0] == waste
    
    """
    if streams is None: streams = bst.main_flowsheet.stream
    isa = isinstance
    return [
        i for i in streams if 
        not i.price and i.isproduct()
        and not isa(i.source, bst.Facility) 
        and i.LHV / i.F_mass > 1000. 
    ]

def get_unaccounted_waste_streams(streams=None, phase=None):
    """
    Return all product streams that do not have a price at the given phase.
    These are assumed to be waste streams that have not been accounted for 
    (i.e., is not sold and needs to be treated in some way).
    
    Parameters
    ----------
    streams : Iterable[Stream], optional
        Stream to filter through.
    phase : str, optional
        Phase of streams. Defaults to 'l' for liquid.
    
    Examples
    --------
    Get all waste streams created (as registerd in the flowsheet):
        
    >>> import biosteam as bst
    >>> bst.main_flowsheet.clear()
    >>> bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    >>> waste = bst.Stream('waste')
    >>> mixer = bst.Mixer('mixer', outs=waste)
    >>> streams = bst.get_unaccounted_waste_streams()
    >>> assert len(streams) == 1 and streams[0] == waste
    
    """
    if streams is None: streams = bst.main_flowsheet.stream
    if phase is None: phase = 'l'
    isa = isinstance
    return [
        i for i in streams if 
        not i.price and i.isproduct()
        and i.phase==phase
        and not isa(i.source, bst.Facility)
    ]

def get_fresh_process_water_streams(streams=None):
    """
    Return all feed water streams without a price.
    
    Parameters
    ----------
    streams : Iterable[Stream], optional
        Stream to filter through.
        
    Examples
    --------
    Get all fresh process water streams created (as registerd in the flowsheet):
        
    >>> import biosteam as bst
    >>> bst.main_flowsheet.clear()
    >>> bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    >>> water = bst.Stream('water', Water=1)
    >>> mixer = bst.Mixer('mixer', ins=water)
    >>> streams = bst.get_fresh_process_water_streams()
    >>> assert len(streams) == 1 and streams[0] == water
    
    """    
    def only_water(stream):
        chemicals = stream.available_chemicals
        return len(chemicals) == 1 and chemicals[0].CAS == '7732-18-5'
    
    if streams is None: streams = bst.main_flowsheet.stream
    return [
        i for i in streams
        if not i.price and i.isfeed() and i.phase=='l' and only_water(i)
    ]
