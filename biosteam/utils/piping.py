# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
This module includes classes and functions concerning Stream objects.
"""
from thermosteam.utils.decorators import registered_franchise
from thermosteam import Stream
from thermosteam.network import (
    AbstractStream, Connection, AbstractInlets, AbstractOutlets, 
    IgnoreDockingWarnings, ignore_docking_warnings, InletPort, OutletPort,
    StreamPorts, AbstractMissingStream
)

__all__ = ('AbstractStream', 'MissingStream', 'MockStream', 'Inlets',
           'Outlets', 'InletPort', 'OutletPort', 'StreamPorts', 
           'as_stream', 'as_upstream', 'as_downstream', 
           'materialize_connections', 'ignore_docking_warnings',
           'IgnoreDockingWarnings', 'Connection')


# %% Utilities

def as_stream(stream):
    isa = isinstance
    if isa(stream, Stream):
        return stream
    elif isa(stream, str):
        return Stream(stream)
    elif stream is None:
        return MissingStream(None, None)
    
def as_upstream(stream, sink):
    stream = as_stream(stream)
    stream._sink = sink
    return stream
        
def as_downstream(stream, source):
    stream = as_stream(stream)
    stream._source = source
    return stream

def materialize_connections(streams):
    for stream in streams:
        if not stream: stream.materialize_connection()


# %% Special Stream objects    
    
class MissingStream(AbstractMissingStream):
    """
    Create a MissingStream object that acts as a dummy in Inlets and Outlets
    objects until replaced by an actual Stream object.
    """
    __slots__ = ()
    T = P = phase = None
    
    def _get_tooltip_string(self, format, full):
        if format not in ('html', 'svg'): return ''
        return '(empty)'
            
    def materialize_connection(self, ID=None):
        """
        Disconnect this missing stream from any unit operations and 
        replace it with a material stream. 
        """
        source = self._source
        sink = self._sink
        if not (source or sink):
            raise RuntimeError("either a source or a sink is required to "
                               "materialize connection")
        material_stream = Stream(ID, thermo=(source or sink).thermo)
        if source: 
            try: 
                source._outs.replace(self, material_stream)
            except: # Must be an auxlet
                material_stream._source = source
        if sink: 
            try:
                sink._ins.replace(self, material_stream)
            except: # Must be an auxlet
                material_stream._sink = sink
        return material_stream
    
    def get_impact(self, key):
        return 0.
    
    def get_CF(self, key):
        return 0.
    
    def reset_cache(self):
        """Does nothing, MissingStream objects do not contain cache."""
    
    def get_data(self):
        return None
    
    def set_data(self, data): pass
    
    def get_total_flow(self, units):
        return 0.
    
    def get_flow(self, units, key=None):
        return 0.
    
    H = Hf = Hnet = LHV = HHV = Hvap = C = F_mol = F_mass = F_vol = cost = price = 0.
    
    def isempty(self):
        return True
    
    def empty(self): pass


@registered_franchise(Stream)
class MockStream:
    """
    Create a MockStream object that acts as a dummy Stream object that is always
    empty. A MockStream behaves similar to a MissingStream, with the following 
    differences: they have an ID and are registered in the flowsheet as streams,
    they cannot serve as inlets or outlets and cannot be materialized, you can
    set life-cycle characterization factors.
    
    """
    __slots__ = ('_ID', 'characterization_factors', 'price')
    line = 'Stream'
    
    def __init__(self, ID):
        self._ID = ID
        self.characterization_factors = {}
        self.price = 0.
        self._register(ID)

    def _register(self, ID):
        replace_ticket_number = isinstance(ID, int)
        if replace_ticket_number: 
            Stream.ticket_numbers[Stream.ticket_name] = ID
        if ID == "" or replace_ticket_number: 
            registry = Stream.registry
            data = registry.data
            ID = Stream._take_ticket()
            while ID in data: ID = Stream._take_ticket()
            registry.register(ID, self)
        elif ID:
            Stream.registry.register_safely(ID, self) 
        else:
            self._ID = Stream._take_unregistered_ticket()
    
    source = sink = None
    
    def isfeed(self):
        return False
    
    def isproduct(self):
        return False
    
    get_CF = Stream.get_CF
    set_CF = Stream.set_CF
    get_impact = MissingStream.get_impact
    reset_cache = MissingStream.reset_cache
    get_data = MissingStream.get_data
    set_data = MissingStream.set_data
    get_total_flow = MissingStream.get_total_flow
    get_flow = MissingStream.get_flow
    H = MissingStream.H
    Hf = MissingStream.Hf
    Hnet = MissingStream.Hnet
    LHV = MissingStream.LHV
    HHV = MissingStream.HHV
    Hvap = MissingStream.Hvap
    C = MissingStream.C
    F_mol = MissingStream.F_mol
    F_mass = MissingStream.F_mass
    F_vol = MissingStream.F_vol
    isempty = MissingStream.isempty
    empty = MissingStream.empty
    cost = MissingStream.cost
    
    __str__ = Stream.__str__
    __repr__ = Stream.__repr__

# %% Unit operation inlets and outlets

class Inlets(AbstractInlets):
    Stream = Stream
    MissingStream = MissingStream


class Outlets(AbstractOutlets):
    Stream = Stream
    MissingStream = MissingStream
