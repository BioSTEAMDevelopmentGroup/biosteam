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
from thermosteam.utils.piping import (
    AbstractStream, Connection, AbstractInlets, AbstractOutlets, Inlet, Outlet, 
    IgnoreDockingWarnings, ignore_docking_warnings,
)

__all__ = ('AbstractStream', 'MissingStream', 'MockStream', 'Inlets',
           'Outlets', 'InletPort', 'OutletPort', 'StreamPorts', 
           'as_stream', 'as_upstream', 'as_downstream', 
           'materialize_connections', 'ignore_docking_warnings',
           'IgnoreDockingWarnings', 'SuperpositionInlet', 'SuperpositionOutlet',
           'Connection', 'Inlet', 'Outlet')


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
    
class MissingStream:
    """
    Create a MissingStream object that acts as a dummy in Inlets and Outlets
    objects until replaced by an actual Stream object.
    """
    __slots__ = ()
    line = 'Stream'
    ID = 'missing stream'
    T = P = phase = None
    
    def __init__(self, source=None, sink=None):
        self._source = source
        self._sink = sink
    
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
    
    def __bool__(self):
        return False

    def __repr__(self):
        return f'<{type(self).__name__}>'
    
    def __str__(self):
        return self.ID
    
    def show(self):
        print(self._basic_info())


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


# %% System inlets and outlets

class InletPort:
    __slots__ = ('sink', 'index')
    
    @classmethod
    def from_inlet(cls, inlet):
        sink = inlet.sink
        if not sink: raise ValueError(f'stream {inlet} is not an inlet to any unit')
        index = sink.ins.index(inlet)
        return cls(sink, index)
    
    def __init__(self, sink, index):
        self.sink = sink
        self.index = index
      
    def __eq__(self, other):
        return self.sink is other.sink and self.index == other.index  
      
    def _sorting_key(self):
        return (self.sink.ID[1:], self.sink.ID, self.index)
        
    def get_stream(self):
        return self.sink.ins[self.index]
    
    def set_stream(self, stream, stacklevel):
        self.sink.ins._set_stream(self.index, stream, stacklevel+1)
    
    def __str__(self):
        return f"{self.index}-{self.sink}"
    
    def __repr__(self):
        return f"{type(self).__name__}({self.sink}, {self.index})"


class OutletPort:
    __slots__ = ('source', 'index')
    
    @classmethod
    def from_outlet(cls, outlet):
        source = outlet.source
        if not source: raise ValueError(f'stream {outlet} is not an outlet to any unit')
        index = source.outs.index(outlet)
        return cls(source, index)
    
    def __init__(self, source, index):
        self.source = source
        self.index = index
    
    def __eq__(self, other):
        return self.source is other.source and self.index == other.index
    
    def _sorting_key(self):
        return (self.source.ID[1:], self.source.ID[0], self.index)
    
    def get_stream(self):
        return self.source.outs[self.index]
    
    def set_stream(self, stream, stacklevel):
        self.source.outs._set_stream(self.index, stream, stacklevel+1)
    
    def __str__(self):
        return f"{self.source}-{self.index}"
    
    def __repr__(self):
        return f"{type(self).__name__}({self.source}, {self.index})"


class StreamPorts:
    __slots__ = ('_ports',)
    
    @classmethod
    def from_inlets(cls, inlets, sort=None):
        return cls([InletPort.from_inlet(i) for i in inlets], sort)
    
    @classmethod
    def from_outlets(cls, outlets, sort=None):
        return cls([OutletPort.from_outlet(i) for i in outlets], sort)
    
    def __init__(self, ports, sort=None):
        if sort: ports = sorted(ports, key=lambda x: x._sorting_key())
        self._ports = tuple(ports)    
    
    def __bool__(self):
        return bool(self._ports)
        
    def __iter__(self):
        for i in self._ports: yield i.get_stream()
    
    def __len__(self):
        return len(self._ports)
    
    def __getitem__(self, index):
        if isinstance(index, slice):
            return self.__class__(self._ports[index])
        else:
            return self._ports[index].get_stream()
    
    def __setitem__(self, index, item):
        isa = isinstance
        if isa(index, int):
            self._set_stream(index, item, 2)
        elif isa(index, slice):
            self._set_streams(index, item, 2)
        else:
            raise IndexError("Only intergers and slices are valid "
                            f"indices for '{type(self).__name__}' objects")
          
    def _set_stream(self, int, stream, stacklevel):
        self._ports[int].set_stream(stream, stacklevel+1)
    
    def _set_streams(self, slice, streams, stacklevel):
        ports = self._ports[slice]
        stacklevel += 1
        if len(streams) == len(ports):
            for i, j in zip(ports, streams): i.set_stream(j, stacklevel)
        else:
            raise IndexError("number of inlets must match the size of slice")
    
    def __repr__ (self):
        ports = ', '.join([str(i) for i in self._ports])
        return f"[{ports}]"


# %% Auxiliary piping

def superposition_property(name):
    @property
    def p(self):
        return getattr(self.port.get_stream(), name)
    @p.setter
    def p(self, value):
        setattr(self.port.get_stream(), name, value)
        
    return p

def _superposition(cls, parent, port):
    excluded = set([*cls.__dict__, port, '_' + port, 'port'])
    for name in parent.__dict__:
        if name in excluded: continue
        setattr(cls, name, superposition_property(name))
    return cls

def superposition(parent, port):
    return lambda cls: _superposition(cls, parent, port)


@superposition(Stream, 'sink')
class SuperpositionInlet(Stream): 
    """Create a SuperpositionInlet that references an inlet from another 
    unit operation."""
    __slots__ = ()
    
    def __init__(self, port, sink=None):
        self.port = port
        self._sink = sink
      

@superposition(Stream, 'source')
class SuperpositionOutlet(Stream):
    """Create a SuperpositionOutlet that references an outlet from another 
    unit operation."""
    __slots__ = ()
    
    def __init__(self, port, source=None):
        self.port = port
        self._source = source
    