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
from warnings import warn
from typing import NamedTuple

__all__ = ('MissingStream', 'MockStream', 'StreamSequence',
           'Inlets', 'Outlets', 'InletPort', 'OutletPort', 'StreamPorts', 
           'as_stream', 'as_upstream', 'as_downstream', 
           'materialize_connections', 'ignore_docking_warnings',
           'IgnoreDockingWarnings', 'Connection', 'InletPipe', 'OutletPipe',
           'mark_disjunction', 'unmark_disjunction')


# %% Utilities for docking

def n_missing(ub, N):
    if ub < N: raise RuntimeError(f"at most {ub} streams expected; {N} given")
    return ub - N

DOCKING_WARNINGS = True

def ignore_docking_warnings(f):
    def g(*args, **kwargs):
        global DOCKING_WARNINGS
        warn = DOCKING_WARNINGS
        DOCKING_WARNINGS = False
        try:
            return f(*args, **kwargs)
        finally:
            DOCKING_WARNINGS = warn
    g.__name__ = f.__name__
    return g

class IgnoreDockingWarnings:
    """
    Ignore docking warnings within a context by using an IgnoreDockingWarnings object.
    
    """
    __slots__ = ('original_value',)
    
    def __enter__(self): 
        global DOCKING_WARNINGS
        self.original_value = DOCKING_WARNINGS
        DOCKING_WARNINGS = False
        return 
    
    def __exit__(self, type, exception, traceback):
        global DOCKING_WARNINGS
        DOCKING_WARNINGS = self.original_value
        if exception: raise exception
        
# %% Configuration bookkeeping

class Connection(NamedTuple):
    source: object
    source_index: int
    stream: object
    sink_index: int
    sink: object
    
    @ignore_docking_warnings
    def reconnect(self):
        # Does not attempt to connect auxiliaries with owners (which should not be possible)
        source = self.source
        sink = self.sink
        if source:
            if not (sink and getattr(sink, '_owner', None) is source):
                source.outs[self.source_index] = self.stream
        else:
            self.stream.disconnect_source()
        if sink:
            if not (source and getattr(source, '_owner', None) is sink):
                sink.ins[self.sink_index] = self.stream
        else:
            self.stream.disconnect_sink()


# %% Piping notation

class InletPipe:
    """
    Create an InletPipe object that connects an stream to a unit using -pipe-
    notation:
    
    Parameters
    ----------
    stream : Stream
    index : int
        
    Examples
    --------
    First create a stream and a Mixer:
    
    >>> from biosteam import Stream, Mixer, settings
    >>> settings.set_thermo(['Water'])
    >>> stream = Stream('s1')
    >>> unit = Mixer('M1', outs=('out'))
    
    Inlet pipes are created using -pipe- notation:
        
    >>> stream-1
    <InletPipe: s1-1>
    
    Use pipe notation to connect the stream:
    
    >>> stream-1-unit
    <Mixer: M1>
    >>> unit.show()
    Mixer: M1
    ins...
    [0] missing stream
    [1] s1
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow: 0
    outs...
    [0] out
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow: 0
    
    """
    __slots__ = ('stream', 'index')
    def __init__(self, stream, index):
        self.stream = stream
        self.index = index

    # Forward pipping
    def __sub__(self, unit):
        unit.ins[self.index] = self.stream
        return unit
    
    # Backward pipping
    __pow__ = __sub__
    
    def __repr__(self):
        return '<' + type(self).__name__ + ': ' + self.stream.ID + '-' + str(self.index) + '>'


class OutletPipe:
    """
    Create an OutletPipe object that connects an stream to a unit using -pipe-
    notation:
    
    Parameters
    ----------
    stream : Stream
    index : int
    
    Examples
    --------
    First create a stream and a Mixer:
    
    >>> from biosteam import Stream, Mixer, settings
    >>> settings.set_thermo(['Water'])
    >>> stream = Stream('s1')
    >>> unit = Mixer('M1')
    
    Outlet pipes are created using -pipe- notation:
        
    >>> 1**stream
    <OutletPipe: 1-s1>
    
    Use -pipe- notation to connect the stream:
    
    >>> unit**0**stream
    <Mixer: M1>
    >>> unit.show()
    Mixer: M1
    ins...
    [0] missing stream
    [1] missing stream
    outs...
    [0] s1
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow: 0
    
    """
    __slots__ = ('stream', 'index')
    def __init__(self, stream, index):
        self.stream = stream
        self.index = index

    # Forward pipping
    def __rsub__(self, unit):
        unit.outs[self.index] = self.stream
        return unit
    
    # Backward pipping
    __rpow__ = __rsub__
    
    def __repr__(self):
        return '<' + type(self).__name__ + ': ' + str(self.index) + '-' + self.stream.ID + '>'


# %% System inlets and outlets

class InletPort:
    __slots__ = ('sink', 'index')
    
    @classmethod
    def from_inlet(cls, inlet):
        if hasattr(inlet, 'port'): return inlet.port
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
        if hasattr(outlet, 'port'): return outlet.port
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
    
    def __add__(self, other):
        return [i.get_stream() for i in self._ports] + other
    
    def __radd__(self, other):
        return other + [i.get_stream() for i in self._ports]
    
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

# %% Stream piping notation and connections

def _stream_sub(self, index):
    if isinstance(index, int):
        return InletPipe(self, index)
    elif hasattr(index, 'separate_out'):
        new = self.copy()
        new.separate_out(index)   
        return new
    else:
        return NotImplemented

def _stream_rsub(self, index):
    if isinstance(index, int):
        return OutletPipe(self, index)
    else:
        return NotImplemented
    
def _stream_pow(self, index):
    if isinstance(index, int):
        return InletPipe(self, index)
    else:
        return NotImplemented

def _stream_rpow(self, index):
    if isinstance(index, int):
        return OutletPipe(self, index)
    else:
        return NotImplemented

def _stream_get_connection(self, junction=None):
    if junction is None: junction = True
    source = self._source
    sink = self._sink
    if not junction and getattr(source, 'line', None) == 'Junction':
        stream = source._ins[0]
        source = stream._source
    else:
        stream = self
    try:
        source_index = source._outs.index(stream) if source else None
    except ValueError: # Auxiliary streams are not in inlets nor outlets  
        try:
            auxname = sink.ID
            source_index = source._auxout_index[auxname]
        except:
            source_index = -1
    if not junction and getattr(source, 'line', None) == 'Junction':
        stream = sink._outs[0]
        sink = stream._sink
    else:
        stream = self
    try:
        sink_index = sink._ins.index(stream) if sink else None
    except ValueError: # Auxiliary streams are not in inlets nor outlets
        try:
            auxname = source.ID
            sink_index = sink._auxin_index[auxname]
        except:
            sink_index = -1
    return Connection(source, source_index, self, sink_index, sink)

Stream.__sub__ = _stream_sub
Stream.__rsub__ = _stream_rsub
Stream.__pow__ = _stream_pow
Stream.__rpow__ = _stream_rpow
Stream.get_connection = _stream_get_connection

# %% Special network configuration

#: Outlet ports of streams marked as disjunctions.
disjunctions = []

def mark_disjunction(stream):
    """Mark stream as deterministic after a linear flowsheet simulation.
    In other words, the stream flow rates and thermal condition will not change
    when solving recycle loops. This will prevent the system from
    forming recycle loops about this stream."""
    port = OutletPort.from_outlet(stream)
    if port not in disjunctions:
        disjunctions.append(port)

def unmark_disjunction(stream):
    port = OutletPort.from_outlet(stream)
    if port in disjunctions:
        disjunctions.remove(port)


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
    __slots__ = ('_source', '_sink')
    line = 'Stream'
    ID = 'missing stream'
    get_connection = Stream.get_connection
    disconnect_source = Stream.disconnect_source
    disconnect_sink = Stream.disconnect_sink
    disconnect = Stream.disconnect
    
    def __init__(self, source=None, sink=None):
        self._source = source
        self._sink = sink
    
    @property
    def source(self):
        return self._source
    @property
    def sink(self):
        return self._sink

    def __bool__(self):
        return False

    def __repr__(self):
        return f'<{type(self).__name__}>'

    def __str__(self):
        return self.ID

    def show(self):
        print(self._basic_info())

    _basic_info = Stream._basic_info
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

class StreamSequence:
    """
    Abstract class for a sequence of streams for a unit.
    
    Abstract methods:
        * _dock(self, stream) -> Stream
        * _redock(self, stream) -> Stream
        * _undock(self) -> None
        * _load_missing_stream(self)
    
    """
    __slots__ = ('_size', '_streams', '_fixed_size')
    MissingStream = MissingStream
    Stream = Stream
    
    def __init_subclass__(cls):
        cls.stream_types = (cls.Stream, cls.MissingStream)
    
    def __init__(self, size, streams, thermo, fixed_size, stacklevel):
        self._size = size
        self._fixed_size = fixed_size
        Stream = self.Stream
        stream_types = self.stream_types
        dock = self._dock
        redock = self._redock
        if streams == ():
            self._streams = [dock(Stream('-', thermo=thermo)) for i in range(size)]
        else:
            isa = isinstance
            if fixed_size:
                self._initialize_missing_streams()
                if streams is not None:
                    if isa(streams, str):
                        self._streams[0] = dock(Stream(streams or '-', thermo=thermo))
                    elif isa(streams, stream_types):
                        self._streams[0] = redock(streams, stacklevel)
                    else:
                        N = len(streams)
                        n_missing(size, N) # Make sure size is not too big
                        loaded_streams = []
                        for i in streams:
                            if isa(i, Stream):
                                s = redock(i, stacklevel)
                            elif i is None:
                                s = self._create_missing_stream()
                            else:
                                s = Stream(i or '-', thermo=thermo)
                                dock(s)
                            loaded_streams.append(s)
                        self._streams[:N] = loaded_streams
            elif streams is not None:
                if isa(streams, str):
                    self._streams = [dock(Stream(streams or '-', thermo=thermo))]
                elif isinstance(streams, stream_types):
                    self._streams = [redock(streams, stacklevel)]
                else:
                    self._streams = loaded_streams = []
                    for i in streams:
                        if isa(i, Stream):
                            s = redock(i, stacklevel)
                        elif i is None:
                            s = self._create_missing_stream()
                        else:
                            s = Stream(i or '-', thermo=thermo)
                            dock(s)
                        loaded_streams.append(s)
            else:
                self._initialize_missing_streams()
        
    def _create_missing_stream(self):
        return self.MissingStream(None, None)
        
    def _create_N_missing_streams(self, N):
        return [self._create_missing_stream() for i in range(N)]
    
    def _initialize_missing_streams(self):
        #: All input streams
        self._streams = self._create_N_missing_streams(self._size)
        
    def __add__(self, other):
        return self._streams + other
    def __radd__(self, other):
        return other + self._streams
    
    # DO NOT DELETE: These should be implemented by child class
    # def _dock(self, stream): return stream
    # def _redock(self, stream, stacklevel): return stream
    # def _undock(self, stream): pass

    def _set_streams(self, slice, streams, stacklevel):
        streams = [self._as_stream(i) for i in streams]
        all_streams = self._streams
        for stream in all_streams[slice]: self._undock(stream)
        all_streams[slice] = streams
        stacklevel += 1
        for stream in all_streams: self._redock(stream, stacklevel)
        if self._fixed_size:
            size = self._size
            N_streams = len(all_streams)
            if N_streams < size:
                N_missing = n_missing(size, N_streams)
                if N_missing:
                    all_streams[N_streams: size] = self._create_N_missing_streams(N_missing)
       
    def _as_stream(self, stream):
        if stream is None:
            stream = self._create_missing_stream()
        elif not isinstance(stream, (Stream, MissingStream)):
            raise TypeError(
                f"'{type(self).__name__}' object can only contain "
                f"streams; not '{type(stream).__name__}' objects"
            )
        return stream
       
    @property
    def size(self):
        return self._streams.__len__()
    
    def __len__(self):
        return self._streams.__len__()
    
    def __bool__(self):
        return bool(self._streams)
    
    def _set_stream(self, int, stream, stacklevel):
        stream = self._as_stream(stream)
        try:
            old_stream = self._streams[int]
        except IndexError as e:
            if int >= self.size and not self._fixed_size: 
                self._streams.append(self._redock(stream, stacklevel+1))
            else:
                raise e from None
        else:
            self._undock(old_stream)
            self._streams[int] = self._redock(stream, stacklevel+1)
    
    def empty(self):
        for i in self._streams: self._undock(i)
        self._initialize_missing_streams()
    
    def insert(self, index, stream):
        if self._fixed_size: 
            raise RuntimeError(f"size of '{type(self).__name__}' object is fixed")
        self._undock(stream)
        self._dock(stream)
        self._streams.insert(index, stream)
    
    def append(self, stream):
        if self._fixed_size: 
            raise RuntimeError(f"size of '{type(self).__name__}' object is fixed")
        self._undock(stream)
        self._dock(stream)
        self._streams.append(stream)
    
    def extend(self, streams):
        if self._fixed_size: 
            raise RuntimeError(f"size of '{type(self).__name__}' object is fixed")
        for i in streams:
            self._undock(i)
            self._dock(i)
            self._streams.append(i)
    
    def replace(self, stream, other_stream):
        index = self.index(stream)
        self[index] = other_stream

    def index(self, stream):
        try:
            return self._streams.index(stream)
        except Exception as e:
            imol = stream._imol
            for n, i in enumerate(self._streams):
                if i and i._imol is imol:
                    return n
            else:
                raise e from None

    def pop(self, index):
        streams = self._streams
        if self._fixed_size:
            stream = streams[index]
            missing_stream = self._create_missing_stream()
            self.replace(stream, missing_stream)
        else:
            stream = streams.pop(index)
        return stream

    def remove(self, stream):
        self.replace(stream, self._create_missing_stream())
        
    def clear(self):
        if self._fixed_size:
            self._initialize_missing_streams()
        else:
            for i in self._streams: self._undock(i)
            self._streams.clear()
    
    def reverse(self):
        self.streams.reverse()
    
    def __iter__(self):
        return iter(self._streams)
    
    def __getitem__(self, index):
        return self._streams[index]
    
    def __setitem__(self, index, item):
        isa = isinstance
        if isa(index, int):
            self._set_stream(index, item, 2)
        elif isa(index, slice):
            self._set_streams(index, item, 2)
        else:
            raise IndexError("Only intergers and slices are valid "
                             f"indices for '{type(self).__name__}' objects")
    
    def __repr__(self):
        return repr(self._streams)


class Inlets(StreamSequence):
    """Create an Inlets object which serves as input streams for a unit."""
    __slots__ = ('_sink',)
    Stream = Stream
    MissingStream = MissingStream
    
    def __init__(self, sink, size, streams, thermo, fixed_size, stacklevel):
        self._sink = sink
        super().__init__(size, streams, thermo, fixed_size, stacklevel)
    
    @property
    def sink(self):
        return self._sink
    
    def _create_missing_stream(self):
        return self.MissingStream(None, self._sink)
    
    def _dock(self, stream): 
        stream._sink = self._sink
        return stream

    def _redock(self, stream, stacklevel): 
        sink = stream._sink
        if sink:
            ins = sink._ins
            if ins is not self:
                if stream in ins:
                    ins.remove(stream)
                    stream._sink = new_sink = self._sink
                    if (DOCKING_WARNINGS 
                        and sink.ID and new_sink.ID
                        and sink.ID != new_sink.ID):
                        warn(f"undocked inlet {stream} from {sink}; "
                             f"{stream} is now docked at {new_sink}", 
                             RuntimeWarning, stacklevel + 1)
                else:
                    stream._sink = self._sink
        else:
            stream._sink = self._sink
        return stream
    
    def _undock(self, stream): 
        stream._sink = None
    
        
class Outlets(StreamSequence):
    """Create an Outlets object which serves as output streams for units."""
    __slots__ = ('_source',)
    Stream = Stream
    MissingStream = MissingStream
    
    def __init__(self, source, size, streams, thermo, fixed_size, stacklevel):
        self._source = source
        super().__init__(size, streams, thermo, fixed_size, stacklevel)
    
    @property
    def source(self):
        return self._source
    
    def _create_missing_stream(self):
        return self.MissingStream(self._source, None)
    
    def _dock(self, stream): 
        stream._source = self._source
        return stream

    def _redock(self, stream, stacklevel): 
        source = stream._source
        if source:
            outs = source._outs
            if outs is not self:
                if stream in outs:
                    outs.remove(stream)
                    stream._source = new_source = self._source
                    if (DOCKING_WARNINGS 
                        and source.ID and new_source.ID
                        and source.ID != new_source.ID):
                        warn(f"undocked outlet {stream} from {source}; "
                             f"{stream} is now docked at {new_source}", 
                             RuntimeWarning, stacklevel + 1)
                else:
                    stream._source = self._source
        else:
            stream._source = self._source
        return stream
    
    def _undock(self, stream): 
        stream._source = None
