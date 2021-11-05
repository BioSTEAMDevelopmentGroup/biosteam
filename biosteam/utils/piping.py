# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
This module includes classes and functions concerning Stream objects.
"""
from thermosteam.utils.decorators import registered_franchise
from thermosteam import Stream
from collections import namedtuple
from warnings import warn
__all__ = ('MissingStream', 'MockStream', 'Inlets', 'Outlets', 'Sink', 'Source',
           'InletPort', 'OutletPort', 'StreamPorts', 'Connection', 
           'as_stream', 'as_upstream', 'as_downstream', 
           'materialize_connections', 'ignore_docking_warnings')

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

# %% Utilities

def pipe_info(source, sink):
    """Return stream information header."""
    # First line
    if source is None:
        source = ''
    else:
        source = f' from {repr(source)}'
    if sink is None:
        sink = ''
    else:
        sink = f' to {repr(sink)}'
    return f"{source}{sink}"

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


# %% Dummy Stream object

class MissingStream:
    """
    Create a MissingStream object that acts as a dummy in Inlets and Outlets
    objects until replaced by an actual Stream object.
    """
    __slots__ = ('_source', '_sink')
    line = 'Stream'
    
    disconnect = Stream.disconnect
    disconnect_source = Stream.disconnect_source
    disconnect_sink = Stream.disconnect_sink
    isfeed = Stream.isfeed
    isproduct = Stream.isproduct
    
    def __init__(self, source, sink):
        self._source = source
        self._sink = sink
    
    def get_connection(self):
        self = self.materialize_connection()
        return self.get_connection()
    
    def materialize_connection(self, ID=""):
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
        if source: source._outs.replace(self, material_stream)
        if sink: sink._ins.replace(self, material_stream)
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
    
    @property
    def link(self):
        return None
    @property
    def H(self):
        return 0.
    @property
    def Hf(self):
        return 0.
    @property
    def Hnet(self):
        return 0.
    @property
    def LHV(self):
        return 0.
    @property
    def HHV(self):
        return 0.
    @property
    def Hvap(self):
        return 0.
    @property
    def C(self):
        return 0.
    @property
    def F_mol(self):
        return 0.
    @property
    def F_mass(self):
        return 0.
    @property
    def F_vol(self):
        return 0.
    @property
    def cost(self):
        return 0.
    @property
    def price(self):
        return 0.
    
    def isempty(self):
        return True
    
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
        return 'missing stream'

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
    
    @property
    def ID(self):
        """Unique identification (str). If set as '', it will choose a default ID."""
        return self._ID

    @ID.setter
    def ID(self, ID):
        Stream._register(self, ID)
    
    get_CF = Stream.get_CF
    set_CF = Stream.set_CF
    get_impact = MissingStream.get_impact
    reset_cache = MissingStream.reset_cache
    get_data = MissingStream.get_data
    set_data = MissingStream.set_data
    get_total_flow = MissingStream.get_total_flow
    get_flow = MissingStream.get_flow
    link = MissingStream.link
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
    cost = MissingStream.cost
    
    __str__ = Stream.__str__
    __repr__ = Stream.__repr__

# %% Utilities

def n_missing(ub, N):
    if ub < N: raise RuntimeError(f"size of streams exceeds {ub}")
    return ub - N

# %% List objects for input and output streams

class StreamSequence:
    """
    Abstract class for a sequence of streams for a Unit object.
    
    Abstract methods:
        * _dock(self, stream) -> Stream
        * _redock(self, stream) -> Stream
        * _undock(self) -> None
        * _load_missing_stream(self)
    
    """
    __slots__ = ('_size', '_streams', '_fixed_size')
        
    def __init__(self, size, streams, thermo, fixed_size, stacklevel):
        self._size = size
        self._fixed_size = fixed_size
        dock = self._dock
        redock = self._redock
        if streams == ():
            self._streams = [dock(Stream(thermo=thermo)) for i in range(size)]
        else:
            isa = isinstance
            stream_types = (Stream, MissingStream)
            if fixed_size:
                self._initialize_missing_streams()
                if streams is not None:
                    if isa(streams, str):
                        self._streams[0] = dock(Stream(streams, thermo=thermo))
                    elif isa(streams, stream_types):
                        self._streams[0] = redock(streams, stacklevel)
                    else:
                        N = len(streams)
                        n_missing(size, N) # Make sure size is not too big
                        self._streams[:N] = [redock(i, stacklevel+1) if isa(i, stream_types)
                                             else dock(Stream(i, thermo=thermo)) for i in streams]
            elif streams is not None:
                if isa(streams, str):
                    self._streams = [dock(Stream(streams, thermo=thermo))]
                elif isa(streams, stream_types):
                    self._streams = [redock(streams, stacklevel)]
                else:
                    self._streams = loaded_streams = []
                    for i in streams:
                        if isa(i, stream_types):
                            s = redock(i, stacklevel)
                        elif i is None:
                            s = self._create_missing_stream()
                        else:
                            s = Stream(i, thermo=thermo)
                            dock(s)
                        loaded_streams.append(s)
            else:
                self._initialize_missing_streams()
        
    def _create_missing_stream(self):
        return MissingStream(None, None)
        
    def _create_N_missing_streams(self, N):
        return [self._create_missing_stream() for i in range(N)]
    
    def _initialize_missing_streams(self):
        #: All input streams
        self._streams = self._create_N_missing_streams(self._size)
        
    def __add__(self, other):
        return self._streams + other
    def __radd__(self, other):
        return other + self._streams
    
    def _dock(self, stream): return stream
    def _redock(self, stream, stacklevel): return stream
    def _undock(self, stream): pass
        
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
        if stream:
            if not isinstance(stream, Stream):
                raise TypeError(
                    f"'{type(self).__name__}' object can only contain "
                    f"'Stream' objects; not '{type(stream).__name__}'"
                )
        elif not isinstance(stream, MissingStream):
            stream = self._create_missing_stream()
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
        self._undock(self._streams[int])
        self._redock(stream, stacklevel+1)
        self._streams[int] = stream
    
    def empty(self):
        for i in self._streams: self._undock(i)
        self._initialize_missing_streams()
    
    def append(self, stream):
        if self._fixed_size: 
            raise RuntimeError(f"size of '{type(self).__name__}' object is fixed")
        self._undock(stream)
        self._dock(stream)
        self._streams.append(stream)
    
    def replace(self, stream, other_stream):
        index = self.index(stream)
        self[index] = other_stream

    def index(self, stream):
        return self._streams.index(stream)

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
        self._undock(stream)
        missing_stream = self._create_missing_stream()
        self.replace(stream, missing_stream)
        
    def clear(self):
        if self._fixed_size:
            self._initialize_missing_streams()
        else:
            self._streams.clear()
    
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
    """Create an Inlets object which serves as input streams for a Unit object."""
    __slots__ = ('_sink', '_fixed_size')
    
    def __init__(self, sink, size, streams, thermo, fixed_size, stacklevel):
        self._sink = sink
        super().__init__(size, streams, thermo, fixed_size, stacklevel)
    
    @property
    def sink(self):
        return self._sink
    
    def _create_missing_stream(self):
        return MissingStream(None, self._sink)
    
    def _dock(self, stream): 
        stream._sink = self._sink
        return stream

    def _redock(self, stream, stacklevel): 
        sink = stream._sink
        if sink:
            ins = sink._ins
            if ins is not self:
                ins.remove(stream)
                stream._sink = new_sink = self._sink
                if (DOCKING_WARNINGS 
                    and sink._ID and new_sink._ID
                    and sink._ID != new_sink._ID):
                    warn(f"undocked inlet stream {stream} from unit {sink}; "
                         f"{stream} is now docked at {new_sink}", 
                         RuntimeWarning, stacklevel + 1)
        else:
            stream._sink = self._sink
        return stream
    
    def _undock(self, stream): 
        stream._sink = None
    
        
class Outlets(StreamSequence):
    """Create an Outlets object which serves as output streams for a Unit object."""
    __slots__ = ('_source',)
    
    def __init__(self, source, size, streams, thermo, fixed_size, stacklevel):
        self._source = source
        super().__init__(size, streams, thermo, fixed_size, stacklevel)
    
    @property
    def source(self):
        return self._source
    
    def _create_missing_stream(self):
        return MissingStream(self._source, None)
    
    def _dock(self, stream): 
        stream._source = self._source
        return stream

    def _redock(self, stream, stacklevel): 
        source = stream._source
        if source:
            outs = source._outs
            if outs is not self:
                outs.remove(stream)
                stream._source = new_source = self._source
                if (DOCKING_WARNINGS 
                    and source._ID and new_source._ID
                    and source._ID != new_source._ID):
                    warn(f"undocked outlet stream {stream} from unit {source}; "
                         f"{stream} is now docked at {new_source}", 
                         RuntimeWarning, stacklevel + 1)
        else:
            stream._source = self._source
        return stream
    
    def _undock(self, stream): 
        stream._source = None


# %% Sink and Source object for piping notation

class Sink:
    """
    Create a Sink object that connects a stream to a unit using piping
    notation:
    
    Parameters
    ----------
    stream : Stream
    index : int
        
    Examples
    --------
    First create a stream and a Mixer:
    
    .. code-block:: python
    
        >>> from biosteam import Stream, Mixer, settings
        >>> settings.set_thermo(['Water'])
        >>> stream = Stream('s1')
        >>> unit = Mixer('M1', outs=('out'))
    
    Sink objects are created using -pipe- notation:
        
    .. code-block:: python
    
        >>> stream-1
        <Sink: s1-1>
    
    Use pipe notation to create a sink and connect the stream:
    
    .. code-block:: python
    
        >>> stream-1-unit # The last unit is returned to continue piping; just ignore this
        <Mixer: M1>
        >>> unit.show()
        Mixer: M1
        ins...
        [0] missing stream
        [1] s1
            phase: 'l', T: 298.15 K, P: 101325 Pa
            flow:  0
        outs...
        [0] out
            phase: 'l', T: 298.15 K, P: 101325 Pa
            flow:  0
    
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


class Source:
    """
    Create a Source object that connects a stream to a unit using piping
    notation:
    
    Parameters
    ----------
    stream : Stream
    index : int
    
    Examples
    --------
    First create a stream and a Mixer:
    
    .. code-block:: python
    
        >>> from biosteam import Stream, Mixer, settings
        >>> settings.set_thermo(['Water'])
        >>> stream = Stream('s1')
        >>> unit = Mixer('M1')
    
    Source objects are created using -pipe- notation:
        
    .. code-block:: python
    
        >>> 1**stream
        <Source: 1-s1>
    
    Use -pipe- notation to create a source and connect the stream:
    
    .. code-block:: python
    
        >>> unit**0**stream # First unit is returned to continue backwards piping; just ignore this
        <Mixer: M1>
        >>> unit.show()
        Mixer: M1
        ins...
        [0] missing stream
        [1] missing stream
        outs...
        [0] s1
            phase: 'l', T: 298.15 K, P: 101325 Pa
            flow:  0
    
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

# %% System pipping

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
    def from_outlets(cls, inlets, sort=None):
        return cls([OutletPort.from_outlet(i) for i in inlets], sort)
    
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


# %% Configuration bookkeeping

Connection = namedtuple('Connection', 
                        ('source', 'source_index', 'stream', 'sink_index', 'sink'),
                        module=__name__)


# %% Stream pipping
        
def __sub__(self, index):
    if isinstance(index, int):
        return Sink(self, index)
    elif isinstance(index, Stream):
        raise TypeError("unsupported operand type(s) for -: "
                        f"'{type(self).__name__}' and '{type(index).__name__}'")
    return index.__rsub__(self)

def __rsub__(self, index):
    if isinstance(index, int):
        return Source(self, index)
    elif isinstance(index, Stream):
        raise TypeError("unsupported operand type(s) for -: "
                       f"'{type(index).__name__}' and '{type(self).__name__}'")
    return index.__sub__(self)

def get_connection(self):
    source = self._source
    source_index = source._outs.index(self) if source else None
    sink = self._sink
    sink_index = sink._ins.index(self) if sink else None
    return Connection(source, source_index, self, sink_index, sink)

Stream.get_connection = get_connection
MissingStream.__pow__ = MissingStream.__sub__ = Stream.__pow__ = Stream.__sub__ = __sub__  # Forward pipping
MissingStream.__rpow__ = MissingStream.__rsub__ = Stream.__rpow__ = Stream.__rsub__ = __rsub__ # Backward pipping    
Stream._basic_info = lambda self: (f"{type(self).__name__}: {self.ID or ''}"
                                   f"{pipe_info(self._source, self._sink)}\n")