# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 16:47:33 2018

This module includes classes and functions relating Stream objects.

@author: Yoel Cortes-Pena
"""
from thermosteam import Stream, MultiStream

__all__ = ('MissingStream', 'Ins', 'Outs', 'Sink', 'Source',
           'as_stream', 'as_upstream', 'as_downstream')

isa = isinstance

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
    if isa(stream, Stream):
        return stream
    elif isa(stream, str):
        return Stream(stream)
    elif stream is None:
        return MissingStream()
    
def as_upstream(stream, sink):
    stream = as_stream(stream)
    stream._sink = sink
    return stream
        
def as_downstream(stream, source):
    stream = as_stream(stream)
    stream._source = source
    return stream


# %% Dummy Stream object

class MissingStream:
    """Create a MissingStream object that acts as a dummy in Ins and Outs objects until replaced by an actual Stream object."""
    __slots__ = ('_source', '_sink')
    
    def __bool__(self):
        return False
    
    def __getattr__(self, key):
        raise AttributeError(str(self))

    def __repr__(self):
        return f'<MissingStream>'
    
    def __str__(self):
        return 'missing stream'


# %% Utilities

def n_missing(ub, N):
    assert ub >= N, f"size of streams exceeds {ub}"
    return ub - N

MissingStream = MissingStream()

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
        
    def __init__(self, size, streams, thermo, fixed_size):
        self._size = size
        self._fixed_size = fixed_size
        dock = self._dock
        redock = self._redock
        if streams == ():
                self._streams = [dock(Stream(thermo=thermo)) for i in range(size)]
        else:
            if fixed_size:
                self._streams = [MissingStream] * size #: All input streams
                if streams:
                    if isa(streams, str):
                        self._streams[0] = dock(Stream(streams, thermo=thermo))
                    elif isa(streams, (Stream, MultiStream)):
                        self._streams[0] = redock(streams)
                    else:
                        N = len(streams)
                        n_missing(size, N) # Assert size is not too big
                        self._streams[:N] = [redock(i) if isa(i, Stream)
                                             else dock(Stream(i, thermo=thermo)) for i in streams]
            else:
                if streams:
                    if isa(streams, str):
                        self._streams = [dock(Stream(streams, thermo=thermo))]
                    elif isa(streams, (Stream, MultiStream)):
                        self._streams = [redock(streams)]
                    else:
                        self._streams = [redock(i) if isa(i, Stream)
                                         else dock(Stream(i, thermo=thermo)) for i in streams]
                else:
                    self._streams = size * [MissingStream]
        
    def __add__(self, other):
        return self._streams + other
    def __radd__(self, other):
        return other + self._streams
    
    def _dock(self, stream): return stream
    def _redock(self, stream): return stream
    def _undock(self, stream): pass
        
    def _set_streams(self, slice, streams):
        all_streams = self._streams
        for stream in all_streams[slice]: self._undock(stream)
        all_streams[slice] = streams
        for stream in all_streams:
            if stream: self._redock(stream)
        if self._fixed_size:
            size = self._size
            N_streams = len(all_streams)
            if N_streams < size:
                N_missing = n_missing(size, N_streams)
                if N_missing:
                    all_streams[N_streams: size] = (MissingStream,) * N_missing   
            
    @property
    def size(self):
        return self._streams.__len__()
    
    def __len__(self):
        return self._streams.__len__()
    
    def _set_stream(self, int, stream):
        self._undock(self._streams[int])
        if stream: self._redock(stream)
        self._streams[int] = stream

    def index(self, stream):
        return self._streams.index(stream)

    def pop(self, index):
        streams = self._streams
        if self._fixed_size:
            stream = streams[index]
            streams[index] = MissingStream
        else:
            stream = streams.pop(index)
        return stream

    def remove(self, stream):
        streams = self._streams
        self._undock(stream)
        if self._fixed_size:
            index = streams.index(stream)
            streams[index] = MissingStream
        else:
            streams.remove(stream)
        
    def clear(self):
        if self._fixed_size:
            self._streams = [MissingStream] * self.size
        else:
            self._streams.clear()
    
    def __iter__(self):
        yield from self._streams
    
    def __getitem__(self, index):
        return self._streams[index]
            
    def __setitem__(self, index, item):
        if isa(index, int):
            if item:
                assert isa(item, Stream), (
                    f"'{type(self).__name__}' object can only contain "
                    f"'Stream' objects; not '{type(item).__name__}'")
            else:
                item = MissingStream
            self._set_stream(index, item)
        elif isa(index, slice):
            streams = []
            for i in item:
                if i:
                    assert isa(i, Stream), (
                        f"'{type(self).__name__}' object can only contain "
                        f"'Stream' objects; not '{type(i).__name__}'")
                    streams.append(i)
                else:
                    streams.append(MissingStream)
            self._set_streams(index, item)
        else:
            raise TypeError(f"Only intergers and slices are valid indices for '{type(self).__name__}' objects")
    
    def __repr__(self):
        return repr(self._streams)


class Ins(StreamSequence):
    """Create an Ins object which serves as input streams for a Unit object."""
    __slots__ = ('_sink', '_fixed_size')
    
    def __init__(self, sink, size, streams, thermo, fixed_size=True):
        self._sink = sink
        super().__init__(size, streams, thermo, fixed_size)
    
    @property
    def sink(self):
        return self._sink
    
    def _dock(self, stream): 
        stream._sink = self._sink
        return stream

    def _redock(self, stream): 
        sink = stream._sink
        if sink:
            ins = sink._ins
            if ins is not self:
                ins.remove(stream)
                stream._sink = self._sink
        else:
            stream._sink = self._sink
        return stream
    
    def _undock(self, stream): 
        stream._sink = None
    
        
class Outs(StreamSequence):
    """Create an Outs object which serves as output streams for a Unit object."""
    __slots__ = ('_source',)
    
    def __init__(self, source, size, streams, thermo, fixed_size=True):
        self._source = source
        super().__init__(size, streams, thermo, fixed_size)
    
    @property
    def source(self):
        return self._source
    
    def _dock(self, stream): 
        stream._source = self._source
        return stream

    def _redock(self, stream): 
        source = stream._source
        if source:
            outs = source._outs
            if outs is not self:
                # Remove from source
                outs.remove(stream)
                stream._source = self._source
        else:
            stream._source = self._source
        return stream
    
    def _undock(self, stream): 
        stream._source = None



# %% Sink and Source object for piping notation

class Sink:
    """
    Create a Sink object that connects a stream to a unit using piping notation:
    
    Parameters
    ----------
    stream : Stream
    index : int
        
    Examples
    --------
    First create a stream and a Mixer:
    
    .. code-block:: python
    
        >>> stream = Stream('s1')
        >>> unit = Mixer('M1')
    
    Sink objects are created using -pipe- notation:
        
    .. code-block:: python
    
        >>> stream-1
        <Sink: s1-1>
    
    Use pipe notation to create a sink and connect the stream:
    
    .. code-block:: python
    
        >>> stream-1-unit
        >>> M1.show()
        
        Mixer: M1
        ins...
        [0] Missing stream
        [1] s1
            phase: 'l', T: 298.15 K, P: 101325 Pa
            flow:  0
        outs...
        [0] d27
            phase: 'l', T: 298.15 K, P: 101325 Pa
            flow:  0
    
    """
    __slots__ = ('stream', 'index')
    def __init__(self, stream, index):
        self.stream = stream or MissingStream
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
    Create a Source object that connects a stream to a unit using piping notation:
    
    Parameters
    ----------
    stream : Stream
    index : int
    
    Examples
    --------
    First create a stream and a Mixer:
    
    .. code-block:: python
    
        >>> stream = Stream('s1')
        >>> unit = Mixer('M1')
    
    Source objects are created using -pipe- notation:
        
    .. code-block:: python
    
        >>> 1**stream
        <Source: 1-s1>
    
    Use -pipe- notation to create a source and connect the stream:
    
    .. code-block:: python
    
        >>> unit**0**stream
        >>> M1.show()
        
        Mixer: M1
        ins...
        [0] Missing stream
        [1] Missing stream
        outs...
        [0] s1
            phase: 'l', T: 298.15 K, P: 101325 Pa
            flow:  0
    
    """
    __slots__ = ('stream', 'index')
    def __init__(self, stream, index):
        self.stream = stream or MissingStream
        self.index = index

    # Forward pipping
    def __rsub__(self, unit):
        unit.outs[self.index] = self.stream
        return unit
    
    # Backward pipping
    __rpow__ = __rsub__
    
    def __repr__(self):
        return '<' + type(self).__name__ + ': ' + str(self.index) + '-' + self.stream.ID + '>'

# %% Pipping
        
def __sub__(self, index):
    if isinstance(index, int):
        return Sink(self, index)
    elif isinstance(index, Stream):
        raise TypeError("unsupported operand type(s) for -: "
                        f"'{type(self)}' and '{type(index)}'")
    return index.__rsub__(self)

def __rsub__(self, index):
    if isinstance(index, int):
        return Source(self, index)
    elif isinstance(index, Stream):
        raise TypeError("unsupported operand type(s) for -: "
                        "'{type(self)}' and '{type(index)}'")
    return index.__sub__(self)

Stream.__pow__ = Stream.__sub__ = __sub__  # Forward pipping
Stream.__rpow__ = Stream.__rsub__ = __rsub__ # Backward pipping    
Stream.sink = property(lambda self: self._sink)
Stream.source = property(lambda self: self._source)
Stream._basic_info = lambda self: (f"{type(self).__name__}: {self.ID or ''}"
                                   f"{pipe_info(self._source, self._sink)}\n")
