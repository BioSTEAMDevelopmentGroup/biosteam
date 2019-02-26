# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 16:47:33 2018

This module includes classes and functions relating Stream objects.

@author: Yoel Cortes-Pena
"""

__all__ = ('MissingStream', 'Ins', 'Outs',
           'Sink', 'Source', 'missing_stream')

# For later use in this module
def missing_method(self, *args, **kwargs):
    raise TypeError(f'Method not supported for {type(self).__name__} object.')


# %% Dummy Stream object

class MissingStream:
    """Create a MissingStream object that acts as a dummy in Ins and Outs objects until replaced by an actual Stream object."""
    ID = 'Missing Stream'
    _missing_stream = None
    __slots__ = ('_sink', '_source')
    
    def __new__(cls):
        s = cls._missing_stream
        if not s:
            cls._missing_stream = s = super().__new__(cls)
        return s
    
    def __getattr__(self, key):
        raise TypeError(f'{self.ID} object')

    def _info(self, **show_units):
        return self.ID

    def __repr__(self):
        return f'<{type(self).__name__}>'

    def show(self, **show_units):
        print(self._info(**show_units))

missing_stream = MissingStream()

# %% List objects for input and output streams

class Ins(list):
    """Create a Ins object which serves as a list of input streams for a Unit object."""
    __slots__ = ('sink',)

    def __init__(self, sink, streams=()):
        if streams is self:
            streams = tuple(self)
        elif not hasattr(streams, '__iter__'):
            streams = (streams,)
        self.sink = sink #: Unit where inputs are attached
        self._clear_sink(sink)
        super().__init__(streams)
        self._fix_sink(sink)
        
    def _clear_sink(self, sink):
        """Remove sink from streams."""
        for s in self:
            if s._sink[0] is sink:
                s._sink = (None, None)
    
    def _fix_sink(self, sink):
        """Set sink for all streams."""
        i = 0
        for s in self:
            s._sink = (sink, i)
            i += 1
    
    def __setitem__(self, index, stream):
        sink = self.sink
        if stream is missing_stream:
            pass
        elif isinstance(index, int):
            s_old = self[index]
            if s_old._sink[0] is sink:
                s_old._sink = (None, None)
            stream._sink = (sink, index)
        elif isinstance(index, slice):
            self._clear_sink(sink)
            self._fix_sink(sink)
        else:
            raise TypeError(f'Only intergers and slices are valid indeces for {type(self).__name__} objects')
        super().__setitem__(index, stream)
           
    def clear(self):
        self._clear_sink(self.sink)
        super().clear()
    
    def extend(self, iterable):
        index = len(self)
        sink = self.sink
        # Add sink to new streams
        for s in iterable:
            s._sink = (sink, index)
            index += 1
        super().extend(iterable)
        
    def append(self, stream):
        index = len(self)
        sink = self.sink
        # Add sink to new stream
        stream._sink = (sink, index)
        super().append(stream)
    
    pop = insert = remove = reverse = sort = missing_method
    

class Outs(list):
    """Create a Outs object which serves as a list of output streams for a Unit object."""
    __slots__ = ('source',)
    
    def __init__(self, source, streams=()):
        if streams is self:
            streams = tuple(self)
        elif not hasattr(streams, '__iter__'):
            streams = (streams,)
        self.source = source #: Unit where Outputs is attached
        self._clear_source(source)
        super().__init__(streams)
        self._fix_source(source)
            
    def _clear_source(self, source):
        """Remove source from streams."""
        for s in self:
            if s._source[0] is source:
                s._source = (None, None)
    
    def _fix_source(self, source):
        """Set source for all streams."""
        i = 0
        for s in self:
            s._source = (source, i)
            i += 1
            
    def __setitem__(self, index, stream):
        source = self.source
        if stream is missing_stream:
            pass
        elif isinstance(index, int):
            s_old = self[index]
            if s_old._source[0] is source:
                s_old._source = (None, None)
            stream._source = (source, index)
        elif isinstance(index, slice):
            self._clear_source(source)
            self._fix_source(source)
        else:
            raise TypeError(f'Only intergers and slices are valid indeces for {type(self).__name__} objects')
        super().__setitem__(index, stream)
           
    def clear(self):
        self._clear_source(self.source)
        super().clear()
    
    def extend(self, iterable):
        index = len(self)
        source = self.source
        # Add source to new streams
        for s in iterable:
            s._source = (source, index)
            index += 1
        super().extend(iterable)
        
    def append(self, stream):
        index = len(self)
        source = self.source
        # Add sink to new stream
        stream._source = (source, index)
        super().append(stream)

    pop = insert = remove = reverse = sort = missing_method

# %% Sink and Source object for piping notation

class Sink:
    """Create a Sink object that connects a stream to a unit using piping notation:
    
    **Parameters**
    
        stream: [Stream]
        
        index: [int]
        
    **Example**
    
    First create a stream and a Mixer:
    
    .. code-block:: python
    
        >>> stream = Stream('s1')
        >>> unit = Mixer('M1')
    
    Sink objects are created using piping notation:
        
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
            flows:  0
        outs...
        [0] d27
            phase: 'l', T: 298.15 K, P: 101325 Pa
            flows:  0
    
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




