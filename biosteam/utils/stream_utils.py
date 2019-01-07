# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 16:47:33 2018

This module includes classes and functions relating Stream objects.

@author: Yoel Cortes-Pena
"""

__all__ = ('get_streams', 'MissingStream', 'Ins', 'Outs', 'Sink', 'Source')

# For later use in this module
def missing_method(self, *args, **kwargs):
    raise TypeError(f'Method not supported for {type(self).__name__} objects.')


# %% Functions relating streams

def get_streams(units):
    """Return tuple of Stream objects from an iterable of Unit objects."""
    streams = []
    for u in units:
        streams += u._ins + u._outs
    return streams


# %% Dummy Stream object

class MissingStream:
    """Create a MissingStream object that acts as a dummy in Ins and Outs objects until replaced by an actual Stream object."""
    ID = 'Missing Stream'
    
    def __getattr__(self, key):
        raise TypeError(f'{self.ID} object')

    def _info(self, **show_units):
        return self.ID

    def __repr__(self):
        return f'<{type(self).__name__}>'

    def show(self, **show_units):
        print(self._info(**show_units))


# %% List objects for input and output streams

class Ins(list):
    """Create a Ins object which serves as a list of input streams for a Unit object."""
    __slots__ = ['sink']

    def __init__(self, sink, streams=None):
        self.sink = sink #: Unit where inputs are attached
        if streams:
            if not hasattr(streams, '__iter__'):
                streams = (streams,)
            ID = sink.ID
            self._clear_sink(ID)
            super().__init__(streams)
            self._fix_sink(ID)
        
    def _clear_sink(self, ID):
        """Remove sink from streams."""
        for s in self:
            if s._sink[0] is ID:
                s._sink = (None, None)
    
    def _fix_sink(self, ID):
        """Set sink for all streams."""
        i = 0
        for s in self:
            s._sink = (ID, i)
            i += 1
    
    def __setitem__(self, index, stream):
        ID = self.sink.ID
        if isinstance(stream, MissingStream):
            pass
        elif isinstance(index, int):
            s_old = self[index]
            if s_old._sink[0] is ID:
                s_old._sink = (None, None)
            stream._sink = (ID, index)
        elif isinstance(index, slice):
            self._clear_sink(ID)
            self._set_sink(ID)
        else:
            raise TypeError(f'Only intergers and slices are valid indeces for {type(self).__name__} objects')
        super().__setitem__(index, stream)
           
    def clear(self):
        self._clear_sink(self.sink.ID)
        super().clear()
    
    def extend(self, iterable):
        index = len(self)
        ID = self.sink.ID
        # Add sink to new streams
        for s in iterable:
            s._sink = (ID, index)
            index += 1
        super().extend(iterable)
        
    def append(self, stream):
        index = len(self)
        ID = self.sink.ID
        # Add sink to new stream
        stream._sink = (ID, index)
        super().append(stream)
    
    pop = insert = remove = reverse = sort = missing_method
    

class Outs(list):
    """Create a Outs object which serves as a list of output streams for a Unit object."""
    __slots__ = ['source']
    
    def __init__(self, unit, streams=None):
        self.source = unit #: Unit where Outputs is attached
        if streams:
            if not hasattr(streams, '__iter__'):
                streams = (streams,)
            ID = unit.ID
            self._clear_source(ID)
            super().__init__(streams)
            self._fix_source(ID)
            
    def _clear_source(self, ID):
        """Remove source from streams."""
        for s in self:
            if s._source[0] is ID:
                s._source = (None, None)
    
    def _fix_source(self, ID):
        """Set source for all streams."""
        i = 0
        for s in self:
            s._source = (ID, i)
            i += 1
            
    def __setitem__(self, index, stream):
        ID = self.source.ID
        if isinstance(stream, MissingStream):
            pass
        elif isinstance(index, int):
            s_old = self[index]
            if s_old._source[0] is ID:
                s_old._source = (None, None)
            stream._source = (ID, index)
        elif isinstance(index, slice):
            self._clear_source(ID)
            self._fix_source(ID)
        else:
            raise TypeError(f'Only intergers and slices are valid indeces for {type(self).__name__} objects')
        super().__setitem__(index, stream)
           
    def clear(self):
        self._clear_source(self.source.ID)
        super().clear()
    
    def extend(self, iterable):
        index = len(self)
        ID = self.source.ID
        # Add source to new streams
        for s in iterable:
            s._source = (ID, index)
            index += 1
        super().extend(iterable)
        
    def append(self, stream):
        index = len(self)
        ID = self.source.ID
        # Add sink to new stream
        stream._source = (ID, index)
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




