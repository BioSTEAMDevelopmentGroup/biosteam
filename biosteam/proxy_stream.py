# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 21:55:57 2019

@author: yoelr
"""
from .stream import Stream
from .mixed_stream import MixedStream

__all__ = ('ProxyStream', 'ProxyMixedStream')

#%% Proxy streams for efficiency

class ProxyStream(Stream):
    """Create a ProxyStream object that shares data with another stream. ProxyStream objects work exactly the same as Stream objects.
    
    **Parameters**
    
        **ID:** [str] ID of proxy stream
    
    """
    _data = {}
    _link = None
    
    def __new__(cls, ID):
        self = object.__new__(cls)
        cls._data[self] = [(None, None), (None, None), '']
        self.ID = ID
        return self
    
    def __init__(self, ID): pass

    def copy_like(self, stream):
        if not self._link is stream: super().copy_like(stream)
    
    @classmethod
    def asproxy(cls, stream):
        cls._data[stream] = [(None, None), (None, None), stream.ID]
        stream.__class__ = cls
        return stream
    
    @property
    def link(self):
        return self._link
    
    @link.setter
    def link(self, stream):
        if not self._link is stream:
            if hasattr(self, 'price'):
                stream.price = self.price
            self.__dict__ = stream.__dict__
            self._link = stream
        
    @property
    def _ID(self):
        return self._data[self][2]
    
    @_ID.setter
    def _ID(self, ID):
        self._data[self][2] = ID
    
    @property
    def _source(self):
        return self._data[self][0]
    
    @_source.setter
    def _source(self, source):
        self._data[self][0] = source
    
    @property
    def _sink(self):
        return self._data[self][1]
    
    @_sink.setter
    def _sink(self, sink):
        self._data[self][1] = sink
        
    def disable_phases(self, phase):
        super().disable_phases(phase)
        self.__class__ = ProxyStream
        
    def enable_phases(self):
        super().enable_phases()
        self.__class__ = ProxyMixedStream

    def _info(self, *args, **kwargs):
        linkinfo = f'\n link: {self.link}'
        if self.__dict__: return super()._info(*args, **kwargs) + linkinfo
        else: return (self._info_header() + linkinfo)
        

class ProxyMixedStream(MixedStream, ProxyStream):
    """Create a ProxyMixedStream object that shares a dictionary with a stream. ProxyMixedStream objects work exactly the same as MixedStream objects.
    
    **Parameters**
    
        **ID:** [str]
    
    """
    
    
    
    
    