# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 21:55:57 2019

@author: yoelr
"""
from ._stream import Stream
from ._mixed_stream import MixedStream

__all__ = ('ProxyStream', 'ProxyMixedStream')

#%% Proxy streams for efficiency

class ProxyStream(Stream):
    """Create a ProxyStream object that shares data with another stream. ProxyStream objects work exactly the same as Stream objects.
    
    **Parameters**
    
        **ID:** [str] ID of proxy stream
        
        **link:** [Stream] proxy stream will share data from this link.
    
    """
    _data = {}
    _link = None
    __slots__ = ()
    
    def __init__(self, ID='', link=None):
        ProxyStream._data[self] = [None, None, None]
        self.ID = ID
        self.link = link

    def copylike(self, stream):
        if not self._link is stream: super().copylike(stream)
    copylike.__doc__ = Stream.copylike.__doc__
    
    @classmethod
    def asproxy(cls, stream):
        """Cast the stream as a proxy object."""
        cls._data[stream] = [stream._source, stream._sink, stream.ID, stream.price]
        stream.__class__ = cls
        return stream
    
    @property
    def link(self):
        """Stream that the proxy stream is sharing data with."""
        return self._link
    
    @link.setter
    def link(self, stream):
        if not self._link is stream:
            self.__dict__ = stream.__dict__
            stream._link = stream
        elif stream and not isinstance(stream, Stream):
            raise TypeError(f"link must be a Stream object or None, not a '{type(stream).__name__}' object.")
       
    @property
    def price(self):
        return self._data[self][3]
    @price.setter
    def price(self, price): self._data[self][3] = price
        
    @property
    def _ID(self):
        return self._data[self][2]
    @_ID.setter
    def _ID(self, ID):
        self._data[self][2] = ID
    
    @property
    def _source(self): return self._data[self][0]
    @_source.setter
    def _source(self, source): self._data[self][0] = source
    
    @property
    def _sink(self):return self._data[self][1]
    @_sink.setter
    def _sink(self, sink): self._data[self][1] = sink
        
    @property
    def P(self):
        if hasattr(self, '_P'): return self._P
        else: return self.__dict__['P']
    @P.setter
    def P(self, P): self._P = P
    
    def disable_phases(self, phase):
        super().disable_phases(phase)
        self.__class__ = ProxyStream
    disable_phases.__doc__ = Stream.disable_phases.__doc__    
    
    def enable_phases(self):
        super().enable_phases()
        self.__class__ = ProxyMixedStream
    enable_phases.__doc__ = Stream.enable_phases.__doc__

    def _info(self, *args, **kwargs):
        linkinfo = f'\n link: {self.link}'
        if self.__dict__: return super()._info(*args, **kwargs) + linkinfo
        else: return (self._info_header() + linkinfo)
        

class ProxyMixedStream(MixedStream, ProxyStream):
    """Create a ProxyMixedStream object that shares a dictionary with a stream. ProxyMixedStream objects work exactly the same as MixedStream objects.
    
    **Parameters**
    
        **ID:** [str]
    
    """
    
    
    
    
    