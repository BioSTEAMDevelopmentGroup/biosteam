# -*- coding: utf-8 -*-
"""
Created on Sat Jun 15 21:01:26 2019

@author: yoelr
"""
from ..._unit import Outs, Stream, missing_stream, metaUnit, Unit, _try_method

__all__ = ('static',)

def _init_outs(self, outs):
    """Initialize output streams."""
    if outs is None:
        self._outs = Outs(self, (Stream.proxy(None) for i in self._N_outs))
    elif not outs:
        self._outs = Outs(self, (Stream.proxy('') for i in self._ins))
    elif isinstance(outs, Stream):
        self._outs = Outs(self, (outs,))
    elif isinstance(outs, str):
        self._outs = Outs(self, (Stream.proxy(outs),))
    else:
        self._outs = Outs(self, (o if isinstance(o, Stream)
                                 else Stream.proxy(o)
                                 for o in outs))
        
def _link_streams(self):
    """Link product to feed."""
    try: self._outs[0].link = self._ins[0]
    except Exception as Error:
        if missing_stream in (self._ins + self._outs):
            raise AttributeError(f'missing stream object in {repr(self)}')
   
def simulate(self):
    """Run rigourous simulation and determine all design requirements and costs."""
    self._link_streams()
    _try_method(self._run)
    _try_method(self._summary)

def _proxy_run(self):
    o = self._outs[0]
    i = self._ins[0]
    o.T = i.T
    o.P = i.P
    o._phase = i._phase
        
def _info(self, T, P, flow, fraction):
    """Information on unit."""
    self._link_streams()
    return Unit._info(self, T, P, flow, fraction)
    
class static(metaUnit):
    """Unit subclasses that implement the static metaclass expect only one input stream (the `feed`) and only one output stream (the `product`). The static metaclass makes the `feed` and `product` share data on material flow rates. Therefore `feed`.mol and `product`.mol will always be equal. Also, the `product`'s T, P, and phase will be equal to the `feed`'s upon simulation if no `_run` method is implemented."""
    def __new__(mcl, name, bases, dct):
        try:
            base, = bases
        except:
            TypeError('cannot create {mcl.__name__} subclass from more than one super class')
           
        cls = super().__new__(mcl, name, bases, dct)
        if not isinstance(base, mcl):
            cls._init_outs = _init_outs
            cls._N_ins = cls._N_outs = 1
            cls._link_streams = _link_streams
            cls.simulate = simulate
            if '_run' not in dct: cls._run = _proxy_run
        return cls
        