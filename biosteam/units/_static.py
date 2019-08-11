# -*- coding: utf-8 -*-
"""
Created on Sat Jun 15 21:01:26 2019

@author: yoelr
"""
from .._unit import Outs, Stream, MissingStream, Unit

__all__ = ('Static',)

class Static(Unit):
    """A Static Unit object expects only one input stream (the `feed`) and only one output stream (the `product`). The static metaclass makes the `feed` and `product` share data on material flow rates. Therefore `feed`.mol and `product`.mol will always be equal. For Static subclasses, the `product`'s T, P, and phase will be equal to the `feed`'s upon simulation if no `_run` method is implemented."""
    _N_ins = _N_outs = 1
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
            
    def link_streams(self):
        """Link product to feed."""
        try: self._outs[0].link = self._ins[0]
        except Exception as Error:
            if MissingStream in (self._ins + self._outs):
                raise AttributeError(f'missing stream object in {repr(self)}')
       
    def simulate(self):
        """Run rigourous simulation and determine all design requirements and costs."""
        self.link_streams()
        super().simulate()
    
    def _run(self):
        o = self._outs[0]
        i = self._ins[0]
        o.T = i.T
        o.P = i.P
        o._phase = i._phase
            
    def _info(self, T, P, flow, fraction):
        """Information on unit."""
        self.link_streams()
        return super()._info(T, P, flow, fraction)
    
    
    