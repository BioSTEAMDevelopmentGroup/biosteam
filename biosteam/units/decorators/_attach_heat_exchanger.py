# -*- coding: utf-8 -*-
"""
Created on Wed May  1 19:05:53 2019

@author: yoelr
"""
from .._hx import HXutility as _HXutility

__all__ = ('attach_heat_exchanger',)

def attach_heat_exchanger(cls):
    if cls._N_heat_utilities != 0:
        raise ValueError("number of heat utilities (_N_heat_utilities) must be 0 to attach heat exchanger")
    
    class NewUnit(cls):
        line = cls.line
            
        def _init(self):
            super()._init()
            self._heat_exchanger = he = _HXutility(None, None) 
            self._heat_utilities = he._heat_utilities
            he._ins = self._ins
            he._outs = self._outs
        
        def _design(self):
            super()._design()
            self._heat_exchanger._design()
            
        def _cost(self):
            super()._cost()
            he = self._heat_exchanger
            he._cost()
            self._results['Cost']['Heat exchanger'] = he._results['Cost']['Heat exchanger'] 
        
        _init.__doc__ = cls.__doc__
        _design.__doc__ = cls.__doc__
        _cost.__doc__ = cls.__doc__
    
    NewUnit.__name__ = cls.__name__
    return NewUnit
 