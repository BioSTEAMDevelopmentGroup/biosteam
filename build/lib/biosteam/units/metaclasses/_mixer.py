# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 14:26:41 2018

@author: yoelr
"""
from ... import Stream
from ..._unit import metaUnit

__all__ = ('mixer', 'run_mixer')

class mixer(metaUnit):
    """Create a mixer Unit class which behaves like a mixer regarding mass and energy balances."""
    def __new__(mcl, name, bases, dct):
        if '_run' in dct:
            raise TypeError(f"cannot use {mcl.__name__} metaclass with an implemented '_run' method")
        dct['_run'] = run_mixer
        dct['_N_ins'] = 2
        dct['_N_outs'] = 1
        return super().__new__(mcl, name, bases, dct)

run_mixer = lambda self: Stream.sum(self.outs[0], self.ins)
    
    
