# -*- coding: utf-8 -*-
"""
Created on Sat Jul 13 02:24:35 2019

@author: yoelr
"""
from ._unit import Unit

__all__ = ('ProcessSpecification',)

class ProcessSpecification(Unit):
    _N_ins = _N_outs = 1
    _power_utility = None
    def __init__(self, run, ins=None, outs=(), thermo=None):
        super().__init__(run.__name__, ins, outs, thermo)
        self.run = run
        
    @property
    def run(self):
        return self._run
    @run.setter
    def run(self, run):
        assert callable(run), "run must be a function"
        self._run = run
        
ProcessSpecification._graphics.node['shape'] = 'octagon'