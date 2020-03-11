# -*- coding: utf-8 -*-
"""
Created on Sat Jul 13 02:24:35 2019

@author: yoelr
"""
from ._unit import Unit
from .utils import colors

__all__ = ('ProcessSpecification',)

class ProcessSpecification(Unit):
    _N_ins = _N_outs = 1
    power_utility = None
    results = None
    heat_utilities = ()
    
    def __init__(self, run, ins=None, outs=(), thermo=None):
        self._load_thermo(thermo)
        self._init_ins(ins)
        self._init_outs(outs)
        self._assert_compatible_property_package()
        self._register(run.__name__)
        self.run = run
        
    @property
    def run(self):
        return self._run
    @run.setter
    def run(self, run):
        assert callable(run), "run must be a function"
        self._run = run
        
orange = colors.orange_tint.tint(50)
orange_tint = orange.tint(75)
node = ProcessSpecification._graphics.node
node['fillcolor'] = orange_tint.HEX + ':' + orange.HEX
node['shape'] = 'note'
node['margin'] = '0.2'

