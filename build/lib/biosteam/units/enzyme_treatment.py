# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 22:07:01 2018

@author: yoelr
"""
from biosteam import Unit
from biosteam.units.hx import HXutility

class EnzymeTreatment(Unit):
    _N_outs = 1
    _N_heat_utilities = 1
    kwargs = {'T': 298.15}  # operating temperature (K)
    
    def _init(self):
        self._heat_exchanger = he = HXutility(self.ID+' hx', None) 
        self.heat_utilities = he.heat_utilities
        he._ins = self._ins
        he._outs = self._outs
        
    def _setup(self):
        T = self.kwargs['T']
        self.outs[0].T = T
    
    def _run(self):
        feed = self.ins[0]
        out = self.outs[0]
        out.mol = self._mol_in
        out.phase = feed.phase
        out.P = feed.P
    
    def _operation(self):
        self._heat_exchanger._operation()
    
    def _design(self):
        self._heat_exchanger._design()
        
    def _cost(self):
        self._heat_exchanger._cost()
    