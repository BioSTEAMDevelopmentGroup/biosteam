# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 22:07:01 2018

@author: yoelr
"""
from ._tank import MixTank
from ._hx import HXutility

class EnzymeTreatment(MixTank):
    """Create an EnzymeTreatment unit that is cost as a MixTank with a heat exchanger."""
    _N_outs = 1
    _kwargs = {'T': 298.15}  # operating temperature (K)
    
    #: Residence time (hr)
    _tau = 1
    
    def _init(self):
        self.T = self._kwargs['T']
        self._heat_exchanger = he = HXutility(None, None, T=self.T) 
        self._heat_utilities = he._heat_utilities
        he._ins = self._ins
        he._outs = self._outs
    
    def _run(self):
        feed = self.ins[0]
        out = self.outs[0]
        out._mol[:] = self._mol_in
        out.phase = feed.phase
        out.P = feed.P
        out.T = self.T
        
    def _design(self):
        super()._design()
        self._heat_exchanger._design()
        
    def _cost(self):
        super()._cost()
        he = self._heat_exchanger
        he._cost()
        self._Cost['Heat exchanger'] = he._Cost['Heat exchanger'] 
    