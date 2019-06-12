# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 22:07:01 2018

@author: yoelr
"""
from ._tank import MixTank
from .decorators import attach_heat_exchanger

@attach_heat_exchanger
class EnzymeTreatment(MixTank):
    """Create an EnzymeTreatment unit that is cost as a MixTank with a heat exchanger."""
    _N_outs = 1
    _kwargs = {'T': 298.15}  # operating temperature (K)
    
    #: Residence time (hr)
    _tau = 1
    
    def _setup(self):
        self.outs[0].T = self._kwargs['T']
    
    def _run(self):
        feed = self.ins[0]
        out = self.outs[0]
        out._mol[:] = self._mol_in
        out.phase = feed.phase
        out.P = feed.P
    