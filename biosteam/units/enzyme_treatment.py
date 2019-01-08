# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 22:07:01 2018

@author: yoelr
"""
from biosteam import Unit
from biosteam.units.hx import HX

class EnzymeTreatment(Unit):
    _N_outs = 1
    _N_heat_util = 1
    kwargs = {'T': 298.15}  # operating temperature (K)
    
    def setup(self):
        T = self.kwargs['T']
        self.outs[0].T = T
        
    def run(self):
        feed = self.ins[0]
        out = self.outs[0]
        out.mol = self._mol_in
        out.phase = feed.phase
        out.P = feed.P
    _simple_run = run
    
    operation = HX.operation