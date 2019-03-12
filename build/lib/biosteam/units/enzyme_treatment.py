# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 22:07:01 2018

@author: yoelr
"""
from biosteam import Unit
from .hx import HXutility

class EnzymeTreatment(Unit):
    _N_outs = 1
    _N_heat_utilities = 1
    _has_power_utility = True
    kwargs = {'T': 298.15}  # operating temperature (K)
    
    #: Fraction of filled tank to total tank volume
    V_wf = 0.8
    
    #: Residence time (hr)
    tau = 1
    
    #: Electricity rate (kW/m3)
    electricity_rate = 0.591
    
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
    
    def _cost(self):
        he = self._heat_exchanger
        he._operation()
        he._design()
        he._cost()
        
        r = self.results
        Cost = r['Cost']
        Cost['Heat exchanger'] = he.results['Cost']['Heat exchanger']
        V = self.tau * self.ins[0].volnet / self.V_wf
        self.power_utility(self.electricity_rate*V)
        Cost['Tank'] = 12080 * V**0.525 * self.CEPCI/525.4 # Same equation as MixTank
    