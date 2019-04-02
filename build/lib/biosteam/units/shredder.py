# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 10:57:12 2019

@author: yoelr
"""
from biosteam import Unit

class Shredder(Unit):
    _has_power_utility = True
    _has_proxystream = True
    _N_outs = 1
    
    #: Original Price (USD)
    C_0 = 2.5e6
    
    #: Original flow rate (kg/hr)
    V_0 = 500e3 
    
    #: Scaling exponent
    exp = 0.60 
    
    #: Electricity rate (kW/(kg/hr))
    electricity_rate = 0.006
    
    #: Original CEPCI 
    CEPCI_0 = 567.3
    
    def _run(self):
        self.outs[0].copylike(self.ins[0])
    
    def _cost(self):
        # Size factor
        massflow = self.ins[0].massnet
        S = massflow/self.V_0
        self._results['Cost']['Shredder'] = (self.CEPCI/self.CEPCI_0) * self.C_0*S**self.exp
        self._power_utility(massflow*self.electricity_rate)