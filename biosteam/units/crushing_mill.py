# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 22:14:01 2018

@author: yoelr
"""
from .. import Unit
from . import Splitter


class CrushingMill(Unit):
    _kwargs = Splitter._kwargs
    _run = Splitter._run
    _has_power_utility = True
    
    #: Original Price (USD)
    C_0 = 1.5e6
    
    #: Original flow rate (kg/hr)
    V_0 = 335e3 
    
    #: Scaling exponent
    exp = 0.60 
    
    #: Electricity rate (kW/(kg/hr))
    electricity_rate = 0.006
    
    #: Original CEPCI 
    CEPCI_0 = 541.7
    
    def _cost(self):
        # Size factor
        massflow = self.ins[0].massnet
        S = massflow/self.V_0
        self._results['Cost']['Crushing mill'] = (self.CEPCI/self.CEPCI_0) * self.C_0*S**self.exp
        self._power_utility(massflow*self.electricity_rate)
