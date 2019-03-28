# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 16:40:50 2019

@author: yoelr
"""
from biosteam import Unit, np

class ScrewFeeder(Unit):
    length = 30 #: ft
    _N_outs = 1
    _has_power_utility = True
    _has_linked_streams = True
    
    bounds = {'Volumetric flow': (400, 10000)}
    
    def _run(self):
        self.outs[0].copy_like(self.ins[0])
    
    def _cost(self):
        feed = self.ins[0]
        r = self.results
        Design = r['Design']
        volbounds = self.bounds['Volumetric flow']
        volnet = feed.volnet*35.315 # ft3/hr
        N = volnet/volbounds[-1]
        Design['N'] = N = np.ceil(N)
        volnet /= N 
        massnet = feed.massnet*0.0006124/N #lb/s
        power = N * 0.0146*massnet**0.85*self.length # hp
        power *= 0.7457 # kW
        self.power_utility(power)
        r['Cost']['Conveying belt and motor'] = N * self.CEPCI/567 * 813*volnet**0.38