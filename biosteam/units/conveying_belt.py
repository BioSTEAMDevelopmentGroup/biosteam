# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 11:10:49 2019

@author: yoelr
"""
from biosteam import Unit, np

class ConveyingBelt(Unit):
    length = 50 #: ft
    height = 30 #: ft
    _N_outs = 1
    _has_power_utility = True
    
    bounds = {'Volumetric flow': (120, 500)}
    
    def _run(self):
        self.outs[0].copy_like(self.ins[0])
    
    def _cost(self):
        feed = self.ins[0]
        r = self.results
        Design = r['Design']
        volbounds = self.bounds['Volumetric flow']
        volnet = feed.volnet*35.315 # ft3/hr
        N_belts = volnet/volbounds[-1]
        Design['N_belts'] = N_belts = np.ceil(N_belts)
        volnet /= N_belts 
        massnet = feed.massnet*0.0006124/N_belts #lb/s
        power = N_belts * 0.00058*massnet**0.82*self.length + self.height*0.00182*massnet# hp
        power *= 0.7457 # kW
        self.power_utility(power)
        r['Cost']['Conveying belt and motor'] = N_belts * self.CEPCI/567 * 813*volnet**0.38
        

