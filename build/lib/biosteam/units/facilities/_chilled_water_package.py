# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 18:28:55 2019

@author: yoelr
"""
from . import Facility
from ...utils.stream_link_options import static
from ..decorators import cost
from ... import HeatUtility
import numpy as np

__all__ = ('ChilledWaterPackage',)

@static
@cost('Duty', S=-14*4184000, kW=3400*0.7457, cost=1375e3, CE=551, n=0.7, BM=1.5)
class ChilledWaterPackage(Facility):
    """Create a chilled water package that is cost based on flow rate of chilled water.
    
    Parameters
    ----------
    Duty : float
        Chilled water duty (kJ/hr)
    
    """
    _N_heat_utilities = 1
    _units = {'Duty': 'kJ/hr'}
    def __init__(self, ID=''):
        thermo = HeatUtility.cooling_agents['Chilled water'].thermo
        super().__init__(ID, 'return_chilled_water', 'chilled_water', thermo=thermo)
        self.chilled_water_utilities = set()
        
    def _design(self):
        cwu = self.chilled_water_utilities
        if not cwu:
            for u in self.system.units:
                if u is self: continue
                for hu in u.heat_utilities:
                    if hu.ID == 'Chilled water': cwu.add(hu)
        self.design_results['Duty'] = duty = sum([i.duty for i in cwu])
        hu = self.heat_utilities[0]
        hu(duty, 330)
        used = self._ins[0]
        used.mol[0] = sum([i.flow for i in cwu])
        used.T = np.array([i.outlet_utility_stream.T for i in cwu]).mean()
        self.outs[0].T = hu.cooling_agents['Chilled water'].T
        