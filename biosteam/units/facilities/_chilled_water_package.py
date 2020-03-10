# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 18:28:55 2019

@author: yoelr
"""
from . import Facility
from ..decorators import cost
from ... import HeatUtility
import numpy as np

__all__ = ('ChilledWaterPackage',)

@cost('Duty', S=-14*4184000, kW=3400*0.7457, cost=1375e3, CE=551, n=0.7, BM=1.5)
class ChilledWaterPackage(Facility):
    """
    Create a chilled water package with capital cost and power based on the flow rate 
    of chilled water as in [1]_.
    
    
    Parameters
    ----------
    ID : str, optional
        Unit ID.
    
    References
    ----------
    .. [1] Humbird, D., Davis, R., Tao, L., Kinchin, C., Hsu, D., Aden, A.,
        Dudgeon, D. (2011). Process Design and Economics for Biochemical 
        Conversion of Lignocellulosic Biomass to Ethanol: Dilute-Acid 
        Pretreatment and Enzymatic Hydrolysis of Corn Stover
        (No. NREL/TP-5100-47764, 1013269). https://doi.org/10.2172/1013269
    
    """
    _N_heat_utilities = 1
    _units = {'Duty': 'kJ/hr'}
    def __init__(self, ID=''):
        chilled_water = HeatUtility.get_cooling_agent('chilled_water')
        super().__init__(ID,
                         ins='recirculated_chilled_water',
                         outs=chilled_water.to_stream(),
                         thermo=chilled_water.thermo)
        self.chilled_water_utilities = set()
        
    def _design(self):
        cwu = self.chilled_water_utilities
        if not cwu:
            for u in self.system.units:
                if u is self: continue
                for hu in u.heat_utilities:
                    if hu.ID == 'chilled_water': cwu.add(hu)
        self.design_results['Duty'] = duty = sum([i.duty for i in cwu])
        hu = self.heat_utilities[0]
        hu(duty, 330)
        used = self.ins[0]
        used.mol[0] = sum([i.flow for i in cwu])
        used.T = np.array([i.outlet_utility_stream.T for i in cwu]).mean()
        