# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
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
    ticket_name = 'CWP'
    network_priority = 0
    _N_heat_utilities = 2
    _units = {'Duty': 'kJ/hr'}
    def __init__(self, ID='', agent=None):
        self.agent = chilled_water = agent or HeatUtility.get_cooling_agent('chilled_water')
        super().__init__(ID,
                         ins='recirculated_chilled_water',
                         outs=chilled_water.to_stream(),
                         thermo=chilled_water.thermo)
        
    def _load_chilled_water_utilities(self):
        self.chilled_water_utilities = cwu = set()
        ID = self.agent.ID
        for u in self.other_units:
            if u is self: continue
            for hu in u.heat_utilities:
                agent = hu.agent 
                if agent and agent.ID == ID: cwu.add(hu)
        
    def _design(self):
        self._load_chilled_water_utilities()
        cwu = self.chilled_water_utilities
        self.design_results['Duty'] = duty = sum([i.duty for i in cwu])        
        hu_cooling, hu_chilled = self.heat_utilities
        hu_chilled.mix_from(cwu)
        hu_chilled.reverse()
        hu_cooling(duty, 330)
        used = self.ins[0]
        used.mol[0] = sum([i.flow for i in cwu])
        used.T = np.array([i.outlet_utility_stream.T for i in cwu]).mean()
        