# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import biosteam as bst
cost = bst.decorators.cost

__all__ = ('CoolingTower',)


# %%

@cost('Flow rate', 'Cooling water pump',
      S=609624, kW=1021, cost=283671, CE=551, n=0.8, BM=3.1)
@cost('Flow rate', 'Cooling tower',
      S=557183, kW=1598, cost=1375e3, CE=551, n=0.7, BM=1.5)
class CoolingTower(bst.Facility):
    """
    Create a cooling tower with capital cost and power based on the flow rate 
    of cooling water as in [1]_.
    
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
    ticket_name = 'CT'
    network_priority = 0
    _units = {'Flow rate': 'kmol/hr'}
    _N_ins = 3
    _N_outs = 3
    evaporation = 0.01
    blowdown = 0.001
    def __init__(self, ID='', agent=None):
        self.agent = cooling_water = agent or bst.settings.get_cooling_agent('cooling_water')
        self.makeup_water = makeup_water = cooling_water.to_stream('cooling_tower_makeup_water')
        return_cooling_water = cooling_water.to_stream()
        cooling_water = return_cooling_water.flow_proxy()
        cooling_tower_chemicals = return_cooling_water.copy('cooling_tower_chemicals')
        cooling_tower_chemicals.price = 3.
        blowdown = makeup_water.copy('cooling_tower_blowdown')
        evaporation = makeup_water.copy('cooling_tower_evaporation')
        super().__init__(ID, (return_cooling_water, makeup_water, cooling_tower_chemicals),
                         (cooling_water, blowdown, evaporation), thermo=cooling_water.thermo)
        self.cooling_water_utilities = []
    
    @property
    def return_cooling_water(self):
        return self.ins[0]
    
    @property
    def cooling_water(self):
        return self.outs[0]
    
    @property
    def blowdown_water(self):
        return self.outs[1]
    
    @property
    def evaporation_water(self):
        return self.outs[2]
    
    def _run(self): pass
        
    def _load_utility_agents(self):
        cwu = self.cooling_water_utilities
        ID = self.agent.ID
        cwu.clear()
        for u in self.other_units:
            for hu in u.heat_utilities:
                agent = hu.agent
                if agent and agent.ID == ID: cwu.append(hu)
        
    def _design(self):
        self._load_utility_agents()
        cwu = self.cooling_water_utilities
        return_cooling_water, makeup_water, cooling_tower_chemicals = self._ins
        hu = self.create_heat_utility()
        self._load_utility_agents()
        hu.mix_from(cwu)            
        return_cooling_water.imol['7732-18-5'] = self.design_results['Flow rate'] = hu.flow
        cooling_tower_chemicals.imass['Water'] = 2 * return_cooling_water.F_mol / 4.4e+05
        self.return_cooling_water.T = hu.inlet_utility_stream.T
        self.evaporation_water.mol[0] = evaporation = hu.flow  * self.evaporation
        self.blowdown_water.mol[0] = blowdown = hu.flow  * self.blowdown
        self.makeup_water.mol[0] = evaporation + blowdown
        self.evaporation_water.T = self.cooling_water.T = self.blowdown_water.T = hu.outlet_utility_stream.T
        self.evaporation_water.phase = 'g'
        hu.reverse()

