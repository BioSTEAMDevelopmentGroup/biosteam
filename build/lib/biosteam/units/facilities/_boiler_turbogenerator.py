# -*- coding: utf-8 -*-
"""
Created on Sun Feb 17 16:36:47 2019

@author: yoelr
"""
from ... import HeatUtility, Unit
from . import Facility
from ..decorators import cost

__all__ = ('BoilerTurbogenerator',)

#: TODO add reference of NREL

@cost('Work', 'Turbogenerator',
      CE=551, S=42200, kW=0, cost=9500e3, n=0.60, BM=1.8)
@cost('Flow rate', 'Hot process water softener system', 
      CE=551, cost=78e3, S=235803, n=0.6, BM=1.8)
@cost('Flow rate', 'Amine addition pkg', 
      CE=551, cost=40e3, S=235803, n=0.6, BM=1.8)
@cost('Flow rate', 'Deaerator',
      CE=551, cost=305e3, S=235803, n=0.6, BM=3.0)
@cost('Flow rate', 'Boiler',
      CE=551, cost=28550e3, kW=1000, S=238686, n=0.6, BM=1.8)
class BoilerTurbogenerator(Facility):
    """
    Create a BoilerTurbogenerator object that will calculate electricity
    generation from burning the feed. It also takes into account how much
    steam is being produced, and the required cooling utility of the turbo
    generator. No emissions or mass balances are taken into account. All 
    capital cost correlations are based on [1]_.
    
    Parameters
    ----------
    ins : stream sequence
        [0] Liquid/solid feed that will be burned.
        
        [1] Gas feed that will be burned.
        
        [2] Make-up water. 
    outs : stream sequence
        [0] Total emissions produced.
        
        [1] Blowdown water.
    boiler_efficiency : float
        Fraction of heat transfered to steam.
    turbo_generator_efficiency : float
        Fraction of steam heat converted to electricity.
    agent : UtilityAgent
        Steam produced.
    
    References
    ----------
    .. [1] Humbird, D., Davis, R., Tao, L., Kinchin, C., Hsu, D., Aden, A.,
        Dudgeon, D. (2011). Process Design and Economics for Biochemical 
        Conversion of Lignocellulosic Biomass to Ethanol: Dilute-Acid 
        Pretreatment and Enzymatic Hydrolysis of Corn Stover
        (No. NREL/TP-5100-47764, 1013269). https://doi.org/10.2172/1013269
    
    """
    network_priority = 0
    
    # TODO: Make this a whole system instead of approximating duty per mol values
    duty_over_mol = 40000 # Superheat steam with 40000 kJ/kmol
    boiler_blowdown = 0.03
    RO_rejection = 0
    _N_ins = 3
    _N_outs = 2
    _N_heat_utilities = 2
    _units = {'Flow rate': 'kg/hr',
              'Work': 'kW'}
    
    def __init__(self, ID='', ins=None, outs=(), *,
                 boiler_efficiency=0.80,
                 turbogenerator_efficiency=0.85,
                 side_steam=None,
                 agent=HeatUtility.get_heating_agent('low_pressure_steam')):
        Unit.__init__(self, ID, ins, outs)
        self.agent = agent
        self.makeup_water = makeup_water = agent.to_stream('boiler_makeup_water')
        loss = makeup_water.flow_proxy()
        loss.ID = 'rejected_water_and_blowdown'
        self.ins[-1] = makeup_water
        self.outs[-1] = loss
        self.boiler_efficiency = boiler_efficiency
        self.turbogenerator_efficiency = turbogenerator_efficiency
        self.steam_utilities = set()
        self.steam_demand = agent.to_stream()
        self.side_steam = side_steam
    
    def _run(self): pass
    
    def _design(self):
        B_eff = self.boiler_efficiency
        TG_eff = self.turbogenerator_efficiency
        steam = self.steam_demand
        steam_utilities = self.steam_utilities
        agent_ID = self.agent.ID
        if not steam_utilities:
            for u in self.system.units:
                if u is self: continue
                for hu in u.heat_utilities:
                    if hu.ID == agent_ID:
                        steam_utilities.add(hu)
        steam.imol['7732-18-5'] = steam_mol = sum([i.flow for i in steam_utilities])
        duty_over_mol = self.duty_over_mol
        feed_solids, feed_gas, _ = self.ins
        emissions, _ = self.outs
        hu_cooling, hu_steam = self.heat_utilities
        H_steam = steam_mol * duty_over_mol
        if self.side_steam: 
            H_steam += self.side_steam.H
        LHV_feed = 0
        for feed in (feed_solids, feed_gas):
            if feed: LHV_feed -= feed.LHV
        
        # This is simply for the mass balance (no special purpose)
        # TODO: In reality, this should be CO2
        emissions_mol = emissions.mol
        emissions_mol[:] = 0
        if feed_solids:
            emissions_mol += feed_solids.mol
        if feed_gas: 
            emissions_mol += feed_gas.mol
        combustion_rxns = self.chemicals.get_combustion_reactions()
        combustion_rxns(emissions_mol)
        emissions.imol['O2'] = 0
        emissions.T = 373.15
        emissions.P = 101325
        emissions.phase = 'g'
        H_content = LHV_feed*B_eff - emissions.H
        
        #: [float] Total steam produced by the boiler (kmol/hr)
        self.total_steam = H_content / duty_over_mol 
        
        self.makeup_water.imol['7732-18-5'] = (
            self.total_steam * self.boiler_blowdown * 1/(1-self.RO_rejection)
        )
        # Heat available for the turbogenerator
        H_electricity = H_content - H_steam
        
        Design = self.design_results
        Design['Flow rate'] = self.total_steam * 18.01528
        
        if H_electricity < 0:
            H_steam = H_content
            cooling = electricity = 0
        else:
            electricity = H_electricity * TG_eff
            cooling = electricity - H_electricity
        hu_cooling(cooling, steam.T)
        hu_steam.agent = self.agent
        hu_steam.cost = -sum([i.cost for i in steam_utilities])
        Design['Work'] = electricity/3600

    def _end_decorated_cost_(self):
        self.power_utility(self.power_utility.rate - self.design_results['Work'])

# Simulation of ethanol production from sugarcane
# in Brazil: economic study of an autonomous
# distillery
# Marina O.S. Diasa,b, Marcelo P. Cunhaa, Charles D.F. Jesusa, Mirna I.G.
# Scandiffioa, Carlos E.V. Rossella,b, Rubens Maciel Filhob, Antonio Bonomia
# a CTBE – Bioethanol Science and Technology National Laboratory, PO Box 6170 –
# CEP 13083-970, Campinas – SP, Brazil, marina.dias@bioetanol.org.br
# b School of Chemical Engineering, University of Campinas, PO Box 6066 – CEP
# 13083-970, Campinas – SP, Brazil
# LHV = 7565 # Bagasse lower heating value (kJ/kg)