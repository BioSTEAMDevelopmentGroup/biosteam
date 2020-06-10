# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from ... import HeatUtility
from . import Facility
from ..decorators import cost
import flexsolve as flx
import biosteam as bst

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
        
        [3] Natural gas to satisfy steam and power requirement.
        
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
    boiler_blowdown = 0.03
    RO_rejection = 0
    _N_ins = 4
    _N_outs = 2
    _N_heat_utilities = 0
    _units = {'Flow rate': 'kg/hr',
              'Work': 'kW'}
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 boiler_efficiency=0.80,
                 turbogenerator_efficiency=0.85,
                 side_steam=None,
                 agent=HeatUtility.get_heating_agent('low_pressure_steam'),
                 other_agents = ()):
        Facility.__init__(self, ID, ins, outs, thermo)
        self.agent = agent
        self.makeup_water = makeup_water = agent.to_stream('boiler_makeup_water')
        loss = makeup_water.flow_proxy()
        loss.ID = 'rejected_water_and_blowdown'
        self.ins[-1] = bst.Stream(None, thermo=self.thermo, phase='g')
        self.ins[-2] = makeup_water
        self.outs[-1] = loss
        self.boiler_efficiency = boiler_efficiency
        self.turbogenerator_efficiency = turbogenerator_efficiency
        self.steam_utilities = set()
        self.power_utilities = set()
        self.steam_demand = agent.to_stream()
        self.side_steam = side_steam
        self.other_agents = other_agents
    
    def _run(self): pass
    
    def _load_utility_agents(self):
        steam_utilities = self.steam_utilities
        steam_utilities.clear()
        agent = self.agent
        units = self.system.units.copy()
        units.remove(self)
        for agent in (*self.other_agents, agent):
            for u in units:
                for hu in u.heat_utilities:
                    if hu.agent is agent:
                        steam_utilities.add(hu)
        self.electricity_demand = sum([u.power_utility.rate for u in units if u.power_utility])

    def _design(self):
        B_eff = self.boiler_efficiency
        TG_eff = self.turbogenerator_efficiency
        steam_demand = self.steam_demand
        Design = self.design_results
        self._load_utility_agents()
        mol_steam = sum([i.flow for i in self.steam_utilities])
        feed_solids, feed_gas, makeup_water, feed_CH4 = self.ins
        emissions, _ = self.outs
        H_steam =  sum([i.duty for i in self.steam_utilities])
        side_steam = self.side_steam
        if side_steam: 
            H_steam += side_steam.H
            mol_steam += side_steam.F_mol
        steam_demand.imol['7732-18-5'] = mol_steam 
        if not mol_steam:
            duty_over_mol = self.agent.Hvap
        else:
            duty_over_mol = H_steam / mol_steam
        emissions_mol = emissions.mol
        emissions.T = self.agent.T
        emissions.P = 101325
        emissions.phase = 'g'
        combustion_rxns = self.chemicals.get_combustion_reactions()
        non_empty_feeds = [i for i in (feed_solids, feed_gas) if not i.isempty()]
        
        def calculate_excess_electricity_at_natual_gas_flow(natural_gas_flow):
            if natural_gas_flow:
                natural_gas_flow = abs(natural_gas_flow)
                feed_CH4.imol['CH4'] = natural_gas_flow
            else:
                feed_CH4.empty()
            H_combustion = feed_CH4.H - feed_CH4.HHV
            emissions_mol[:] = feed_CH4.mol
            for feed in non_empty_feeds:
                H_combustion += feed.H - feed.HHV
                emissions_mol[:] += feed.mol
            
            combustion_rxns.force_reaction(emissions_mol)
            emissions.imol['O2'] = 0
            H_content = B_eff * H_combustion - emissions.H
            
            #: [float] Total steam demand produced by the boiler (kmol/hr)
            self.total_steam = H_content / duty_over_mol 
            Design['Flow rate'] = flow_rate = self.total_steam * 18.01528
            
            # Heat available for the turbogenerator
            H_electricity = H_content - H_steam
            
            if H_electricity < 0:
                self.cooling_duty = electricity = 0
            else:
                electricity = H_electricity * TG_eff
                self.cooling_duty = electricity - H_electricity
            
            Design['Work'] = work = electricity/3600
            boiler = self.cost_items['Boiler']
            rate_boiler = boiler.kW * flow_rate / boiler.S
            return work - self.electricity_demand - rate_boiler
        
        excess_electricity = calculate_excess_electricity_at_natual_gas_flow(0)
        if excess_electricity < 0:
            f = calculate_excess_electricity_at_natual_gas_flow
            flx.aitken_secant(f, 10, 20, xtol=1, ytol=1)
        
        hu_cooling = bst.HeatUtility()
        hu_cooling(self.cooling_duty, steam_demand.T)
        hus_heating = bst.HeatUtility.sum_by_agent(tuple(self.steam_utilities))
        for hu in hus_heating: hu.reverse()
        self.heat_utilities = (*hus_heating, hu_cooling)
        
        makeup_water.imol['7732-18-5'] = (
                self.total_steam * self.boiler_blowdown * 1/(1-self.RO_rejection)
        )

    def _cost(self):
        self._decorated_cost()
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
