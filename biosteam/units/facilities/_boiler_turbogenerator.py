# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
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
import thermosteam as tmo
from warnings import warn

__all__ = ('BoilerTurbogenerator', 'Boiler')

#: TODO add reference of NREL

@cost('Flow rate', 'Hot process water softener system', 
      CE=551, cost=78e3, S=235803, n=0.6, BM=1.8)
@cost('Flow rate', 'Amine addition pkg', 
      CE=551, cost=40e3, S=235803, n=0.6, BM=1.8)
@cost('Flow rate', 'Deaerator',
      CE=551, cost=305e3, S=235803, n=0.6, BM=3.0)
@cost('Flow rate', 'Boiler',
      CE=551, cost=28550e3, kW=1371, S=238686, n=0.6, BM=1.8)
@cost('Ash disposal', 'Baghouse bags',
      CE=551, cost=466183. / 4363., n=1.0, lifetime=5)
class Boiler(Facility):
    """
    Create a Boiler object that will produce steam from burning feeds.
    All capital cost correlations are based on [1]_.
    
    Parameters
    ----------
    ins : stream sequence
        [0] Liquid/solid feed that will be burned.
        
        [1] Gas feed that will be burned.
        
        [2] Make-up water. 
        
        [3] Natural gas to satisfy steam and power requirement.
        
        [4] Lime for flue gas desulfurization.
        
        [5] Boiler chemicals.
        
    outs : stream sequence
        [0] Total emissions produced.
        
        [1] Blowdown water.
        
        [2] Ash disposal.
        
        [3] Remainder feed solids.
        
        [4] Remainder feed gas.
        
    boiler_efficiency : float
        Fraction of heat transfered to steam.
    agent : UtilityAgent, optional
        Steam produced. Defaults to low pressure steam.
    other_agents = () : Iterable[UtilityAgent]
        Other steams produced.
    
    Notes
    -----
    The flow rate of natural gas, lime, and boiler chemicals (streams 3-5)
    is set by the BoilerTurbogenerator object during simulation.
    
    References
    ----------
    .. [1] Humbird, D., Davis, R., Tao, L., Kinchin, C., Hsu, D., Aden, A.,
        Dudgeon, D. (2011). Process Design and Economics for Biochemical 
        Conversion of Lignocellulosic Biomass to Ethanol: Dilute-Acid 
        Pretreatment and Enzymatic Hydrolysis of Corn Stover
        (No. NREL/TP-5100-47764, 1013269). https://doi.org/10.2172/1013269
    
    """
    ticket_name = 'BT'
    network_priority = 0
    boiler_blowdown = 0.03
    RO_rejection = 0
    _N_ins = 6
    _N_outs = 5
    _N_heat_utilities = 0
    _units = {'Flow rate': 'kg/hr',
              'Work': 'kW',
              'Ash disposal': 'kg/hr'}
    
    def __init__(self, ID='', ins=None, 
                 outs=('emissions',
                       'rejected_water_and_blowdown',
                       'ash_disposal',
                       'remainder_feed_solids',
                       'remainder_feed_gas'),
                 thermo=None, *,
                 boiler_efficiency=0.80,
                 side_steam=None,
                 agent=None,
                 other_agents = ()):
        Facility.__init__(self, ID, ins, outs, thermo)
        self.agent = agent = agent or HeatUtility.get_heating_agent('low_pressure_steam')
        self.define_utility('Natural gas', self.natural_gas)
        self.define_utility('Ash disposal', self.ash_disposal)
        self.boiler_efficiency = boiler_efficiency
        self.steam_utilities = set()
        self.power_utilities = set()
        self.steam_demand = agent.to_stream()
        self.side_steam = side_steam
        self.other_agents = other_agents
        chemicals = self.chemicals
        if 'SO2' in chemicals:
            CAS_lime = '1305-62-0'
            if CAS_lime in chemicals or 'Ca(OH)2' in chemicals:
                if 'Ca(OH)2' not in chemicals:
                    chemicals.set_synonym(CAS_lime, 'Ca(OH)2')
                self.desulfurization_reaction =  tmo.Reaction(
                    'SO2 + Ca(OH)2 + 0.5 O2 -> CaSO4 + H2O', 'SO2', 0.92
                )
                self._ID_lime = 'Ca(OH)2'
                return
            CAS_lime = '1305-78-8'
            if CAS_lime in chemicals or 'CaO' in chemicals:
                if 'CaO' not in chemicals:
                    chemicals.set_synonym(CAS_lime, 'CaO')
                self.desulfurization_reaction =  tmo.Reaction(
                    'SO2 + CaO + 0.5 O2 -> CaSO4', 'SO2', 0.92
                )
                self._ID_lime = 'CaO'
                return
    
    @property
    def makeup_water(self):
        """[Stream] Makeup water due to boiler blowdown."""
        return self.ins[2]
    
    @property
    def natural_gas(self):
        """[Stream] Natural gas to satisfy steam and electricity requirements."""
        return self.ins[3]
    
    @property
    def ash_disposal(self):
        """[Stream] Ash disposal."""
        return self.outs[2]
    
    def _run(self): pass

    def _load_utility_agents(self):
        steam_utilities = self.steam_utilities
        steam_utilities.clear()
        agent = self.agent
        units = self.other_units
        for agent in (*self.other_agents, agent):
            ID = agent.ID
            for u in units:
                for hu in u.heat_utilities:
                    agent = hu.agent
                    if agent and agent.ID == ID:
                        steam_utilities.add(hu)

    def _design(self):
        B_eff = self.boiler_efficiency
        steam_demand = self.steam_demand
        Design = self.design_results
        self._load_utility_agents()
        mol_steam = sum([i.flow for i in self.steam_utilities])
        feed_solids, feed_gas, makeup_water, feed_CH4, lime, chems = self.ins
        emissions, blowdown_water, ash_disposal, remainder_feed_solids, remainder_feed_gas = self.outs
        chemicals = emissions.chemicals
        if not lime.price:
            lime.price = 0.19937504680689402
        if not chems.price:
            chems.price = 4.995862254032183
        H_steam =  sum([i.duty for i in self.steam_utilities])
        side_steam = self.side_steam
        if side_steam: 
            H_steam += side_steam.H
            mol_steam += side_steam.F_mol
        steam_demand.imol['7732-18-5'] = mol_steam 
        duty_over_mol = 39000 # kJ / mol-superheated steam 
        emissions_mol = emissions.mol
        emissions.T = self.agent.T
        emissions.P = 101325
        emissions.phase = 'g'
        self.combustion_reactions = combustion_rxns = chemicals.get_combustion_reactions()
        non_empty_feeds = [i for i in (feed_solids, feed_gas) if not i.isempty()]
        
        def calculate_excess_heat_at_natual_gas_flow(natural_gas_flow):
            if natural_gas_flow:
                natural_gas_flow = abs(natural_gas_flow)
                feed_CH4.imol['CH4'] = natural_gas_flow
            else:
                feed_CH4.empty()
            H_combustion = feed_CH4.H + feed_CH4.HHV
            emissions_mol[:] = feed_CH4.mol
            for feed in non_empty_feeds:
                H_combustion += feed.H + feed.HHV
                emissions_mol[:] += feed.mol
            
            combustion_rxns.force_reaction(emissions_mol)
            emissions.imol['O2'] = 0
            
            H_content = B_eff * H_combustion - emissions.H
            
            #: [float] Total steam produced by the boiler (kmol/hr)
            self.total_steam = H_content / duty_over_mol 
            Design['Flow rate'] = self.total_steam * 18.01528
            
            # Excess heat available
            H_excess = H_content - H_steam
            return H_excess
        
        remainder_feed_gas.empty()
        remainder_feed_solids.empty()
        gas_fraction_burned = solids_fraction_burned = 1.
        excess_heat = calculate_excess_heat_at_natual_gas_flow(0)
        if excess_heat < -1:
            f = calculate_excess_heat_at_natual_gas_flow
            lb = 0.
            ub = 100.
            while f(ub) < 0.: 
                lb = ub
                ub *= 2
            flx.IQ_interpolation(f, lb, ub, xtol=1, ytol=1)
        elif excess_heat > 1:
            def calculate_excess_heat_with_diverted_gas(fraction_burned):
                H_combustion = fraction_burned*(feed_gas.H + feed_gas.HHV) + (feed_solids.H + feed_solids.HHV)
                emissions_mol[:] = fraction_burned*feed_gas.mol + feed_solids.mol
                combustion_rxns.force_reaction(emissions_mol)
                emissions.imol['O2'] = 0
                H_content = B_eff * H_combustion - emissions.H
                
                #: [float] Total steam produced by the boiler (kmol/hr)
                self.total_steam = H_content / duty_over_mol 
                Design['Flow rate'] = self.total_steam * 18.01528
                
                # Excess heat available
                H_excess = H_content - H_steam
                return H_excess
            if f(0.) > 1.: 
                gas_fraction_burned = 0.
            else:
                f = calculate_excess_heat_with_diverted_gas
                gas_fraction_burned = flx.IQ_interpolation(f, 0., 1., xtol=1, ytol=1)
            
            if gas_fraction_burned == 0.:
                def calculate_excess_heat_with_diverted_solids(fraction_burned):
                    H_combustion = fraction_burned * (feed_solids.H + feed_solids.HHV)
                    emissions_mol[:] = fraction_burned * feed_solids.mol
                    combustion_rxns.force_reaction(emissions_mol)
                    emissions.imol['O2'] = 0
                    H_content = B_eff * H_combustion - emissions.H
                    
                    #: [float] Total steam produced by the boiler (kmol/hr)
                    self.total_steam = H_content / duty_over_mol 
                    Design['Flow rate'] = self.total_steam * 18.01528
                    
                    # Excess heat available
                    H_excess = H_content - H_steam
                    return H_excess
                f = calculate_excess_heat_with_diverted_solids
                solids_fraction_burned = flx.IQ_interpolation(f, 0., 1., xtol=1, ytol=1)
            
        remainder_feed_solids.mol = (1. - solids_fraction_burned) * feed_solids.mol
        remainder_feed_gas.mol = (1. - gas_fraction_burned) * feed_gas.mol
        hus_heating = bst.HeatUtility.sum_by_agent(tuple(self.steam_utilities))
        for hu in hus_heating: hu.reverse()
        self.heat_utilities = tuple(hus_heating)
        water_index = chemicals.index('7732-18-5')
        blowdown_water.mol[water_index] = makeup_water.mol[water_index] = (
                self.total_steam * self.boiler_blowdown * 1/(1-self.RO_rejection)
        )
        ash_IDs = [i.ID for i in self.chemicals if not i.formula]
        emissions_mol = emissions.mol
        if 'SO2' in chemicals: 
            ash_IDs.append('CaSO4')
            lime_index = chemicals.index(self._ID_lime)
            self.desulfurization_reaction.force_reaction(emissions)
            # FGD lime scaled based on SO2 generated,	
            # 20% stoichiometetric excess based on P52 of ref [1]
            
            lime.mol[lime_index] = lime_mol = max(0, - emissions_mol[lime_index] * 1.2)
            emissions_mol[emissions_mol < 0.] = 0.
        # About 0.4536 kg/hr of boiler chemicals are needed per 234484 kg/hr steam produced
        chems.imol['Ash'] = boiler_chems = 1.9345e-06 * Design['Flow rate']
        ash_disposal.empty()
        ash_disposal.copy_flow(emissions, IDs=tuple(ash_IDs), remove=True)
        ash_disposal.imol['Ash'] += boiler_chems
        dry_ash = ash_disposal.F_mass
        ash_disposal.imass['Water'] = moisture = dry_ash * 0.3 # ~20% moisture
        if 'SO2' in chemicals:
            if self._ID_lime == '1305-62-0': # Ca(OH)2
                lime.imol['Water'] = 4 * lime_mol # Its a slurry
            else: # CaO
                lime.imol['Water'] = 5 * lime_mol 
        

@cost('Work', 'Turbogenerator',
      CE=551, S=42200, kW=0, cost=9500e3, n=0.60, BM=1.8)
@cost('Flow rate', 'Hot process water softener system', 
      CE=551, cost=78e3, S=235803, n=0.6, BM=1.8)
@cost('Flow rate', 'Amine addition pkg', 
      CE=551, cost=40e3, S=235803, n=0.6, BM=1.8)
@cost('Flow rate', 'Deaerator',
      CE=551, cost=305e3, S=235803, n=0.6, BM=3.0)
@cost('Flow rate', 'Boiler',
      CE=551, cost=28550e3, kW=1371, S=238686, n=0.6, BM=1.8)
@cost('Ash disposal', 'Baghouse bags',
      CE=551, cost=466183. / 4363., n=1.0, lifetime=5)
class BoilerTurbogenerator(Facility):
    """
    Create a BoilerTurbogenerator object that will calculate electricity
    generation from burning the feed. It also takes into account how much
    steam is being produced, and the required cooling utility of the turbo
    generator. All capital cost correlations are based on [1]_.
    
    Parameters
    ----------
    ins : stream sequence
        [0] Liquid/solid feed that will be burned.
        
        [1] Gas feed that will be burned.
        
        [2] Make-up water. 
        
        [3] Natural gas to satisfy steam and power requirement.
        
        [4] Lime for flue gas desulfurization.
        
        [5] Boiler chemicals.
        
    outs : stream sequence
        [0] Total emissions produced.
        
        [1] Blowdown water.
        
        [2] Ash disposal.
        
    boiler_efficiency : float
        Fraction of heat transfered to steam.
    turbo_generator_efficiency : float
        Fraction of steam heat converted to electricity.
    agent : UtilityAgent, optional
        Steam produced. Defaults to low pressure steam.
    other_agents = () : Iterable[UtilityAgent]
        Other steams produced.
    natural_gas_price = 0.218 : float
        Price of natural gas [USD/kg].
    ash_disposal_price = -0.0318 : float
        Price of disposing ash [USD/kg].
        
    Notes
    -----
    The flow rate of natural gas, lime, and boiler chemicals (streams 3-5)
    is set by the BoilerTurbogenerator object during simulation.
    
    References
    ----------
    .. [1] Humbird, D., Davis, R., Tao, L., Kinchin, C., Hsu, D., Aden, A.,
        Dudgeon, D. (2011). Process Design and Economics for Biochemical 
        Conversion of Lignocellulosic Biomass to Ethanol: Dilute-Acid 
        Pretreatment and Enzymatic Hydrolysis of Corn Stover
        (No. NREL/TP-5100-47764, 1013269). https://doi.org/10.2172/1013269
    
    """
    ticket_name = 'BT'
    network_priority = 0
    boiler_blowdown = 0.03
    RO_rejection = 0
    _N_ins = 6
    _N_outs = 3
    _N_heat_utilities = 0
    _units = {'Flow rate': 'kg/hr',
              'Work': 'kW',
              'Ash disposal': 'kg/hr'}
    
    def __init__(self, ID='', ins=None, 
                 outs=('emissions',
                       'rejected_water_and_blowdown',
                       'ash_disposal'),
                 thermo=None, *,
                 boiler_efficiency=0.80,
                 turbogenerator_efficiency=0.85,
                 side_steam=None,
                 agent=None,
                 other_agents = (),
                 natural_gas_price=0.218,
                 ash_disposal_price=-0.0318):
        Facility.__init__(self, ID, ins, outs, thermo)
        self.agent = agent = agent or HeatUtility.get_heating_agent('low_pressure_steam')
        self.define_utility('Natural gas', self.natural_gas)
        self.define_utility('Ash disposal', self.ash_disposal)
        self.boiler_efficiency = boiler_efficiency
        self.turbogenerator_efficiency = turbogenerator_efficiency
        self.steam_utilities = set()
        self.power_utilities = set()
        self.steam_demand = agent.to_stream()
        self.side_steam = side_steam
        self.other_agents = other_agents
        self._load_components()
        
    def _load_components(self):
        chemicals = self.chemicals
        if 'SO2' in chemicals:
            CAS_lime = '1305-62-0'
            if CAS_lime in chemicals or 'Ca(OH)2' in chemicals:
                if 'Ca(OH)2' not in chemicals:
                    chemicals.set_synonym(CAS_lime, 'Ca(OH)2')
                self.desulfurization_reaction =  tmo.Reaction(
                    'SO2 + Ca(OH)2 + 0.5 O2 -> CaSO4 + H2O', 'SO2', 0.92, chemicals
                )
                self._ID_lime = 'Ca(OH)2'
                return
            CAS_lime = '1305-78-8'
            if CAS_lime in chemicals or 'CaO' in chemicals:
                if 'CaO' not in chemicals:
                    chemicals.set_synonym(CAS_lime, 'CaO')
                self.desulfurization_reaction =  tmo.Reaction(
                    'SO2 + CaO + 0.5 O2 -> CaSO4', 'SO2', 0.92, chemicals
                )
                self._ID_lime = 'CaO'
                return
    
    @property
    def makeup_water(self):
        """[Stream] Makeup water due to boiler blowdown."""
        return self.ins[2]
    
    @property
    def natural_gas(self):
        """[Stream] Natural gas to satisfy steam and electricity requirements."""
        return self.ins[3]
    
    @property
    def ash_disposal(self):
        """[Stream] Ash disposal."""
        return self.outs[2]
    
    def _run(self): pass

    def _load_utility_agents(self):
        steam_utilities = self.steam_utilities
        steam_utilities.clear()
        agent = self.agent
        units = self.other_units
        for agent in (*self.other_agents, agent):
            ID = agent.ID
            for u in units:
                for hu in u.heat_utilities:
                    agent = hu.agent
                    if agent and agent.ID == ID:
                        steam_utilities.add(hu)
        self.electricity_demand = sum([u.power_utility.consumption for u in units])

    def _design(self):
        B_eff = self.boiler_efficiency
        TG_eff = self.turbogenerator_efficiency
        steam_demand = self.steam_demand
        Design = self.design_results
        chemicals = self.chemicals
        self._load_utility_agents()
        mol_steam = sum([i.flow for i in self.steam_utilities])
        feed_solids, feed_gas, makeup_water, feed_CH4, lime, chems = self.ins
        emissions, blowdown_water, ash_disposal = self.outs
        if not lime.price:
            lime.price = 0.19937504680689402
        if not chems.price:
            chems.price = 4.995862254032183
        H_steam =  sum([i.duty for i in self.steam_utilities])
        side_steam = self.side_steam
        if side_steam: 
            H_steam += side_steam.H
            mol_steam += side_steam.F_mol
        steam_demand.imol['7732-18-5'] = mol_steam 
        duty_over_mol = 39000 # kJ / mol-superheated steam 
        emissions_mol = emissions.mol
        emissions.T = self.agent.T
        emissions.P = 101325
        emissions.phase = 'g'
        self.combustion_reactions = combustion_rxns = chemicals.get_combustion_reactions()
        non_empty_feeds = [i for i in (feed_solids, feed_gas) if not i.isempty()]
        
        def calculate_excess_electricity_at_natual_gas_flow(natural_gas_flow):
            if natural_gas_flow:
                natural_gas_flow = abs(natural_gas_flow)
                feed_CH4.imol['CH4'] = natural_gas_flow
            else:
                feed_CH4.empty()
            H_combustion = feed_CH4.H + feed_CH4.HHV
            emissions_mol[:] = feed_CH4.mol
            for feed in non_empty_feeds:
                H_combustion += feed.H + feed.HHV
                emissions_mol[:] += feed.mol
            
            combustion_rxns.force_reaction(emissions_mol)
            emissions.imol['O2'] = 0
            
            H_content = B_eff * H_combustion - emissions.H
            
            #: [float] Total steam produced by the boiler (kmol/hr)
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
            lb = 0.
            ub = - excess_electricity * 3600 / feed_CH4.chemicals.CH4.LHV
            while f(ub) < 0.: 
                lb = ub
                ub *= 2
            flx.IQ_interpolation(f, lb, ub, xtol=1, ytol=1)
        
        hu_cooling = bst.HeatUtility()
        hu_cooling(self.cooling_duty, steam_demand.T)
        hus_heating = bst.HeatUtility.sum_by_agent(tuple(self.steam_utilities))
        for hu in hus_heating: hu.reverse()
        self.heat_utilities = (*hus_heating, hu_cooling)
        water_index = chemicals.index('7732-18-5')
        blowdown_water.mol[water_index] = makeup_water.mol[water_index] = (
                self.total_steam * self.boiler_blowdown * 1/(1-self.RO_rejection)
        )
        ash_IDs = [i.ID for i in self.chemicals if not i.formula]
        emissions_mol = emissions.mol
        if 'SO2' in chemicals: 
            ash_IDs.append('CaSO4')
            lime_index = emissions.chemicals.index(self._ID_lime)
            sulfur_index = emissions.chemicals.index('CaSO4')
            self.desulfurization_reaction.force_reaction(emissions)
            # FGD lime scaled based on SO2 generated,	
            # 20% stoichiometetric excess based on P52 of ref [1]
            
            lime.mol[lime_index] = lime_mol = max(0, emissions_mol[sulfur_index] * 1.2)
            emissions_mol[emissions_mol < 0.] = 0.
        else:
            lime.empty()
        # About 0.4536 kg/hr of boiler chemicals are needed per 234484 kg/hr steam produced
        chems.imol['Ash'] = boiler_chems = 1.9345e-06 * Design['Flow rate']
        ash_disposal.empty()
        ash_disposal.copy_flow(emissions, IDs=tuple(ash_IDs), remove=True)
        ash_disposal.imol['Ash'] += boiler_chems
        dry_ash = ash_disposal.F_mass
        ash_disposal.imass['Water'] = moisture = dry_ash * 0.3 # ~20% moisture
        Design['Ash disposal'] = dry_ash + moisture
        if 'SO2' in chemicals:
            if self._ID_lime == '1305-62-0': # Ca(OH)2
                lime.imol['Water'] = 4 * lime_mol # Its a slurry
            else: # CaO
                lime.imol['Water'] = 5 * lime_mol 
        
    def _cost(self):
        self._decorated_cost()
        self.power_utility.production = self.design_results['Work']

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
