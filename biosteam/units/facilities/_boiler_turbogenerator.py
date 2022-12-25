# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
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
        Fraction of heat transferred to steam.
    turbo_generator_efficiency : float
        Fraction of steam heat converted to electricity.
    agent : UtilityAgent, optional
        Steam produced. Defaults to low pressure steam.
    other_agents = () : Iterable[UtilityAgent]
        Other steams produced.
    natural_gas_price : float
        Price of natural gas [USD/kg]. Same as `bst.stream_utility_prices['Natural gas']`,
        defaults to 0.218.
    ash_disposal_price : float
        Price of disposing ash [USD/kg]. Same as `bst.stream_utility_prices['Ash disposal']`,
        defaults to -0.0318.
    satisfy_system_electricity_demand : bool
        Whether to purchase natural gas to satisfy system electricity demand
        if there is not enough heat from the feed waste and gas.
        If True, will purchase natural gas to satisfy system heat and electricity demand
        (even if there is not enough heat from the feed wastes and gas);
        if False, will only purchase natural gas to satisfy system heat demand
        (i.e., electricity will be purchased from the grid if there is not
         enough heat from the feed wastes and gas).
        In either case, if there is excess heat from the feed wastes and gas,
        electricity will still be produced
        (i.e., this arg only affects the calculation of natural gas flow).
        
        
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
    default_boiler_efficiency_basis = 'LHV and sensible heat'
    _N_ins = 6
    _N_outs = 3
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
                 natural_gas_price=None,
                 ash_disposal_price=None,
                 T_emissions=None,
                 satisfy_system_electricity_demand=True,
                 boiler_efficiency_basis=None,
        ):
        if boiler_efficiency_basis is None:
            boiler_efficiency_basis = self.default_boiler_efficiency_basis
        self.boiler_efficiency_basis = boiler_efficiency_basis
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
        self.T_emissions = self.agent.T if T_emissions is None else T_emissions # Assume no heat integration
        if natural_gas_price is not None: self.natural_gas_price = natural_gas_price
        if ash_disposal_price is not None: self.ash_disposal_price = ash_disposal_price
        self.satisfy_system_electricity_demand = satisfy_system_electricity_demand
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
    def blowdown_water(self):
        return self.outs[1]
    
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
    
    @property
    def natural_gas_price(self):
        """[Float] Price of natural gas, same as `bst.stream_utility_prices['Natural gas']`."""
        return bst.stream_utility_prices['Natural gas']
    
    @natural_gas_price.setter
    def natural_gas_price(self, new_price):
        bst.stream_utility_prices['Natural gas'] = new_price
    
    @property
    def ash_disposal_price(self):
        """[Float] Price of ash disposal, same as `bst.stream_utility_prices['Ash disposal']`."""
        return bst.stream_utility_prices['Ash disposal']
    
    @ash_disposal_price.setter
    def ash_disposal_price(self, ash_disposal_price):
        bst.stream_utility_prices['Ash disposal'] = ash_disposal_price
    
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
        feed_CH4.phase = 'g'
        feed_CH4.set_property('T', 60, 'degF')
        feed_CH4.set_property('P', 14.73, 'psi')
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
        emissions.T = self.T_emissions
        emissions.P = 101325
        emissions.phase = 'g'
        self.combustion_reactions = combustion_rxns = chemicals.get_combustion_reactions()
        non_empty_feeds = [i for i in (feed_solids, feed_gas) if not i.isempty()]
        boiler_efficiency_basis = self.boiler_efficiency_basis
        electricity_demand = self.electricity_demand if self.satisfy_system_electricity_demand else 0
        def calculate_excess_electricity_at_natual_gas_flow(natural_gas_flow):
            emissions_mol[:] = 0.
            if natural_gas_flow:
                natural_gas_flow = abs(natural_gas_flow)
                feed_CH4.imol['CH4'] = natural_gas_flow
            else:
                feed_CH4.empty()
            for feed in non_empty_feeds: emissions_mol[:] += feed.mol
            combustion_rxns.force_reaction(emissions_mol)
            emissions.imol['O2'] = 0
            if boiler_efficiency_basis == 'LHV and sensible heat':
                H_combustion = feed_CH4.HHV - emissions.H
                for feed in non_empty_feeds: H_combustion += feed.H + feed.HHV
            elif boiler_efficiency_basis == 'LHV':
                H_combustion = feed_CH4.LHV
                for feed in non_empty_feeds: H_combustion += feed.LHV
            elif boiler_efficiency_basis == 'HHV':
                H_combustion = feed_CH4.HHV
                for feed in non_empty_feeds: H_combustion += feed.HHV
            else:
                raise ValueError(
                    f"invalid boiler efficiency basis {boiler_efficiency_basis}; "
                    f"valid values include 'LHV and sensible heat', 'LHV', or 'HHV'"
                )
            H_content = B_eff * H_combustion
            #: [float] Total steam produced by the boiler (kmol/hr)
            self.total_steam = H_content / duty_over_mol 
            Design['Flow rate'] = flow_rate = self.total_steam * 18.01528
            
            # Heat available for the turbogenerator
            H_electricity = H_content - H_steam
            
            electricity = H_electricity * TG_eff  # Electricity produced
            self.cooling_duty = electricity - H_electricity
            
            Design['Work'] = work = electricity/3600
            boiler = self.cost_items['Boiler']
            rate_boiler = boiler.kW * flow_rate / boiler.S
            return work - electricity_demand - rate_boiler
        
        self._excess_electricity_without_natural_gas = excess_electricity = calculate_excess_electricity_at_natual_gas_flow(0)
        if excess_electricity < 0:
            f = calculate_excess_electricity_at_natual_gas_flow
            lb = 0.
            ub = - excess_electricity * 3600 / feed_CH4.chemicals.CH4.LHV
            while f(ub) < 0.: 
                lb = ub
                ub *= 2
            flx.IQ_interpolation(f, lb, ub, xtol=1, ytol=1)
        
        hu_cooling = bst.HeatUtility()
        try:
            hu_cooling(self.cooling_duty, steam_demand.T)
        except:
            breakpoint()
        hus_heating = bst.HeatUtility.sum_by_agent(tuple(self.steam_utilities))
        for hu in hus_heating: hu.reverse()
        self.heat_utilities = [*hus_heating, hu_cooling]
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
            # 20% stoichiometric excess based on P52 of ref [1]
            
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