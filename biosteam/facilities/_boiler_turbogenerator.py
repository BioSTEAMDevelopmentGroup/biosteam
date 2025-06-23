# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import flexsolve as flx
import biosteam as bst
import thermosteam as tmo
cost = bst.decorators.cost

__all__ = ('BoilerTurbogenerator', 'Boiler')

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
class BoilerTurbogenerator(bst.Facility):
    """
    Create a BoilerTurbogenerator object that will calculate electricity
    generation from burning the feed. It also takes into account how much
    steam is being produced, and the required cooling utility of the turbo
    generator. All capital cost correlations are based on [1]_.
    
    Parameters
    ----------
    ins : 
        * [0] Liquid/solid feed that will be burned.
        * [1] Gas feed that will be burned.
        * [2] Make-up water. 
        * [3] Natural gas/fuel to satisfy steam and electricity demand.
        * [4] Lime for flue gas desulfurization.
        * [5] Boiler chemicals.
        * [6] Air or oxygen-rich gas.
    outs : 
        * [0] Total emissions produced.
        * [1] Blowdown water.
        * [2] Ash disposal.
    boiler_efficiency : float, optional
        Fraction of heat transferred to steam. Defaults to 0.8.
    turbo_generator_efficiency : float, optional
        Fraction of steam heat converted to electricity. Defaults to 0.85.
    agent : UtilityAgent, optional
        Steam produced. Defaults to low pressure steam.
    other_agents = () : Iterable[UtilityAgent], optional
        Other steams produced. Defaults to all other heating agents.
    fuel_source : str, optional
        Name fuel used to satisfy steam and electricity demand. Defaults to 'CH4'.
    fuel_price : float, optional
        Price of natural gas [USD/kg]. Same as `bst.stream_utility_prices['Natural gas']`,
        which defaults to 0.218.
    ash_disposal_price : float, optional
        Price of disposing ash [USD/kg]. Same as `bst.stream_utility_prices['Ash disposal']`,
        which defaults to -0.0318.
    satisfy_system_electricity_demand : bool, optional
        Whether to purchase natural gas to satisfy system electricity demand
        if there is not enough heat from process feeds (i.e., inlets 0 and 1).
        If True, natural gas is purchased to satisfy system heat and electricity demand
        when there is not enough heat from the feed and gas.
        If False, natural gas is only purchased to satisfy system heat demand
        and electricity will be purchased from the grid if there is not
        enough heat from the feeds.
        In either case, if there is excess heat from the process feeds,
        electricity will still be produced.
    boiler_efficiency_basis : str, optional
        Basis of boiler efficiency. Defaults to 'LHV' (i.e., lower heating value).
        'HHV' (i.e., higher heating value) is also a valid basis. 
        
    Examples
    --------
    Create a boiler-turbogenerator system that uses sugarcane bagasse to 
    produce steam for a distillation unit and any excess steam for surplus electricity:
    
    >>> import biosteam as bst
    >>> from biorefineries import cane
    >>> chemicals = cane.create_sugarcane_chemicals()
    >>> chemicals.define_group(
    ...     name='Fiber',
    ...     IDs=['Cellulose', 'Hemicellulose', 'Lignin'],
    ...     composition=[0.4704 , 0.2775, 0.2520],
    ...     wt=True, # Composition is given as weight
    ... )
    >>> bst.settings.set_thermo(chemicals)
    >>> dilute_ethanol = bst.Stream('dilute_ethanol', Water=1390, Ethanol=590)
    >>> bagasse = bst.Stream('bagasse', Water=0.4, Fiber=0.6, total_flow=8e4, units='kg/hr')
    >>> with bst.System('sys') as sys:
    ...     D1 = bst.BinaryDistillation('D1', ins=dilute_ethanol, Lr=0.999, Hr=0.89, k=1.25, LHK=('Ethanol', 'Water'))
    ...     BT = bst.BoilerTurbogenerator('BT')
    ...     BT.ins[0] = bagasse
    >>> sys.simulate()
    >>> BT.show()
    BoilerTurbogenerator: BT
    ins...
    [0] bagasse  
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Water          1.78e+03
                        Cellulose      139
                        Hemicellulose  101
                        Lignin         79.5
    [1] -  
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow: 0
    [2] -  
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Water  488
    [3] -  
        phase: 'g', T: 288.71 K, P: 101560 Pa
        flow: 0
    [4] -  
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow: 0
    [5] -  
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Ash  0.567
    [6] -  
        phase: 'g', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): O2  9.85e+03
                        N2  4.23e+04
    outs...
    [0] emissions  
        phase: 'g', T: 394.15 K, P: 101325 Pa
        flow (kmol/hr): Water  3.19e+03
                        CO2    1.98e+03
                        O2     7.84e+03
                        N2     4.23e+04
    [1] rejected_water_and_blowdown  
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Water  488
    [2] ash_disposal  
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Water  0.00944
                        Ash    0.567
    
    >>> BT.results() # Steam and electricity are produced, so costs are negative
    Boiler turbogenerator                                      Units        BT
    Electricity           Power                                   kW -1.31e+05
                          Cost                                USD/hr -1.02e+04
    Low pressure steam    Duty                                 kJ/hr -7.32e+07
                          Flow                               kmol/hr -1.89e+03
                          Cost                                USD/hr      -450
    Cooling water         Duty                                 kJ/hr -8.42e+07
                          Flow                               kmol/hr  5.75e+04
                          Cost                                USD/hr      28.1
    Fuel (inlet)          Flow                                 kg/hr         0
                          Cost                                USD/hr         0
    Ash disposal (outlet) Flow                                 kg/hr     0.737
                          Cost                                USD/hr    0.0234
    Design                Work                                    kW  1.33e+05
                          Flow rate                            kg/hr  2.93e+05
                          Ash disposal                         kg/hr     0.737
    Purchase cost         Baghouse bags                          USD      81.1
                          Boiler                                 USD  3.33e+07
                          Deaerator                              USD  3.58e+05
                          Amine addition pkg                     USD  4.69e+04
                          Hot process water softener system      USD  9.16e+04
                          Turbogenerator                         USD  1.94e+07
    Total purchase cost                                          USD  5.32e+07
    Utility cost                                              USD/hr -1.07e+04
    
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
    network_priority = 1
    boiler_blowdown = 0.03
    # Reverse osmosis (RO) typically rejects 25% of water, but the boiler-feed water is assumed to come after RO.
    # Setting this parameter to a fraction more than zero effectively assumes that this unit accounts for reverse osmosis.
    RO_rejection = 0 
    _N_ins = 7
    _N_outs = 3
    _units = {'Flow rate': 'kg/hr',
              'Work': 'kW',
              'Ash disposal': 'kg/hr'}
    
    def __init__(self, ID='', ins=None, 
                 outs=('emissions',
                       'rejected_water_and_blowdown',
                       'ash_disposal'),
                 thermo=None, *,
                 boiler_efficiency=None,
                 turbogenerator_efficiency=None,
                 side_steam=None,
                 agent=None,
                 other_agents=None,
                 fuel_price=None,
                 natural_gas_price=None,
                 ash_disposal_price=None,
                 T_emissions=None,
                 CO2_emissions_concentration=None,
                 satisfy_system_electricity_demand=None,
                 boiler_efficiency_basis=None,
                 fuel_source=None,
                 oxygen_rich_gas_composition=None,
        ):
        if boiler_efficiency_basis is None: boiler_efficiency_basis = 'LHV'
        if boiler_efficiency is None: boiler_efficiency = 0.80
        if turbogenerator_efficiency is None: turbogenerator_efficiency = 0.85
        if satisfy_system_electricity_demand is None: satisfy_system_electricity_demand = True
        if fuel_source is None: fuel_source = 'CH4'
        if oxygen_rich_gas_composition is None: oxygen_rich_gas_composition = dict(O2=21, N2=79, phase='g', units='kg/hr')
        if CO2_emissions_concentration is None: CO2_emissions_concentration = 0.055  # Usually between 4 - 7 for biomass and natural gas (https://www.sciencedirect.com/science/article/pii/S0957582021005127)
        bst.Facility.__init__(self, ID, ins, outs, thermo)
        settings = bst.settings
        self.boiler_efficiency_basis = boiler_efficiency_basis
        self.agent = agent = agent or settings.get_heating_agent('low_pressure_steam')
        self.fuel_source = fuel_source
        self.define_utility('Fuel', self.fuel)
        self.define_utility('Ash disposal', self.ash_disposal)
        self.boiler_efficiency = boiler_efficiency
        self.turbogenerator_efficiency = turbogenerator_efficiency
        self.steam_utilities = []
        self.steam_demand = agent.to_stream()
        self.side_steam = side_steam
        self.other_agents = [i for i in settings.heating_agents if i is not agent] if other_agents is None else other_agents
        self.CO2_emissions_concentration = CO2_emissions_concentration
        self.oxygen_rich_gas_composition = oxygen_rich_gas_composition 
        if T_emissions is not None: 
            raise ValueError('setting T_emissions is not yet implemented')
            # T_emissions should dictate the efficiency of either the boiler or the turbogenerator.
            # Note that T_emissions is tied to the energy balance.
            self.T_emissions = T_emissions
        if natural_gas_price is not None: 
            self.fuel_price = natural_gas_price
        elif fuel_price is not None: 
            self.fuel_price = fuel_price
        if ash_disposal_price is not None: self.ash_disposal_price = ash_disposal_price
        self.satisfy_system_electricity_demand = satisfy_system_electricity_demand
      
    def _get_desulfurization_rxn_and_coreactant(self):
        try:
            return self.desulfurization_reaction, self._ID_lime
        except:
            chemicals = self.chemicals
            CAS_lime = '1305-62-0'
            has_common_name = 'Ca(OH)2' in chemicals
            if CAS_lime in chemicals or has_common_name:
                if not has_common_name: chemicals.set_synonym(CAS_lime, 'Ca(OH)2')
                self.desulfurization_reaction = rxn = tmo.Reaction(
                    'SO2 + Ca(OH)2 + 0.5 O2 -> CaSO4 + H2O', 'SO2', 0.92, chemicals
                )
                self._ID_lime = ID = 'Ca(OH)2'
                return rxn, ID
            CAS_lime = '1305-78-8'
            has_common_name = 'CaO' in chemicals
            if CAS_lime in chemicals or has_common_name:
                if not has_common_name: chemicals.set_synonym(CAS_lime, 'CaO')
                self.desulfurization_reaction = rxn = tmo.Reaction(
                    'SO2 + CaO + 0.5 O2 -> CaSO4', 'SO2', 0.92, chemicals
                )
                self._ID_lime = ID = 'CaO'
                return rxn, ID
        raise RuntimeError(
            "lime is required for boiler, but no chemical 'CaO' or 'Ca(OH)2' "
            "available in thermodynamic property package"
        )
        
    @property
    def emissions(self):
        return self.outs[0]
    
    @property
    def blowdown_water(self):
        return self.outs[1]
    
    @property
    def makeup_water(self):
        """[Stream] Makeup water due to boiler blowdown."""
        return self.ins[2]
    
    @property
    def fuel(self):
        """[Stream] Fuel used to satisfy steam and electricity requirements."""
        return self.ins[3]
    natural_gas = fuel
    
    @property
    def air(self):
        """[Stream] Air or oxygen rich gas used to supply oxygen for combustion."""
        return self.ins[6]
    oxygen_rich_gas = air
    
    @property
    def ash_disposal(self):
        """[Stream] Ash disposal."""
        return self.outs[2]
    
    @property
    def fuel_price(self):
        """[Float] Price of natural gas, same as `bst.stream_utility_prices['Natural gas']`."""
        return bst.stream_utility_prices['Fuel']
    @fuel_price.setter
    def fuel_price(self, new_price):
        bst.stream_utility_prices['Fuel'] = new_price
    natural_gas_price = fuel_price
    
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
                        steam_utilities.append(hu)
        self.electricity_demand = sum([u.power_utility.consumption for u in units])
    
    def _design(self):
        B_eff = self.boiler_efficiency
        TG_eff = self.turbogenerator_efficiency
        steam_demand = self.steam_demand
        Design = self.design_results
        chemicals = self.chemicals
        self._load_utility_agents()
        mol_steam = sum([i.flow for i in self.steam_utilities])
        feed_solids, feed_gas, makeup_water, fuel, lime, chems, oxygen_rich_gas = self.ins
        oxygen_rich_gas.empty()
        if self.fuel_source == 'CH4':
            fuel.phase = 'g'
            fuel.set_property('T', 60, 'degF')
            fuel.set_property('P', 14.73, 'psi')
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
        self.combustion_reactions = combustion_rxns = chemicals.get_combustion_reactions()
        non_empty_feeds = [i for i in (feed_solids, feed_gas) if not i.isempty()]
        boiler_efficiency_basis = self.boiler_efficiency_basis
        fuel_source = self.fuel_source
        def calculate_excess_electricity_at_natual_gas_flow(fuel_flow):
            if fuel_flow:
                fuel_flow = abs(fuel_flow)
                fuel.imol[fuel_source] = fuel_flow
            else:
                fuel.empty()
            if boiler_efficiency_basis == 'LHV':
                H_combustion = fuel.LHV
                for feed in non_empty_feeds: H_combustion += feed.LHV
            elif boiler_efficiency_basis == 'HHV':
                H_combustion = fuel.HHV
                for feed in non_empty_feeds: H_combustion += feed.HHV
            else:
                raise ValueError(
                    f"invalid boiler efficiency basis {boiler_efficiency_basis}; "
                    f"valid values include 'LHV', or 'HHV'"
                )
            self.H_content = H_content = B_eff * H_combustion 
            self.H_loss_to_emissions = H_combustion - H_content
            H_electricity = H_content - H_steam # Heat available for the turbogenerator
            electricity = H_electricity * TG_eff  # Electricity produced
            self.cooling_duty = electricity - H_electricity
            Design['Work'] = work = electricity / 3600
            duty_over_mol = 39000 # kJ / mol-superheated steam 
            self.total_steam = H_content / duty_over_mol #: [float] Total steam produced by the boiler (kmol/hr)
            Design['Flow rate'] = flow_rate = self.total_steam * 18.01528
            
            if self.satisfy_system_electricity_demand:
                boiler = self.cost_items['Boiler']
                rate_boiler = boiler.kW * flow_rate / boiler.S
                return work - self.electricity_demand - rate_boiler
            else:
                return work
        
        self._excess_electricity_without_fuel = excess_electricity = calculate_excess_electricity_at_natual_gas_flow(0)
        if excess_electricity < 0:
            f = calculate_excess_electricity_at_natual_gas_flow
            lb = 0.
            fuel.imol[fuel_source] = 1
            ub = - excess_electricity * 3600 / fuel.LHV
            while f(ub) < 0.: 
                lb = ub
                ub *= 2
            flx.IQ_interpolation(f, lb, ub, xtol=1, ytol=1)
        
        if self.cooling_duty > 0.: 
            # In the event that no electricity is produced and the solver
            # solution for natural gas is slightly below the requirement for steam
            # (this would lead to a positive duty).
            self.cooling_duty = 0.
            Design['Work'] = 0.
        
        emissions.T = 298.15 # Will be updated later with the energy balance
        emissions.P = 101325
        emissions.phase = 'g'
        emissions_mol = emissions.mol
        emissions_mol[:] = fuel.mol
        for feed in non_empty_feeds: emissions_mol[:] += feed.mol
        combustion_rxns.force_reaction(emissions_mol)
        O2_consumption = -emissions.imol['O2']
        oxygen_rich_gas.reset_flow(**self.oxygen_rich_gas_composition)
        z_O2 = oxygen_rich_gas.imol['O2'] / oxygen_rich_gas.F_mol
        oxygen_rich_gas.F_mol = O2_consumption / z_O2
        emissions_mol += oxygen_rich_gas.mol
        F_emissions = emissions.F_mass
        z_CO2 = emissions.imass['CO2'] / F_emissions
        z_CO2_target = self.CO2_emissions_concentration
        if z_CO2 > z_CO2_target:
            F_emissions_new = z_CO2 * F_emissions / z_CO2_target
            dF_emissions = F_emissions_new - F_emissions
            oxygen_rich_gas.F_mass = F_mass_O2_new = oxygen_rich_gas.F_mass + dF_emissions
            emissions_mol += oxygen_rich_gas.mol * (dF_emissions / F_mass_O2_new)
        emissions.H += self.H_loss_to_emissions
        hu_cooling = bst.HeatUtility()
        hu_cooling(self.cooling_duty, steam_demand.T)
        hus_heating = bst.HeatUtility.sum_by_agent(self.steam_utilities)
        for hu in hus_heating: hu.reverse()
        self.heat_utilities = [*hus_heating, hu_cooling]
        water_index = chemicals.index('7732-18-5')
        makeup_water.mol[water_index] = blowdown_water.mol[water_index] = (
                self.total_steam * self.boiler_blowdown * 1 / (1 - self.RO_rejection)   
        )
        ash_IDs = [i.ID for i in self.chemicals if not i.formula]
        emissions_mol = emissions.mol
        SO2_produced = 0
        if 'SO2' in chemicals:
            SO2_produced += emissions.imol['SO2']
        if 'CaSO4' in chemicals:
            SO2_produced += emissions.imol['CaSO4']
            ash_IDs.append('CaSO4')
        if SO2_produced: 
            rxn, ID_lime = self._get_desulfurization_rxn_and_coreactant()
            # FGD lime scaled based on SO2 generated,	
            # 20% stoichiometric excess based on P52 of ref [1]
            rxn.force_reaction(emissions)
            lime.imol[ID_lime] = lime_mol = SO2_produced * 1.2
            emissions_mol.remove_negatives()
        else:
            lime.empty()
        # About 0.4536 kg/hr of boiler chemicals are needed per 234484 kg/hr steam produced
        chems.imol['Ash'] = boiler_chems = 1.9345e-06 * Design['Flow rate']
        ash_disposal.empty()
        ash_disposal.copy_flow(emissions, IDs=tuple(ash_IDs), remove=True)
        ash_disposal.imol['Ash'] += boiler_chems
        dry_ash = ash_disposal.F_mass
        moisture = min(emissions.imass['Water'], dry_ash * 0.3) # ~20% moisture
        ash_disposal.imass['Water'] = moisture
        emissions.imass['Water'] -= moisture
        Design['Ash disposal'] = dry_ash + moisture
        if SO2_produced:
            if ID_lime == 'Ca(OH)2': # Ca(OH)2
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
class Boiler(bst.Facility):
    """
    Create a Boiler object that will calculate steam
    generation from burning the feed.
    All capital cost correlations are based on [1]_.
    
    Parameters
    ----------
    ins : 
        * [0] Liquid/solid feed that will be burned.
        * [1] Gas feed that will be burned.
        * [2] Make-up water. 
        * [3] Natural gas/fuel to satisfy steam and electricity demand.
        * [4] Lime for flue gas desulfurization.
        * [5] Boiler chemicals.
        * [6] Air or oxygen-rich gas.
    outs : 
        * [0] Total emissions produced.
        * [1] Blowdown water.
        * [2] Ash disposal.
    boiler_efficiency : float, optional
        Fraction of heat transferred to steam. Defaults to 0.8.
    agent : UtilityAgent, optional
        Steam produced. Defaults to low pressure steam.
    other_agents = () : Iterable[UtilityAgent], optional
        Other steams produced. Defaults to all other heating agents.
    fuel_source : str, optional
        Name fuel used to satisfy steam and electricity demand. Defaults to 'CH4'.
    fuel_price : float, optional
        Price of natural gas [USD/kg]. Same as `bst.stream_utility_prices['Natural gas']`,
        which defaults to 0.218.
    ash_disposal_price : float, optional
        Price of disposing ash [USD/kg]. Same as `bst.stream_utility_prices['Ash disposal']`,
        which defaults to -0.0318.
    satisfy_system_steam_demand : bool, optional
        Whether to purchase natural gas to satisfy system steam demand
        if there is not enough heat from process feeds (i.e., inlets 0 and 1).
        If True, natural gas is purchased to satisfy system heat demand
        when there is not enough heat from the process feeds.
        In either case, if there is excess heat from the process feeds,
        excess steam will be produced.
    boiler_efficiency_basis : str, optional
        Basis of boiler efficiency. Defaults to 'LHV' (i.e., lower heating value).
        'HHV' (i.e., higher heating value) is also a valid basis. 
        
    Examples
    --------
    Create a boiler-turbogenerator system that uses sugarcane bagasse to 
    produce steam for a distillation unit and any excess steam for surplus electricity:
    
    >>> import biosteam as bst
    >>> from biorefineries import cane
    >>> chemicals = cane.create_sugarcane_chemicals()
    >>> chemicals.define_group(
    ...     name='Fiber',
    ...     IDs=['Cellulose', 'Hemicellulose', 'Lignin'],
    ...     composition=[0.4704 , 0.2775, 0.2520],
    ...     wt=True, # Composition is given as weight
    ... )
    >>> bst.settings.set_thermo(chemicals)
    >>> dilute_ethanol = bst.Stream('dilute_ethanol', Water=1390, Ethanol=590)
    >>> bagasse = bst.Stream('bagasse', Water=0.4, Fiber=0.6, total_flow=8e4, units='kg/hr')
    >>> with bst.System('sys') as sys:
    ...     D1 = bst.BinaryDistillation('D1', ins=dilute_ethanol, Lr=0.999, Hr=0.89, k=1.25, LHK=('Ethanol', 'Water'))
    ...     B = bst.Boiler('B')
    ...     B.ins[0] = bagasse
    >>> sys.simulate()
    >>> B.show()
    Boiler: B
    ins...
    [0] bagasse  
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Water          1.78e+03
                        Cellulose      139
                        Hemicellulose  101
                        Lignin         79.5
    [1] -  
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow: 0
    [2] -  
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Water  488
    [3] -  
        phase: 'g', T: 288.71 K, P: 101560 Pa
        flow: 0
    [4] -  
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow: 0
    [5] -  
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Ash  0.567
    [6] -  
        phase: 'g', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): O2  9.85e+03
                        N2  4.23e+04
    outs...
    [0] emissions  
        phase: 'g', T: 394.15 K, P: 101325 Pa
        flow (kmol/hr): Water  3.19e+03
                        CO2    1.98e+03
                        O2     7.84e+03
                        N2     4.23e+04
    [1] rejected_water_and_blowdown  
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Water  488
    [2] ash_disposal  
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Water  0.00944
                        Ash    0.567
    
    >>> B.results() # Steam is produced, so utility cost is negative
    Boiler                                                     Units         B
    Electricity           Power                                   kW  1.68e+03
                          Cost                                USD/hr       132
    Low pressure steam    Duty                                 kJ/hr -7.32e+07
                          Flow                               kmol/hr -1.89e+03
                          Cost                                USD/hr      -450
    Fuel (inlet)          Flow                                 kg/hr         0
                          Cost                                USD/hr         0
    Ash disposal (outlet) Flow                                 kg/hr     0.737
                          Cost                                USD/hr    0.0234
    Design                Flow rate                            kg/hr  2.93e+05
                          Ash disposal                         kg/hr     0.737
    Purchase cost         Baghouse bags                          USD      81.1
                          Boiler                                 USD  3.33e+07
                          Deaerator                              USD  3.58e+05
                          Amine addition pkg                     USD  4.69e+04
                          Hot process water softener system      USD  9.16e+04
    Total purchase cost                                          USD  3.38e+07
    Utility cost                                              USD/hr      -319
    
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
    ticket_name = 'B'
    network_priority = 1
    boiler_blowdown = 0.03
    # Reverse osmosis (RO) typically rejects 25% of water, but the boiler-feed water is assumed to come after RO.
    # Setting this parameter to a fraction more than zero effectively assumes that this unit accounts for reverse osmosis.
    RO_rejection = 0 
    _N_ins = 7
    _N_outs = 3
    _units = {'Flow rate': 'kg/hr',
              'Work': 'kW',
              'Ash disposal': 'kg/hr'}
    
    def __init__(self, ID='', ins=None, 
                 outs=('emissions',
                       'rejected_water_and_blowdown',
                       'ash_disposal'),
                 thermo=None, *,
                 boiler_efficiency=None,
                 side_steam=None,
                 agent=None,
                 other_agents=None,
                 fuel_price=None,
                 natural_gas_price=None,
                 ash_disposal_price=None,
                 T_emissions=None,
                 CO2_emissions_concentration=None,
                 satisfy_system_steam_demand=None,
                 boiler_efficiency_basis=None,
                 fuel_source=None,
                 oxygen_rich_gas_composition=None,
        ):
        if boiler_efficiency_basis is None: boiler_efficiency_basis = 'LHV'
        if boiler_efficiency is None: boiler_efficiency = 0.80
        if satisfy_system_steam_demand is None: satisfy_system_electricity_demand = True
        if fuel_source is None: fuel_source = 'CH4'
        if oxygen_rich_gas_composition is None: oxygen_rich_gas_composition = dict(O2=21, N2=79, phase='g', units='kg/hr')
        if CO2_emissions_concentration is None: CO2_emissions_concentration = 0.055  # Usually between 4 - 7 for biomass and natural gas (https://www.sciencedirect.com/science/article/pii/S0957582021005127)
        bst.Facility.__init__(self, ID, ins, outs, thermo)
        settings = bst.settings
        self.boiler_efficiency_basis = boiler_efficiency_basis
        self.agent = agent = agent or settings.get_heating_agent('low_pressure_steam')
        self.fuel_source = fuel_source
        self.define_utility('Fuel', self.fuel)
        self.define_utility('Ash disposal', self.ash_disposal)
        self.boiler_efficiency = boiler_efficiency
        self.steam_utilities = []
        self.steam_demand = agent.to_stream()
        self.side_steam = side_steam
        self.other_agents = [i for i in settings.heating_agents if i is not agent] if other_agents is None else other_agents
        self.CO2_emissions_concentration = CO2_emissions_concentration
        self.oxygen_rich_gas_composition = oxygen_rich_gas_composition 
        if T_emissions is not None: 
            raise ValueError('setting T_emissions is not yet implemented')
            # T_emissions should dictate the efficiency of either the boiler or the turbogenerator.
            # Note that T_emissions is tied to the energy balance.
            self.T_emissions = T_emissions
        if natural_gas_price is not None: 
            self.fuel_price = natural_gas_price
        elif fuel_price is not None: 
            self.fuel_price = fuel_price
        if ash_disposal_price is not None: self.ash_disposal_price = ash_disposal_price
        self.satisfy_system_steam_demand = satisfy_system_steam_demand
      
    def _get_desulfurization_rxn_and_coreactant(self):
        try:
            return self.desulfurization_reaction, self._ID_lime
        except:
            chemicals = self.chemicals
            CAS_lime = '1305-62-0'
            has_common_name = 'Ca(OH)2' in chemicals
            if CAS_lime in chemicals or has_common_name:
                if not has_common_name: chemicals.set_synonym(CAS_lime, 'Ca(OH)2')
                self.desulfurization_reaction = rxn = tmo.Reaction(
                    'SO2 + Ca(OH)2 + 0.5 O2 -> CaSO4 + H2O', 'SO2', 0.92, chemicals
                )
                self._ID_lime = ID = 'Ca(OH)2'
                return rxn, ID
            CAS_lime = '1305-78-8'
            has_common_name = 'CaO' in chemicals
            if CAS_lime in chemicals or has_common_name:
                if not has_common_name: chemicals.set_synonym(CAS_lime, 'CaO')
                self.desulfurization_reaction = rxn = tmo.Reaction(
                    'SO2 + CaO + 0.5 O2 -> CaSO4', 'SO2', 0.92, chemicals
                )
                self._ID_lime = ID = 'CaO'
                return rxn, ID
        raise RuntimeError(
            "lime is required for boiler, but no chemical 'CaO' or 'Ca(OH)2' "
            "available in thermodynamic property package"
        )
        
    @property
    def emissions(self):
        return self.outs[0]
    
    @property
    def blowdown_water(self):
        return self.outs[1]
    
    @property
    def makeup_water(self):
        """[Stream] Makeup water due to boiler blowdown."""
        return self.ins[2]
    
    @property
    def fuel(self):
        """[Stream] Fuel used to satisfy steam and electricity requirements."""
        return self.ins[3]
    natural_gas = fuel
    
    @property
    def air(self):
        """[Stream] Air or oxygen rich gas used to supply oxygen for combustion."""
        return self.ins[6]
    oxygen_rich_gas = air
    
    @property
    def ash_disposal(self):
        """[Stream] Ash disposal."""
        return self.outs[2]
    
    @property
    def fuel_price(self):
        """[Float] Price of natural gas, same as `bst.stream_utility_prices['Natural gas']`."""
        return bst.stream_utility_prices['Fuel']
    @fuel_price.setter
    def fuel_price(self, new_price):
        bst.stream_utility_prices['Fuel'] = new_price
    natural_gas_price = fuel_price
    
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
                        steam_utilities.append(hu)
    
    def _design(self):
        B_eff = self.boiler_efficiency
        steam_demand = self.steam_demand
        Design = self.design_results
        chemicals = self.chemicals
        self._load_utility_agents()
        mol_steam = sum([i.flow for i in self.steam_utilities])
        feed_solids, feed_gas, makeup_water, fuel, lime, chems, oxygen_rich_gas = self.ins
        oxygen_rich_gas.empty()
        if self.fuel_source == 'CH4':
            fuel.phase = 'g'
            fuel.set_property('T', 60, 'degF')
            fuel.set_property('P', 14.73, 'psi')
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
        self.combustion_reactions = combustion_rxns = chemicals.get_combustion_reactions()
        non_empty_feeds = [i for i in (feed_solids, feed_gas) if not i.isempty()]
        boiler_efficiency_basis = self.boiler_efficiency_basis
        fuel_source = self.fuel_source
        def calculate_excess_steam_at_natual_gas_flow(fuel_flow):
            if fuel_flow:
                fuel_flow = abs(fuel_flow)
                fuel.imol[fuel_source] = fuel_flow
            else:
                fuel.empty()
            if boiler_efficiency_basis == 'LHV':
                H_combustion = fuel.LHV
                for feed in non_empty_feeds: H_combustion += feed.LHV
            elif boiler_efficiency_basis == 'HHV':
                H_combustion = fuel.HHV
                for feed in non_empty_feeds: H_combustion += feed.HHV
            else:
                raise ValueError(
                    f"invalid boiler efficiency basis {boiler_efficiency_basis}; "
                    f"valid values include 'LHV', or 'HHV'"
                )
            self.H_content = H_content = B_eff * H_combustion 
            self.H_loss_to_emissions = H_combustion - H_content
            H_excess = H_content - H_steam # Heat available for the turbogenerator
            duty_over_mol = 39000 # kJ / mol-superheated steam 
            self.total_steam = H_content / duty_over_mol #: [float] Total steam produced by the boiler (kmol/hr)
            Design['Flow rate'] = flow_rate = self.total_steam * 18.01528
            return H_excess
        
        H_excess = calculate_excess_steam_at_natual_gas_flow(0)
        if H_excess < 0:
            f = calculate_excess_steam_at_natual_gas_flow
            lb = 0.
            fuel.imol[fuel_source] = 1
            ub = - H_excess / fuel.LHV
            while f(ub) < 0.: 
                lb = ub
                ub *= 2
            flx.IQ_interpolation(f, lb, ub, xtol=1, ytol=1)
        else:
            pass
            # raise NotImplementedError('BioSTEAM Boiler cannot yet produce excess steam')
        emissions.T = 298.15 # Will be updated later with the energy balance
        emissions.P = 101325
        emissions.phase = 'g'
        emissions_mol = emissions.mol
        emissions_mol[:] = fuel.mol
        for feed in non_empty_feeds: emissions_mol[:] += feed.mol
        combustion_rxns.force_reaction(emissions_mol)
        O2_consumption = -emissions.imol['O2']
        oxygen_rich_gas.reset_flow(**self.oxygen_rich_gas_composition)
        z_O2 = oxygen_rich_gas.imol['O2'] / oxygen_rich_gas.F_mol
        oxygen_rich_gas.F_mol = O2_consumption / z_O2
        emissions_mol += oxygen_rich_gas.mol
        F_emissions = emissions.F_mass
        z_CO2 = emissions.imass['CO2'] / F_emissions
        z_CO2_target = self.CO2_emissions_concentration
        if z_CO2 > z_CO2_target:
            F_emissions_new = z_CO2 * F_emissions / z_CO2_target
            dF_emissions = F_emissions_new - F_emissions
            oxygen_rich_gas.F_mass = F_mass_O2_new = oxygen_rich_gas.F_mass + dF_emissions
            emissions_mol += oxygen_rich_gas.mol * (dF_emissions / F_mass_O2_new)
        emissions.H += self.H_loss_to_emissions
        hus_heating = bst.HeatUtility.sum_by_agent(self.steam_utilities)
        for hu in hus_heating: hu.reverse()
        self.heat_utilities = hus_heating
        water_index = chemicals.index('7732-18-5')
        makeup_water.mol[water_index] = blowdown_water.mol[water_index] = (
                self.total_steam * self.boiler_blowdown * 1 / (1 - self.RO_rejection)   
        )
        ash_IDs = [i.ID for i in self.chemicals if not i.formula]
        emissions_mol = emissions.mol
        SO2_produced = 0
        if 'SO2' in chemicals:
            SO2_produced += emissions.imol['SO2']
        if 'CaSO4' in chemicals:
            SO2_produced += emissions.imol['CaSO4']
            ash_IDs.append('CaSO4')
        if SO2_produced: 
            rxn, ID_lime = self._get_desulfurization_rxn_and_coreactant()
            # FGD lime scaled based on SO2 generated,	
            # 20% stoichiometric excess based on P52 of ref [1]
            rxn.force_reaction(emissions)
            lime.imol[ID_lime] = lime_mol = SO2_produced * 1.2
            emissions_mol.remove_negatives()
        else:
            lime.empty()
        # About 0.4536 kg/hr of boiler chemicals are needed per 234484 kg/hr steam produced
        chems.imol['Ash'] = boiler_chems = 1.9345e-06 * Design['Flow rate']
        ash_disposal.empty()
        ash_disposal.copy_flow(emissions, IDs=tuple(ash_IDs), remove=True)
        ash_disposal.imol['Ash'] += boiler_chems
        dry_ash = ash_disposal.F_mass
        moisture = min(emissions.imass['Water'], dry_ash * 0.3) # ~20% moisture
        ash_disposal.imass['Water'] = moisture
        emissions.imass['Water'] -= moisture
        Design['Ash disposal'] = dry_ash + moisture
        if SO2_produced:
            if ID_lime == 'Ca(OH)2': # Ca(OH)2
                lime.imol['Water'] = 4 * lime_mol # Its a slurry
            else: # CaO
                lime.imol['Water'] = 5 * lime_mol 
        