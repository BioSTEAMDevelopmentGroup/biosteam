# -*- coding: utf-8 -*-
"""
"""
import biosteam as bst
from .._unit import design_input
from math import sqrt, pi
import numpy as np

__all__ = ('FluidizedCatalyticCracking',)

class FluidizedCatalyticCracking(bst.Unit):
    """
    Create a reactor for catalytic cracking.
    
    Parameters
    ----------
    ins :
        * [0] liquid feed (upressurized)
        * [1] make-up catalyst
        * [2] air (unpressurized)
        * [3] steam (superheated)
    outs :
        * [0] liquid product
        * [1] spent catalyst
        * [2] flue gas
    
    Other Parameters
    ----------------
    ...
    
    Examples
    --------
    ...
    
    References
    ----------
    ...
    
    """
    _N_ins = 4
    _N_outs = 3
    auxiliary_unit_names = (
        'pump', 'feed_preheater', 'air_compressor', 'condenser', 
        'riser', 'stripper', 'regenerator', 'reactor',
    )
    vessel_material: str = design_input(default='Stainless steel 316')
    feed_pressure = design_input(
        'Feed pressure', 'Pa', lb=0, default=2 * 101325,
        doc="Feed pressure [Pa]. Can operate at a low pressure, but must "
            "accomodate pressure drop. Defaults to 2 atm.",
    )
    feed_vapor_fraction = design_input(
        'Feed vapor fraction', lb=0, ub=1, default=0,
    )
    catalyst_to_feed_ratio = design_input(
        'Catalyst to feed ratio', lb=0, default=5,
        doc="Catalyst to feed ratio by weight. Defaults to 5.",
    )
    spent_catalyst_to_feed_ratio = design_input(
        'Spent catalyst to feed ratio', lb=0,
        default=1e-6, # TODO: Find value; value is abitrary for now
        doc="Spent catalyst to feed ratio by weight. Defaults to 1e-6.",
    )
    product_loss = design_input(
        'Product loss', lb=0, ub=1, default=0.5e-2,
        doc="Fraction of product lost (combusted in the regenerator). "
            "Defaults to 0.005.",
    )
    riser_product_residence_time = design_input(
        'Riser product residence time', 'hr', lb=0,
        default=2.1 / 3600, # 1.8 to 2.4 seconds preferred with VGO
    )
    riser_length_to_diameter = design_input(
        'Riser length to diameter ratio', lb=0, default=20, # Minimally preferred
    )
    reactor_vessel_dissengagement_height = design_input(
        'Reactor vessel disengagement height', 'm', lb=0, default=4.572, # 15 ft minimally
    )
    reactor_vessel_exit_velocity = design_input(
        'Reactor vessel exit velocity', 'm/hr', lb=0, default=2194.56,
    )
    stripper_steam_rate = design_input(
        'Stripper steam rate', lb=0, default=0.25, # 0.2-0.3 wt % steam per catalyst
        doc="Steam to catalyst ratio in the stripper by weight. Defaults to 0.25.",
    )
    stripper_catalyst_residence_time = design_input(
        'Stripper catalyst residence time', 'hr', lb=0,
        default=75 / 3600, # typically 60 to 90 seconds
    )
    stripper_length_to_diameter = design_input(
        'Stripper length to diameter ratio', lb=0, default=1,
    )
    CO2_concentration = design_input(
        'CO2 concentration', lb=0, ub=1, default=0.055,
        doc="Maximum CO2 mass fraction in the flue gas; air flow is increased "
            "to meet this concentration. Defaults to 0.055.",
    )
    regenerator_catalyst_residence_time = design_input(
        'Regenerator catalyst residence time', 'hr', lb=0,
        default=7.5 / 60, # 5 to 10 min for 1-stage regenerator
    )
    regenerator_length_to_diameter = design_input(
        'Regenerator length to diameter ratio', lb=0, default=0.5, # Typically 0.3 to 0.7
    )
    regenerator_pressure = design_input(
        'Regenerator pressure', 'Pa', lb=0, default=282685, # 40 psig
    )

    def _init(self, reaction, bulk_reactant=None):
        self.reaction = reaction
        self.bulk_reactant = bulk_reactant

    @property
    def feed(self): return self.ins[0]
    
    @property
    def makeup_catalyst(self): return self.ins[1]
    
    @property
    def air(self): return self.ins[2]
    
    @property
    def steam(self): return self.ins[3]
    
    @property
    def product(self): return self.outs[0]
    
    @property
    def spent_catalyst(self): return self.outs[1]
    
    @property
    def flue_gas(self): return self.outs[2]
    
    def _setup(self):
        pump = self.auxiliary(
            'pump', bst.Pump, ins=self.feed, P=self.feed_pressure,
        )
        self.auxiliary(
            'feed_preheater', bst.HXutility, ins=pump-0, V=self.feed_vapor_fraction, rigorous=True,
        )
        self.auxiliary(
            'air_compressor', bst.IsentropicCompressor, ins=self.air, P=self.regenerator_pressure,
        )

    def _estimate_reactor_pressure_drop(self):
        return 40530.0

    def _estimate_regenerator_pressure_drop(self):
        return 40530.0

    def _run(self):
        feed, fresh_catalyst, air, steam = self.ins
        product, discarded_catalyst, flue_gas = self.outs
        self.pump.run()
        self.feed_preheater.run()
        F_mass_feed = feed.F_mass 
        fresh_catalyst.imass['Catalyst'] = discarded_catalyst.imass['Catalyst'] = F_mass_feed * self.spent_catalyst_to_feed_ratio
        self.catalyst_recirculation = catalyst_recirculation = F_mass_feed * self.catalyst_to_feed_ratio # kg / h
        steam.copy_like(bst.settings.get_heating_agent('low_pressure_steam'))
        steam.imass['Water'] = catalyst_recirculation * self.stripper_steam_rate
        product.P = self.feed_pressure - self._estimate_reactor_pressure_drop()
        flue_gas.P = discarded_catalyst.P = self.regenerator_pressure - self._estimate_regenerator_pressure_drop()
        product.mol = feed.mol + steam.mol
        if self.bulk_reactant:
            if self.reaction.X != 1:
                raise RuntimeError('reaction extent must be 100% with bulk reactants')
            model_reactant, reactants = self.bulk_reactant
            F_reactant_mass = product.imass[reactants].sum()
            product.imol[reactants] = 0
            product.imass[model_reactant] = F_reactant_mass
        product.phase = 'g'
        flue_gas.phase = 'g'
        self.reaction(product)
        # Product loss will dictate temperature of recirculated catalyst
        product.split_to(flue_gas, product, self.product_loss, energy_balance=False)
        product_loss = flue_gas.mol.copy()
        combustion = self.chemicals.get_combustion_reactions()
        combustion.force_reaction(flue_gas)
        O2 = -flue_gas.imass['O2']
        N2 = 0.79 / 0.21 * O2
        air.imass['O2', 'N2'] = [O2, N2]
        flue_gas.mol += air.mol
        F_emissions = flue_gas.F_mass
        z_CO2 = flue_gas.imass['CO2'] / F_emissions
        z_CO2_target = self.CO2_concentration
        if z_CO2 > z_CO2_target:
            F_emissions_new = z_CO2 * F_emissions / z_CO2_target
            dF_emissions = F_emissions_new - F_emissions
            air.F_mass = F_mass_O2_new = air.F_mass + dF_emissions
            flue_gas.mol += air.mol * (dF_emissions / F_mass_O2_new)
        self.air_compressor.run()
        regenerated_catalyst = bst.Stream(
            None, 
            Catalyst=catalyst_recirculation,
            units='kg/hr'
        )
        spent_catalyst = bst.Stream(
            None, 
            Catalyst=catalyst_recirculation, 
            units='kg/hr'
        )
        spent_catalyst.mol += product_loss
        compressed_air = self.air_compressor-0
        heated_feed = self.feed_preheater-0
        regenerator_outlets = [flue_gas, regenerated_catalyst, discarded_catalyst]
        regenerator_inlets = [compressed_air, spent_catalyst, fresh_catalyst]
        reactor_outlets = [product, spent_catalyst]
        reactor_inlets = [heated_feed, steam, regenerated_catalyst]
        feed_temperature = self.feed_preheater.outs[0].T
        for i in regenerator_outlets + reactor_outlets: i.T = feed_temperature # Initial guess
        dTs = np.ones(2)
        while np.abs(dTs).sum() > 0.1:
            # Energy balance correction equations
            # dT_regenerator * (C_flue_gas + C_catalyst_regen + C_discarded_catalyst) - dT_reactor * C_catalyst_spent = DeltaH_regen
            # dT_reactor * (C_catalyst_spent + C_product) - dT_regenerator * C_catalyst_regen = DeltaH_reactor
            A = np.array(
                [[sum([i.C for i in regenerator_outlets]), -spent_catalyst.C],
                 [-regenerated_catalyst.C, sum([i.C for i in reactor_outlets])]]
            )
            b = np.array(
                [sum([i.Hnet for i in regenerator_inlets]) - sum([i.Hnet for i in regenerator_outlets]),
                 sum([i.Hnet for i in reactor_inlets]) - sum([i.Hnet for i in reactor_outlets])]
            )
            dT_regenerator, dT_reactor = dTs = np.linalg.solve(A, b)
            for i in regenerator_outlets: i.T += dT_regenerator
            for i in reactor_outlets: i.T += dT_reactor

    def _design(self):
        # volume = L * pi * D**2 / 4
        # volume = D * L2D * pi * D**2 / 4
        # volume = D**3 * L2D * pi / 4
        # D = (volume * 4 / (pi * L2D)) ** (1 / 3)
        
        # Riser
        self.riser_volume = riser_volume = self.riser_product_residence_time * self.product.F_vol # m3
        L2D = self.riser_length_to_diameter
        self.riser_diameter = riser_diameter = (riser_volume * 4 / (pi * L2D)) ** (1 / 3)
        self.riser_length = riser_length = riser_diameter * L2D
        self.riser = bst.AuxiliaryPressureVessel(
            self.feed_pressure, riser_diameter, riser_length,
            pressure_units='Pa', length_units='m', 
            orientation='Vertical', material=self.vessel_material, 
        )
        
        # Reactor
        self.reactor_length = reactor_length = self.reactor_vessel_dissengagement_height
        self.reactor_diameter = reactor_diameter = sqrt(
            4 * self.product.F_vol / self.reactor_vessel_exit_velocity / pi
        )
        self.reactor = bst.AuxiliaryPressureVessel(
            self.feed_pressure, reactor_diameter, reactor_length,
            pressure_units='Pa', length_units='m', 
            orientation='Vertical', material=self.vessel_material, 
        )
        
        # Stripper
        catalyst_recirculation_volume = self.catalyst_recirculation / self.chemicals.Catalyst.rho(self.product.T, self.product.P)
        self.stripper_volume = stripper_volume = (
            self.stripper_catalyst_residence_time
            * catalyst_recirculation_volume
        )
        L2D = self.stripper_length_to_diameter
        self.stripper_diameter = stripper_diameter = (stripper_volume * 4 / (pi * L2D)) ** (1 / 3)
        self.stripper_length = stripper_length = stripper_diameter * L2D
        self.stripper = bst.AuxiliaryPressureVessel(
            self.feed_pressure, stripper_diameter, stripper_length,
            pressure_units='Pa', length_units='m', 
            orientation='Vertical', material=self.vessel_material, 
        )
        
        # Regenerator
        self.regenerator_volume = regenerator_volume = (
            self.regenerator_catalyst_residence_time
            * catalyst_recirculation_volume
        )
        L2D = self.regenerator_length_to_diameter
        self.regenerator_diameter = regenerator_diameter = (regenerator_volume * 4 / (pi * L2D)) ** (1 / 3)
        self.regenerator_length = regenerator_length = regenerator_diameter * L2D
        self.regenerator = bst.AuxiliaryPressureVessel(
            self.regenerator_pressure, regenerator_diameter, regenerator_length,
            pressure_units='Pa', length_units='m', 
            orientation='Vertical', material=self.vessel_material, 
        )
        
        for i in self.auxiliary_units: i._design()
    
    def _cost(self):
        for i in self.auxiliary_units: i._cost()


FCC = FluidizedCatalyticCracking