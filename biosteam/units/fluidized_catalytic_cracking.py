# -*- coding: utf-8 -*-
"""
"""
import biosteam as bst
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

    def _init(
            self,
            reaction,
            bulk_reactant=None,
            vessel_material=None,
            feed_pressure=None,
            feed_vapor_fraction=None,
            catalyst_to_feed_ratio=None,
            spent_catalyst_to_feed_ratio=None,
            product_loss=None,
            riser_product_residence_time=None,
            # riser_entrance_superficial_velocity=None,
            riser_length_to_diameter=None,
            reactor_vessel_dissengagement_height=None, 
            reactor_vessel_exit_velocity=None,
            stripper_steam_rate=None,
            stripper_catalyst_residence_time=None,
            stripper_length_to_diameter=None,
            # regenerator_exit_velocity=None,
            CO2_concentration=None,
            regenerator_catalyst_residence_time=None,
            regenerator_length_to_diameter=None,
            regenerator_pressure=None,
        ):
        if vessel_material is None:
            vessel_material = 'Stainless steel 316'
        if feed_pressure is None:
            feed_pressure = 2 * 101325 # Can operate at a low pressure, but must accomodate pressure drop
        if feed_vapor_fraction is None:
            feed_vapor_fraction = 0
        if catalyst_to_feed_ratio is None:
            catalyst_to_feed_ratio = 5 # by weight
        if spent_catalyst_to_feed_ratio is None:
            # TODO: Find value
            # Mitias recommends (rule of thumb) 
            spent_catalyst_to_feed_ratio = 1e-6 # by weight; Value is abitrary for now 
        if product_loss is None:
            product_loss = 0.5e-2
        if riser_product_residence_time is None:
            riser_product_residence_time = 2.1 / 3600  # hr; 1.8 to 2.4 seconds preferred with VGO 
        # if riser_entrance_superficial_velocity is None:
        #     riser_entrance_superficial_velocity = 21945.6 # m / h; minimally 15 to 20 ft / s
        if riser_length_to_diameter is None:
            riser_length_to_diameter = 20 # Minimally preferred
        if reactor_vessel_dissengagement_height is None:
            reactor_vessel_dissengagement_height = 4.572 # m; 15 ft minimally
        if reactor_vessel_exit_velocity is None:
            reactor_vessel_exit_velocity = 2194.56 # m / s; 
        if stripper_steam_rate is None:
            stripper_steam_rate = 0.25 # wt %; 0.2-0.3 wt % steam per catalyst
        if stripper_catalyst_residence_time is None:
            stripper_catalyst_residence_time = 75 / 3600 # hr; typically 60 to 90 seconds
        if stripper_length_to_diameter is None:
            stripper_length_to_diameter = 1
        if CO2_concentration is None: 
            CO2_concentration = 0.055
        # if regenerator_exit_velocity is None:
        #     regenerator_exit_velocity = 2743.2 # m / h; between 2 to 3 ft / s
        if regenerator_catalyst_residence_time is None:
            regenerator_catalyst_residence_time = 7.5 / 60 # h; 5 to 10 min for 1-stage regenerator
        if regenerator_length_to_diameter is None:
            regenerator_length_to_diameter = 0.5 # Typically 0.3 to 0.7 
        if regenerator_pressure is None:
            regenerator_pressure = 282685 # 40 psig
        self.reaction = reaction
        self.bulk_reactant = bulk_reactant
        self.vessel_material = vessel_material
        self.feed_pressure = feed_pressure
        self.catalyst_to_feed_ratio = catalyst_to_feed_ratio
        self.spent_catalyst_to_feed_ratio = spent_catalyst_to_feed_ratio
        self.product_loss = product_loss
        self.riser_product_residence_time = riser_product_residence_time
        # self.riser_entrance_superficial_velocity = riser_entrance_superficial_velocity
        self.riser_length_to_diameter = riser_length_to_diameter
        self.reactor_vessel_dissengagement_height = reactor_vessel_dissengagement_height
        self.reactor_vessel_exit_velocity = reactor_vessel_exit_velocity
        self.stripper_steam_rate = stripper_steam_rate
        self.stripper_catalyst_residence_time = stripper_catalyst_residence_time
        self.stripper_length_to_diameter = stripper_length_to_diameter
        self.CO2_concentration = CO2_concentration
        # self.regenerator_exit_velocity = regenerator_exit_velocity
        self.regenerator_catalyst_residence_time = regenerator_catalyst_residence_time
        self.regenerator_length_to_diameter = regenerator_length_to_diameter
        self.regenerator_pressure = regenerator_pressure
        self.feed_vapor_fraction = feed_vapor_fraction
        
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