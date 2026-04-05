# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
This module contains unit operations for drying solids.

.. contents:: :local:
    
.. autoclass:: biosteam.units.drying.DrumDryer

References
----------
.. [1] Kwiatkowski, J. R.; McAloon, A. J.; Taylor, F.; Johnston, D. B. 
    Modeling the Process and Costs of Fuel Ethanol Production by the Corn 
    Dry-Grind Process. Industrial Crops and Products 2006, 23 (3), 288–296.
    https://doi.org/10.1016/j.indcrop.2005.08.004.

"""
import flexsolve as flx
import numpy as np
from thermosteam import separations as sep
from .design_tools import (
    CEPCI_by_year,
    cylinder_diameter_from_volume, 
    cylinder_area,
)
from .decorators import cost
from .._unit import Unit
from math import exp, log
from thermosteam import separations
from warnings import warn
from ..exceptions import DesignWarning

__all__ = ('SprayDryer', 'DrumDryer', 'ThermalOxidizer')

@cost('Evaporation rate', CE=567, units='lb/hr', # lb=30, ub=3000, 
      BM=2.06, f=lambda W: exp(8.5133 + 0.9847*(logW:=log(W)) - 0.0561 * logW * logW)
)
class SprayDryer(Unit):
    _units = {'Evaporation rate': 'lb/hr'}
    _N_ins = 1
    _N_outs = 2
    
    def _init(self, moisture_content=0.90):
        self.moisture_content = moisture_content
        
    def _run(self):
        feed = self.ins[0]
        water, solids = self.outs
        solids.copy_like(feed)
        water.copy_flow(solids, 'Water', remove=True)
        separations.adjust_moisture_content(solids, water, self.moisture_content)
        water.phase = 'g'
        
    def _design(self):
        water, solids = self.outs
        self.design_results['Evaporation rate'] = water.get_total_flow('lb/hr')


# TODO: The drum dryer is carbon steel. Add material factors later
@cost('Peripheral drum area', CE=CEPCI_by_year[2007], ub=7854.0, BM=2.06,
      S=1235.35, units='m2', n=0.6, cost=0.52 * 2268000., kW=938.866)
class DrumDryer(Unit):
    """
    Create a drum dryer that dries solids by passing hot air 
    (heated by burning natural gas).
    
    Parameters
    ----------
    ins : 
        * [0] Wet solids.
        * [1] Dry gas.
        * [2] Natural gas.
    outs : 
        * [0] Dried solids
        * [1] Hot gas
        * [2] Emissions
    split : dict[str, float]
        Component splits to hot gas (stream [1]).
    RH : float, optional
        Relative humidity of hot gas as a fraction. Defaults to 0.80.
    H : float, optional
        Specific evaporation rate [kg/hr/m3]. Defaults to 20. 
    length_to_diameter : float, optional
        Note that the drum is horizontal. Defaults to 25.
    T : float, optional
        Operating temperature [K]. Defaults to 343.15.
    moisture_content : float
        Moisture content of solids [wt / wt]. Defaults to 0.10.
        
    Notes
    -----
    The flow rate for gas in the inlet is calculated to meet the `RH` specification
    (i.e. relative humidity of hot gas depending on moisture_ID evaporated). The flow 
    rate of inlet natural gas is also altered to meet the heat demand.
    
    The default parameter values are based on heuristics for drying 
    dried distillers grains with solubles (DDGS).
    
    Examples
    --------
    >>> import biosteam as bst
    >>> from biorefineries import corn as c
    >>> bst.settings.set_thermo(c.create_chemicals())
    >>> feed = bst.Stream('feed', phase='l', T=352.33, P=101325,
    ...     Water=0.6749, Ethanol=5.041e-06, Ash=0.01978, Yeast=0.008452, 
    ...     CaO=0.0001446, TriOlein=0.02702, H2SO4=0.001205, Fiber=0.1508, 
    ...     SolubleProtein=0.04805, InsolubleProtein=0.06967, 
    ...     total_flow=32720, units='kg/hr',
    ... )
    >>> dryer = bst.DrumDryer('D610', 
    ...     (feed, 'dryer_air', 'natural_gas'), 
    ...     ('dryed_solids', 'hot_air', 'emissions'),
    ...     moisture_content=0.10, split=dict(Ethanol=1.0)
    ... )
    >>> dryer.simulate()
    >>> dryer.show('cwt100')
    DrumDryer: D610
    ins...
    [0] feed
        phase: 'l', T: 352.33 K, P: 101325 Pa
        composition (%): Water             67.5
                         Ethanol           0.000504
                         Ash               1.98
                         Yeast             0.845
                         CaO               0.0145
                         TriOlein          2.7
                         H2SO4             0.12
                         Fiber             15.1
                         SolubleProtein    4.8
                         InsolubleProtein  6.97
                         ----------------  3.27e+04 kg/hr
    [1] dryer_air
        phase: 'g', T: 298.15 K, P: 1.01325e+06 Pa
        composition (%): O2  21
                         N2  79
                         --  1.32e+06 kg/hr
    [2] natural_gas
        phase: 'g', T: 298.15 K, P: 101325 Pa
        composition (%): CH4  100
                         ---  2.45e+03 kg/hr
    outs...
    [0] dryed_solids
        phase: 'l', T: 343.15 K, P: 101325 Pa
        composition (%): Water             10
                         Ash               5.48
                         Yeast             2.34
                         CaO               0.04
                         TriOlein          7.48
                         H2SO4             0.334
                         Fiber             41.7
                         SolubleProtein    13.3
                         InsolubleProtein  19.3
                         ----------------  1.18e+04 kg/hr
    [1] hot_air
        phase: 'g', T: 343.15 K, P: 1.01325e+06 Pa
        composition (%): Water    1.56
                         Ethanol  1.23e-05
                         O2       20.7
                         N2       77.8
                         -------  1.34e+06 kg/hr
    [2] emissions
        phase: 'g', T: 373.15 K, P: 101325 Pa
        composition (%): Water  45
                         CO2    55
                         -----  1.22e+04 kg/hr
                        
    >>> dryer.results()
    Drum dryer                                 Units     D610
    Electricity         Power                     kW      845
                        Cost                  USD/hr       66
    Natural gas (inlet) Flow                   kg/hr 2.45e+03
                        Cost                  USD/hr      534
    Design              Evaporation            kg/hr 2.09e+04
                        Volume                       1.05e+03
                        Diameter                   m     3.76
                        Length                             94
                        Peripheral drum area      m2 1.11e+03
    Purchase cost       Drum dryer               USD  1.2e+06
    Total purchase cost                          USD  1.2e+06
    Utility cost                              USD/hr      600
    
    """
    # auxiliary_unit_names = ('heat_exchanger',)
    _units = {'Evaporation': 'kg/hr',
              'Peripheral drum area': 'm2',
              'Diameter': 'm'}
    _N_ins = 3
    _N_outs = 3
    
    @property
    def isplit(self):
        """[ChemicalIndexer] Componentwise split of feed to 0th outlet stream."""
        return self._isplit
    @property
    def split(self):
        """[Array] Componentwise split of feed to 0th outlet stream."""
        return self._isplit.data
    
    @property
    def natural_gas(self):
        """[Stream] Natural gas to satisfy steam and electricity requirements."""
        return self.ins[2]
    
    def _init(self, split, RH=0.80, H=20., length_to_diameter=25, T=343.15, P=10*101325,
              moisture_content=0.15, utility_agent='Natural gas', gas_composition=None,
              moisture_ID=None):
        self._isplit = self.chemicals.isplit(split)
        self.define_utility('Natural gas', self.natural_gas)
        self.P = P
        self.T = T
        self.RH = RH
        self.H = H
        self.gas_composition = [('N2', 0.79), ('O2', 0.21)] if gas_composition is None else gas_composition
        self.length_to_diameter = length_to_diameter
        self.moisture_content = moisture_content
        self.utility_agent = utility_agent
        self.moisture_ID = 'Water' if moisture_ID is None else moisture_ID
        
    @property
    def utility_agent(self):
        return self._utility_agent
    
    @utility_agent.setter
    def utility_agent(self, utility_agent):
        if utility_agent not in ('Natural gas', 'Steam'):
            raise ValueError(f"utility agent must be either 'Steam' or 'Natural gas'; not '{utility_agent}'")
        self._utility_agent = utility_agent

    def _get_moisture_vapor_pressure(self, T):
        chemical = self.thermo.chemicals[self.moisture_ID]
        return chemical.Psat(T)

    def _convert_air_mol_to_mass(self, n_air, gas_composition):
        mol_weight_air = 0.
        for chem, x in gas_composition:
            chem_mol_weight = self.thermo.chemicals[chem].MW
            mol_weight_air += x / chem_mol_weight
        return n_air / mol_weight_air

    def _run(self):
        wet_solids, dry_gas, natural_gas = self.ins
        dry_solids, hot_air, emissions = self.outs
        wet_solids.split_to(hot_air, dry_solids, self.split)
        sep.adjust_moisture_content(dry_solids, hot_air, self.moisture_content, self.moisture_ID)
        hot_air.P = dry_gas.P = self.P
        emissions.phase = dry_gas.phase = natural_gas.phase = hot_air.phase = 'g'
        design_results = self.design_results
        design_results['Evaporation'] = hot_air.F_mass

        # Calculate n_moisture_ID and n_evap_compounds
        n_moisture = hot_air.imol[self.moisture_ID]

        # Calculate y_moisture_ID
        Psat = self._get_moisture_vapor_pressure(self.T)
        P = self.P
        if Psat > P: 
            warn(
                f'saturated pressure of {self.moisture_ID} ({int(Psat)} Pa) '
                f'is greater than operating pressure ({int(P)} Pa)',
                category=DesignWarning,
            )
            Psat = P
        y_moisture = self.RH * Psat / P
        
        # Calculate total gas flow (molar basis)
        n_other = n_moisture * (1 - y_moisture) / y_moisture
        dry_gas.reset_flow(**dict(self.gas_composition), total_flow=n_other) # In kmol / hr
        hot_air.mol += dry_gas.mol
        dry_solids.T = hot_air.T = self.T
        emissions.T = self.T + 30.
        natural_gas.empty()
        emissions.empty()
        if self.utility_agent == 'Natural gas':
            LHV = self.chemicals.CH4.LHV
            def f(CH4):
                CO2 = CH4    
                H2O = 2. * CH4
                natural_gas.imol['CH4'] = CH4
                emissions.imol['CO2', 'H2O'] = [CO2, H2O]    
                duty = (dry_solids.H + hot_air.H + emissions.H) - (wet_solids.H + dry_gas.H + natural_gas.H)
                CH4 = duty / LHV
                return CH4
            flx.wegstein(f, 0., 1e-3)
        
    def _design(self):
        length_to_diameter = self.length_to_diameter
        design_results = self.design_results
        design_results['Volume'] = volume = design_results['Evaporation'] / self.H 
        design_results['Diameter'] = diameter = cylinder_diameter_from_volume(volume, length_to_diameter)
        design_results['Length'] = length = diameter * length_to_diameter
        design_results['Peripheral drum area'] = cylinder_area(diameter, length)
        if self.utility_agent == 'Steam':
            self.add_heat_utility(self.H_out - self.H_in, self.T)
        
        
class ThermalOxidizer(Unit):
    """
    Create a ThermalOxidizer that burns any remaining combustibles.
    
    Parameters
    ----------
    ins : 
        [0] Feed gas
        [1] Air
        [2] Natural gas
    outs : 
        Emissions.
    tau : float, optional
        Residence time [hr]. Defaults to 0.00014 (0.5 seconds).
    duty_per_kg : float, optional
        Duty per kg of feed. Defaults to 105858 kJ / kg.
    V_wf : float, optional
        Fraction of working volume. Defaults to 0.95.
    
    Notes
    -----
    Adiabatic operation is assumed. Simulation and cost is based on [1]_.

    """
    _N_ins = 3
    _N_outs = 1
    max_volume = 20. # m3
    _F_BM_default = {'Vessels': 2.06} # Assume same as dryer
    
    @property
    def natural_gas(self):
        """[Stream] Natural gas to satisfy steam and electricity requirements."""
        return self.ins[2]
    
    def _init(self, tau=0.00014, duty_per_kg=61.):
        self.define_utility('Natural gas', self.natural_gas)
        self.tau = tau
        self.duty_per_kg = duty_per_kg

    def _run(self):
        feed, air, ng = self.ins
        ng.imol['CH4'] = self.duty_per_kg * feed.F_mass / self.chemicals.CH4.LHV
        ng.phase = 'g'
        air.phase = 'g'
        emissions, = self.outs
        ng.P = air.P = emissions.P = feed.P
        emissions.phase = 'g'
        ng_burned = ng.copy()
        combustion_rxns = self.chemicals.get_combustion_reactions()
        # Enough oxygen must be present in air to burn natural gas
        combustion_rxns.force_reaction(ng_burned)
        O2 = max(-ng_burned.imol['O2'], 0.)
        air.imol['N2', 'O2'] = [0.79/0.21 * O2, O2]
        # Enough oxygen must be present in air to burn feed as well
        emissions.mix_from(self.ins)
        dummy_emissions = emissions.copy()
        combustion_rxns.force_reaction(dummy_emissions)
        O2 = max(-dummy_emissions.imol['O2'], 0.) # Missing oxygen
        air.imol['N2', 'O2'] += [0.79/0.21 * O2, O2]
        emissions.mix_from(self.ins)
        # Account for temperature raise
        combustion_rxns.adiabatic_reaction(emissions)
        
    def _design(self):
        design_results = self.design_results
        volume = self.tau * self.ins[0].F_vol
        V_max = self.max_volume
        design_results['Number of vessels'] = N = np.ceil(volume / V_max)
        design_results['Vessel volume'] = volume / N
        design_results['Total volume'] = volume
        
    def _cost(self):
        design_results = self.design_results
        N = design_results['Number of vessels']
        vessel_volume = design_results['Vessel volume']
        C = self.baseline_purchase_costs
        C['Vessels'] = N * 918300. * (vessel_volume / 13.18)**0.6