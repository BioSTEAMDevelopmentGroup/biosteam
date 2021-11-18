# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
This module contains unit operations for drying solids.

.. contents:: :local:
    
Unit operations
---------------
.. autoclass:: biosteam.units.drying.DrumDryer


"""
import biosteam as bst
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

__all__ = ('DrumDryer', 'ThermalOxidizer')

# TODO: The drum dryer is carbon steel. Add material factors later
@cost('Peripheral drum area', CE=CEPCI_by_year[2007], ub=7854.0, BM=2.06,
      S=1235.35, units='m2', n=0.6, cost=0.52 * 2268000., kW=938.866)
class DrumDryer(Unit):
    """
    Create a drum dryer that dries solids by passing hot air 
    (heated by burning natural gas).
    
    Parameters
    ----------
    ins : stream sequence
        [0] Wet solids.
        [1] Air.
        [2] Natural gas.
    outs : stream sequence
        [0] Dried solids
        [1] Hot air
        [2] Emissions
    split : dict[str, float]
        Component splits to hot air (stream [1]).
    R : float, optional
        Flow of hot air over evaporation. Defaults to 1.4 wt gas / wt evap.
    H : float, optional
        Specific evaporation rate [kg/hr/m3]. Defaults to 20. 
    length_to_diameter : float, optional
        Note that the drum is horizontal. Defaults to 25.
    T : float, optional
        Operating temperature [K]. Defaults to 343.15.
    moisture_content : float
        Moisutre content of solids [wt / wt]. Defaults to 0.10.
        
    Notes
    -----
    The flow rate for air in the inlet is varied to meet the `R` specification
    (i.e. flow of hot air over flow rate evaporated). The flow rate of inlet natural
    gas is also altered to meet the heat demand.
    
    The default parameter values are based on heuristics for drying 
    dried distillers grains with solubles (DDGS).
    
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
        return self._isplit._data
    
    @property
    def natural_gas(self):
        """[Stream] Natural gas to satisfy steam and electricity requirements."""
        return self.ins[2]
    
    def __init__(self, ID="", ins=None, outs=(), thermo=None, *,
                 split, R=1.4, H=20., length_to_diameter=25, T=343.15,
                 moisture_content=0.15, utility_agent='Natural gas'):
        super().__init__(ID, ins, outs, thermo)
        self._isplit = self.chemicals.isplit(split)
        self.define_utility('Natural gas', self.natural_gas)
        self.T = T
        self.R = R
        self.H = H
        self.length_to_diameter = length_to_diameter
        self.moisture_content = moisture_content
        self.utility_agent = utility_agent
        
    @property
    def utility_agent(self):
        return self._utility_agent
    
    @utility_agent.setter
    def utility_agent(self, utility_agent):
        if utility_agent == 'Natural gas':
            pass
        elif utility_agent == 'Steam':
            self.heat_utilities = (bst.HeatUtility(),)
        else:
            raise ValueError(f"utility agent must be either 'Steam' or 'Natural gas'; not '{utility_agent}'")
        self._utility_agent = utility_agent
        
    def _run(self):
        wet_solids, air, natural_gas = self.ins
        dry_solids, hot_air, emissions = self.outs
        wet_solids.split_to(hot_air, dry_solids, self.split)
        sep.adjust_moisture_content(dry_solids, hot_air, self.moisture_content)
        emissions.phase = air.phase = natural_gas.phase = hot_air.phase = 'g'
        design_results = self.design_results
        design_results['Evaporation'] = evaporation = hot_air.F_mass
        air.imass['N2', 'O2'] = np.array([0.78, 0.32]) * self.R * evaporation
        hot_air.mol += air.mol
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
                duty = (dry_solids.H + hot_air.H + emissions.H) - (wet_solids.H + air.H + natural_gas.H)
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
            self.heat_utilities[0](self.H_out - self.H_in, self.T)
        
        
class ThermalOxidizer(Unit):
    """
    Create a ThermalOxidizer that burns any remaining combustibles.
    
    Parameters
    ----------
    ins : stream sequence
        [0] Feed gas
        [1] Air
        [2] Natural gas
    outs : stream
        Emissions.
    tau : float, optional
        Residence time [hr]. Defaults to 0.00014 (0.5 seconds).
    duty_per_kg : float, optional
        Duty per kg of feed. Defaults to 105858 kJ / kg.
    V_wf : float, optional
        Fraction of working volume. Defaults to 0.95.
    
    Notes
    -----
    Adiabatic operation is assummed. Simulation and cost is based on [1]_.
    
    References
    ----------
    .. [1] Kwiatkowski, J. R.; McAloon, A. J.; Taylor, F.; Johnston, D. B. 
        Modeling the Process and Costs of Fuel Ethanol Production by the Corn 
        Dry-Grind Process. Industrial Crops and Products 2006, 23 (3), 288–296.
        https://doi.org/10.1016/j.indcrop.2005.08.004.

    """
    _N_ins = 3
    _N_outs = 1
    max_volume = 20. # m3
    _F_BM_default = {'Vessels': 2.06} # Assume same as dryer
    
    @property
    def natural_gas(self):
        """[Stream] Natural gas to satisfy steam and electricity requirements."""
        return self.ins[2]
    
    def __init__(self, *args, tau=0.00014, duty_per_kg=61., V_wf=0.95, 
                 **kwargs):
        Unit.__init__(self, *args, **kwargs)
        self.define_utility('Natural gas', self.natural_gas)
        self.tau = tau
        self.duty_per_kg = duty_per_kg
        self.V_wf = V_wf

    def _run(self):
        feed, air, ng = self.ins
        ng.imol['CH4'] = self.duty_per_kg * feed.F_mass / self.chemicals.CH4.LHV
        emissions, = self.outs
        emissions.mix_from([feed, ng])
        combustion_rxns = self.chemicals.get_combustion_reactions()
        combustion_rxns.force_reaction(emissions)
        O2 = max(-emissions.imol['O2'], 0.)
        air.imol['N2', 'O2'] = [0.78/0.32 * O2, O2]
        emissions.mix_from(self.ins)
        combustion_rxns.adiabatic_reaction(emissions)
        
    def _design(self):
        design_results = self.design_results
        volume = self.tau * self.outs[0].F_vol / self.V_wf
        V_max = self.max_volume
        design_results['Number of vessels'] = N = np.ceil(volume / V_max)
        design_results['Vessel volume'] = volume / N
        design_results['Total volume'] = volume
        
    def _cost(self):
        design_results = self.design_results
        total_volume = design_results['Total volume']
        N = design_results['Number of vessels']
        vessel_volume = design_results['Vessel volume']
        C = self.baseline_purchase_costs
        C['Vessels'] = N * 918300. * (vessel_volume / 13.18)**0.6