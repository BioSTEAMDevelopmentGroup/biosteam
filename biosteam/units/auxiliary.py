# -*- coding: utf-8 -*-
"""
This module contains functions for adding auxliary unit operations.
"""
from . import design_tools as design
import biosteam as bst

__all__ = ('VacuumSystem',)

class VacuumSystem: # Can function as an auxiliary unit operation
    __slots__ = (
        'power_utility',
        'heat_utilities', 
        'baseline_purchase_costs',
        'purchase_costs',
        'installed_costs',
        'F_M', 'F_D', 'F_P', 'F_BM',
    )
    def __init__(self, F_mass, F_vol, P, volume, vacuum_system_preference):
        self.baseline_purchase_costs = self.purchase_costs = self.installed_costs = capex = {} # Assume all costs the same
        self.F_M = self.F_D = self.F_P = self.F_BM = {} # No factors
        vacuum_results = design.compute_vacuum_system_power_and_cost(
            F_mass, F_vol, P, volume, self.vacuum_system_preference
        )
        name = vacuum_results['Name']
        capex[name] = vacuum_results['Cost']
        heating_agent = vacuum_results['Heating agent']
        if heating_agent: # Steam turbine
            vacuum_steam = bst.HeatUtility()
            vacuum_steam.set_utility_by_flow_rate(heating_agent, vacuum_results['Steam flow rate'])
            if vacuum_results['Condenser']: 
                vacuum_cooling_water = bst.HeatUtility()
                vacuum_cooling_water(-vacuum_steam.unit_duty, 373.15)
        self.power_utility = bst.PowerUtility(vacuum_results['Work'])

def run_auxiliary_heat_exchanger_loop(unit):
    if hasattr(unit, 'heat_exchanger'):
        pass
    else: # Create heat exchanger and add them
        pass