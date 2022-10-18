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
        'owner',
    )
    def __init__(self, unit=None, vacuum_system_preference=None, 
                 F_mass=None, F_vol=None, P_suction=None, vessel_volume=None):
        if unit:
            # Deduce arguments if unit given
            if F_mass is None or F_vol is None:
                # If vapor is condensed, assume vacuum system is after condenser
                vents = [i for i in unit.outs if 'g' in i.phase]
                F_mass = 0.
                F_vol = 0.
                for vapor in vents:
                    hx = vapor.sink
                    if isinstance(hx, bst.HX):
                        index = hx.ins.index(vapor)
                        vapor = hx.outs[index]
                        if 'g' not in vapor.phase: continue
                    if isinstance(vapor, bst.MultiStream):
                        F_mass += vapor['g'].F_mass
                        F_vol += vapor['g'].F_vol
                    else:
                        F_mass += vapor.F_mass
                        F_vol += vapor.F_vol
            if P_suction is None:
                P_suction = getattr(unit, 'P', None) or getattr(unit, 'P_suction', None)
                if P_suction is None: P_suction = min([i.P for i in unit.outs])
            if vessel_volume is None:
                if 'Total volume' in unit.design_results:
                    vessel_volume = unit.get_design_result('Total volume', 'm3')
                else:
                    raise ValueError("'Total volume' was not found in design results; "
                                     "'vesse_volume' parameter could not be deduced")
        else:
            # In case user does not supply all arguments
            params = [('F_mass', F_mass), 
                      ('F_vol', F_vol),
                      ('P_suction', P_suction), 
                      ('vessel_volume', vessel_volume)]
            for name, value in params:
                if value is None: 
                    raise ValueError(
                        f"missing argument '{name}'; cannot deduce '{name}'"
                         "when no unit is specified"
                    )
        self.baseline_purchase_costs = self.purchase_costs = self.installed_costs = capex = {} # Assume all costs the same
        self.F_M = self.F_D = self.F_P = self.F_BM = {} # No factors
        vacuum_results = design.compute_vacuum_system_power_and_cost(
            F_mass, F_vol, P_suction, vessel_volume, vacuum_system_preference
        )
        name = vacuum_results['Name']
        capex[name] = vacuum_results['Cost']
        heating_agent = vacuum_results['Heating agent']
        self.heat_utilities = heat_utilities = []
        if heating_agent: # Steam turbine
            vacuum_steam = bst.HeatUtility()
            heat_utilities.append(vacuum_steam)
            vacuum_steam.set_utility_by_flow_rate(heating_agent, vacuum_results['Steam flow rate'])
            if vacuum_results['Condenser']: 
                vacuum_cooling_water = bst.HeatUtility()
                vacuum_cooling_water.append(vacuum_steam)
                vacuum_cooling_water(-vacuum_steam.unit_duty, 373.15)
        self.power_utility = bst.PowerUtility(vacuum_results['Work'])

    def _load_costs(self): pass # Necessary as an auxiliary unit operation