# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
.. contents:: :local:

.. autoclass:: biosteam.units.vacuum_system.VacuumSystem

"""
from . import design_tools as design
from .decorators import auxiliary
from typing import Optional, Union
import biosteam as bst
from .._unit import Unit

__all__ = ('VacuumSystem',)
       
@auxiliary
class VacuumSystem: 
    """
    Create an auxiliary vacuum system for a unit operation. 
    
    Parameters
    ----------
    unit : 
        Main unit operation containing the auxiliary vacuum system. Missing parameters
        are estimated through this unit operation, if given.
    F_mass : 
        Vapor mass flow rate entering vacuum system from vessel [kg/hr] (not including inleakage).
        Defaults to the vapor volumetric flow rate exiting `unit`.
    F_vol : 
        Vapor volumetric flow rate entering vacuum system from vessel [m3/hr] (not including inleakage).
        Defaults to the vapor mass flow rate exiting `unit`.
    P_suction : 
        Suction pressure [Pa]. Defaults to `unit.P` or the minimum outlet pressure.
    vessel_volume : 
        Vacuum volume [m3]. Defaults to `unit.design_results['Total volume']`.
    vacuum_system_preference : 
        Name(s) of preferred vacuum systems. Valid names include 'Liquid-ring pump', 
        'Steam-jet ejector', and 'Dry-vacuum pump'.
    
    Notes
    -----
    The vacuum system is sized/costed based on the vapor flow rate through the 
    vacuum system, which includes the inleakage into the vessel through fixtures 
    and cracks. The inleakage is a function of the suction pressure and the 
    total vessel volume as in [1]_. 
    
    BioSTEAM's CSTR, Flash, and Distillation columns automatically include a 
    vacuum system when needed.
    
    Examples
    --------
    Create vessel with a vacuum system:
        
    >>> import biosteam as bst
    >>> class VacuumVessel(bst.Unit):
    ...     auxiliary_unit_names = ('vacuum_system',) # Mark attributes as auxiliary
    ...     _units = {'Total volume': 'm3'} # This is needed for the vacuum system
    ...     P = 1000 # Pa
    ...     tau = 4 # hr
    ...
    ...     def _run(self):
    ...         self.outs[0].P = 1000 # Pa
    ...     
    ...     def _design(self):
    ...         self.design_results['Total volume'] = self.feed.F_vol * self.tau
    ...         self.vacuum_system = bst.VacuumSystem(self)
    ...
    >>> bst.settings.set_thermo(['Water'])
    >>> feed = bst.Stream('feed', Water=1e6)
    >>> V1 = VacuumVessel('V1', ins=feed)
    >>> V1.simulate()
    >>> V1.results()
    Vacuum vessel                                                Units       V1
    Medium pressure steam Duty                                   kJ/hr 2.09e+07
                          Flow                                 kmol/hr      575
                          Cost                                  USD/hr      159
    Design                Total volume                              m3 7.23e+04
    Purchase cost         Vacuum system - Steam-jet ejecto...      USD 8.18e+04
    Total purchase cost                                            USD 8.18e+04
    Utility cost                                                USD/hr      159
    
    For simplicity, this example does not include the cost of the vessel, but 
    vessel costs should be included for techno-economic analysis.
    
    References
    ----------
    .. [1] Seider, W. D.; Lewin, D. R.; Seader, J. D.; Widagdo, S.; Gani, R.;
        Ng, M. K. Cost Accounting and Capital Cost Estimation.
        In Product and Process Design Principles; Wiley, 2017; pp 426–485.
    
    """
    __slots__ = (
        'power_utility',
        'heat_utilities', 
        'baseline_purchase_costs',
        'purchase_costs',
        'installed_costs',
        'F_M', 'F_D', 'F_P', 'F_BM',
        'owner',
    )
    
    def __init__(self,
            unit: Optional[Unit]=None, 
            vacuum_system_preference: Optional[str]=None, 
            F_mass: Optional[float]=None,
            F_vol: Optional[float]=None, 
            P_suction: Optional[float]=None, 
            vessel_volume: Optional[Union[float, list[float]]]=None
        ):
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
                                     "'vessel_volume' parameter could not be deduced")
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
                heat_utilities.append(vacuum_steam)
                vacuum_cooling_water(-vacuum_steam.unit_duty, 373.15)
        self.power_utility = bst.PowerUtility(vacuum_results['Work'])