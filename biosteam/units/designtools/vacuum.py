# -*- coding: utf-8 -*-
"""
Created on Fri Feb 22 17:31:50 2019

@author: yoelr
"""
import numpy as np
from biosteam.utils import checkbounds
from biosteam.exceptions import DesignError

__all__ = ('vacuum_system', '_calc_MotorEfficiency', '_calc_BreakEfficiency')

ln = np.log

# %% Data

# System types of vacuum systems
# Volumetric flowrate ranges, (cfm) and lower limit of suction (torr)
_steamjet_ejectors = {
    'One stage':               ((10, 1000000), 100),
    'Two stage':               ((10, 1000000),  15),
    'Three stage':             ((10, 1000000),   2)}
_liquid_ring = {
    'One stage water seal':    (( 3,   18000),  50),
    'Two stage water seal':    (( 3,   18000),  25),
    'Oil seal':                (( 3,   18000),  10)}
_dry_vacuum = {
    'Three-stage rotary lobe': ((60,     240), 1.5),
    'Three-stage claw':        ((60,     270), 0.3),
    'Screw compressor':        ((50,    1400), 0.1)}

_default_vacuum_systems = {'Liquid-ring pump': _liquid_ring,
                           'Steam-jet ejector': _steamjet_ejectors,
                           'Dry-vacuum pump': _dry_vacuum}
            
_air_density = 1.2041 # kg/m3 dry air

# %% Calculate vacuum system requirements

def vacuum_system(massflow:'kg/hr', volflow:'m3/hr',
                  P_suction:'Pa', vol:'m3',
                  vacuum_system_preference=None):
    """Return dictionary of results
    
    **Parameters**
    
        **massflow:** [float] Vapor mass flow rate entering vacuum system from vessel (not including inleakage) 
        
        **volflow:** [float] Vapor volumetric flow rate entering vacuum system from vessel (not including inleakage)
        
        **P_suction:** [float] Suction pressure
        
        **vol:** [float] Vacuum volume
        
        **vacuum_system_preference:** tuple[str] or [str] Name(s) of preferred vacuum systems
    
    **Return**
    
        **power:** [float] (kW)
        
        **cost:** [float] (USD)
    
    """
    vacuum_systems = _prefered_VacuumSystems(vacuum_system_preference)
    P_suction *= 7.5006e-3 # to torr
    if vol:
        vol *= 35.315 # to ft3
        air_massflow = _calc_AirInleakage(vol, P_suction) # lb/hr
        air_volflow = 0.26697*air_massflow/_air_density # cfm
    else:
        air_massflow = 0
        air_volflow = 0
    volflow_cfm = 0.5886*volflow + air_volflow
    if volflow_cfm < 3.01:
        factor = 3.01/volflow_cfm
        volflow_cfm = 3.01
    else:
        factor = 1
    massflow_kgph = (massflow + 0.4536*air_massflow)*factor # kg/hr
    massflow_lbph = 2.205*massflow_kgph
    power = _power(massflow_kgph, P_suction)
    vacuum_sys, grade = _select_VacuumSystem(vacuum_systems, volflow_cfm, P_suction)
    cost = _cost(vacuum_sys, grade, massflow_lbph, volflow_cfm, P_suction)
    return power, cost

# %% Supporting functions

def _calc_BreakEfficiency(q:'gpm'):
    if q < 50: q = 50
    elif q > 5000: q = 5000
    return -0.316 + 0.24015*ln(q) - 0.01199*ln(q)**2

def _calc_MotorEfficiency(Pb):
    if Pb < 1: Pb = 1
    elif Pb > 1500: Pb = 1500
    return 0.8 + 0.0319*ln(Pb) - 0.00182*ln(Pb)**2

def _prefered_VacuumSystems(preference):
    if preference is None:
        vacuum_systems = _default_vacuum_systems.values()
    else:
        defaults = _default_vacuum_systems.keys()
        if isinstance(preference, str):
            if preference not in defaults:
                raise ValueError(f"Preference have at least one of the following: {defaults}.")
            preference = (preference,)
        else:
            for name in preference:
                if name not in defaults:
                    raise ValueError(f"Preference have at least one of the following: {defaults}.")
        
        vacuum_systems = [_default_vacuum_systems[name] for name in preference]
    return vacuum_systems


def _available_VacuumSystems(flow:'cfm', P_suction:'torr'):
    types = []
    for vacuumtype, vacuum_sys in _default_vacuum_systems.items():
        for grade, flowrange_minsuction in vacuum_sys.items():
            flowrange, minsuction = flowrange_minsuction
            if checkbounds(flow, flowrange) and P_suction > minsuction:
                types.append((vacuumtype, grade))
    return types

def _select_VacuumSystem(vacuum_systems, volflow, P_suction):
    for vacuum_sys in vacuum_systems:
        for grade, flowrange_minsuction in vacuum_sys.items():
            flowrange, minsuction = flowrange_minsuction
            if checkbounds(volflow, flowrange) and P_suction > minsuction:
                return (vacuum_sys, grade)
    raise DesignError('No vacuum system available at current flow  and suction pressure.')

def _lazycalc_AirInleakage(V:'ft^3', P:'Pa') -> 'lb/hr':
    if P > 11999.013: k = 0.2
    elif P > 4132.993: k = 0.15
    elif P > 2799.77: k = 0.10
    elif P > 133.322: k = 0.051
    else: raise ValueError('Cannot calculate air inleakage at pressures lower than 133.32 Pascal.')
    return k*V**0.667

def _calc_AirInleakage(V:'ft^3', P:'torr') -> 'lb/hr':
    return 5 + (0.0298 + 0.03088*ln(P) - 5.733e-4*ln(P)**2)*V**0.66
    
def _power(massflow:'kg/hr',  P_suction:'torr') -> 'kW':
    SF = massflow/P_suction # Size factor
    if SF < 0.2: SF = 0.2
    elif SF > 16: SF = 16
    Pb = 21.4*SF**0.924 # Break power assuming Liquid-ring (NASH)
    mu_m = _calc_MotorEfficiency(Pb)
    return Pb/mu_m

def _cost(vacuum_sys, grade, massflow:'lb/hr', volflow:'cfm', P_suction:'torr') -> 'USD':
    # Costing for different vacuums ASSUMING THAT ALL ARE CARBON STEEL!!!!
    if vacuum_sys is _steamjet_ejectors:
        S = massflow/P_suction
        Cp = 1915*(S**(0.41))
        if grade == 'One stage':
            Cs = 1
        elif grade == 'Two stage':
            Cs = 1.8
        elif grade == 'Three stage':
            Cs = 2.1
        Cost = Cp*Cs
    elif vacuum_sys is _liquid_ring:
        S = volflow
        Cp = 8250*(S**(0.37))
        if grade == 'One stage water seal':
            Cs = 1
        elif grade == 'Two stage water seal':
            Cs = 1.8
        elif grade == 'Oil seal':
            Cs = 2.1
        Cost = Cp*Cs
    elif vacuum_sys is _dry_vacuum:
        S = volflow
        if grade == 'Three-stage rotary lobe':
            Cp = 8075*(S**(0.41))
        elif grade == 'Three-stage claw':
            Cp = 9785*(S**(0.36))
        elif grade == 'Screw compressor':
            Cp = 10875*(S**(0.38))
        Cost = Cp
    return Cost
    
