# -*- coding: utf-8 -*-
"""
Functional algorithms for the design and purchase cost estimation of
vacuum systems.

References
----------
.. [1] Seider, W. D.; Lewin, D. R.; Seader, J. D.; Widagdo, S.; Gani, R.;
    Ng, M. K. Cost Accounting and Capital Cost Estimation.
    In Product and Process Design Principles; Wiley, 2017; pp 426–485.


"""
from numpy import log as ln
from . import mechanical as mch
from biosteam.utils import checkbounds
from biosteam.exceptions import DesignError
from flexsolve import njitable
import biosteam as bst

__all__ = ('compute_vacuum_system_power_and_cost',)


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

def compute_vacuum_system_power_and_cost(
        F_mass, F_vol, P_suction, vessel_volume,
        vacuum_system_preference=None):
    """
    Return power (kW) and cost (USD) of vacuum system.
    
    Parameters
    ----------
    F_mass : float
        Vapor mass flow rate entering vacuum system from vessel in kg/hr (not including inleakage).
    F_vol : float
        Vapor volumetric flow rate entering vacuum system from vessel in m3/hr (not including inleakage).
    P_suction : float
        Suction pressure in Pa
    vessel_volume : float
        Vacuum volume in m3
    vacuum_system_preference : 'Liquid-ring pump', 'Steam-jet ejector', or 'Dry-vacuum pump'
        Name(s) of preferred vacuum systems
    
    """
    P_suction *= 7.5006e-3 # to torr
    if vessel_volume:
        vessel_volume *= 35.315 # to ft3
        F_mass_air = calculate_air_inleakage(vessel_volume, P_suction) # lb/hr
        F_vol_air = 0.26697*F_mass_air/_air_density # cfm
    else:
        F_vol_air = F_mass_air = 0
    F_vol_cfm = 0.5886*F_vol + F_vol_air
    if F_vol_cfm < 3.01:
        factor = 3.01/F_vol_cfm
        F_vol_cfm = 3.01
    else:
        factor = 1
    F_mass_kgph = (F_mass + 0.4536*F_mass_air)*factor # kg/hr
    F_mass_lbph = 2.205 * F_mass_kgph
    power = calculate_vacuum_power(F_mass_kgph, P_suction)
    vacuum_systems = get_prefered_vacuum_systems(vacuum_system_preference)
    vacuum_sys, grade = select_vacuum_system(vacuum_systems, F_vol_cfm, P_suction)
    base_cost = calculate_vacuum_cost(vacuum_sys, grade, F_mass_lbph, F_vol_cfm, P_suction)
    cost = bst.CE / 567  * base_cost
    return power, cost


# %% Supporting functions

def get_prefered_vacuum_systems(preference):
    if preference is None:
        vacuum_systems = _default_vacuum_systems.values()
    else:
        defaults = _default_vacuum_systems.keys()
        if isinstance(preference, str):
            if preference not in defaults:
                raise ValueError(f"preference have at least one of the following: {defaults}")
            preference = (preference,)
        else:
            for name in preference:
                if name not in defaults:
                    raise ValueError(f"preference have at least one of the following: {defaults}")
        
        vacuum_systems = [_default_vacuum_systems[name] for name in preference]
    return vacuum_systems


def get_available_vacuum_systems(F_vol_cfm, P_suction):
    """
    Return available vacuum type and grade
    
    Parameters
    ----------
    F_vol_cfm : float
        Vapor volumetric flow rate entering vacuum system from vessel in cfm (including inleakage).
    P_suction : float
        Suction pressure in Torr
    """
    types = []
    for vacuumtype, vacuum_sys in _default_vacuum_systems.items():
        for grade, flowrange_minsuction in vacuum_sys.items():
            flowrange, minsuction = flowrange_minsuction
            if checkbounds(F_vol_cfm, flowrange) and P_suction > minsuction:
                types.append((vacuumtype, grade))
    return types

def select_vacuum_system(vacuum_systems, F_vol_cfm, P_suction):
    """
    Return a heuristic vacuum type and grade
    
    Parameters
    ----------
    F_vol_cfm : float
        Vapor volumetric flow rate entering vacuum system from vessel in cfm (including inleakage).
    P_suction : float
        Suction pressure in Torr
    """
    for vacuum_sys in vacuum_systems:
        for grade, flowrange_minsuction in vacuum_sys.items():
            flowrange, minsuction = flowrange_minsuction
            if checkbounds(F_vol_cfm, flowrange) and P_suction > minsuction:
                return (vacuum_sys, grade)
<<<<<<< HEAD
    raise DesignError('no vacuum system available at current flow and suction pressure')

=======
    for vacuum_sys in vacuum_systems:
        for grade, flowrange_minsuction in vacuum_sys.items():
            flowrange, minsuction = flowrange_minsuction
            if F_vol_cfm < flowrange[-1] and P_suction > minsuction:
                return (vacuum_sys, grade)
    raise DesignError('no vacuum system available at current flow and suction pressure')

@njitable(cache=True)
>>>>>>> cd2c5013aaf9b5bc94bb764b52fd37db183472f1
def calculate_heuristic_air_inleakage(V, P):
    """
    Return air in-leakage in kg/hr through a heuristic calculation.
    
    Parameters
    ----------
    V : float
        Vacuum volume in m3
    P : float
        Suction pressure in Pa
    
    """
    if P > 11999.013: k = 0.2
    elif P > 4132.993: k = 0.15
    elif P > 2799.77: k = 0.10
    elif P > 133.322: k = 0.051
    else: raise ValueError('cannot calculate air inleakage at pressures lower than 133.32 Pascal')
    return k*V**0.667

<<<<<<< HEAD
@njitable
=======
@njitable(cache=True)
>>>>>>> cd2c5013aaf9b5bc94bb764b52fd37db183472f1
def calculate_air_inleakage(V, P):
    """
    Return air in-leakage in kg/hr.
    
    Parameters
    ----------
    V : float
        Vacuum volume in m3
    P : float
        Suction pressure in Torr
    """
    return 5 + (0.0298 + 0.03088*ln(P) - 5.733e-4*ln(P)**2)*V**0.66

def calculate_vacuum_power(F_mass,  P_suction):
    """
    Return vacuum power (after accounting for motor efficiency) in kW.
    
    Parameters
    ----------
    F_mass : float
        Total mass flow rate entering vacuum system in kg/hr (including inleakage).
    P_suction : float
        Suction pressure in Torr
    
    """
    SF = F_mass/P_suction # Size factor
    if SF < 0.2: SF = 0.2
    elif SF > 16: SF = 16
    Pb = 21.4*SF**0.924 # Break power assuming Liquid-ring (NASH)
    mu_m = mch.motor_efficiency(Pb)
    return Pb/mu_m

def calculate_vacuum_cost(vacuum_sys, grade, F_mass_lbph, F_vol_cfm, P_suction):
    # Costing for different vacuums ASSUMING THAT ALL ARE CARBON STEEL!!!!
    if vacuum_sys is _steamjet_ejectors:
        S = F_mass_lbph/P_suction
        Cp = 1915*(S**(0.41))
        if grade == 'One stage':
            Cs = 1
        elif grade == 'Two stage':
            Cs = 1.8
        elif grade == 'Three stage':
            Cs = 2.1
        Cost = Cp*Cs
    elif vacuum_sys is _liquid_ring:
        S = F_vol_cfm
        Cp = 8250*(S**(0.37))
        if grade == 'One stage water seal':
            Cs = 1
        elif grade == 'Two stage water seal':
            Cs = 1.8
        elif grade == 'Oil seal':
            Cs = 2.1
        Cost = Cp*Cs
    elif vacuum_sys is _dry_vacuum:
        S = F_vol_cfm
        if grade == 'Three-stage rotary lobe':
            Cp = 8075*(S**(0.41))
        elif grade == 'Three-stage claw':
            Cp = 9785*(S**(0.36))
        elif grade == 'Screw compressor':
            Cp = 10875*(S**(0.38))
        Cost = Cp
    return Cost
    
