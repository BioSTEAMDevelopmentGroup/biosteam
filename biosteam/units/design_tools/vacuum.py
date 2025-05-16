# -*- coding: utf-8 -*-
"""
Functional algorithms for the design and purchase cost estimation of
vacuum systems.

References
----------
.. [1] Seider, W. D.; Lewin, D. R.; Seader, J. D.; Widagdo, S.; Gani, R.;
    Ng, M. K. Cost Accounting and Capital Cost Estimation.
    In Product and Process Design Principles; Wiley, 2017; pp 426–485.
.. [2] Amazon. Robinair (15115) VacuMaster Single Stage Vacuum Pump - Single-Stage, 1.5 CFM. 
    https://www.amazon.com/Robinair-15115-VacuMaster-Single-Vacuum/dp/B005CO9FDW?ref_=ast_sto_dp.
    Accessed on 09/28/2023.
.. [3] Amazon. Robinair (15300) VacuMaster Economy Vacuum Pump - 2-Stage, 3 CFM.
    https://www.amazon.com/Robinair-15300-VacuMaster-Economy-Vacuum/dp/B000O1E5UQ?ref_=ast_sto_dp
    Accessed on 09/28/2023.

"""
from numpy import log as ln
from . import mechanical as mch
from math import ceil
from biosteam.utils import checkbounds
from biosteam.exceptions import DesignError
from numba import njit
import biosteam as bst

__all__ = ('compute_vacuum_system_power_and_cost',)


# %% Data

# System types of vacuum systems
# Volumetric flowrate ranges, (cfm) and lower limit of suction (torr)
# Rotary vane pumps based on ref. [2], [3]
_steamjet_ejectors = {
    'One stage':               ((10, 1000000), 100),
    'Two stage':               ((10, 1000000),  15),
    'Three stage':             ((10, 1000000),   2)}
_liquid_ring = {
    'One stage water seal':    (( 3,   18000),  50),
    'Two stage water seal':    (( 3,   18000),  25),
    'Oil seal':                (( 3,   18000),  10)}
_dry_vacuum = {
    'Three stage rotary lobe': ((60,     240), 1.5),
    'Three stage claw':        ((60,     270), 0.3),
    'Screw compressor':        ((50,    1400), 0.1)}
_rotary_vane = {
    'One stage':               ((0,     1.51), 0.115),
    'Two stage':               ((1.5,   3.01), 0.035)}

_default_vacuum_systems = {'Liquid-ring pump': _liquid_ring,
                           'Steam-jet ejector': _steamjet_ejectors,
                           'Dry-vacuum pump': _dry_vacuum,
                           'Rotary-vane pump': _rotary_vane}

_default_rotary_vane_work_cost = {
    'One stage': (1/5 * 0.7457, 127*1.08), # hp to kW; 2023 USD (including tax & shipping)
    'Two stage': (1/3 * 0.7457, 248*1.08)
}

_air_density = 1.2041 # kg/m3 dry air

# %% Calculate vacuum system requirements

def compute_vacuum_system_power_and_cost(
        F_mass, F_vol, P_suction, vessel_volume,
        vacuum_system_preference=None):
    """
    Return a dictionary of vacuum system requirements.
    
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
    # if F_vol_cfm < 3.01:
        # factor = 3.01/F_vol_cfm
        # F_vol_cfm = 3.01
    # else:
    #     factor = 1
    factor = 1
    F_mass_kgph = (F_mass + 0.4536*F_mass_air)*factor # kg/hr
    F_mass_lbph = 2.205 * F_mass_kgph
    vacuum_systems = get_preferred_vacuum_systems(vacuum_system_preference)
    name, grade, N = select_vacuum_system(vacuum_systems, F_vol_cfm, P_suction, bool(vacuum_system_preference))
    base_cost = calculate_vacuum_cost(name, grade, F_mass_lbph, F_vol_cfm, P_suction)
    cost = bst.CE / 567.  * base_cost
    if name == 'Steam-jet ejector':
        for agent in bst.HeatUtility.heating_agents:
            if agent.P > 689475.: break # 100 psig
        steam = 0.41631 * F_mass_kgph # [kmol/hr] 7.5 weight steam/ weight gas
        work = 0.
        has_condenser = grade != 'One stage'
    elif name == 'Rotary-vane pump':
        has_condenser = False
        agent = None
        steam = 0.
        N = 1
        work = _default_rotary_vane_work_cost[grade][0]
    else:
        has_condenser = False
        agent = None
        steam = 0.
        work = calculate_vacuum_power(F_mass_kgph, P_suction)
    return {'Work': work, 
            'Cost': N * cost, 
            'Name': f"{name}, {grade.lower()}",
            'In parallel': N,
            'Condenser': has_condenser,
            'Steam flow rate': steam, 
            'Heating agent': agent}


# %% Supporting functions

def get_preferred_vacuum_systems(preference):
    defaults = _default_vacuum_systems
    if preference is None:
        return defaults
    else:
        if isinstance(preference, str):
            if preference not in defaults:
                raise ValueError(f"preference must be in the following list: {list(defaults)}")
            preference = (preference, *defaults)
        else:
            for name in preference:
                if name not in defaults:
                    raise ValueError(f"preference must be in the following list: {list(defaults)}")
        return {name: defaults[name] for name in preference}

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

def select_vacuum_system(vacuum_systems, F_vol_cfm, P_suction, ignore_F_lb=False):
    """
    Return a heuristic vacuum type and grade
    
    Parameters
    ----------
    F_vol_cfm : float
        Vapor volumetric flow rate entering vacuum system from vessel in cfm (including inleakage).
    P_suction : float
        Suction pressure in Torr
    """
    for name, vacuum_sys in vacuum_systems.items():
        for grade, flowrange_minsuction in vacuum_sys.items():
            flowrange, minsuction = flowrange_minsuction
            if ignore_F_lb:
                if F_vol_cfm < flowrange[-1] and P_suction > minsuction:
                    return (name, grade, 1)
            elif checkbounds(F_vol_cfm, flowrange) and P_suction > minsuction:
                    return (name, grade, 1)
    for name, vacuum_sys in vacuum_systems.items(): # Flow rate too large
        for grade, flowrange_minsuction in vacuum_sys.items():
            flowrange, minsuction = flowrange_minsuction
            if P_suction > minsuction:
                N = ceil(F_vol_cfm / flowrange[-1])
                return (name, grade, N)
    raise DesignError(f'no vacuum system available at current flow ({F_vol_cfm:.2f} cfm) and suction pressure ({P_suction:.2f} Torr)')

@njit(cache=True)
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

@njit(cache=True)
def calculate_air_inleakage(V, P):
    """
    Return air in-leakage in kg/hr.
    
    Parameters
    ----------
    V : float
        Vacuum volume in ft3
    P : float
        Suction pressure in Torr
    """
    lnP = ln(P) 
    return 5. + (0.0298 + 0.03088*lnP - 5.733e-4*lnP*lnP)*V**0.66

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
    elif SF > 16.: SF = 16.
    Pb = 21.4*SF**0.924 # Break power assuming Liquid-ring (NASH)
    mu_m = mch.motor_efficiency(Pb)
    return Pb/mu_m

def calculate_vacuum_cost(vacuum_sys, grade, F_mass_lbph, F_vol_cfm, P_suction):
    # Costing for different vacuums ASSUMING THAT ALL ARE CARBON STEEL!!!!
    if vacuum_sys ==  'Steam-jet ejector':
        S = F_mass_lbph/P_suction
        Cp = 1915*S**0.41
        if grade == 'One stage':
            Cs = 1
        elif grade == 'Two stage':
            Cs = 1.8 * 1.6 # 2 stage + 1 condenser
        elif grade == 'Three stage':
            Cs = 2.1 * 2.3 # 3 stage + 2 condenser
        Cost = Cp*Cs
    elif vacuum_sys == 'Liquid-ring pump':
        S = F_vol_cfm
        Cp = 8250*S**0.37
        if grade == 'One stage water seal':
            Cs = 1
        elif grade == 'Two stage water seal':
            Cs = 1.8
        elif grade == 'Oil seal':
            Cs = 2.1
        Cost = Cp*Cs
    elif vacuum_sys == 'Dry-vacuum pump':
        S = F_vol_cfm
        if grade == 'Three stage rotary lobe':
            Cp = 8075*S**0.41
        elif grade == 'Three stage claw':
            Cp = 9785*S**0.36
        elif grade == 'Screw compressor':
            Cp = 10875*S**0.38
        Cost = Cp
    elif vacuum_sys == 'Rotary-vane pump':
        Cost = _default_rotary_vane_work_cost[grade][1] / 708. * 567.   # !!! 708 is the 2021 CEPCI, need to be updated to 2023 CEPCI
    return Cost
    

    
    