# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
General functional algorithms for the design of pumps and motors.

"""
import biosteam as bst
import numpy as np
from math import log, exp
from numba import njit, objmode

__all__ = ('brake_efficiency', 'motor_efficiency', 'pump_efficiency',
           'nearest_NEMA_motor_size', 'electric_motor_cost',
           )

#: [tuple] All NEMA motor sizes in increasing order, in horsepower.
nema_sizes_hp = (0.25, 0.3333333333333333, 0.5, 0.75, 1.0, 1.5, 2.0,
                 3.0, 4.0, 5.0, 5.5, 7.5, 10.0, 15.0, 20.0, 25.0, 30.0,
                 40.0, 50.0, 60.0, 75.0, 100.0, 125.0, 150.0, 175.0, 200.0, 
                 250.0, 300.0, 350.0, 400.0, 450.0, 500.0)


@njit(cache=True)
def brake_efficiency(q):
    """Return brake efficiency given flow rate in gpm."""
    if q < 50: q = 50
    elif q > 5000: q = 5000
    logq = np.log(q)
    return -0.316 + 0.24015*logq - 0.01199*logq*logq

@njit(cache=True)
def motor_efficiency(Pb):
    """Return motor efficiency given brake power in hp (used for estimating pump efficiency)."""
    if Pb < 1: Pb = 1
    elif Pb > 1500: Pb = 1500
    logPb = np.log(Pb)
    return 0.8 + 0.0319*logPb - 0.00182*logPb*logPb

def electric_motor_cost(Pc):
    """
    Return the baseline purchase cost of an electric motor given
    the shaft power in hp.
    """
    return _electric_motor_cost(Pc, bst.CE) 

@njit(cache=True)
def _electric_motor_cost(Pc, CE):
    lnp = log(Pc)
    lnp2 = lnp*lnp
    lnp3 = lnp2*lnp
    lnp4 = lnp3*lnp
    return exp(5.9332 + 0.16829*lnp
               - 0.110056*lnp2 + 0.071413*lnp3
               - 0.0063788*lnp4) * CE / 567

@njit(cache=True)
def noncondensing_steam_turbine_cost(Pc):
    """
    Return the baseline purchase cost of a noncondensing steam turbine given
    the shaft power in hp.
    """
    return 10660. * Pc ** 0.41

@njit(cache=True)
def condensing_steam_turbine_cost(Pc):
    """
    Return the baseline purchase cost of a condensing steam turbine given
    the shaft power in hp.
    """
    return 28350. * Pc ** 0.405

@njit(cache=True)
def gas_turbine_cost(Pc):
    """
    Return the baseline purchase cost of a gas turbine given
    the shaft power in hp.
    """
    return 2835. * Pc ** 0.76

@njit(cache=True)
def internal_combustion_engine_cost(Pc):
    """
    Return the baseline purchase cost of an internal combustion engine given
    the shaft power in hp.
    """
    return 1588. * Pc ** 0.75

@njit(cache=True)
def pump_efficiency(q, p):
    """Return pump efficiency.
    
    Parameters
    ----------
    q : float
        Volumetric flow rate in gpm.
    p : float
        Power in hp.
    """
    mup = brake_efficiency(q)
    mum = motor_efficiency(p/mup)
    return mup*mum

def nearest_NEMA_motor_size(power):
    for nearest_power in nema_sizes_hp:
        if nearest_power >= power: return nearest_power
    return nearest_power

@njit(cache=True)    
def calculate_NPSH(P_suction, P_vapor, rho_liq):
    """Return NPSH in ft given suction and vapor pressure in Pa and density in kg/m^3."""
    # Note: NPSH = (P_suction - P_vapor)/(rho_liq*gravity)
    # Taking into account units, NPSH will be equal to return value
    return 0.334438*(P_suction - P_vapor)/rho_liq