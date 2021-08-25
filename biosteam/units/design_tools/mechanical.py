# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
General functional algorithms for the design of pumps and motors.

"""
import numpy as np
from numba import njit

__all__ = ('brake_efficiency', 'motor_efficiency', 'pump_efficiency',
           'nearest_NEMA_motor_size')

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
    """Return motor efficiency given brake power in hp."""
    if Pb < 1: Pb = 1
    elif Pb > 1500: Pb = 1500
    logPb = np.log(Pb)
    return 0.8 + 0.0319*logPb - 0.00182*logPb*logPb

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
    raise ValueError(f'no NEMA motor size bigger than {power} hp')

@njit(cache=True)    
def calculate_NPSH(P_suction, P_vapor, rho_liq):
    """Return NPSH in ft given suction and vapor pressure in Pa and density in kg/m^3."""
    # Note: NPSH = (P_suction - P_vapor)/(rho_liq*gravity)
    # Taking into account units, NPSH will be equal to return value
    return 0.334438*(P_suction - P_vapor)/rho_liq