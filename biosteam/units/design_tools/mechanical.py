# -*- coding: utf-8 -*-
"""
General functional algorithms for the design of pumps and motors.

"""
import numpy as np
from flexsolve import njitable
from fluids.pump import nema_sizes_hp

__all__ = ('brake_efficiency', 'motor_efficiency', 'pump_efficiency',
           'nearest_NEMA_motor_size')

@njitable
def brake_efficiency(q):
    """Return brake efficiency given flow rate in gpm."""
    if q < 50: q = 50
    elif q > 5000: q = 5000
    return -0.316 + 0.24015*np.log(q) - 0.01199*np.log(q)**2

@njitable
def motor_efficiency(Pb):
    """Return motor efficiency given brake power in hp."""
    if Pb < 1: Pb = 1
    elif Pb > 1500: Pb = 1500
    return 0.8 + 0.0319*np.log(Pb) - 0.00182*np.log(Pb)**2

@njitable
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

@njitable    
def calculate_NPSH(P_suction, P_vapor, rho_liq):
    """Return NPSH in ft given suction and vapor pressure in Pa and density in kg/m^3."""
    # Note: NPSH = (P_suction - P_vapor)/(rho_liq*gravity)
    # Taking into account units, NPSH will be equal to return value
    return 0.334438*(P_suction - P_vapor)/rho_liq