# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
General functional algorithms for the design and purchase cost estimation
of distillation columns.

References
----------
.. [1] Seider, W. D., Lewin,  D. R., Seader, J. D., Widagdo, S., Gani, R.,
    & Ng, M. K. (2017). Product and Process Design Principles. Wiley.
    Cost Accounting and Capital Cost Estimation (Chapter 16)
.. [2] M. Duss, R. Taylor. (2018)
    Predict Distillation Tray Efficiency. AICHE 
.. [3] Green, D. W. Distillation. In Perry’s Chemical Engineers’
    674 Handbook, 9 ed.; McGraw-Hill Education, 2018.

"""
import numpy as np
from . import utils
from numba import njit
import biosteam as bst
from warnings import warn

__all__ = ('compute_purchase_cost_of_trays',
           'compute_empty_tower_cost',
           'compute_plaform_ladder_cost',
           'compute_tower_weight',
           'compute_tower_wall_thickness',
           'compute_tray_base_purchase_cost',
           'compute_n_trays_factor',
           'compute_murphree_stage_efficiency',
           'compute_flow_parameter',
           'compute_max_capacity_parameter',
           'compute_max_vapor_velocity',
           'compute_downcomer_area_fraction',
           'compute_tower_diameter',
           'compute_tower_height')

@njit(cache=True)
def minimum_thickness_from_diameter(D):
    return 0.03125 * D + 0.125

@njit(cache=True)
def compute_purchase_cost_of_trays(N_T, Di):
    """
    Return total cost of all trays at BioSTEAM's CEPCI.
    
    Parameters
    ----------
    N_T : int
        Number of trays.
    Di : float
        Inner diameter [ft].
    
    Notes
    -----
    The purchase cost is given by [1]_. See source code for details.
    The purchase cost is scaled according to BioSTEAM's Chemical
    Plant Cost Index, `biosteam.CE`.
    
    """
    F_CE = bst.CE/500.
    C_BT = compute_tray_base_purchase_cost(Di)
    F_NT = compute_n_trays_factor(N_T)
    return N_T * F_CE * F_NT * C_BT

@njit(cache=True)
def compute_empty_tower_cost(W):
    """
    Return the cost [C_V; in USD] of an empty tower vessel at BioSTEAM's CEPCI.
    
    Parameters
    ----------
    W : float
        Weight [lb].
    
    
    Notes
    -----
    The purchase cost is given by [1]_. See source code for details.
    
    """
    logW = np.log(W)
    return bst.CE/500. * np.exp(7.2756 + 0.18255*logW + 0.02297*logW*logW)

@njit(cache=True)
def compute_plaform_ladder_cost(Di, L):
    """
    Return the cost [C_PL; in USD] of platforms and ladders at BioSTEAM's CEPCI.
    
    Parameters
    ----------
    Di: float
        Inner diameter [ft].
    L: float
        Legnth [ft].
    
    Notes
    -----
    The purchase cost is given by [1]_. See source code for details.
    
    """
    return bst.CE/500. * 300.9*Di**0.63316*L**0.80161

@njit(cache=True)
def compute_tower_weight(Di, L, tv, rho_M):
    """
    Return the weight [W; in lb] of the tower assuming 2:1 elliptical head.
    
    Parameters
    ----------
    Di : float
        Inner diameter [ft].
    L :  float
        Legnth [ft].
    tv : float
        Shell thickness [in].
    rho_M: floa
        Density of material [lb/in^3].
    
    Notes
    -----
    The tower weight is given by [1]_. See source code for details.
    
    """
    Di = Di*12.
    L = L*12.
    return np.pi*(Di+tv)*(L+0.8*Di)*tv*rho_M

@njit(cache=True)
def compute_tower_wall_thickness(Po, Di, L, S=15000., E=None, M=29.5):
    """
    Return the wall thinkness [tv; in inches] designed to withstand the
    internal pressure and the wind/earthquake load at the bottom.
    
    Parameters
    ----------
    Po : float
        Operating internal pressure [psi].
    Di : float
        Internal diameter [ft].
    L : float
        Height [ft].
    S : float
        Maximum stress [psi].
    E : float
        Fractional weld efficiency
    M : float
        Elasticity [psi].
    
    Notes
    -----
    The wall thickness is given by [1]_. See source code for details.
    
    Warning
    -------
    This function is only applicable to positive internal pressures (no vacuums).
    Vacuum pressure vessels may require stiffening rings and higher vessel thickness.
    
    """
    # TODO: Incorporate temperature for choosing S and M
    Di = Di*12. # ft to in
    L = L*12.
    
    E_check = E is None
    if E_check:
        # Assume carbon steel with thickness more than 1.25 in
        E = 1.0 
    
    # Get design pressure, which should be higher than operating pressure.
    Po_gauge = Po - 14.69
    if Po_gauge < 5.:
        Pd = 10.
    elif Po_gauge < 1000.:
        logPo = np.log(Po)
        Pd = np.exp(0.60608 + 0.91615*logPo) + 0.0015655*logPo*logPo
    else:
        Pd = 1.1*Po_gauge
    
    # Calculate thinkess according to ASME pressure-vessel code.
    ts = Pd*Di/(2.*S*E-1.2*Pd)
    
    if E_check:
        # Weld efficiency of 0.85 for low thickness carbon steel
        if ts < 1.25:
            E = 0.85
            ts = Pd*Di/(2.*S*E-1.2*Pd)
    
    # Add corrosion allowence
    ts += 1/8.
    
    # Minimum thickness for vessel rigidity may be larger
    Di_ft = Di/12.
    ts_min = minimum_thickness_from_diameter(Di_ft) if Di_ft > 4. else 0.25
    if ts < ts_min:
        ts = ts_min
    
    # Calculate thickness to withstand wind/earthquake load
    Do = Di + ts
    tw = 0.22*(Do + 18.)*L*L/(S*Do*Do)
    tv = tw if tw > ts else ts
    
    # Vessels are fabricated from metal plates with small increments
    if tv < 0.5:
        tv = utils.approx2step(tv, 3/16, 1/16)
    elif tv < 2.:
        tv = utils.approx2step(tv, 0.5, 1/8)
    elif tv < 3.:
        tv = utils.approx2step(tv, 2., 1/4)
    return tv

@njit(cache=True)
def compute_tray_base_purchase_cost(Di):
    """Return the base cost of a tray [C_BT; USD] at a CE of 500.
    
    Parameters
    ----------
    Di : float
        Inner diameter [ft].
    
    Notes
    -----
    The purchase cost is given by [1]_. See source code for details.
    
    """
    return 412.6985 * np.exp(0.1482*Di)

@njit(cache=True)
def compute_n_trays_factor(N_T):
    """
    Return the cost factor for number of trays, F_NT.
    
    Parameters
    ----------
    N_T: Number of trays
    
    Notes
    -----
    The cost factor is given by [1]_. See source code for details.
    
    """
    if N_T < 20.:
        F_NT = 2.25/1.0414**N_T
    else:
        F_NT = 1.
    return F_NT

@njit(cache=True)
def compute_murphree_stage_efficiency(mu, alpha, L, V):
    """
    Return the sectional murphree efficiency, E_mv.
    
    Parameters
    ----------
    mu: float
        Viscosity [mPa*s]
    alpha: float
        Relative volatility.    
    L: float
        Liquid flow rate by mol.
    V: float
        Vapor flow rate by mol.
    
    Notes
    -----
    The efficiency is given by [2]_. See source code for details.
    
    """
    S = alpha*V/L # Stripping factor
    e = 0.503*mu**(-0.226)*(S if S > 1. else 1./S)**(-0.08 )
    if e < 1.: return e
    else: return 1.

@njit(cache=True)
def compute_flow_parameter(L, V, rho_V, rho_L):
    """
    Return the flow parameter, F_LV.
    
    Parameters
    ----------
    L : float
        Liquid flow rate by mass.
    V : float
        Vapor flow rate by mass.
    rho_V : float
        Vapor density.
    rho_L : float
        Liquid density.
    
    Notes
    -----
    The flow parameter is given by [3]_. See source code for details.
    
    """
    return L/V*(rho_V/rho_L)**0.5

@njit(cache=True)
def compute_max_capacity_parameter(TS, F_LV):
    """Return the maximum capacity parameter before flooding [C_sbf; in m/s].
    
    Parameters
    ----------
    TS : float
        Tray spacing [mm].
    F_LV : float
        Flow parameter.
    
    Notes
    -----
    The max capacity parameter is given by [3]_. See source code for details.
    
    """
    return 0.0105 + 8.127e-4*TS**0.755*np.exp(-1.463*F_LV**0.842)

@njit(cache=True)
def compute_max_vapor_velocity(C_sbf, sigma, rho_L, rho_V, F_F, A_ha):
    """
    Return the maximum allowable vapor velocity
    through the net area of flow before flooding [U_f; in m/s].
    
    Parameters
    ----------
    C_sbf : 
        Maximum Capacity Parameter (m/s)
    sigma : 
        Liquid surface tension (dyn/cm)
    rho_L : 
        Liquid density
    rho_V : 
        Vapor density
    F_F : 
        Foaming factor
    A_ha : 
        Ratio of open area, A_h, to active area, A_a
    
    Notes
    -----
    The max vapor velocity is given by [3]_. See source code for details.
    
    """
    F_ST = (sigma/20.)**0.2 # Surface tension factor
    
    # Working area factor
    if A_ha >= 0.1 and A_ha <= 1.:
        F_HA = 1.
    elif A_ha >= 0.06:
        F_HA = 5.*A_ha + 0.5
    else:
        raise ValueError("ratio of open to active area, 'A', must be between 0.06 and 1") 
    
    return C_sbf * F_HA * F_ST * np.sqrt((rho_L-rho_V)/rho_V)

@njit(cache=True)
def compute_downcomer_area_fraction(F_LV):
    """
    Return the ratio of downcomer area to net (total) area, `A_dn`.
    
    Parameters
    ----------
    F_LV : float
        Flow parameter.

    Notes
    -----
    The fraction of downcomer area is given by [3]_. See source code for details.

    """
    if F_LV < 0.1:
        A_dn = 0.1
    elif F_LV < 1.:
        A_dn = 0.1 + (F_LV-0.1)/9.
    else:
        A_dn = 0.2
    return A_dn

@njit(cache=True)
def compute_tower_diameter(V_vol, U_f, f, A_dn):
    """Return tower diameter [D_T; in meter].
    
    Parameters
    ----------
    V_vol : float
        Vapor volumetric flow rate [m^3/s].
    U_f : float
        Maximum vapor velocity before flooding [m/s].
    f : float
        Ratio of actual velocity to `U_f`.
    A_dn : float
        Ratio of downcomer area to net (total) area.
    
    Notes
    -----
    The tower diameter is given by [3]_. See source code for details.
    
    """
    Di = np.sqrt(4.*V_vol/(f*U_f*np.pi*(1.-A_dn)))
    if Di < 0.914:
        # Make sure diameter is not too small
        Di = 0.914
    return Di

@njit(cache=True)
def compute_tower_height(TS, N_stages: int, top=True, bot=True):
    """
    Return the height of a tower [H; in meter].
    
    Parameters
    ----------
    TS : float
        Tray spacing [mm].
    N_stages : float
        Number of stages.
    
    Notes
    -----
    The tower height is given by [3]_. See source code for details.
    
    """
    # 3 m bottoms surge capacity, 1.25 m above top tray to remove entrained liquid
    H = TS*N_stages/1000.
    if top:
        H += 1.2672
    if bot:
        H += 3.
    return H 