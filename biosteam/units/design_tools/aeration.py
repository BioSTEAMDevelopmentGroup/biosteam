# -*- coding: utf-8 -*-
"""
"""

__all__ = (
    'H_O2',
    'C_O2_L',
    'ka_L',
    'log_mean_driving_force',
)

from math import exp, log

def H_O2(T): 
    """
    Return Henry's law constant [mol O2 / kg / bar] for O2 given the temperature [K].
    
    Data from NIST Standard Reference Database 69: NIST Chemistry WebBook:
    https://webbook.nist.gov/cgi/cbook.cgi?ID=C7782447&Mask=10#Notes
    
    """
    k_H = 0.0013
    A = 1500
    dTinv = (1. / T - 1. / 298.15)
    return k_H * exp(A * dTinv)

def C_O2_L(T, P_O2): 
    """O2 concentration [mol / kg] in the liquid given the temperature [K] and oxygen
    partial pressure of O2 in the gas [bar]."""
    return P_O2 * H_O2(T)

def ka_L(P, V, U):
    """
    Return the lumped mass transfer coefficient (k_L) and the
    mean bubble specific interfacial area (a) given the
    gassed power input (P; W), the total volume (V; m3), and the
    gas superficial velocity in the reactor (m/s).
    
    Correlation from:
    Gregory T. Benz. Bioreactor Design for Chemical Engineers. Benz Technology International, Inc.
    https://www.academia.edu/19636928/Bioreactor_Design_for_Chemical_Engineers

    """
    return 0.002 * (P / V) ** 0.7 * U ** 0.2

def log_mean_driving_force(C_sat_out, C_sat_in, C_out, C_in=None):
    """
    Return the driving force for mass transfer. In small vessels (<1 m tall) 
    where both liquid concentration and saturation are almost constant, the simple
    form is adequate [1]_. In tall vessels, the log-mean driving force should be 
    used for more accuracy, since both the local concentration and the saturation
    concentration are different in the top and bottom of a bioreactor [1]_. 
    
    Parameters
    ----------
    C_sat_out : float
        Saturated concentration entering the bioreactor.
    C_sat_out : float
        Saturated concentration exiting the bioreactor.
    C_out : float, optional
        Outlet concentration. 
    C_in : float
        Inlet concentration. Defaults to the outlet concentration, which
        assumes perfect mixing and oxygen uptake rate is the same everywhere.
        
    References
    ----------
    [1] Benz, G. T. Bioreactor Design for Chemical Engineers. AICHE CEP 2011, 21–26.

    """
    if C_in is None: C_in = C_out # Assume perfect mixing and oxygen uptake rate is the same everywhere.
    dC_out = C_sat_out - C_out
    dC_in = C_sat_in - C_in
    if dC_out < 1e-9: # This is not theoretically possible, but may be an assumption
        dC_out = 1e-9 # Assume it near zero
    return (dC_out - dC_in) / log(dC_out / dC_in)