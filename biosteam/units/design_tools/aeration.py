# -*- coding: utf-8 -*-
"""
"""

__all__ = (
    'H_O2',
    'C_O2_L',
    'ka_L'
)

from math import exp

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

