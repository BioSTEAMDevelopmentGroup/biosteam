# -*- coding: utf-8 -*-
"""
"""

__all__ = (
    'C_L',
    'C_O2_L',
    'kLa',
    'P_at_kLa',
    'log_mean_driving_force',
)

from math import exp, log

kLa_coefficients = {
    # Garcia-Ochoa, F.; Gomez, E. Bioreactor Scale-up and Oxygen Transfer Rate
    # in Microbial Processes: An Overview. Biotechnology Advances 2009, 27 (2),
    # 153–176. https://doi.org/10.1016/j.biotechadv.2008.10.006.
    # Name: (a, b, c) # kLa = a * (P/V) ^ b * U^c 
    'Figueiredo & Calderbank': (0.026, 0.6, 0.8),
    "Van't Riet": (0.026, 0.4, 0.5),
}

#: Henry's law coefficients. Data from NIST Standard Reference Database 69: NIST Chemistry WebBook:
#: https://webbook.nist.gov/cgi/cbook.cgi?ID=C7782447&Mask=10#Notes
H_coefficients = { # k_H, A
    'O2': (0.0013, 1500),
    'CO2': (0.035, 2400),
    'H2': (0.00078, 500), 
}

def Henrys_law_constant(T, k_H, A):
    """
    Return Henry's law constant [mol / kg / bar] for a chemical given the temperature 
    [K] and the coefficients.
    
    """
    dTinv = (1. / T - 1. / 298.15)
    return k_H * exp(A * dTinv)
    
def C_L(T, Py, chemical):
    """Chemical concentration [mol / kg] in the liquid given the temperature [K] and its
    partial pressure in the gas [bar]."""
    return Py * Henrys_law_constant(T, *H_coefficients[chemical])

def C_O2_L(T, P_O2): # Commonly used, so here for convinience
    """O2 concentration [mol / kg] in the liquid given the temperature [K] and oxygen
    partial pressure of O2 in the gas [bar]."""
    return P_O2 * Henrys_law_constant(T, *H_coefficients['O2'])

def kLa(P, V, U, coefficients=None):
    """
    Return the lumped mass transfer coefficient and the
    mean bubble specific interfacial area (k_L*a; 1/s) given the
    gassed power input (P; W), the total volume (V; m3), and the
    gas superficial velocity in the reactor (m/s).
    
    Other Parameters
    ----------------
    coefficients : Iterable[float]|str, optional
        Name of author or an iterable in the form of a, b, c. The
        
    Notes
    ------
    The correlation is kLa = a * (P/V) ^ b * U ^ c.
        
    Correlation from:
    Van’t Riet, K. Review of Measuring Methods and Results in Nonviscous 
    Gas-Liquid Mass Transfer in Stirred Vessels. Ind. Eng. Chem. Proc. Des. 
    Dev. 1979, 18 (3), 357–364. https://doi.org/10.1021/i260071a001.


    """
    if coefficients is None:
        a, b, c = kLa_coefficients['Figueiredo & Calderbank']
    elif isinstance(coefficients, str):
        a, b, c = kLa_coefficients[coefficients]
    else:
        a, b, c = coefficients
    return a * (P / V) ** b * U ** c

def P_at_kLa(kLa, V, U, coefficients=None):
    """
    Return the gassed power input (P; W) given the lumped mass transfer 
    coefficient and the mean bubble specific interfacial area (k_L*a; 1/s),
    the total volume (V; m3), and the gas superficial velocity in the reactor (m/s).
    
    Other Parameters
    ----------------
    coefficients : Iterable[float]|str, optional
        Name of author or an iterable in the form of a, b, c. The
        
    Notes
    ------
    The correlation is kLa = a * (P/V) ^ b * U ^ c.
    
    Correlation from:
    Van’t Riet, K. Review of Measuring Methods and Results in Nonviscous 
    Gas-Liquid Mass Transfer in Stirred Vessels. Ind. Eng. Chem. Proc. Des. 
    Dev. 1979, 18 (3), 357–364. https://doi.org/10.1021/i260071a001.

    """
    if coefficients is None:
        a, b, c = kLa_coefficients['Figueiredo & Calderbank']
    elif isinstance(coefficients, str):
        a, b, c = kLa_coefficients[coefficients]
    else:
        a, b, c = coefficients
    return (kLa / (a * U ** c)) ** (1 / b) * V

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