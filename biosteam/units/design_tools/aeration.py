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

def get_efficiency_factor():
    """
    Return the efficiency factor for mass transfer in
    bubble column reactors.    
    Take present that the units of the efficiency factor are [m^-1].
    This values are 
        - 0.016-0.028 m^-1 for coarse bubble diffusers
        - 0.07-0.1 m^-1 for fine bubble diffusers
    """
    K = None
    return K

def kLa_corr_bubcol(U_g, Q, mu, k, n,correlation = "Deshpande", 
                    system_like = None):
    """
    Returns the kLa coefficient given a certain correlation used, these are for bubble columns. 
    Various correlations are available in literature. 
    Mainly we are going to cover the correlations from the following sources:
    
    *correlation = "Deshpande"
    -  Deshpande, S. S., Kar, K., Pressler, J., Tebeka, I., Martins, B., Rosenfeld, D., & Biggs, J. (2019).
        Mass transfer estimation for bubble column scale up. Chemical Engineering Science, 205, 350–357.
        https://doi.org/10.1016/j.ces.2019.05.011
        This follows a power law correlation, which is fairly used in the literature. The correlation is given by:
            kla = alpha U_g ^ beta
        parameters:
        - U_g = superficial gas velocity, [m/s]
        - alpha = K (H(T) * M_l)/(Rhat T * rho_l) 
            - K = efficiency change per unit change in liquid height above sparger, [m^-1]
            - H(T)  = Henry’s law constant for oxygen in water @ T (in K) [Pa]
            - M_l = molar mass of water (0.018 kg/mol), [kg/mol]
            - Rhat = ideal gas constant (8.314 J/mol-K), [J/mol-K]
            - T = temperature [K]
            - rho_l = mass density of water (1000 kg/m3), [kg/m3]
        - beta = 1    
    *correlation ="DeJesus"
    - De Jesus, S. S., Moreira Neto, J., & Maciel Filho, R. (2017). 
        Hydrodynamics and mass transfer in bubble column, conventional airlift,
        stirred airlift and stirred tank bioreactors, using viscous fluid: A comparative study. 
        Biochemical Engineering Journal, 118, 70–81. https://doi.org/10.1016/j.bej.2016.11.019

        This correlation is based on:
            kla = a Q ^ c * mu ^ d  -> for Newtonian fluids
            kla = e Q ^ h * k ^ i * n ^ j -> for non-Newtonian fluids
        
        parameters:
        - Q = specific air flow rate [vmm]
        - mu = Dynamic viscosity [Pa.s]
        - k = Consistency index [Pa.s^n]
        - gamma_AV = Average Shear rate [s^-1]
        - n = flow behavior index
        - a, c, d = constants

        In this case, the case studies give the following values that might be used in the case 
        in which the fluid behaves similarly:
        - Glycerol:
            - a = 1.27e-2
            - c = 0.51
            - d = -0.12
        - Xanthan:
            - e = 5.80e-3
            - h = 0.66
            - i = -0.50
            - j = -0.80
    """
    if correlation == "Deshpande":
        K = get_efficiency_factor()
        H_at_293 = 3.9e-9 #[Pa]
        M_l = 0.018 #[kg/mol]
        Rhat = 8.314 #[J/mol-K]
        T = 293.15 #[K]
        rho_l = 1000
        alpha = K * (H_at_293 * M_l)/(Rhat * T * rho_l)
        beta = 1
        return alpha * U_g ** beta
    elif correlation == "DeJesus":
        if system_like == "Glycerol":
            a = 1.27e-2
            c = 0.51
            d = -0.12
            return a * Q ** c * mu ** d
        elif system_like == "Xanthan":
            e = 5.80e-3
            h = 0.66
            i = -0.50
            j = -0.80
            return e * Q ** h * k ** i * n ** j


    
    else:
        raise ValueError(f"Correlation {correlation} not available. Please choose a valid correlation")
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