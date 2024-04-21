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

from math import exp, log, sqrt

kLa_coefficients = {
    # Garcia-Ochoa, F.; Gomez, E. Bioreactor Scale-up and Oxygen Transfer Rate
    # in Microbial Processes: An Overview. Biotechnology Advances 2009, 27 (2),
    # 153–176. https://doi.org/10.1016/j.biotechadv.2008.10.006.
    # Name: (a, b, c) # kLa = a * (P/V) ** b * U**c 
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
    The correlation is kLa = a * (P/V) ** b * U ** c.
        
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
    Take present that the units of the efficiency factor are [m**-1].
    This values are 
        - 0.016-0.028 m^-1 for coarse bubble diffusers
        - 0.07-0.1 m^-1 for fine bubble diffusers
    """
    K = None
    return K

    
def kla_bubcol_Deshpande(K, M_l, Rhat, T, rho_l, U_g):
    """
    Returns the KLa coefficient for a bubble column reactor based on the Deshpande et al. (2019) correlation.
    ------
    Parameters:
    K: float
        Efficiency change per unit change in liquid height above sparger, [m^-1]
    M_l: float
        Molar mass, e.g. water =  0.018 kg/mol, [kg/mol]
    Rhat: float
        Ideal gas constant, [J/mol-K]
    T: float
        Temperature, [K]
    rho_l: float
        Mass density of the liquid, e.g. 1000 kg/m^3 for water [kg/m^3]
    U_g: float
        Superficial gas velocity, [m/s]
    ----------
    References
    Deshpande, S. S., Kar, K., Pressler, J., Tebeka, I., Martins, B., Rosenfeld, D., & Biggs, J. (2019).
        Mass transfer estimation for bubble column scale up. Chemical Engineering Science, 205, 350–357.
        https://doi.org/10.1016/j.ces.2019.05.011
    """
    H_at_T = Henrys_law_constant(T)
    alpha = K * (H_at_T * M_l)/(Rhat * T * rho_l)
    beta = 1
    return alpha * U_g ** beta

def kla_bubcol_DeJesus(Q, mu, k, n, system_like = None, kargs = None):
    """
    Returns the KLa coefficient for a bubble column reactor based on the De Jesus et al. (2017) correlation.
    ------
    Parameters:
    Q: float
        Specific air flow rate, [vmm]
    mu: float
        Dynamic viscosity, [Pa.s]
    k: float
        Consistency index, [Pa.s^n]
    n: float
        Flow behavior index
    system_like: str
        Name of the system, e.g. "Glycerol" or "Xanthan"
    ----------
    References
    De Jesus, S. S., Moreira Neto, J., & Maciel Filho, R. (2017). 
        Hydrodynamics and mass transfer in bubble column, conventional airlift,
        stirred airlift and stirred tank bioreactors, using viscous fluid: A comparative study. 
        Biochemical Engineering Journal, 118, 70–81. https://doi.org/10.1016/j.bej.2016.11.019
    """
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
        if kargs is not None and system_like is None:
            raise ValueError("To use the De Jesus et al. (2017) correlation without\n a system_like, please provide the necessary kargs")
        else:
            if system_like == 'Newtonian':
                a = kargs['a']
                c = kargs['c']
                d = kargs['d']
                return a * Q ** c * mu ** d
            elif system_like == 'Non-Newtonian':
                e = kargs['e']
                h = kargs['h']
                i = kargs['i']
                j = kargs['j']
                return e * Q ** h * k ** i * n ** j
            else:
                raise ValueError("Please provide a valid system_like (in this case 'Newtonian' or 'Non-Newtonian')")

def kla_bubcol_Akita_Yoshida(D, mu_l, rho_l, D_l, g, sigma_l, epsilon_g):
    """
    Returns the KLa coefficient for a bubble column reactor based on the Akita & Yoshida (1973) correlation.

    ------
    Parameters:
    D: float
        Diameter of the column, [m]
    mu_l: float
        Viscosity of the liquid, [Cp]
    rho_l: float    
        Density of the liquid, [kg/m^3]
    D_l: float
        Diffusivity ofliquid, [m^2/s]
    g: float
        Gravitational acceleration, [m/s^2] 
    sigma_l: float
        The surface tension ofliquid [N/m]
    epsilon_g: float
        Overall gashold-up
    ------
    References:
    Akita, K., & Yoshida, F. (1973).
    Gas Holdup and Volumetric Mass Transfer Coefficient in Bubble Columns. 
    Effects of Liquid Properties. Industrial & Engineering Chemistry Process Design and Development, 12(1), 76–80. 
    https://doi.org/10.1021/i260045a015
    ------
    Notes:
    This correlation is valid under the following range of values:
    D = 0.15 m -> Diameter of the column
    H = 4 m -> Liquid height in the column
    V_g = 0 - 0.33 m/s -> Superficial gasvelocity [m/s]
    And the conditions are:
    ->  O2-Water, Aqueous Solutions of Glycerol, Glycol, Methanol, Sodium Sulphite
    """

    kla = (D_l/(D**2)) * 0.6 *(mu_l/(rho_l * D_l))**(0.5) * ((g*rho_l*D**2)/sigma_l)**(0.62) * ((g * D**3 * rho_l**2)/mu_l)**(0.31) * epsilon_g**(1.1)
    return kla


def kla_bubcol_Posarac_Tekic(D, mu_l, rho_l, D_l, g, sigma_l, epsilon_g):

    """
    Returns the KLa coefficient for a bubble column reactor based on the  Pošarac & Tekić (1987).

    ------
    Parameters:
    D: float
        Diameter of the column, [m]
    mu_l: float
        Viscosity of the liquid, [Cp]
    rho_l: float    
        Density of the liquid, [kg/m^3]
    D_l: float
        Diffusivity ofliquid, [m^2/s]
    g: float
        Gravitational acceleration, [m/s^2] 
    sigma_l: float
        The surface tension ofliquid [N/m]
    epsilon_g: float
        Overall gashold-up

    ---
    References:
    Pošarac, D., & Tekić, M. N. (1987). 
    Gas holdup and volumetric mass transfer coefficient in bubble columns with dilute alcohol solutions. 
    AIChE Journal, 33(3), 497–499. https://doi.org/10.1002/aic.690330316

    ------
    Notes:
    This correlation is valid under the following range of values:
    D = 0.1 m -> Diameter of the column
    H = 2.5 m -> Liquid height in the column
    V_g = 0.008 - 0.08 m/s -> Superficial gasvelocity [m/s]
    And the conditions are:
    ->  Air-Water, Aqueous Solution of Methanol, Ethanol, i-Propanol, nButanol
    """
    kla = (D_l/(D**2)) * 0.961 *(mu_l/(rho_l * D_l))**(0.5) * ((g*rho_l*D**2)/sigma_l)**(0.62) * ((g * D**3 * rho_l**2)/mu_l**2)**(0.31) * epsilon_g**(1.1)

    return kla


def kla_bubcol_Seno(D, mu_l, rho_l, D_l, g, sigma_l, epsilon_g, V_g, V_l):
    """
    Returns the KLa coefficient for a bubble column reactor based on the Seno et al. (1990) correlation.
    ------
    Parameters:
    D: float
        Diameter of the column, [m]
    mu_l: float
        Viscosity of the liquid, [Cp]
    rho_l: float    
        Density of the liquid, [kg/m^3]
    D_l: float
        Diffusivity ofliquid, [m^2/s]
    g: float
        Gravitational acceleration, [m/s^2] 
    sigma_l: float
        The surface tension ofliquid [N/m]
    epsilon_g: float
        Overall gashold-up
    V_g: float
        Superficial gas velocity, [m/s]
    V_l: float
        Superficial liquid velocity, [m/s]
    ------
    References:
    Seno, T., Uchida, S., & Tsuyutani, S. (1990).
        Mass transfer in countercurrent and cocurrent bubble columns. Chemical Engineering & Technology, 13(1), 113–118.
        https://doi.org/10.1002/ceat.270130115
    
    -----
    Notes:
    This correlation is valid under the following range of values:
    D = 0.0464 m -> Diameter of the column
    H = 1.36 m -> Liquid height in the column
    V_g = 0.005 - 0.04 m/s -> Superficial gasvelocity [m/s]
    V_l = 0.005 - 0.1 m/s -> Superficial liquid velocity [m/s]
    rho_l = 995 - 1043 kg/m^3 -> Density of the liquid [kg/m^3]
    mu_l = 0.653- 1.31 Cp -> Viscosity of the liquid [Cp]
    sigma_l = 0.0348-0.0728 N/m -> The surface tension ofliquid [N/m]
    D_l = 1.68-3.24 × 10^-9 m^2/s
    Conditions:
    O2- Water, Aqueous Solution of Butanol, Polyoxyethylene sorbitan monolaurate with silicone oil

    """
    kla = (D_l/(D**2)) * 0.6 * (V_g/(V_g+V_l))**-0.39 *(mu_l/(rho_l * D_l))**(0.5) * ((g * rho_l**2 * D**2)/sigma_l)**(0.62) * ((g * D**3 * rho_l**2)/mu_l**2)**(0.31) * epsilon_g**(1.1)
    return kla

def kla_bubcol_Suh(D, mu_l, rho_l, D_l, g, sigma_l, epsilon_g, V_g, V_l):
    """
    Returns the KLa coefficient for a bubble column reactor based on the Suh et al. (1991) correlation.
    -------
    Parameters:
    Parameters:
    D: float
        Diameter of the column, [m]
    mu_l: float
        Viscosity of the liquid, [Cp]
    rho_l: float    
        Density of the liquid, [kg/m^3]
    D_l: float
        Diffusivity ofliquid, [m^2/s]
    g: float
        Gravitational acceleration, [m/s^2] 
    sigma_l: float
        The surface tension ofliquid [N/m]
    epsilon_g: float
        Overall gashold-up
    V_g: float
        Superficial gas velocity, [m/s]
    V_l: float
        Superficial liquid velocity, [m/s]


    -------
    References:
    Suh, I. ‐S., Schumpe, A., Deckwer, W. ‐D., & Kulicke, W. ‐M. (1991). 
        Gas‐liquid mass transfer in the bubble column with viscoelastic liquid. 
        The Canadian Journal of Chemical Engineering, 69(2), 506–512. https://doi.org/10.1002/cjce.5450690215

    ------
    Notes:
    This correlation is valid under the following range of values:
    D: 0.15 m -> Diameter of the column [m]
    H: 2.9 m -> Liquid height in the column [m]
    V_g: 0.005 - 0.04 m/s -> Superficial gasvelocity [m/s]
    rho_l: 1001 - 1264 kg/m^3 -> Density of the liquid [kg/m^3]
    sigma_l: 0.0656-0.0746 N/m -> The surface tension of the liquid [N/m]
    D_l: 0.318-2.5 * 10^-9 m^2/s -> Diffusivity of liquid [m^2/s]

    And the conditions are:
    Air-Aqueous Sucrose Solution
    """

    kla = (D_l/(D**2)) * 0.018 * (mu_l/(rho_l * D_l))**(0.5) * ((g*rho_l*D**2)/sigma_l)**(0.75) * (g*D**3*rho_l**2/(mu_l**2))**(0.39)*(V_g/(sqrt(g*D)))

    return kla

def kla_bubcol_Shah(D, V_l, g, V_g, rho_l, mu_l):
    """
    Returns the KLa coefficient for a bubble column reactor based on the Shah et al. (2012) correlation.
    -------
    Parameters:
    D: float
        Diameter of the column, [m]
    V_l:
        Superficial liquid velocity, [m/s]
    g: 
        Gravitational acceleration, [m/s^2]
    V_g:
        Superficial gas velocity, [m/s]
    rho_l:
        Density of the liquid, [kg/m^3]
    mu_l:   
        Viscosity of the liquid, [Cp]
    -------
    References:
    Shah, M., Kiss, A. A., Zondervan, E., Van Der Schaaf, J., & De Haan, A. B. (2012).
        Gas Holdup, Axial Dispersion, and Mass Transfer Studies in Bubble Columns. 
        Industrial & Engineering Chemistry Research, 51(43), 14268–14278. https://doi.org/10.1021/ie301227t
    -------
    Notes:
    This correlation is valid under the following range of values:
    D: 0.29 m -> Diameter of the column [m]
    H: 0.2 m -> Liquid height in the column [m]
    V_g: 0.021 - 0.105 m/s -> Superficial gasvelocity [m/s]
    V_l: 0.0005 - 0.002 m/s -> Superficial liquid velocity [m/s]
    mu_l: 1 - 50 Cp -> Viscosity of the liquid [Cp]
    """

    kla = (V_l/D) * 0.0029 * (V_g**2 / (g * D))**0.301 * (V_l**2/(g * D))**-0.511 * (g*rho_l**3 * D**3 / (mu_l**2))**0.12
    return kla

def kla_stirred_Labik(N, D, v_s, P0):
    """
    Returns the KLa coefficient for a stirred tank reactor based on the Labik et al. (2017) correlation.
    ------
    Parameters:
    N: float
        impeller frequency [s^-1]
    D: float
        impeller diameter [m]
    v_s: float
        gas superficial velocity [m/s]
    P0: float
        impeller power number (P_u/rhh N^3 D^5) [-]
    ------
    References:
    Labík, L., Moucha, T., Petříček, R., Rejl, J. F., Valenz, L., & Haidl, J. (2017). 
        Volumetric mass transfer coefficient in viscous liquid in mechanically agitated fermenters. 
        Measurement and correlation. Chemical Engineering Science, 170, 451–463. https://doi.org/10.1016/j.ces.2017.04.006
    """
    kla = 0.295 * ( N * D)**2.083 * v_s**0.461 * P0 ** 0.737
    return kla




def kla_stirred_Galaction(aeration_type: str, V, C_x, v_s, P_a = None, P = None, organism_type = None):
    """
    Return the kla for a stirred tank reactor based on the Galaction et al. (2004) correlation.
    ------
    Parameters:
    aeration_type: str
        Type of aireation, in this case either 'Surface' or 'Submerged'
    P_a: float
        power consumption for mixing of aerated broths, [W]
    P: float
        power consumption for mixing of non-aerated broths, [W]
    V: float
        volume of the medium, [m**3]
    C_x: float
        Biomass concentration, [g/l dry weight]
    v_s: float
        superficial air velocity, [m/s]
    ---
    References:
    Galaction, A.-I., Cascaval, D., Oniscu, C., & Turnea, M. (2004). 
        Prediction of oxygen mass transfer coefficients in stirred bioreactors for bacteria, yeasts and fungus broths.
        Biochemical Engineering Journal, 20(1), 85–94. https://doi.org/10.1016/j.bej.2004.02.005
    """
    def kla_submerged(alpha, beta, gamma, delta, C_x, P_a, V, v_s):
        """
        Auxiliary function to calculate kla for submerged aeration.
        """
        return alpha * C_x ** beta * (P_a / V) ** gamma * v_s ** delta
    def kla_surface(alpha, beta, P, V):
        """
        Auxiliary function to calculate kla for surface aeration.
        """
        return 1/exp(alpha) * (P/V) ** beta
    
    if aeration_type == 'Submerged':
        
        if organism_type == 'P.shermanii':
            alpha = 6.586
            beta = -0.282
            gamma = 0.0286
            delta = 0.429
            return kla_submerged(alpha, beta, gamma, delta, C_x, P_a, V, v_s)
        elif organism_type == 'S.cerevisiae':
            alpha = 52.44
            beta = -0.702
            gamma = -0.0762
            delta = 0.514
            return kla_submerged(alpha, beta, gamma, delta, C_x, P_a, V, v_s)
        elif organism_type == 'P.chrysogenum_pellet':
            alpha = 0.193
            beta = - 0.269
            gamma = 0.0288
            delta = 0.257
            return kla_submerged(alpha, beta, gamma, delta, C_x, P_a, V, v_s)
        elif organism_type == 'P.chrysogenum_mycelia':
            alpha = 33.59
            beta = -1.012
            gamma = -0.0463
            delta = 0.94
            return kla_submerged(alpha, beta, gamma, delta, C_x, P_a, V, v_s)
        
        else:
            raise ValueError("Please provide a valid organism_type (in this case 'P.shermanii', 'S.cerevisiae', 'P.chrysogenum_pellet' or 'P.chrysogenum_mycelia')")

    elif aeration_type == 'Surface':
        if organism_type == 'P.shermanii':
            alpha = 10.35 + 0.63 *log(C_x)
            beta = 0.694 - 0.069 * log(C_x)
            return kla_surface(alpha, beta, P, V)
        elif organism_type == 'S.cerevisiae':
            alpha = 11.64 + 0.307 * log(C_x)
            beta = 1.14 - 0.187 * log(C_x)
            return kla_surface(alpha, beta, P, V)
        elif organism_type == 'P.chrysogenum_pellet':
            alpha = 12.77 + 0.227 * log(C_x)
            beta = 0.62 - 0.109 * log(C_x)
            return kla_surface(alpha, beta, P, V)
        elif organism_type == 'P.chrysogenum_mycelia':
            alpha = 11.46 + 0.729 * log(C_x)
            beta = 0.093 + 0.057 * log(C_x)
            return kla_surface(alpha, beta, P, V)
        else: 
            raise ValueError("Please provide a valid organism_type (in this case 'P.shermanii', 'S.cerevisiae', 'P.chrysogenum_pellet' or 'P.chrysogenum_mycelia')")
    else:
        raise ValueError("Please provide a valid aeration_type (in this case 'Surface' or 'Submerged')")




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
    The correlation is kLa = a * (P/V) ** b * U ** c.
    
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