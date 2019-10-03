# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 23:04:44 2019

@author: yoelr
"""
from .thermo_model import TDependentModel
from .model_handler import TDependentModelHandler
from .dippr import ModelEQ101
from .utils import sample_dict
from math import log, exp
import numpy as np
import os

folder = os.path.join(os.path.dirname(__file__), 'Vapor Pressure')

WagnerMcGarry = sample_dict(folder, "Wagner Original McGarry.tsv",
                            excluded=(0,))
Wagner = sample_dict(folder, "Wagner Collection Poling.tsv",
                     excluded=(0,))
Antoine = sample_dict(folder, "Antoine Collection Poling.tsv",
                      excluded=(0,))
AntoineExtended = sample_dict(folder, "Antoine Extended Collection Poling.tsv",
                              excluded=(0,))
Perrys2_8 = sample_dict(folder, "Table 2-8 Vapor Pressure of Inorganic and Organic Liquids.tsv",
                        excluded=(0,))
VDI_PPDS_3 = sample_dict(folder, "VDI PPDS Boiling temperatures at different pressures.tsv",
                         excluded=(0,))

def AntoineModel(A, B, C, Tmin, Tmax):
    r'''Return a TDependentModel object that calculates vapor pressure (Pa) of a chemical at arbitrary temperatures (K) using the Antoine equation. Parameters `A`, `B`, and `C` are chemical-dependent. Parameters can be found in numerous sources; however units of the coefficients used vary. Originally proposed by Antoine (1888) [2]_.

    .. math::
        \log_{\text{base}} P^{\text{sat}} = A - \frac{B}{T+C}

    Parameters
    ----------
    A, B, C : floats
        Regressed coefficients for Antoine equation for a chemical

    Notes
    -----
    Assumes coefficients are for calculating vapor pressure in Pascal. 
    Coefficients should be consistent with input temperatures in Kelvin;
    however, if both the given temperature and units are specific to degrees
    Celcius, the result will still be correct.
    
    References
    ----------
    .. [1] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
    .. [2] Antoine, C. 1888. Tensions des Vapeurs: Nouvelle Relation Entre les 
       Tensions et les Tempé. Compt.Rend. 107:681-684.
    .. [3] Yaws, Carl L. The Yaws Handbook of Vapor Pressure: Antoine 
       Coefficients. 1 edition. Houston, Tex: Gulf Publishing Company, 2007.
    '''
    def Antoine(T):
        return 10.0**(A - B / (T + C))
    return TDependentModel(Antoine, Tmin, Tmax)

def AntoineExtendedModel(Tc, to, A, B, C, n, E, F, Tmin, Tmax):
    r'''Return a TDependentModel object that calculates vapor pressure (Pa) of a chemical at arbitrary temperatures (K) using the TRC Extended Antoine equation. Parameters are chemical dependent, and said to be from the  Thermodynamics Research Center (TRC) at Texas A&M. Coefficients for various chemicals can be found in [1]_.

    .. math::
        \log_{10} P^{sat} = A - \frac{B}{T + C} + 0.43429x^n + Ex^8 + Fx^{12}
        
        x = \max \left(\frac{T-t_o-273.15}{T_c}, 0 \right)

    Parameters
    ----------
    A, B, C, n, E, F : floats
        Regressed coefficients for the Antoine Extended (TRC) equation,
        specific for each chemical, [-]
    
    Notes
    -----
    Assumes coefficients are for calculating vapor pressure in Pascal. 
    Coefficients should be consistent with input temperatures in Kelvin;

    References
    ----------
    .. [1] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
    
    '''
    max_ = max
    def Antoine_extended(T):
        x = max_((T - to - 273.15) / Tc, 0.0)
        return 10.0**(A - B / (T + C) + 0.43429 * x**n + E * x**8 + F * x**12)
    return TDependentModel(Antoine_extended, Tmin, Tmax)

def WagnerMcGarryModel(A, B, C, D, Tc, Pc, Tmin, Tmax):
    r'''Return a TDependentModel object that calculates vapor pressure (Pa) at arbitrary temperatures (K) using the Wagner equation (3, 6 form).

    .. math::
        \ln P^{sat}= \ln P_c + \frac{a\tau + b \tau^{1.5} + c\tau^3 + d\tau^6}
        {T_r}
        
        \tau = 1 - \frac{T}{T_c}

    Parameters
    ----------
    Tc : float
        Critical temperature, [K]
    Pc : float
        Critical pressure, [Pa]
    a, b, c, d : floats
        Parameters for wagner equation. Specific to each chemical. [-]

    References
    ----------
    .. [1] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
    .. [2] McGarry, Jack. "Correlation and Prediction of the Vapor Pressures of
       Pure Liquids over Large Pressure Ranges." Industrial & Engineering
       Chemistry Process Design and Development 22, no. 2 (April 1, 1983):
       313-22. doi:10.1021/i200021a023.
    '''
    def Wagner_McGarry(T):
        Tr = T / Tc
        tau = 1.0 - Tr
        return Pc * exp((A * tau + B * tau**1.5 + C * tau**3 + D * tau**6) / Tr)
    return TDependentModel(Wagner_McGarry, Tmin, Tmax)

def WagnerModel(Tc, Pc, A, B, C, D, Tmin, Tmax):
    r'''Return a TDependentModel object that calculates vapor pressure (Pa) at arbitrary temperatures (K) using the Wagner equation (2.5, 5 form).
    
    .. math::
        \ln P^{sat}= \ln P_c + \frac{a\tau + b \tau^{1.5} + c\tau^{2.5}
        + d\tau^5} {T_r}

        \tau = 1 - \frac{T}{T_c}

    Parameters
    ----------
    Tc : float
        Critical temperature, [K]
    Pc : float
        Critical pressure, [Pa]
    a, b, c, d : floats
        Parameters for wagner equation. Specific to each chemical. [-]

    References
    ----------
    .. [1] Wagner, W. "New Vapour Pressure Measurements for Argon and Nitrogen and
       a New Method for Establishing Rational Vapour Pressure Equations."
       Cryogenics 13, no. 8 (August 1973): 470-82. doi:10.1016/0011-2275(73)90003-9
    .. [2] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
    
    '''
    def Wagner(T):
        Tr = T / Tc
        τ = 1.0 - Tr
        return Pc * exp((A*τ + B*τ**1.5 + C*τ**2.5 + D*τ**5) / Tr)
    return TDependentModel(Wagner, Tmin, Tmax)

def BoilingCriticalRelationModel(Tb, Tc, Pc, Tmin, Tmax):
    r'''Return a TDependentModel object that calculates vapor pressure (Pa) of a fluid at arbitrary temperatures (K) using a CSP relationship as in [1]_.

    .. math::
        \ln P^{sat}_r = h\left( 1 - \frac{1}{T_r}\right)

        h = T_{br} \frac{\ln(P_c/101325)}{1-T_{br}}

    Parameters
    ----------
    Tb : float
        Boiling temperature of fluid [K]
    Tc : float
        Critical temperature of fluid [K]
    Pc : float
        Critical pressure of fluid [Pa]

    Notes
    -----
    Units are Pa. Formulation makes intuitive sense; a logarithmic form of
    interpolation.

    References
    ----------
    .. [1] Reid, Robert C..; Prausnitz, John M.;; Poling, Bruce E.
       The Properties of Gases and Liquids. McGraw-Hill Companies, 1987.
    '''
    Tbr = Tb / Tc
    h = Tbr * log(Pc / 101325.0) / (1 - Tbr)
    def Boiling_critical_relation(T):
        return exp(h * (1 - Tc / T)) * Pc
    return TDependentModel(Boiling_critical_relation, Tmin, Tmax)

def LeeKeslerModel(Tc, Pc, ω, Tmin, Tmax):
    r'''Return a TDependentModel object that calculates vapor pressure (Pa) of a fluid at arbitrary temperatures (K) using a CSP relationship by [1]_.

    .. math::
        \ln P^{sat}_r = f^{(0)} + \omega f^{(1)}

        f^{(0)} = 5.92714-\frac{6.09648}{T_r}-1.28862\ln T_r + 0.169347T_r^6

        f^{(1)} = 15.2518-\frac{15.6875}{T_r} - 13.4721 \ln T_r + 0.43577T_r^6

    Parameters
    ----------
    Tc : float
        Critical temperature of fluid [K]
    Pc : float
        Critical pressure of fluid [Pa]
    ω : float
        Acentric factor [-]

    Notes
    -----
    This equation appears in [1]_ in expanded form.
    The reduced pressure form of the equation ensures predicted vapor pressure 
    cannot surpass the critical pressure.

    References
    ----------
    .. [1] Lee, Byung Ik, and Michael G. Kesler. "A Generalized Thermodynamic
       Correlation Based on Three-Parameter Corresponding States." AIChE Journal
       21, no. 3 (1975): 510-527. doi:10.1002/aic.690210313.
    .. [2] Reid, Robert C..; Prausnitz, John M.;; Poling, Bruce E.
       The Properties of Gases and Liquids. McGraw-Hill Companies, 1987.
    '''
    def Lee_Kesler(T):
        Tr = T / Tc
        Tra = Tr**6
        logTr = log(Tr)
        f0 = 5.92714 - 6.09648 / Tr - 1.28862 * logTr + 0.169347 * Tra
        f1 = 15.2518 - 15.6875 / Tr - 13.4721 * logTr + 0.43577 * Tra
        return exp(f0 + ω * f1) * Pc
    return TDependentModel(Lee_Kesler, Tmin, Tmax)

def AmbroseWaltonModel(Tc, Pc, ω, Tmin, Tmax):
    r'''Return a TDependentModel object that calculates vapor pressure (Pa) of a fluid at arbitrary temperatures (K) using a CSP relationship by [1]_.

    .. math::
        \ln P_r=f^{(0)}+\omega f^{(1)}+\omega^2f^{(2)}

        f^{(0)}=\frac{-5.97616\tau + 1.29874\tau^{1.5}- 0.60394\tau^{2.5}
        -1.06841\tau^5}{T_r}

        f^{(1)}=\frac{-5.03365\tau + 1.11505\tau^{1.5}- 5.41217\tau^{2.5}
        -7.46628\tau^5}{T_r}

        f^{(2)}=\frac{-0.64771\tau + 2.41539\tau^{1.5}- 4.26979\tau^{2.5}
        +3.25259\tau^5}{T_r}

        \tau = 1-T_{r}

    Parameters
    ----------
    Tc : float
        Critical temperature of fluid [K]
    Pc : float
        Critical pressure of fluid [Pa]
    ω : float
        Acentric factor [-]

    Notes
    -----
    Somewhat more accurate than the Lee Kesler formulation.

    References
    ----------
    .. [1] Ambrose, D., and J. Walton. "Vapour Pressures up to Their Critical
       Temperatures of Normal Alkanes and 1-Alkanols." Pure and Applied
       Chemistry 61, no. 8 (1989): 1395-1403. doi:10.1351/pac198961081395.
    .. [2] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
    '''
    def Ambrose_Walton(T):
        Tr = T / Tc
        τ = 1 - Tr
        τa = τ**1.5
        τb = τ**2.5
        τc = τ**5
        f0 = -5.97616 * τ + 1.29874 * τa - 0.60394 * τb - 1.06841 * τc
        f1 = -5.03365 * τ + 1.11505 * τa - 5.41217 * τb - 7.46628 * τc
        f2 = -0.64771 * τ + 2.41539 * τa - 4.26979 * τb + 3.25259 * τc
        return Pc * exp((f0 + f1 * ω + f2 * ω**2) / Tr)
    return TDependentModel(Ambrose_Walton, Tmin, Tmax)

def SanjariModel(Tc, Pc, ω, Tmin, Tmax):
    r'''Return a TDependentModel object that calculates vapor pressure (Pa) of a fluid at arbitrary temperatures (K) using a CSP relationship by [1]_. Although developed for refrigerants, this model should have some general predictive ability.

    .. math::
        P^{sat} = P_c\exp(f^{(0)} + \omega f^{(1)} + \omega^2 f^{(2)})

        f^{(0)} = a_1 + \frac{a_2}{T_r} + a_3\ln T_r + a_4 T_r^{1.9}

        f^{(1)} = a_5 + \frac{a_6}{T_r} + a_7\ln T_r + a_8 T_r^{1.9}

        f^{(2)} = a_9 + \frac{a_{10}}{T_r} + a_{11}\ln T_r + a_{12} T_r^{1.9}

    Parameters
    ----------
    Tc : float
        Critical temperature of fluid [K]
    Pc : float
        Critical pressure of fluid [Pa]
    ω : float
        Acentric factor [-]

    Notes
    -----
    a[1-12] are as follows:
    6.83377, -5.76051, 0.90654, -1.16906,
    5.32034, -28.1460, -58.0352, 23.57466,
    18.19967, 16.33839, 65.6995, -35.9739.

    For a claimed fluid not included in the regression, R128, the claimed AARD
    was 0.428%. A re-calculation using 200 data points from 125.45 K to
    343.90225 K evenly spaced by 1.09775 K as generated by NIST Webbook April
    2016 produced an AARD of 0.644%. It is likely that the author's regression
    used more precision in its coefficients than was shown here. Nevertheless,
    the function is reproduced as shown in [1]_.

    For Tc=808 K, Pc=1100000 Pa, ω=1.1571, this function actually declines
    after 770 K.

    References
    ----------
    .. [1] Sanjari, Ehsan, Mehrdad Honarmand, Hamidreza Badihi, and Ali
       Ghaheri. "An Accurate Generalized Model for Predict Vapor Pressure of
       Refrigerants." International Journal of Refrigeration 36, no. 4
       (June 2013): 1327-32. doi:10.1016/j.ijrefrig.2013.01.007.
    '''
    def Sanjari(T):
        Tr = T / Tc
        logTr = log(Tr)
        Ta = Tr**1.9
        f0 = 6.83377 + -5.76051 / Tr + 0.90654 * logTr + -1.16906 * Ta
        f1 = 5.32034 + -28.1460 / Tr + -58.0352 * logTr + 23.57466 * Ta
        f2 = 18.19967 + 16.33839 / Tr + 65.6995 * logTr + -35.9739 * Ta
        return Pc * exp(f0 + f1 * ω + f2 * ω**2)
    return TDependentModel(Sanjari, Tmin, Tmax)

def EOSModel(EOS, Tmin, Tmax):
    r'''Return a TDependentModel object of the saturated pressure of an EOS model.'''
    def EOS(T): return EOS.Psat(T)
    return TDependentModel(EOS, Tmin, Tmax)

def EdalatModel(Tc, Pc, ω, Tmin, Tmax):
    r'''Return a TDependentModel object that calculates vapor pressure (Pa) of a fluid at arbitrary temperatures (K) using a CSP relationship by [1]_. Claimed to have a higher accuracy than the Lee-Kesler CSP relationship.

    .. math::
        \ln(P^{sat}/P_c) = \frac{a\tau + b\tau^{1.5} + c\tau^3 + d\tau^6}
        {1-\tau}
        
        a = -6.1559 - 4.0855\omega
        
        b = 1.5737 - 1.0540\omega - 4.4365\times 10^{-3} d
        
        c = -0.8747 - 7.8874\omega
        
        d = \frac{1}{-0.4893 - 0.9912\omega + 3.1551\omega^2}
        
        \tau = 1 - \frac{T}{T_c}
        
    Parameters
    ----------
    Tc : float
        Critical temperature of fluid [K]
    Pc : float
        Critical pressure of fluid [Pa]
    ω : float
        Acentric factor [-]

    Notes
    -----
    [1]_ found an average error of 6.06% on 94 compounds and 1106 data points.

    References
    ----------
    .. [1] Edalat, M., R. B. Bozar-Jomehri, and G. A. Mansoori. "Generalized 
       Equation Predicts Vapor Pressure of Hydrocarbons." Oil and Gas Journal; 
       91:5 (February 1, 1993).
    '''
    def Edalat(T):
        τ = 1.0 - T / Tc
        a = -6.1559 - 4.0855 * ω
        c = -0.8747 - 7.8874 * ω
        d = 1.0 / (-0.4893 - 0.9912 * ω + 3.1551 * ω**2)
        b = 1.5737 - 1.0540 * ω - 4.4365E-3 * d
        lnPr = (a * τ + b * τ**1.5 + c * τ**3.0 + d * τ**6.0) / (1.0 - τ)
        return exp(lnPr) * Pc
    return TDependentModel(Edalat, Tmin, Tmax)

def VaporPressure(CAS, Tb=None, Tc=None, Pc=None, omega=None, eos=None, models=None):
    if models is None:
        models = []
    if CAS in WagnerMcGarry:
        A, B, C, D, Pc, Tc, Tmin = WagnerMcGarry[CAS]
        Tmax = Tc
        models.append(WagnerMcGarryModel(A, B, C, D, Tc, Pc, Tmin, Tmax))
    if CAS in Wagner:
        A, B, C, D, Tc, Pc, Tmin, Tmax = Wagner[CAS]
        # Some Tmin values are missing; Arbitrary choice of 0.1 lower limit
        if np.isnan(Tmin):
            Tmin = Tmax * 0.1
        models.append(WagnerModel(A, B, C, D, Tc, Pc, Tmin, Tmax))
    if CAS in AntoineExtended:
        A, B, C, Tc, to, n, E, F, Tmin, Tmax = AntoineExtended[CAS]
        models.append(AntoineExtendedModel(A, B, C, Tc, to, n, E, F, Tmin, Tmax))
    if CAS in Antoine:
        A, B, C, Tmin, Tmax = Antoine[CAS]
        models.append(AntoineModel(A, B, C, Tmin, Tmax))
    if CAS in Perrys2_8:
        C1, C2, C3, C4, C5, Tmin, Tmax = Perrys2_8[CAS]
        models.append(ModelEQ101(C1, C2, C3, C4, C5, Tmin=Tmin, Tmax=Tmax))
    if CAS in VDI_PPDS_3:
        Tm, Tc, Pc, A, B, C, D = VDI_PPDS_3[CAS]
        models.append(WagnerModel(Tc, Pc, A, B, C, D, Tmin, Tmax))
    if all((Tb, Pc, Pc)):
        models.append(BoilingCriticalRelationModel(Tb, Tc, Pc))
    if all((Tc, Pc, omega)):
        models.append(LeeKeslerModel(Tc, Pc, omega))
        models.append(AmbroseWaltonModel(Tc, Pc, omega))
        models.append(SanjariModel(Tc, Pc, omega))
        models.append(EdalatModel(Tc, Pc, omega))
    if eos is not None:
        Tmin = 0.01
        Tmax = Tc
        if eos:
            models.append(EOSModel(eos, Tmin, Tmax))
    if not models:
        raise ValueError(f"no vapor pressure models available for CAS {CAS}")
    acting_model = models[0]
    return TDependentModelHandler(models, acting_model)