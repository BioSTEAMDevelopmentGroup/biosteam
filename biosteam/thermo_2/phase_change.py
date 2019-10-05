# -*- coding: utf-8 -*-
'''Chemical Engineering Design Library (ChEDL). Utilities for process modeling.
Copyright (C) 2016, Caleb Bell <Caleb.Andrew.Bell@gmail.com>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.'''

from __future__ import division

__all__ = ['Yaws_data', 'Tb_methods', 'Tb', 'Tm_ON_data', 'Tm_methods', 'Tm', 
           'Clapeyron', 'Pitzer', 'SMK', 'MK', 'Velasco', 'Riedel', 'Chen', 
           'Liu', 'Vetere', 'GharagheiziHvap_data', 'CRCHvap_data', 
           'Perrys2_150', 'VDI_PPDS_4', 'Alibakhshi_Cs', 'Watson', 
           'enthalpy_vaporization_methods', 'EnthalpyVaporization', 
           'CRCHfus_data', 'Hfus', 'GharagheiziHsub_data', 'Hsub', 'Tliquidus']

import os
import numpy as np
import pandas as pd
from numba import njit
from .thermo_model import TDependentModel, InterpolatedTDependentModel
from .model_handler import TDependentModelHandler
from .utils import log
from .utils import R, pi, N_A
from .miscdata import CRC_organic_data, CRC_inorganic_data
from .miscdata import _VDISaturationDict, VDI_tabular_data
from .utils import  mixing_simple, none_and_length_check
from .vapor_pressure import VaporPressure
from .coolprop import has_CoolProp, PropsSI, coolprop_dict, coolprop_fluids
from .dippr import ModelEQ106

folder = os.path.join(os.path.dirname(__file__), 'Phase Change')


Yaws_data = pd.read_csv(os.path.join(folder,
'Yaws Boiling Points.tsv'), sep='\t', index_col=0)

Tm_ON_data = pd.read_csv(os.path.join(folder, 'OpenNotebook Melting Points.tsv'),
                         sep='\t', index_col=0)

GharagheiziHvap_data = pd.read_csv(os.path.join(folder, 'Ghazerati Appendix Vaporization Enthalpy.tsv'),
                                   sep='\t', index_col=0)

CRCHvap_data = pd.read_csv(os.path.join(folder, 'CRC Handbook Heat of Vaporization.tsv'),
                           sep='\t', index_col=0)

CRCHfus_data = pd.read_csv(os.path.join(folder, 'CRC Handbook Heat of Fusion.tsv'),
                                    sep='\t', index_col=0)

GharagheiziHsub_data = pd.read_csv(os.path.join(folder, 'Ghazerati Appendix Sublimation Enthalpy.tsv'),
                                    sep='\t', index_col=0)

Perrys2_150 = pd.read_csv(os.path.join(folder, 'Table 2-150 Heats of Vaporization of Inorganic and Organic Liquids.tsv'),
                          sep='\t', index_col=0)
_Perrys2_150_values = Perrys2_150.values

VDI_PPDS_4 = pd.read_csv(os.path.join(folder, 'VDI PPDS Enthalpies of vaporization.tsv'),
                          sep='\t', index_col=0)
_VDI_PPDS_4_values = VDI_PPDS_4.values

Alibakhshi_Cs = pd.read_csv(os.path.join(folder, 'Alibakhshi one-coefficient enthalpy of vaporization.tsv'),
                          sep='\t', index_col=0)


### Boiling Point at 1 atm

CRC_ORG = 'CRC_ORG'
CRC_INORG = 'CRC_INORG'
YAWS = 'YAWS'
PSAT_DEFINITION = 'PSAT_DEFINITION'
NONE = 'NONE'

Tb_methods = [CRC_INORG, CRC_ORG, YAWS, PSAT_DEFINITION]


def Tb(CASRN, AvailableMethods=False, Method=None, IgnoreMethods=[PSAT_DEFINITION]):
    r'''This function handles the retrieval of a chemical's boiling
    point. Lookup is based on CASRNs. Will automatically select a data
    source to use if no Method is provided; returns None if the data is not
    available.

    Prefered sources are 'CRC Physical Constants, organic' for organic
    chemicals, and 'CRC Physical Constants, inorganic' for inorganic
    chemicals. Function has data for approximately 13000 chemicals.

    Parameters
    ----------
    CASRN : string
        CASRN [-]

    Returns
    -------
    Tb : float
        Boiling temperature, [K]
    methods : list, only returned if AvailableMethods == True
        List of methods which can be used to obtain Tb with the given inputs

    Other Parameters
    ----------------
    Method : string, optional
        A string for the method name to use, as defined by constants in
        Tb_methods
    AvailableMethods : bool, optional
        If True, function will determine which methods can be used to obtain
        Tb for the desired chemical, and will return methods instead of Tb
    IgnoreMethods : list, optional
        A list of methods to ignore in obtaining the full list of methods,
        useful for for performance reasons and ignoring inaccurate methods

    Notes
    -----
    A total of four methods are available for this function. They are:

        * 'CRC_ORG', a compillation of data on organics
          as published in [1]_.
        * 'CRC_INORG', a compillation of data on
          inorganic as published in [1]_.
        * 'YAWS', a large compillation of data from a
          variety of sources; no data points are sourced in the work of [2]_.
        * 'PSAT_DEFINITION', calculation of boiling point from a
          vapor pressure calculation. This is normally off by a fraction of a
          degree even in the best cases. Listed in IgnoreMethods by default
          for performance reasons.

    Examples
    --------
    >>> Tb('7732-18-5')
    373.124

    References
    ----------
    .. [1] Haynes, W.M., Thomas J. Bruno, and David R. Lide. CRC Handbook of
       Chemistry and Physics, 95E. Boca Raton, FL: CRC press, 2014.
    .. [2] Yaws, Carl L. Thermophysical Properties of Chemicals and
       Hydrocarbons, Second Edition. Amsterdam Boston: Gulf Professional
       Publishing, 2014.
    '''
    def list_methods():
        methods = []
        if CASRN in CRC_inorganic_data.index and not np.isnan(CRC_inorganic_data.at[CASRN, 'Tb']):
            methods.append(CRC_INORG)
        if CASRN in CRC_organic_data.index and not np.isnan(CRC_organic_data.at[CASRN, 'Tb']):
            methods.append(CRC_ORG)
        if CASRN in Yaws_data.index:
            methods.append(YAWS)
        if PSAT_DEFINITION not in IgnoreMethods:
            try:
                # For some chemicals, vapor pressure range will exclude Tb
                VaporPressure(CASRN=CASRN).solve_prop(101325.)
                methods.append(PSAT_DEFINITION)
            except:  # pragma: no cover
                pass
        if IgnoreMethods:
            for Method in IgnoreMethods:
                if Method in methods:
                    methods.remove(Method)
        methods.append(NONE)
        return methods
    if AvailableMethods:
        return list_methods()
    if not Method:
        Method = list_methods()[0]

    if Method == CRC_INORG:
        return float(CRC_inorganic_data.at[CASRN, 'Tb'])
    elif Method == CRC_ORG:
        return float(CRC_organic_data.at[CASRN, 'Tb'])
    elif Method == YAWS:
        return float(Yaws_data.at[CASRN, 'Tb'])
    elif Method == PSAT_DEFINITION:
        return VaporPressure(CASRN=CASRN).solve_prop(101325.)
    elif Method == NONE:
        return None
    else:
        raise Exception('Failure in in function')


### Melting Point


OPEN_NTBKM = 'OPEN_NTBKM'

Tm_methods = [OPEN_NTBKM, CRC_INORG, CRC_ORG]


def Tm(CASRN, AvailableMethods=False, Method=None, IgnoreMethods=[]):
    r'''This function handles the retrieval of a chemical's melting
    point. Lookup is based on CASRNs. Will automatically select a data
    source to use if no Method is provided; returns None if the data is not
    available.

    Prefered sources are 'Open Notebook Melting Points', with backup sources
    'CRC Physical Constants, organic' for organic chemicals, and
    'CRC Physical Constants, inorganic' for inorganic chemicals. Function has
    data for approximately 14000 chemicals.

    Parameters
    ----------
    CASRN : string
        CASRN [-]

    Returns
    -------
    Tm : float
        Melting temperature, [K]
    methods : list, only returned if AvailableMethods == True
        List of methods which can be used to obtain Tm with the given inputs

    Other Parameters
    ----------------
    Method : string, optional
        A string for the method name to use, as defined by constants in
        Tm_methods
    AvailableMethods : bool, optional
        If True, function will determine which methods can be used to obtain
        Tm for the desired chemical, and will return methods instead of Tm
    IgnoreMethods : list, optional
        A list of methods to ignore in obtaining the full list of methods

    Notes
    -----
    A total of three sources are available for this function. They are:

        * 'OPEN_NTBKM, a compillation of data on organics
          as published in [1]_ as Open Notebook Melting Points; Averaged 
          (median) values were used when
          multiple points were available. For more information on this
          invaluable and excellent collection, see
          http://onswebservices.wikispaces.com/meltingpoint.
        * 'CRC_ORG', a compillation of data on organics
          as published in [2]_.
        * 'CRC_INORG', a compillation of data on
          inorganic as published in [2]_.

    Examples
    --------
    >>> Tm(CASRN='7732-18-5')
    273.15

    References
    ----------
    .. [1] Bradley, Jean-Claude, Antony Williams, and Andrew Lang.
       "Jean-Claude Bradley Open Melting Point Dataset", May 20, 2014.
       https://figshare.com/articles/Jean_Claude_Bradley_Open_Melting_Point_Datset/1031637.
    .. [2] Haynes, W.M., Thomas J. Bruno, and David R. Lide. CRC Handbook of
       Chemistry and Physics, 95E. Boca Raton, FL: CRC press, 2014.
    '''
    def list_methods():
        methods = []
        if CASRN in Tm_ON_data.index:
            methods.append(OPEN_NTBKM)
        if CASRN in CRC_inorganic_data.index and not np.isnan(CRC_inorganic_data.at[CASRN, 'Tm']):
            methods.append(CRC_INORG)
        if CASRN in CRC_organic_data.index and not np.isnan(CRC_organic_data.at[CASRN, 'Tm']):
            methods.append(CRC_ORG)
        if IgnoreMethods:
            for Method in IgnoreMethods:
                if Method in methods:
                    methods.remove(Method)
        methods.append(NONE)
        return methods
    if AvailableMethods:
        return list_methods()
    if not Method:
        Method = list_methods()[0]

    if Method == OPEN_NTBKM:
        return float(Tm_ON_data.at[CASRN, 'Tm'])
    elif Method == CRC_INORG:
        return float(CRC_inorganic_data.at[CASRN, 'Tm'])
    elif Method == CRC_ORG:
        return float(CRC_organic_data.at[CASRN, 'Tm'])
    elif Method == NONE:
        return None
    else:
        raise Exception('Failure in in function')


### Enthalpy of Vaporization at T

def Clapeyron(T, Tc, Pc, dZ, Psat):
    r'''Calculates enthalpy of vaporization at arbitrary temperatures using the
    Clapeyron equation.

    The enthalpy of vaporization is given by:

    .. math::
        \Delta H_{vap} = RT \Delta Z \frac{\ln (P_c/Psat)}{(1-T_{r})}

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    Tc : float
        Critical temperature of fluid [K]
    Pc : float
        Critical pressure of fluid [Pa]
    dZ : function(T)
        Change in compressibility factor between liquid and gas, []
    Psat : function(T)
        Saturation pressure of fluid [Pa], optional

    Returns
    -------
    Hvap : float
        Enthalpy of vaporization, [J/mol]

    Notes
    -----
    No original source is available for this equation.
    [1]_ claims this equation overpredicts enthalpy by several percent.
    Under Tr = 0.8, dZ = 1 is a reasonable assumption.
    This equation is most accurate at the normal boiling point.

    Internal units are bar.

    WARNING: I believe it possible that the adjustment for pressure may be incorrect

    Examples
    --------
    Problem from Perry's examples.

    >>> Clapeyron(T=294.0, Tc=466.0, Pc=5.55E6)
    26512.354585061985

    References
    ----------
    .. [1] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
    '''
    Tr = T/Tc
    return R*T*dZ(T)*log(Pc/Psat(T))/(1. - Tr)

@njit
def Pitzer(T, Tc, omega):
    r'''Calculates enthalpy of vaporization at arbitrary temperatures using a
    fit by [2]_ to the work of Pitzer [1]_; requires a chemical's critical
    temperature and acentric factor.

    The enthalpy of vaporization is given by:

    .. math::
        \frac{\Delta_{vap} H}{RT_c}=7.08(1-T_r)^{0.354}+10.95\omega(1-T_r)^{0.456}

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    Tc : float
        Critical temperature of fluid [K]
    omega : float
        Acentric factor [-]

    Returns
    -------
    Hvap : float
        Enthalpy of vaporization, [J/mol]

    Notes
    -----
    This equation is listed in [3]_, page 2-487 as method #2 for estimating
    Hvap. This cites [2]_.

    The recommended range is 0.6 to 1 Tr. Users should expect up to 5% error.
    T must be under Tc, or an exception is raised.

    The original article has been reviewed and found to have a set of tabulated
    values which could be used instead of the fit function to provide additional
    accuracy.

    Examples
    --------
    Example as in [3]_, p2-487; exp: 37.51 kJ/mol

    >>> Pitzer(452, 645.6, 0.35017)
    36696.736640106414

    References
    ----------
    .. [1] Pitzer, Kenneth S. "The Volumetric and Thermodynamic Properties of
       Fluids. I. Theoretical Basis and Virial Coefficients."
       Journal of the American Chemical Society 77, no. 13 (July 1, 1955):
       3427-33. doi:10.1021/ja01618a001
    .. [2] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
    .. [3] Green, Don, and Robert Perry. Perry's Chemical Engineers' Handbook,
       Eighth Edition. McGraw-Hill Professional, 2007.
    '''
    Tr = T/Tc
    return R*Tc * (7.08*(1. - Tr)**0.354 + 10.95*omega*(1. - Tr)**0.456)

@njit
def SMK(T, Tc, omega):
    r'''Calculates enthalpy of vaporization at arbitrary temperatures using a
    the work of [1]_; requires a chemical's critical temperature and
    acentric factor.

    The enthalpy of vaporization is given by:

    .. math::
         \frac{\Delta H_{vap}} {RT_c} =
         \left( \frac{\Delta H_{vap}} {RT_c} \right)^{(R1)} + \left(
         \frac{\omega - \omega^{(R1)}} {\omega^{(R2)} - \omega^{(R1)}} \right)
         \left[\left( \frac{\Delta H_{vap}} {RT_c} \right)^{(R2)} - \left(
         \frac{\Delta H_{vap}} {RT_c} \right)^{(R1)} \right]

        \left( \frac{\Delta H_{vap}} {RT_c} \right)^{(R1)}
        = 6.537 \tau^{1/3} - 2.467 \tau^{5/6} - 77.251 \tau^{1.208} +
        59.634 \tau + 36.009 \tau^2 - 14.606 \tau^3

        \left( \frac{\Delta H_{vap}} {RT_c} \right)^{(R2)} - \left(
        \frac{\Delta H_{vap}} {RT_c} \right)^{(R1)}=-0.133 \tau^{1/3} - 28.215
        \tau^{5/6} - 82.958 \tau^{1.208} + 99.00 \tau  + 19.105 \tau^2 -2.796 \tau^3

        \tau = 1-T/T_c

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    Tc : float
        Critical temperature of fluid [K]
    omega : float
        Acentric factor [-]

    Returns
    -------
    Hvap : float
        Enthalpy of vaporization, [J/mol]

    Notes
    -----
    The original article has been reviewed and found to have coefficients with
    slightly more precision. Additionally, the form of the equation is slightly
    different, but numerically equivalent.

    The refence fluids are:

    :math:`\omega_0` = benzene = 0.212

    :math:`\omega_1` = carbazole = 0.461

    A sample problem in the article has been verified. The numerical result
    presented by the author requires high numerical accuracy to obtain.

    Examples
    --------
    Problem in [1]_:

    >>> SMK(553.15, 751.35, 0.302)
    39866.17647797959

    References
    ----------
    .. [1] Sivaraman, Alwarappa, Joe W. Magee, and Riki Kobayashi. "Generalized
       Correlation of Latent Heats of Vaporization of Coal-Liquid Model Compounds
       between Their Freezing Points and Critical Points." Industrial &
       Engineering Chemistry Fundamentals 23, no. 1 (February 1, 1984): 97-100.
       doi:10.1021/i100013a017.
    '''
    omegaR1, omegaR2 = 0.212, 0.461
    A10 = 6.536924
    A20 = -2.466698
    A30 = -77.52141
    B10 = 59.63435
    B20 = 36.09887
    B30 = -14.60567

    A11 = -0.132584
    A21 = -28.21525
    A31 = -82.95820
    B11 = 99.00008
    B21 = 19.10458
    B31 = -2.795660

    tau = 1. - T/Tc
    L0 = A10*tau**(1/3.) + A20*tau**(5/6.) + A30*tau**(1-1/8. + 1/3.) + \
        B10*tau + B20*tau**2 + B30*tau**3

    L1 = A11*tau**(1/3.) + A21*tau**(5/6.0) + A31*tau**(1-1/8. + 1/3.) + \
        B11*tau + B21*tau**2 + B31*tau**3

    domega = (omega - omegaR1)/(omegaR2 - omegaR1)
    return R*Tc*(L0 + domega*L1)

@njit
def MK(T, Tc, omega):
    r'''Calculates enthalpy of vaporization at arbitrary temperatures using a
    the work of [1]_; requires a chemical's critical temperature and
    acentric factor.

    The enthalpy of vaporization is given by:

    .. math::
        \Delta H_{vap} =  \Delta H_{vap}^{(0)} + \omega \Delta H_{vap}^{(1)} + \omega^2 \Delta H_{vap}^{(2)}

        \frac{\Delta H_{vap}^{(i)}}{RT_c} = b^{(j)} \tau^{1/3} + b_2^{(j)} \tau^{5/6}
        + b_3^{(j)} \tau^{1.2083} + b_4^{(j)}\tau + b_5^{(j)} \tau^2 + b_6^{(j)} \tau^3

        \tau = 1-T/T_c

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    Tc : float
        Critical temperature of fluid [K]
    omega : float
        Acentric factor [-]

    Returns
    -------
    Hvap : float
        Enthalpy of vaporization, [J/mol]

    Notes
    -----
    The original article has been reviewed. A total of 18 coefficients are used:

    WARNING: The correlation has been implemented as described in the article,
    but its results seem different and with some error.
    Its results match with other functions however.

    Has poor behavior for low-temperature use.

    Examples
    --------
    Problem in article for SMK function.

    >>> MK(553.15, 751.35, 0.302)
    38727.993546377205

    References
    ----------
    .. [1] Morgan, David L., and Riki Kobayashi. "Extension of Pitzer CSP
       Models for Vapor Pressures and Heats of Vaporization to Long-Chain
       Hydrocarbons." Fluid Phase Equilibria 94 (March 15, 1994): 51-87.
       doi:10.1016/0378-3812(94)87051-9.
    '''
    bs = [[5.2804, 0.080022, 7.2543],
          [12.8650, 273.23, -346.45],
          [1.1710, 465.08, -610.48],
          [-13.1160, -638.51, 839.89],
          [0.4858, -145.12, 160.05],
          [-1.0880, 74.049, -50.711]]

    tau = 1. - T/Tc
    H0 = (bs[0][0]*tau**(0.3333) + bs[1][0]*tau**(0.8333) + bs[2][0]*tau**(1.2083) +
    bs[3][0]*tau + bs[4][0]*tau**(2) + bs[5][0]*tau**(3))*R*Tc

    H1 = (bs[0][1]*tau**(0.3333) + bs[1][1]*tau**(0.8333) + bs[2][1]*tau**(1.2083) +
    bs[3][1]*tau + bs[4][1]*tau**(2) + bs[5][1]*tau**(3))*R*Tc

    H2 = (bs[0][2]*tau**(0.3333) + bs[1][2]*tau**(0.8333) + bs[2][2]*tau**(1.2083) +
    bs[3][2]*tau + bs[4][2]*tau**(2) + bs[5][2]*tau**(3))*R*Tc

    return H0 + omega*H1 + omega**2*H2

@njit
def Velasco(T, P, Tc, omega):
    r'''Calculates enthalpy of vaporization at arbitrary temperatures using a
    the work of [1]_; requires a chemical's critical temperature and
    acentric factor.

    The enthalpy of vaporization is given by:

    .. math::
        \Delta_{vap} H = RT_c(7.2729 + 10.4962\omega + 0.6061\omega^2)(1-T_r)^{0.38}

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    Tc : float
        Critical temperature of fluid [K]
    omega : float
        Acentric factor [-]

    Returns
    -------
    Hvap : float
        Enthalpy of vaporization, [J/mol]

    Notes
    -----
    The original article has been reviewed. It is regressed from enthalpy of
    vaporization values at 0.7Tr, from 121 fluids in REFPROP 9.1.
    A value in the article was read to be similar, but slightly too low from
    that calculated here.

    Examples
    --------
    From graph, in [1]_ for perfluoro-n-heptane.

    >>> Velasco(333.2, 476.0, 0.5559)
    33299.41734936356

    References
    ----------
    .. [1] Velasco, S., M. J. Santos, and J. A. White. "Extended Corresponding
       States Expressions for the Changes in Enthalpy, Compressibility Factor
       and Constant-Volume Heat Capacity at Vaporization." The Journal of
       Chemical Thermodynamics 85 (June 2015): 68-76.
       doi:10.1016/j.jct.2015.01.011.
    '''
    return (7.2729 + 10.4962*omega + 0.6061*omega**2)*(1-T/Tc)**0.38*R*Tc


### Enthalpy of Vaporization at Normal Boiling Point.

@njit
def Riedel(Tb, Tc, Pc):
    r'''Calculates enthalpy of vaporization at the boiling point, using the
    Ridel [1]_ CSP method. Required information are critical temperature
    and pressure, and boiling point. Equation taken from [2]_ and [3]_.

    The enthalpy of vaporization is given by:

    .. math::
        \Delta_{vap} H=1.093 T_b R\frac{\ln P_c-1.013}{0.930-T_{br}}

    Parameters
    ----------
    Tb : float
        Boiling temperature of fluid [K]
    Tc : float
        Critical temperature of fluid [K]
    Pc : float
        Critical pressure of fluid [Pa]

    Returns
    -------
    Hvap : float
        Enthalpy of vaporization at the normal boiling point, [J/mol]

    Notes
    -----
    This equation has no example calculation in any source. The source has not
    been verified. It is equation 4-144 in Perry's. Perry's also claims that
    errors seldom surpass 5%.

    [2]_ is the source of example work here, showing a calculation at 0.0%
    error.

    Internal units of pressure are bar.

    Examples
    --------
    Pyridine, 0.0% err vs. exp: 35090 J/mol; from Poling [2]_.

    >>> Riedel(388.4, 620.0, 56.3E5)
    35089.78989646058

    References
    ----------
    .. [1] Riedel, L. "Eine Neue Universelle Dampfdruckformel Untersuchungen
       Uber Eine Erweiterung Des Theorems Der Ubereinstimmenden Zustande. Teil
       I." Chemie Ingenieur Technik 26, no. 2 (February 1, 1954): 83-89.
       doi:10.1002/cite.330260206.
    .. [2] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
    .. [3] Green, Don, and Robert Perry. Perry's Chemical Engineers' Handbook,
       Eighth Edition. McGraw-Hill Professional, 2007.
    '''
    Pc = Pc/1E5  # Pa to bar
    Tbr = Tb/Tc
    return 1.093*Tb*R*(log(Pc) - 1.013)/(0.93 - Tbr)

@njit
def Chen(Tb, Tc, Pc):
    r'''Calculates enthalpy of vaporization using the Chen [1]_ correlation
    and a chemical's critical temperature, pressure and boiling point.

    The enthalpy of vaporization is given by:

    .. math::
        \Delta H_{vb} = RT_b \frac{3.978 T_r - 3.958 + 1.555 \ln P_c}{1.07 - T_r}

    Parameters
    ----------
    Tb : float
        Boiling temperature of the fluid [K]
    Tc : float
        Critical temperature of fluid [K]
    Pc : float
        Critical pressure of fluid [Pa]

    Returns
    -------
    Hvap : float
        Enthalpy of vaporization, [J/mol]

    Notes
    -----
    The formulation presented in the original article is similar, but uses
    units of atm and calorie instead. The form in [2]_ has adjusted for this.
    A method for estimating enthalpy of vaporization at other conditions
    has also been developed, but the article is unclear on its implementation.
    Based on the Pitzer correlation.

    Internal units: bar and K

    Examples
    --------
    Same problem as in Perry's examples.

    >>> Chen(294.0, 466.0, 5.55E6)
    26705.893506174052

    References
    ----------
    .. [1] Chen, N. H. "Generalized Correlation for Latent Heat of Vaporization."
       Journal of Chemical & Engineering Data 10, no. 2 (April 1, 1965): 207-10.
       doi:10.1021/je60025a047
    .. [2] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
    '''
    Tbr = Tb/Tc
    Pc = Pc/1E5  # Pa to bar
    return R*Tb*(3.978*Tbr - 3.958 + 1.555*log(Pc))/(1.07 - Tbr)

@njit
def Liu(Tb, Tc, Pc):
    r'''Calculates enthalpy of vaporization at the normal boiling point using
    the Liu [1]_ correlation, and a chemical's critical temperature, pressure
    and boiling point.

    The enthalpy of vaporization is given by:

    .. math::
        \Delta H_{vap} = RT_b \left[ \frac{T_b}{220}\right]^{0.0627} \frac{
        (1-T_{br})^{0.38} \ln(P_c/P_A)}{1-T_{br} + 0.38 T_{br} \ln T_{br}}

    Parameters
    ----------
    Tb : float
        Boiling temperature of the fluid [K]
    Tc : float
        Critical temperature of fluid [K]
    Pc : float
        Critical pressure of fluid [Pa]

    Returns
    -------
    Hvap : float
        Enthalpy of vaporization, [J/mol]

    Notes
    -----
    This formulation can be adjusted for lower boiling points, due to the use
    of a rationalized pressure relationship. The formulation is taken from
    the original article.

    A correction for alcohols and organic acids based on carbon number,
    which only modifies the boiling point, is available but not implemented.

    No sample calculations are available in the article.

    Internal units: Pa and K

    Examples
    --------
    Same problem as in Perry's examples

    >>> Liu(294.0, 466.0, 5.55E6)
    26378.566319606754

    References
    ----------
    .. [1] LIU, ZHI-YONG. "Estimation of Heat of Vaporization of Pure Liquid at
       Its Normal Boiling Temperature." Chemical Engineering Communications
       184, no. 1 (February 1, 2001): 221-28. doi:10.1080/00986440108912849.
    '''
    Tbr = Tb/Tc
    return R*Tb*(Tb/220.)**0.0627*(1. - Tbr)**0.38*log(Pc/101325.) \
        / (1 - Tbr + 0.38*Tbr*log(Tbr))

@njit
def Vetere(Tb, Tc, Pc, F=1):
    r'''Calculates enthalpy of vaporization at the boiling point, using the
    Vetere [1]_ CSP method. Required information are critical temperature
    and pressure, and boiling point. Equation taken from [2]_.

    The enthalpy of vaporization is given by:

    .. math::
        \frac {\Delta H_{vap}}{RT_b} = \frac{\tau_b^{0.38}
        \left[ \ln P_c - 0.513 + \frac{0.5066}{P_cT_{br}^2}\right]}
        {\tau_b + F(1-\tau_b^{0.38})\ln T_{br}}

    Parameters
    ----------
    Tb : float
        Boiling temperature of fluid [K]
    Tc : float
        Critical temperature of fluid [K]
    Pc : float
        Critical pressure of fluid [Pa]
    F : float, optional
        Constant for a fluid, [-]

    Returns
    -------
    Hvap : float
        Enthalpy of vaporization at the boiling point, [J/mol]

    Notes
    -----
    The equation cannot be found in the original source. It is believed that a
    second article is its source, or that DIPPR staff have altered the formulation.

    Internal units of pressure are bar.

    Examples
    --------
    Example as in [2]_, p2-487; exp: 25.73

    >>> Vetere(294.0, 466.0, 5.55E6)
    26363.430021286465

    References
    ----------
    .. [1] Vetere, Alessandro. "Methods to Predict the Vaporization Enthalpies
       at the Normal Boiling Temperature of Pure Compounds Revisited."
       Fluid Phase Equilibria 106, no. 1-2 (May 1, 1995): 1–10.
       doi:10.1016/0378-3812(94)02627-D.
    .. [2] Green, Don, and Robert Perry. Perry's Chemical Engineers' Handbook,
       Eighth Edition. McGraw-Hill Professional, 2007.
    '''
    Tbr = Tb/Tc
    taub = 1-Tb/Tc
    Pc = Pc/1E5
    term = taub**0.38*(log(Pc)-0.513 + 0.5066/Pc/Tbr**2) / (taub + F*(1-taub**0.38)*log(Tbr))
    return R*Tb*term


### Enthalpy of Vaporization at STP.


### Enthalpy of Vaporization adjusted for T
@njit
def Watson(T, Hvap_ref, T_Ref, Tc, exponent=0.38):
    '''
    Adjusts enthalpy of vaporization of enthalpy for another temperature, for one temperature.
    '''
    Tr = T/Tc
    Trefr = T_Ref/Tc
    H2 = Hvap_ref*((1-Tr)/(1-Trefr))**exponent
    return H2

@njit
def Alibakhshi(T, Tc, C):
    return (4.5*pi*N_A)**(1./3.)*4.2E-7*(Tc-6.) - R/2.*T*log(T) + C*T

@njit
def VDI_PPDS(T, Tc, A, B, C, D, E):
    tau = 1. - T/Tc
    return R*Tc*(A*tau**(1/3.) + B*tau**(2/3.) + C*tau + D*tau**2 + E*tau**6)


def EnthalpyVaporization(CAS, Tb=None, Tc=None, Pc=None, omega=None,
                         similarity_variable=None, Psat=None, Zg=None):
    models = []
    # if has_CoolProp and self.CASRN in coolprop_dict:
    #     methods.append(COOLPROP)
    #     self.CP_f = coolprop_fluids[self.CASRN]
    #     Tmins.append(self.CP_f.Tt); Tmaxs.append(self.CP_f.Tc)
    if CAS in Perrys2_150.index:
        _, Tc, C1, C2, C3, C4, Perrys2_150_Tmin, Perrys2_150_Tmax = _Perrys2_150_values[Perrys2_150.index.get_loc(CAS)].tolist()
        models.append(ModelEQ106(Tc, C1, C2, C3, C4, Tmin=Perrys2_150_Tmin, Tmax=Perrys2_150_Tmax))
    if CAS in VDI_PPDS_4.index:
        _,  MW, Tc, A, B, C, D, E = _VDI_PPDS_4_values[VDI_PPDS_4.index.get_loc(CAS)].tolist()
        models.append(TDependentModel(VDI_PPDS, 0, Tc, args=(Tc, A, B, C, D, E)))
    if all((Tb, Tc, Pc)):
        models.extend([TDependentModel(f, 0, Tc, args=(Tb, Tc, Pc))
                       for f in (Riedel, Chen, Vetere, Liu)])
    if all((Tc, Pc)):
        models.append(TDependentModel(Clapeyron, 0, Tc, args=(Tc, Pc, Zg, Psat)))
    if all((Tc, omega)):
        models.extend([TDependentModel(f, 0, Tc, args=(Tc, omega)) for f in (MK, SMK, Velasco, Pitzer)])
    if CAS in _VDISaturationDict:
        Ts, Hvaps = VDI_tabular_data(CAS, 'Hvap')
        models.append(InterpolatedTDependentModel(Ts, Hvaps, Ts[0], Ts[-1]))
    if CAS in Alibakhshi_Cs.index and Tc:
        C = float(Alibakhshi_Cs.at[CAS, 'C'])
        models.append(TDependentModel(Alibakhshi, 0, Tc, args=(Tc, C)))
    if CAS in CRCHvap_data.index and not np.isnan(CRCHvap_data.at[CAS, 'HvapTb']):
        Tb = float(CRCHvap_data.at[CAS, 'Tb'])
        Hvap = float(CRCHvap_data.at[CAS, 'HvapTb'])
        models.append(TDependentModel(Watson, 0, 10e6, args=(Hvap, Tb, Tc, 0.38)))
    if CAS in CRCHvap_data.index and not np.isnan(CRCHvap_data.at[CAS, 'Hvap298']):
        Hvap = float(CRCHvap_data.at[CAS, 'Hvap298'])
        models.append(TDependentModel(Watson, 0, 10e6, args=(Hvap, 298., Tc, 0.38)))
    if CAS in GharagheiziHvap_data.index:
        Hvap = float(GharagheiziHvap_data.at[CAS, 'Hvap298'])
        models.append(TDependentModel(Watson, 0, 10e6, args=(Hvap, 298., Tc, 0.38)))
    return TDependentModelHandler(models)


### Heat of Fusion


def Hfus(T=298.15, P=101325, MW=None, AvailableMethods=False, Method=None, CASRN=''):  # pragma: no cover
    '''This function handles the calculation of a chemical's enthalpy of fusion.
    Generally this, is used by the chemical class, as all parameters are passed.
    Calling the function directly works okay.

    Enthalpy of fusion is a weak function of pressure, and its effects are
    neglected.

    This API is considered experimental, and is expected to be removed in a
    future release in favor of a more complete object-oriented interface.

    '''
    def list_methods():
        methods = []
        if CASRN in CRCHfus_data.index:
            methods.append('CRC, at melting point')
        methods.append('None')
        return methods
    if AvailableMethods:
        return list_methods()
    if not Method:
        Method = list_methods()[0]
    # This is the calculate, given the method section
    if Method == 'CRC, at melting point':
        _Hfus = CRCHfus_data.at[CASRN, 'Hfus']
    elif Method == 'None' or not MW:
        _Hfus = None
    else:
        raise Exception('Failure in in function')
    return _Hfus
#print Hfus(CASRN='75-07-0')
#['CRC, at melting point', 'None']


### Heat of Sublimation


def Hsub(T=298.15, P=101325, AvailableMethods=False, Method=None, CASRN=''):  # pragma: no cover
    '''This function handles the calculation of a chemical's enthalpy of sublimation.
    Generally this, is used by the chemical class, as all parameters are passed.


    This API is considered experimental, and is expected to be removed in a
    future release in favor of a more complete object-oriented interface.
    '''
    def list_methods():
        methods = []
#        if Hfus(T=T, P=P, MW=MW, CASRN=CASRN) and Hvap(T=T, P=P, MW=MW, CASRN=CASRN):
#            methods.append('Hfus + Hvap')
        if CASRN in GharagheiziHsub_data.index:
            methods.append('Ghazerati Appendix, at 298K')
        methods.append('None')
        return methods
    if AvailableMethods:
        return list_methods()
    if not Method:
        Method = list_methods()[0]
    # This is the calculate, given the method section
#    if Method == 'Hfus + Hvap':
#        p1 = Hfus(T=T, P=P, MW=MW, CASRN=CASRN)
#        p2 = Hvap(T=T, P=P, MW=MW, CASRN=CASRN)
#        if p1 and p2:
#            _Hsub = p1 + p2
#        else:
#            _Hsub = None
    if Method == 'Ghazerati Appendix, at 298K':
        _Hsub = float(GharagheiziHsub_data.at[CASRN, 'Hsub'])
    elif Method == 'None':
        return None
    else:
        raise Exception('Failure in in function')
    return _Hsub


#print Hsub(CASRN='101-81-5')


### Liquidus line for mixtures
def Tliquidus(Tms=None, ws=None, xs=None, CASRNs=None, AvailableMethods=False,
              Method=None):  # pragma: no cover
    '''This function handles the retrival of a mixtures's liquidus point.

    This API is considered experimental, and is expected to be removed in a
    future release in favor of a more complete object-oriented interface.

    >>> Tliquidus(Tms=[250.0, 350.0], xs=[0.5, 0.5])
    350.0
    >>> Tliquidus(Tms=[250, 350], xs=[0.5, 0.5], Method='Simple')
    300.0
    >>> Tliquidus(Tms=[250, 350], xs=[0.5, 0.5], AvailableMethods=True)
    ['Maximum', 'Simple', 'None']
    '''
    def list_methods():
        methods = []
        if none_and_length_check([Tms]):
            methods.append('Maximum')
            methods.append('Simple')
        methods.append('None')
        return methods
    if AvailableMethods:
        return list_methods()
    if not Method:
        Method = list_methods()[0]
    # This is the calculate, given the method section
    if Method == 'Maximum':
        _Tliq = max(Tms)
    elif Method == 'Simple':
        _Tliq = mixing_simple(xs, Tms)
    elif Method == 'None':
        return None
    else:
        raise Exception('Failure in in function')
    return _Tliq

#print Tliquidus(Tms=[250, 350], xs=[0.5, 0.5], AvailableMethods=True)
