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

__all__ = ['isobaric_expansion', 'isothermal_compressibility', 
'Cp_minus_Cv', 'speed_of_sound', 'Joule_Thomson',
'phase_identification_parameter', 'phase_identification_parameter_phase',
'isentropic_exponent', 'Vm_to_rho', 'rho_to_Vm', 
'Z', 'B_to_Z', 'B_from_Z', 'Z_from_virial_density_form', 
'Z_from_virial_pressure_form', 'zs_to_ws', 'ws_to_zs', 'zs_to_Vfs', 
'Vfs_to_zs', 'none_and_length_check', 'normalize', 'mixing_simple', 
'mixing_logarithmic', 'has_matplotlib', 'to_num', 'CAS2int', 
'int2CAS', 'Parachor', 'SG_to_API', 'API_to_SG', 'SG',
'phase_select_property', 'allclose_variable', 'horner', 
'polylog2']

from cmath import sqrt as csqrt
from bisect import bisect_left
import numpy as np
from scipy.optimize import brenth
from scipy.misc import derivative
from scipy.integrate import quad
from scipy.interpolate import interp1d, interp2d
import os
import pandas as pd

try:
    import matplotlib.pyplot as plt
    has_matplotlib = True
except:
    has_matplotlib = False

from math import (acos, acosh, asin, asinh, atan, atan2, atanh, ceil, copysign,
                  cos, cosh, degrees, e,  exp, fabs, 
                  factorial, floor, fmod, frexp, fsum, hypot, isinf, 
                  isnan, ldexp, log, log10, log1p, modf, pi, pow, 
                  radians, sin, sinh, sqrt, tan, tanh, trunc) # Not supported in Python 2.6: expm1, erf, erfc,gamma lgamma

__all__.extend(['acos', 'acosh', 'asin', 'asinh', 'atan', 'atan2', 'atanh', 
'ceil', 'copysign', 'cos', 'cosh', 'degrees', 'e', 'exp', 
'fabs', 'factorial', 'floor', 'fmod', 'frexp', 'fsum', 
'hypot', 'isinf', 'isnan', 'ldexp',  'log', 'log10', 'log1p', 'modf', 
'pi', 'pow', 'radians', 'sin', 'sinh', 'sqrt', 'tan', 'tanh', 'trunc'])
__all__.extend(['R', 'k', 'N_A', 'calorie', 'epsilon_0']) # 'expm1', 'erf', 'erfc',  'lgamma', 'gamma', 
    
# Obtained from SciPy 0.19 (2014 CODATA)
# Included here so calculations are consistent across SciPy versions
R = 8.3144598
k = 1.38064852e-23
N_A = 6.022140857e+23
calorie = 4.184
epsilon_0 = 8.854187817620389e-12

def sample_dict(folder, name, excluded=()):
    dirname = os.path.join(folder, name)
    table = pd.read_csv(dirname, sep='\t', index_col=0)
    columns = [j for i,j in enumerate(table.columns) if i not in excluded]
    data = np.array(table[columns])
    return {j:data[i, :] for (i, j) in enumerate(table.index)}

def to_num(values):
    r'''Legacy function to turn a list of strings into either floats
    (if numeric), stripped strings (if not) or None if the string is empty.
    Accepts any numeric formatting the float function does.

    Parameters
    ----------
    values : list
        list of strings

    Returns
    -------
    values : list
        list of floats, strings, and None values [-]

    Examples
    --------
    >>> to_num(['1', '1.1', '1E5', '0xB4', ''])
    [1.0, 1.1, 100000.0, '0xB4', None]
    '''
    for i in range(len(values)):
        try:
            values[i] = float(values[i])
        except:
            if values[i] == '':
                values[i] = None
            else:
                values[i] = values[i].strip()
                pass
    return values


def CAS2int(i):
    r'''Converts CAS number of a compounds from a string to an int. This is
    helpful when storing large amounts of CAS numbers, as their strings take up
    more memory than their numerical representational. All CAS numbers fit into
    64 bit ints.

    Parameters
    ----------
    CASRN : string
        CASRN [-]

    Returns
    -------
    CASRN : int
        CASRN [-]

    Notes
    -----
    Accomplishes conversion by removing dashes only, and then converting to an
    int. An incorrect CAS number will change without exception.

    Examples
    --------
    >>> CAS2int('7704-34-9')
    7704349
    '''
    return int(i.replace('-', ''))


def int2CAS(i):
    r'''Converts CAS number of a compounds from an int to an string. This is
    helpful when dealing with int CAS numbers.

    Parameters
    ----------
    CASRN : int
        CASRN [-]

    Returns
    -------
    CASRN : string
        CASRN [-]

    Notes
    -----
    Handles CAS numbers with an unspecified number of digits. Does not work on
    floats.

    Examples
    --------
    >>> int2CAS(7704349)
    '7704-34-9'
    '''
    i = str(i)
    return i[:-3]+'-'+i[-3:-1]+'-'+i[-1]


def Parachor(MW, rhol, rhog, sigma):
    r'''Calculate Parachor for a pure species, using its density in the
    liquid and gas phases, surface tension, and molecular weight.

    .. math::
        P = \frac{\sigma^{0.25} MW}{\rho_L - \rho_V}
    
    Parameters
    ----------
    MW : float
        Molecular weight, [g/mol]
    rhol : float
        Liquid density [kg/m^3]
    rhog : float
        Gas density [kg/m^3]
    sigma : float
        Surface tension, [N/m]

    Returns
    -------
    P : float
        Parachor, [N^0.25*m^2.75/mol]

    Notes
    -----
    To convert the output of this function to units of [mN^0.25*m^2.75/kmol], 
    multiply by 5623.4132519.
    
    Values in group contribution tables for Parachor are often listed as 
    dimensionless, in which they are multiplied by 5623413 and the appropriate
    units to make them dimensionless.
    
    Examples
    --------
    Calculating Parachor from a known surface tension for methyl isobutyl 
    ketone at 293.15 K
    
    >>> Parachor(100.15888, 800.8088185536124, 4.97865317223119, 0.02672166960656005)
    5.088443542210164e-05
    
    Converting to the `dimensionless` form:
    
    >>> 5623413*5.088443542210164e-05
    286.14419565030687
    
    Compared to 274.9 according to a group contribution method described in
    [3]_.

    References
    ----------
    .. [1] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
    .. [2] Green, Don, and Robert Perry. Perry's Chemical Engineers' Handbook,
       8E. McGraw-Hill Professional, 2007.
    .. [3] Danner, Ronald P, and Design Institute for Physical Property Data.
       Manual for Predicting Chemical Process Design Data. New York, N.Y, 1982.
    '''
    rhol, rhog = rhol*1000., rhog*1000. # Convert kg/m^3 to g/m^3
    return sigma**0.25*MW/(rhol-rhog) # (N/m)**0.25*g/mol/(g/m^3)

def SG_to_API(SG):
    r'''Calculates specific gravity of a liquid given its API, as shown in
    [1]_.

    .. math::
        \text{API gravity} = \frac{141.5}{\text{SG}} - 131.5

    Parameters
    ----------
    SG : float
        Specific gravity of the fluid at 60 degrees Farenheight [-]

    Returns
    -------
    API : float
        API of the fluid [-]

    Notes
    -----
    Defined only at 60 degrees Fahrenheit.

    Examples
    --------
    >>> SG_to_API(0.7365)
    60.62491513917175

    References
    ----------
    .. [1] API Technical Data Book: General Properties & Characterization.
    American Petroleum Institute, 7E, 2005.
    '''
    return 141.5/SG - 131.5


def API_to_SG(API):
    r'''Calculates API of a liquid given its specific gravity, as shown in
    [1]_.

    .. math::
        \text{SG at}~60^\circ\text{F} = \frac{141.5}{\text{API gravity} +131.5}

    Parameters
    ----------
    API : float
        API of the fluid [-]

    Returns
    -------
    SG : float
        Specific gravity of the fluid at 60 degrees Farenheight [-]

    Notes
    -----
    Defined only at 60 degrees Fahrenheit.

    Examples
    --------
    >>> API_to_SG(60.62)
    0.7365188423901728

    References
    ----------
    .. [1] API Technical Data Book: General Properties & Characterization.
    American Petroleum Institute, 7E, 2005.
    '''
    return 141.5/(API + 131.5)


def SG(rho, rho_ref=999.0170824078306):
    r'''Calculates the specific gravity of a substance with respect to another
    substance; by default, this is water at 15.6 °C (60 °F). For gases, 
    normally the reference density is 1.2 kg/m^3, that of dry air. However, in 
    general specific gravity should always be specified with respect to the
    temperature and pressure of its reference fluid. This can vary widely.
    
    .. math::
        SG = \frac{\rho}{\rho_{ref}}

    Parameters
    ----------
    rho : float
        Density of the substance, [kg/m^3]
    rho_ref : float, optional
        Density of the reference substance, [kg/m^3]

    Returns
    -------
    SG : float
        Specific gravity of the substance with respect to the reference 
        density, [-]

    Notes
    -----
    Another common reference point is water at 4°C (rho_ref=999.9748691393087).
    Specific gravity is often used by consumers instead of density.
    The reference for solids is normally the same as for liquids - water.
    
    Examples
    --------
    >>> SG(860)
    0.8608461408159591
    '''
    return rho/rho_ref


def isobaric_expansion(V, dV_dT):
    r'''Calculate the isobaric coefficient of a thermal expansion, given its 
    molar volume at a certain `T` and `P`, and its derivative of molar volume
    with respect to `T`.

    .. math::
        \beta = \frac{1}{V}\left(\frac{\partial V}{\partial T} \right)_P

    Parameters
    ----------
    V : float
        Molar volume at `T` and `P`, [m^3/mol]
    dV_dT : float
        Derivative of molar volume with respect to `T`, [m^3/mol/K]

    Returns
    -------
    beta : float
        Isobaric coefficient of a thermal expansion, [1/K]
        
    Notes
    -----
    For an ideal gas, this expression simplified to:
    
    .. math::
        \beta = \frac{1}{T}

    Examples
    --------
    Calculated for hexane from the PR EOS at 299 K and 1 MPa (liquid):
    
    >>> isobaric_expansion(0.000130229900873546, 1.58875261849113e-7)
    0.0012199599384121608

    References
    ----------
    .. [1] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
    '''
    return dV_dT/V


def isothermal_compressibility(V, dV_dP):
    r'''Calculate the isothermal coefficient of a compressibility, given its 
    molar volume at a certain `T` and `P`, and its derivative of molar volume
    with respect to `P`.

    .. math::
        \kappa = -\frac{1}{V}\left(\frac{\partial V}{\partial P} \right)_T

    Parameters
    ----------
    V : float
        Molar volume at `T` and `P`, [m^3/mol]
    dV_dP : float
        Derivative of molar volume with respect to `P`, [m^3/mol/Pa]

    Returns
    -------
    kappa : float
        Isothermal coefficient of a compressibility, [1/Pa]
        
    Notes
    -----
    For an ideal gas, this expression simplified to:
    
    .. math::
        \kappa = \frac{1}{P}

    Examples
    --------
    Calculated for hexane from the PR EOS at 299 K and 1 MPa (liquid):
    
    >>> isothermal_compressibility(0.000130229900873546, -2.72902118209903e-13)
    2.095541165119158e-09

    References
    ----------
    .. [1] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
    '''
    return -dV_dP/V


def phase_identification_parameter(V, dP_dT, dP_dV, d2P_dV2, d2P_dVdT):
    r'''Calculate the Phase Identification Parameter developed in [1]_ for
    the accurate and efficient determination of whether a fluid is a liquid or
    a gas based on the results of an equation of state. For supercritical 
    conditions, this provides a good method for choosing which property 
    correlations to use.
    
    .. math::
        \Pi = V \left[\frac{\frac{\partial^2 P}{\partial V \partial T}}
        {\frac{\partial P }{\partial T}}- \frac{\frac{\partial^2 P}{\partial 
        V^2}}{\frac{\partial P}{\partial V}} \right]

    Parameters
    ----------
    V : float
        Molar volume at `T` and `P`, [m^3/mol]
    dP_dT : float
        Derivative of `P` with respect to `T`, [Pa/K]
    dP_dV : float
        Derivative of `P` with respect to `V`, [Pa*mol/m^3]
    d2P_dV2 : float
        Second derivative of `P` with respect to `V`, [Pa*mol^2/m^6]
    d2P_dVdT : float
        Second derivative of `P` with respect to both `V` and `T`, [Pa*mol/m^3/K]

    Returns
    -------
    PIP : float
        Phase Identification Parameter, [-]
        
    Notes
    -----
    Heuristics were used by process simulators before the invent of this 
    parameter. 
    
    The criteria for liquid is Pi > 1; for vapor, Pi <= 1.
    
    There is also a solid phase mechanism available. For solids, the Solid  
    Phase Identification Parameter is greater than 1, like liquids; however,  
    unlike liquids, d2P_dVdT is always >0; it is < 0 for liquids and gases.

    Examples
    --------
    Calculated for hexane from the PR EOS at 299 K and 1 MPa (liquid):
    
    >>> phase_identification_parameter(0.000130229900874, 582169.397484, 
    ... -3.66431747236e+12, 4.48067893805e+17, -20518995218.2)
    11.33428990564796

    References
    ----------
    .. [1] Venkatarathnam, G., and L. R. Oellrich. "Identification of the Phase
       of a Fluid Using Partial Derivatives of Pressure, Volume, and 
       Temperature without Reference to Saturation Properties: Applications in 
       Phase Equilibria Calculations." Fluid Phase Equilibria 301, no. 2 
       (February 25, 2011): 225-33. doi:10.1016/j.fluid.2010.12.001.
    .. [2] Jayanti, Pranava Chaitanya, and G. Venkatarathnam. "Identification
       of the Phase of a Substance from the Derivatives of Pressure, Volume and
       Temperature, without Prior Knowledge of Saturation Properties: Extension
       to Solid Phase." Fluid Phase Equilibria 425 (October 15, 2016): 269-277.
       doi:10.1016/j.fluid.2016.06.001.
    '''
    return V*(d2P_dVdT/dP_dT - d2P_dV2/dP_dV)


def phase_identification_parameter_phase(d2P_dVdT, V=None, dP_dT=None, dP_dV=None, d2P_dV2=None):
    r'''Uses the Phase Identification Parameter concept developed in [1]_ and 
    [2]_ to determine if a chemical is a solid, liquid, or vapor given the 
    appropriate thermodynamic conditions.

    The criteria for liquid is PIP > 1; for vapor, PIP <= 1.

    For solids, PIP(solid) is defined to be d2P_dVdT. If it is larger than 0, 
    the species is a solid. It is less than 0 for all liquids and gases.

    Parameters
    ----------
    d2P_dVdT : float
        Second derivative of `P` with respect to both `V` and `T`, [Pa*mol/m^3/K]
    V : float, optional
        Molar volume at `T` and `P`, [m^3/mol]
    dP_dT : float, optional
        Derivative of `P` with respect to `T`, [Pa/K]
    dP_dV : float, optional
        Derivative of `P` with respect to `V`, [Pa*mol/m^3]
    d2P_dV2 : float, optionsl
        Second derivative of `P` with respect to `V`, [Pa*mol^2/m^6]

    Returns
    -------
    phase : str
        Either 's', 'l' or 'g'
    
    Notes
    -----
    The criteria for being a solid phase is checked first, which only
    requires d2P_dVdT. All other inputs are optional for this reason.
    However, an exception will be raised if the other inputs become 
    needed to determine if a species is a liquid or a gas.
        
    Examples
    --------
    Calculated for hexane from the PR EOS at 299 K and 1 MPa (liquid):
    
    >>> phase_identification_parameter_phase(-20518995218.2, 0.000130229900874, 
    ... 582169.397484, -3.66431747236e+12, 4.48067893805e+17)
    'l'

    References
    ----------
    .. [1] Venkatarathnam, G., and L. R. Oellrich. "Identification of the Phase
       of a Fluid Using Partial Derivatives of Pressure, Volume, and 
       Temperature without Reference to Saturation Properties: Applications in 
       Phase Equilibria Calculations." Fluid Phase Equilibria 301, no. 2 
       (February 25, 2011): 225-33. doi:10.1016/j.fluid.2010.12.001.
    .. [2] Jayanti, Pranava Chaitanya, and G. Venkatarathnam. "Identification
       of the Phase of a Substance from the Derivatives of Pressure, Volume and
       Temperature, without Prior Knowledge of Saturation Properties: Extension
       to Solid Phase." Fluid Phase Equilibria 425 (October 15, 2016): 269-277.
       doi:10.1016/j.fluid.2016.06.001.
    '''
    if d2P_dVdT > 0:
        return 's'
    else:
        PIP = phase_identification_parameter(V=V, dP_dT=dP_dT, dP_dV=dP_dV, 
                                             d2P_dV2=d2P_dV2, d2P_dVdT=d2P_dVdT)
        return 'l' if PIP > 1 else 'g'


def Cp_minus_Cv(T, dP_dT, dP_dV):
    r'''Calculate the difference between a real gas's constant-pressure heat
    capacity and constant-volume heat capacity, as given in [1]_, [2]_, and
    [3]_. The required derivatives should be calculated with an equation of
    state.

    .. math::
        C_p - C_v = -T\left(\frac{\partial P}{\partial T}\right)_V^2/
        \left(\frac{\partial P}{\partial V}\right)_T

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    dP_dT : float
        Derivative of `P` with respect to `T`, [Pa/K]
    dP_dV : float
        Derivative of `P` with respect to `V`, [Pa*mol/m^3]

    Returns
    -------
    Cp_minus_Cv : float
        Cp - Cv for a real gas, [J/mol/K]
        
    Notes
    -----
    Equivalent expressions are:
    
    .. math::
        C_p - C_v= -T\left(\frac{\partial V}{\partial T}\right)_P^2/\left(
        \frac{\partial V}{\partial P}\right)_T
        
        C_p - C_v = T\left(\frac{\partial P}{\partial T}\right)
        \left(\frac{\partial V}{\partial T}\right)

    Note that these are not second derivatives, only first derivatives, some
    of which are squared.

    Examples
    --------
    Calculated for hexane from the PR EOS at 299 K and 1 MPa (liquid):
    
    >>> Cp_minus_Cv(299, 582232.475794113, -3665180614672.253)
    27.654681381642394

    References
    ----------
    .. [1] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
    .. [2] Walas, Stanley M. Phase Equilibria in Chemical Engineering. 
       Butterworth-Heinemann, 1985.
    .. [3] Gmehling, Jurgen, Barbel Kolbe, Michael Kleiber, and Jurgen Rarey.
       Chemical Thermodynamics for Process Simulation. 1st edition. Weinheim: 
       Wiley-VCH, 2012.
    '''
    return -T*dP_dT**2/dP_dV
    
    
def speed_of_sound(V, dP_dV, Cp, Cv, MW=None):
    r'''Calculate a real fluid's speed of sound. The required derivatives should 
    be calculated with an equation of state, and `Cp` and `Cv` are both the
    real fluid versions. Expression is given in [1]_ and [2]_; a unit conversion
    is further performed to obtain a result in m/s. If MW is not provided the 
    result is returned in units of m*kg^0.5/s/mol^0.5.

    .. math::
        w = \left[-V^2 \left(\frac{\partial P}{\partial V}\right)_T \frac{C_p}
        {C_v}\right]^{1/2}
        
    Parameters
    ----------
    V : float
        Molar volume of fluid, [m^3/mol]
    dP_dV : float
        Derivative of `P` with respect to `V`, [Pa*mol/m^3]
    Cp : float
        Real fluid heat capacity at constant pressure, [J/mol/K]
    Cv : float
        Real fluid heat capacity at constant volume, [J/mol/K]
    MW : float, optional
        Molecular weight, [g/mol]

    Returns
    -------
    w : float
        Speed of sound for a real gas, [m/s or m*kg^0.5/s/mol^0.5 or MW missing]
        
    Notes
    -----
    An alternate expression based on molar density is as follows:
    
    .. math::
       w = \left[\left(\frac{\partial P}{\partial \rho}\right)_T \frac{C_p}
       {C_v}\right]^{1/2}

    The form with the unit conversion performed inside it is as follows:
    
    .. math::
        w = \left[-V^2 \frac{1000}{MW}\left(\frac{\partial P}{\partial V}
        \right)_T \frac{C_p}{C_v}\right]^{1/2}
    
    Examples
    --------
    Example from [2]_:
    
    >>> speed_of_sound(V=0.00229754, dP_dV=-3.5459e+08, Cp=153.235, Cv=132.435, MW=67.152)
    179.5868138460819

    References
    ----------
    .. [1] Gmehling, Jurgen, Barbel Kolbe, Michael Kleiber, and Jurgen Rarey.
       Chemical Thermodynamics for Process Simulation. 1st edition. Weinheim: 
       Wiley-VCH, 2012.
    .. [2] Pratt, R. M. "Thermodynamic Properties Involving Derivatives: Using 
       the Peng-Robinson Equation of State." Chemical Engineering Education 35,
       no. 2 (March 1, 2001): 112-115. 
    '''
    if not MW:
        return (-V**2*dP_dV*Cp/Cv)**0.5
    else:
        return (-V**2*1000./MW*dP_dV*Cp/Cv)**0.5


def Joule_Thomson(T, V, Cp, dV_dT=None, beta=None):
    r'''Calculate a real fluid's Joule Thomson coefficient. The required 
    derivative should be calculated with an equation of state, and `Cp` is the
    real fluid versions. This can either be calculated with `dV_dT` directly, 
    or with `beta` if it is already known.

    .. math::
        \mu_{JT} = \left(\frac{\partial T}{\partial P}\right)_H = \frac{1}{C_p}
        \left[T \left(\frac{\partial V}{\partial T}\right)_P - V\right]
        = \frac{V}{C_p}\left(\beta T-1\right)
        
    Parameters
    ----------
    T : float
        Temperature of fluid, [K]
    V : float
        Molar volume of fluid, [m^3/mol]
    Cp : float
        Real fluid heat capacity at constant pressure, [J/mol/K]
    dV_dT : float, optional
        Derivative of `V` with respect to `T`, [m^3/mol/K]
    beta : float, optional
        Isobaric coefficient of a thermal expansion, [1/K]

    Returns
    -------
    mu_JT : float
        Joule-Thomson coefficient [K/Pa]
            
    Examples
    --------
    Example from [2]_:
    
    >>> Joule_Thomson(T=390, V=0.00229754, Cp=153.235, dV_dT=1.226396e-05)
    1.621956080529905e-05

    References
    ----------
    .. [1] Walas, Stanley M. Phase Equilibria in Chemical Engineering. 
       Butterworth-Heinemann, 1985.
    .. [2] Pratt, R. M. "Thermodynamic Properties Involving Derivatives: Using 
       the Peng-Robinson Equation of State." Chemical Engineering Education 35,
       no. 2 (March 1, 2001): 112-115. 
    '''
    if dV_dT:
        return (T*dV_dT - V)/Cp
    elif beta:
        return V/Cp*(beta*T - 1.)
    else:
        raise Exception('Either dV_dT or beta is needed')


def isentropic_exponent(Cp, Cv):
    r'''Calculate the isentropic coefficient of a gas, given its constant-
    pressure and constant-volume heat capacity.

    .. math::
        k = \frac{C_p}{C_v}

    Parameters
    ----------
    Cp : float
        Gas heat capacity at constant pressure, [J/mol/K]
    Cv : float
        Gas heat capacity at constant volume, [J/mol/K]

    Returns
    -------
    k : float
        Isentropic exponent, [-]

    Examples
    --------
    >>> isentropic_exponent(33.6, 25.27)
    1.329639889196676

    References
    ----------
    .. [1] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
    '''
    return Cp/Cv


def Vm_to_rho(Vm, MW):
    r'''Calculate the density of a chemical, given its molar volume and
    molecular weight.

    .. math::
        \rho = \frac{MW}{1000\cdot VM}

    Parameters
    ----------
    Vm : float
        Molar volume, [m^3/mol]
    MW : float
        Molecular weight, [g/mol]

    Returns
    -------
    rho : float
        Density, [kg/m^3]

    Examples
    --------
    >>> Vm_to_rho(0.000132, 86.18)
    652.8787878787879

    References
    ----------
    .. [1] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
    '''
    return (Vm)**-1*MW/1000.


def rho_to_Vm(rho, MW):
    r'''Calculate the molar volume of a chemical, given its density and
    molecular weight.

    .. math::
        V_m = \left(\frac{1000 \rho}{MW}\right)^{-1}

    Parameters
    ----------
    rho : float
        Density, [kg/m^3]
    MW : float
        Molecular weight, [g/mol]

    Returns
    -------
    Vm : float
        Molar volume, [m^3/mol]

    Examples
    --------
    >>> rho_to_Vm(652.9, 86.18)
    0.00013199571144126206

    References
    ----------
    .. [1] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
    '''
    return (rho*1000./MW)**-1


def Z(T, P, V):
    r'''Calculates the compressibility factor of a gas, given its
    temperature, pressure, and molar volume.

    .. math::
        Z = \frac{PV}{RT}

    Parameters
    ----------
    T : float
        Temperature, [K]
    P : float
        Pressure [Pa]
    V : float
        Molar volume, [m^3/mol]

    Returns
    -------
    Z : float
        Compressibility factor, [-]

    Examples
    --------
    >>> Z(600, P=1E6, V=0.00463)
    0.9281019876560912

    References
    ----------
    .. [1] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
    '''
    return V*P/T/R


def B_to_Z(B, T, P):
    r'''Calculates the compressibility factor of a gas, given its
    second virial coefficient.

    .. math::
        Z = 1 + \frac{BP}{RT}

    Parameters
    ----------
    B : float
        Second virial coefficient, [m^3/mol]
    T : float
        Temperature, [K]
    P : float
        Pressure [Pa]

    Returns
    -------
    Z : float
        Compressibility factor, [-]

    Notes
    -----
    Other forms of the virial coefficient exist.

    Examples
    --------
    >>> B_to_Z(-0.0015, 300, 1E5)
    0.9398638020957176

    References
    ----------
    .. [1] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
    '''
    return 1. + B*P/R/T


def B_from_Z(Z, T, P):
    r'''Calculates the second virial coefficient of a pure species, given the
    compressibility factor of the gas.

    .. math::
        B = \frac{RT(Z-1)}{P}

    Parameters
    ----------
    Z : float
        Compressibility factor, [-]
    T : float
        Temperature, [K]
    P : float
        Pressure [Pa]

    Returns
    -------
    B : float
        Second virial coefficient, [m^3/mol]

    Notes
    -----
    Other forms of the virial coefficient exist.

    Examples
    --------
    >>> B_from_Z(0.94, 300, 1E5)
    -0.0014966027640000014

    References
    ----------
    .. [1] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
    '''
    return (Z - 1)*R*T/P


def Z_from_virial_density_form(T, P, *args):
    r'''Calculates the compressibility factor of a gas given its temperature,
    pressure, and molar density-form virial coefficients. Any number of
    coefficients is supported.

    .. math::
        Z = \frac{PV}{RT} = 1 + \frac{B}{V} + \frac{C}{V^2} + \frac{D}{V^3}
        + \frac{E}{V^4} \dots

    Parameters
    ----------
    T : float
        Temperature, [K]
    P : float
        Pressure, [Pa]
    B to Z : float, optional
        Virial coefficients, [various]

    Returns
    -------
    Z : float
        Compressibility factor at T, P, and with given virial coefficients, [-]

    Notes
    -----
    For use with B or with B and C or with B and C and D, optimized equations 
    are used to obtain the compressibility factor directly.
    If more coefficients are provided, uses numpy's roots function to solve 
    this equation. This takes substantially longer as the solution is 
    numerical.
    
    If no virial coefficients are given, returns 1, as per the ideal gas law.
    
    The units of each virial coefficient are as follows, where for B, n=1, and
    C, n=2, and so on.
    
    .. math::
        \left(\frac{\text{m}^3}{\text{mol}}\right)^n

    Examples
    --------
    >>> Z_from_virial_density_form(300, 122057.233762653, 1E-4, 1E-5, 1E-6, 1E-7)
    1.2843496002100001

    References
    ----------
    .. [1] Prausnitz, John M., Rudiger N. Lichtenthaler, and Edmundo Gomes de 
       Azevedo. Molecular Thermodynamics of Fluid-Phase Equilibria. 3rd 
       edition. Upper Saddle River, N.J: Prentice Hall, 1998.
    .. [2] Walas, Stanley M. Phase Equilibria in Chemical Engineering. 
       Butterworth-Heinemann, 1985.
    '''
    l = len(args)
    if l == 1:
        return 1/2. + (4*args[0]*P + R*T)**0.5/(2*(R*T)**0.5)
#        return ((R*T*(4*args[0]*P + R*T))**0.5 + R*T)/(2*P)
    if l == 2:
        B, C = args
        # A small imaginary part is ignored
        return (P*(-(3*B*R*T/P + R**2*T**2/P**2)/(3*(-1/2 + csqrt(3)*1j/2)*(-9*B*R**2*T**2/(2*P**2) - 27*C*R*T/(2*P) + csqrt(-4*(3*B*R*T/P + R**2*T**2/P**2)**(3+0j) + (-9*B*R**2*T**2/P**2 - 27*C*R*T/P - 2*R**3*T**3/P**3)**(2+0j))/2 - R**3*T**3/P**3)**(1/3.+0j)) - (-1/2 + csqrt(3)*1j/2)*(-9*B*R**2*T**2/(2*P**2) - 27*C*R*T/(2*P) + csqrt(-4*(3*B*R*T/P + R**2*T**2/P**2)**(3+0j) + (-9*B*R**2*T**2/P**2 - 27*C*R*T/P - 2*R**3*T**3/P**3)**(2+0j))/2 - R**3*T**3/P**3)**(1/3.+0j)/3 + R*T/(3*P))/(R*T)).real
    if l == 3:
        # Huge mess. Ideally sympy could optimize a function for quick python 
        # execution. Derived with kate's text highlighting
        B, C, D = args
        P2 = P**2 
        RT = R*T
        BRT = B*RT
        T2 = T**2
        R2 = R**2
        RT23 = 3*R2*T2
        mCRT = -C*RT
        P2256 = 256*P2
        
        RT23P2256 = RT23/(P2256)
        big1 = (D*RT/P - (-BRT/P - RT23/(8*P2))**2/12 - RT*(mCRT/(4*P) - RT*(BRT/(16*P) + RT23P2256)/P)/P)
        big3 = (-BRT/P - RT23/(8*P2))
        big4 = (mCRT/P - RT*(BRT/(2*P) + R2*T2/(8*P2))/P)
        big5 = big3*(-D*RT/P + RT*(mCRT/(4*P) - RT*(BRT/(16*P) + RT23P2256)/P)/P)
        big2 = 2*big1/(3*(big3**3/216 - big5/6 + big4**2/16 + csqrt(big1**3/27 + (-big3**3/108 + big5/3 - big4**2/8)**2/4))**(1/3))
        big7 = 2*BRT/(3*P) - big2 + 2*(big3**3/216 - big5/6 + big4**2/16 + csqrt(big1**3/27 + (-big3**3/108 + big5/3 - big4**2/8)**2/4))**(1/3) + R2*T2/(4*P2)
        return (P*(((csqrt(big7)/2 + csqrt(4*BRT/(3*P) - (-2*C*RT/P - 2*RT*(BRT/(2*P) + R2*T2/(8*P2))/P)/csqrt(big7) + big2 - 2*(big3**3/216 - big5/6 + big4**2/16 + csqrt(big1**3/27 + (-big3**3/108 + big5/3 - big4**2/8)**2/4))**(1/3) + R2*T2/(2*P2))/2 + RT/(4*P))))/R/T).real

    args = list(args)
    args.reverse()
    args.extend([1, -P/R/T])
    solns = np.roots(args)
    rho = [i for i in solns if not i.imag and i.real > 0][0].real # Quicker than indexing where imag ==0
    return P/rho/R/T


def Z_from_virial_pressure_form(P, *args):
    r'''Calculates the compressibility factor of a gas given its pressure, and 
    pressure-form virial coefficients. Any number of coefficients is supported.

    .. math::
        Z = \frac{Pv}{RT} = 1 + B'P + C'P^2 + D'P^3 + E'P^4 \dots

    Parameters
    ----------
    P : float
        Pressure, [Pa]
    B to Z : float, optional
        Pressure form Virial coefficients, [various]

    Returns
    -------
    Z : float
        Compressibility factor at P, and with given virial coefficients, [-]

    Notes
    -----
    Note that although this function does not require a temperature input, it  
    is still dependent on it because the coefficients themselves normally are
    regressed in terms of temperature.
    
    The use of this form is less common than the density form. Its coefficients
    are normally indicated with the "'" suffix.
    
    If no virial coefficients are given, returns 1, as per the ideal gas law.
    
    The units of each virial coefficient are as follows, where for B, n=1, and
    C, n=2, and so on.
    
    .. math::
        \left(\frac{1}{\text{Pa}}\right)^n

    Examples
    --------
    >>> Z_from_virial_pressure_form(102919.99946855308, 4.032286555169439e-09, 1.6197059494442215e-13, 6.483855042486911e-19)
    1.00283753944
    
    References
    ----------
    .. [1] Prausnitz, John M., Rudiger N. Lichtenthaler, and Edmundo Gomes de 
       Azevedo. Molecular Thermodynamics of Fluid-Phase Equilibria. 3rd 
       edition. Upper Saddle River, N.J: Prentice Hall, 1998.
    .. [2] Walas, Stanley M. Phase Equilibria in Chemical Engineering. 
       Butterworth-Heinemann, 1985.
    '''
    return 1 + P*sum([coeff*P**i for i, coeff in enumerate(args)])


def zs_to_ws(zs, MWs):
    r'''Converts a list of mole fractions to mass fractions. Requires molecular
    weights for all species.

    .. math::
        w_i = \frac{z_i MW_i}{MW_{avg}}

        MW_{avg} = \sum_i z_i MW_i

    Parameters
    ----------
    zs : iterable
        Mole fractions [-]
    MWs : iterable
        Molecular weights [g/mol]

    Returns
    -------
    ws : iterable
        Mass fractions [-]

    Notes
    -----
    Does not check that the sums add to one. Does not check that inputs are of
    the same length.

    Examples
    --------
    >>> zs_to_ws([0.5, 0.5], [10, 20])
    [0.3333333333333333, 0.6666666666666666]
    '''
    Mavg = sum(zi*MWi for zi, MWi in zip(zs, MWs))
    ws = [zi*MWi/Mavg for zi, MWi in zip(zs, MWs)]
    return ws


def ws_to_zs(ws, MWs):
    r'''Converts a list of mass fractions to mole fractions. Requires molecular
    weights for all species.

    .. math::
        z_i = \frac{\frac{w_i}{MW_i}}{\sum_i \frac{w_i}{MW_i}}

    Parameters
    ----------
    ws : iterable
        Mass fractions [-]
    MWs : iterable
        Molecular weights [g/mol]

    Returns
    -------
    zs : iterable
        Mole fractions [-]

    Notes
    -----
    Does not check that the sums add to one. Does not check that inputs are of
    the same length.

    Examples
    --------
    >>> ws_to_zs([0.3333333333333333, 0.6666666666666666], [10, 20])
    [0.5, 0.5]
    '''
    tot = sum(w/MW for w, MW in zip(ws, MWs))
    zs = [w/MW/tot for w, MW in zip(ws, MWs)]
    return zs


def zs_to_Vfs(zs, Vms):
    r'''Converts a list of mole fractions to volume fractions. Requires molar
    volumes for all species.

    .. math::
        \text{Vf}_i = \frac{z_i V_{m,i}}{\sum_i z_i V_{m,i}}

    Parameters
    ----------
    zs : iterable
        Mole fractions [-]
    VMs : iterable
        Molar volumes of species [m^3/mol]

    Returns
    -------
    Vfs : list
        Molar volume fractions [-]

    Notes
    -----
    Does not check that the sums add to one. Does not check that inputs are of
    the same length.

    Molar volumes are specified in terms of pure components only. Function
    works with any phase.

    Examples
    --------
    Acetone and benzene example

    >>> zs_to_Vfs([0.637, 0.363], [8.0234e-05, 9.543e-05])
    [0.5960229712956298, 0.4039770287043703]
    '''
    vol_is = [zi*Vmi for zi, Vmi in zip(zs, Vms)]
    tot = sum(vol_is)
    return [vol_i/tot for vol_i in vol_is]


def Vfs_to_zs(Vfs, Vms):
    r'''Converts a list of mass fractions to mole fractions. Requires molecular
    weights for all species.

    .. math::
        z_i = \frac{\frac{\text{Vf}_i}{V_{m,i}}}{\sum_i
        \frac{\text{Vf}_i}{V_{m,i}}}

    Parameters
    ----------
    Vfs : iterable
        Molar volume fractions [-]
    VMs : iterable
        Molar volumes of species [m^3/mol]

    Returns
    -------
    zs : list
        Mole fractions [-]

    Notes
    -----
    Does not check that the sums add to one. Does not check that inputs are of
    the same length.

    Molar volumes are specified in terms of pure components only. Function
    works with any phase.

    Examples
    --------
    Acetone and benzene example

    >>> Vfs_to_zs([0.596, 0.404], [8.0234e-05, 9.543e-05])
    [0.6369779395901142, 0.3630220604098858]
    '''
    mols_i = [Vfi/Vmi for Vfi, Vmi in zip(Vfs, Vms)]
    mols = sum(mols_i)
    return [mol_i/mols for mol_i in mols_i]


def none_and_length_check(all_inputs, length=None):
    r'''Checks inputs for suitability of use by a mixing rule which requires
    all inputs to be of the same length and non-None. A number of variations
    were attempted for this function; this was found to be the quickest.

    Parameters
    ----------
    all_inputs : array-like of array-like
        list of all the lists of inputs, [-]
    length : int, optional
        Length of the desired inputs, [-]

    Returns
    -------
    False/True : bool
        Returns True only if all inputs are the same length (or length `length`)
        and none of the inputs contain None [-]

    Notes
    -----
    Does not check for nan values.

    Examples
    --------
    >>> none_and_length_check(([1, 1], [1, 1], [1, 30], [10,0]), length=2)
    True
    '''
    if not length:
        length = len(all_inputs[0])
    for things in all_inputs:
        if None in things or len(things) != length:
            return False
    return True


def allclose_variable(a, b, limits, rtols=None, atols=None):
    '''Returns True if two arrays are element-wise equal within several 
    different tolerances. Tolerance values are always positive, usually
    very small. Based on numpy's allclose function.
    
    Only atols or rtols needs to be specified; both are used if given.
    
    Parameters
    ----------
    a, b : array_like
        Input arrays to compare.
    limits : array_like
        Fractions of elements allowed to not match to within each tolerance.
    rtols : array_like
        The relative tolerance parameters.
    atols : float
        The absolute tolerance parameters.

    Returns
    -------
    allclose : bool
        Returns True if the two arrays are equal within the given
        tolerances; False otherwise.
            
    Examples
    --------
    10 random similar variables, all of them matching to within 1E-5, allowing 
    up to half to match up to 1E-6.
    
    >>> x = [2.7244322249597719e-08, 3.0105683900110473e-10, 2.7244124924802327e-08, 3.0105259397637556e-10, 2.7243929226310193e-08, 3.0104990272770901e-10, 2.7243666849384451e-08, 3.0104101821236015e-10, 2.7243433745917367e-08, 3.0103707421519949e-10]
    >>> y = [2.7244328304561904e-08, 3.0105753470546008e-10, 2.724412872417824e-08,  3.0105303055834564e-10, 2.7243914341030203e-08, 3.0104819238021998e-10, 2.7243684057561379e-08, 3.0104299541023674e-10, 2.7243436694839306e-08, 3.010374130526363e-10]
    >>> allclose_variable(x, y, limits=[.0, .5], rtols=[1E-5, 1E-6])
    True
    '''
    l = float(len(a))
    if rtols is None and atols is None:
        raise Exception('Either absolute errors or relative errors must be supplied.')
    elif rtols is None:
        rtols = [0 for i in atols]
    elif atols is None:
        atols = [0 for i in rtols]
    
    for atol, rtol, lim in zip(atols, rtols, limits):
        matches = np.count_nonzero(np.isclose(a, b, rtol=rtol, atol=atol))
        if 1-matches/l > lim:
            return False
    return True


def horner(coeffs, x):
    r'''Simple function to calculate the value of a polynomial at a specific
    value of `x`, using the Horner evaluation scheme

    Parameters
    ----------
    coeffs : array-like
        Coefficients, where coeffs[0] is multiplied by the largest power of x,
        and coeffs[-1] is added to the sum with no multiplication.
    x : float
        Value to evaluate the polynomial at

    Returns
    -------
    y : float
        Evaluated result

    Notes
    -----
    Efficient. Faster than numpy.polyval.

    Examples
    --------
    >>> horner([1,2,3], 3)
    18
    '''
    tot = 0
    for c in coeffs:
        tot = tot * x + c
    return tot


def polylog2(x):
    r'''Simple function to calculate PolyLog(2, x) from ranges 0 <= x <= 1,
    with relative error guaranteed to be < 1E-7 from 0 to 0.99999. This
    is a Pade approximation, with three coefficient sets with splits at 0.7 
    and 0.99. An exception is raised if x is under 0 or above 1. 
    

    Parameters
    ----------
    x : float
        Value to evaluate PolyLog(2, x) T

    Returns
    -------
    y : float
        Evaluated result

    Notes
    -----
    Efficient (2-4 microseconds). No implementation of this function exists in 
    SciPy. Derived with mpmath's pade approximation.
    Required for the entropy integral of 
    :obj:`thermo.heat_capacity.Zabransky_quasi_polynomial`.

    Examples
    --------
    >>> polylog2(0.5)
    0.5822405264516294
    '''
    if 0 <= x <= 0.7:
        p = [0.06184590404457956, -0.7460693871557973, 2.2435704485433376, -2.1944070385048526, 0.3382265629285811, 0.2791966558569478]
        q = [-0.005308735283483908, 0.1823421262956287, -1.2364596896290079, 2.9897802200092296, -2.9365321202088004, 1.0]
        offset = 0.26
    elif 0.7 < x <= 0.99:
        p = [7543860.817140365, -10254250.429758755, -4186383.973408412, 7724476.972409749, -3130743.609030545, 600806.068543299, -62981.15051292659, 3696.7937385473397, -114.06795167646395, 1.4406337969700391]
        q = [-1262997.3422452002, 10684514.56076485, -16931658.916668657, 10275996.02842749, -3079141.9506451315, 511164.4690136096, -49254.56172495263, 2738.0399260270983, -81.36790509581284, 1.0]
        offset = 0.95
    elif 0.99 < x <= 1:
        p = [8.548256176424551e+34, 1.8485781239087334e+35, -2.1706889553798647e+34, 8.318563643438321e+32, -1.559802348661511e+31, 1.698939241177209e+29, -1.180285031647229e+27, 5.531049937687143e+24, -1.8085903366375877e+22, 4.203276811951035e+19, -6.98211620300421e+16, 82281997048841.92, -67157299796.61345, 36084814.54808544, -11478.108105137717, 1.6370226052761176]
        q = [-1.9763570499484274e+35, 1.4813997374958851e+35, -1.4773854824041134e+34, 5.38853721252814e+32, -9.882387315028929e+30, 1.0635231532999732e+29, -7.334629044071992e+26, 3.420655574477631e+24, -1.1147787784365177e+22, 2.584530363912858e+19, -4.285376337404043e+16, 50430830490687.56, -41115254924.43107, 22072284.971253656, -7015.799744041691, 1.0]
        offset = 0.999
    else:
        raise Exception('Approximation is valid between 0 and 1 only.')
    x = x - offset
    return horner(p, x)/horner(q, x)


def normalize(values):
    r'''Simple function which normalizes a series of values to be from 0 to 1,
    and for their sum to add to 1.

    .. math::
        x = \frac{x}{sum_i x_i}

    Parameters
    ----------
    values : array-like
        array of values

    Returns
    -------
    fractions : array-like
        Array of values from 0 to 1

    Notes
    -----
    Does not work on negative values.

    Examples
    --------
    >>> normalize([3, 2, 1])
    [0.5, 0.3333333333333333, 0.16666666666666666]
    '''
    tot = sum(values)
    return [i/tot for i in values]


def mixing_simple(fracs, props):
    r'''Simple function calculates a property based on weighted averages of
    properties. Weights could be mole fractions, volume fractions, mass
    fractions, or anything else.

    .. math::
        y = \sum_i \text{frac}_i \cdot \text{prop}_i

    Parameters
    ----------
    fracs : array-like
        Fractions of a mixture
    props: array-like
        Properties

    Returns
    -------
    prop : value
        Calculated property

    Notes
    -----
    Returns None if any fractions or properties are missing or are not of the
    same length.

    Examples
    --------
    >>> mixing_simple([0.1, 0.9], [0.01, 0.02])
    0.019000000000000003
    '''
    if not none_and_length_check([fracs, props]):
        return None
    result = sum(frac*prop for frac, prop in zip(fracs, props))
    return result


def mixing_logarithmic(fracs, props):
    r'''Simple function calculates a property based on weighted averages of
    logarithmic properties.

    .. math::
        y = \sum_i \text{frac}_i \cdot \log(\text{prop}_i)

    Parameters
    ----------
    fracs : array-like
        Fractions of a mixture
    props: array-like
        Properties

    Returns
    -------
    prop : value
        Calculated property

    Notes
    -----
    Does not work on negative values.
    Returns None if any fractions or properties are missing or are not of the
    same length.

    Examples
    --------
    >>> mixing_logarithmic([0.1, 0.9], [0.01, 0.02])
    0.01866065983073615
    '''
    if not none_and_length_check([fracs, props]):
        return None
    return exp(sum(frac*log(prop) for frac, prop in zip(fracs, props)))


def phase_select_property(phase=None, s=None, l=None, g=None, V_over_F=None):
    r'''Determines which phase's property should be set as a default, given
    the phase a chemical is, and the property values of various phases. For the
    case of liquid-gas phase, returns None. If the property is not available
    for the current phase, or if the current phase is not known, returns None.

    Parameters
    ----------
    phase : str
        One of {'s', 'l', 'g', 'two-phase'}
    s : float
        Solid-phase property
    l : float
        Liquid-phase property
    g : float
        Gas-phase property
    V_over_F : float
        Vapor phase fraction

    Returns
    -------
    prop : float
        The selected/calculated property for the relevant phase

    Notes
    -----
    Could calculate mole-fraction weighted properties for the two phase regime.
    Could also implement equilibria with solid phases.

    Examples
    --------
    >>> phase_select_property(phase='g', l=1560.14, g=3312.)
    3312.0
    '''
    if phase == 's':
        return s
    elif phase == 'l':
        return l
    elif phase == 'g':
        return g
    elif phase == 'two-phase':
        return None  #TODO: all two-phase properties?
    elif phase is None:
        return None
    else:
        raise Exception('Property not recognized')


class MultiCheb1D(object):
    '''Simple class to store set of coefficients for multiple chebyshev 
    approximations and perform calculations from them.
    '''
    def __init__(self, points, coeffs):
        self.points = points
        self.coeffs = coeffs
        self.N = len(points)-1
        
    def __call__(self, x):
        coeffs = self.coeffs[bisect_left(self.points, x)]
        return coeffs(x)
#        return self.chebval(x, coeffs)
                
    @staticmethod
    def chebval(x, c):
        # copied from numpy's source, slightly optimized
        # https://github.com/numpy/numpy/blob/v1.13.0/numpy/polynomial/chebyshev.py#L1093-L1177
        x2 = 2.*x
        c0 = c[-2]
        c1 = c[-1]
        for i in range(3, len(c) + 1):
            tmp = c0
            c0 = c[-i] - c1
            c1 = tmp + c1*x2
        return c0 + c1*x
