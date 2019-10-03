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

# __all__ = ['Poling_data', 'TRC_gas_data', '_PerryI', 'CRC_standard_data', 
#            'LastovkaShawModel', 'Lastovka_Shaw_integral', 
#            'Lastovka_Shaw_integral_over_T', 'TRCCp', 
#            'TRCCp_integral', 'TRCCp_integral_over_T', 
#            'heat_capacity_gas_methods', 'HeatCapacityGas', 
#            'Rowlinson_Poling', 'Rowlinson_Bondi', 'Dadgostar_Shaw', 
#            'Zabransky_quasi_polynomial', 'Zabransky_quasi_polynomial_integral',
#            'Zabransky_quasi_polynomial_integral_over_T', 'Zabransky_cubic', 
#            'Zabransky_cubic_integral', 'Zabransky_cubic_integral_over_T',
#            'Zabransky_quasipolynomial', 'Zabransky_spline',
#            'ZABRANSKY_TO_DICT', 'heat_capacity_liquid_methods', 
#            'HeatCapacityLiquid', 'Lastovka_solid', 'Lastovka_solid_integral', 
#            'Lastovka_solid_integral_over_T', 'heat_capacity_solid_methods', 
#            'HeatCapacitySolid']
import os
from io import open
from .utils import log, exp, polylog2
from cmath import log as clog, exp as cexp
import numpy as np
# import pandas as pd

# from scipy.integrate import quad
from .thermo_model import TDependentModel, ConstantTDependentModel, InterpolatedTDependentModel
from .model_handler import TDependentModelHandler
from .utils import R, calorie
from .utils import to_num
from .miscdata import _VDISaturationDict, VDI_tabular_data
# from .electrochem import (LaliberteHeatCapacityModel,
#                           _Laliberte_Heat_Capacity_ParametersDict)
from .utils import sample_dict
# from .coolprop import has_CoolProp, coolprop_dict, CoolProp_T_dependent_property,\
#                       coolprop_fluids, PropsSI

folder = os.path.join(os.path.dirname(__file__), 'Heat Capacity')
Poling_data = sample_dict(folder, 'PolingDatabank.tsv', (0,))
TRC_gas_data = sample_dict(folder, 'TRC Thermodynamics of Organic Compounds in the Gas State.tsv', (0,))


_PerryI = {}
with open(os.path.join(folder, 'Perrys Table 2-151.tsv'), encoding='utf-8') as f:
    '''Read in a dict of heat capacities of irnorganic and elemental solids.
    These are in section 2, table 151 in:
    Green, Don, and Robert Perry. Perry's Chemical Engineers' Handbook,
    Eighth Edition. McGraw-Hill Professional, 2007.

    Formula:
    Cp(Cal/mol/K) = Const + Lin*T + Quadinv/T^2 + Quadinv*T^2

    Phases: c, gls, l, g.
    '''
    next(f)
    for line in f:
        values = to_num(line.strip('\n').split('\t'))
        (CASRN, _formula, _phase, _subphase, Const, Lin, Quadinv, Quad, Tmin,
         Tmax, err) = values
        if Lin is None:
            Lin = 0
        if Quadinv is None:
            Quadinv = 0
        if Quad is None:
            Quad = 0
        if CASRN in _PerryI and CASRN:
            a = _PerryI[CASRN]
            a.update({_phase: {"Formula": _formula, "Phase": _phase,
                               "Subphase": _subphase, "Const": Const,
                               "Lin": Lin, "Quadinv": Quadinv, "Quad": Quad,
                               "Tmin": Tmin, "Tmax": Tmax, "Error": err}})
            _PerryI[CASRN] = a
        else:
            _PerryI[CASRN] = {_phase: {"Formula": _formula, "Phase": _phase,
                                       "Subphase": _subphase, "Const": Const,
                                       "Lin": Lin, "Quadinv": Quadinv,
                                       "Quad": Quad, "Tmin": Tmin,
                                       "Tmax": Tmax, "Error": err}}


#    '''Read in a dict of 2481 thermodynamic property sets of different phases from:
#        Haynes, W.M., Thomas J. Bruno, and David R. Lide. CRC Handbook of
#        Chemistry and Physics. [Boca Raton, FL]: CRC press, 2014.
#        Warning: 11 duplicated chemicals are present and currently clobbered.
CRC_standard_data = sample_dict(folder, 'CRC Standard Thermodynamic Properties of Chemical Substances.tsv')



### Heat capacities of gases

def LastovkaShawModel(MW, similarity_variable, cyclic_aliphatic=False, *, Tmin, Tmax):
    r'''Return a TDependentModel object that calculates ideal-gas constant-pressure heat capacity (kJ/kmol/K) at arbitrary temperatures (K) with the similarity variable concept and method as shown in [1]_.

    .. math::
        C_p = MW \cdot C_p^{mass}
        C_p^{mass} = \left(A_2 + \frac{A_1 - A_2}{1 + \exp(\frac{\alpha-A_3}{A_4})}\right)
        + (B_{11} + B_{12}\alpha)\left(-\frac{(C_{11} + C_{12}\alpha)}{T}\right)^2
        \frac{\exp(-(C_{11} + C_{12}\alpha)/T)}{[1-\exp(-(C_{11}+C_{12}\alpha)/T)]^2}\\
        + (B_{21} + B_{22}\alpha)\left(-\frac{(C_{21} + C_{22}\alpha)}{T}\right)^2
        \frac{\exp(-(C_{21} + C_{22}\alpha)/T)}{[1-\exp(-(C_{21}+C_{22}\alpha)/T)]^2}

    Parameters
    ----------
    similarity_variable : float
        similarity variable as defined in [1]_, [mol/g]

    Notes
    -----
    A1 = 0.58, A2 = 1.25, A3 = 0.17338003, A4 = 0.014, B11 = 0.73917383,
    B12 = 8.88308889, C11 = 1188.28051, C12 = 1813.04613, B21 = 0.0483019,
    B22 = 4.35656721, C21 = 2897.01927, C22 = 5987.80407.

    References
    ----------
    .. [1] Lastovka, Vaclav, and John M. Shaw. "Predictive Correlations for
       Ideal Gas Heat Capacities of Pure Hydrocarbons and Petroleum Fractions."
       Fluid Phase Equilibria 356 (October 25, 2013): 338-370.
       doi:10.1016/j.fluid.2013.07.023.
    '''
    a = similarity_variable
    if cyclic_aliphatic:
        A1 = -0.1793547
        A2 = 3.86944439
        first = A1 + A2*a
    else:
        A1 = 0.58
        A2 = 1.25
        A3 = 0.17338003 # 803 instead of 8003 in another paper
        A4 = 0.014
        first = A2 + (A1-A2)/(1. + exp((a - A3)/A4))
        # Personal communication confirms the change
    B11 = 0.73917383
    B12 = 8.88308889
    C11 = 1188.28051
    C12 = 1813.04613
    B21 = 0.0483019
    B22 = 4.35656721
    C21 = 2897.01927
    C22 = 5987.80407
    def Lastovka_Shaw(T):
        return MW * (first + (B11 + B12*a)*((C11+C12*a)/T)**2*exp(-(C11 + C12*a)/T)/(1.-exp(-(C11+C12*a)/T))**2
                     + (B21 + B22*a)*((C21+C22*a)/T)**2*exp(-(C21 + C22*a)/T)/(1.-exp(-(C21+C22*a)/T))**2)
    return TDependentModel(Lastovka_Shaw, Tmin, Tmax,
                           integrate=LastovkaShawIntegralFunction(MW, similarity_variable, cyclic_aliphatic),
                           integrate_over_T=LastovkaShawIntegralOverTFunction(MW, similarity_variable, cyclic_aliphatic))

def LastovkaShawIntegralFunction(MW, similarity_variable, cyclic_aliphatic=False):
    r'''Return a function that calculates the integral of ideal-gas constant-pressure heat capacitiy (kJ/kmol/K) with the similarity variable concept and method as shown in [1]_.

    Parameters
    ----------
    similarity_variable : float
        similarity variable as defined in [1]_, [mol/g]

    Notes
    -----
    The model is for predicting specific enthalpy, not molar enthalpy like most other methods!
    Integral was computed with SymPy.

    References
    ----------
    .. [1] Lastovka, Vaclav, and John M. Shaw. "Predictive Correlations for
       Ideal Gas Heat Capacities of Pure Hydrocarbons and Petroleum Fractions."
       Fluid Phase Equilibria 356 (October 25, 2013): 338-370.
       doi:10.1016/j.fluid.2013.07.023.
    '''
    a = similarity_variable
    if cyclic_aliphatic:
        A1 = -0.1793547
        A2 = 3.86944439
        first = A1 + A2*a
    else:
        A1 = 0.58
        A2 = 1.25
        A3 = 0.17338003 # 803 instead of 8003 in another paper
        A4 = 0.014
        first = A2 + (A1-A2)/(1.+exp((a-A3)/A4)) # One reference says exp((a-A3)/A4)
        # Personal communication confirms the change

    B11 = 0.73917383
    B12 = 8.88308889
    C11 = 1188.28051
    C12 = 1813.04613
    B21 = 0.0483019
    B22 = 4.35656721
    C21 = 2897.01927
    C22 = 5987.80407
    def Lastovka_Shaw_integral(T):
        return MW * (T*first - (B11 + B12*a)*(-C11 - C12*a)**2/(-C11 - C12*a + (C11 
        + C12*a)*exp((-C11 - C12*a)/T)) - (B21 + B22*a)*(-C21 - C22*a)**2/(-C21 
        - C22*a + (C21 + C22*a)*exp((-C21 - C22*a)/T)))
    return Lastovka_Shaw_integral

def LastovkaShawIntegralOverTFunction(MW, similarity_variable, cyclic_aliphatic=False):
    r'''Return a function that calculates the integral over temperature of ideal-gas constant-pressure 
    heat capacitiy (kJ/kmol/K) with the similarity variable concept and method as shown in
    [1]_.

    Parameters
    ----------
    T : float
        Temperature of gas [K]
    similarity_variable : float
        similarity variable as defined in [1]_, [mol/g]

    Notes
    -----
    Integral was computed with SymPy.

    References
    ----------
    .. [1] Lastovka, Vaclav, and John M. Shaw. "Predictive Correlations for
       Ideal Gas Heat Capacities of Pure Hydrocarbons and Petroleum Fractions."
       Fluid Phase Equilibria 356 (October 25, 2013): 338-370.
       doi:10.1016/j.fluid.2013.07.023.
    '''
    a = similarity_variable
    if cyclic_aliphatic:
        A1 = -0.1793547
        A2 = 3.86944439
        first = A1 + A2*a
    else:
        A1 = 0.58
        A2 = 1.25
        A3 = 0.17338003 # 803 instead of 8003 in another paper
        A4 = 0.014
        first = A2 + (A1-A2)/(1. + exp((a - A3)/A4))

    a2 = a*a
    B11 = 0.73917383
    B12 = 8.88308889
    C11 = 1188.28051
    C12 = 1813.04613
    B21 = 0.0483019
    B22 = 4.35656721
    C21 = 2897.01927
    C22 = 5987.80407
    def Lastovka_Shaw_integral_over_T(T):
        S = (first*clog(T) + (-B11 - B12*a)*clog(cexp((-C11 - C12*a)/T) - 1.) 
            + (-B11*C11 - B11*C12*a - B12*C11*a - B12*C12*a2)/(T*cexp((-C11
            - C12*a)/T) - T) - (B11*C11 + B11*C12*a + B12*C11*a + B12*C12*a2)/T)
        S += ((-B21 - B22*a)*clog(cexp((-C21 - C22*a)/T) - 1.) + (-B21*C21 - B21*C22*a
            - B22*C21*a - B22*C22*a2)/(T*cexp((-C21 - C22*a)/T) - T) - (B21*C21
            + B21*C22*a + B22*C21*a + B22*C22*a**2)/T)
        return MW * S.real
    # There is a non-real component, but it is only a function of similariy 
    # variable and so will always cancel out.
    return Lastovka_Shaw_integral_over_T


def TRCCpModel(a0, a1, a2, a3, a4, a5, a6, a7, *, Tmin, Tmax):
    r'''Return a TDependentModel object that calculates ideal gas heat capacity (J/mol/K) using the model developed in [1]_.
    
    .. math::
        C_p = R\left(a_0 + (a_1/T^2) \exp(-a_2/T) + a_3 y^2
        + (a_4 - a_5/(T-a_7)^2 )y^j \right)

        y = \frac{T-a_7}{T+a_6} \text{ for } T > a_7 \text{ otherwise } 0

    Parameters
    ----------
    a1-a7 : float
        Coefficients

    Notes
    -----
    j is set to 8. Analytical integrals are available for this expression.

    References
    ----------
    .. [1] Kabo, G. J., and G. N. Roganov. Thermodynamics of Organic Compounds
       in the Gas State, Volume II: V. 2. College Station, Tex: CRC Press, 1994.
    '''
    def TRCCp(T):
        if T <= a7:
            y = 0.
        else:
            y = (T - a7)/(T + a6)
        return R*(a0 + (a1/T**2)*exp(-a2/T) + a3*y**2 + (a4 - a5/(T-a7)**2 )*y**8.)
    return TDependentModel(TRCCp, Tmin, Tmax,
                           integrate=TRCCpIntegralFunction(a0, a1, a2, a3, a4, a5, a6, a7),
                           integrate_over_T=TRCCpIntegralOverTFunction(a0, a1, a2, a3, a4, a5, a6, a7))


def TRCCpIntegralFunction(a0, a1, a2, a3, a4, a5, a6, a7):
    r'''Return function that integrates ideal gas heat capacity using the model developed in [1]_.
    Best used as a delta only.

    The difference in enthalpy with respect to 0 K is given by:

    .. math::
        \frac{H(T) - H^{ref}}{RT} = a_0 + a_1x(a_2)/(a_2T) + I/T + h(T)/T
        
        h(T) = (a_5 + a_7)\left[(2a_3 + 8a_4)\ln(1-y)+ \left\{a_3\left(1 + 
        \frac{1}{1-y}\right) + a_4\left(7 + \frac{1}{1-y}\right)\right\}y
        + a_4\left\{3y^2 + (5/3)y^3 + y^4 + (3/5)y^5 + (1/3)y^6\right\} 
        + (1/7)\left\{a_4 - \frac{a_5}{(a_6+a_7)^2}\right\}y^7\right]
        
        h(T) = 0 \text{ for } T \le a_7

        y = \frac{T-a_7}{T+a_6} \text{ for } T > a_7 \text{ otherwise } 0

    Parameters
    ----------
    a1-a7 : float
        Coefficients

    Notes
    -----
    Analytical integral as provided in [1]_ and verified with numerical
    integration. 

    References
    ----------
    .. [1] Kabo, G. J., and G. N. Roganov. Thermodynamics of Organic Compounds
       in the Gas State, Volume II: V. 2. College Station, Tex: CRC Press, 1994.
    '''
    first = a6 + a7
    def TRCCp_integral(Ta, Tb):
        T = Ta
        y = 0. if T <= a7 else (T - a7)/(T + a6)
        y2 = y*y
        y4 = y2*y2
        if T <= a7:
            h = 0.0
        else:
            second = (2.*a3 + 8.*a4)*log(1. - y)
            third = (a3*(1. + 1./(1. - y)) + a4*(7. + 1./(1. - y)))*y
            fourth = a4*(3.*y2 + 5./3.*y*y2 + y4 + 0.6*y4*y + 1/3.*y4*y2)
            fifth = 1/7.*(a4 - a5/((a6 + a7)**2))*y4*y2*y
            h = first*(second + third + fourth + fifth)
        H1 = (a0 + a1*exp(-a2/T)/(a2*T) + h/T)*R*T
        
        T = Tb
        y = 0. if T <= a7 else (T - a7)/(T + a6)
        y2 = y*y
        y4 = y2*y2
        if T <= a7:
            h = 0.0
        else:
            second = (2.*a3 + 8.*a4)*log(1. - y)
            third = (a3*(1. + 1./(1. - y)) + a4*(7. + 1./(1. - y)))*y
            fourth = a4*(3.*y2 + 5./3.*y*y2 + y4 + 0.6*y4*y + 1/3.*y4*y2)
            fifth = 1/7.*(a4 - a5/((a6 + a7)**2))*y4*y2*y
            h = first*(second + third + fourth + fifth)
        H2 = (a0 + a1*exp(-a2/T)/(a2*T) + h/T)*R*T
        
        return H2 - H1
    return TRCCp_integral


def TRCCpIntegralOverTFunction(a0, a1, a2, a3, a4, a5, a6, a7):
    r'''Return function that integrates ideal gas heat capacity over T (J/mol/K) using the model developed in 
    [1]_. Best used as a delta only.

    The difference in ideal-gas entropy with respect to 0 K is given by:

    .. math::
        \frac{S^\circ}{R} = J + a_0\ln T + \frac{a_1}{a_2^2}\left(1
        + \frac{a_2}{T}\right)x(a_2) + s(T)

        s(T) = \left[\left\{a_3 + \left(\frac{a_4 a_7^2 - a_5}{a_6^2}\right)
        \left(\frac{a_7}{a_6}\right)^4\right\}\left(\frac{a_7}{a_6}\right)^2
        \ln z + (a_3 + a_4)\ln\left(\frac{T+a_6}{a_6+a_7}\right)
        +\sum_{i=1}^7 \left\{\left(\frac{a_4 a_7^2 - a_5}{a_6^2}\right)\left(
        \frac{-a_7}{a_6}\right)^{6-i} - a_4\right\}\frac{y^i}{i}
        - \left\{\frac{a_3}{a_6}(a_6 + a_7) + \frac{a_5 y^6}{7a_7(a_6+a_7)}
        \right\}y\right]

        s(T) = 0 \text{ for } T \le a_7
        
        z = \frac{T}{T+a_6} \cdot \frac{a_7 + a_6}{a_7}

        y = \frac{T-a_7}{T+a_6} \text{ for } T > a_7 \text{ otherwise } 0

    Parameters
    ----------
    a1-a7 : float
        Coefficients

    Notes
    -----
    Analytical integral as provided in [1]_ and verified with numerical
    integration. 

    References
    ----------
    .. [1] Kabo, G. J., and G. N. Roganov. Thermodynamics of Organic Compounds
       in the Gas State, Volume II: V. 2. College Station, Tex: CRC Press, 1994.
    '''
    # Possible optimizations: pre-cache as much as possible.
    # If this were replaced by a cache, much of this would not need to be computed.
    def TRCCp_integral_over_T(Ta, Tb):
        T = Ta
        y = 0. if T <= a7 else (T - a7)/(T + a6)
        z = T/(T + a6)*(a7 + a6)/a7
        if T <= a7:
            s = 0.
        else:
            a72 = a7*a7
            a62 = a6*a6
            a7_a6 = a7/a6 # a7/a6
            a7_a6_2 = a7_a6*a7_a6
            a7_a6_4 = a7_a6_2*a7_a6_2
            x1 = (a4*a72 - a5)/a62 # part of third, sum
            first = (a3 + ((a4*a72 - a5)/a62)*a7_a6_4)*a7_a6_2*log(z)
            second = (a3 + a4)*log((T + a6)/(a6 + a7))
            fourth = -(a3/a6*(a6 + a7) + a5*y**6/(7.*a7*(a6 + a7)))*y
            third = np.array([(x1*(-a7_a6)**(6-i) - a4)*y**i/i for i in np.arange(1, 8)]).sum()
            s = first + second + third + fourth
        H1_ = a0*log(T) + a1/(a2*a2)*(1. + a2/T)*exp(-a2/T) + s
        
        T = Tb
        y = 0. if T <= a7 else (T - a7)/(T + a6)
        z = T/(T + a6)*(a7 + a6)/a7
        if T <= a7:
            s = 0.
        else:
            a72 = a7*a7
            a62 = a6*a6
            a7_a6 = a7/a6 # a7/a6
            a7_a6_2 = a7_a6*a7_a6
            a7_a6_4 = a7_a6_2*a7_a6_2
            x1 = (a4*a72 - a5)/a62 # part of third, sum
            first = (a3 + ((a4*a72 - a5)/a62)*a7_a6_4)*a7_a6_2*log(z)
            second = (a3 + a4)*log((T + a6)/(a6 + a7))
            fourth = -(a3/a6*(a6 + a7) + a5*y**6/(7.*a7*(a6 + a7)))*y
            third = np.array([(x1*(-a7_a6)**(6-i) - a4)*y**i/i for i in np.arange(1, 8)]).sum()
            s = first + second + third + fourth
        H2_ = a0*log(T) + a1/(a2*a2)*(1. + a2/T)*exp(-a2/T) + s
        return R*(H2_ - H1_)
        
    return TRCCp_integral_over_T
    
def PolingModel(A, B, C, D, E, Tmin, Tmax):
    R_ = R
    def Poling(T):
        return R_*(A + B*T + C*T**2 + D*T**3 + E*T**4)
    def Poling_integrate(Ta, Tb):
        H2 = (((((0.2*E)*Tb + 0.25*D)*Tb + C/3.)*Tb + 0.5*B)*Tb + A)*Tb
        H1 = (((((0.2*E)*Ta + 0.25*D)*Ta + C/3.)*Ta + 0.5*B)*Ta + A)*Ta
        return R_*(H2 - H1)
    def Poling_integrate_over_T(Ta, Tb):
        S2 = ((((0.25*E)*Tb + D/3.)*Tb + 0.5*C)*Tb + B)*Tb
        S1 = ((((0.25*E)*Ta + D/3.)*Ta + 0.5*C)*Ta + B)*Ta
        return R_*(S2-S1 + A*log(Tb/Ta))
    return TDependentModel(Poling, Tmin, Tmax, integrate=Poling_integrate,
                           integrate_over_T=Poling_integrate_over_T)

# TODO: Make cool prop model use Cp_f
# def CoolPropCpModel(Cp_f, CAS, Tmin, Tmax):
#     def coolprop(T): return PropsSI('Cp0molar', 'T', T,'P', 101325.0, CAS)
#     return TDependentModel(coolprop, Tmin, Tmax, compile=False)


# Heat capacity gas methods:
# TRCIG = 'TRC Thermodynamics of Organic Compounds in the Gas State (1994)'
# POLING = 'Poling et al. (2001)'
# POLING_CONST = 'Poling et al. (2001) constant'
# CRCSTD = 'CRC Standard Thermodynamic Properties of Chemical Substances'
# VDI_TABULAR = 'VDI Heat Atlas'
# LASTOVKA_SHAW = 'Lastovka and Shaw (2013)'
# COOLPROP = 'CoolProp'

def HeatCapacityGas(CAS, MW=None, similarity_variable=None, cyclic_aliphatic=None, models=None):
    if not models: models = []
    if CAS in TRC_gas_data:
        Tmin, Tmax, a0, a1, a2, a3, a4, a5, a6, a7, _, _, _ = TRC_gas_data[CAS]
        models.append(TRCCpModel(a0, a1, a2, a3, a4, a5, a6, a7, Tmin=Tmin, Tmax=Tmax))
    if CAS in Poling_data:
        Tmin, Tmax, A, B, C, D, E, Cpg, Cpl = Poling_data[CAS]
        if not np.isnan(a0):
            models.append(PolingModel(A, B, C, D, E, Tmin, Tmax))
        if not np.isnan(Cpg):
            models.append(ConstantTDependentModel(Cpg, Tmin, Tmax))
    if CAS in CRC_standard_data:
        Cpg = CRC_standard_data[CAS][-1]
        if not np.isnan(Cpg):
            models.append(ConstantTDependentModel(Cpg, Tmin=0, Tmax=1e6))
    if CAS in _VDISaturationDict:
        # NOTE: VDI data is for the saturation curve, i.e. at increasing
        # pressure; it is normally substantially higher than the ideal gas
        # value
        Ts, Cpgs = VDI_tabular_data(CAS, 'Cp (g)')
        models.append(InterpolatedTDependentModel(Ts, Cpgs, Tmin=Ts[0], Tmax=Ts[-1]))
    # TODO: Make cool prop model use Cp_f
    # if has_CoolProp and CAS in coolprop_dict:
    #     CP_f = coolprop_fluids[CAS]
    #     Tmin = CP_f.Tt; Tmax = CP_f.Tc
    #     models.append(CoolPropCpModel(CP_f, CAS, Tmin, Tmax))
    if MW and similarity_variable:
        models.append(LastovkaShawModel(MW, similarity_variable, cyclic_aliphatic, Tmin=Tmin, Tmax=Tmax))
    if not models:
        raise ValueError(f"no gas heat capacity models available for CAS {CAS}")
    active_model = models[0]
    return TDependentModelHandler(models, active_model)
    
### Heat capacities of liquids

def RowlinsonPolingModel(Tc, ω, Cpgm, Tmin, Tmax):
    r'''Return a TDependentModel object that calculates liquid constant-pressure heat capacity (J/mol/K) at arbitrary temperatures (K) with the [1]_ CSP method. This equation is not terrible accurate.

    .. math::
        \frac{Cp^{L} - Cp^{g}}{R} = 1.586 + \frac{0.49}{1-T_r} +
        \omega\left[ 4.2775 + \frac{6.3(1-T_r)^{1/3}}{T_r} + \frac{0.4355}{1-T_r}\right]

    Parameters
    ----------
    Tc : float
        Critical temperature of fluid [K]
    ω : float
        Acentric factor for fluid, [-]
    Cpgm : fuction(T)
        Should return constant-pressure gas heat capacity, [J/mol/K]

    Notes
    -----
    Poling compared 212 substances, and found error at 298K larger than 10%
    for 18 of them, mostly associating. Of the other 194 compounds, AARD is 2.5%.

    References
    ----------
    .. [1] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
    '''
    def Rowlinson_Poling(T):
        Tr = T/Tc
        return Cpgm(T) + R*(1.586 + 0.49/(1.-Tr) + ω*(4.2775 + 6.3*(1-Tr)**(1/3.)/Tr + 0.4355/(1.-Tr)))
    return TDependentModel(Rowlinson_Poling, Tmin, Tmax)


def RowlinsonBondiModel(Tc, ω, Cpgm, Tmin, Tmax):
    r'''Return a TDependentModel object that calculates liquid constant-pressure heat capacity (J/mol/K) at arbitrary temperatures (K) with the CSP method shown in [1]_.

    .. math::
        \frac{Cp^L - Cp^{ig}}{R} = 1.45 + 0.45(1-T_r)^{-1} + 0.25\omega
        [17.11 + 25.2(1-T_r)^{1/3}T_r^{-1} + 1.742(1-T_r)^{-1}]

    Parameters
    ----------
    Tc : float
        Critical temperature of fluid [K]
    ω : float
        Acentric factor for fluid, [-]
    Cpgm : fuction(T)
        Should return constant-pressure gas heat capacity, [J/mol/K]

    Notes
    -----
    Less accurate than `Rowlinson_Poling`.

    References
    ----------
    .. [1] Poling, Bruce E. The Properties of Gases and Liquids. 5th edition.
       New York: McGraw-Hill Professional, 2000.
    .. [2] Gesellschaft, V. D. I., ed. VDI Heat Atlas. 2nd edition.
       Berlin; New York:: Springer, 2010.
    .. [3] J.S. Rowlinson, Liquids and Liquid Mixtures, 2nd Ed.,
       Butterworth, London (1969).
    '''
    def Rowlinson_Bondi(T):
        Tr = T/Tc
        return Cpgm(T) + R*(1.45 + 0.45/(1.-Tr) + 0.25*ω*(17.11 + 25.2*(1-Tr)**(1/3.)/Tr + 1.742/(1.-Tr)))
    return TDependentModel(Rowlinson_Bondi, Tmin, Tmax)


def DadgostarShawModel(similarity_variable, MW, Tmin, Tmax):
    r'''Return a TDependentModel object that calculates liquid constant-pressure heat capacity (J/mol/K) at arbitrary temperatures (K) with the similarity variable concept and method as shown in [1]_.

    .. math::
        C_{p}^{mol} = MW \cdot C_{p}^{mass}
        C_{p}^{mass} = 24.5(a_{11}\alpha + a_{12}\alpha^2)+ (a_{21}\alpha
        + a_{22}\alpha^2)T +(a_{31}\alpha + a_{32}\alpha^2)T^2

    Parameters
    ----------
    similarity_variable : float
        similarity variable as defined in [1]_, [mol/g]

    Notes
    -----
    Many restrictions on its use.

    a11 = -0.3416; a12 = 2.2671; a21 = 0.1064; a22 = -0.3874l;
    a31 = -9.8231E-05; a32 = 4.182E-04

    References
    ----------
    .. [1] Dadgostar, Nafiseh, and John M. Shaw. "A Predictive Correlation for
       the Constant-Pressure Specific Heat Capacity of Pure and Ill-Defined
       Liquid Hydrocarbons." Fluid Phase Equilibria 313 (January 15, 2012):
       211-226. doi:10.1016/j.fluid.2011.09.015.
    '''
    a = similarity_variable
    a2 = a*a
    a11 = -0.3416
    a12 = 2.2671
    a21 = 0.1064
    a22 = -0.3874
    a31 = -9.8231E-05
    a32 = 4.182E-04
    # Didn't seem to improve the comparison; sum of errors on some
    # points included went from 65.5  to 286.
    # Author probably used more precision in their calculation.
    #    theta = 151.8675
    #    constant = 3*R*(theta/T)**2*exp(theta/T)/(exp(theta/T)-1)**2
    constant = 24.5
    first = constant*(a11*a + a12*a2)
    second = a21*a + a22*a2
    third = a31*a + a32*a2
    def Dadgostar_Shaw(T):
        return (first + second*T + third*T**2) * MW
    def Dadgostar_Shaw_integral(Ta, Tb):
        T = Ta
        T2 = T*T
        H2 = T2*T/3.*third + T2*0.5*second + T*first
        T = Tb
        T2 = T*T
        H1 = T2*T/3.*third + T2*0.5*second + T*first
        return MW*(H2 - H1)
    def Dadgostar_Shaw_integral_over_T(Ta, Tb):
        T = Ta
        S1 = T*T*0.5*third + T*second + first*log(T)
        T = Tb
        S2 = T*T*0.5*third + T*second + first*log(T)
        return MW*(S2 - S1)
    return TDependentModel(Dadgostar_Shaw, Tmin, Tmax,
                           integrate=Dadgostar_Shaw_integral,
                           integrate_over_T=Dadgostar_Shaw_integral_over_T)
        
zabransky_dict_sat_s = {}
zabransky_dict_sat_p = {}
zabransky_dict_const_s = {}
zabransky_dict_const_p = {}
zabransky_dict_iso_s = {}
zabransky_dict_iso_p = {}

# C means average heat capacity values, from less rigorous experiments
# sat means heat capacity along the saturation line
# p means constant-pressure values, 
# True means it has a splie, False otherwise
type_to_zabransky_dict = {('C', True): zabransky_dict_const_s, 
                          ('C', False):   zabransky_dict_const_p,
                          ('sat', True):  zabransky_dict_sat_s,
                          ('sat', False): zabransky_dict_sat_p,
                          ('p', True):    zabransky_dict_iso_s,
                          ('p', False):   zabransky_dict_iso_p}

with open(os.path.join(folder, 'Zabransky.tsv'), encoding='utf-8') as f:
    next(f)
    for line in f:
        values = to_num(line.strip('\n').split('\t'))
        (CAS, name, Type, uncertainty, Tmin, Tmax, a1s, a2s, a3s, a4s, a1p, a2p, a3p, a4p, a5p, a6p, Tc) = values
        spline = bool(a1s) # False if Quasypolynomial, True if spline
        d = type_to_zabransky_dict[(Type, spline)]
        if spline:
            if CAS not in d:
                d[CAS] = [(a1s, a2s, a3s, a4s, Tmin, Tmax)]
            else:
                d[CAS].append((a1s, a2s, a3s, a4s, Tmin, Tmax))
        else:
            # No duplicates for quasipolynomials
            d[CAS] = (a1p, a2p, a3p, a4p, a5p, a6p, Tc, Tmin, Tmax)

def ZabranskyQuasiPolynomialModel(Tc, a1, a2, a3, a4, a5, a6, Tmin, Tmax):
    r'''Calculates liquid heat capacity using the model developed in [1]_.

    .. math::
        \frac{C}{R}=A_1\ln(1-T_r) + \frac{A_2}{1-T_r}
        + \sum_{j=0}^m A_{j+3} T_r^j

    Parameters
    ----------
    Tc : float
        Critical temperature of fluid, [K]
    a1-a6 : float
        Coefficients

    Notes
    -----
    Used only for isobaric heat capacities, not saturation heat capacities.
    Designed for reasonable extrapolation behavior caused by using the reduced
    critical temperature. Used by the authors of [1]_ when critical temperature
    was available for the fluid. Analytical integrals are available for this expression.

    References
    ----------
    .. [1] Zabransky, M., V. Ruzicka Jr, V. Majer, and Eugene S. Domalski.
       Heat Capacity of Liquids: Critical Review and Recommended Values.
       2 Volume Set. Washington, D.C.: Amer Inst of Physics, 1996.
    '''
    Tc2 = Tc*Tc
    Tc3 = Tc2*Tc
    def Zabransky_quasi_polynomial(T):
        Tr = T/Tc
        return R*(a1*log(1-Tr) + a2/(1-Tr) + a3 + a4*Tr + a5*Tr**2 + a6*Tr**3)
    def Zabransky_quasi_polynomial_integral(Ta, Tb):
        T = Ta
        term = T - Tc
        H1 = R*(T*(T*(T*(T*a6/(4.*Tc3) + a5/(3.*Tc2)) + a4/(2.*Tc)) - a1 + a3) 
                  + T*a1*log(1. - T/Tc) - 0.5*Tc*(a1 + a2)*log(term*term))
        T = Tb
        term = Tb - Tc
        H2 = R*(T*(T*(T*(T*a6/(4.*Tc3) + a5/(3.*Tc2)) + a4/(2.*Tc)) - a1 + a3) 
                  + T*a1*log(1. - T/Tc) - 0.5*Tc*(a1 + a2)*log(term*term))
        return H2 - H1
    def Zabransky_quasi_polynomial_integral_over_T(Ta, Tb):
        T = Ta
        term = T - Tc
        logT = log(T)
        H1 = R*(a3*logT -a1*polylog2(T/Tc) - a2*(-logT + 0.5*log(term*term))
                  + T*(T*(T*a6/(3.*Tc3) + a5/(2.*Tc2)) + a4/Tc))
        T = Tb
        term = T - Tc
        logT = log(T)
        H2 = R*(a3*logT -a1*polylog2(T/Tc) - a2*(-logT + 0.5*log(term*term))
                  + T*(T*(T*a6/(3.*Tc3) + a5/(2.*Tc2)) + a4/Tc))
        return H2 - H1
    return TDependentModel(Zabransky_quasi_polynomial, Tmin, Tmax,
                           integrate=Zabransky_quasi_polynomial_integral,
                           integrate_over_T=Zabransky_quasi_polynomial_integral_over_T)


def ZabranskyCubicModel(a1, a2, a3, a4, Tmin, Tmax):
    r'''Return a TDependentModel object that calculates liquid heat capacity using the model developed in [1]_.

    .. math::
        \frac{C}{R}=\sum_{j=0}^3 A_{j+1} \left(\frac{T}{100}\right)^j

    Parameters
    ----------
    T : float
        Temperature [K]
    a1-a4 : float
        Coefficients

    Notes
    -----
    Most often form used in [1]_.
    Analytical integrals are available for this expression.

    References
    ----------
    .. [1] Zabransky, M., V. Ruzicka Jr, V. Majer, and Eugene S. Domalski.
       Heat Capacity of Liquids: Critical Review and Recommended Values.
       2 Volume Set. Washington, D.C.: Amer Inst of Physics, 1996.
    '''
    def Zabransky_cubic(T):
        T = T/100.
        return R*(((a4*T + a3)*T + a2)*T + a1)
    def Zabransky_cubic_integral(Ta, Tb):
        T = Ta/100.
        H1 = 100*R*T*(T*(T*(T*a4*0.25 + a3/3.) + a2*0.5) + a1)
        T = Tb/100.
        H2 = 100*R*T*(T*(T*(T*a4*0.25 + a3/3.) + a2*0.5) + a1)
        return H2 - H1
    def Zabransky_cubic_integral_over_T(Ta, Tb):
        T = Ta/100.
        H1 = R*(T*(T*(T*a4/3 + a3/2) + a2) + a1*log(T))
        T = Tb/100.
        H2 = R*(T*(T*(T*a4/3 + a3/2) + a2) + a1*log(T))
        return H2 - H1
    return TDependentModel(Zabransky_cubic, Tmin, Tmax,
                           integrate=Zabransky_cubic_integral,
                           integrate_over_T=Zabransky_cubic_integral_over_T)



# Heat capacity liquid methods:
ZABRANSKY_SPLINE = 'Zabransky spline, averaged heat capacity'
ZABRANSKY_QUASIPOLYNOMIAL = 'Zabransky quasipolynomial, averaged heat capacity'
ZABRANSKY_SPLINE_C = 'Zabransky spline, constant-pressure'
ZABRANSKY_QUASIPOLYNOMIAL_C = 'Zabransky quasipolynomial, constant-pressure'
ZABRANSKY_SPLINE_SAT = 'Zabransky spline, saturation'
ZABRANSKY_QUASIPOLYNOMIAL_SAT = 'Zabransky quasipolynomial, saturation'
ROWLINSON_POLING = 'Rowlinson and Poling (2001)'
ROWLINSON_BONDI = 'Rowlinson and Bondi (1969)'
DADGOSTAR_SHAW = 'Dadgostar and Shaw (2011)'


ZABRANSKY_TO_DICT = {ZABRANSKY_SPLINE: zabransky_dict_const_s,
                     ZABRANSKY_QUASIPOLYNOMIAL: zabransky_dict_const_p,
                     ZABRANSKY_SPLINE_C: zabransky_dict_iso_s,
                     ZABRANSKY_QUASIPOLYNOMIAL_C: zabransky_dict_iso_p,
                     ZABRANSKY_SPLINE_SAT: zabransky_dict_sat_s,
                     ZABRANSKY_QUASIPOLYNOMIAL_SAT: zabransky_dict_sat_p}

def HeatCapacityLiquid(CAS, Tc=None, omega=None, MW=None, similarity_variable=None, Cpgm=None):
    models = []
    if CAS in zabransky_dict_const_s:
        for args in zabransky_dict_const_s[CAS]:
            models.append(ZabranskyCubicModel(*args))
    if CAS in zabransky_dict_const_p:
        models.append(ZabranskyQuasiPolynomialModel(*zabransky_dict_const_p[CAS]))
    if CAS in zabransky_dict_iso_s:
        for args in zabransky_dict_iso_s[CAS]:
            models.append(ZabranskyCubicModel(*args))
    if CAS in zabransky_dict_iso_p:
        models.append(ZabranskyQuasiPolynomialModel(*zabransky_dict_iso_p[CAS]))
    if CAS in Poling_data:
        Tmin, Tmax, A, B, C, D, E, Cpg, Cpl = Poling_data[CAS]
        if not np.isnan(Cpg):
            models.append(ConstantTDependentModel(Cpl, Tmin, Tmax))
    if CAS in CRC_standard_data:
        Cpl = CRC_standard_data[CAS][-5]
        if not np.isnan(Cpl):
            models.append(ConstantTDependentModel(Cpl, Tmin=0, Tmax=1e6))
    # Saturation functions
    if CAS in zabransky_dict_sat_s:
        for args in zabransky_dict_sat_s[CAS]:
            models.append(ZabranskyCubicModel(*args))
    if CAS in zabransky_dict_sat_p:
        models.append(ZabranskyQuasiPolynomialModel(*zabransky_dict_sat_p[CAS]))
    if CAS in _VDISaturationDict:
        # NOTE: VDI data is for the saturation curve, i.e. at increasing
        # pressure; it is normally substantially higher than the ideal gas
        # value
        Ts, Cpls = VDI_tabular_data(CAS, 'Cp (l)')
        models.append(InterpolatedTDependentModel(Ts, Cpls, Ts[0], Ts[-1]))
    if Tc and omega and Cpgm:
        models.append(RowlinsonBondiModel(Tc, omega, Cpgm, Tmin, Tmax))
        models.append(RowlinsonPolingModel(Tc, omega, Cpgm, Tmin, Tmax))
    # TODO: Implement coolprop
    # if has_CoolProp and CAS in coolprop_dict:
    #     methods.append(COOLPROP)
    #     self.CP_f = coolprop_fluids[CAS]
    #     Tmins.append(self.CP_f.Tt); Tmaxs.append(self.CP_f.Tc)
    if MW and similarity_variable:
        models.append(DadgostarShawModel(similarity_variable, MW, Tmin, Tmax))
    active_model = models[0]
    return TDependentModelHandler(models, active_model)

### Solid

def LastovkaSolidModel(similarity_variable, MW, Tmin, Tmax):
    r'''Calculate solid constant-pressure heat capacitiy with the similarity
    variable concept and method as shown in [1]_.

    .. math::
        C_p = 3(A_1\alpha + A_2\alpha^2)R\left(\frac{\theta}{T}\right)^2
        \frac{\exp(\theta/T)}{[\exp(\theta/T)-1]^2}
        + (C_1\alpha + C_2\alpha^2)T + (D_1\alpha + D_2\alpha^2)T^2

    Parameters
    ----------
    similarity_variable : float
        similarity variable as defined in [1]_, [mol/g]

    Notes
    -----
    Many restrictions on its use. Trained on data with MW from 12.24 g/mol
    to 402.4 g/mol, C mass fractions from 61.3% to 95.2%,
    H mass fractions from 3.73% to 15.2%, N mass fractions from 0 to 15.4%,
    O mass fractions from 0 to 18.8%, and S mass fractions from 0 to 29.6%.
    Recommended for organic compounds with low mass fractions of hetero-atoms
    and especially when molar mass exceeds 200 g/mol. This model does not show
    and effects of phase transition but should not be used passed the triple
    point.

    Original model is in terms of J/g/K. Note that the model s for predicting
    mass heat capacity, not molar heat capacity like most other methods!

    A1 = 0.013183; A2 = 0.249381; theta = 151.8675; C1 = 0.026526;
    C2 = -0.024942; D1 = 0.000025; D2 = -0.000123.

    Examples
    --------
    >>> Lastovka_solid(300, 0.2139)
    1682.063629222013

    References
    ----------
    .. [1] Laštovka, Václav, Michal Fulem, Mildred Becerra, and John M. Shaw.
       "A Similarity Variable for Estimating the Heat Capacity of Solid Organic
       Compounds: Part II. Application: Heat Capacity Calculation for
       Ill-Defined Organic Solids." Fluid Phase Equilibria 268, no. 1-2
       (June 25, 2008): 134-41. doi:10.1016/j.fluid.2008.03.018.
    '''
    A1 = 0.013183
    A2 = 0.249381
    theta = 151.8675
    C1 = 0.026526
    C2 = -0.024942
    D1 = 0.000025
    D2 = -0.000123
    sim2 = similarity_variable*similarity_variable
    def Lastovka_solid(T):
        return MW * (3*(A1*similarity_variable + A2*similarity_variable**2)*R*(theta/T)**2
                     * exp(theta/T)/(exp(theta/T)-1)**2
                     + (C1*similarity_variable + C2*similarity_variable**2)*T
                     + (D1*similarity_variable + D2*similarity_variable**2)*T**2)
    def Lastovka_solid_integral(T):
        return MW * (T**3*(1000.*D1*similarity_variable/3. 
            + 1000.*D2*sim2/3.) + T*T*(500.*C1*similarity_variable 
            + 500.*C2*sim2)
            + (3000.*A1*R*similarity_variable*theta
            + 3000.*A2*R*sim2*theta)/(exp(theta/T) - 1.))
    def Lastovka_solid_integral_over_T(T):
        exp_theta_T = exp(theta/T)
        return MW * (-3000.*R*similarity_variable*(A1 + A2*similarity_variable)*log(exp_theta_T - 1.) 
            + T**2*(500.*D1*similarity_variable + 500.*D2*sim2)
            + T*(1000.*C1*similarity_variable + 1000.*C2*sim2)
            + (3000.*A1*R*similarity_variable*theta 
            + 3000.*A2*R*sim2*theta)/(T*exp_theta_T - T) 
            + (3000.*A1*R*similarity_variable*theta 
            + 3000.*A2*R*sim2*theta)/T)
    return TDependentModel(Lastovka_solid, Tmin, Tmax,
                           integrate=Lastovka_solid_integral,
                           integrate_over_T=Lastovka_solid_integral_over_T)
    
def Perry151Model(Const, Lin, Quad, Quadinv, Tmin, Tmax):
    def Perry_151(T):
        return (Const + Lin*T + Quad/T**2 +Quadinv*T**2)*calorie
    def Perry_151_integral(Ta, Tb):
        H1 = (Const*Ta + 0.5*Lin*Ta**2 - Quadinv/Ta + Quad*Ta**3/3.)
        H2 = (Const*Tb + 0.5*Lin*Tb**2 - Quadinv/Tb + Quad*Tb**3/3.)
        return (H2 - H1)*calorie
    def Perry_151_integral_over_T(Ta, Tb):
        S1 = Const*log(Ta) + Lin*Ta - Quadinv/(2.*Ta**2) + 0.5*Quad*Ta**2
        S2 = Const*log(Tb) + Lin*Tb - Quadinv/(2.*Tb**2) + 0.5*Quad*Tb**2
        return (S2 - S1)*calorie
    return TDependentModel(Perry_151, Tmin, Tmax,
                           integrate=Perry_151_integral,
                           integrate_over_T=Perry_151_integral_over_T)

# Heat capacity solid methods
LASTOVKA_S = 'Lastovka, Fulem, Becerra and Shaw (2008)'
PERRY151 = '''Perry's Table 2-151'''

def HeatCapacitySolid(CAS, similarity_variable, MW):
    models = []
    Tmin = 0
    Tmax = 2000
    if CAS in _PerryI:
        vals = _PerryI[CAS]
        if 'c' in vals:
            c = vals['c']
            Tmin = c['Tmin']
            Tmax = c['Tmax']
            models.append(Perry151Model(c['Const'], c['Lin'], c['Quad'], c['Quadinv'], Tmin, Tmax))
    if CAS in CRC_standard_data:
        Cpc = CRC_standard_data[CAS][3]
        if not np.isnan(Cpc):
            models.append(ConstantTDependentModel(float(Cpc), 200, 350))
    if similarity_variable and MW:
        models.append(LastovkaSolidModel(similarity_variable, MW, Tmin, Tmax))
    active_model = models[0]
    return TDependentModelHandler(models, active_model)
    

