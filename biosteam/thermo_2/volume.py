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

__all__ = ['COSTALD_data', 'SNM0_data', 'Perry_l_data', 'CRC_inorg_l_data', 
'VDI_PPDS_2',
'CRC_inorg_l_const_data', 'CRC_inorg_s_const_data', 'CRC_virial_data', 
'Yen_Woods_saturation', 'Rackett', 'Yamada_Gunn', 'Townsend_Hales', 
'Bhirud_normal', 'COSTALD', 'Campbell_Thodos', 'SNM0', 'CRC_inorganic', 
'volume_liquid_methods', 'volume_liquid_methods_P', 'VolumeLiquid', 
'COSTALD_compressed', 'Amgat', 'Rackett_mixture', 'COSTALD_mixture', 
'ideal_gas', 'volume_gas_methods', 'VolumeGas', 
'volume_gas_mixture_methods', 'volume_solid_mixture_methods', 'Goodman',
 'volume_solid_methods', 'VolumeSolid',
'VolumeLiquidMixture', 'VolumeGasMixture', 'VolumeSolidMixture']

import os
import numpy as np
from scipy.interpolate import interp1d
import pandas as pd
from numba import njit
from .utils import R
from .utils import log, exp
from .utils import Vm_to_rho, rho_to_Vm, mixing_simple, none_and_length_check
from .virial import BVirial_Pitzer_Curl, BVirial_Abbott, BVirial_Tsonopoulos, BVirial_Tsonopoulos_extended
from .miscdata import _VDISaturationDict, VDI_tabular_data
from .dippr import ModelEQ105
from .electrochem import _Laliberte_Density_ParametersDict, LaliberteDensityModel
from .coolprop import has_CoolProp, PropsSI, PhaseSI, coolprop_fluids, coolprop_dict, CoolProp_T_dependent_property
from .thermo_model import TPDependentModel, TInterpolatedTPDependentModel, ConstantTPDependentModel
from .model_handler import TPDependentModelHandler
from .eos import PR78

folder = os.path.join(os.path.dirname(__file__), 'Density')

COSTALD_data = pd.read_csv(os.path.join(folder, 'COSTALD Parameters.tsv'),
                           sep='\t', index_col=0)

SNM0_data = pd.read_csv(os.path.join(folder, 'Mchaweh SN0 deltas.tsv'),
                        sep='\t', index_col=0)

Perry_l_data = pd.read_csv(os.path.join(folder, 'Perry Parameters 105.tsv'),
                           sep='\t', index_col=0)
_Perry_l_data_values = Perry_l_data.values

VDI_PPDS_2 = pd.read_csv(os.path.join(folder, 'VDI PPDS Density of Saturated Liquids.tsv'),
                          sep='\t', index_col=0)
_VDI_PPDS_2_values = VDI_PPDS_2.values


CRC_inorg_l_data = pd.read_csv(os.path.join(folder, 'CRC Inorganics densties of molten compounds and salts.tsv'),
                               sep='\t', index_col=0)
_CRC_inorg_l_data_values = CRC_inorg_l_data.values

CRC_inorg_l_const_data = pd.read_csv(os.path.join(folder, 'CRC Liquid Inorganic Constant Densities.tsv'),
                                     sep='\t', index_col=0)

CRC_inorg_s_const_data = pd.read_csv(os.path.join(folder, 'CRC Solid Inorganic Constant Densities.tsv'),
                                     sep='\t', index_col=0)

CRC_virial_data = pd.read_csv(os.path.join(folder, 'CRC Virial polynomials.tsv'),
                              sep='\t', index_col=0)
_CRC_virial_data_values = CRC_virial_data.values

### Critical-properties based

def YenWoodsModel(Tc, Vc, Zc, Tmin, Tmax, Pmin=100., Pmax=10e6):
    r'''Return a TPDependentModel object that calculates saturation liquid volume [m^3/mol],
    using the Yen and Woods [1]_ CSP method and a chemical's critical properties.

    The molar volume of a liquid is given by:

    .. math::
        Vc/Vs = 1 + A(1-T_r)^{1/3} + B(1-T_r)^{2/3} + D(1-T_r)^{4/3}

        D = 0.93-B

        A = 17.4425 - 214.578Z_c + 989.625Z_c^2 - 1522.06Z_c^3

        B = -3.28257 + 13.6377Z_c + 107.4844Z_c^2-384.211Z_c^3
        \text{ if } Zc \le 0.26

        B = 60.2091 - 402.063Z_c + 501.0 Z_c^2 + 641.0 Z_c^3
        \text{ if } Zc \ge 0.26


    Parameters
    ----------
    Tc : float
        Critical temperature of fluid [K]
    Vc : float
        Critical volume of fluid [m^3/mol]
    Zc : float
        Critical compressibility of fluid, [-]

    Notes
    -----
    Original equation was in terms of density, but it is converted here.

    No example has been found, nor are there points in the article. However,
    it is believed correct. For compressed liquids with the Yen-Woods method,
    see the `YenWoods_compressed` function.

    References
    ----------
    .. [1] Yen, Lewis C., and S. S. Woods. "A Generalized Equation for Computer
       Calculation of Liquid Densities." AIChE Journal 12, no. 1 (1966):
       95-99. doi:10.1002/aic.690120119
    '''
    A = 17.4425 - 214.578*Zc + 989.625*Zc**2 - 1522.06*Zc**3
    if Zc <= 0.26:
        B = -3.28257 + 13.6377*Zc + 107.4844*Zc**2 - 384.211*Zc**3
    else:
        B = 60.2091 - 402.063*Zc + 501.0*Zc**2 + 641.0*Zc**3
    D = 0.93 - B
    def Yen_Woods_saturation(T, P):
        Tr = T/Tc
        return Vc/(1 + A*(1-Tr)**(1/3.) + B*(1-Tr)**(2/3.) + D*(1-Tr)**(4./3.))
    return TPDependentModel(Yen_Woods_saturation, Tmin, Tmax, Pmin, Pmax)

def RackettModel(Tc, Pc, Zc, Tmin, Tmax, Pmin=100., Pmax=10e6):
    r'''Return a TPDependentModel object that calculates saturation liquid volume [m^3/mol], using Rackett CSP method and
    critical properties.

    The molar volume of a liquid is given by:

    .. math::
        V_s = \frac{RT_c}{P_c}{Z_c}^{[1+(1-{T/T_c})^{2/7} ]}

    Units are all currently in m^3/mol - this can be changed to kg/m^3

    Parameters
    ----------
    Tc : float
        Critical temperature of fluid [K]
    Pc : float
        Critical pressure of fluid [Pa]
    Zc : float
        Critical compressibility of fluid, [-]

    Notes
    -----
    Units are dependent on gas constant R, imported from scipy
    According to Reid et. al, underpredicts volume for compounds with Zc < 0.22

    Examples
    --------
    Propane, example from the API Handbook

    >>> Vm_to_rho(Rackett(272.03889, 369.83, 4248000.0, 0.2763), 44.09562)
    531.3223212651092

    References
    ----------
    .. [1] Rackett, Harold G. "Equation of State for Saturated Liquids."
       Journal of Chemical & Engineering Data 15, no. 4 (1970): 514-517.
       doi:10.1021/je60047a012
    '''
    def Rackett(T, P):
        return R*Tc/Pc*Zc**(1 + (1 - T/Tc)**(2/7.))
    return TPDependentModel(Rackett, Tmin, Tmax, Pmin, Pmax)


def YamadaGunnModel(Tc, Pc, omega, Tmin, Tmax, Pmin=100., Pmax=10e6):
    r'''Return a TPDependentModel object that calculates saturation liquid volume [m^3/mol], using Yamada and Gunn CSP method
    and a chemical's critical properties and acentric factor.

    The molar volume of a liquid is given by:

    .. math::
        V_s = \frac{RT_c}{P_c}{(0.29056-0.08775\omega)}^{[1+(1-{T/T_c})^{2/7}]}

    Units are in m^3/mol.

    Parameters
    ----------
    Tc : float
        Critical temperature of fluid [K]
    Pc : float
        Critical pressure of fluid [Pa]
    omega : float
        Acentric factor for fluid, [-]

    Notes
    -----
    This equation is an improvement on the Rackett equation.
    This is often presented as the Rackett equation.
    The acentric factor is used here, instead of the critical compressibility
    A variant using a reference fluid also exists

    References
    ----------
    .. [1] Gunn, R. D., and Tomoyoshi Yamada. "A Corresponding States
        Correlation of Saturated Liquid Volumes." AIChE Journal 17, no. 6
        (1971): 1341-45. doi:10.1002/aic.690170613
    .. [2] Yamada, Tomoyoshi, and Robert D. Gunn. "Saturated Liquid Molar
        Volumes. Rackett Equation." Journal of Chemical & Engineering Data 18,
        no. 2 (1973): 234-36. doi:10.1021/je60057a006
    '''
    constant = R*Tc/Pc*(0.29056 - 0.08775*omega)
    def Yamada_Gunn(T, P):
        return constant**(1 + (1 - T/Tc)**(2/7.))
    return TPDependentModel(Yamada_Gunn, Tmin, Tmax, Pmin, Pmax)


def TownsendHalesModel(Tc, Vc, omega, Tmin, Tmax, Pmin=100., Pmax=10e6):
    r'''Return a TPDependentModel object that calculates saturation liquid volume [m^3/mol], using the Townsend and Hales CSP method as modified from the original Riedel equation. Uses chemical critical volume and temperature, as well as acentric factor

    The density of a liquid is given by:

    .. math::
        Vs = V_c/\left(1+0.85(1-T_r)+(1.692+0.986\omega)(1-T_r)^{1/3}\right)

    Parameters
    ----------
    Tc : float
        Critical temperature of fluid [K]
    Vc : float
        Critical volume of fluid [m^3/mol]
    omega : float
        Acentric factor for fluid, [-]
        
    Notes
    -----
    The requirement for critical volume and acentric factor requires all data.
    References
    ----------
    .. [1] Hales, J. L, and R Townsend. "Liquid Densities from 293 to 490 K of
       Nine Aromatic Hydrocarbons." The Journal of Chemical Thermodynamics
       4, no. 5 (1972): 763-72. doi:10.1016/0021-9614(72)90050-X
    '''
    def Townsend_Hales(T, P):
        Tr = T/Tc
        return Vc/(1 + 0.85*(1-Tr) + (1.692 + 0.986*omega)*(1-Tr)**(1/3.))
    return TPDependentModel(Townsend_Hales, Tmin, Tmax, Pmin, Pmax)


Bhirud_normal_Trs = [0.98, 0.982, 0.984, 0.986, 0.988, 0.99, 0.992, 0.994,
            0.996, 0.998, 0.999, 1]
Bhirud_normal_lnU0s = [-1.6198, -1.604, -1.59, -1.578, -1.564, -1.548, -1.533,
              -1.515, -1.489, -1.454, -1.425, -1.243]
Bhirud_normal_lnU1 = [-0.4626, -0.459, -0.451, -0.441, -0.428, -0.412, -0.392,
              -0.367, -0.337, -0.302, -0.283, -0.2629]
Bhirud_normal_lnU0_interp = interp1d(Bhirud_normal_Trs, Bhirud_normal_lnU0s, kind='cubic')
Bhirud_normal_lnU1_interp = interp1d(Bhirud_normal_Trs, Bhirud_normal_lnU1, kind='cubic')


def BhirudNormalModel(Tc, Pc, omega, Tmin, Tmax, Pmin=100., Pmax=10e6):
    r'''Return a TPDependentModel object that calculates saturation liquid volume [m^3/mol] using the Bhirud [1]_ CSP method.
    Uses Critical temperature and pressure and acentric factor.

    The density of a liquid is given by:

    .. math::
        &\ln \frac{P_c}{\rho RT} = \ln U^{(0)} + \omega\ln U^{(1)}

        &\ln U^{(0)} = 1.396 44 - 24.076T_r+ 102.615T_r^2
        -255.719T_r^3+355.805T_r^4-256.671T_r^5 + 75.1088T_r^6

        &\ln U^{(1)} = 13.4412 - 135.7437 T_r + 533.380T_r^2-
        1091.453T_r^3+1231.43T_r^4 - 728.227T_r^5 + 176.737T_r^6

    Parameters
    ----------
    Tc : float
        Critical temperature of fluid [K]
    Pc : float
        Critical pressure of fluid [Pa]
    omega : float
        Acentric factor for fluid, [-]

    Notes
    -----
    Claimed inadequate by others.

    An interpolation table for ln U values are used from Tr = 0.98 - 1.000.
    Has terrible behavior at low reduced temperatures.

    References
    ----------
    .. [1] Bhirud, Vasant L. "Saturated Liquid Densities of Normal Fluids."
       AIChE Journal 24, no. 6 (November 1, 1978): 1127-31.
       doi:10.1002/aic.690240630
    '''
    def Bhirud_normal(T, P):
        Tr = T/Tc
        if Tr <= 0.98:
            lnU0 = 1.39644 - 24.076*Tr + 102.615*Tr**2 - 255.719*Tr**3 \
                + 355.805*Tr**4 - 256.671*Tr**5 + 75.1088*Tr**6
            lnU1 = 13.4412 - 135.7437*Tr + 533.380*Tr**2-1091.453*Tr**3 \
                + 1231.43*Tr**4 - 728.227*Tr**5 + 176.737*Tr**6
        elif Tr > 1:
            raise Exception('Critical phase, correlation does not apply')
        else:
            lnU0 = Bhirud_normal_lnU0_interp(Tr)
            lnU1 = Bhirud_normal_lnU1_interp(Tr)
    
        Unonpolar = exp(lnU0 + omega*lnU1)
        return (Unonpolar*R*T)/Pc
    return TPDependentModel(Bhirud_normal, Tmin, Tmax, Pmin, Pmax)


def CostaldModel(Tc, Vc, omega, Tmin, Tmax, Pmin=100., Pmax=10e6):
    r'''Return a TPDependentModel object that calculates saturation liquid volume [m^3/mol] using the Costald CSP method.

    A popular and accurate estimation method. If possible, fit parameters are
    used; alternatively critical properties work well.

    The density of a liquid is given by:

    .. math::
        V_s=V^*V^{(0)}[1-\omega_{SRK}V^{(\delta)}]

        V^{(0)}=1-1.52816(1-T_r)^{1/3}+1.43907(1-T_r)^{2/3}
        - 0.81446(1-T_r)+0.190454(1-T_r)^{4/3}

        V^{(\delta)}=\frac{-0.296123+0.386914T_r-0.0427258T_r^2-0.0480645T_r^3}
        {T_r-1.00001}

    Units are that of critical or fit constant volume.

    Parameters
    ----------
    Tc : float
        Critical temperature of fluid [K]
    Vc : float
        Critical volume of fluid [m^3/mol].
        This parameter is alternatively a fit parameter
    omega : float
        (ideally SRK) Acentric factor for fluid, [-]
        This parameter is alternatively a fit parameter.

    Notes
    -----
    196 constants are fit to this function in [1]_.
    Range: 0.25 < Tr < 0.95, often said to be to 1.0

    This function has been checked with the API handbook example problem.

    References
    ----------
    .. [1] Hankinson, Risdon W., and George H. Thomson. "A New Correlation for
       Saturated Densities of Liquids and Their Mixtures." AIChE Journal
       25, no. 4 (1979): 653-663. doi:10.1002/aic.690250412
    '''
    def Costald(T, P):
        Tr = T/Tc
        V_delta = (-0.296123 + 0.386914*Tr - 0.0427258*Tr**2
            - 0.0480645*Tr**3)/(Tr - 1.00001)
        V_0 = 1 - 1.52816*(1-Tr)**(1/3.) + 1.43907*(1-Tr)**(2/3.) \
            - 0.81446*(1-Tr) + 0.190454*(1-Tr)**(4/3.)
        return Vc*V_0*(1-omega*V_delta)
    return TPDependentModel(Costald, Tmin, Tmax, Pmin, Pmax)


def CampbellThodosModel(Tb, Tc, Pc, M, dipole=None, hydroxyl=False, Pmin=100., Pmax=10e6, *, Tmin, Tmax,):
    r'''Return a TPDependentModel object that calculates saturation liquid volume [m^3/mol] using the Campbell-Thodos [1]_
    CSP method.

    An old and uncommon estimation method.

    .. math::
        V_s = \frac{RT_c}{P_c}{Z_{RA}}^{[1+(1-T_r)^{2/7}]}

        Z_{RA} = \alpha + \beta(1-T_r)

        \alpha = 0.3883-0.0179s

        s = T_{br} \frac{\ln P_c}{(1-T_{br})}

        \beta = 0.00318s-0.0211+0.625\Lambda^{1.35}

        \Lambda = \frac{P_c^{1/3}} { M^{1/2} T_c^{5/6}}

    For polar compounds:

    .. math::
        \theta = P_c \mu^2/T_c^2

        \alpha = 0.3883 - 0.0179s - 130540\theta^{2.41}

        \beta = 0.00318s - 0.0211 + 0.625\Lambda^{1.35} + 9.74\times
        10^6 \theta^{3.38}

    Polar Combounds with hydroxyl groups (water, alcohols)

    .. math::
        \alpha = \left[0.690T_{br} -0.3342 + \frac{5.79\times 10^{-10}}
        {T_{br}^{32.75}}\right] P_c^{0.145}

        \beta = 0.00318s - 0.0211 + 0.625 \Lambda^{1.35} + 5.90\Theta^{0.835}

    Parameters
    ----------
    Tb : float
        Boiling temperature of the fluid [K]
    Tc : float
        Critical temperature of fluid [K]
    Pc : float
        Critical pressure of fluid [Pa]
    M : float
        Molecular weight of the fluid [g/mol]
    dipole : float, optional
        Dipole moment of the fluid [debye]
    hydroxyl : bool, optional
        Swith to use the hydroxyl variant for polar fluids
        
    Notes
    -----
    If a dipole is provided, the polar chemical method is used.
    The paper is an excellent read.
    Pc is internally converted to atm.

    References
    ----------
    .. [1] Campbell, Scott W., and George Thodos. "Prediction of Saturated
       Liquid Densities and Critical Volumes for Polar and Nonpolar
       Substances." Journal of Chemical & Engineering Data 30, no. 1
       (January 1, 1985): 102-11. doi:10.1021/je00039a032.
    '''
    Pc = Pc/101325.
    def Campbell_Thodos(T, P):
        Tr = T/Tc
        Tbr = Tb/Tc
        s = Tbr * log(Pc)/(1-Tbr)
        Lambda = Pc**(1/3.)/(M**0.5*Tc**(5/6.))
        alpha = 0.3883 - 0.0179*s
        beta = 0.00318*s - 0.0211 + 0.625*Lambda**(1.35)
        if dipole:
            theta = Pc*dipole**2/Tc**2
            alpha -= 130540 * theta**2.41
            beta += 9.74E6 * theta**3.38
        if hydroxyl:
            beta = 0.00318*s - 0.0211 + 0.625*Lambda**(1.35) + 5.90*theta**0.835
            alpha = (0.69*Tbr - 0.3342 + 5.79E-10/Tbr**32.75)*Pc**0.145
        Zra = alpha + beta*(1-Tr)
        Vs = R*Tc/(Pc*101325)*Zra**(1+(1-Tr)**(2/7.))
        return Vs
    return TPDependentModel(Campbell_Thodos, Tmin, Tmax, Pmin, Pmax)


def SNM0Model(Tc, Vc, omega, delta_SRK=None, *, Tmin, Tmax, Pmin=100., Pmax=10e6):
    r'''Return a TPDependentModel object that calculates saturation liquid volume [m^3/mol] using the Mchaweh, Moshfeghian
    model [1]_. Designed for simple calculations.

    .. math::
        V_s = V_c/(1+1.169\tau^{1/3}+1.818\tau^{2/3}-2.658\tau+2.161\tau^{4/3}

        \tau = 1-\frac{(T/T_c)}{\alpha_{SRK}}

        \alpha_{SRK} = [1 + m(1-\sqrt{T/T_C}]^2

        m = 0.480+1.574\omega-0.176\omega^2

    If the fit parameter `delta_SRK` is provided, the following is used:

    .. math::
        V_s = V_C/(1+1.169\tau^{1/3}+1.818\tau^{2/3}-2.658\tau+2.161\tau^{4/3})
        /\left[1+\delta_{SRK}(\alpha_{SRK}-1)^{1/3}\right]

    Parameters
    ----------
    Tc : float
        Critical temperature of fluid [K]
    Vc : float
        Critical volume of fluid [m^3/mol]
    omega : float
        Acentric factor for fluid, [-]
    delta_SRK : float, optional
        Fitting parameter [-]
    
    Notes
    -----
    73 fit parameters have been gathered from the article.

    References
    ----------
    .. [1] Mchaweh, A., A. Alsaygh, Kh. Nasrifar, and M. Moshfeghian.
       "A Simplified Method for Calculating Saturated Liquid Densities."
       Fluid Phase Equilibria 224, no. 2 (October 1, 2004): 157-67.
       doi:10.1016/j.fluid.2004.06.054
    '''
    m = 0.480 + 1.574*omega - 0.176*omega*omega
    def SNM0(T, P):
        Tr = T/Tc
        alpha_SRK = (1. + m*(1. - Tr**0.5))**2
        tau = 1. - Tr/alpha_SRK
        rho0 = 1. + 1.169*tau**(1/3.) + 1.818*tau**(2/3.) - 2.658*tau + 2.161*tau**(4/3.)
        V0 = 1./rho0
        if not delta_SRK:
            return Vc*V0
        else:
            return Vc*V0/(1. + delta_SRK*(alpha_SRK - 1.)**(1/3.))
    return TPDependentModel(SNM0, Tmin, Tmax, Pmin, Pmax)


def CRCInorganicModel(rho0, k, Tm, MW, Tmin, Tmax, Pmin=0, Pmax=1e15):
    r'''Return a TPDependentModel object that calculates saturation liquid volume [m^3/mol] of a molten element or salt at temperature
    above the melting point. Some coefficients are given nearly up to the
    boiling point.

    The mass density (g/ml) of the inorganic liquid is given by:

    .. math::
        \rho = \rho_{0} - k(T-T_m)

    Parameters
    ----------
    rho0 : float
        Mass density of the liquid at Tm, [kg/m^3]
    k : float
        Linear temperature dependence of the mass density, [kg/m^3/K]
    Tm : float
        The normal melting point, used in the correlation [K]
    MW : float
        Molecular weight [g/mol]
    Notes
    -----
    This linear form is useful only in small temperature ranges.
    Coefficients for one compound could be used to predict the temperature
    dependence of density of a similar compound.

    References
    ----------
    .. [1] Haynes, W.M., Thomas J. Bruno, and David R. Lide. CRC Handbook of
        Chemistry and Physics, 95E. [Boca Raton, FL]: CRC press, 2014.
    '''
    f = MW/1000
    def CRC_inorganic(T, P): return (rho0 - k*(T-Tm))*f
    return TPDependentModel(CRC_inorganic, Tmin, Tmax)

def VDI_PPDSModel(A, B, C, D, Tc, rhoc, MW, Tmin, Tmax, Pmin=100., Pmax=10e6):
    k = MW / 1e9
    Vc = rhoc * k
    A *= k
    B *= k
    C *= k
    D *= k
    def VDI_PPDS(T, P):
        tau = 1. - T/Tc
        return Vc + A*tau**0.35 + B*tau**(2/3.) + C*tau + D*tau**(4/3.)
    return TPDependentModel(VDI_PPDS, Tmin, Tmax, Pmin, Pmax)


COOLPROP = 'COOLPROP'
PERRYDIPPR = 'PERRYDIPPR'
VDI_PPDS = 'VDI_PPDS'
MMSNM0 = 'MMSNM0'
MMSNM0FIT = 'MMSNM0FIT'
VDI_TABULAR = 'VDI_TABULAR'
HTCOSTALD = 'HTCOSTALD'
HTCOSTALDFIT = 'HTCOSTALDFIT'
COSTALD_COMPRESSED = 'COSTALD_COMPRESSED'
RACKETT = 'RACKETT'
RACKETTFIT = 'RACKETTFIT'
YEN_WOODS_SAT = 'YEN_WOODS_SAT'
YAMADA_GUNN = 'YAMADA_GUNN'
BHIRUD_NORMAL = 'BHIRUD_NORMAL'
TOWNSEND_HALES = 'TOWNSEND_HALES'
CAMPBELL_THODOS = 'CAMPBELL_THODOS'
EOS = 'EOS'


CRC_INORG_L = 'CRC_INORG_L'
CRC_INORG_L_CONST = 'CRC_INORG_L_CONST'

volume_liquid_methods = [PERRYDIPPR, VDI_PPDS, COOLPROP, MMSNM0FIT, VDI_TABULAR,
                         HTCOSTALDFIT, RACKETTFIT, CRC_INORG_L,
                         CRC_INORG_L_CONST, MMSNM0, HTCOSTALD,
                         YEN_WOODS_SAT, RACKETT, YAMADA_GUNN,
                         BHIRUD_NORMAL, TOWNSEND_HALES, CAMPBELL_THODOS]
'''Holds all low-pressure methods available for the VolumeLiquid class, for use
in iterating over them.'''

volume_liquid_methods_P = [COOLPROP, COSTALD_COMPRESSED, EOS]
'''Holds all high-pressure methods available for the VolumeLiquid class, for
use in iterating over them.'''


def VolumeLiquid(CAS, MW=None, Tb=None, Tc=None, Pc=None, Vc=None, Zc=None,
                 omega=None, dipole=None, Psat=None, eos=None):
    models = []
    # if has_CoolProp and self.CASRN in coolprop_dict:
    #     methods.append(COOLPROP); methods_P.append(COOLPROP)
    #     self.CP_f = coolprop_fluids[self.CASRN]
    #     Tmins.append(self.CP_f.Tt); Tmaxs.append(self.CP_f.Tc)
    Tmin = 50
    Tmax = 1000
    if CAS in CRC_inorg_l_data.index:
        _, MW, rho, k, Tm, Tmax = _CRC_inorg_l_data_values[CRC_inorg_l_data.index.get_loc(CAS)].tolist()
        Tmin = Tm
        models.append(CRCInorganicModel(rho, k, Tm, MW, Tmin, Tmax))
    if CAS in Perry_l_data.index:
        _, C1, C2, C3, C4, Tmin, Tmax = _Perry_l_data_values[Perry_l_data.index.get_loc(CAS)].tolist()
        models.append(ModelEQ105(C1, C2, C3, C4, Tmin, Tmax, vol=True))
    if Tc and Pc and CAS in COSTALD_data.index and not np.isnan(COSTALD_data.at[CAS, 'Z_RA']):
        Zc_ = float(COSTALD_data.at[CAS, 'Z_RA'])
        models.append(RackettModel(Tc, Pc, Zc_, Tmin, Tmax))
        # Roughly data at STP; not guaranteed however; not used for Trange
    if Tc and CAS in COSTALD_data.index:
        Vc = float(COSTALD_data.at[CAS, 'Vchar'])
        omega = float(COSTALD_data.at[CAS, 'omega_SRK'])
        Tmin = 0
        Tmax = Tc
        models.append(CostaldModel(Tc, Vc, omega, Tmin, Tmax))
    if CAS in VDI_PPDS_2.index:
        _, MW, Tc_, rhoc, A, B, C, D = _VDI_PPDS_2_values[VDI_PPDS_2.index.get_loc(CAS)].tolist()
        Tmin = 0. # TODO: Add better Tmin
        models.append(VDI_PPDSModel(A, B, C, D, Tc_, rhoc, MW, Tmin, Tc_))
    if CAS in _VDISaturationDict:
        Ts, Vls = VDI_tabular_data(CAS, 'Volume (l)')
        models.append(TInterpolatedTPDependentModel(Ts, Vls, Tmin=Ts[0], Tmax=Ts[-1]))
    if all((Tc, Pc, omega)):
        models.insert(0, CostaldCompressedModel(Psat, Tc, Pc, omega, models[0].evaluate, 50, 500))
        if eos: models.insert(0, EOS)
    if all((Tc, Vc, Zc)):
        models.append(YenWoodsModel(Tc, Vc, Zc, 0, Tc))
    if all((Tc, Pc, Zc)):
        models.append(RackettModel(Tc, Pc, Zc, 0, Tc))
    if all((Tc, Pc, omega)):
        models.append(YamadaGunnModel(Tc, Pc, omega, 0, Tc))
        models.append(BhirudNormalModel(Tc, Pc, omega, 0, Tc))
    if all((Tc, Vc, omega)):
        models.append(TownsendHalesModel(Tc, Vc, omega, 0, Tc))
        models.append(RackettModel(Tc, Vc, omega, 0, Tc))
        if CAS in SNM0_data.index:
            SNM0_delta_SRK = float(SNM0_data.at[CAS, 'delta_SRK'])
            models.append(SNM0Model(Tc, Vc, omega, SNM0_delta_SRK))
        models.append(SNM0Model(Tc, Vc, omega, Tmin=0, Tmax=Tc))
    if all((Tc, Vc, omega, Tb, MW)):
        models.append(CampbellThodosModel(Tb, Tc, Pc, MW, Tmin=0, Tmax=Tc))
    if CAS in CRC_inorg_l_const_data.index:
        Vl = float(CRC_inorg_l_const_data.at[CAS, 'Vm'])
        models.append(ConstantTPDependentModel(Vl, Tmin, Tmax))
    return TPDependentModelHandler(models)
        


def CostaldCompressedModel(Psat, Tc, Pc, omega, Vs, Tmin, Tmax, Pmin=0, Pmax=1e15):
    r'''Return a TPDependentModel object that calculates saturation liquid volume [m^3/mol] using the COSTALD [1]_ CSP method and a chemical's critical properties.

    The molar volume of a liquid is given by:

    .. math::
        V = V_s\left( 1 - C \ln \frac{B + P}{B + P^{sat}}\right)

        \frac{B}{P_c} = -1 + a\tau^{1/3} + b\tau^{2/3} + d\tau + e\tau^{4/3}

        e = \exp(f + g\omega_{SRK} + h \omega_{SRK}^2)

        C = j + k \omega_{SRK}

    Parameters
    ----------
    Psat : float
        Saturation pressure of the fluid [Pa]
    Tc : float
        Critical temperature of fluid [K]
    Pc : float
        Critical pressure of fluid [Pa]
    omega : float
        (ideally SRK) Acentric factor for fluid, [-]
        This parameter is alternatively a fit parameter.
    Vs : float
        Saturation liquid volume, [m^3/mol]

    Returns
    -------
    V_dense : float
        High-pressure liquid volume, [m^3/mol]

    Notes
    -----
    Original equation was in terms of density, but it is converted here.

    The example is from DIPPR, and exactly correct.
    This is DIPPR Procedure 4C: Method for Estimating the Density of Pure
    Organic Liquids under Pressure.

    References
    ----------
    .. [1]  Thomson, G. H., K. R. Brobst, and R. W. Hankinson. "An Improved
       Correlation for Densities of Compressed Liquids and Liquid Mixtures."
       AIChE Journal 28, no. 4 (July 1, 1982): 671-76. doi:10.1002/aic.690280420
    '''
    a = -9.070217
    b = 62.45326
    d = -135.1102
    f = 4.79594
    g = 0.250047
    h = 1.14188
    j = 0.0861488
    k = 0.0344483
    e = exp(f + g*omega + h*omega**2)
    C = j + k*omega
    def Costald_compressed(T, P):
        tau = 1 - T/Tc
        B = Pc*(-1 + a*tau**(1/3.) + b*tau**(2/3.) + d*tau + e*tau**(4/3.))
        return Vs(T, P)*(1 - C*log((B + P)/(B + Psat(T))))
    return TPDependentModel(Costald_compressed, Tmin, Tmax, Pmin, Pmax, compile=False)


### Gases
@njit
def ideal_gas(T, P):
    r'''Calculates ideal gas molar volume.
    The molar volume of an ideal gas is given by:

    .. math::
        V = \frac{RT}{P}

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    P : float
        Pressure of fluid [Pa]

    Returns
    -------
    V : float
        Gas volume, [m^3/mol]

    Examples
    --------
    >>> ideal_gas(298.15, 101325.)
    0.02446539540458919
    '''
    return R*T/P

Tmin = 0; Tmax = 10e6
Pmin = 0; Pmax = 10e15

ideal_gas_model = TPDependentModel(ideal_gas,
                                   Tmin, Tmax,
                                   Pmin, Pmax,
                                   isaccurate=False)

def TsonopoulosExtendedModel(Tc, Pc, omega, a=0, b=0, species_type='', dipole=0, order=0):
    def Tsonopoulos_extended(T, P):
        return ideal_gas(T, P) + BVirial_Tsonopoulos_extended(T, Tc, Pc, omega, a, b,
                                                              species_type, dipole, order)
    return TPDependentModel(Tsonopoulos_extended, Tmin, Tmax, Pmin, Pmax)
 
def TsonopoulosModel(Tc, Pc, omega):
    def Tsonopoulos(T, P):
        return ideal_gas(T, P) + BVirial_Tsonopoulos(T, Tc, Pc, omega)
    return TPDependentModel(Tsonopoulos, Tmin, Tmax, Pmin, Pmax)

def AbbottModel(Tc, Pc, omega):
    def Abbott(T, P):
        return ideal_gas(T, P) + BVirial_Abbott(T, Tc, Pc, omega)
    return TPDependentModel(Abbott, Tmin, Tmax, Pmin, Pmax)

def PitzerCurlModel(Tc, Pc, omega):
    def Pitzer_Curl(T, P):
        return ideal_gas(T, P) + BVirial_Pitzer_Curl(T, Tc, Pc, omega)
    return TPDependentModel(Pitzer_Curl, Tmin, Tmax, Pmin, Pmax)
    
def CRCVirialModel(a1, a2, a3, a4, a5):
    def CRC_virial(T, P):
        t = 298.15/T - 1.
        return ideal_gas(T, P) + (a1 + a2*t + a3*t**2 + a4*t**3 + a5*t**4)/1e6
    return TPDependentModel(CRC_virial, Tmin, Tmax, Pmin, Pmax)

def EOSModel(eos):
    def EOS(T, P): return eos[0].to_TP(T, P).V_g
    return TPDependentModel(EOS, Tmin, Tmax, Pmin, Pmax)

#PR = 'PR'
CRC_VIRIAL = 'CRC_VIRIAL'
TSONOPOULOS_EXTENDED = 'TSONOPOULOS_EXTENDED'
TSONOPOULOS = 'TSONOPOULOS'
ABBOTT = 'ABBOTT'
PITZER_CURL = 'PITZER_CURL'
IDEAL = 'IDEAL'
NONE = 'NONE'
volume_gas_methods = [COOLPROP, EOS, CRC_VIRIAL, TSONOPOULOS_EXTENDED, TSONOPOULOS,
                      ABBOTT, PITZER_CURL, IDEAL]
'''Holds all methods available for the VolumeGas class, for use in
iterating over them.'''

def VolumeGas(CAS, Tc=None, Pc=None, omega=None, eos=None):
    models = []
    # no point in getting Tmin, Tmax
    if all((Tc, Pc, omega)):
        if eos: models.append(EOSModel(eos))
        models.extend((TsonopoulosExtendedModel(Tc, Pc, omega),
                       TsonopoulosModel(Tc, Pc, omega),
                       AbbottModel(Tc, Pc, omega),
                       PitzerCurlModel(Tc, Pc, omega)))
    if CAS in CRC_virial_data.index:
        a1, a2, a3, a4, a5 = _CRC_virial_data_values[CRC_virial_data.index.get_loc(CAS)].tolist()[1:]
        models.append(CRCVirialModel(a1, a2, a3, a4, a5))
    models.append(ideal_gas_model)
    return TPDependentModelHandler(models)
    # if has_CoolProp and CAS in coolprop_dict:
    #     methods_P.append(COOLPROP)
    #     self.CP_f = coolprop_fluids[self.CASRN]


### Solids

def GoodmanModel(Tt, Vl):
    r'''Calculates solid volume at T using the simple relationship
    by a member of the DIPPR.

    The molar volume of a solid is given by:

    .. math::
        \frac{1}{V_m} = \left( 1.28 - 0.16 \frac{T}{T_t}\right)
        \frac{1}{{Vm}_L(T_t)}

    Parameters
    ----------
    T : float
        Temperature of fluid [K]
    Tt : float
        Triple temperature of fluid [K]
    Vl : float
        Liquid volume, [m^3/mol]

    Returns
    -------
    Vs : float
        Solid volume, [m^3/mol]

    Notes
    -----
    Works to the next solid transition temperature or to approximately 0.3Tt.

    References
    ----------
    .. [1] Goodman, Benjamin T., W. Vincent Wilding, John L. Oscarson, and
       Richard L. Rowley. "A Note on the Relationship between Organic Solid
       Density and Liquid Density at the Triple Point." Journal of Chemical &
       Engineering Data 49, no. 6 (2004): 1512-14. doi:10.1021/je034220e.
    '''
    def Goodman(T): return (1.28 - 0.16*(T/Tt))*Vl
    return TPDependentModel(Goodman, Tmin, Tmax, Pmin, Pmax)


GOODMAN = 'GOODMAN'
CRC_INORG_S = 'CRC_INORG_S'
volume_solid_methods = [GOODMAN, CRC_INORG_S]
'''Holds all methods available for the VolumeSolid class, for use in
iterating over them.'''


def VolumeSolid(CAS):
    r'''Class for dealing with solid molar volume as a function of temperature.
    Consists of one constant value source, and one simple estimator based on
    liquid molar volume.

    Parameters
    ----------
    CASRN : str, optional
        CAS number
    MW : float, optional
        Molecular weight, [g/mol]
    Tt : float, optional
        Triple temperature
    Vml_Tt : float, optional
        Liquid molar volume at the triple point

    Notes
    -----
    A string holding each method's name is assigned to the following variables
    in this module, intended as the most convenient way to refer to a method.
    To iterate over all methods, use the list stored in
    :obj:`volume_solid_methods`.

    **CRC_INORG_S**:
        Constant values in [1]_, for 1872 chemicals.
    **GOODMAN**:
        Simple method using the liquid molar volume. Good up to 0.3*Tt.
        See :obj:`Goodman` for details.

    See Also
    --------
    Goodman

    References
    ----------
    .. [1] Haynes, W.M., Thomas J. Bruno, and David R. Lide. CRC Handbook of
       Chemistry and Physics. [Boca Raton, FL]: CRC press, 2014.
    '''
    models = []
    if CAS in CRC_inorg_s_const_data.index:
        CRC_INORG_S_Vm = float(CRC_inorg_s_const_data.at[CAS, 'Vm'])
        models.append(ConstantTPDependentModel(CRC_INORG_S_Vm, Tmin, Tmax, Pmin, Pmax))
    return TPDependentModelHandler(models)