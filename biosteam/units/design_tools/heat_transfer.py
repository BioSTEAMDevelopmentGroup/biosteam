# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
General functional algorithms for the design of heat exchangers.

"""
from numpy import log as ln
from numba import njit

__all__ = ('counter_current_heat_exchange',
           'heuristic_overall_heat_transfer_coefficient',
           'heuristic_pressure_drop',
           'heuristic_tubeside_and_shellside_pressure_drops',
           'order_streams',
           'compute_Fahkeri_LMTD_correction_factor',
           'compute_heat_transfer_area',
           'compute_LMTD')

# %% Functional heat exchanger

def heat_exchange_to_condition(s_in, s_out, T=None, phase=None, 
                               H_lim=None, heating=None):
    """
    Set the outlet stream condition and return duty required to achieve 
    said condition.
    """
    H_lim_given = H_lim is not None
    H_in = s_in.H
    if H_lim_given:
        if heating:
            if H_in > H_lim: return 0.
        else:
            if H_in < H_lim: return 0.
    if phase:
        s_out.T = T
        s_out.phase = phase
        if H_lim_given:
            if heating:
                if s_out.H > H_lim: s_out.H = H_lim
            else:
                if s_out.H < H_lim: s_out.H = H_lim
    elif len(s_out.vle_chemicals) == 1 and heating is not None:
        bp = s_in.bubble_point_at_P()
        tol = 1e-3
        dT = T - bp.T
        s_out.T = T
        if dT < -tol:
            s_out.phase = 'l'
            if H_lim_given:
                if heating:
                    if s_out.H > H_lim: s_out.H = H_lim
                else:
                    if s_out.H < H_lim: s_out.vle(H=H_lim , P=s_out.P)
        elif dT > tol:
            s_out.phase = 'g'
            if H_lim_given:
                if heating:
                    if s_out.H > H_lim: s_out.vle(H=H_lim , P=s_out.P)
                else:
                    if s_out.H < H_lim: s_out.H = H_lim
        else:
            s_out.phase = 'g' if heating else 'l'
    else:
        s_out.vle(T=T, P=s_out.P)
        if H_lim_given:
            if heating:
                if s_out.H > H_lim: s_out.vle(H=H_lim , P=s_out.P)
            else:
                if s_out.H < H_lim: s_out.vle(H=H_lim , P=s_out.P)
    
    # Sanity check
    if heating:
        assert s_out.T > s_in.T - 0.01, "failed to solve temperature"
        if s_out.T < s_in.T: s_out.T = s_in.T
    else:
        assert s_out.T <= s_in.T + 0.01, "failed to solve temperature"
        if s_out.T > s_in.T: s_out.T = s_in.T
    
    return s_out.H - H_in
    

def counter_current_heat_exchange(s0_in, s1_in, s0_out, s1_out,
                                  dT, T_lim0=None, T_lim1=None,
                                  phase0=None, phase1=None,
                                  H_lim0=None, H_lim1=None):
    """
    Allow outlet streams to exchange heat until either the given temperature 
    limits or the pinch temperature and return the total heat transfer 
    [Q; in kJ/hr].
    """
    # Counter current heat exchange setup:
    # First find the hot inlet, cold inlet, hot outlet and cold outlet streams
    # along with the maximum temperature approaches for the hotside and the 
    # cold side.
    if s0_in.T > s1_in.T:
        s_hot_in = s0_in
        s_cold_in = s1_in
        s_hot_out = s0_out
        s_cold_out = s1_out
        T_lim_coldside = T_lim0
        T_lim_hotside = T_lim1
        H_lim_coldside = H_lim0
        H_lim_hotside = H_lim1
        phase_coldside = phase0
        phase_hotside = phase1
    else:
        s_cold_in = s0_in
        s_hot_in = s1_in
        s_cold_out = s0_out
        s_hot_out = s1_out
        T_lim_hotside = T_lim0
        T_lim_coldside = T_lim1
        H_lim_hotside = H_lim0
        H_lim_coldside = H_lim1
        phase_hotside = phase0
        phase_coldside = phase1
    
    if (s_hot_in.T - s_cold_in.T) <= dT: return 0. # No heat exchange
    
    T_pinch_coldside = s_cold_in.T + dT
    if T_lim_coldside:
        if T_lim_coldside > s_hot_in.T:
            return 0. # No heat exchange
        else:
            T_lim_coldside = max(T_pinch_coldside, T_lim_coldside) 
    else:
        T_lim_coldside = T_pinch_coldside
    
    T_pinch_hotside = s_hot_in.T - dT
    if T_lim_hotside:
        if T_lim_hotside < s_cold_in.T:
            return 0. # No heat exchange
        else:
            T_lim_hotside = min(T_pinch_hotside, T_lim_hotside) 
    else:
        T_lim_hotside = T_pinch_hotside
       
    # Find which side reaches the pinch first by selecting the side that needs
    # the least heat transfer to reach the pinch.
       
    # Pinch on the cold side
    Q_hot_stream = heat_exchange_to_condition(s_hot_in, s_hot_out, 
                                              T_lim_coldside, phase_coldside,
                                              H_lim_coldside, heating=False)
    
    # Pinch on the hot side
    Q_cold_stream = heat_exchange_to_condition(s_cold_in, s_cold_out, 
                                               T_lim_hotside, phase_hotside,
                                               H_lim_hotside, heating=True)
    
    if Q_hot_stream == Q_cold_stream == 0.:
        s0_out.copy_like(s0_in)
        s1_in.copy_like(s1_out)
        return 0.
    
    if Q_hot_stream > 0 or Q_cold_stream < 0:
        # Sanity check
        if Q_hot_stream / s_hot_in.C < 0.1 or Q_cold_stream / s_cold_in.C > -0.1:
            s0_out.copy_like(s0_in)
            s1_in.copy_like(s1_out)
            return 0.
        raise RuntimeError('inlet stream not in vapor-liquid equilibrium')
    
    if Q_cold_stream < -Q_hot_stream:
        # Pinch on the hot side            
        Q = Q_cold_stream
        if phase_coldside:
            s_hot_out.H = s_hot_in.H - Q
        else:
            s_hot_out.vle(H=s_hot_in.H - Q, P=s_hot_out.P)
    else:
        # Pinch on the cold side
        Q = Q_hot_stream
        if phase_hotside:
            s_cold_out.H = s_cold_in.H - Q
        else:
            s_cold_out.vle(H=s_cold_in.H - Q, P=s_cold_out.P)
    
    return abs(Q)


# %% Heuristics

def heuristic_overall_heat_transfer_coefficient(ci, hi, co, ho):
    """
    Return a heuristic estimate of the overall heat transfer coefficient
    [U; in kW/m^2/K]. Assume `U` is 1.0 kW/m^2/K if heat exchange is between
    a condensing fluid and a vaporizing fluid and 1.0 kW/m^2/K otherwise.
    
    Parameters
    ----------
    ci : Stream
        Cold inlet stream.
    hi : Stream
        Hot inlet stream.
    co : Stream
        Cold outlet stream.
    ho : Stream
        Hot outlet stream.
    
    Returns
    -------
    U : float
        overall heat transfer coefficient [kW/m^2/K].
    
    """
    # TODO: Base U on Table 18.5, Warren D. Seider et. al. Product and Process Design Principles. (2016)
    cip, hip, cop, hop = ci.phase, hi.phase, co.phase, ho.phase
    phases = cip + hip + cop + hop
    if 'g' in phases:
        if ('g' in hip and 'l' in hop) and ('l' in cip and 'g' in cop):
            return 1.0
        else:
            return 0.5
    else:
        return 0.5

def heuristic_pressure_drop(inlet_phase, outlet_phase):
    """
    Return a heuristic estimate of the pressure drop [dP; in psi]. If the fluid 
    changes phase, `dP` is 1.5 psi. If the fluid remains a liquid, `dP` is 5 psi.
    If the fluid remains a gas, `dP` is 3 psi.
    
    Parameters
    ----------
    inlet_phase: str
    outlet_phase : str
        
    Returns
    -------
    dP : float
        Pressure drop [psi].
    
    """
    if ('l' in inlet_phase and 'g' in outlet_phase) or ('g' in inlet_phase and 'l' in outlet_phase):
        # Latent fluid (boiling or condensing)
        dP = 1.5
    elif 'l' in inlet_phase:
        # Sensible liquid
        dP = 5.
    elif 'g' in outlet_phase:
        # Sensible vapor
        dP = 3.
    return dP

def heuristic_tubeside_and_shellside_pressure_drops(ci, hi, co, ho,
                                                    tubeside_iscooling=True):
    """
    Return an estimate of tubeside and shellside pressure drops.
    
    Parameters
    ----------
    ci : Stream
        Cold inlet stream.
    hi : Stream
        Hot inlet stream.
    co : Stream
        Cold outlet stream.
    ho : Stream
        Hot outlet stream.
    tubeside_iscooling : bool
        True of tubeside fluid is cooling.
    
    Returns
    -------  
    dP_tube : float
        Tubeside pressure drop (psi)
    dP_shell : float
        Shellside pressure drop (psi)
    
    """
    dP_c = heuristic_pressure_drop(ci.phase, co.phase)
    dP_h = heuristic_pressure_drop(hi.phase, ho.phase)
    if tubeside_iscooling:
        dP_tube = dP_h
        dP_shell = dP_c
    else:
        dP_tube = dP_c
        dP_shell = dP_h
    return dP_tube, dP_shell

# %% General functions

def order_streams(in_a, in_b, out_a, out_b):
    """
    Return cold and hot inlet and outlet streams.
    
    Parameters
    ----------
    in_a : Stream
        Inlet a.
    in_b : Stream
        Inlet b.
    out_a : Stream
        Outlet a.
    out_b : Stream
        Outlet b.
    
    Returns
    -------
    ci : Stream
        Cold inlet.
    hi : Stream
        Hot inlet.
    co : Stream
        Cold outlet.
    ho : Stream
        Hot outlet.
    
    """
    if in_a.T < in_b.T:
        return in_a, in_b, out_a, out_b
    else:
        return in_b, in_a, out_b, out_a

# %% Computational functions

@njit(cache=True)
def compute_fallback_Fahkeri_LMTD_correction_factor(P, N_shells):
    """Return LMTF correction factor using the fallback equation for
    `compute_Fahkeri_LMTD_correction_factor` when logarithms cannot be computed."""
    # A, J, and K are dummy variables
    A = N_shells - N_shells*P
    W = A/(A + P)
    if 0.999 < W < 1.001:
        Ft = 1
    else:
        J = W/(1. - W)
        K = (J + 2**-0.5)/(J - 2**-0.5)
        if K <= 1:
            Ft = 1
        else:
            Ft = (2**0.5*J)/ln(K)
    return Ft

@njit(cache=True)
def compute_Fahkeri_LMTD_correction_factor(Tci, Thi, Tco, Tho, N_shells):
    r"""
    Return the log-mean temperature difference correction factor `Ft` 
    for a shell-and-tube heat exchanger with one or an even number of tube 
    passes, and a given number of shell passes, with the expression given in 
    [1]_ and also shown in [2]_.
    
    .. math::
        F_t=\frac{S\ln W}{\ln \frac{1+W-S+SW}{1+W+S-SW}}
    
        S = \frac{\sqrt{R^2+1}}{R-1}
        
        W = \left(\frac{1-PR}{1-P}\right)^{1/N}
        
        R = \frac{T_{in}-T_{out}}{t_{out}-t_{in}}
        
        P = \frac{t_{out}-t_{in}}{T_{in}-t_{in}}
        
    If R = 1 and logarithms cannot be evaluated:
        
    .. math::
        W' = \frac{N-NP}{N-NP+P}
        
        F_t = \frac{\sqrt{2}\frac{1-W'}{W'}}{\ln\frac{\frac{W'}{1-W'}+\frac{1}
        {\sqrt{2}}}{\frac{W'}{1-W'}-\frac{1}{\sqrt{2}}}}
        
    Parameters
    ----------
    Tci : float
        Inlet temperature of cold fluid, [K]
    Thi : float
        Inlet temperature of hot fluid, [K]
    Tco : float
        Outlet temperature of cold fluid, [K]        
    Tho : float
        Outlet temperature of hot fluid, [K]
    shells : int, optional
        Number of shell-side passes, [-]
    
    Returns
    -------
    Ft : float
        Log-mean temperature difference correction factor, [-]
    
    Notes
    -----
    This expression is symmetric - the same result is calculated if the cold
    side values are swapped with the hot side values. It also does not 
    depend on the units of the temperature given.
    
    Examples
    --------
    compute_Fahkeri_LMTD_correction_factor(Tci=15, Tco=85, Thi=130, Tho=110, N_shells=1)
    0.9438358829645933
    
    References
    ----------
    .. [1] Fakheri, Ahmad. "A General Expression for the Determination of the 
       Log Mean Temperature Correction Factor for Shell and Tube Heat 
       Exchangers." Journal of Heat Transfer 125, no. 3 (May 20, 2003): 527-30.
       doi:10.1115/1.1571078.
    .. [2] Hall, Stephen. Rules of Thumb for Chemical Engineers, Fifth Edition.
       Oxford; Waltham, MA: Butterworth-Heinemann, 2012.
    
    """
    if (Tco - Tci) < 0.01:
        R = 1
    else:
        R = (Thi - Tho)/(Tco - Tci)
    P = (Tco - Tci)/(Thi - Tci)
    if 0.999 < R < 1.001:
        Ft = compute_fallback_Fahkeri_LMTD_correction_factor(P, N_shells)
    else:
        W = ((1. - P*R)/(1. - P))**(1./N_shells)
        S = (R*R + 1.)**0.5/(R - 1.)
        K = (1. + W - S + S*W)/(1. + W + S - S*W)
        if K <= 0.001 or 0.999 < K < 1.001:
            Ft = compute_fallback_Fahkeri_LMTD_correction_factor(P, N_shells)
        else:
            Ft = S*ln(W)/ln(K)
    if Ft > 1.0:
        Ft = 1.0
    elif Ft < 0.5:
        # Bad design, probably a heat exchanger network operating
        # too close to the pinch. Fahkeri may not be valid, so give
        # a conservative estimate of the correction factor.
        Ft = 0.5
    return Ft

@njit(cache=True)
def compute_heat_transfer_area(LMTD, U, Q, ft):
    """
    Return required heat transfer area by LMTD correction factor method.
    
    Parameters
    ----------
    LMTD : float
        Log mean temperature difference
    U : float
        Heat transfer coefficient
    Q : float
        Duty
    
    """
    return Q/(U*LMTD*ft)   

@njit(cache=True)
def compute_LMTD(Thi, Tho, Tci, Tco, counterflow=True):
    r'''
    Return the log-mean temperature difference of an ideal counterflow
    or co-current heat exchanger.

    .. math::
        \Delta T_{LMTD}=\frac{\Delta T_1-\Delta T_2}{\ln(\Delta T_1/\Delta T_2)}

        \text{For countercurrent:      } \\
        \Delta T_1=T_{h,i}-T_{c,o}\\
        \Delta T_2=T_{h,o}-T_{c,i}

        \text{Parallel Flow Only:} \\
        {\Delta T_1=T_{h,i}-T_{c,i}}\\
        {\Delta T_2=T_{h,o}-T_{c,o}}

    Parameters
    ----------
    Thi : float
        Inlet temperature of hot fluid [K].
    Tho : float
        Outlet temperature of hot fluid [K].
    Tci : float
        Inlet temperature of cold fluid [K].
    Tco : float
        Outlet temperature of cold fluid [K].
    counterflow : bool, optional
        Whether the exchanger is counterflow or co-current.

    Returns
    -------
    LMTD : float
        Log-mean temperature difference [K]

    Notes
    -----
    Any consistent set of units produces a consistent output.

    Examples
    --------
    >>> compute_LMTD(100., 60., 30., 40.2)
    43.200409294131525
    >>> compute_LMTD(100., 60., 30., 40.2, counterflow=False)
    39.75251118049003
    
    '''
    if counterflow:
        dTF1 = Thi-Tco
        dTF2 = Tho-Tci
    else:
        dTF1 = Thi-Tci
        dTF2 = Tho-Tco
    dTF21 = dTF2 - dTF1
    if abs(dTF21) < 1e-8:
        LMTD = dTF1
    else:
        LMTD = dTF21/ln(dTF2 / dTF1)
    return LMTD