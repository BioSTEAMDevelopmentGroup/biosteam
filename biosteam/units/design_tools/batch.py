# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
General functional algorithms for batch design.

"""
from math import ceil
from flexsolve import wegstein

__all__ = ('size_batch',)

def size_batch(F_vol, tau_reaction, tau_cleaning, V_wf, 
               V_max=None, N_reactors=None, loading_time=None) -> dict:
    r"""
    Solve for batch reactor volume, cycle time, and loading time.
    
    Parameters
    ----------
    F_vol : float
        Volumetric flow rate.
    tau_reaction : float
        Reaction time.
    tau_cleaning : float
        Cleaning in place time.
    V_wf : float
        Fraction of working volume.
    V_max : int
        Maximum volume of a reactor.
    N_reactors : int
        Number of reactors.
    loading_time :
        Loading time of batch reactor. If not given, it will assume each vessel is constantly
        being filled.
        
    Returns
    -------
    dict
        * 'Reactor volume': float
        * 'Batch time': float
        * 'Loading time': float
       
    Notes
    -----
    By assuming no downtime, the total volume of all reactors is:
        
    .. math::
        
        V_T = F_{vol}(\tau_{reaction} + \tau_{cleaning} + \tau_{loading})
        
    where :math:`V_T` is the total volume required, :math:`F_{vol}` is the 
    volumetric flow rate of the feed, :math:`\tau_{reaction}` is the 
    reaction time, :math:`\tau_{cleaning}` is the cleaning and unloading time, 
    and :math:`\tau_{loading}` is the time required to load a vessel. This
    equation makes the conservative assumption that no reaction takes place 
    when the tank is being filled.
    
    The number of reactors is:
    
    .. math::
        
        N_{reactors} = ceil(\frac{V_T}{V_{max} \cdot V_{wf}})
    
    where :math:`N_{reactors}` is the number of reactor vessels, 
    :math:`V_{max}` is the maximum reactor volume, and 
    :math:`V_{wf}` is the fraction of working volume in a reactor.
    
    The working volume of an individual reactor is:
    
    .. math::
    
        V_{i,working} = \frac{V_T}{N_{reactors}}
    
    If there is no upstream storage and a vessel is constantly being filled,
    the time required to load a reactor is:
        
    .. math::

        \tau_{loading} = \frac{V_{i,working}}{F_{vol}}
    
    Note that the actual volume of a reactor is:
        
    .. math::

        V_i = \frac{V_{i,working}}{f}
    
    Plugging in and solving for the total volume, :math:`V_{T}`, we have: 
    
    .. math::

        V_T = F_{vol}\frac{\tau_{reaction} + \tau_{cleaning}}{1 - \frac{1}{N_{reactors}}}
    
    If the number of reactors is specified but not the loading time, :math:`V_T` is first calculated, 
    then :math:`V_{i, working}`, :math:`\tau_{loading}`, and :math:`V_i`. 
    
    If neither the number of reactors nor the loading time are specified, we solve the equations iteratively
    until :math:`\tau_{loading}` converges.
    
    If the loading time is given but not the number of reactors, 
    then the equation for :math:`tau_{loading}` does not apply and 
    we first compute :math:`N_reactors` given :math:`V_{max}.
    
    Units of measure may vary so long as they are consistent. The loading time
    can be considered the cycle time in this scenario.
    
    Examples
    --------
    Size batch given a maximum reactor volume of 1,000 m3 and zero loading time:
    
    >>> from biosteam.units.design_tools import size_batch
    >>> F_vol = 1e3; tau_reaction = 30; tau_cleaning = 3; V_wf = 0.95
    >>> size_batch(
    ...     F_vol, tau_reaction, tau_cleaning, V_wf, 
    ...     V_max=1e3, loading_time=0
    ... )
    {'Reactor volume': 992.48,
     'Batch time': 33,
     'Loading time': 0,
     'Number of reactors': 35}
    
    Size batch given 35 reactors and zero loading time:
    
    >>> size_batch(
    ...     F_vol, tau_reaction, tau_cleaning, V_wf, 
    ...     N_reactors=35, loading_time=0
    ... )
    {'Reactor volume': 992.48,
     'Batch time': 33,
     'Loading time': 0,
     'Number of reactors': 35}
    
    Size batch given a maximum reactor volume of 1,000 m3 and assume
    the constant loading:
    
    >>> size_batch(
    ...     F_vol, tau_reaction, tau_cleaning, V_wf, 
    ...     V_max=1000,
    ... )    
    {'Reactor volume': 992.4812030075188,
     'Batch time': 33.94285714285714,
     'Loading time': 0.9428571428571428,
     'Number of reactors': 36}
    
    Size batch given 36 reactors and assume
    the constant loading:
    
    >>> size_batch(
    ...     F_vol, tau_reaction, tau_cleaning, V_wf, 
    ...     N_reactors=36
    ... )    
    {'Reactor volume': 992.4812030075188,
     'Batch time': 33.94285714285714,
     'Loading time': 0.9428571428571428,
     'Number of reactors': 36}
    
    
    """
    if sum([i is None for i in [N_reactors, V_max]]) != 1:
        raise ValueError('must pass either `N_reactors` or `V_max`')
    if loading_time is None:
        if N_reactors is None:
            # Solve iteratively
            def f(tau_loading):
                N_reactors = F_vol * (tau_reaction + tau_cleaning + tau_loading) / (V_max * V_wf)
                V_T = F_vol * (tau_reaction + tau_cleaning) / (1 - 1 / N_reactors)
                V_i = V_T/N_reactors
                tau_loading = V_i/F_vol
                return tau_loading
            
            tau_loading = wegstein(f, 0.5, checkconvergence=False, checkiter=False)
            N_reactors = ceil(F_vol * (tau_reaction + tau_cleaning + tau_loading) / (V_max * V_wf))
            V_T = F_vol * (tau_reaction + tau_cleaning + tau_loading)
            V_i = V_T / N_reactors
        else:
            # Total volume of all reactors, assuming no downtime
            V_T = F_vol * (tau_reaction + tau_cleaning) / (1 - 1 / N_reactors)
            
            # Volume of an individual reactor
            V_i = V_T / N_reactors
            
            # Time required to load a reactor
            tau_loading = V_i / F_vol
    else:
        tau_loading = loading_time
        
        if N_reactors is None:
            # Number of reactors
            N_reactors = ceil(F_vol * (tau_reaction + tau_cleaning + tau_loading) / (V_max * V_wf))

        # Total volume of all reactors
        V_T = F_vol * (tau_reaction + tau_cleaning + tau_loading)
        
        # Volume of an individual reactor
        V_i = V_T / N_reactors
    
    # Total batch time
    tau_batch = tau_reaction + tau_cleaning + tau_loading
    
    # Account for excess volume
    V_i /= V_wf 
    
    return {'Reactor volume': V_i,
            'Batch time': tau_batch,
            'Loading time': tau_loading,
            'Number of reactors': N_reactors}
        
        