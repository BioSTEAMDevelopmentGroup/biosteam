# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
General functional algorithms for batch design.

"""

__all__ = ('size_batch',)

def size_batch(F_vol, tau_reaction, tau_cleaning, N_reactors, V_wf) -> dict:
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
    N_reactors : int
        Number of reactors.
    V_wf : float
        Fraction of working volume.
    
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
        
    where :math:`V_T` is the total volume of all reactors, :math:`F_{vol}` is the 
    volumetric flow rate of the feed, :math:`\tau_{reaction}` is the 
    reaction time, :math:`\tau_{cleaning}` is the cleaning and unloading time, 
    and :math:`\tau_{loading}` is the time required to load a vessel. This
    equation makes the conservative assumption that no reaction takes place 
    when the tank is being filled.
    
    The working volume of an individual reactor is:
    
    .. math::
    
        V_{i,working} = \frac{V_T}{N_{reactors}}
    
    where :math:`N_{reactors}` is the number of reactor vessels.
    
    The time required to load a reactor (assuming no downtime) is:
        
    .. math::

        \tau_{loading} = \frac{V_{i,working}}{F_{vol}}
    
    Note that the the actual volume of a reactor is:
        
    .. math::

        V_i = \frac{V_{i,working}}{f}
        
    where f is the fraction of working volume in a reactor.
    
    Plugging in and solving for the total volume, :math:`V_{T}`, we have: 
    
    .. math::

        V_T = F_{vol}\frac{\tau_{reaction} + \tau_{cleaning}}{1 - \frac{1}{N_{reactors}}}
    
    Using this equation, :math:`V_T` is first calculated, then :math:`V_{i, working}`, 
    :math:`\tau_{loading}`, and :math:`V_i`.
    
    Units of measure may vary so long as they are consistent. The loading time
    can be considered the cycle time in this scenario.
    
    """
    # Total volume of all reactors, assuming no downtime
    V_T = F_vol * (tau_reaction + tau_cleaning) / (1 - 1 / N_reactors)
    
    # Volume of an individual reactor
    V_i = V_T/N_reactors
    
    # Time required to load a reactor
    tau_loading = V_i/F_vol
    
    # Total batch time
    tau_batch = tau_reaction + tau_cleaning + tau_loading
    
    # Account for excess volume
    V_i /= V_wf 
    
    return {'Reactor volume': V_i,
            'Batch time': tau_batch,
            'Loading time': tau_loading}
        
        