# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
General functional algorithms for batch design.

"""

__all__ = ('size_batch',)

def size_batch(F_vol, tau_reaction, tau_cleaning, N_reactors, V_wf) -> dict:
    """
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
    Units of measure may vary so long as they are consistent. The loading time
    can be considered the cycle time.
        
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
        
        