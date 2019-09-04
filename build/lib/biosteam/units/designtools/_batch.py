# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 22:56:55 2018

@author: yoelr
"""

__all__ = ('size_batch',)

def size_batch(v_0, tau, tau_0, N_reactors, V_wf) -> dict:
    """Solve for batch reactor volume, cycle time, and loading time
    
    Parameters
    ----------
    v_0 
        Flow rate
    tau
        Reaction time
    tau_0
        Turnaround time
    N_reactors
        Number of reactors
    V_wf
        Fraction of working volume
    
    Returns
    -------
    dict[str: float]
        * 'Reactor volume'
        * 'Cycle time'
        * 'Loading time'
        
    """
    V_T = v_0*(tau + tau_0)/(1-1/N_reactors) # Reacting volume
    V_i = V_T/N_reactors
    t_L = V_i/v_0
    t_B = tau + tau_0 + t_L
    V_i /= V_wf
    
    return {'Reactor volume': V_i,
            'Cycle time': t_B, 
            'Loading time': t_L}
        
        