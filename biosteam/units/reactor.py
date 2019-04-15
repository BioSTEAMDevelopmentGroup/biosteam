# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 22:56:55 2018

@author: yoelr
"""
from .. import Unit, np
from . import Tank


class Reactor(Unit):
    """Abstract class for reactors."""
    

# %% Reactor classes
    
class BatchReactor(Reactor):
    
    _N_heat_utilities = 1
    
    # Exponential scaling
    C_0 = None #: Original Price
    V_0 = None #: Original Volume
    exp = None #: Scaling exponent
    
    @classmethod
    def _solve(cls, 
               v_0: 'Flow rate',
               tau: 'Reaction time',
               tau_0: 'Cleaning and unloading time',
               N_reactors: 'Loading time per volume',
               V_wf: 'Fraction of working volume') -> 'Results [dict]':
        
        V_T = v_0*(tau + tau_0)/(1-1/N_reactors) # Reacting volume
        V_i = V_T/N_reactors
        t_L = V_i/v_0
        t_B = tau + tau_0 + t_L
        V_i /= V_wf
        
        return {'Reactor volume': V_i,
                'Cycle time': t_B, 
                'Loading time': t_L}
    
    
# %%    
