# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 22:56:55 2018

@author: yoelr
"""
from biosteam import Unit, np
from . import Tank


class Reactor(Unit):
    """Abstract class for reactors."""
    _N_heat_util = 1
    _tau = 1
    
    # Residence time property
    tau = Tank.tau 
    
    def _calc_utility(self):
        self.heat_utilities[0](self.Hnet, self.outs[0].T)

# %% Reactor classes
    
class BatchReactor(Reactor):
    
    C_0 = None #: Original Price
    V_0 = None #: Original Volume
    exp = None #: Scaling exponent
    
    @classmethod
    def _solve(cls, 
               v_0: 'Flow rate',
               tau: 'Reaction time',
               tau_0: 'Cleaning and unloading time',
               f: 'Loading time per volume',
               V_wf: 'Fraction of working volume') -> 'Results [dict]':
        V_T = cls._calc_TotalVolume(v_0, tau, tau_0, f)
        t_L = f*V_T
        N = cls._calc_Nreactors(v_0, f)
        t_B = cls._calc_CycleTime(t_L, tau, tau_0)
        T_D = N*t_L - t_B
        V_T *= (t_B + T_D)/t_B # Actual working volume required 
        V_i = V_T/N
        V_i /= V_wf
        
        return {'Number of reactors': N,
                'Reactor volume': V_i,
                'Cycle time': t_B, 
                'Loading time': t_L,
                'Total dead time': T_D}
        
    @staticmethod
    def _calc_TotalVolume(v_0: 'Flow rate',
                          tau: 'Reaction time',
                          tau_0: 'Cleaning and unloading time',
                          f: 'Loading time per volume') -> 'V_T':
        return (tau + tau_0)/(1/v_0 - f)
        
    @staticmethod
    def _calc_CycleTime(t_L: 'Loading time',
                        tau: 'Reaction time',
                        tau_0: 'Cleaning and unloading time') -> 't_B':
        return t_L + tau + tau_0
    
    @staticmethod
    def _calc_Nreactors(v_0: 'Flow rate',
                        f: 'Loading time per volume') -> 'V_T':
        return np.ceil(1/(f*v_0))
    
    
# %%    
