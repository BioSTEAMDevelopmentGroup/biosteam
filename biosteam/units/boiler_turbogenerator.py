# -*- coding: utf-8 -*-
"""
Created on Sun Feb 17 16:36:47 2019

@author: yoelr
"""
from biosteam import Unit

#: TODO add reference of NREL

class BoilerTurbogenerator(Unit):
    """Create a Boiler_TurboGenerator object that will calculate electricity generation from burning the feed. It also takes into account how much steam is being produced, and the required cooling utility of the turbo generator. No emissions or mass balances are taken into account.
    
    **Parameters**
    
        **boiler_efficiency:** [float] Fraction of heat transfered to steam.
        **turbo_generator_efficiency:** [float] Fraction of heat converted to electricity.
    
    **ins**
    
        [0] Feed that will be burned.
    
    **outs**
    
        [0] Steam extracted from the boiler
        
    .. Note::
        
        The steam (outs[0]) must be specified with the right flow rate, temperature, pressure, and phase before simulation.
    
    """
    _N_heat_utilities = 1
    _has_power_utility = True
    
    #: Base CEPCI 
    CEPCI_0 = 567
    
    #: Base flow rate (kg/hr)
    F_0 = 143445.52660170314
    
    #: Base capital cost price (USD)    
    C_0 = 36e6
    
    #: Scaling exponent for costing
    exp = 0.6
    
    kwargs = {'boiler_efficiency': 0.80,
              'turbo_generator_efficiency': 0.85}
    
    def _cost(self):
        B_eff, TG_eff = self.kwargs.values()
        feed = self.ins[0]
        steam, emission = self.outs
        
        H_steam = steam.H
        H_content = feed.Hc * B_eff
        
        emission.mol[:] = feed.mol
        emission.phase = 'g'
        H_electicity = H_content - H_steam - emission.H
        electricity = H_electicity * TG_eff
        cooling = electricity - H_electicity
        hu = self.heat_utilities[0]
        hu(cooling, steam.T)
        self.power_utility(-electricity/3600)
        flow = self.ins[0].massnet
        r = self.results
        r['Cost']['Purchase price'] = self.CEPCI/self.CEPCI_0 * self.C_0*(flow/self.F_0)**self.exp
    

# Water_index = steam.get_index('7732-18-5')
# steam.mass[Water_index] = S_ext
# steam.P = S_P
# steam.T = steam.bubble_point()[0]