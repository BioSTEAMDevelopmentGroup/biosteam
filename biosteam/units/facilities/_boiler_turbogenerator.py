# -*- coding: utf-8 -*-
"""
Created on Sun Feb 17 16:36:47 2019

@author: yoelr
"""
from ... import Unit
from ..decorators import cost

__all__ = ('BoilerTurbogenerator',)

#: TODO add reference of NREL

@cost('Heat', CE=567.3, cost=36e6, S=104e7, n=0.95)
class BoilerTurbogenerator(Unit):
    """Create a Boiler_TurboGenerator object that will calculate electricity generation from burning the feed. It also takes into account how much steam is being produced, and the required cooling utility of the turbo generator. No emissions or mass balances are taken into account.
    
    **Parameters**
    
        **boiler_efficiency:** [float] Fraction of heat transfered to steam.
    
        **turbo_generator_efficiency:** [float] Fraction of steam heat converted to electricity.

        **steam:** [Stream] Steam product from the boiler.
    
    **ins**
    
        [0] Feed that will be burned.
    
    **outs**
    
        [0] Emission (burned feed)
        
    
    """
    _N_heat_utilities = 2
    _has_power_utility = True
    _units = {'Heat': 'kJ/hr'}
    _kwargs = {'boiler_efficiency': 0.80,
               'turbo_generator_efficiency': 0.85,
               'steam': None}
    
    def _design(self):
        B_eff, TG_eff, steam = self._kwargs.values()
        feed = self.ins[0]
        emission, = self.outs
        hu_cooling, hu_steam = self._heat_utilities
        lps = hu_steam.heating_agents['Low pressure steam']
        H_steam = steam.H
        self._Design['Heat'] = feed_Hc = feed.Hc
        
        # This is simply for the mass balance (no special purpose)
        emission.mol[:] = feed.mol # TODO: In reality, this should be CO2
        emission.T = steam.T
        emission.phase = 'g'
        
        feed_massnet = feed.massnet
        H_content = feed_Hc*B_eff
        H_electricity = H_content - H_steam - feed_massnet*1130 # 50 percent of bagasse is water, so remove latent heat of vaporization
        if H_electricity < 0:
            H_steam = H_content
            cooling = electricity = 0
        else:
            electricity = H_electricity * TG_eff
            cooling = electricity - H_electricity
        hu_cooling(cooling, steam.T)
        hu_steam.duty = -H_steam
        hu_steam.ID = 'Low pressure steam'
        hu_steam.flow = molnet = -steam.molnet
        hu_steam.cost = molnet*lps.price_kmol + H_steam*lps.price_kJ
        self._power_utility(-electricity/3600)
        

# Simulation of ethanol production from sugarcane
# in Brazil: economic study of an autonomous
# distillery
# Marina O.S. Diasa,b, Marcelo P. Cunhaa, Charles D.F. Jesusa, Mirna I.G.
# Scandiffioa, Carlos E.V. Rossella,b, Rubens Maciel Filhob, Antonio Bonomia
# a CTBE – Bioethanol Science and Technology National Laboratory, PO Box 6170 –
# CEP 13083-970, Campinas – SP, Brazil, marina.dias@bioetanol.org.br
# b School of Chemical Engineering, University of Campinas, PO Box 6066 – CEP
# 13083-970, Campinas – SP, Brazil
# LHV = 7565 # Bagasse lower heating value (kJ/kg)