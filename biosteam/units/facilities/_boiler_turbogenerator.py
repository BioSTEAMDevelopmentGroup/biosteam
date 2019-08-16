# -*- coding: utf-8 -*-
"""
Created on Sun Feb 17 16:36:47 2019

@author: yoelr
"""
from ... import Stream, Species, HeatUtility, Unit
from . import Facility
from ..decorators import cost

__all__ = ('BoilerTurbogenerator',)

lps = HeatUtility.heating_agents['Low pressure steam']

#: TODO add reference of NREL

@cost('Work', 'Turbogenerator',
      CE=551, S=42200, kW=0, cost=9500e3, n=0.60, BM=1.8)
@cost('Flow rate', 'Hot process water softener system', 
      CE=551, cost=78e3, S=235803, n=0.6, BM=1.8)
@cost('Flow rate', 'Amine addition pkg', 
      CE=551, cost=40e3, S=235803, n=0.6, BM=1.8)
@cost('Flow rate', 'Deaerator',
      CE=551, cost=305e3, S=235803, n=0.6, BM=3.0)
@cost('Flow rate', 'Boiler',
      CE=551, cost=28550e3, kW=2752, S=238686, n=0.6, BM=1.8)
class BoilerTurbogenerator(Facility):
    """Create a Boiler_TurboGenerator object that will calculate electricity generation from burning the feed. It also takes into account how much steam is being produced, and the required cooling utility of the turbo generator. No emissions or mass balances are taken into account.
    
    **Parameters**
    
        **boiler_efficiency:** [float] Fraction of heat transfered to steam.
    
        **turbo_generator_efficiency:** [float] Fraction of steam heat converted to electricity.
    
    **ins**
    
        [0] Feed that will be burned.
    
    **outs**
    
        [0] Emission (burned feed)
        
    
    """
    _N_heat_utilities = 2
    _has_power_utility = True
    _units = {'Flow rate': 'kg/hr',
              'Work': 'kW'}
    
    def __init__(self, ID='', ins=None, outs=(), *,
                 boiler_efficiency=0.80,
                 turbogenerator_efficiency=0.85):
        Unit.__init__(self, ID, ins, outs)
        self.boiler_efficiency = boiler_efficiency
        self.turbogenerator_efficiency = turbogenerator_efficiency
        self.steam_utilities = set()
        self.steam_demand = Stream(None, species=Species('Water'), P=lps.P, T=lps.T, phase='g')
    
    def _run(self): pass
    
    def _design(self):
        B_eff = self.boiler_efficiency
        TG_eff = self.turbogenerator_efficiency
        steam = self.steam_demand
        steam_utilities = self.steam_utilities
        if not steam_utilities:
            for u in self.system.units:
                if u is self: continue
                for hu in u._heat_utilities:
                    if hu.ID == 'Low pressure steam':
                        steam_utilities.add(hu)
        steam._mol[0] = sum([i.flow for i in steam_utilities])
        feed = self._ins[0]
        emission = self._outs[0]
        hu_cooling, hu_steam = self._heat_utilities
        H_steam = steam.H
        feed_Hc = feed.Hc
        
        # This is simply for the mass balance (no special purpose)
        emission.mol[:] = feed.mol # TODO: In reality, this should be CO2
        
        # A percent of bagasse is water, so remove latent heat of vaporization
        feed_massnet = feed.massnet
        moisture_content = feed._mass[feed.index('7732-18-5')]/feed_massnet
        H_content = feed_Hc*B_eff - feed_massnet*2260*moisture_content
        
        # Heat available for the turbogenerator
        H_electricity = H_content - H_steam 
        
        #: [float] Total steam produced by the boiler (kmol/hr)
        self.total_steam = H_content/70000 # Superheat steam with 70000 kJ/kmol
        Design = self._Design
        Design['Flow rate'] = self.total_steam * 18.01528
        
        if H_electricity < 0:
            H_steam = H_content
            cooling = electricity = 0
        else:
            electricity = H_electricity * TG_eff
            cooling = electricity - H_electricity
        hu_cooling(cooling, steam.T)
        hu_steam.ID = 'Low pressure steam'
        hu_steam.cost = -sum([i.cost for i in steam_utilities])
        Design['Work'] = electricity/3600        

    def _end_decorated_cost_(self):
        self._power_utility(self._power_utility.rate - self._Design['Work'])

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