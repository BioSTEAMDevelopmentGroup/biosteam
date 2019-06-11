# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 00:55:03 2019

Mechanical equipment

@author: yoelr
"""

from ... import Unit, Mixer, Stream
from ..metaclasses import splitter
from .. import decorators
from .designtools import CEPCI_by_year as CE

# %% Decorators 

_massflow_units = {'Flow rate': 'kg/hr'}
def _design_kg_hr(self):
    self._results['Design']['Flow rate'] = self._ins[0].massnet

def cost_kg_hr(*, cost, exp, S, kW=0, CE=CE[2009], N=1, basis='Flow rate'):
    def decorator(cls):
        cls._design = _design_kg_hr
        cls._units = _massflow_units
        return decorators.cost(basis, cost=cost, exp=exp,
                               CE=CE, S=S, kW=kW, N=N)(cls)
    return decorator

# %% Units

@cost_kg_hr(cost=13329690, exp=0.6, S=94697, kW=511.321)
class FeedStockHandling(Unit): BM = 1.7
    
@cost_kg_hr(cost=6000, exp=0.5, S=136260)
class SulfuricAcidMixer(Unit): BM = 1.0
    
@cost_kg_hr(cost=19812400, exp=0.6, S=83333, kW=5290)
class PretreatmentReactorSystem(Unit):
    BM = 1.5
    def _init(self):
        self._water_mass = Stream.indices('Water')
    def _design(self):
        feed = self._ins[0]
        self._results['Design']['Flow rate'] = feed.massnet - feed.mass[self._water_index]
        
    