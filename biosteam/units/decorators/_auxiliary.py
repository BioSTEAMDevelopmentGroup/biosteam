# -*- coding: utf-8 -*-
"""
The auxiliary decorator allows a class to function as an auxiliary unit.
"""
__all__ = ('auxiliary',)


def _setup(self):
    results = (self.baseline_purchase_costs, self.purchase_costs, 
               self.installed_costs, self.F_M, self.F_D, self.F_P,
               self.F_BM)
    for i in results: i.clear()
    for i in self.heat_utilities: i.empty()
    self.heat_utilities.clear()
    self.power_utility.empty()
    
def auxiliary(cls):
    """Make class able to function as an auxiliary unit. The class should 
    compute all results during initialization."""
    if not hasattr(cls, '_setup'): cls._setup = _setup
    return cls