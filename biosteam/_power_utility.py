# -*- coding: utf-8 -*-
"""
Created on Sun Nov 11 11:20:42 2018

@author: yoelr
"""
from . import _Q
from ._utils import DisplayUnits

__all__ = ('PowerUtility',)

class PowerUtility:
    """Create an PowerUtility object that can calculates the cost of power.
    
    **__call__()**
    
       Calculate utility requirements given the essential parameters.
        
        **Parameters**
        
            rate: [float] Power requirement (kW)
            
    **Class Parameters**
    
        **price:** ($/kW-hr)
    
    **Examples**
    
        Create a PowerUtility object:
        
        .. code-block:: python
        
           >>> pu = PowerUtility()
           >>> pu
           <PowerUtility: None>
           
        Call object to calculate cost:
            
        .. code-block:: python
        
           >>> pu(rate=500)
           >>> pu
           <PowerUtility: 500 kW, 30 USD/hr>
           
        Results are accessible:
            
        .. code-block:: python
        
           >>> pu.rate, pu.cost
           (500, 30.)
           
        See the object with different units:
            
        .. code-block:: python
        
           >>> pu.show(rate='BTU/s', cost='USD/yr')
           PowerUtility: rate=474 BTU/s, cost=2.63e+05 USD/yr
    
    """
    _units = dict(rate='kW', cost='USD/hr')
    
    #: [DisplayUnits] Units of measure for IPython display
    display_units = DisplayUnits(**_units)
    
    __slots__ = ('rate', 'cost')
    price = 0.0782 #: USD/kWhr
    
    def __init__(self):
        self.rate = 0
        self.cost = 0
    
    def __call__(self, rate:'kW'):
        """Calculate cost and save. 
        
        **Parameters**
        
            rate: [float] Power requirement (kW)
        
        """
        self.rate = rate
        self.cost = self.price * rate
    
    # Representation
    def _info_units(self):
        # Get units of measure
        units = self._units
        rate_units, cost_units = self.display_units
        rate = _Q(self.rate, units['rate']).to(rate_units).magnitude
        cost = _Q(self.cost, units['cost']).to(cost_units).magnitude
        return rate, cost, rate_units, cost_units
        
    def _info(self):
        if self.rate:
            rate, cost, rate_units, cost_units = self._info_units()
            return (f'{type(self).__name__}: \n'
                   +f' rate: {rate:.3g} {rate_units}\n'
                   +f' cost: {cost:.3g} {cost_units}')
        else:
            return (f'{type(self).__name__}: None')

    def _ipython_display_(self):
        print(self._info())

    def __repr__(self):
        if self.rate:
            rate, cost, rate_units, cost_units = self._info_units()
            return (f'<{type(self).__name__}: {self.rate:.3g} {rate_units}, {self.cost:.3g} {cost_units}>')
        else:
            return (f'<{type(self).__name__}: None>')
        
        