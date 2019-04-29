# -*- coding: utf-8 -*-
"""
Created on Sun Nov 11 11:20:42 2018

@author: yoelr
"""
from . import _Q

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
    def _info(self, **show_units):
        # Get units of measure
        su = show_units
        units = self._units
        Rate = su.get('rate') or units['rate']
        Cost = su.get('cost') or units['cost']
        if self.rate:
            rate = _Q(self.rate, units['rate']).to(Rate).magnitude
            cost = _Q(self.cost, units['cost']).to(Cost).magnitude
            return (f'{type(self).__name__}: rate={rate:.3g} {Rate}, cost={cost:.3g} {Cost}')
        else:
            return (f'{type(self).__name__}: None')

    def show(self, **show_units):
        print(self._info(**show_units))

    def __repr__(self):
        units = self._units
        if self.rate:
            Rate = units['rate']
            Cost = units['cost']
            return (f'<{type(self).__name__}: {self.rate:.3g} {Rate}, {self.cost:.3g} {Cost}>')
        else:
            return (f'<{type(self).__name__}: None>')
        
        