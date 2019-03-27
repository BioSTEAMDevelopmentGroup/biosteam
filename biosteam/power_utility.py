# -*- coding: utf-8 -*-
"""
Created on Sun Nov 11 11:20:42 2018

@author: yoelr
"""
from biosteam import Q_

__all__ = ('PowerUtility',)

class PowerUtility:
    """Create an PowerUtility object that can calculates the cost of power.
    
    **__call__()**
    
       Calculate utility requirements given the essential parameters.
        
        **Parameters**
        
            power: [float] Power requirement (kW)
            
    **Class Parameters**
    
        **price:** ($/kW-hr)
    
    **Examples**
    
        Create a PowerUtility object:
        
        .. code-block:: python
        
           >>> from biosteam import PowerUtility
           >>> pu = PowerUtility()
           >>> pu
           <PowerUtility: None>
           
        Call object to calculate cost:
            
        .. code-block:: python
        
           >>> pu(power=500)
           >>> pu
           <PowerUtility: 500 kW, 30 USD/hr>
           
        Results are accessible:
            
        .. code-block:: python
        
           >>> pu.power, pu.cost
           (500, 30.)
           
        See the object with different units:
            
        .. code-block:: python
        
           >>> pu.show(power='BTU/s', cost='USD/yr')
           PowerUtility: power=474 BTU/s, post=2.63e+05 USD/yr
    
    """
    _units = dict(power='kW', cost='USD/hr')
    __slots__ = ('power', 'cost')
    price = 0.0782 #: USD/kWhr
    
    def __init__(self):
        self.power = 0
        self.cost = 0
    
    def __call__(self, power:'kW'):
        """Return dictionary of utility requirements given the essential parameters.
        
        **Parameters**
        
            power: [float] Power requirement (kW)
        
        """
        self.power = power
        self.cost = self.price * power
    
    # Representation
    def _info(self, **show_units):
        # Get units of measure
        su = show_units
        units = self._units
        Power = su.get('power') or units['power']
        Cost = su.get('cost') or units['cost']
        if self.power:
            power = Q_(self.power, units['power']).to(Power).magnitude
            cost = Q_(self.cost, units['cost']).to(Cost).magnitude
            return (f'{type(self).__name__}: power={power:.3g} {Power}, cost={cost:.3g} {Cost}')
        else:
            return (f'{type(self).__name__}: None')

    def show(self, **show_units):
        print(self._info(**show_units))

    def __repr__(self):
        units = self._units
        if self.power:
            Power = units['power']
            Cost = units['cost']
            return (f'<{type(self).__name__}: {self.power:.3g} {Power}, {self.cost:.3g} {Cost}>')
        else:
            return (f'<{type(self).__name__}: None>')
        
        