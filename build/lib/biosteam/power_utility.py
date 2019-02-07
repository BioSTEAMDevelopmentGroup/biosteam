# -*- coding: utf-8 -*-
"""
Created on Sun Nov 11 11:20:42 2018

@author: yoelr
"""
from bookkeep import SmartBook, UnitManager, Q_

class PowerUtility:
    """Create an PowerUtility object that can calculates the cost of power.
    
    **__call__()**
    
        Return dictionary of utility requirements given the essential parameters.
        
        **Parameters**
        
            power: [float] Power requirement (kW)
            
        **Returns**
            * 'Power': Unit duty requirement (kW)
            * 'Cost': Cost of utility (USD/hr)
    
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
           {'Power': 500 (kW),
            'Cost': 30 (USD/hr)}
           
        PowerUtility objects retain results:
            
        .. code-block:: python
        
           >>> pu
           <PowerUtility: 500 kW, 30 USD/hr>
           >>> pu.results
           {'Power': 500 (kW),
            'Cost': 30 (USD/hr)}
           
        See the object with different units:
            
        .. code-block:: python
        
           >>> pu.show(Power='BTU/s', Cost='USD/yr')
           PowerUtility: Power=474 BTU/s, Cost=2.63e+05 USD/yr
    
    """
    _units = UnitManager([], Power='kW', Cost='USD/hr')
    __slots__ = ('results',)
    price = 0.0782
    
    def __init__(self, source):
        self.results = SmartBook(self._units, source=source)
    
    def __call__(self, power:'kW'):
        """Return dictionary of utility requirements given the essential parameters.
        
        **Parameters**
        
            power: [float] Power requirement (kW)
            
        **Returns**
            * 'Power': Unit duty requirement (kW)
            * 'Cost': Cost of utility (USD/hr)
        
        """
        results = self.results
        results['Power'] = power
        results['Cost'] = self.price * power
        return results
    
    # Representation
    def _info(self, **show_units):
        # Get units of measure
        su = show_units
        r = self.results
        Power = su.get('Power') or su.get('power') or r.units['Power']
        Cost = su.get('Cost') or su.get('cost') or r.units['Cost']
        if r:
            power = Q_(r['Power'], r.units['Power']).to(Power).magnitude
            cost = Q_(r['Cost'], r.units['Cost']).to(Cost).magnitude
            return (f'{type(self).__name__}: Power={power:.3g} {Power}, Cost={cost:.3g} {Cost}')
        else:
            return (f'{type(self).__name__}: None')

    def show(self, **show_units):
        print(self._info(**show_units))

    def __repr__(self):
        r = self.results
        if r:
            power = r['Power']
            Power = r.units['Power']
            cost = r['Cost']
            Cost = r.units['Cost']
            return (f'<{type(self).__name__}: {power:.3g} {Power}, {cost:.3g} {Cost}>')
        else:
            return (f'<{type(self).__name__}: None>')