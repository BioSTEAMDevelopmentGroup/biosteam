# -*- coding: utf-8 -*-
"""
Created on Sun Nov 11 11:20:42 2018

@author: yoelr
"""
from thermosteam.base.units_of_measure import DisplayUnits, convert

__all__ = ('PowerUtility',)

class PowerUtility:
    """Create an PowerUtility object that, when called, calculates the cost of power (kW) and saves the power and cost.
            
    Examples
    --------
    Create a PowerUtility object:
    
    .. code-block:: python
    
       >>> pu = PowerUtility()
       >>> pu
       <PowerUtility: 0 kW, 0 USD/hr>
       
    Call object to calculate cost:
        
    .. code-block:: python
    
       >>> pu(rate=500)
       >>> pu
       <PowerUtility: 500 kW, 39.1 USD/hr>
       
    Results are accessible:
        
    .. code-block:: python
    
       >>> pu.rate, pu.cost
       (500, 39.1)
    
    """
    __slots__ = ('rate', 'cost')
    
    #: [DisplayUnits] Units of measure for IPython display
    display_units = DisplayUnits(rate='kW', cost='USD/hr')
    
    #: [float] USD/kWhr
    price = 0.0782
    
    def __init__(self):
        #: Power requirement (kW)
        self.rate = 0
        
        #: Cost (USD/hr)
        self.cost = 0
    
    def __bool__(self):
        return bool(self.rate)
    
    def __call__(self, rate):
        """Calculate cost and save results. 
        
        Parameters
        ----------
        rate : float
               Power requirement (kW)
        
        """
        self.rate = rate
        self.cost = self.price * rate
    
    def show(self, rate=None, cost=None):
        # Get units of measure
        display_units = self.display_units
        rate_units = rate or display_units.rate
        cost_units = cost or display_units.cost
        rate = convert(self.rate, 'kW', rate_units)
        cost = convert(self.cost, 'USD/hr', cost_units)
        return (f'<{type(self).__name__}: {rate:.3g} {rate_units}, {cost:.3g} {cost_units}>')
    _ipython_display = show    
    
    def __repr__(self):
        return (f'<{type(self).__name__}: {self.rate:.3g} kW, {self.cost:.3g} USD/hr>')    