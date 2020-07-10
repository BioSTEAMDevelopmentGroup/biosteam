# -*- coding: utf-8 -*-
<<<<<<< HEAD
"""
Created on Sun Nov 11 11:20:42 2018

@author: yoelr
"""
from thermosteam.base.units_of_measure import DisplayUnits, convert
=======
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from thermosteam.units_of_measure import DisplayUnits, convert
>>>>>>> cd2c5013aaf9b5bc94bb764b52fd37db183472f1

__all__ = ('PowerUtility',)

class PowerUtility:
    """
<<<<<<< HEAD
    Create an PowerUtility object that, when called, calculates the cost of power (kW) and saves the power and cost.
=======
    Create an PowerUtility object that, when called, calculates the cost of 
    power [kW] and saves the rate and cost.
>>>>>>> cd2c5013aaf9b5bc94bb764b52fd37db183472f1
            
    Examples
    --------
    Create a PowerUtility object:
    
<<<<<<< HEAD
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
=======
    >>> pu = PowerUtility()
    >>> pu
    <PowerUtility: 0 kW, 0 USD/hr>
       
    Call object to calculate cost:
        
    >>> pu(rate=500.)
    >>> pu
    <PowerUtility: 500 kW, 39.1 USD/hr>
       
    Results are accessible:
        
    >>> pu.rate, pu.cost
    (500.0, 39.1)
    
    PowerUtility objects have `consumption` and `production` attributes
    which are updated when setting the rate with the assumption that
    a positive rate means no production (only consumption) and a negative
    rate means no consumption (only consumption).
    
    >>> pu.consumption, pu.production
    (500.0, 0.0)
    
    It is possible to have both consumption and production by setting these
    attributes individually (instead of setting rate)
    
    >>> pu.production = 100. 
    >>> pu.rate
    400.0
    
    Notice how the rate is equal the consumption minus the production.
    
    """
    __slots__ = ('consumption', 'production')
>>>>>>> cd2c5013aaf9b5bc94bb764b52fd37db183472f1
    
    #: [DisplayUnits] Units of measure for IPython display
    display_units = DisplayUnits(rate='kW', cost='USD/hr')
    
    #: [float] USD/kWhr
    price = 0.0782
    
    def __init__(self):
<<<<<<< HEAD
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
=======
        #: Electricity consumption [kW]
        self.consumption = 0.
        
        #: Electricity production [kW]
        self.production = 0.
    
    @property
    def rate(self):
        """Power requirement [kW]."""
        return self.consumption - self.production
    @rate.setter
    def rate(self, rate):
        rate = float(rate)
        if rate >= 0.:
            self.consumption = rate
            self.production = 0.
        else:
            self.consumption = 0.
            self.production = -rate
    
    @property
    def cost(self):
        """Cost [USD/hr]"""
        return self.price * self.rate
    
    def __bool__(self):
        return bool(self.consumption or self.production)
    
    def __call__(self, rate):
        """Set rate in kW."""
        self.rate = rate
    
    def mix_from(self, power_utilities):
        """Mix in requirements of power utilities."""
        self.consumption = sum([i.consumption for i in power_utilities])
        self.production = sum([i.production for i in power_utilities])
    
    @classmethod
    def sum(cls, power_utilities):
        """
        Return a PowerUtility object that represents the sum of power utilities.
        """
        power_utility = cls()
        power_utility.mix_from(power_utilities)
        return power_utility
>>>>>>> cd2c5013aaf9b5bc94bb764b52fd37db183472f1
    
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