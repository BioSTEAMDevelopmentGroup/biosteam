# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 13:20:52 2019

Module for estimating utility price based on:

How to Estimate Utility Costs. Engineering Practice. CHEMICAL ENGINEERING. WWW.CHE.COM, APRIL 2006. Gael D. Ulrich and Palligarnai T. VasudevanUniversity of New Hampshire

https://terpconnect.umd.edu/~nsw/chbe446/HowToEstimateUtilityCosts-UlrichVasudevan2006.pdf
    
Utility estimates are often complicated because they depend on both inflation and energy costs. This simplified approach offers a two-factor utility-cost equation and the relevant coefficients for a number of utilities. 

@author: yoelr
"""
import numpy as np

ln = np.log


# %% Classes to calculate utilities

__all__ = ('Utility', 'Electricity', 'CoolingWater')

class Utility:
    """Abstract class for utilities."""
    __slots__ = ()
    ID = 'Utility'
    fuelprice = 12 #: $/GJ
    CEPCI = 567.5 #: Chemical Engineering Plant Cost Index
    units = {}
        
    @property
    def price(self):
        """Price of utility."""
        return self._a * self.CEPCI + self._b * self.fuelprice
    
    def __repr__(self):
        units = self.units.get('price')
        if units: units = ' ' + units 
        return f"<{self.ID}: price={self.price:.2g}{units}>"
    
    def _info(self):
        units = self.units.get('price')
        if units: units = ' ' + units
        info =(f"{self.ID}:\n"
             + f" price: {self.price:.2g}{units}")  
        attr2show = self._attr2show
        if attr2show:
            for attr in attr2show:
                value = getattr(self, attr)
                units = self.units.get(attr, '')
                if units:
                    units = ' ' + units
                    info += f"\n {attr}: {value:.2g}{units}"
                else:
                    info += f"\n {attr}: {value}{units}"
        return info
    
    def show(self):
        print(self._info())
    

class Electricity(Utility):
    __slots__ = ('purchased', 'grassroots')
    ID = 'Electricity'
    _units = 'USD/kWhr'
    _attr2show = __slots__
    units = {'price': 'USD/kWhr'}
    
    def __init__(self, purchased=True, grassroots=False):
        self.purchased = purchased
        self.grassroots = grassroots
    
    @property
    def _a(self):
        if self.purchased:
            return 1.3e-4
        elif self.grassroots:
            return 1.1e-4
        else:
            return 1.3e-4
    
    @property
    def _b(self):
        if self.purchased:
            return 0.01
        else:
            return 0.011

class CoolingWater(Utility):
    
    __slots__ = ('grassroots', '_flow', '_b')
    ID = 'Cooling water'
    _attr2show = ('flow', 'grassroots')
    units = {'price': 'USD/m^3',
             'flow': 'm^3/s'}
    _flowbounds = (0.01, 10)
    
    def __init__(self, flow, grassroots=True):
        self.grassroots = grassroots
        self.flow = flow
        self._b = 0.003
    
    @property
    def _a(self):
        if self.grassroots:
            return 7e-5 + 2.5e-5 * self._flow**-1 
        else:
            return 1e-4 + 3.0e-5 * self._flow**-1
    
    @property
    def flow(self):
        """Volumetric flow rate (m3/s)."""
        return self._flow
    
    @flow.setter
    def flow(self, flow):
        lb, ub = self._flowbounds
        if not lb <= flow <= ub:
            raise ValueError("flow must be between 0.01 and 10 m3/s")
        self._flow = flow
        
    
    
    
    
    
    
    
    
    
    
    
    