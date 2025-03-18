# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from ._feature import Feature

__all__ = ('Indicator', 'indicator', 'Metric', 'metric')

class Indicator(Feature):
    """
    Create an Indicator object that serves as an argument for Model objects.
    
    Parameters
    ----------
    name : str
        Name of indicator.
    units : str
        Indicator units of measure.
    getter : function
        Should take no arguments and return the indicator value.
    element : str
        Element corresponding to indicator
    
    """
    __slots__ = ('getter', 'last_value')
    distribution = None
    def __init__(self, name, getter, units=None, element=None):
        if name is None and hasattr(getter, '__name__'): name = getter.__name__
        super().__init__(name, units, element)
        self.getter = getter
        
    def __call__(self):
        self.last_value = self.getter()
        return self.last_value
    
    def get(self):
        """Return value of indicator. This method used cached values."""
        try:
            return self.last_value
        except:
            self.last_value = self.getter()
            return self.last_value
    
    def difference(self):
        """Return the difference between the current indicator value and the last one 
        evaluated by calling this object."""
        return self.getter() - self.last_value
    
    
def indicator(getter=None, name=None, units=None, element='Biorefinery'):
    if not getter: return lambda getter: indicator(getter, name, units, element)
    return Indicator(name, getter, units, element)

metric = indicator
Metric = Indicator