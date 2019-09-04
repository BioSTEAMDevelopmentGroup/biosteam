# -*- coding: utf-8 -*-
"""
Created on Mon Sep  2 02:52:51 2019

@author: yoelr
"""

__all__ = ('Metric',)

class Metric:
    """Create a Metric object that serves as an argument for Model objects.
    
    Parameters
    ----------
    name : str
        Name of metric.
    units : str
        Metric units of measure.
    getter : function
        Should take no arguments and return metric value.
    
    """
    __slots__ = ('name', 'units', 'getter')
    
    def __init__(self, name, units, getter):
        self.name = name
        self.units = units
        self.getter = getter
    
    def describe(self):
        """Return description of metric."""
        if self.units:
            return self.name + ' (' + self.units + ')'
        else:
            return self.name
    
    def __repr__(self):
        return f"<Metric: {self.describe()}>"