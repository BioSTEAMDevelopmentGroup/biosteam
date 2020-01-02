# -*- coding: utf-8 -*-
"""
Created on Mon Sep  2 02:52:51 2019

@author: yoelr
"""
from ._variable import Variable
__all__ = ('Metric',)

class Metric(Variable):
    """Create a Metric object that serves as an argument for Model objects.
    
    Parameters
    ----------
    name : str
        Name of metric.
    units : str
        Metric units of measure.
    getter : function
        Should take no arguments and return the metric value.
    element : str
        Element corresponding to metric
    
    """
    __slots__ = ('name', 'units', 'getter', 'element')
    distribution = None
    def __init__(self, name, getter, units=None, element='Biorefinery'):
        self.name = name
        self.units = units
        self.getter = getter
        self.element = element
        
    def __call__(self):
        return self.getter()
    