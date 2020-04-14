# -*- coding: utf-8 -*-
"""
Created on Tue May 14 14:20:53 2019

@author: yoelr
"""
__all__ = ('Parameter',)
from ._variable import Variable

class Parameter(Variable):
    """
    Create a Parameter object that, when called, runs the setter and
    the simulate functions.
    
    Parameters
    ----------
    name : str
        Name of parameter.
    setter : function
        Should set the parameter.
    simulate : function
        Should simulate parameter effects.
    element : object
        Element associated to parameter.
    system : System
        System associated to parameter.
    distribution : chaospy.Dist
        Parameter distribution.
    units : str
        Units of parameter.
    baseline : float
        Baseline value of parameter.
    
    """
    __slots__ = ('name', 'setter', 'simulate', 'element',
                 'system', 'distribution', 'units', 'baseline')
    
    def __init__(self, name, setter, simulate,
                 element, system, distribution,
                 units, baseline):
        self.name = name.replace('_', ' ').capitalize()
        self.setter = setter
        self.simulate = simulate
        self.element = element
        self.system = system
        self.distribution = distribution
        self.units = units
        self.baseline = baseline
    
    def __call__(self, value):
        self.setter(value)
        self.simulate()
    
    def _info(self):
        return (f'{type(self).__name__}: {self.name}\n'
             + (f' element: {self.element_name}' if self.element else '')
             + (f' system: {self.system}\n' if self.system else '')
             + (f' units: {self.units}\n' if self.units else '')
             + (f' distribution: {self.distribution}\n' if self.distribution else ''))
    
                
    
               