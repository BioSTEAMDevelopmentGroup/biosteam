# -*- coding: utf-8 -*-
<<<<<<< HEAD
"""
Created on Tue May 14 14:20:53 2019

@author: yoelr
=======
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
>>>>>>> cd2c5013aaf9b5bc94bb764b52fd37db183472f1
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
<<<<<<< HEAD
    
    """
    __slots__ = ('name', 'setter', 'simulate', 'element',
                 'system', 'distribution', 'units', 'baseline')
    
    def __init__(self, name, setter, simulate,
                 element, system, distribution,
                 units, baseline):
=======
    bounds : tuple[float, float]
        Lower and upper bounds of parameter.
    
    """
    __slots__ = ('name', 'setter', 'simulate', 'element',
                 'system', 'distribution', 'units', 'baseline',
                 'bounds')
    
    def __init__(self, name, setter, simulate,
                 element, system, distribution,
                 units, baseline, bounds):
>>>>>>> cd2c5013aaf9b5bc94bb764b52fd37db183472f1
        self.name = name.replace('_', ' ').capitalize()
        self.setter = setter
        self.simulate = simulate
        self.element = element
        self.system = system
        self.distribution = distribution
        self.units = units
        self.baseline = baseline
<<<<<<< HEAD
=======
        self.bounds = bounds
>>>>>>> cd2c5013aaf9b5bc94bb764b52fd37db183472f1
    
    def __call__(self, value):
        self.setter(value)
        self.simulate()
    
    def _info(self):
        return (f'{type(self).__name__}: {self.name}\n'
             + (f' element: {self.element_name}' if self.element else '')
             + (f' system: {self.system}\n' if self.system else '')
             + (f' units: {self.units}\n' if self.units else '')
             + (f' distribution: {self.distribution}\n' if self.distribution else ''))
    
                
    
               