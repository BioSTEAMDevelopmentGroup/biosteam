# -*- coding: utf-8 -*-
<<<<<<< HEAD
"""
Created on Mon Sep  2 02:52:51 2019

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
from ._variable import Variable
__all__ = ('Metric',)

class Metric(Variable):
<<<<<<< HEAD
    """Create a Metric object that serves as an argument for Model objects.
=======
    """
    Create a Metric object that serves as an argument for Model objects.
>>>>>>> cd2c5013aaf9b5bc94bb764b52fd37db183472f1
    
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
<<<<<<< HEAD
    def __init__(self, name, getter, units=None, element='Biorefinery'):
=======
    def __init__(self, name, getter=None, units=None, element='Biorefinery'):
        if not getter: return lambda getter: Metric(name, getter, units, element)
>>>>>>> cd2c5013aaf9b5bc94bb764b52fd37db183472f1
        self.name = name
        self.units = units
        self.getter = getter
        self.element = element
        
    def __call__(self):
        return self.getter()
    