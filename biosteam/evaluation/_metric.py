# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from ._variable import Variable
__all__ = ('Metric',)

class Metric(Variable):
    """
    Create a Metric object that serves as an argument for Model objects.
    
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
    def __init__(self, name, getter=None, units=None, element='Biorefinery'):
        if not getter: return lambda getter: Metric(name, getter, units, element)
        self.name = name
        self.units = units
        self.getter = getter
        self.element = element
        
    def __call__(self):
        return self.getter()
    