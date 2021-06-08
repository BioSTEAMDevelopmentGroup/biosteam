# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from ._variable import Variable
from ..utils import format_title

__all__ = ('Metric', 'metric')

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
    __slots__ = ('getter',)
    distribution = None
    def __init__(self, name, getter, units=None, element='Biorefinery'):
        if name is None and hasattr(getter, '__name__'): name = format_title(getter.__name__)
        super().__init__(name, units, element)
        self.getter = getter
        
    def __call__(self):
        return self.getter()
    
def metric(getter=None, name=None, units=None, element='Biorefinery'):
    if not getter: return lambda getter: metric(getter, name, units, element)
    return Metric(name, getter, units, element)