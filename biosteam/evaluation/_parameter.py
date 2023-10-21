# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
__all__ = ('Parameter',)

from ._feature import Feature
from ..utils import format_title
import biosteam as bst
from inspect import signature

class Parameter(Feature):
    """
    Create a Parameter object that, when called, runs the setter and
    the simulate functions.
    
    Parameters
    ----------
    name : str
        Name of parameter.
    setter : function
        Should set the parameter.
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
    bounds : tuple[float, float]
        Lower and upper bounds of parameter.
    kind : str
        * 'design': Parameter only affects unit operation design.
        * 'coupled': Parameter affects mass and energy balances.
        * 'isolated': Parameter does not affect the system in any way.
    hook : Callable
        Should return the new parameter value given the sample.
    scale : float, optional
        The sample is multiplied by the scale before setting.
        
    """
    __slots__ = ('setter', 'system', 'distribution', 
                 'baseline', 'bounds', 'kind', 'hook',
                 'description', 'scale')
    
    def __init__(self, name, setter, element, system, distribution,
                 units, baseline, bounds, kind, hook, description, scale):
        if not name: name, *_ = signature(setter).parameters.keys()
        super().__init__(format_title(name), units, element)
        if kind is None: kind = 'isolated'
        self.setter = setter.setter if isinstance(setter, Parameter) else setter
        self.system = system
        self.distribution = distribution
        if not bounds:
            if distribution: bounds = (distribution.lower[0], distribution.upper[0])
        if bounds and baseline is None:
            baseline = 0.5 * (bounds[0] + bounds[1])
        self.baseline = baseline
        self.bounds = bounds
        self.kind = kind
        self.hook = hook
        self.description = description
        self.scale = scale
    
    @classmethod
    def sort_parameters(cls, parameters):
        if not parameters: return
        try:
            system, = set([i.system for i in parameters])
        except:
            raise ValueError('all parameters must have the same system to sort')
        unit_path = system.units
        length = len(unit_path)
        def key(parameter):
            if parameter.kind == 'coupled':
                unit = parameter.unit
                if unit: return unit_path.index(unit) 
            return length
        parameters.sort(key=key) 
    
    @property
    def unit(self):
        """Unit operation directly associated to parameter."""
        element = self.element
        if isinstance(element, bst.Unit): return element
        elif isinstance(element, bst.Stream): return element._sink
    
    @property
    def subsystem(self):
        """Subsystem directly associated to parameter."""
        system = self.system
        if not system:
            return None
        else:
            unit = self.unit
            return system._downstream_system(unit) if unit else system
    
    def simulate(self, **dyn_sim_kwargs):
        """Simulate parameter."""
        kind = self.kind
        if kind in ('design', 'cost'):
            unit = self.unit
            if not unit: raise RuntimeError(f'no unit to run {kind} algorithm')
            unit._reevaluate()
        elif kind == 'coupled':
            subsystem = self.subsystem
            if not subsystem: raise RuntimeError(f'no system to run {kind} algorithm')
            self.subsystem.simulate(**dyn_sim_kwargs)
        elif kind == 'isolated':
            pass
        else:
            raise RuntimeError(f"invalid parameter kind '{kind}'")
    
    def __call__(self, value):
        if self.hook: value = self.hook(value)
        self.setter(value if self.scale is None else value * self.scale)
        self.simulate()
    
    def _info(self):
        return (f'{type(self).__name__}: {self.name}\n'
             + (f' element: {self.element_name}' if self.element else '')
             + (f' system: {self.system}\n' if self.system else '')
             + (f' units: {self.units}\n' if self.units else '')
             + (f' distribution: {self.distribution}\n' if self.distribution else ''))
    