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
from warnings import warn
from inspect import signature

try:
    from chaospy import distributions as shape
except:
    warn('chaospy not installed; cannot use automation features for uncertainty and optimization', RuntimeWarning, stacklevel=2)
else:
    # Fix compatibility with new chaospy version
    import chaospy as cp
    version_components = cp.__version__.split('.')
    CP_MAJOR, CP_MINOR = int(version_components[0]), int(version_components[1])
    CP4 = (CP_MAJOR, CP_MINOR) >= (4, 0)
    if CP4:
        def save_repr_init(f):
            defaults = list(signature(f).parameters.values())[1:]
            defaults = {i.name: i.default for i in defaults}
            def init(self, *args, **kwargs):
                if not hasattr(self, '_repr'):
                    self._repr = params = defaults.copy()
                    for i, j in zip(params, args): params[i] = j
                    params.update(kwargs)
                f(self, *args, **kwargs)
            return init
        
        shapes = cp.distributions
        Distribution = cp.distributions.Distribution
        baseshapes = set([i for i in cp.distributions.baseclass.__dict__.values()
                          if isinstance(i, type) and issubclass(i, Distribution)])
        for i in shapes.__dict__.values():
            if isinstance(i, type) and issubclass(i, Distribution) and i not in baseshapes:
                i.__init__ = save_repr_init(i.__init__)
        del save_repr_init, shapes, baseshapes, Distribution, i
    del version_components, CP_MAJOR, CP_MINOR, CP4

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
    coupled : str
        Whether parameter is coupled to the system's mass and energy balances.
        This allows a ConvergenceModel to predict it's impact on recycle loops.
        Defaults to False.
    hook : Callable
        Should return the new parameter value given the sample.
        
    """
    __slots__ = ('setter', 'system', 'distribution', 
                 'baseline', 'bounds', 'coupled', 'hook',
                 'description', 'active', 'last_value')
    
    def __init__(self, name, setter, element, system, distribution,
                 units, baseline, bounds, coupled, hook, description):
        if not name: name, *_ = signature(setter).parameters.keys()
        super().__init__(format_title(name), units, element)
        self.setter = setter.setter if isinstance(setter, Parameter) else setter
        self.system = system
        if isinstance(distribution, str):
            match distribution.lower():
                case 'triangle' | 'triangular':
                    distribution = shape.Triangle(bounds[0], baseline, bounds[1])
                case 'uniform':
                    distribution = shape.Uniform(*bounds)
                case _:
                    raise ValueError(f"invalid distribution {distribution!r}; distribution must be either 'triangular' or 'uniform'")
        elif bounds is None:
            if distribution: bounds = (distribution.lower[0], distribution.upper[0])
        elif not distribution:
            distribution = shape.Uniform(*bounds)
        if bounds and baseline is None:
            baseline = 0.5 * (bounds[0] + bounds[1])
        self.distribution = distribution
        self.baseline = baseline
        self.bounds = bounds
        self.coupled = coupled
        self.hook = hook
        self.description = description
        self.active = True
        self.last_value = None
    
    @classmethod
    def sort_parameters(cls, parameters):
        if not parameters: return
        try:
            system, = set([i.system for i in parameters])
        except:
            raise ValueError('all parameters must have the same system to sort')
        if system is None: return
        unit_path = system.units
        length = len(unit_path)
        def key(parameter):
            if parameter.coupled:
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
        if self.coupled:
            subsystem = self.subsystem
            if not subsystem: raise RuntimeError('no system to simulate')
            self.subsystem.simulate(**dyn_sim_kwargs)
        else:
            unit = self.unit
            if hasattr(unit, '_reevaluate'): unit._reevaluate()
    
    def __call__(self, value):
        if self.hook: value = self.hook(value)
        self.setter(value)
        self.last_value = value
        self.simulate()
    
    def _info(self):
        return (f'{type(self).__name__}: {self.name}\n'
             + (f' element: {self.element_name}' if self.element else '')
             + (f' system: {self.system}\n' if self.system else '')
             + (f' units: {self.units}\n' if self.units else '')
             + (f' distribution: {self.distribution}\n' if self.distribution else ''))
    