# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from ._block import Block
from ._parameter import Parameter
from .. import Unit
from thermosteam import Stream
import numpy as np
import pandas as pd
from chaospy import J, Uniform
from .evaluation_tools import load_default_parameters

__all__ = ('State',)

# %% functions

def parameter_unit(parameter):
    element = parameter.element
    if isinstance(element, Unit): return element
    elif isinstance(element, Stream): return element._sink

def parameter(system, element, setter, kind, name,
              distribution, units, baseline, bounds):
    if kind == 'coupled':
        block = Block(element, system); simulate = block.simulate
    elif kind == 'isolated':
        block = Block(element, None); simulate = None
    elif kind == 'design':
        block = Block(element, None); simulate = element._summary
    elif kind == 'cost':
        block = Block(element, None); simulate = element._cost
    else:
        raise ValueError(f"kind must be either 'coupled', 'isolated', 'design', or 'cost' (not {kind}).")
    return block.parameter(setter, simulate, name, distribution, 
                           units, baseline, bounds)

class UpdateWithSkipping:
    __slots__ = ('cache', 'parameters')
    def __init__(self, parameters):
        self.parameters = tuple(parameters)
        self.cache = None
    def __call__(self, sample, specification=None):
        try:
            same_arr = self.cache==sample
            for p, x, same in zip(self.parameters, sample, same_arr):
                if same: continue
                p.setter(x)
            if specification: specification()
            for p, x, same in zip(self.parameters, sample, same_arr):
                if same: continue
                if p.system: 
                    p.simulate()
                    break
                else: p.simulate()
            self.cache = sample.copy()
        except Exception as Error:
            self.cache = None
            raise Error

class UpdateWithoutSkipping:
    __slots__ = ('parameters',)
    def __init__(self, parameters):
        self.parameters = tuple(parameters)
    def __call__(self, sample, specification=None):
        for p, x in zip(self.parameters, sample): p.setter(x)
        if specification: specification()
        for p, x in zip(self.parameters, sample): 
            if p.system: 
                p.simulate()
                break
            else: p.simulate()


# %%
    
class State:
    """
    Create a State object that can update the `system` state given
    a sample of parameter states.
    
    Parameters
    ----------
    system : System
    specification=None : Function, optional
        Loads speficications once all parameters are set.
    skip=False : bool, optional
        If True, skip simulation for repeated states.
    parameters=None : Iterable[Parameter], optional
        Parameters to sample from.
    
    """
    __slots__ = ('_system', # [System]
                 '_parameters', # list[Parameter] All parameters.
                 '_update', # [function] Updates system state.
                 '_specification', # [function] Loads speficications once all parameters are set.
                 '_skip') # [bool] If True, skip simulation for repeated states
    
    load_default_parameters = load_default_parameters
    
    def __init__(self, system, specification=None, skip=False, parameters=None):
        self.specification = specification
        if parameters:
            self.set_parameters(parameters)
        else:
            self._parameters = []
        self._system = system
        self._update = None
        self._skip = skip
    
    specification = Unit.specification
        
    def __len__(self):
        return len(self._parameters)
    
    def copy(self):
        """Return copy."""
        copy = self.__new__(type(self))
        copy._parameters = list(self._parameters)
        copy._system = self._system
        copy._update = self._update
        copy._specification = self._specification
        copy._skip = self._skip
        return copy
    
    def get_baseline_sample(self):
        return np.array([i.baseline for i in self.get_parameters()])
    
    def set_parameters(self, parameters):
        """Set parameters."""
        parameters = list(parameters)
        isa = isinstance
        for i in parameters:
            assert isa(i, Parameter), 'all elements must be Parameter objects'
        self._erase()
        self._parameters = parameters
    
    def get_parameters(self):
        """Return parameters."""
        if not self._update: self._load_parameters()
        return tuple(self._parameters)
    
    def get_distribution_summary(self):
        """Return dictionary of shape name-DataFrame pairs."""
        parameters = self.get_parameters()
        if not parameters: return None
        parameters_by_shape = {}
        shape_keys = {}
        for p in parameters:
            distribution = p.distribution
            if not distribution: continue
            shape = type(distribution).__name__
            if shape in parameters_by_shape:
                parameters_by_shape[shape].append(p)
            else:
                parameters_by_shape[shape] = [p]
                shape_keys[shape] = tuple(distribution._repr.keys())
        tables_by_shape = {}
        for shape, parameters in parameters_by_shape.items():
            data = []
            columns = ('Element', 'Name', 'Units', 'Shape', *shape_keys[shape])
            for p in parameters:
                distribution = p.distribution
                element = p.element_name
                name = p.name.replace(element, '')
                units = p.units or ''
                values = distribution._repr.values()
                data.append((element, name, units, shape, *values))
            tables_by_shape[shape] =  pd.DataFrame(data, columns=columns)
        return tables_by_shape    
    
    def parameter(self, setter=None, element=None, kind='isolated', name=None, 
                  distribution=None, units=None, baseline=None, bounds=None):
        """
        Define and register parameter.
        
        Parameters
        ----------    
        setter : function
                 Should set parameter in the element.
        element : Unit or :class:`~thermosteam.Stream`
                  Element in the system being altered.
        kind : {'coupled', 'isolated', 'design', 'cost'}
            * 'coupled': parameter is coupled to the system.
            * 'isolated': parameter does not affect the system but does affect the element (if any).
            * 'design': parameter only affects design and cost of the element.
            * 'cost': parameter only affects cost of the element.
        name : str
               Name of parameter. If None, default to argument name of setter.
        distribution : chaospy.Dist
                       Parameter distribution.
        units : str
                Parameter units of measure
        baseline : float
            Baseline value of parameter.
        bounds : tuple[float, float]
            Lower and upper bounds of parameter.
        
        Notes
        -----
        If kind is 'coupled', account for downstream operations. Otherwise,
        only account for given element. If kind is 'design' or 'cost', 
        element must be a Unit object.
        
        """
        if not setter:
            return lambda setter: self.parameter(setter, element, kind, name,
                                                 distribution, units, baseline,
                                                 bounds)
        p = parameter(self._system, element, setter, kind, name, 
                      distribution, units, baseline, bounds)
        self._parameters.append(p)
        self._erase()
        return setter
    
    def sample(self, N, rule, uniform=False): 
        """
        Return N samples from parameter distribution at given rule.
        
        Parameters
        ----------
        N : int
            Number of samples.
        rule : str
            Sampling rule.
        uniform=False : bool, optional
            Whether to assume a joint uniform distribution across parameter bounds. 
            Otherwise, sample from a joint distribution of all parameters.
        
        Notes
        -----
        Use the following ``rule`` flag for sampling:
        
        +-------+-------------------------------------------------+
        | key   | Description                                     |
        +=======+=================================================+
        | ``C`` | Roots of the first order Chebyshev polynomials. |
        +-------+-------------------------------------------------+
        | ``NC``| Chebyshev nodes adjusted to ensure nested.      |
        +-------+-------------------------------------------------+
        | ``K`` | Korobov lattice.                                |
        +-------+-------------------------------------------------+
        | ``R`` | Classical (Pseudo-)Random samples.              |
        +-------+-------------------------------------------------+
        | ``RG``| Regular spaced grid.                            |
        +-------+-------------------------------------------------+
        | ``NG``| Nested regular spaced grid.                     |
        +-------+-------------------------------------------------+
        | ``L`` | Latin hypercube samples.                        |
        +-------+-------------------------------------------------+
        | ``S`` | Sobol low-discrepancy sequence.                 |
        +-------+-------------------------------------------------+
        | ``H`` | Halton low-discrepancy sequence.                |
        +-------+-------------------------------------------------+
        | ``M`` | Hammersley low-discrepancy sequence.            |
        +-------+-------------------------------------------------+
            
        """
        if not self._update: self._load_parameters()
        parameters = self._parameters
        if uniform:
            distributions = [Uniform(*i.bounds) for i in parameters]
        else:
            distributions = [i.distribution for i in parameters]
        return J(*distributions).sample(N, rule).transpose()
    
    def _erase(self):
        """Erase cached data."""
        self._update = None
    
    def _load_parameters(self):
        """Load parameters."""
        system = self._system
        unit_path = system._unit_path + list(system._facilities)
        length = len(unit_path)
        index =  unit_path.index
        self._parameters.sort(key=lambda x: index(parameter_unit(x)) if x.system else length)
        self._update = (UpdateWithSkipping if self._skip 
                        else UpdateWithoutSkipping)(self._parameters)
    
    def __call__(self, sample):
        """Update state given sample of parameters."""
        if not self._update: self._load_parameters()
        self._update(np.asarray(sample, dtype=float), self._specification)
    
    def _repr(self):
        return f'{type(self).__name__}: {self._system}'
    
    def __repr__(self):
        return '<' + self._repr() + '>'
    
    def _info(self):
        if not self._update: self._load_parameters()
        if not self._parameters: return f'{self._repr()}\n (No parameters)'
        lines = []
        lenghts_block = []
        lastblk = None
        for i in self._parameters:
            blk = i.element_name
            element = len(blk)*' ' if blk==lastblk else blk
            lines.append(f" {element}${i.name}\n")
            lastblk = blk
            lenghts_block.append(len(blk))
        maxlen_block = max(lenghts_block)
        out = f'{self._repr()}\n'
        maxlen_block = max(maxlen_block, 7)
        out += ' Element:' + (maxlen_block - 7)*' ' + ' Parameter:\n'
        for newline, len_ in zip(lines, lenghts_block):
            newline = newline.replace('$', ' '*(maxlen_block-len_) + '  ')
            out += newline
        return out.rstrip('\n ')
    
    def show(self):
        """Return information on metric parameters."""
        print(self._info())
    _ipython_display_ = show
    