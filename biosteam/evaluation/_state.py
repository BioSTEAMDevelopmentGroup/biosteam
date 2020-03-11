# -*- coding: utf-8 -*-
"""
Created on Thu May  9 13:38:57 2019

@author: Guest Group
"""
from ._block import Block
from .. import Unit
from thermosteam import Stream
import numpy as np
import pandas as pd
from chaospy import J
from .evaluation_tools import load_default_parameters

__all__ = ('State',)

# %% functions

def param_unit(param):
    element = param.element
    if isinstance(element, Unit): return element
    elif isinstance(element, Stream): return element._sink

def parameter(system, element, setter, kind, name, distribution, units, baseline):
    if kind == 'coupled':
        return Block(element, system).parameter(setter, name=name,
                    distribution=distribution, units=units, baseline=baseline)
    elif kind == 'isolated':
        return Block(element, None).parameter(setter, name=name,
                    distribution=distribution, units=units, baseline=baseline)
    elif kind == 'design':
        return Block(element, None).parameter(setter, element._summary,
                    name, distribution=distribution, units=units, baseline=baseline)
    elif kind == 'cost':
        if hasattr(element, '_end'):
            def simulate():
                element._cost()
                element._end()
        else:
            simulate = element._cost
        return Block(element, None).parameter(setter, simulate, name,
                    distribution=distribution, units=units, baseline=baseline)
    raise ValueError(f"kind must be either 'coupled', 'isolated', 'design', or 'cost' (not {kind}).")

class UpdateWithSkipping:
    __slots__ = ('cache', 'params')
    def __init__(self, params):
        self.params = tuple(params)
        self.cache = None
    def __call__(self, sample):
        try:
            sim = None
            for p, x, same in zip(self.params, sample, self.cache==sample):
                if same: continue
                p.setter(x)
                if sim: continue
                if p.system: sim = p.simulate 
                else: p.simulate()
            if sim: sim()
            self.cache = sample.copy()
        except Exception as Error:
            self.cache = None
            raise Error

class UpdateWithoutSkipping:
    __slots__ = ('params',)
    def __init__(self, params):
        self.params = tuple(params)
    def __call__(self, sample):
        sim = None
        for p, x in zip(self.params, sample):
            p.setter(x)
            if sim: continue
            if p.system: sim = p.simulate 
            else: p.simulate()
        if sim: sim()


# %%
    
class State:
    """
    Create a State object that can update the `system` state given
    a sample of parameter states.
    
    Parameters
    ----------
    system : System
    skip=False : bool
        If True, skip simulation for repeated states
    
    """
    __slots__ = ('_system', # [System]
                 '_params', # list[Parameter] All parameters.
                 '_update', # [function] Updates system state.
                 '_skip') # [bool] If True, skip simulation for repeated states
    
    load_default_parameters = load_default_parameters
    
    def __init__(self, system, skip=False):
        self._system = system
        self._params = []
        self._update = None
        self._skip = skip
    
    def __len__(self):
        return len(self._params)
    
    def copy(self):
        """Return copy."""
        copy = self.__new__(type(self))
        copy._params = list(self._params)
        copy._system = self._system
        copy._update = self._update
        copy._skip = self._skip
        return copy
    
    def get_baseline_sample(self):
        return np.array([i.baseline for i in self.get_parameters()])
    
    def get_parameters(self):
        """Return parameters."""
        if not self._update: self._loadparams()
        return tuple(self._params)
    
    def get_distribution_summary(self):
        """Return dictionary of shape name-DataFrame pairs."""
        params = self.get_parameters()
        if not params: return None
        params_by_shape = {}
        shape_keys = {}
        for p in params:
            distribution = p.distribution
            if not distribution: continue
            shape = type(distribution).__name__
            if shape in params_by_shape:
                params_by_shape[shape].append(p)
            else:
                params_by_shape[shape] = [p]
                shape_keys[shape] = tuple(distribution._repr.keys())
        tables_by_shape = {}
        for shape, params in params_by_shape.items():
            data = []
            columns = ('Element', 'Name', 'Units', 'Shape', *shape_keys[shape])
            for p in params:
                distribution = p.distribution
                element = p.element_name
                name = p.name.replace(element, '')
                units = p.units or ''
                values = distribution._repr.values()
                data.append((element, name, units, shape, *values))
            tables_by_shape[shape] =  pd.DataFrame(data, columns=columns)
        return tables_by_shape
    
    def parameter(self, setter=None, element=None, kind='isolated',
                  name=None, distribution=None, units=None, baseline=None):
        """Define and register parameter.
        
        Parameters
        ----------    
        setter : function
                 Should set parameter in the element.
        element : Unit or Stream
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
        
        Notes
        -----
        
        If kind is 'coupled', account for downstream operations. Otherwise, only account for given element. If kind is 'design' or 'cost', element must be a Unit object.
        
        """
        if not setter:
            return lambda setter: self.parameter(setter, element, kind, name,
                                                 distribution, units, baseline)
        param = parameter(self._system, element, setter, kind,
                          name, distribution, units, baseline)
        self._params.append(param)
        self._erase()
        return setter
    
    def sample(self, N, rule): 
        """Return N samples from parameter distribution at given rule.
        
        Parameters
        ----------
        N : int
            Number of samples.
        rule : str
               Sampling rule.
        
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
        if not self._update: self._loadparams()
        return J(*[i.distribution for i in self._params]).sample(N, rule).transpose()
    
    def _erase(self):
        """Erase cached data."""
        self._update = None
    
    def _loadparams(self):
        """Load parameters."""
        length = len(self._system._unit_path)
        index =  self._system._unit_path.index
        self._params.sort(key=lambda x: index(param_unit(x)) if x.system else length)
        self._update = (UpdateWithSkipping if self._skip 
                        else UpdateWithoutSkipping)(self._params)
    
    def __call__(self, sample):
        """Update state given sample of parameters."""
        if not self._update: self._loadparams()
        return self._update(np.asarray(sample, dtype=float))
    
    def _repr(self):
        return f'{type(self).__name__}: {self._system}'
    
    def __repr__(self):
        return '<' + self._repr() + '>'
    
    def _info(self):
        if not self._update: self._loadparams()
        if not self._params: return f'{self._repr()}\n (No parameters)'
        lines = []
        lenghts_block = []
        lastblk = None
        for i in self._params:
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
    