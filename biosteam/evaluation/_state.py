# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
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
from .evaluation_tools import load_default_parameters

__all__ = ('State',)

# %% Fix compatibility with new chaospy version

import chaospy as cp
version_components = cp.__version__.split('.')
CP_MAJOR, CP_MINOR = int(version_components[0]), int(version_components[1])
CP4 = (CP_MAJOR, CP_MINOR) >= (4, 0)
if CP4:
    from inspect import signature
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
    del signature, save_repr_init, shapes, baseshapes, Distribution, i
del version_components, CP_MAJOR, CP_MINOR, CP4


# %% Functions

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
        Parameter.check_indices_unique(parameters)
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
        Parameter.check_index_unique(p, self._parameters)
        self._parameters.append(p)
        self._erase()
        return p
    
    def problem(self):
        """
        Return a dictionary of parameter metadata (referred to as "problem") 
        to be used for sampling by ``SALib``.
        
        See Also
        --------
        `SALib basics <https://salib.readthedocs.io/en/latest/basics.html#an-example>`_
        
        """
        params = self.get_parameters()
        return {
            'num_vars': len(params),
            'names': [i.name for i in params],
            'bounds': [i.bounds if i.bounds
                       else (i.distribution.lower[0], i.distribution.upper[0])
                       for i in params]
        }

    def sample(self, N, rule, **kwargs): 
        """
        Return N samples from parameter distribution at given rule.
        
        Parameters
        ----------
        N : int
            Number of samples.
        rule : str
            Sampling rule.
        
        Other Parameters
        ----------------
        kwargs : 
            Key word arguments for sensitivity analysis sampling.
        
        Notes
        -----
        For sampling from a joint distribution of all parameters, use the 
        following ``rule`` flags:
        
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
            
        If sampling for sensitivity analysis, use the following ``rule`` flags:
        
        +------------+--------------------------------------------+
        | key        | Description                                     |
        +============+============================================+
        | ``MORRIS`` | Samples for Morris One-at-A-Time (OAT)     |
        +------------+--------------------------------------------+
        | ``RBD``    | Chebyshev nodes adjusted to ensure nested. |
        +------------+--------------------------------------------+
        | ``FAST``   | Korobov lattice.                           |
        +------------+--------------------------------------------+
        | ``SOBOL``  | Classical (Pseudo-) Random samples.         |
        +------------+--------------------------------------------+
        
        This method relies on the ``SALib`` library for sampling schemes 
        specific to sensitivity analysis.
        
        """
        if not self._update: self._load_parameters()
        parameters = self._parameters
        problem = self.problem()
        rule = rule.upper()
        if rule in ('C', 'NC', 'K', 'R', 'RG', 'NG', 'L', 'S', 'H', 'M'):
            shape = cp.distributions
            distributions = [i.distribution for i in parameters]
            samples = shape.J(*distributions).sample(N, rule).transpose()
        else:
            if rule == 'MORRIS':
                from SALib.sample import morris as sampler
            elif rule in ('FAST', 'EFAST'):
                from SALib.sample import fast_sampler as sampler
            elif rule == 'RBD':
                from SALib.sample import latin as sampler
            elif rule == 'SOBOL':
                from SALib.sample import saltelli as sampler
            else:
                raise ValueError(f"invalid rule '{rule}'")
            samples = sampler.sample(problem, N=N, **kwargs)
        return samples
    
    def _erase(self):
        """Erase cached data."""
        self._update = None
    
    def _load_parameters(self):
        """Load parameters."""
        system = self._system
        unit_path = system.units
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
    