# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from ._parameter import Parameter
from .. import System
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

class UpdateWithSkipping:
    __slots__ = ('cache', 'parameters')
    def __init__(self, parameters):
        self.parameters = tuple(parameters)
        self.cache = None
    def __call__(self, sample, specification=None, thorough=True):
        try:
            same_arr = self.cache==sample
            for p, x, same in zip(self.parameters, sample, same_arr):
                if same: continue
                p.setter(x)
            if specification: specification()
            for p, x, same in zip(self.parameters, sample, same_arr):
                if same: continue
                p.simulate()
                if p.kind == 'coupled': break
            self.cache = sample.copy()
        except Exception as Error:
            self.cache = None
            raise Error


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
    parameters=None : Iterable[Parameter], optional
        Parameters to sample from.
    
    """
    __slots__ = ('_system', # [System]
                 '_parameters', # list[Parameter] All parameters.
                 '_N_parameters_cache', # [int] Number of parameters previously loaded
                 '_specification', # [function] Loads speficications once all parameters are set.
                 '_sample_cache') # [1d array] Last sample evaluated.
    
    load_default_parameters = load_default_parameters
    
    def __init__(self, system, specification=None, parameters=None):
        self.specification = specification
        if parameters:
            self.set_parameters(parameters)
        else:
            self._parameters = []
        self._system = system
        self._specification = specification
        self._N_parameters_cache = None
    
    specification = System.specification
    
    @property
    def system(self):
        return self._system
    
    def __len__(self):
        return len(self._parameters)
    
    def copy(self):
        """Return copy."""
        copy = self.__new__(type(self))
        copy._parameters = list(self._parameters)
        copy._system = self._system
        copy._N_parameters_cache = self._N_parameters_cache
        copy._specification = self._specification
        return copy
    
    def get_baseline_sample(self):
        """Return a pandas Series object of parameter baseline values."""
        parameters = self.get_parameters()
        sample = {}
        for p in parameters:
            baseline = p.baseline
            if baseline is None: raise RuntimeError(f'{p} has no baseline value')  
            sample[p.index] = baseline
        return pd.Series(sample)
    
    def _erase(self):
        """Erase cached data."""
        self._N_parameters_cache = self._sample_cache = None
    
    @property
    def parameters(self):
        """tuple[Parameter] All parameters added to the model."""
        return self.get_parameters()
    @parameters.setter
    def parameters(self, parameters):
        self.set_parameters(parameters)
    
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
        self._load_parameters()
        return tuple(self._parameters)
    
    def get_joint_distribution(self):
        """
        Return a chaospy joint distribution object constituting of all
        parameter objects.
        
        """
        return cp.distributions.J(*[i.distribution for i in self.get_parameters()])
    
    def get_distribution_summary(self, xlfile=None):
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
        if xlfile:
            with pd.ExcelWriter(xlfile) as writer:
                for shape, df in tables_by_shape.items():
                    df.to_excel(writer, sheet_name=shape)
        return tables_by_shape    
    
    def parameter(self, setter=None, element=None, kind='isolated', name=None, 
                  distribution=None, units=None, baseline=None, bounds=None, 
                  hook=None, description=None):
        """
        Define and register parameter.
        
        Parameters
        ----------    
        setter : function
                 Should set parameter in the element.
        element : Unit or :class:`~thermosteam.Stream`
                  Element in the system being altered.
        kind : {'coupled', 'isolated', 'design', 'cost'}, optional
            * 'coupled': parameter is coupled to the system.
            * 'isolated': parameter does not affect the system but does affect the element (if any).
            * 'design': parameter only affects design and/or cost of the element.
        name : str, optional
               Name of parameter. If None, default to argument name of setter.
        distribution : chaospy.Dist
                       Parameter distribution.
        units : str, optional
                Parameter units of measure
        baseline : float, optional
            Baseline value of parameter.
        bounds : tuple[float, float], optional
            Lower and upper bounds of parameter.
        hook : Callable, optional
            Should return the new parameter value given the sample.
        
        Notes
        -----
        If kind is 'coupled', account for downstream operations. Otherwise,
        only account for given element. If kind is 'design' or 'cost', 
        element must be a Unit object.
        
        """
        if not setter:
            return lambda setter: self.parameter(setter, element, kind, name,
                                                 distribution, units, baseline,
                                                 bounds, hook, description)
        p = Parameter(name, setter, element or 'biorefinery',
                      self.system, distribution, units, 
                      baseline, bounds, kind, hook, description)
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
            'bounds': [i.bounds for i in params]
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
        **kwargs :
            Keyword arguments passed to sampler.
        
        Notes
        -----
        This method relies on the ``chaospy`` library for sampling from 
        distributions, and the``SALib`` library for sampling schemes 
        specific to sensitivity analysis.
        
        For sampling from a joint distribution of all parameters, 
        use the following ``rule`` flags:
        
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
        | key        | Description                                |
        +============+============================================+
        | ``MORRIS`` | Samples for Morris One-at-A-Time (OAT)     |
        +------------+--------------------------------------------+
        | ``RBD``    | Chebyshev nodes adjusted to ensure nested. |
        +------------+--------------------------------------------+
        | ``FAST``   | Korobov lattice.                           |
        +------------+--------------------------------------------+
        | ``SOBOL``  | Classical (Pseudo-) Random samples.        |
        +------------+--------------------------------------------+
        
        Note that only distribution bounds (i.e. lower and upper bounds) are 
        taken into account for sensitivity analysis, the type of distribution
        (e.g., triangle vs. uniform) do not affect the sampling. 
        
        """
        rule = rule.upper()
        if rule in ('C', 'NC', 'K', 'R', 'RG', 'NG', 'L', 'S', 'H', 'M'):
            samples = self.get_joint_distribution().sample(N, rule).transpose()
        else:
            if rule == 'MORRIS':
                from SALib.sample import morris as sampler
            elif rule in ('FAST', 'EFAST'):
                from SALib.sample import fast_sampler as sampler
            elif rule == 'RBD':
                from SALib.sample import latin as sampler
            elif rule == 'SOBOL' or rule == 'SALTELLI':
                from SALib.sample import saltelli as sampler
            else:
                raise ValueError(f"invalid rule '{rule}'")
            problem = kwargs.pop('problem') if 'problem' in kwargs else self.problem()
            samples = sampler.sample(problem, N=N, **kwargs)
        return samples
    
    def _load_parameters(self):
        """Load parameters."""
        parameters = self._parameters
        N_parameters = len(parameters)
        if N_parameters != self._N_parameters_cache:
            self._N_parameters_cache = N_parameters
            Parameter.sort_parameters(parameters)
    
    def _update_state(self, sample, thorough=True):
        try:
            if thorough: 
                for f, s in zip(self._parameters, sample): f.setter(s)
                self._specification() if self._specification else self._system.simulate()
            else:
                same_arr = self._sample_cache==sample
                for p, x, same in zip(self._parameters, sample, same_arr):
                    if same: continue
                    p.setter(x)
                if self._specification: self._specification()
                for p, x, same in zip(self._parameters, sample, same_arr):
                    if same: continue
                    p.simulate()
                    if p.kind == 'coupled': break
                self._sample_cache = sample.copy()
        except Exception as Error:
            self._sample_cache = None
            raise Error
    
    def __call__(self, sample, thorough=True):
        """Update state given sample of parameters."""
        self._load_parameters()
        self._update_state(np.asarray(sample, dtype=float), self._specification)
    
    def _repr(self):
        return f'{type(self).__name__}: {self._system}'
    
    def __repr__(self):
        return '<' + self._repr() + '>'
    
    def _info(self):
        if not self._parameters: return f'{self._repr()}\n (No parameters)'
        self._load_parameters()
        lines = []
        lengths_block = []
        lastblk = None
        for i in self._parameters:
            blk = i.element_name
            element = len(blk)*' ' if blk==lastblk else blk
            lines.append(f" {element}${i.name}\n")
            lastblk = blk
            lengths_block.append(len(blk))
        maxlen_block = max(lengths_block)
        out = f'{self._repr()}\n'
        maxlen_block = max(maxlen_block, 7)
        out += ' Element:' + (maxlen_block - 7)*' ' + ' Parameter:\n'
        for newline, len_ in zip(lines, lengths_block):
            newline = newline.replace('$', ' '*(maxlen_block-len_) + '  ')
            out += newline
        return out.rstrip('\n ')
    
    def show(self):
        """Return information on metric parameters."""
        print(self._info())
    _ipython_display_ = show
    