# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-, Yoel Cortes-Pena <yoelcortes@gmail.com>,
#                      Yalin Li <mailto.yalin.li@gmail.com>,
#                      Sarang Bhagwat <sarangb2@gmail.com>
#
# This module implements a filtering feature from the stats module of the QSDsan library:
# QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems
# Copyright (C) 2020-, Yalin Li <mailto.yalin.li@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

from scipy.spatial.distance import cdist
from scipy.optimize import shgo, differential_evolution
import numpy as np
import pandas as pd
from chaospy import distributions as shape
from ._indicator import Indicator
from ._feature import MockFeature
from ._utils import var_indices, var_columns, indices_to_multiindex
from ._prediction import ConvergenceModel
from .._unit import Unit
from biosteam.exceptions import FailedEvaluation
from warnings import warn
from collections.abc import Sized
from biosteam.utils import Timer
from typing import Optional, Callable
from ._parameter import Parameter
from .evaluation_tools import load_default_parameters
import pickle

__all__ = ('Model', 'EasyInputModel')

def replace_nones(values, replacement):
    for i, j in enumerate(values):
        if j is None: values[i] = replacement
    return values

def codify(statement):
    statement = replace_apostrophes(statement)
    statement = replace_newline(statement)
    return statement

def replace_newline(statement):
    statement = statement.replace('\n', ';')
    return statement

def replace_apostrophes(statement):
    statement = statement.replace('’', "'").replace('‘', "'").replace('“', '"').replace('”', '"')
    return statement

def create_function(code, namespace):
    def wrapper_fn(statement):
        def f(x):
            namespace['x'] = x
            exec(codify(statement), namespace)
        return f
    function = wrapper_fn(code)
    return function

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

# %% Simulation of process systems

class Model:
    """
    Create a Model object that allows for evaluation over a sample space.
    
    Parameters
    ----------
    system : System
        Should reflect the model state.
    indicators : tuple[Indicator]
        Indicators to be evaluated by model.
    specification=None : Function, optional
        Loads specifications once all parameters are set. Specification should 
        simulate the system as well.
    parameters=None : Iterable[Parameter], optional
        Parameters to sample from.
    exception_hook : callable(exception, sample)
        Function called after a failed evaluation. The exception hook should 
        return either None or indicator values given the exception and sample.

    """
    __slots__ = (
        '_system', # [System]
        '_parameters', # list[Parameter] All parameters.
        '_optimized_parameters', # list[Parameter] Parameters being optimized.
        '_specification', # [function] Loads specifications once all parameters are set.
        'table',            # [DataFrame] All arguments and results.
        'retry_evaluation', # [bool] Whether to retry evaluation if it fails
        'convergence_model',# [ConvergenceModel] Prediction model for recycle convergence.
        '_indicators',         # tuple[Indicator] Indicators to be evaluated by model.
        '_index',           # list[int] Order of sample evaluation for performance.
        '_samples',         # [array] Argument sample space.
        '_exception_hook',  # [callable(exception, sample)] Should return either None or indicator value given an exception and the sample.
    )
    default_optimizer_options = {
        'shgo': dict(f_tol=1e-3, minimizer_kwargs=dict(f_tol=1e-3)),
        'differential evolution': {'seed': 0, 'popsize': 12, 'tol': 1e-3}
    }
    default_optimizer = 'shgo'
    default_convergence_model = None # Optional[str] Default convergence model
    load_default_parameters = load_default_parameters
    
    @property
    def specification(self) -> Callable:
        """Process specification."""
        return self._specification
    @specification.setter
    def specification(self, specification):
        if specification:
            if callable(specification):
                self._specification = specification
            else:
                raise AttributeError(
                    "specification must be callable or None; "
                   f"not a '{type(specification).__name__}'"
                )
        else:
            self._specification = None
    
    @property
    def system(self):
        return self._system
    
    def __len__(self):
        return len(self._parameters)
    
    def get_baseline_scenario(self, parameters=None, array=False):
        """Return a pandas Series object of parameter baseline values."""
        if parameters is None: parameters = self._parameters
        sample = [] if array else {}
        for p in parameters:
            baseline = p.baseline
            if baseline is None: raise RuntimeError(f'{p} has no baseline value')  
            if array: sample.append(baseline)
            else: sample[p.index] = baseline
        return np.array(sample) if array else pd.Series(sample)
    get_baseline_sample = get_baseline_scenario
    
    @property
    def optimized_parameters(self):
        """tuple[Parameter] All parameters optimized in the model."""
        return tuple(self._optimized_parameters)
    @optimized_parameters.setter
    def optimized_parameters(self, parameters):
        self._optimized_parameters = parameters = list(parameters)
        isa = isinstance
        for i in parameters:
            assert isa(i, Parameter), 'all elements must be Parameter objects'
        Parameter.check_indices_unique(self.features)
    
    @property
    def parameters(self):
        """tuple[Parameter] All parameters added to the model."""
        return self.get_parameters()
    @parameters.setter
    def parameters(self, parameters):
        self.set_parameters(parameters)
    
    def set_parameters(self, parameters):
        """Set parameters."""
        self._parameters = parameters = list(parameters)
        isa = isinstance
        for i in parameters:
            assert isa(i, Parameter), 'all elements must be Parameter objects'
        Parameter.check_indices_unique(self.features)
    
    def parameters_from_df(self, df_or_filename, namespace=None):
        """
        Load a list (from a DataFrame or spreadsheet) of distributions and statements
        to load values for user-selected parameters.
        
        Parameters
        ----------
        df_or_filename : pandas.DataFrame or file path to a spreadsheet with the following column titles.
                         
            * 'Parameter name' [String] Name of the parameter.
            
            * 'Element' [String] 
                        
            * 'Units' [String] 
                        
            * 'Baseline' [float] The baseline value of the parameter.
                        
            * 'Shape' {'Uniform', 'Triangular'} The shape of the parameter distribution.
                        
            * 'Lower' [float] The lower value defining the shape of the parameter distribution.
                        
            * 'Midpoint' [float] The midpoint value defining the shape of a 'Triangular' parameter distribution.
                        
            * 'Upper' [float] The upper value defining the shape of the parameter distribution.
                        
            * 'Load statement' [String] A statement executed to load the value of the parameter. 
            
        namespace : dict, optional
            Dictionary used to update the namespace accessed when executing
            statements to load values into model parameters. Defaults to the
            system's flowsheet dict.
                        
        """
        
        df = df_or_filename
        if type(df) is not pd.DataFrame:
            try: 
                df = pd.read_excel(df_or_filename)
            except:
                df = pd.read_csv(df_or_filename)
                
        if namespace is None: namespace = {}
        namespace = self.system.flowsheet.to_dict() | namespace
        
        param = self.parameter
        
        for i, row in df.iterrows():
            name = row['Parameter name']
            element = row['Element'] # currently only compatible with String elements
            try:
                coupled = row['Coupled']
            except:
                coupled = row['Kind'] == 'coupled'
            units = row['Units']
            baseline = row['Baseline']
            shape_data = row['Shape']
            lower, midpoint, upper = row['Lower'], row['Midpoint'], row['Upper']
            load_statements = row['Load statement']
            
            D = None
            if shape_data.lower() in ['triangular', 'triangle',]:
                D = shape.Triangle(lower, midpoint, upper)
            elif shape_data.lower() in ['uniform',]:
                if not str(midpoint)=='nan':
                    raise ValueError(f"The parameter distribution for {name} ({element}) is 'Uniform' but was associated with a given midpoint value.")
                D = shape.Uniform(lower, upper)
                
            param(name=name, 
                  setter=create_function(load_statements, namespace), 
                  element=element, 
                  coupled=coupled, 
                  units=units,
                  baseline=baseline, 
                  distribution=D)
            
    def get_parameters(self):
        """Return parameters."""
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
    
    def optimized_parameter(self, *args, **kwargs):
        return self.parameter(*args, **kwargs, optimized=True, coupled=True)
    
    def parameter(self, 
            setter: Optional[Callable]=None,
            element: Optional[Unit]=None, 
            coupled: Optional[bool]=None,
            name: Optional[str]=None, 
            distribution: Optional[str|shape.baseclass.Distribution]=None, 
            units: Optional[str]=None, 
            baseline: Optional[float]=None, 
            bounds: Optional[tuple[float, float]]=None, 
            hook: Optional[Callable]=None, 
            description: Optional[str]=None, 
            optimized=False, 
            kind=None, # For backwards compatibility
            safe=False,
        ):
        """
        Define and register parameter.
        
        Parameters
        ---------*    
        setter : 
            Should set parameter in the element.
        element : 
            Element in the system being altered.
        coupled : 
            Whether parameter is coupled to the system's mass and energy balances.
            This allows a ConvergenceModel to predict it's impact on recycle loops.
            Defaults to False.
        name : 
            Name of parameter. If None, default to argument name of setter.
        distribution : 
            Parameter distribution.
        units : 
            Parameter units of measure
        baseline : 
            Baseline value of parameter.
        bounds : 
            Lower and upper bounds of parameter.
        hook : 
            Should return the new parameter value given the sample.
        
        """
        if kind == 'coupled': coupled = True
        elif coupled is None: coupled = False
        if isinstance(setter, Parameter):
            if element is None: element = setter.element
            if coupled is None: coupled = setter.coupled
            if name is None: name = setter.name
            if distribution is None: distribution = setter.distribution
            if units is None: units = setter.units
            if baseline is None: baseline = setter.baseline
            if bounds is None: bounds = setter.bounds
            if hook is None: hook = setter.hook
            if description is None: description = setter.description
            setter = setter.setter
        elif isinstance(setter, MockFeature):
            if element is None: element = setter.element
            if name is None: name = setter.name
            if units is None: units = setter.units
        elif not setter:
            return lambda setter: self.parameter(setter, element, coupled, name,
                                                 distribution, units, baseline,
                                                 bounds, hook, description, optimized)
        p = Parameter(name, setter, element,
                      self.system, distribution, units, 
                      baseline, bounds, coupled, hook, description)
        if safe:
            Parameter.check_index_unique(p, self.features)
        else:
            key = (p.element_name, p.name)
            dct = {(i.element_name, i.name): i for i in self.features}
            if key in dct:
                old_p = dct[key]
                try:
                    if optimized:
                        self._optimized_parameters.remove(old_p)
                    else:
                        self._parameters.remove(old_p)
                except:
                    raise ValueError(
                             "each feature must have a unique element and name; "
                            f"feature with element {repr(p.element)} "
                            f"and name {repr(p.name)} already present"
                        )
        if optimized:
            self._optimized_parameters.append(p)
        else:
            self._parameters.append(p)
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
                from SALib.sample import sobol as sampler
            else:
                raise ValueError(f"invalid rule '{rule}'")
            problem = kwargs.pop('problem') if 'problem' in kwargs else self.problem()
            samples = sampler.sample(problem, N=N, **kwargs)
        return samples
    
    def _objective_function(self, sample, loss, parameters, convergence_model=None, **kwargs):
        for f, s in zip(parameters, sample): 
            f.setter(s if f.scale is None else f.scale * s)
        if convergence_model:
            with convergence_model.practice(sample, parameters):
                self._specification() if self._specification else self._system.simulate(**kwargs)
        return loss()
    
    def _update_state(self, sample, convergence_model=None, **kwargs):
        for i, (f, value) in enumerate(zip(self._parameters, sample)): 
            if f.active: 
                f.setter(value)
                f.last_value = value
            else:
                sample[i] = f.last_value
        if convergence_model:
            with convergence_model.practice(sample):
                return self._specification() if self._specification else self._system.simulate(**kwargs)
        else:
            return self._specification() if self._specification else self._system.simulate(**kwargs)
    
    def _evaluate_sample(self, sample, convergence_model=None, **kwargs):
        state_updated = False
        try:
            self._update_state(sample, convergence_model, **kwargs)
            state_updated = True
            return [i() for i in self.indicators]
        except Exception as exception:
            if self.retry_evaluation and not state_updated:
                self._reset_system()
                try:
                    self._update_state(sample, convergence_model, **kwargs)
                    return [i() for i in self.indicators]
                except Exception as new_exception: 
                    exception = new_exception
            if self._exception_hook: 
                values = self._exception_hook(exception, sample)
                self._reset_system()
                if isinstance(values, Sized) and len(values) == len(self.indicators):
                    return values
                elif values is not None:
                    raise RuntimeError('exception hook must return either None or '
                                       'an array of indicator values for the given sample')
            return self._failed_evaluation()
    
    def __init__(self, system, indicators=None, specification=None, 
                 parameters=None, retry_evaluation=None, exception_hook=None):
        self.specification = specification
        if parameters:
            self.set_parameters(parameters)
        else:
            self._parameters = []
        self._optimized_parameters = []
        self._system = system
        self._specification = specification
        self.indicators = indicators or ()
        self.exception_hook = 'warn' if exception_hook is None else exception_hook 
        self.retry_evaluation = bool(system) if retry_evaluation is None else retry_evaluation
        self.table = None
        self._erase()
        
    def copy(self):
        """Return copy."""
        copy = self.__new__(type(self))
        copy._parameters = self._parameters.copy()
        copy._optimized_parameters = self._optimized_parameters.copy()
        copy._system = self._system
        copy._specification = self._specification
        copy._indicators = self._indicators
        if self.table is None:
            copy._samples = copy.table = None
        else:
            copy.table = self.table.copy()
            copy._samples = self._samples
            copy._index = self._index
        return copy
    
    def _erase(self):
        """Erase cached data."""
        self._samples = None
    
    @property
    def exception_hook(self):
        """
        [callable(exception, sample)] Function called after a failed 
        evaluation. The exception hook should return either None or indicator 
        values given the exception and sample.
        
        """
        return self._exception_hook
    @exception_hook.setter
    def exception_hook(self, exception_hook):
        if not exception_hook or callable(exception_hook):
            self._exception_hook = exception_hook
            return
        if isinstance(exception_hook, str):
            if exception_hook == 'ignore':
                self._exception_hook = lambda exception, sample: None
            elif exception_hook == 'warn':
                self._exception_hook = lambda exception, sample: warn(FailedEvaluation(f"[{type(exception).__name__}] {exception}"), stacklevel=6)
            elif exception_hook == 'raise':
                def raise_exception(exception, sample): raise exception from None
                self._exception_hook = raise_exception
            else:
                raise ValueError(f"invalid exception hook name '{exception_hook}'; "
                                  "valid names include 'ignore', 'warn', and 'raise'")
        else:
            raise ValueError('exception hook must be either a callable, a string, or None')
    
    @property
    def indicators(self):
        """tuple[Indicator] Indicators to be evaluated by model."""
        return tuple(self._indicators)
    @indicators.setter
    def indicators(self, indicators):
        self._indicators = indicators = list(indicators)
        isa = isinstance
        for i in indicators:
            if not isa(i, Indicator):
                raise ValueError(f"indicators must be '{Indicator.__name__}' "
                                 f"objects, not '{type(i).__name__}'")
        Indicator.check_indices_unique(self.features)
    
    # Backwards compatibility
    metrics = indicators
    @property
    def _metrics(self): return self._indicators
    @_metrics.setter
    def _metrics(self, metrics): self._indicators = metrics
    
    @property
    def features(self):
        return (*self._parameters, *self._optimized_parameters, *self._indicators)
    
    def indicator(self, getter=None, name=None, units=None, element=None, safe=False):
        """
        Define and register indicator.
        
        Parameters
        ----------    
        getter : function, optional
                 Should return indicator.
        name : str, optional
               Name of indicator. If None, defaults to the name of the getter.
        units : str, optional
                Indicator units of measure
        element : object, optional
                  Element being evaluated. Works mainly for bookkeeping. 
                  Defaults to 'Biorefinery'.
        
        Notes
        -----
        This method works as a decorator.
        
        """
        if isinstance(getter, Indicator):
            if name is None: name = getter.name
            if units is None: units = getter.units
            if element is None: element = getter.element
            getter = getter.getter
        elif isinstance(getter, MockFeature):
            if element is None: element = getter.element
            if name is None: name = getter.name
            if units is None: units = getter.units
        elif not getter: 
            return lambda getter: self.indicator(getter, name, units, element)
        indicator = Indicator(name, getter, units, element)
        if safe:
            Indicator.check_index_unique(indicator, self.features, safe)
        else:
            key = (indicator.element_name, indicator.name)
            dct = {(i.element_name, i.name): i for i in self.features}
            if key in dct:
                old_indicator = dct[key]
                try:
                    self._indicators.remove(old_indicator)
                except:
                    raise ValueError(
                             "each feature must have a unique element and name; "
                            f"feature with element {repr(indicator.element)} "
                            f"and name {repr(indicator.name)} already present"
                        )
        self._indicators.append(indicator)
        return indicator 
    metric = indicator
    
    def _sample_hook(self, samples, parameters):
        if any([p.hook for p in parameters]):
            return np.array(
                [[(i if p.hook is None else p.hook(i))
                  for p, i in zip(parameters, row)]
                 for row in samples]
            )
        else:
            return samples
    
    def _load_sample_order(self, samples, parameters, distance, algorithm=None):
        """
        Sort simulation order to optimize convergence speed
        by minimizing perturbations to the system between simulations.
        
        """
        N_samples = samples.shape[0]
        if N_samples < 2: 
            self._index = list(range(N_samples))
            return
        if distance is None: distance = 'cityblock'
        length = samples.shape[0]
        columns = [i for i, parameter in enumerate(parameters) if parameter.coupled]
        samples = samples.copy()
        samples = samples[:, columns]
        samples_min = samples.min(axis=0)
        samples_max = samples.max(axis=0)
        samples_diff = samples_max - samples_min
        normalized_samples = (samples - samples_min) / samples_diff
        nearest_arr = cdist(normalized_samples, normalized_samples, metric=distance)
        if algorithm is None: algorithm = 'nearest neighbor'
        if algorithm == 'nearest neighbor':
            nearest_arr = np.argsort(nearest_arr, axis=1)
            remaining = set(range(length))
            self._index = index = [0]
            nearest = nearest_arr[0, 1]
            index.append(nearest)
            remaining.remove(0)
            remaining.remove(nearest)
            N_remaining = length - 2
            N_remaining_last = N_remaining + 1
            while N_remaining:
                assert N_remaining_last != N_remaining, "issue in sorting algorithm"
                N_remaining_last = N_remaining
                for i in range(1, length):
                    next_nearest = nearest_arr[nearest, i]
                    if next_nearest in remaining:
                        nearest = next_nearest
                        remaining.remove(nearest)
                        index.append(nearest)
                        N_remaining -= 1
                        break
        else:
            raise ValueError(f"algorithm {algorithm!r} is not available (yet); "
                              "only 'nearest neighbor' is available")
        
    def load_samples(self, samples=None, sort=None, file=None, 
                     autoload=None, autosave=None, distance=None):
        """
        Load samples for evaluation.
        
        Parameters
        ----------
        samples : numpy.ndarray, dim=2, optional
            All parameter samples to evaluate.
        sort : bool, optional
            Whether to internally sort the samples to optimize convergence speed
            by minimizing perturbations to the system between simulations. The
            optimization problem is equivalent to the travelling salesman problem;
            each scenario of (normalized) parameters represent a point in the path.
            Defaults to True.
        file : str, optional
            File to load/save samples and simulation order to/from.
        autosave : bool, optional
            Whether to save samples and simulation order to file (when not loaded from file).
        autoload : bool, optional
            Whether to load samples and simulation order from file (if possible).
        distance : str, optional
            Distance indicator used for sorting. Defaults to 'cityblock'.
            See scipy.spatial.distance.cdist for options.
        algorithm : str, optional
            Algorithm used for sorting. Defaults to 'nearest neighbor'.
            Note that neirest neighbor is a greedy algorithm that is known to result, 
            on average, in paths 25% longer than the shortest path.
        
        """
        parameters = self._parameters
        if autoload:
            try:
                with open(file, "rb") as f: (self._samples, self._index) = pickle.load(f)
            except FileNotFoundError: pass
            else:
                if (samples is None or (samples.shape == self._samples.shape and (samples == self._samples).all())):
                    indicators = self._indicators
                    empty_indicator_data = np.zeros((len(samples), len(indicators)))
                    self.table = pd.DataFrame(np.hstack((samples, empty_indicator_data)),
                                              columns=var_columns(parameters + indicators),
                                              dtype=float)
                    return 
        if not isinstance(samples, np.ndarray):
            raise TypeError(f'samples must be an ndarray, not a {type(samples).__name__} object')
        if samples.ndim == 1:
            samples = samples[:, np.newaxis]
        elif samples.ndim != 2:
            raise ValueError('samples must be 2 dimensional')
        N_parameters = len(parameters)
        if samples.shape[1] != N_parameters:
            raise ValueError(f'number of parameters in samples ({samples.shape[1]}) must be equal to the number of parameters ({N_parameters})')
        indicators = self._indicators
        samples = self._sample_hook(samples, parameters)
        if sort is None: sort = True
        if sort and any([i.coupled for i in parameters]): 
            self._load_sample_order(samples, parameters, distance)
        else:
            self._index = list(range(samples.shape[0]))
        empty_indicator_data = np.zeros((len(samples), len(indicators)))
        self.table = pd.DataFrame(np.hstack((samples, empty_indicator_data)),
                                  columns=var_columns(parameters + indicators),
                                  dtype=float)
        self._samples = samples
        if autosave:
            obj = (samples, self._index)
            try:
                with open(file, 'wb') as f: pickle.dump(obj, f)
            except FileNotFoundError:
                import os
                head, tail = os.path.split(file)
                os.mkdir(head)
                with open(file, 'wb') as f: pickle.dump(obj, f)
            
    def single_point_sensitivity(self, 
            etol=0.01, array=False, parameters=None, indicators=None, evaluate=None, 
            **kwargs
        ):
        if parameters is None: parameters = self.parameters
        bounds = [i.bounds for i in parameters]
        sample = [i.baseline for i in parameters]
        N_parameters = len(parameters)
        index = range(N_parameters)
        if indicators is None: indicators = self.indicators
        N_indicators = len(indicators)
        values_lb = np.zeros([N_parameters, N_indicators])
        values_ub = values_lb.copy()
        if evaluate is None: evaluate = self._evaluate_sample
        baseline_1 = np.array(evaluate(sample, **kwargs))
        sys = self.system
        if not sys.isdynamic: kwargs['recycle_data'] = sys.get_recycle_data()
        for i in index:
            sample_lb = sample.copy()
            sample_ub = sample.copy()
            lb, ub = bounds[i]
            p = parameters[i]
            hook = p.hook
            if hook is not None:
                lb = hook(lb)
                ub = hook(ub)
            sample_lb[i] = lb
            sample_ub[i] = ub
            values_lb[i, :] = evaluate(sample_lb, **kwargs)
            values_ub[i, :] = evaluate(sample_ub, **kwargs)
        baseline_2 = np.array(evaluate(sample, **kwargs))
        error = np.abs(baseline_2 - baseline_1)
        index, = np.where(error > 1e-6)
        error = error[index]
        relative_error = error / np.maximum.reduce([np.abs(baseline_1[index]), np.abs(baseline_2[index])])
        for i, idx in enumerate(index):
            if relative_error[i] > etol:
                print(RuntimeError(
                    f"inconsistent model; {indicators[idx]} has a value of "
                    f"{baseline_1[idx]} before evaluating sensitivity and "
                    f"{baseline_2[idx]} after"
                ))
        baseline = 0.5 * (baseline_1 + baseline_2)
        if array:
            return baseline, values_lb, values_ub
        else:
            indicator_index = var_columns(indicators)
            baseline = pd.Series(baseline, index=indicator_index)
            df_lb = pd.DataFrame(values_lb, 
                                index=var_columns(parameters),
                                columns=indicator_index)
            df_ub = df_lb.copy()
            df_ub[:] = values_ub
            return baseline, df_lb, df_ub
    
    def load_pickled_results(self, file=None, safe=True):
        table = self.table
        with open(file, "rb") as f:
            number, values, table_index, table_columns = pickle.load(f)
        if (table_index != table.index).any() or (table_columns != table.columns).any():
            if safe: raise ValueError('table layout does not match autoload file')
        del table_index, table_columns
        indicators = self._indicators
        table[var_indices(indicators)] = replace_nones(values, [np.nan] * len(indicators))
    
    def optimize(self, 
            loss, 
            parameters=None, 
            method=None, 
            convergence_model=None, 
            options=None,
        ):
        if parameters is None:
            parameters = self._optimized_parameters
        if method is None:
            method = self.default_optimizer
        else:
            method = method.lower()
        if options is None and method in self.default_optimizer_options: 
            options = self.default_optimizer_options[method]
        if isinstance(convergence_model, str):
            convergence_model = ConvergenceModel(
                system=self.system,
                predictors=parameters,
                model_type=convergence_model,
            )
        elif convergence_model is None:
            convergence_model = ConvergenceModel(
                system=self.system,
                predictors=parameters,
                model_type=self.default_convergence_model,
            )
        objective_function = self._objective_function
        args = (loss, parameters, convergence_model)
        bounds = np.array([p.bounds for p in parameters])
        if method == 'shgo':
            result = shgo(
                objective_function, bounds, args, options=options,
            )
        elif method == 'differential evolution':
            result = differential_evolution(
                objective_function, bounds, args, **options
            )
        else:
            raise ValueError(f'invalid optimization method {method!r}')
        return result, convergence_model
    
    def evaluate(self, notify=0, file=None, autosave=0, autoload=False,
                 convergence_model=None, **kwargs):
        """
        Evaluate indicators over the loaded samples and save values to `table`.
        
        Parameters
        ----------
        notify=0 : int, optional
            If 1 or greater, notify elapsed time after the given number of sample evaluations. 
        file : str, optional
            Name of file to save/load pickled evaluation results.
        autosave : int, optional
            If 1 or greater, save pickled evaluation results after the given number of sample evaluations.
        autoload : bool, optional
            Whether to load pickled evaluation results from file.
        convergence_model : ConvergencePredictionModel, optional
            A prediction model for accelerated system convergence. Defaults
            to no convergence model and the last solution
            as the initial guess for the next scenario. If a string is passed, 
            a ConvergenceModel will be created using that model type.
        kwargs : dict
            Any keyword arguments passed to :func:`biosteam.System.simulate`.
        
        Warning
        -------
        Any changes made to either the model or the samples will not be accounted
        for when autoloading and may lead to misleading results.
        
        """
        samples = self._samples
        if samples is None: raise RuntimeError('must load samples before evaluating')
        evaluate_sample = self._evaluate_sample
        table = self.table
        if isinstance(convergence_model, str):
            convergence_model = ConvergenceModel(
                system=self.system,
                predictors=self.parameters,
                model_type=convergence_model,
            )
        if notify:
            timer = Timer()
            timer.start()
            count = [0]
            def evaluate(sample, convergence_model, count=count, **kwargs):
                count[0] += 1
                values = evaluate_sample(sample, convergence_model, **kwargs)
                if not count[0] % notify:
                    print(f"{count} Elapsed time: {timer.elapsed_time:.0f} sec")
                return values
        else:
            evaluate = evaluate_sample
        N_samples, _ = samples.shape
        if autoload: 
            try:
                with open(file, "rb") as f:
                    number, values, table_index, table_columns = pickle.load(f)
                if (table_index != table.index).any() or (table_columns != table.columns).any():
                    raise ValueError('table layout does not match autoload file')
                del table_index, table_columns
                index = self._index[number:]
            except:
                number = 0
                index = self._index
                values = [None] * N_samples 
            else:
                if notify: count[0] = number
        else:
            number = 0
            index = self._index
            values = [None] * N_samples
        export = 'export_state_to' in kwargs
        layout = table.index, table.columns
        try:
            for number, i in enumerate(index, number + 1): 
                if export: kwargs['sample_id'] = i
                values[i] = evaluate(samples[i], convergence_model, **kwargs)
                if autosave and not number % autosave: 
                    obj = (number, values, *layout)
                    try:
                        with open(file, 'wb') as f: pickle.dump(obj, f)
                    except FileNotFoundError:
                        import os
                        head, tail = os.path.split(file)
                        os.mkdir(head)
                        with open(file, 'wb') as f: pickle.dump(obj, f)
        finally:
            table[var_indices(self._indicators)] = replace_nones(values, [np.nan] * len(self.indicators))
    
    def _reset_system(self):
        if self._system is None: return 
        self._system.empty_outlet_streams()
        self._system.reset_cache()
    
    def _failed_evaluation(self):
        self._reset_system()
        return [np.nan] * len(self.indicators)
    
    def indicators_at_baseline(self):
        """Return indicator values at baseline sample."""
        baseline = self.get_baseline_scenario()
        return self(baseline)
    metrics_at_baseline = indicators_at_baseline
    
    def evaluate_across_coordinate(self,
            name, f_coordinate, coordinate, *, 
            xlfile=None, notify=0, notify_coordinate=True,
            multi_coordinate=False, 
            simulation_independent_coordinate=False,
            f_evaluate=None
        ):
        """
        Evaluate across coordinate and save sample indicators.
        
        Parameters
        ----------
        name : str or tuple[str]
            Name of coordinate
        f_coordinate : function
            Should change state of system given the coordinate.
        coordinate : array
            Coordinate values.
        xlfile : str, optional
            Name of file to save. File must end with ".xlsx"
        rule : str, optional
            Sampling rule. Defaults to 'L'.
        notify_coordinate : bool, optional
            If True, notify elapsed time after each coordinate evaluation. Defaults to True.
        notify : int, optional
            Notify elapsed time after given number of scenario evaluations.
        f_evaluate : callable, optional
            Function to evaluate model. Defaults to evaluate method.
        
        """
        if (isinstance(f_coordinate, Parameter)
            and f_coordinate in self.parameters
            and f_coordinate.active):
            active = f_coordinate.active
            f_coordinate.active = False
            f_coordinate = f_coordinate.setter
            try:
                return self.evaluate_across_coordinate(
                    name, f_coordinate, coordinate,
                    xlfile=xlfile, notify=notify, notify_coordinate=notify_coordinate,
                    multi_coordinate=multi_coordinate,
                    simulation_independent_coordinate=simulation_independent_coordinate,
                    f_evaluate=f_evaluate
                )
            finally:
                f_coordinate.active = active
        N_points = len(coordinate)
        if simulation_independent_coordinate:
            if f_evaluate is not None: 
                raise ValueError(
                    'cannot pass `f_evaluate` if coordinate is '
                    'independent from simulation'
                )
            f_evaluate = self.evaluate
            samples = self._samples
            if samples is None: raise RuntimeError('must load samples before evaluating')
            evaluate_sample = self._evaluate_sample
            table = self.table
            if notify:
                timer = Timer()
                timer.start()
                count = [0]
                def evaluate(sample, count=count):
                    count[0] += 1
                    values = evaluate_sample(sample)
                    if not count[0] % notify:
                        print(f"{count} Elapsed time: {timer.elapsed_time:.0f} sec")
                    return values
            else:
                evaluate = evaluate_sample
            N_samples, _ = samples.shape
            number = 0
            index = self._index
            N_samples, _ = self.table.shape
            indicator_indices = var_indices(self.indicators)
            shape = (N_samples, N_points)
            indicator_data = {i: np.zeros(shape) for i in indicator_indices}
            for number, i in enumerate(index, number + 1): 
                evaluate(samples[i])
                for j, x in enumerate(coordinate):
                    f_coordinate(x)
                    for key, indicator in zip(indicator_data, self.indicators):
                        data = indicator_data[key]
                        try:
                            data[i, j] = indicator()
                        except:
                            data[i, j] = None
        else:
            if f_evaluate is None: f_evaluate = self.evaluate
            
            # Initialize timer
            if notify_coordinate:
                from biosteam.utils import Timer
                timer = Timer()
                timer.start()
                def evaluate(**kwargs):
                    f_evaluate(**kwargs)
                    print(f"[Coordinate {n}] Elapsed time: {timer.elapsed_time:.0f} sec")
            else:
                evaluate = f_evaluate
            indicator_data = None
            for n, x in enumerate(coordinate):
                f_coordinate(*x) if multi_coordinate else f_coordinate(x)
                evaluate(notify=notify)
                if indicator_data is None:
                    # Initialize data containers dynamically in case samples are loaded during evaluation
                    N_samples, _ = self.table.shape
                    indicator_indices = var_indices(self.indicators)
                    shape = (N_samples, N_points)
                    indicator_data = {i: np.zeros(shape) for i in indicator_indices}
                for indicator in indicator_data:
                    indicator_data[indicator][:, n] = self.table[indicator]
        
        if xlfile:
            if multi_coordinate:
                columns = pd.MultiIndex.from_tuples(coordinate,
                                                    names=name)
            else:
                columns = pd.Index(coordinate, name=name)
            
            # Save data to excel
            data = pd.DataFrame(data=np.zeros([N_samples, N_points]),
                                columns=columns)
            
            with pd.ExcelWriter(xlfile) as writer:
                for indicator in self.indicators:
                    data[:] = indicator_data[indicator.index]
                    data.to_excel(writer, sheet_name=indicator.short_description)
        return indicator_data
    
    def spearman(self, parameters=None, indicators=None):
        warn(DeprecationWarning('this method will be deprecated in biosteam 2.25; '
                                'use spearman_r instead'), stacklevel=2)
        return self.spearman_r(parameters, indicators)[0]
    
    def spearman_r(self, parameters=None, indicators=None, filter=None, **kwargs): # pragma: no cover
        """
        Return two DataFrame objects of Spearman's rho and p-values between indicators 
        and parameters.
        
        Parameters
        ----------
        parameters : Iterable[Parameter], optional
            Parameters to be correlated with indicators. Defaults to all parameters.
        indicators : Iterable[Indicator], optional 
            Indicators to be correlated with parameters. Defaults to all indicators.
        
        Other Parameters
        ----------------
        filter : Callable(x, y) -> x, y, or string, optional
            Function that accepts 1d arrays of x and y values and returns 
            filtered x and y values to correlate. May also
            be one of the following strings:
                
            * 'none': no filter
            
            * 'omit nan': all NaN values are ignored in correlation
            
            * 'propagate nan': NaN values return NaN correlation results
            
            * 'raise nan': NaN values will raise a ValueError
            
        **kwargs :
            Keyword arguments passed to :func:`scipy.stats.spearmanr`.
        
        See Also
        --------
        :func:`scipy.stats.spearmanr`
        
        """
        from scipy.stats import spearmanr
        return self._correlation(spearmanr, parameters, indicators, filter, kwargs)
    
    def pearson_r(self, parameters=None, indicators=None, filter=None, **kwargs):
        """
        Return two DataFrame objects of Pearson's rho and p-values between indicators 
        and parameters.
        
        Parameters
        ----------
        parameters : Iterable[Parameter], optional
            Parameters to be correlated with indicators. Defaults to all parameters.
        indicators : Iterable[Indicator], optional 
            Indicators to be correlated with parameters. Defaults to all indicators.
            
        Other Parameters
        ----------------
        filter : Callable(x, y) -> x, y, or string, optional
            Function that accepts 1d arrays of x and y values and returns 
            filtered x and y values to correlate. May also
            be one of the following strings:
                
            * 'none': no filter
            
            * 'omit nan': all NaN values are ignored in correlation
            
            * 'propagate nan': NaN values return NaN correlation results
            
            * 'raise nan': NaN values will raise a ValueError
            
        **kwargs :
            Keyword arguments passed to :func:`scipy.stats.pearsonr`.
        
        See Also
        --------
        :func:`scipy.stats.pearsonr`
        
        """
        from scipy.stats import pearsonr
        return self._correlation(pearsonr, parameters, indicators, filter, kwargs)
    
    def kendall_tau(self, parameters=None, indicators=None, filter=None, **kwargs):
        """
        Return two DataFrame objects of Kendall's tau and p-values between indicators 
        and parameters.
        
        Parameters
        ----------
        parameters : Iterable[Parameter], optional
            Parameters to be correlated with indicators. Defaults to all parameters.
        indicators : Iterable[Indicator], optional 
            Indicators to be correlated with parameters. Defaults to all indicators.
            
        Other Parameters
        ----------------
        filter : Callable(x, y) -> x, y, or string, optional
            Function that accepts 1d arrays of x and y values and returns 
            filtered x and y values to correlate. May also
            be one of the following strings:
            
            * 'none': no filter
            
            * 'omit nan': all NaN values are ignored in correlation
            
            * 'propagate nan': NaN values return NaN correlation results
            
            * 'raise nan': NaN values will raise a ValueError
        
        **kwargs :
            Keyword arguments passed to :func:`scipy.stats.kendalltau`.
        
        See Also
        --------
        :func:`scipy.stats.kendalltau`
        
        """
        from scipy.stats import kendalltau
        return self._correlation(kendalltau, parameters, indicators, filter, kwargs)
    
    def kolmogorov_smirnov_d(self, parameters=None, indicators=None, thresholds=[],
                             filter=None, **kwargs):
        """
        Return two DataFrame objects of Kolmogorov–Smirnov's D and p-values
        with the given thresholds for the indicators.
        
        For one particular parameter, all of the sampled values will be divided into two sets,
        one where the resulting indicator value is smaller than or equal to the threshold,
        and the other where the resulting value is larger than the threshold.
        
        Kolmogorov–Smirnov test will then be performed for these two sets of values
        for the particular parameter.
        
        Parameters
        ----------
        parameters : Iterable[Parameter], optional
            Parameters to be correlated with indicators. Defaults to all parameters.
        indicators : Iterable[Indicator], optional 
            Indicators to be correlated with parameters. Defaults to all indicators.
        thresholds : Iterable[float]
            The thresholds for separating parameters into sets.
            
        Other Parameters
        ----------------
        filter : Callable(x, y) -> x, y, or string, optional
            Function that accepts 1d arrays of x and y values and returns 
            filtered x and y values to correlate. May also
            be one of the following strings:
                
            * 'none': no filter
            
            * 'omit nan': all NaN values are ignored in correlation
            
            * 'propagate nan': NaN values return NaN correlation results
            
            * 'raise nan': NaN values will raise a ValueError
            
        **kwargs :
            Keyword arguments passed to :func:`scipy.stats.kstest`.
        
        See Also
        --------
        :func:`scipy.stats.kstest`
        
        """
        from scipy.stats import kstest
        indicators = indicators or self.indicators
        if len(thresholds) != len(indicators):
            raise ValueError(f'The number of indicators {len(indicators)} must match '
                             f'the number of thresholds ({len(thresholds)}).')
        kwargs['thresholds'] = thresholds
        return self._correlation(kstest, parameters, indicators, filter, kwargs)
    
    def _correlation(self, f, parameters, indicators, filter, kwargs):
        """
        Return two DataFrame objects of statistics and p-values between indicators 
        and parameters.
        
        Parameters
        ----------
        f : Callable
            Function with signature f(x, y) -> stat, p
        parameters : Iterable[Parameter], optional
            Parameters to be correlated with indicators. Defaults to all parameters.
        indicators : Iterable[Indicator], optional 
            Indicators to be correlated with parameters. Defaults to all indicators.
        filter : Callable or string, optional
            Function that accepts 1d arrays of x and y values and returns 
            filtered x and y values to correlate. May also
            be one of the following strings:
                
            * 'none': no filter
            
            * 'omit nan': all NaN values are ignored in correlation
            
            * 'propagate nan': NaN values return NaN correlation results
            
            * 'raise nan': NaN values will raise a ValueError
            
        kwargs : dict
            Keyword arguments passed to `f`.
            
        """
        if not parameters: parameters = self._parameters
        table = self.table
        values = table.values.transpose()
        index = table.columns.get_loc
        parameter_indices = var_indices(parameters)
        parameter_data = [values[index(i)] for i in parameter_indices]
        indicator_indices = var_indices(indicators or self.indicators)
        indicator_data = [values[index(i)] for i in indicator_indices]
                
        if not filter: filter = 'propagate nan'
        if isinstance(filter, str):
            name = filter.lower()
            if name == 'omit nan':
                def corr(x, y):
                    index = ~(np.isnan(x) | np.isnan(y))
                    return f(x[index], y[index], **kwargs)
            elif name == 'propagate nan':
                def corr(x, y):
                    if np.isnan(x).any() or np.isnan(y).any(): 
                        NaN = float('NaN')
                        return NaN, NaN
                    else:
                        return f(x, y, **kwargs)
            elif name == 'raise nan':
                def corr(x, y):
                    if np.isnan(x).any() or np.isnan(y).any():
                        raise ValueError('table entries contain NaN values')
                    return f(x, y)
            elif name == 'none':
                corr = lambda x, y: f(x, y, **kwargs)
            else: #: pragma: no cover
                raise ValueError(
                    f"invalid filter '{filter}'; valid filter names are: "
                    "'omit nan', 'propagate nan', 'raise nan', and 'none'"
                )
        elif callable(filter): 
            corr = lambda x, y: f(*filter(x, y), **kwargs)
        else: #: pragma: no cover
            raise TypeError("filter must be either a string or a callable; "
                            "not a '{type(filter).__name__}' object")

        if 'thresholds' not in kwargs:
            data = np.array([[corr(p.astype(float), m.astype(float)) for m in indicator_data] for p in parameter_data])
        else: # KS test
            thresholds = kwargs.pop('thresholds')
            data = np.array(
                [
                    [corr(parameter_data[n_p][indicator_data[n_m]<=thresholds[n_m]],
                          parameter_data[n_p][indicator_data[n_m]>thresholds[n_m]]) \
                     for n_m in range(len(indicator_data))] \
                for n_p in range(len(parameter_data))])
        
        index = indices_to_multiindex(parameter_indices, ('Element', 'Parameter'))
        columns = indices_to_multiindex(indicator_indices, ('Element', 'Indicator'))        
        return [pd.DataFrame(i, index=index, columns=columns) for i in (data[..., 0], data[..., 1])]
        
    def create_fitted_model(self, parameters, indicators): # pragma: no cover
        from pipeml import FittedModel
        Xdf = self.table[[i.index for i in parameters]]
        ydf_index = indicators.index if isinstance(indicators, Indicator) else [i.index for i in indicators]
        ydf = self.table[ydf_index]
        return FittedModel.from_dfs(Xdf, ydf)
    
    def __call__(self, sample, **kwargs):
        """Return pandas Series of indicator values at given sample."""
        self._update_state(np.asarray(sample, dtype=float), **kwargs)
        return pd.Series({i.index: i() for i in self._indicators})
    
    def _repr(self, m):
        clsname = type(self).__name__
        newline = "\n" + " "*(len(clsname)+2)
        return f'{clsname}: {newline.join([i.describe() for i in self.indicators])}'
    
    def __repr__(self):
        return f'<{type(self).__name__}: {len(self.parameters)}-parameters, {len(self.indicators)}-indicators>'
    
    def _info(self, p, m):
        info = f'{type(self).__name__}:'
        parameters = self._parameters
        if parameters: 
            p_max = len(parameters)
            if p is None: p = p_max
            ptitle = 'parameters: '
            info += '\n' + ptitle
            newline = "\n" + " "*(len(ptitle))
            info += newline.join([
                parameters[i].describe(
                    distribution=False, bounds=False
                ) for i in range(p)
            ])
            if p < p_max: info += newline + "..." 
        else:
            info += '\nparameters: None'
        indicators = self._indicators
        if indicators: 
            m_max = len(indicators)
            if m is None: m = m_max
            mtitle = 'indicators: '
            info += '\n' + mtitle
            newline = "\n" + " "*(len(mtitle))
            info += newline.join([
                indicators[i].describe() for i in range(m)
            ])
            if m < m_max: info += newline + "..." 
        else:
            info += '\nindicators: None'
        return info
    
    def show(self, p=None, m=None):
        """Return information on p-parameters and m-indicators."""
        print(self._info(p, m))
    _ipython_display_ = show
    
EasyInputModel = Model
Model.load_parameter_distributions = Model.parameters_from_df
