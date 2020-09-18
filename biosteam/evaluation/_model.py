# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import numpy as np
import pandas as pd
from ._state import State
from ._metric import Metric
from biosteam import speed_up

__all__ = ('Model',)

# %% Functions

var_indices = lambda vars: [var.index for var in vars]
var_columns = lambda vars, names=None: indices_to_multiindex(var_indices(vars), names)
indices_to_multiindex = lambda indices, names=None: pd.MultiIndex.from_tuples(
                                            indices,
                                            names=names or ('Element', 'Variable'),
                                        )

# %% Grid of simulation blocks

class Model(State):
    """
    Create a Model object that allows for evaluation over a sample space.
    
    Parameters
    ----------
    system : System
        Should reflect the model state.
    metrics : tuple[Metric]
        Metrics to be evaluated by model.
    specification=None : Function, optional
        Loads speficications once all parameters are set.
    skip=False : bool, optional
        If True, skip simulation for repeated states.
    params=None : Iterable[Parameter], optional
        Parameters to sample from.
    
    Examples
    --------
    For a complete example, visit :doc:`../tutorial/Monte_Carlo`. An abridged example is presented here:
    
    >>> from chaospy import distributions as shape
    >>> from biorefineries import lipidcane as lc
    >>> import biosteam as bst
    >>> solve_IRR = lc.lipidcane_tea.solve_IRR
    >>> total_utility_cost = lambda: lc.lipidcane_tea.utility_cost / 10**6 # In 10^6 USD/yr
    >>> metrics = (bst.Metric('Internal rate of return', lc.lipidcane_tea.solve_IRR, '%'),
    ...            bst.Metric('Utility cost', total_utility_cost, '10^6 USD/yr'))
    >>> model = bst.Model(lc.lipidcane_sys, metrics)
    
    The Model object begins with no parameters:
    
    >>> model.show()
    Model: Biorefinery internal rate of return [%]
           Biorefinery utility cost [10^6 USD/yr]
     (No parameters)
     
    Add number of fermentation reactors as a "design" parameter (which doesn't affect mass or energy balances'):
        
    >>> R301 = bst.main_flowsheet.unit.R301 # The Fermentation Unit
    >>> @model.parameter(name='Number of reactors',
    ...                  element=R301, kind='design',
    ...                  distribution=shape.Uniform(4, 10))
    ... def set_N_reactors(N):
    ...     R301.N = round(N)
    
    The decorator uses the function to create a Parameter object and add it to the model:

    >>> parameters = model.get_parameters()
    >>> parameters
    (<Parameter: [Fermentation-R301] Number of reactors>,)

    Calling a Parameter object will update the parameter and results:

    >>> set_N_reactors_parameter = parameters[0]
    >>> set_N_reactors_parameter(5)
    >>> R301.purchase_cost / 1e6
    1.7
    >>> set_N_reactors_parameter(8)
    >>> R301.purchase_cost / 1e6
    2.0
    
    Add the fermentation unit base cost as a "cost" parameter with a triangular distribution (which doesn't affect mass and energy balances nor design requirements'):
    
    >>> reactors_cost_coefficients = R301.cost_items['Reactors']
    >>> mid = reactors_cost_coefficients.n # Most probable at baseline value
    >>> lb = mid - 0.1 # Minimum
    >>> ub = mid + 0.1 # Maximum
    >>> @model.parameter(element=R301, kind='cost',
    ...              distribution=shape.Triangle(lb, mid, ub))
    ... def set_exponential_cost_coefficient(exponential_cost_coefficient):
    ...     reactors_cost_coefficients.n = exponential_cost_coefficient
    
    Note that if the name was not defined, it defaults to the setter's signature:
    
    >>> model.get_parameters()
    (<Parameter: [Fermentation-R301] Number of reactors>,
     <Parameter: [Fermentation-R301] Exponential cost coefficient>)

    Add feedstock price as a "isolated" parameter (which doesn't affect unit operations in any way):
                                                   
    >>> lipidcane = lc.lipidcane # The feedstock stream
    >>> lb = lipidcane.price * 0.9 # Minimum price
    >>> ub = lipidcane.price * 1.1 # Maximum price
    >>> @model.parameter(element=lipidcane, kind='isolated', units='USD/kg',
    ...                  distribution=shape.Uniform(lb, ub))
    ... def set_feed_price(feedstock_price):
    ...     lipidcane.price = feedstock_price                                   
    
    Add lipid fraction as a "coupled" parameter (which affects mass and energy balances):   
    
    >>> from biorefineries.lipidcane.utils import set_lipid_fraction
    >>> # Note that if the setter function is already made,
    >>> # you can pass it as the first argument
    >>> set_lipid_fraction = model.parameter(set_lipid_fraction,
    ...                                      element=lipidcane, kind='coupled',
    ...                                      distribution=shape.Uniform(0.05, 0.10))
    
    Add fermentation efficiency as a "coupled" parameter:

    >>> @model.parameter(element=R301, kind='coupled',
    ...                  distribution=shape.Triangle(0.85, 0.90, 0.95))
    ... def set_fermentation_efficiency(efficiency):
    ...     R301.efficiency = efficiency

    Note that all parameters are stored in the model with highly coupled parameters first:

    >>> model.show()
    Model: Biorefinery internal rate of return [%]
           Biorefinery utility cost [10^6 USD/yr]
     Element:           Parameter:
     Stream-lipidcane   Lipid fraction
     Fermentation-R301  Efficiency
                        Number of reactors
                        Exponential cost coefficient
     Stream-lipidcane   Feedstock price        

    Get dictionary that contain DataFrame objects of parameter distributions:

    >>> df_dct = model.get_distribution_summary()
    >>> df_dct['Uniform']    
                 Element                Name   Units    Shape  lower  upper
    0   Stream-lipidcane      Lipid fraction          Uniform   0.05    0.1
    1  Fermentation-R301  Number of reactors          Uniform      4     10
    2   Stream-lipidcane     Feedstock price  USD/kg  Uniform 0.0311  0.038
    
    Evaluate sample:
        
    >>> model([0.05, 0.85, 8, 100000, 0.040]) # Returns metrics (IRR and utility cost)
    Biorefinery  Internal rate of return [%]   0.102
                 Utility cost [10^6 USD/yr]      -18
    dtype: float64
    
    Sample from a joint distribution, and simulate samples:

    >>> import numpy as np
    >>> np.random.seed(1234) # For consistent results
    >>> N_samples = 10
    >>> rule = 'L' # For Latin-Hypercube sampling
    >>> samples = model.sample(N_samples, rule)
    >>> model.load_samples(samples)
    >>> model.evaluate(jit=False)
    >>> table = model.table # All evaluations are stored as a pandas DataFrame
    >>> table['Biorefinery'] # Only biorefinery metrics
    Variable  Internal rate of return [%]  Utility cost [10^6 USD/yr]
    0                                0.16                       -18.8
    1                               0.135                       -17.9
    2                               0.181                       -25.5
    3                               0.149                       -19.8
    4                               0.162                       -25.1
    5                               0.168                       -22.2
    6                               0.145                       -20.5
    7                               0.177                       -21.6
    8                                0.18                       -24.6
    9                                0.17                       -23.5
    
    Note that coupled parameters are on the left most columns, and are ordered 
    from upstream to downstream (e.g. <Stream: Lipid cane> is upstream from <Fermentation: R301>):

    >>> # Reset for future tests
    >>> bst.process_tools.default_utilities()
    >>> bst.CE = 567.5

    """
    __slots__ = ('_user_table',     # [DataFrame] All arguments and results (made by user).
                 '_table',          # [DataFrame] All arguments and results.
                 '_metrics',        # tuple[Metric] Metrics to be evaluated by model.
                 '_index',          # list[int] Order of sample evaluation for performance.
                 '_samples',        # [array] Argument sample space.
                 '_setters',        # list[function] Cached parameter setters.
                 '_getters',        # list[function] Cached metric getters.
                 '_failed_metrics', # list[np.nan] Cached NaN values for failed evaluations.
                 '_metric_indices', # list[Hashable] Cached metric indices.
    )
    def __init__(self, system, metrics, specification=None, skip=False, parameters=None):
        super().__init__(system, specification, skip, parameters)
        self.metrics = metrics
        self._samples = self._table = self._user_table = None
        
    def copy(self):
        """Return copy."""
        copy = super().copy()
        copy._metrics = self._metrics
        if self._update:
            copy._table = self._table.copy()
            copy._user_table = self._user_table.copy() if self._user_table else None
        else:
            copy._samples = copy._table = self._user_table = None
        return copy
    
    def _erase(self):
        """Erase cached data."""
        self._update = self._samples = None
    
    @property
    def metrics(self):
        """tuple[Metric] Metrics to be evaluated by model."""
        return self._metrics
    @metrics.setter
    def metrics(self, metrics):
        for i in metrics:
            if not isinstance(i, Metric):
                raise ValueError(f"metrics must be '{Metric.__name__}' "
                                 f"objects, not '{type(i).__name__}'")
        self._metrics = tuple(metrics)
    
    @property
    def table(self):
        """[DataFrame] Table of the sample space and results."""
        return self._table if self._user_table is None else self._user_table
    @table.setter
    def table(self, table):
        self._user_table = table
    
    def load_samples(self, samples):
        """Load samples for evaluation
        
        Parameters
        ----------
        samples : numpy.ndarray, dim=2
        
        """
        if not self._update: self._load_parameters()
        parameters = self._parameters
        N_parameters = len(parameters)
        if not isinstance(samples, np.ndarray):
            raise TypeError(f'samples must be an ndarray, not a {type(samples).__name__} object')
        if samples.ndim == 1:
            samples = samples[:, np.newaxis]
        elif samples.ndim != 2:
            raise ValueError('samples must be 2 dimensional')
        if samples.shape[1] != N_parameters:
            raise ValueError(f'number of parameters in samples ({samples.shape[1]}) must be equal to the number of parameters ({len(N_parameters)})')
        key = lambda x: samples[x][i]
        N_samples = len(samples)
        index = list(range(N_samples))
        for i in range(N_parameters-1,  -1, -1):
            if not parameters[i].system: break
            index.sort(key=key)
        self._index = index
        metrics = self._metrics
        N_metrics = len(metrics)
        empty_metric_data = np.zeros((N_samples, N_metrics))
        self._table = pd.DataFrame(np.hstack((samples, empty_metric_data)),
                                   columns=var_columns(tuple(parameters) + metrics))
        self._samples = samples
        self._setters = [i.setter for i in parameters]
        self._getters = [i.getter for i in metrics]
        self._failed_metrics = N_metrics * [np.nan]
        self._metric_indices = var_indices(metrics)
        
    def evaluate(self, thorough=True, jit=True):
        """
        Evaluate metric over the argument sample space and save values to ``table``.
        
        Parameters
        ----------
        thorough : bool
            If True, simulate the whole system with each sample.
            If False, simulate only the affected parts of the system.
        jit : bool
            Whether to JIT compile functions with Numba to speed up simulation.
        
        """
        if jit: speed_up()
        samples = self._samples
        if samples is None: raise RuntimeError('must load samples before evaluating')
        evaluate_sample = self._evaluate_sample_thorough if thorough else self._evaluate_sample_smart
        self._table[self._metric_indices] = [evaluate_sample(samples[i]) for i in self._index]
    
    def _evaluate_sample_thorough(self, sample):
        for f, s in zip(self._setters, sample): f(s)
        if self._specification: self._specification()
        try:
            self._system.simulate()
            return [i() for i in self._getters]
        except:
            self._reset_system()
            try:
                self._system.simulate()
                return [i() for i in self._getters]
            except:
                return self._failed_evaluation()
    
    def _evaluate_sample_smart(self, sample):
        try:
            self._update(sample, self._specification)
            return [i() for i in self._getters]
        except:
            self._reset_system()
            try: 
                self._update(sample, self._specification)
                return [i() for i in self._getters]
            except:
                return self._failed_evaluation()
    
    def _reset_system(self):
        self._system.empty_recycles()
        self._system.reset_cache()
    
    def _failed_evaluation(self):
        self._reset_system()
        return self._failed_metrics
    
    def metrics_at_baseline(self):
        """Return metric values at baseline sample."""
        baseline = self.get_baseline_sample()
        return self(baseline)
    
    def evaluate_across_coordinate(self, name, f_coordinate, coordinate,
                                   *, xlfile=None, notify=True,
                                   multi_coordinate=False):
        """
        Evaluate across coordinate and save sample metrics.
        
        Parameters
        ----------
        name : str or tuple[str]
            Name of coordinate
        f_coordinate : function
            Should change state of system given the coordinate.
        coordinte : array
            Coordinate values.
        xlfile : str, optional
            Name of file to save. File must end with ".xlsx"
        rule='L' : str
            Sampling rule.
        notify=True : bool, optional
            If True, notify elapsed time after each coordinate evaluation.
        
        """
        if self._samples is None: raise RuntimeError('must load samples before evaluating')
        table = self._table
        N_samples, N_parameters = table.shape
        N_points = len(coordinate)
        
        # Initialize data containers
        metric_indices = self._metric_indices
        shape = (N_samples, N_points)
        metric_data = {i: np.zeros(shape) for i in metric_indices}
        
        # Initialize timer
        if notify:
            from biosteam.utils import TicToc
            timer = TicToc()
            timer.tic()
            def evaluate():
                self.evaluate()
                print(f"[{n}] Elapsed time: {timer.elapsed_time:.0f} sec")
        else:
            evaluate = self.evaluate
        
        for n, x in enumerate(coordinate):
            f_coordinate(*x) if multi_coordinate else f_coordinate(x)
            evaluate()
            for metric in metric_data:
                metric_data[metric][:, n] = table[metric]
        
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
                for i, metric in zip(metric_indices, metric_data):
                    data[:] = metric_data[metric]
                    data.to_excel(writer, sheet_name=i[1])
        return metric_data
    
    def spearman(self, parameters=(), metrics=()):
        """
        Return DataFrame of Spearman's rank correlation for metrics vs parameters.
        
        Parameters
        ----------
        parameters=() : Iterable[Parameter], defaults to all parameters
            Parameters to be correlated with metrics .
        metrics=() : Iterable[Metric], defaults to all metrics
            Metrics to be correlated with parameters.
        
        """
        from scipy.stats import spearmanr
        if not parameters: parameters = self._parameters
        table = self.table
        parameter_indices = var_indices(parameters)
        parameter_data = [table[i] for i in parameter_indices]
        metric_indices = var_indices(metrics) if metrics else self._metric_indices 
        metric_data = [table[i] for i in metric_indices]
        rhos = [[spearmanr(p, m)[0] for m in metric_data] for p in parameter_data]
        return pd.DataFrame(rhos, 
                            index=indices_to_multiindex(
                                parameter_indices, 
                                ('Element', 'Parameter')
                            )
               )
    
    def create_fitted_model(self, parameters, metrics):
        from pipeml import FittedModel
        Xdf = self.table[[i.index for i in parameters]]
        ydf_index = metrics.index if isinstance(metrics, Metric) else [i.index for i in metrics]
        ydf = self.table[ydf_index]
        return FittedModel.from_dfs(Xdf, ydf)
    
    def __call__(self, sample):
        """Return pandas Series of metric values at given sample."""
        super().__call__(sample)
        return pd.Series({i.index: i.getter() for i in self._metrics})
    
    def _repr(self):
        clsname = type(self).__name__
        newline = "\n" + " "*(len(clsname)+2)
        return f'{clsname}: {newline.join([i.describe() for i in self.metrics])}'
        
    def __repr__(self):
        return f'<{type(self).__name__}: {", ".join([i.name for i in self.metrics])}>'
    
        