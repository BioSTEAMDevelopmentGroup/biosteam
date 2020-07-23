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
from pipeml import FittedModel
from biosteam.utils import TicToc
from scipy.stats import spearmanr
from biosteam import speed_up

__all__ = ('Model',)

# %% Functions

var_indices = lambda vars: [var.index for var in vars]
var_columns = lambda vars: pd.MultiIndex.from_tuples(
                                 var_indices(vars),
                                 names=('Element', 'Variable'))
            

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
    :doc:`../tutorial/Monte_Carlo`
    
    """
    __slots__ = ('_table',   # [DataFrame] All arguments and results
                 '_metrics', # tuple[Metric] Metrics to be evaluated by model
                 '_index',   # list[int] Order of sample evaluation for performance
                 '_samples', # [array] Argument sample space
                 '_setters', # list[function] Cached parameter setters
                 '_getters', # list[function] Cached metric getters
                 '_failed_metrics', # list[np.nan] Cached NaN values for failed evaluations
                 '_metric_indices', # list[Hashable] Cached metric indices.
    )
    def __init__(self, system, metrics, specification=None, skip=False, params=None):
        super().__init__(system, specification, skip, params)
        self.metrics = metrics
        self._samples = self._table = None
        
    def copy(self):
        """Return copy."""
        copy = super().copy()
        copy._metrics = self._metrics
        if self._update:
            copy._table = self._table.copy()
        else:
            copy._samples = copy._table = None
        return copy
    
    def _erase(self):
        """Erase cached data."""
        self._update = self._table = self._samples = None
    
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
        """[DataFrame] Table of the sample space with results in the final column."""
        return self._table
    
    def load_samples(self, samples):
        """Load samples for evaluation
        
        Parameters
        ----------
        samples : numpy.ndarray, dim=2
        
        """
        if not self._update: self._load_params()
        params = self._params
        paramlen = len(params)
        if not isinstance(samples, np.ndarray):
            raise TypeError(f'samples must be an ndarray, not a {type(samples).__name__} object')
        if samples.ndim == 1:
            samples = samples[:, np.newaxis]
        elif samples.ndim != 2:
            raise ValueError('samples must be 2 dimensional')
        if samples.shape[1] != paramlen:
            raise ValueError(f'number of parameters in samples ({samples.shape[1]}) must be equal to the number of parameters ({len(params)})')
        key = lambda x: samples[x][i]
        N_samples = len(samples)
        index = list(range(N_samples))
        for i in range(paramlen-1,  -1, -1):
            if not params[i].system: break
            index.sort(key=key)
        self._index = index
        metrics = self._metrics
        N_metrics = len(metrics)
        empty_metric_data = np.zeros((N_samples, N_metrics))
        self._table = pd.DataFrame(np.hstack((samples, empty_metric_data)),
                                   columns=var_columns(tuple(params) + metrics))
        self._samples = samples
        self._setters = [i.setter for i in params]
        self._getters = [i.getter for i in metrics]
        self._failed_metrics = N_metrics * [np.nan]
        self._metric_indices = var_indices(metrics)
        
    def evaluate(self, thorough=True):
        """
        Evaluate metric over the argument sample space and save values to ``table``.
        
        Parameters
        ----------
        thorough : bool
            If True, simulate the whole system with each sample.
            If False, simulate only the affected parts of the system.
        
        """
        speed_up()
        samples = self._samples
        if samples is None:
            raise ValueError('must load samples or distribution before evaluating')
        evaluate_sample = self._evaluate_sample_thorough if thorough else self._evaluate_sample_smart
        self.table[self._metric_indices] = [evaluate_sample(samples[i]) for i in self._index]
    
    def _evaluate_sample_thorough(self, sample):
        for f, s in zip(self._setters, sample): f(s)
        if self._specification: self._specification()
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
            return self._failed_evaluation()
    
    def _failed_evaluation(self):
        self._system.empty_recycles()
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
        table = self.table
        N_samples, N_parameters = table.shape
        N_points = len(coordinate)
        
        # Initialize data containers
        metric_indices = self._metric_indices
        metric_data = {i: np.zeros([N_samples, N_points]) for i in metric_indices}
        
        # Initialize timer
        if notify:
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
                    data.to_excel(writer, sheet_name=i.name)
        return metric_data
    
    def spearman(self, metrics=()):
        """
        Return DataFrame of Spearman's rank correlation for metrics vs parameters.
        
        Parameters
        ----------
        metrics=() : Iterable[Metric], defaults to all metrics
            Metrics to be correlated with parameters.
        
        """
        data = self.table
        param_cols = list(data)
        all_metric_names = var_indices(self.metrics)
        if not metrics: 
            metric_names = all_metric_names
        else:
            metric_names = var_indices(metrics)
        params = list(self._params)
        
        for i in all_metric_names: param_cols.remove(i)
        param_descriptions = [i.describe() for i in params]
        allrhos = []
        for name in metric_names:
            rhos = []
            metric = data[name]
            for col in param_cols:
                rho, p = spearmanr(data[col], metric)
                rhos.append(rho)
            allrhos.append(rhos)
        allrhos = np.array(allrhos).transpose()
        return pd.DataFrame(allrhos, index=param_descriptions,
                            columns=metric_names)
    
    def create_fitted_model(self, parameters, metrics):
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
    
        