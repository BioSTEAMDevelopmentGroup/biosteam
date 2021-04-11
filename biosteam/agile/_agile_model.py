# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from ..evaluation._state import State
from ..evaluation._metric import Metric
from ..evaluation._parameter import Parameter
from ..evaluation._model import Model
from ..evaluation._utils import var_indices, var_columns, indices_to_multiindex
from ..utils import format_title
import chaospy as cp
import numpy as np
import pandas as pd

__all__ = ('AgileModel',)


# %% Agile model


class AgileModel:
    
    __slots__ = (
        'table',          # [DataFrame] All arguments and results.
        '_models',        # list[Model] Models for each system scenario.
        '_index',          # list[int] Order of sample evaluation for performance.
        '_agile_scenario',  # [ScenarioPad] Compiles system scenarios.
        '_parameters',    # list[Parameter] All parameters.
        '_metrics',       # tuple[Metric] Metrics to be evaluated by model.
        '_samples',       # [array] Argument sample space.
    )
    
    def __init__(self, models, agile_scenario, parameters=None, metrics=None):
        self._models = tuple(models)
        self._agile_scenario = agile_scenario
        self._parameters = parameters or []
        self._metrics = metrics or []
    
    @property
    def parameters(self):
        return tuple(self._parameters)
    
    @property
    def models(self):
        return self._models
    
    @property
    def agile_scenario(self):
        return self._agile_scenario
    
    @property
    def system(self):
        return None
    
    def get_all_metrics(self):
        return sum([i._metrics for i in self._models], self._metrics)
    
    def get_all_parameters(self):
        models = self._models
        for i in models: i._load_parameters()
        return sum([i._parameters for i in models], self._parameters)
    
    def get_full_samples(self):
        models = self._models
        samples = [i._samples for i in models]
        samples.append(self._samples)
        return np.vstack(samples)
    
    def _erase(self):
        """Erase cached data."""
        self._samples = None
    
    metrics = Model.metrics
    metric = Model.metric
    parameter = Model.parameter
    sample = Model.sample
    _sorted_samples = Model._sorted_samples
    
    def get_joint_distribution(self):
        """
        Return a chaospy joint distribution object constituting of all
        parameter objects.
        
        """
        return cp.distributions.J(*[i.distribution for i in self.get_all_parameters()])

    def load_samples(self, samples):
        """
        Load samples for evaluation.
        
        Parameters
        ----------
        samples : numpy.ndarray, dim=2
        
        """
        metrics = self.get_all_metrics()
        parameters = self.get_all_parameters()
        Parameter.check_indices_unique(parameters)
        N_parameters = len(parameters)
        samples = self.get_full_samples()
        if not isinstance(samples, np.ndarray):
            raise TypeError(f'samples must be an ndarray, not a {type(samples).__name__} object')
        if samples.ndim == 1:
            samples = samples[:, np.newaxis]
        elif samples.ndim != 2:
            raise ValueError('samples must be 2 dimensional')
        if samples.shape[1] != N_parameters:
            raise ValueError(f'number of parameters in samples ({samples.shape[1]}) must be equal to the number of parameters ({len(N_parameters)})')
        samples = self._sorted_samples(samples, parameters)
        empty_metric_data = np.zeros((len(samples), len(metrics)))
        self.table = pd.DataFrame(np.hstack((samples, empty_metric_data)),
                                  columns=var_columns(parameters + metrics))
        self._samples = samples
        
    def evaluate(self, thorough=True, jit=True, notify=False):
        """
        Evaluate metrics over the loaded samples and save values to `table`.
        
        Parameters
        ----------
        jit : bool
            Whether to JIT compile functions with Numba to speed up simulation.
        notify=False : bool, optional
            If True, notify elapsed time after each sample evaluation. 
        
        """
        samples = self._samples
        if samples is None: raise RuntimeError('must load samples before evaluating')
        evaluate_sample = self._evaluate_sample
        table = self.table
        if notify:
            from biosteam.utils import TicToc
            timer = TicToc()
            timer.tic()
            def evaluate(sample, count=[0]):
                count[0] += 1
                values = evaluate_sample(sample, thorough)
                print(f"{count} Elapsed time: {timer.elapsed_time:.0f} sec")
                return values
        else:
            evaluate = evaluate_sample
        parameters = self.get_all_parameters()
        metrics = self.get_all_metrics()
        table[var_indices(metrics)] = [evaluate(samples[i], parameters, metrics) 
                                       for i in enumerate(self._index)]
    
    def _evaluate_sample(self, sample, parameters, metrics):
        Nf = len(self._parameters)
        for f, s in zip(parameters[:Nf], sample[:Nf]): f.setter(s)
        agile_scenario = self._agile_scenario
        scenarios = []
        metrics = []
        for i in self._models:
            N = Nf
            Nf = N + len(i._parameters)
            for f, s in zip(parameters[N:Nf], sample[N:Nf]): f.setter(s)
            i.simulate()
            scenario = agile_scenario.create_scenario(i.system)
            scenarios.append(scenario)
            metrics.extend([i() for i in i._metrics])
        metrics.extend([i() for i in self._metrics])
        agile_scenario.compile_scenarios(scenarios)
        return metrics
    
    def problem(self):
        """
        Return a dictionary of parameter metadata (referred to as "problem") 
        to be used for sampling by ``SALib``.
        
        See Also
        --------
        `SALib basics <https://salib.readthedocs.io/en/latest/basics.html#an-example>`_
        
        """
        params = self.get_all_parameters()
        return {
            'num_vars': len(params),
            'names': [i.name for i in params],
            'bounds': [i.bounds if i.bounds
                       else (i.distribution.lower[0], i.distribution.upper[0])
                       for i in params]
        }
    
    spearman_r = Model.spearman_r
    pearson_r = Model.pearson_r
    kendall_tau = Model. kendall_tau
    kolmogorov_smirnov_d = Model.kolmogorov_smirnov_d
    _correlation = Model._correlation
    _repr = Model._repr    
    __repr__ = Model.__repr__
    show = Model.show
    _ipython_display_ = Model._ipython_display_
    