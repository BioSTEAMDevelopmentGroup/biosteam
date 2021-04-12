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
from biosteam import speed_up
import chaospy as cp
import numpy as np
import pandas as pd

__all__ = ('AgileModel',)


# %% Agile model


class AgileModel:
    
    __slots__ = (
        'table',          # [DataFrame] All arguments and results.
        '_agile_scenario',  # [ScenarioPad] Compiles system scenarios.
        '_parameters',    # list[Parameter] All parameters.
        '_metrics',       # tuple[Metric] Metrics to be evaluated by model.
        '_samples',       # [array] Argument sample space.
    )
    
    def __init__(self, agile_scenario):
        self._agile_scenario = agile_scenario
        self._parameters = {}
        self._metrics = {}
    
    @property
    def agile_scenario(self):
        return self._agile_scenario
    
    def _erase(self):
        """Erase cached data."""
        self._samples = None
    
    def parameter(self, setter=None, system=None, element=None, kind='isolated',
                  name=None, distribution=None, units=None, baseline=None, bounds=None):
        """
        Define and register parameter.
        
        Parameters
        ----------    
        setter : function
                 Should set parameter in the element.
        system : System, optional
            System which the parameter pertains to.
        element : Unit or :class:`~thermosteam.Stream`
                  Element in the system being altered.
        kind : {'coupled', 'isolated', 'design', 'cost'}
            * 'coupled': parameter is coupled to the system.
            * 'isolated': parameter does not affect the system but does affect the element (if any).
            * 'design': parameter only affects design and/or cost of the element.
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
            return lambda setter: self.parameter(setter, system, element, kind,
                                                 name, distribution, units, 
                                                 baseline, bounds)
        p = Parameter(name, setter, element, system, distribution,
                      units, baseline, bounds, kind)
        parameters = self._parameters
        all_parameters = sum([list(i) for i in parameters.values()], [])
        Parameter.check_index_unique(p, all_parameters)
        if system in parameters:
            system_parameters = parameters[system]
        else:
            parameters[system] = system_parameters = []
        system_parameters.append(p)
        self._erase()
        return p
    
    def metric(self, getter=None, name=None, units=None, system=None, element='Biorefinery'):
        """
        Define and register metric.
        
        Parameters
        ----------    
        getter : function, optional
                 Should return metric.
        name : str, optional
               Name of metric. If None, defaults to the name of the getter.
        units : str, optional
                Metric units of measure
        system : System, optional
            System which the metric pertains to.
        element : object, optional
                  Element being evaluated. Works mainly for bookkeeping. 
                  Defaults to 'Biorefinery'.
        
        Notes
        -----
        This method works as a decorator.
        
        """
        if not getter: return lambda getter: self.metric(getter, name, units, system, element)
        if not name and hasattr(getter, '__name__'):
            name = format_title(getter.__name__)
        metric = Metric(name, getter, units, element)
        metrics = self._metrics
        all_metrics = sum([list(i) for i in metrics.values()], [])
        Metric.check_index_unique(metric, all_metrics)
        if system in metrics:
            system_metrics = metrics[system]
        else:
            metrics[system] = system_metrics = []
        system_metrics.append(metric)
        return metric 
    
    sample = Model.sample
    get_joint_distribution = Model.get_joint_distribution
    load_samples = Model.load_samples

    def _load_sample_order(self, samples, parameters):
        key = lambda x: samples[x, i]
        index = list(range(samples.shape[0]))
        for i in range(self._N_parameters_cache-1,  -1, -1):
            index.sort(key=key)
        self._index = index
        
    def load_samples(self, samples):
        """Load samples for evaluation
        
        Parameters
        ----------
        samples : numpy.ndarray, dim=2
        
        """
        self._load_parameters()
        parameters = self._parameters
        Parameter.check_indices_unique(parameters)
        N_parameters = self._N_parameters_cache
        if not isinstance(samples, np.ndarray):
            raise TypeError(f'samples must be an ndarray, not a {type(samples).__name__} object')
        if samples.ndim == 1:
            samples = samples[:, np.newaxis]
        elif samples.ndim != 2:
            raise ValueError('samples must be 2 dimensional')
        if samples.shape[1] != N_parameters:
            raise ValueError(f'number of parameters in samples ({samples.shape[1]}) must be equal to the number of parameters ({len(N_parameters)})')
        metrics = self._metrics
        self._load_sample_order(samples, parameters)
        empty_metric_data = np.zeros((len(samples), len(metrics)))
        self.table = pd.DataFrame(np.hstack((samples, empty_metric_data)),
                                  columns=var_columns(parameters + metrics))
        self._samples = samples
        
    def evaluate(self, thorough=True, jit=True, notify=False):
        """
        Evaluate metrics over the loaded samples and save values to `table`.
        
        Parameters
        ----------
        thorough : bool
            If True, simulate the whole system with each sample.
            If False, simulate only the affected parts of the system and
            skip simulation for repeated states.
        jit : bool
            Whether to JIT compile functions with Numba to speed up simulation.
        notify=False : bool, optional
            If True, notify elapsed time after each sample evaluation. 
        
        """
        if jit: speed_up()
        samples = self._samples
        if samples is None: raise RuntimeError('must load samples before evaluating')
        evaluate_sample = self._evaluate_sample
        table = self.table
        if notify:
            from biosteam.utils import TicToc
            timer = TicToc()
            timer.tic()
            def evaluate(sample, thorough, count=[0]):
                count[0] += 1
                values = evaluate_sample(sample, thorough)
                print(f"{count} Elapsed time: {timer.elapsed_time:.0f} sec")
                return values
        else:
            evaluate = evaluate_sample
        index = self._index
        values = [None] * len(index)
        for i in index:
            values[i] = evaluate(samples[i], thorough)
        table[var_indices(self.metrics)] = values
        
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
    
    problem = Model.problem
    spearman_r = Model.spearman_r
    pearson_r = Model.pearson_r
    kendall_tau = Model. kendall_tau
    kolmogorov_smirnov_d = Model.kolmogorov_smirnov_d
    _correlation = Model._correlation
    _repr = Model._repr    
    __repr__ = Model.__repr__
    show = Model.show
    _ipython_display_ = Model._ipython_display_
    