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


def all_variables(vdct):
    return sum([tuple(i) for i in vdct.values()], ())

def get_system_variables(vdct, system):
    if system in vdct:
        system_variables = vdct[system]
    else:
        vdct[system] = system_variables = []
    return system_variables

# %% Agile model


class AgileModel:
    
    __slots__ = (
        'table',              # [DataFrame] All arguments and results.
        '_agile_scenario',    # [ScenarioPad] Compiles system scenarios.
        '_N_parameters_cache', # list[int] Number of parameters cache.
        '_system_parameters', # dict[System, Parameter] All parameters.
        '_system_metrics',    # dict[System, Metric] Metrics to be evaluated by model.
        '_samples',           # [array] Argument sample space.
        '_index',             # list[int] Sample evaluation order.
    )
    
    def __init__(self, agile_scenario):
        self._agile_scenario = agile_scenario
        self._system_parameters = {None: []}
        self._system_metrics = {None: []}
        self._N_parameters_cache = None
    
    @property
    def agile_scenario(self):
        return self._agile_scenario
    
    @property
    def metrics(self):
        metrics = self._system_metrics
        return all_variables(metrics)
    
    @property
    def parameters(self):
        return self.get_parameters()
        
    def get_parameters(self):
        parameters = self._system_parameters
        return all_variables(parameters)
    
    def get_system_metrics(self, system):
        return get_system_variables(self._system_metrics, system)
    
    def get_system_parameters(self, system):
        return get_system_variables(self._system_parameters, system)
    
    def _erase(self):
        """Erase cached data."""
        self._samples = None
    
    def _load_parameters(self):
        """Load parameters."""
        sp = self._system_parameters
        systems = list(sp)
        systems.sort(key=lambda x: len(x.units) if x else 0)
        self._system_parameters = sp = {i: sp[i] for i in systems}
        parameters_lists = tuple(sp.values())
        N_parameters = sum([len(i) for i in parameters_lists])
        if N_parameters != self._N_parameters_cache:
            self._N_parameters_cache = N_parameters
            for sys, parameters in zip(systems, parameters_lists):
                if sys: Parameter.sort_parameters(parameters)
    
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
        Parameter.check_index_unique(p, self.parameters)
        self.get_system_parameters(system).append(p)
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
        Metric.check_index_unique(metric, self.metrics)
        self.get_system_metrics(system).append(metric)
        return metric 
    
    sample = Model.sample

    def _load_sample_order(self, samples, parameters):
        key = lambda x: samples[x, i]
        index = list(range(samples.shape[0]))
        for i in range(self._N_parameters_cache - 1, -1, -1):
            if not parameters[i].system: continue
            index.sort(key=key)
        self._index = index
        
    def load_samples(self, samples):
        """
        Load samples for evaluation
        
        Parameters
        ----------
        samples : numpy.ndarray, dim=2
        
        """
        self._load_parameters()
        parameters = all_variables(self._system_parameters)
        N_parameters = self._N_parameters_cache
        if not isinstance(samples, np.ndarray):
            raise TypeError(f'samples must be an ndarray, not a {type(samples).__name__} object')
        if samples.ndim == 1:
            samples = samples[:, np.newaxis]
        elif samples.ndim != 2:
            raise ValueError('samples must be 2 dimensional')
        if samples.shape[1] != N_parameters:
            raise ValueError(f'number of parameters in samples ({samples.shape[1]}) must be equal to the number of parameters ({len(N_parameters)})')
        metrics = self.metrics
        self._load_sample_order(samples, parameters)
        empty_metric_data = np.zeros((len(samples), len(metrics)))
        self.table = pd.DataFrame(np.hstack((samples, empty_metric_data)),
                                  columns=var_columns(parameters + metrics))
        self._samples = samples
        
    def evaluate(self, jit=True, notify=False):
        """
        Evaluate metrics over the loaded samples and save values to `table`.
        
        Parameters
        ----------
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
            def evaluate(sample, count=[0]):
                count[0] += 1
                values = evaluate_sample(sample)
                print(f"{count} Elapsed time: {timer.elapsed_time:.0f} sec")
                return values
        else:
            evaluate = evaluate_sample
        index = self._index
        values = [None] * len(index)
        for i in index: values[i] = evaluate(samples[i])
        table[var_indices(self.metrics)] = values
        
    def _evaluate_sample(self, sample):
        scenarios = []
        all_metrics = [] 
        agile_scenario = self._agile_scenario
        system_parameters = self._system_parameters
        system_metrics = self._system_metrics
        general_parameters = system_parameters[None]
        general_metrics = system_metrics[None]
        Nf = len(general_parameters)
        for f, s in zip(general_parameters, sample[:Nf]): f.setter(s)
        for system in system_parameters:
            if not system: continue
            parameters = system_parameters[system]
            metrics = system_metrics[system]
            N = Nf
            Nf += len(parameters)
            for f, s in zip(parameters, sample[N:Nf]): f.setter(s)
            system.simulate()
            scenario = agile_scenario.create_scenario(system)
            scenarios.append(scenario)
            all_metrics.extend([i() for i in metrics])
        all_metrics.extend([i() for i in general_metrics])
        agile_scenario.compile_scenarios(scenarios)
    
    get_joint_distribution = Model.get_joint_distribution
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
    