# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import numpy as np
from thermosteam import Stream
from warnings import warn
from ._feature import Feature
from ._parameter import Parameter
from .._system import System, JointRecycleData

class RecycleResponse(Feature):
    __slots__ = (
        'predictors',
    )
    T_max = 1.0
    T_min = 0.5
    P_max = 1.0
    P_min = 0.5
    flow_max = 10.0
    flow_min = 0.1
    def __init__(self, element, name):
        self.element = element
        self.name = name
        match name:
            case 'T':
                units = 'K'
            case 'P':
                units = 'Pa'
            case _:
                units = 'kmol/hr'
        self.units = units
        self.predictors = [] # Indices of predictors
    
    def get_predictors(self, samples):
        return samples[..., self.predictors]
    
    def set(self, value):
        old_value = self.get()
        name = self.name
        recycle = self.element
        match name:
            case 'T':
                min_value = self.T_min * old_value
                max_value = self.T_max * old_value
                if value < min_value:
                    value = min_value
                elif value > max_value:
                    value = max_value
                recycle.T = value
            case 'P':
                min_value = self.P_min * old_value
                max_value = self.P_max * old_value
                if value < min_value:
                    value = min_value
                elif value > max_value:
                    value = max_value
                recycle.P = value
            case _:
                min_value = self.flow_min * old_value
                max_value = self.flow_max * old_value
                if value < min_value:
                    value = min_value
                elif value > max_value:
                    value = max_value
                recycle.imol[name] = value
    
    def get(self):
        name = self.name
        recycle = self.element
        match name:
            case 'T':
                return recycle.T
            case 'P':
                return recycle.P
            case _:
                return recycle.imol[name]

class RecycleModel:
    __slots__ = (
        'system', 
        'predictors', 
        'responses', 
        'data',
    )
    def __init__(self, system: System, predictors: tuple[Parameter]):
        self.system = system
        self.predictors = predictors
        self.load_responses()
        
    def evaluate_parameters(self, sample, default=None, **kwargs):
        for p, value in zip(self.predictors, sample):
            if p.scale is not None: value *= p.scale
            p.setter(value)
        try:
            self.system.simulate(design_and_cost=False, **kwargs)
            return self.system.get_recycle_data()
        except:
            self.system.empty_recycles()
            return default
        
    def load_responses(self): 
        """
        Select responses and their respective predictors through single point 
        sensitivity. Also store the simulation data for fitting later.
        """        
        predictors = self.predictors
        system = self.system
        bounds = [i.bounds for i in predictors]
        sample = [i.baseline for i in predictors]
        N_predictors = len(predictors)
        index = range(N_predictors)
        evaluate = self.evaluate_parameters
        baseline_1 = evaluate(sample)
        values = []
        values_at_bounds = []
        samples = [sample]
        for i in index:
            sample_lb = sample.copy()
            sample_ub = sample.copy()
            lb, ub = bounds[i]
            p = predictors[i]
            hook = p.hook
            if hook is not None:
                lb = hook(lb)
                ub = hook(ub)
            sample_lb[i] = lb
            sample_ub[i] = ub
            samples.append(sample_lb)
            samples.append(sample_ub)
            values_lb = evaluate(sample_lb, default=baseline_1, recycle_data=baseline_1)
            values_ub = evaluate(sample_ub, default=baseline_1, recycle_data=baseline_1)
            values.append(values_lb)
            values.append(values_ub)
            values_at_bounds.append(
                (values_lb, values_ub)
            )
        baseline_2 = evaluate(sample, recycle_data=baseline_1)
        names = baseline_1.get_names()
        assert names == baseline_2.get_names()
        arr1 = baseline_1.to_array()
        arr2 = baseline_2.to_array()
        error = np.abs(arr1 - arr2)
        index, = np.where(error > system.molar_tolerance)
        error = error[index]
        names = [names[i] for i in index]
        relative_error = error / np.maximum.reduce([np.abs(arr1[index]), np.abs(arr2[index])])
        tol = 0.01 + system.relative_molar_tolerance
        bad_index = [i for i, bad in enumerate(relative_error > tol) if bad]
        if bad_index:
            bad_names = [names[i] for i in bad_index]
            bad_names = ', '.join(bad_names)
            warn(
               f"inconsistent model; recycle loops on [{bad_names}] do not "
                "match at baseline before and after single point "
               f"sensitivity analysis ({100 * relative_error[bad_index]} % error)",
               RuntimeWarning
            )
        self.responses: dict[tuple[Stream, str], RecycleResponse] = self.create_sensitive_reponses(
            baseline_1, values_at_bounds, tol
        )
        self.data = {i: [] for i in ['samples', *self.responses]}
        self.extend_data(samples, values)
        
    def create_sensitive_reponses(self, 
            baseline: JointRecycleData,
            bounds: tuple[JointRecycleData], 
            tol: float
        ):
        responses = {}
        baseline = baseline.to_dict()
        for p, (lb, ub) in enumerate(bounds):
            recycle_data = (lb.to_dict(), baseline, ub.to_dict())
            keys = set()
            for i in recycle_data: keys.update(i)
            for key in keys:
                values = [dct[key] for dct in recycle_data if key in dct]
                mean = np.mean(values)
                if any([(i - mean) / mean > tol for i in values]):
                    if key in responses:
                        response = responses[key]
                    else:
                        responses[key] = response = RecycleResponse(*key)
                    response.predictors.append(p)
        return responses
        
    def append_data(self, sample, recycle_data=None):
        data = self.data
        data['samples'].append(sample)
        if recycle_data is None: recycle_data = self.system.get_recycle_data()
        for key, value in recycle_data.to_tuples(): 
            if key in data: data[key].append(value)
            
    def extend_data(self, samples, recycle_data):
        for args in zip(samples, recycle_data): self.append_data(*args)
        
    def predict(self, sample):
        data = self.data
        samples = np.array(data['samples'])
        inverse = np.linalg.inv
        for key, response in self.responses.items():
            y = data[key]
            X = response.get_predictors(samples)
            Xt = X.transpose()
            betas = inverse(Xt @ X) @ Xt @ y
            x = response.get_predictors(sample)
            response.set(np.dot(betas, x))