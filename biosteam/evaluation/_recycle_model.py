# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from numba import njit
import numpy as np
from thermosteam import Stream
from warnings import warn
from ._parameter import Parameter
from .._system import System, JointRecycleData

@njit(cache=True)
def R2(y_true, y_pred):
    dy = y_true - y_pred
    dy_m = y_true - np.mean(y_true)
    return np.dot(dy, dy) / np.dot(dy_m, dy_m)

@njit(cache=True)
def fit_and_predict_linear_model(X, y, x):
    Xt = X.transpose()
    betas = np.linalg.inv(Xt @ X) @ Xt @ y
    return np.dot(betas, x)

def recycle_response(recycle, name):
    if name == 'T':
        response = RecycleTemperature(recycle)
    elif name == 'P':
        response = RecyclePressure(recycle)
    else:
        response = RecycleFlow(recycle, name)
    return response
    
class RecycleResponse:
    __slots__ = (
        'recycle', 'predictors',
    )
    def __init__(self, recycle):
        self.recycle = recycle
        self.predictors = [] # Indices of predictors
    
    def filter_value(self, value):
        old_value = self.get()
        min_value = self.min * old_value
        max_value = self.max * old_value
        if value < min_value:
            value = min_value
        elif value > max_value:
            value = max_value
        return value
    
    def __repr__(self):
        return f"{type(self).__name__}({self.recycle})"


class RecycleTemperature(RecycleResponse):
    __slots__ = ()
    units = 'K'
    max = 1.5
    min = 0.5
    
    def get(self):
        return self.recycle.T
    
    def set(self, value):
        self.recycle.T = self.filter_value(value)
    

class RecyclePressure(RecycleResponse):
    __slots__ = ()
    units = 'Pa'
    max = 1.5
    min = 0.5
    
    def set(self, value):
        self.recycle.P = self.filter_value(value)
        
    def get(self):
        return self.recycle.P


class RecycleFlow(RecycleResponse):
    __slots__ = (
        'name',
    )
    units = 'kmol/hr'
    max = 10.0
    min = 0.1
    
    def __init__(self, recycle, name):
        self.recycle = recycle
        self.name = name
        self.predictors = [] # Indices of predictors
        
    def set(self, value):
        self.recycle.imol[self.name] = self.filter_value(value)
    
    def get(self):
        return self.recycle.imol[self.name]

    def __repr__(self):
        return f"{type(self).__name__}({self.recycle}, {self.name!r})"


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
        self.responses: dict[tuple[Stream, str], RecycleResponse] = {}
        self.load_responses()
        
    def R2(self):
        return {
            'null': self.R2_null(),
            'predicted': self.R2_predicted(),
            'data': self.R2_data(),
        }
        
    def R2_null(self):
        pass
    
    def R2_predicted(self):
        pass
    
    def R2_data(self):
        pass
    
    def practice(self, sample):
        """
        Predict and set recycle responses given the sample, then append actual
        simulation result to data.
        
        Must be used in a with-statement as follows:
            
        ```python
        with recycle_model.practice(sample):
            recycle_model.system.simulate() # Or other simulation code.
            
        ```
        
        This method effectively does the same as running:
            
        ```python
        recycle_model.predict(sample)
        recycle_model.system.simulate() # Or other simulation code.
        recycle_model.append_data(sample)
        ```
        
        Warning
        -------
        Does nothing if not used in with statement containing simulation.
        
        """
        return self
        
    def __enter__(self, sample):
        data = self.data
        null_responses = data['null']
        predicted_responses = data['predicted']
        samples = np.array(data['samples'])
        for key, response in self.responses.items():
            null_responses[key].append(response.get())
            index = response.predictors
            prediction = fit_and_predict_linear_model(
                X=samples[:, index],
                y=data[key], 
                x=sample[index],
            )
            response.set(prediction)
            predicted_responses[key].append(prediction)
        data['samples'].append(sample)
    
    def __exit__(self, type, exception, traceback):
        data = self.data
        if exception: 
            samples = data['samples']
            del samples[-1]
            raise exception
        for key, response in self.responses.items():
            data[key].append(response.get())
        
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
            p = predictors[i]
            if p.kind != 'coupled': 
                values_at_bounds.append(
                    (baseline_1, baseline_1)
                )
                continue
            sample_lb = sample.copy()
            sample_ub = sample.copy()
            lb, ub = bounds[i]
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
        responses = self.add_sensitive_reponses(
            baseline_1, values_at_bounds, tol
        )
        self.data = data = {i: [] for i in ['samples', *responses]}
        data['null'] = {i: [] for i in responses}
        data['predicted'] = {i: [] for i in responses}
        self.extend_data(samples, values)
        
    def add_sensitive_reponses(self, 
            baseline: JointRecycleData,
            bounds: tuple[JointRecycleData], 
            tol: float
        ):
        responses = self.responses
        baseline_dct = baseline.to_dict()
        for p, (lb, ub) in enumerate(bounds):
            if lb is baseline is ub: continue
            recycle_data = (lb.to_dict(), baseline_dct, ub.to_dict())
            keys = set()
            for i in recycle_data: keys.update(i)
            for key in keys:
                values = [dct[key] for dct in recycle_data if key in dct]
                mean = np.mean(values)
                if any([(i - mean) / mean > tol for i in values]):
                    if key in responses:
                        response = responses[key]
                    else:
                        responses[key] = response = recycle_response(*key)
                    response.predictors.append(p)
        return responses
        
    def append_data(self, sample, recycle_data=None):
        data = self.data
        data['samples'].append(sample)
        if recycle_data is None: 
            for key, response in self.responses.items():
                data[key].append(response.get())
        else:
            for key, value in recycle_data.to_tuples(): 
                if key in data: data[key].append(value)
            
    def extend_data(self, samples, recycle_data):
        for args in zip(samples, recycle_data): self.append_data(*args)
        
    def predict(self, sample):
        data = self.data
        samples = np.array(data['samples'])
        for key, response in self.responses.items():
            index = response.predictors
            prediction = fit_and_predict_linear_model(
                X=samples[:, index],
                y=data[key], 
                x=sample[index],
            )
            response.set(prediction)