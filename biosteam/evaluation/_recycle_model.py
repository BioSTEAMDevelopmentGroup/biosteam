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
from typing import Optional, Callable
from scipy.spatial.distance import cdist
from ._parameter import Parameter
from .._system import System, JointRecycleData

@njit(cache=True)
def pearson_correlation_coefficient(X, y):
    n_predictors = X.shape[1]
    result = np.zeros(n_predictors)
    my = y.mean()
    ym = y - my
    ssym = np.dot(ym, ym)
    for i in range(n_predictors):
        x = X[:, i]
        mx = x.mean()
        xm = x - mx
        r_num = np.dot(xm, ym)
        r_den = np.sqrt(np.dot(xm, xm) * ssym)
        result[i] = r_num / r_den
    return result

@njit(cache=True)
def R2(y_actual, y_predicted):
    dy = y_actual - y_predicted
    dy_m = y_actual - np.mean(y_actual)
    SSR = np.dot(dy, dy)
    SST = np.dot(dy_m, dy_m)
    if SST == 0.:
        if SSR == 0: 
            return 1.
        else:
            return -np.inf
    else:
        return 1 - SSR / SST

@njit(cache=True)
def fit_linear_model(X, y):
    Xt = X.transpose()
    return np.linalg.inv(Xt @ X) @ Xt @ y

class LinearRegressor:
    __slots__ = (
        'coefficients',
    )
    def __init__(self):
        self.coefficients = None
        
    def fit(self, X, y):
        self.coefficients = fit_linear_model(X, y)
    
    def predict(self, x):
        return np.dot(self.coefficients, x[0])

    def __repr__(self):
        return f"{type(self).__name__}()"
    
    
def recycle_response(recycle, name, model):
    if name == 'T':
        response = RecycleTemperature(recycle, model)
    elif name == 'P':
        response = RecyclePressure(recycle, model)
    else:
        response = RecycleFlow(recycle, name, model)
    return response
    

class RecycleResponse:
    __slots__ = (
        'recycle', 'model', 'predictors',
    )
    def __init__(self, recycle, model):
        self.recycle = recycle
        self.model = model
        self.predictors = [] # Indices of predictors
    
    def fit(self, X, y):
        self.model.fit(X[:, self.predictors], y)

    def predict(self, x):
        return float(
            self.model.predict(x[None, self.predictors])
        )
    
    def predict_locally(self, X, y, x, x_min, x_range, distance, weight):
        index = self.predictors
        x = x[None, index]
        X = X[:, index]
        x_min = x_min[index]
        x_range = x_range[index]
        X_norm = (X - x_min) / x_range
        x_norm = (x - x_min) / x_range
        distances = cdist(x_norm, X_norm, metric=distance)
        exact_match = distances == 0
        if exact_match.any(): return y[exact_match[0]].mean()
        w = np.sqrt(weight(distances))
        X = X * w.transpose()
        y = y * w[0]
        self.model.fit(X, y)
        return float(
            self.model.predict(x)
        )
    
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
        return f"{type(self).__name__}({self.recycle}, {self.model})"


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
    
    def __init__(self, recycle, name, model):
        self.recycle = recycle
        self.name = name
        self.model = model
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
        'case_study',
        'model_type',
        'recess',
        'distance',
        'fitted',
        'predictors_lb',
        'predictors_range',
        'weight',
        'nfits',
        'local_weighted',
    )
    default_model_type = LinearRegressor
    def __init__(self, 
            system: System,
            predictors: tuple[Parameter],
            model_type: Optional[str|Callable]=None,
            recess: Optional[int]=None, 
            distance: Optional[str]=None,
            weight:  Optional[Callable]=None,
            nfits: Optional[int]=None,
            local_weighted: Optional[int]=None,
        ):
        self.system = system
        self.predictors = predictors
        self.responses: dict[tuple[Stream, str], RecycleResponse] = {}
        if model_type is None: model_type = self.default_model_type
        if isinstance(model_type, str):
            model_type = model_type.lower()
            if model_type == 'linear regressor':
                model_type = LinearRegressor
            elif model_type == 'linear svr': # linear support vector machine regression
                from sklearn.svm import LinearSVR
                from sklearn.pipeline import make_pipeline
                from sklearn.preprocessing import StandardScaler
                model_type = lambda: make_pipeline(
                    StandardScaler(),
                    LinearSVR()
                )
                if nfits is None: nfits = 4
            elif model_type == 'svr':
                from sklearn.svm import SVR
                from sklearn.pipeline import make_pipeline
                from sklearn.preprocessing import StandardScaler
                model_type = lambda: make_pipeline(
                    StandardScaler(),
                    SVR()
                )
                if nfits is None: nfits = 2
            else:
                raise ValueError('unknown model type {model_type!r}')
        if recess is None: 
            if model_type is LinearRegressor:
                recess = 0
            else:
                recess = 2 * sum([i.kind == 'coupled' for i in predictors]) ** 2
        if local_weighted is None:
            if model_type is LinearRegressor:
                local_weighted = True
            else:
                local_weighted = False
        if local_weighted:
            if distance is None: distance = 'cityblock'
            if recess: raise ValueError('local weighted recycle model cannot recess')
            if nfits: raise ValueError('local weighted recycle model must fit every time; cannot pass nfits argument')
            if weight is None: weight = lambda x: 1. / x
        self.model_type = model_type
        self.recess = recess
        self.distance = distance
        self.fitted = False
        self.load_responses(model_type)
        self.weight = weight
        self.nfits = nfits
        self.local_weighted = local_weighted
        
    def fitted_responses(self):
        data = self.data
        responses = self.responses
        fitted = {i: [] for i, j in responses.items()}
        samples = np.array(data['samples'])
        if not self.fitted:
            if self.local_weighted: self.fit()
            else: return fitted
        for i, sample in enumerate(samples):
            for key, response in responses.items():
                fitted[key].append(
                    response.predict(sample)
                )
        return fitted
    
    def R2(self, last=None):
        null, null_dct = self.R2_null(last)
        predicted, predicted_dct = self.R2_predicted(last)
        fitted, fitted_dct = self.R2_fitted()
        return (
            {'null': null,
             'predicted': predicted,
             'fitted': fitted},
            {'null': null_dct,
             'predicted': predicted_dct,
             'fitted': fitted_dct},   
        )
            
    def _R2(self, dataset, last=None):
        results = {}
        data = self.data
        actual = data['actual']
        fitted = dataset == 'fitted'
        if fitted and not self.fitted:
            if self.local_weighted: self.fit()
            else: return np.nan, {}
        case_studies = data['case studies']
        predicted = self.fitted_responses() if fitted else data[dataset]
        for key, response in self.responses.items():
            (recycle, name) = key
            name = f"{recycle}.{name}"
            y_actual = np.array(actual[key]) 
            if not fitted: y_actual = y_actual[case_studies]
            y_predicted = np.array(predicted[key])
            if last is not None:
                y_actual = y_actual[-last:]
                y_predicted = y_predicted[-last:]
            results[name] = R2(
                y_actual,
                y_predicted,
            )
        mean = sum(results.values()) / len(results)
        lb = min(results.values())
        ub = max(results.values())
        return (lb, mean, ub), results
    
    def R2_null(self, last=None):
        return self._R2('null', last)
        
    def R2_predicted(self, last=None):
        return self._R2('predicted', last)
    
    def R2_fitted(self):
        return self._R2('fitted')
    
    def practice(self, case_study):
        """
        Predict and set recycle responses given the sample, then append actual
        simulation result to data.
        
        Must be used in a with-statement as follows:
            
        ```python
        with recycle_model.practice(case_study): # the case study is an unsimulated sample
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
        self.case_study = case_study
        return self
        
    def __enter__(self):
        data = self.data
        null_responses = data['null']
        predicted = data['predicted']
        actual = data['actual']
        sample_list = data['samples']
        samples = np.array(sample_list)
        case_study = self.case_study
        n_samples = len(sample_list)
        if self.local_weighted:
            for key, response in self.responses.items():
                null_responses[key].append(response.get())
                prediction = response.predict_locally(
                    samples, np.array(actual[key]), case_study, 
                    self.predictors_lb, self.predictors_range, self.distance,
                    self.weight,
                )
                response.set(prediction)
                predicted[key].append(prediction)
        elif (not n_samples % (self.recess + 1)  # Recess is over
              and self.nfits is not None
              and self.fitted < self.nfits):
            self.fitted += 1
            for key, response in self.responses.items():
                null_responses[key].append(response.get())
                response.fit(samples, np.array(actual[key]))
                prediction = response.predict(case_study) 
                response.set(prediction)
                predicted[key].append(prediction)
        else:
            for key, response in self.responses.items():
                null = response.get()
                null_responses[key].append(null)
                if self.fitted:
                    prediction = response.predict(case_study) 
                    response.set(prediction)
                else:
                    prediction = null
                predicted[key].append(prediction)
        data['case studies'].append(n_samples)
        sample_list.append(case_study)
    
    def __exit__(self, type, exception, traceback, total=[]):
        data = self.data
        samples = data['samples']
        if exception: 
            del samples[-1]
            raise exception
        actual = data['actual']
        for key, response in self.responses.items():
            actual[key].append(response.get())
        systems = self.system.subsystems
        n = len(systems)
        total.append(sum([i._iter for i in systems]))
        print(
            sum(total) / len(total) / n
        )
        print(self.R2(10)[0])
        del self.case_study
        
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
        
    def load_responses(self, model_type): 
        """
        Select responses and their respective predictors through single point 
        sensitivity. Also store the simulation data for fitting later.
        """        
        predictors = self.predictors
        system = self.system
        bounds = [i.bounds for i in predictors]
        sample = [i.baseline for i in predictors]
        N_predictors = len(predictors)
        self.predictors_lb = predictors_lb = np.zeros(N_predictors)
        self.predictors_range = predictors_range = np.zeros(N_predictors)
        for i, (lb, ub) in enumerate(bounds):
            predictors_lb[i] = lb
            predictors_range[i] = ub - lb
        index = range(len(predictors))
        evaluate = self.evaluate_parameters
        baseline_1 = evaluate(sample)
        values = []
        values_at_bounds = []
        samples = [sample]
        for i, p in enumerate(predictors):
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
        arr1 = baseline_1.to_array()
        arr2 = baseline_2.to_array()
        error = np.abs(arr1 - arr2)
        index, = np.where(error > system.molar_tolerance)
        error = error[index]
        relative_error = error / np.maximum.reduce([np.abs(arr1[index]), np.abs(arr2[index])])
        tol = 0.01 + system.relative_molar_tolerance
        bad_index = [index[i] for i, bad in enumerate(relative_error > tol) if bad]
        if bad_index:
            keys = baseline_1.get_keys()
            names = baseline_1.get_names()
            bad_keys = set([keys[i] for i in bad_index])
            bad_names = [names[i] for i in bad_index]
            bad_names = ', '.join(bad_names)
            warn(
               f"inconsistent model; recycle loops on [{bad_names}] do not "
                "match at baseline before and after single point "
               f"sensitivity analysis ({100 * relative_error[bad_index]} % error)",
               RuntimeWarning
            )
        else:
            bad_keys = set()
        responses = self.add_sensitive_reponses(
            baseline_1, values_at_bounds, tol, model_type, bad_keys
        )
        self.data = data = {'samples': [], 'case studies': []}
        for name in ('actual', 'predicted', 'null'): 
            data[name] = {key: [] for key in responses}
        self.extend_data(samples, values)
        
    def add_sensitive_reponses(self, 
            baseline: JointRecycleData,
            bounds: tuple[JointRecycleData], 
            tol: float,
            model_type: object,
            exclude_responses: set[tuple[Stream, str]],
        ):
        responses = self.responses
        baseline_dct = baseline.to_dict()
        for p, (lb, ub) in enumerate(bounds):
            if lb is baseline is ub: continue
            recycle_data = (lb.to_dict(), baseline_dct, ub.to_dict())
            keys = set()
            for i in recycle_data: keys.update(i)
            for key in keys:
                if key in exclude_responses: continue
                values = [dct[key] for dct in recycle_data if key in dct]
                mean = np.mean(values)
                if any([(i - mean) / mean > tol for i in values]):
                    if key in responses:
                        response = responses[key]
                    else:
                        responses[key] = response = recycle_response(*key, model_type())
                    response.predictors.append(p)
        return responses
        
    def append_data(self, sample, recycle_data=None):
        data = self.data
        actual = data['actual']
        data['samples'].append(sample)
        if recycle_data is None: 
            for key, response in self.responses.items():
                actual[key].append(response.get())
        else:
            for key, value in recycle_data.to_tuples(): 
                if key in actual: actual[key].append(value)
            
    def extend_data(self, samples, recycle_data):
        for args in zip(samples, recycle_data): self.append_data(*args)
    
    def fit(self):
        data = self.data
        actual = data['actual']
        samples = np.array(data['samples'])
        for key, response in self.responses.items():
            response.fit(samples, np.array(actual[key]))    
    
    def predict(self, sample):
        for key, response in self.responses.items():
            prediction = response.predict(sample)
            response.set(prediction)