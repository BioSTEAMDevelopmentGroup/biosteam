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
from itertools import combinations, product
from ._parameter import Parameter
from .._system import System, JointRecycleData

__all__ = (
    'GenericResponse',
    'ConvergenceModel',
    'NullConvergenceModel',
    'Average',
    'LinearRegressor',
    'InterceptLinearRegressor',
)

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


class InterceptLinearRegressor:
    __slots__ = (
        'coefficients',
    )
    def __init__(self):
        self.coefficients = None
        
    def fit(self, X, y):
        m, n = X.shape
        Xi = np.ones([m, n + 1])
        Xi[:, 1:] = X 
        self.coefficients = fit_linear_model(Xi, y)
    
    def predict(self, x):
        xi = np.ones(x.shape[1] + 1)
        xi[1:] = x[0]
        return np.dot(self.coefficients, xi)

    def __repr__(self):
        return f"{type(self).__name__}()"


class Average:
    __slots__ = (
        'mean',
    )
    def __init__(self):
        self.mean = None
        
    def fit(self, X, y):
        self.mean = np.mean(y)
    
    def predict(self, x):
        return self.mean

    def __repr__(self):
        return f"{type(self).__name__}()"
    
    
fast_fit_model_types = set([
    LinearRegressor, 
    InterceptLinearRegressor,
    Average
])
linear_model_types = set([
    LinearRegressor, 
    InterceptLinearRegressor, 
])
    
def recycle_response(recycle, name, model):
    if name == 'T':
        response = RecycleTemperature(recycle, model)
    elif name == 'P':
        response = RecyclePressure(recycle, model)
    else:
        response = RecycleFlow(recycle, name, model)
    return response
    

class Response:
    __slots__ = (
        'element', 'model', 'predictors',
        'max', 'min',
    )
    def __init__(self, element, model):
        self.element = element
        self.model = model
        self.predictors = [] # Indices of predictors
        self.max = self.min = None
    
    def __hash__(self):
        return hash((self.element, self.name))
    
    def __eq__(self, other):
        return (self.element, self.name) == other
    
    def fit(self, X, y):
        self.model.fit(X[:, self.predictors], y)

    def predict(self, x):
        return self.filter_value(
            float(self.model.predict(x[None, self.predictors]))
        )
    
    def predict_locally(
            self, x, X, y, distance, weight
        ):
        index = self.predictors
        x = x[None, index]
        X = X[:, index]
        model = self.model
        distances = cdist(x, X, metric=distance)
        exact_match = distances == 0
        if exact_match.any(): return y[exact_match[0]].mean()
        w = weight(distances)
        X = X * w.transpose()
        y = y * w[0]
        self.model.fit(X, y)
        yi = model.predict(x)
        return self.filter_value(float(yi))
    
    def update_limits(self, value):
        if self.max is None:
            self.max = value
            self.min = value
        elif value > self.max:
            self.max = value
        elif value < self.min:
            self.min = value
    
    def filter_value(self, value):
        if value < self.min:
            return self.min
        elif value > self.max:
            return self.max
        else:
            return value
    
    def __str__(self):
        return f"{self.element}.{self.name}"
    
    def __repr__(self):
        return f"{type(self).__name__}({self.element!r})"


class GenericResponse(Response):
    __slots__ = (
        'name',
        'units',
    )
    
    def __init__(self, 
            element, name, model=None, units=None, max=None, min=None,
        ):
        self.name = name
        self.units = units
        self.max = max
        self.min = min
        super().__init__(element, model)
        
    def get(self):
        return getattr(self.element, self.name)
    
    def set(self, value):
        setattr(self.element, self.name, self.filter_value(value))
    
    def __repr__(self):
        return f"{type(self).__name__}({self.element!r}, {self.name!r})"
        

class RecycleTemperature(Response):
    __slots__ = (
        'element',
    )
    name = 'T'
    units = 'K'
    
    def get(self):
        return self.element.T
    
    def set(self, value):
        self.element.T = value
    

class RecyclePressure(Response):
    __slots__ = (
        'element',
    )
    name = 'P'
    units = 'Pa'
    
    def set(self, value):
        self.element.P = value
        
    def get(self):
        return self.element.P


class RecycleFlow(Response):
    __slots__ = (
        'element',
        'name',
    )
    units = 'kmol/hr'
    
    def __init__(self, element, name, model):
        self.name = name
        super().__init__(element, model)
        
    def set(self, value):
        self.element.imol[self.name] = self.filter_value(value)
    
    def get(self):
        return self.element.imol[self.name]

    def __repr__(self):
        return f"{type(self).__name__}({self.element!r}, {self.name!r})"
    

class ConvergenceModel:
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
        'sample_length',
        'sample_min',
        'sample_range',
        'weight',
        'nfits',
        'local_weighted',
        'normalization',
        'interaction_pairs',
        'predictor_index',
        'save_prediction',
    )
    absolute_response_tolerance = 0.001
    relative_response_tolerance = 0.01
    def __init__(self, 
            predictors: tuple[Parameter],
            model_type: Optional[Callable]=None,
            recess: Optional[int]=None, 
            distance: Optional[str]=None,
            weight:  Optional[Callable]=None,
            nfits: Optional[int]=None,
            local_weighted: Optional[int]=None,
            system: Optional[System] = None,
            responses: Optional[set[Response]]=None,
            interaction_pairs: Optional[bool] = None,
            normalization: Optional[bool] = None,
            load_responses: Optional[bool] = None,
            save_prediction: Optional[bool] = True,
        ):
        if system is None:
            systems = set([i.system for i in predictors])
            try:
                system, = systems
            except:
                if systems:
                    raise ValueError('predictors must share the same system')
                else:
                    raise ValueError('no system available')
        self.system = system
        self.responses = set() if responses is None else responses
        if model_type is None: 
            if len(predictors) > 10:
                model_type = LinearRegressor
            else:
                model_type = InterceptLinearRegressor
        if isinstance(model_type, str):
            model_type = model_type.lower()
            if model_type == 'linear regressor':
                model_type = LinearRegressor
            elif model_type == 'intercept linear regressor':
                model_type = InterceptLinearRegressor
            elif model_type == 'linear svr': # linear support vector machine regression
                from sklearn.svm import LinearSVR
                from sklearn.pipeline import make_pipeline
                from sklearn.preprocessing import StandardScaler
                model_type = lambda: make_pipeline(
                    StandardScaler(),
                    LinearSVR()
                )
                if nfits is None: nfits = 100
            elif model_type == 'svr':
                from sklearn.svm import SVR
                from sklearn.pipeline import make_pipeline
                from sklearn.preprocessing import StandardScaler
                model_type = lambda: make_pipeline(
                    StandardScaler(),
                    SVR()
                )
                if nfits is None: nfits = 10
            elif model_type == 'average':
                model_type = Average
            else:
                raise ValueError('unknown model type {model_type!r}')
        if recess is None: 
            if model_type in fast_fit_model_types:
                recess = 0
            else:
                recess = 5 * sum([i.coupled for i in predictors])
        if interaction_pairs is None: interaction_pairs = False
        if local_weighted is None:
            if model_type in fast_fit_model_types:
                local_weighted = True
            else:
                local_weighted = False
        if local_weighted:
            if normalization is None: normalization = model_type in linear_model_types
            if distance is None: distance = 'cityblock'
            if recess: raise ValueError('local weighted recycle model cannot recess')
            if nfits: raise ValueError('local weighted recycle model must fit every time; cannot pass nfits argument')
            if weight is None: weight = lambda x: np.exp(-x / 2)
        self.model_type = model_type
        self.recess = recess
        self.distance = distance
        self.fitted = 0
        self.weight = weight
        self.nfits = nfits
        self.local_weighted = local_weighted
        self.interaction_pairs = interaction_pairs
        self.normalization = normalization
        self.load_predictors(predictors)
        self.save_prediction = save_prediction
        if load_responses is None or load_responses: self.load_responses()
        
    def with_interactions(self, sample):
        n = len(sample)
        new_sample = np.zeros(self.sample_length)
        new_sample[:n] = sample
        for i, (a, b) in self.interaction_pairs:
            new_sample[i] = sample[a] * sample[b]
        return new_sample
    
    def normalize_sample(self, sample):
        return (sample - self.sample_min) / self.sample_range
    
    def reframe_sample(self, sample, predictors=None):
        if predictors is not None and predictors != self.predictors:
            index = {j: i for i, j in enumerate(self.predictors)}
            sample = self.data['samples'][-1]
            for p, s in zip(predictors, sample): sample[index[p]] = s
        if self.predictor_index is not None:
            sample = np.asarray(sample)[self.predictor_index]
        if self.normalization:
            sample = self.normalize_sample(sample)
        if self.interaction_pairs:
            sample = self.with_interactions(sample)
        return sample
        
    def fitted_responses(self):
        data = self.data
        responses = self.responses
        fitted = {i: [] for i in responses}
        samples = np.array(data['samples'])
        try: self.fit()
        except: return {}
        for i, sample in enumerate(samples):
            for response in responses:
                fitted[response].append(
                    response.predict(sample)
                )
        return fitted
    
    def R2(self, last=None):
        predicted, predicted_dct = self.R2_predicted(last)
        fitted, fitted_dct = self.R2_fitted()
        return (
            {'predicted': predicted,
             'fitted': fitted},
            {'predicted': predicted_dct,
             'fitted': fitted_dct},   
        )
            
    def _R2(self, dataset, last=None):
        results = {}
        data = self.data
        actual = data['actual']
        fitted = dataset == 'fitted'
        predicted = self.fitted_responses() if fitted else data[dataset]
        for response in self.responses:
            name = str(response)
            y_actual = np.array(actual[response])
            y_predicted = np.array(predicted[response])
            if last is None: 
                index = len(y_predicted)
            else:
                index = last
            y_actual = y_actual[-index:]
            y_predicted = y_predicted[-index:]
            results[name] = R2(
                y_actual,
                y_predicted,
            )
        mean = sum(results.values()) / len(results)
        lb = min(results.values())
        ub = max(results.values())
        return {'min': lb, 'mean': mean, 'max': ub}, results
        
    def R2_predicted(self, last=None):
        return self._R2('predicted', last)
    
    def R2_fitted(self):
        return self._R2('fitted')
    
    def practice(self, case_study, predictors=None):
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
        recycle_data = recycle_model.evaluate_system_convergence() # Or other simulation code.
        recycle_model.append_data(sample, recycle_data)
        ```
        
        Warning
        -------
        Does nothing if not used in with statement containing simulation.
        
        """
        self.case_study = self.reframe_sample(case_study, predictors)
        return self
        
    def __enter__(self):
        data = self.data
        if self.save_prediction: predicted = data['predicted']
        actual = data['actual']
        sample_list = data['samples']
        n_samples = len(sample_list)
        samples = np.array(sample_list)
        case_study = self.case_study
        if self.save_prediction: 
            if self.local_weighted:
                for response in self.responses:
                    prediction = response.predict_locally(
                        case_study, samples, np.array(actual[response]),
                        self.distance, self.weight,
                    )
                    response.set(prediction)
                    predicted[response].append(prediction)
            elif (not n_samples % (self.recess + 1)  # Recess is over
                  and (self.nfits is None or self.fitted < self.nfits)):
                self.fitted += 1
                for response in self.responses:
                    response.fit(samples, np.array(actual[response]))
                    prediction = response.predict(case_study) 
                    response.set(prediction)
                    predicted[response].append(prediction)
            elif self.fitted:
                for response in self.responses:
                    prediction = response.predict(case_study) 
                    response.set(prediction)
                    predicted[response].append(prediction)
        elif self.local_weighted:
            for response in self.responses:
                response.set(
                    response.predict_locally(
                        case_study, samples, np.array(actual[response]),
                        self.distance, self.weight,
                    )
                )
        elif (not n_samples % (self.recess + 1)  # Recess is over
              and (self.nfits is None or self.fitted < self.nfits)):
            self.fitted += 1
            for response in self.responses:
                response.fit(samples, np.array(actual[response]))
                response.set(response.predict(case_study))
        elif self.fitted:
            for response in self.responses:
                response.set(response.predict(case_study))
        sample_list.append(case_study)
    
    def __exit__(self, type, exception, traceback, total=[]):
        del self.case_study
        data = self.data
        if exception and self.fitted:
            del data['samples'][-1]
            raise exception
        actual = data['actual']
        for response in self.responses:
            value = response.get()
            response.update_limits(value)
            actual[response].append(value)
        
    def evaluate_system_convergence(self, sample, default=None, **kwargs):
        system = self.system
        for p, value in zip(self.predictors, sample): 
            p.setter(value)
            p.last_value = value
        try:
            system.simulate(design_and_cost=False, **kwargs)
        except Exception as error:
            if default is None: raise error
            # TODO: py3.11, add eror context
            #     raise RuntimeError(
            #         'baseline simulation failed; '
            #         'could not load responses for convergence prediction model'
            #     )
            system.empty_recycles()
            recycles_data = default
        else:
            recycles_data = system.get_recycle_data()
        return recycles_data
       
    def load_predictors(self, predictors):
        predictor_index = np.array([i.coupled for i in predictors])
        self.predictor_index = predictor_index = None if predictor_index.all() else np.where(predictor_index)[0]
        if predictor_index is not None: predictors = [predictors[i] for i in predictor_index]
        self.predictors = predictors
        bounds = [i.bounds for i in predictors]
        n = len(predictors)
        if self.interaction_pairs:
            self.interaction_pairs = [*enumerate(combinations(range(n), 2), n)]
            self.sample_length = int(n * (n - 1) / 2) + n
        else:
            self.interaction_pairs = ()
        if self.normalization:
            self.sample_min = sample_min = np.zeros(n)
            self.sample_range = sample_range = np.zeros(n)
            for i, (lb, ub) in enumerate(bounds):
                sample_min[i] = lb
                sample_range[i] = ub - lb
    
    def load_responses(self): 
        """
        Select material responses and their respective predictors through single point 
        sensitivity. Also store the simulation data for fitting later.
        """
        predictors = self.predictors
        responses = self.responses
        hooks = [i.hook for i in predictors]
        all_bounds = [i.bounds for i in predictors]
        sample = [i.baseline for i in predictors]
        evaluate = self.evaluate_system_convergence        
        baseline_1 = evaluate(sample)
        values = []
        values_at_bounds = []
        samples = [sample]
        for i, p in enumerate(predictors):
            hook = hooks[i]
            bounds = all_bounds[i]
            for limit in bounds:
                new_sample = sample.copy()
                if hook is not None: limit = hook(limit)
                new_sample[i] = limit
                samples.append(new_sample)
                values.append(
                    evaluate(new_sample, default=baseline_1, recycle_data=baseline_1)
                )
            values_at_bounds.append(values[-2:])
        if self.interaction_pairs:
            # [(lb, lb), (lb, ub), (ub, lb), (ub, ub)]
            pairwise_bounds_index = [*enumerate(product(range(2), range(2)))]
            for _, pair in self.interaction_pairs:
                new_samples = [sample.copy() for i in range(4)]
                for i, bounds_index in pairwise_bounds_index:
                    for b, p in zip(bounds_index, pair):
                        hook = hooks[p]
                        bounds = all_bounds[p]
                        limit = bounds[b]
                        new_sample = new_samples[i]
                        if hook is not None: limit = hook(limit)
                        new_sample[p] = limit
                for new_sample in new_samples:
                    values.append(
                        evaluate(new_sample, default=baseline_1, recycle_data=baseline_1)
                    ) 
                    samples.append(new_sample)
        baseline_2 = evaluate(sample, recycle_data=baseline_1)
        arr1 = baseline_1.to_array()
        arr2 = baseline_2.to_array()
        error = np.abs(arr1 - arr2)
        index, = np.where(error > self.absolute_response_tolerance)
        error = error[index]
        relative_error = error / np.maximum.reduce([np.abs(arr1[index]), np.abs(arr2[index])])
        rtol = self.relative_response_tolerance
        bad_index = [i for i, bad in enumerate(relative_error > rtol) if bad]
        if bad_index:
            keys = baseline_1.get_keys()
            names = baseline_1.get_names()
            relative_error = relative_error[bad_index]
            bad_index = [index[i] for i in bad_index]
            bad_keys = set([keys[i] for i in bad_index])
            bad_names = [names[i] for i in bad_index]
            bad_names = ', '.join(bad_names)
            warn(
                f"inconsistent model; recycle loops on [{bad_names}] do not "
                "match at baseline before and after single point "
                f"sensitivity analysis ({100 * relative_error} % error)",
                RuntimeWarning
            )
        else:
            bad_keys = set()
        self.add_sensitive_reponses(
            baseline_1, values_at_bounds, bad_keys
        )
        self.data = data = {'samples': []}
        for name in ('actual', 'predicted'): 
            data[name] = {key: [] for key in responses}
        predictor_index = self.predictor_index
        self.predictor_index = None
        self.extend_data(samples, values)
        self.predictor_index = predictor_index
        
    def add_sensitive_reponses(self, 
            baseline: JointRecycleData,
            bounds: tuple[JointRecycleData],
            exclude_responses: set[tuple[Stream, str]],
        ):
        model_type = self.model_type
        responses = self.responses
        baseline = baseline.to_dict()
        responses_dct = {}
        rtol = self.relative_response_tolerance
        for p, (lb, ub) in enumerate(bounds):
            if lb is baseline is ub: continue
            recycle_data = (lb.to_dict(), baseline, ub.to_dict())
            keys = set()
            for i in recycle_data: keys.update(i)
            for key in keys:
                if key in exclude_responses: continue
                values = [dct.get(key, 0.) for dct in recycle_data]
                mean = np.mean(values)
                if any([(i - mean) > mean * rtol for i in values]):
                    if key in responses_dct:
                        response = responses_dct[key]
                    else:
                        if isinstance(key, Response):
                            response = key
                        else:
                            response = recycle_response(*key, model_type())
                        responses_dct[key] = response
                        responses.add(response)
                    response.predictors.append(p)
        interaction_pairs = self.interaction_pairs
        for response in responses:
            if interaction_pairs:
                predictors_set = set(response.predictors)
                response.predictors.extend(
                    [i for i, (a, b) in interaction_pairs
                     if a in predictors_set and b in predictors_set]
                )
            if response.model is None: 
                response.model = model_type()
        
    def append_data(self, sample, recycle_data=None):
        data = self.data
        actual = data['actual']
        data['samples'].append(self.reframe_sample(sample))
        dct = recycle_data.to_dict()
        for response in self.responses:
            value = dct.get(response, 0.)
            response.update_limits(value)
            actual[response].append(value)
            
    def extend_data(self, samples, recycle_data):
        for args in zip(samples, recycle_data): self.append_data(*args)
    
    def fit(self):
        data = self.data
        actual = data['actual']
        samples = np.array(data['samples'])
        for response in self.responses:
            response.fit(samples, np.array(actual[response]))    
    
    def predict(self, sample):
        for response in self.responses:
            prediction = response.predict(sample)
            response.set(prediction)
            

class NullConvergenceModel:
    __slots__ = (
        'system',
        'predictors', 
        'predictor_index',
        'responses', 
        'data',
        'case_study',
        'interaction_pairs',
        'normalization',
    )
    absolute_response_tolerance = ConvergenceModel.absolute_response_tolerance
    relative_response_tolerance = ConvergenceModel.relative_response_tolerance
    practice = ConvergenceModel.practice
    evaluate_system_convergence = ConvergenceModel.evaluate_system_convergence
    load_responses = ConvergenceModel.load_responses
    load_predictors = ConvergenceModel.load_predictors
    add_sensitive_reponses = ConvergenceModel.add_sensitive_reponses
    append_data = ConvergenceModel.append_data
    extend_data = ConvergenceModel.extend_data
    reframe_sample = ConvergenceModel.reframe_sample
    
    def __init__(self, 
            predictors: tuple[Parameter],
            system: Optional[System] = None,
            responses: Optional[set[Response]]=None,
            load_responses: Optional[bool] = None,
        ):
        if system is None:
            systems = set([i.system for i in predictors])
            try:
                system, = systems
            except:
                if systems:
                    raise ValueError('predictors do not share the same system')
                else:
                    raise ValueError('no system available')
        self.system = system
        self.interaction_pairs = None
        self.normalization = None
        self.load_predictors(predictors)
        self.responses = set() if responses is None else responses
        if load_responses is None or load_responses: self.load_responses()
    
    def model_type(self): return None
    
    def R2(self, last=None):
        results = {}
        data = self.data
        actual = data['actual']
        null_responses = data['predicted']
        for response in self.responses:
            name = str(response)
            y_actual = np.array(actual[response]) 
            y_predicted = np.array(null_responses[response])
            if last is None: last = len(y_predicted)
            y_actual = y_actual[-last:]
            y_predicted = y_predicted[-last:]
            results[name] = R2(
                y_actual,
                y_predicted,
            )
        mean = sum(results.values()) / len(results)
        lb = min(results.values())
        ub = max(results.values())
        return {'min': lb, 'mean': mean, 'max': ub}, results
        
    def __enter__(self):
        data = self.data
        null_responses = data['predicted']
        sample_list = data['samples']
        sample_list.append(self.case_study)
        for response in self.responses:
            null_responses[response].append(response.get())
    
    def __exit__(self, type, exception, traceback, total=[]):
        del self.case_study
        data = self.data
        if exception: 
            data = self.data
            if exception:
                del data['samples'][-1]
                null_responses = data['predicted']
                for response in self.responses:
                    del null_responses[response][-1]
                raise exception
        actual = data['actual']
        for response in self.responses:
            actual[response].append(response.get())
    