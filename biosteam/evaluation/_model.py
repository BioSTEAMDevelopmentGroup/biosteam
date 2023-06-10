# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>,
#                          Yalin Li <mailto.yalin.li@gmail.com>
#
# This module implements a filtering feature from the stats module of the QSDsan library:
# QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems
# Copyright (C) 2020-, Yalin Li <mailto.yalin.li@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
from scipy.spatial.distance import cdist
import numpy as np
import pandas as pd
from ._state import State
from ._metric import Metric
from ._feature import MockFeature
from ._parameter import Parameter
from ._utils import var_indices, var_columns, indices_to_multiindex
from biosteam.exceptions import FailedEvaluation
from warnings import warn
from collections.abc import Sized
from biosteam.utils import TicToc
import pickle

__all__ = ('Model',)

def replace_nones(values, replacement):
    for i, j in enumerate(values):
        if j is None: values[i] = replacement
    return values

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
        Loads specifications once all parameters are set. Specification should 
        simulate the system as well.
    params=None : Iterable[Parameter], optional
        Parameters to sample from.
    exception_hook : callable(exception, sample)
        Function called after a failed evaluation. The exception hook should 
        return either None or metric values given the exception and sample.

    """
    __slots__ = (
        'table',            # [DataFrame] All arguments and results.
        'retry_evaluation', # [bool] Whether to retry evaluation if it fails
        '_metrics',         # tuple[Metric] Metrics to be evaluated by model.
        '_index',           # list[int] Order of sample evaluation for performance.
        '_samples',         # [array] Argument sample space.
        '_exception_hook',  # [callable(exception, sample)] Should return either None or metric value given an exception and the sample.
    )
    def __init__(self, system, metrics=None, specification=None, 
                 parameters=None, retry_evaluation=True, exception_hook='warn'):
        super().__init__(system, specification, parameters)
        self.metrics = metrics or ()
        self.exception_hook = exception_hook
        self.retry_evaluation = retry_evaluation
        self.table = None
        self._erase()
        
    def copy(self):
        """Return copy."""
        copy = super().copy()
        copy._metrics = self._metrics
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
        evaluation. The exception hook should return either None or metric 
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
                def raise_exception(exception, sample): raise exception
                self._exception_hook = raise_exception
            else:
                raise ValueError(f"invalid exception hook name '{exception_hook}'; "
                                  "valid names include 'ignore', 'warn', and 'raise'")
        else:
            raise ValueError('exception hook must be either a callable, a string, or None')
    
    @property
    def metrics(self):
        """tuple[Metric] Metrics to be evaluated by model."""
        return tuple(self._metrics)
    @metrics.setter
    def metrics(self, metrics):
        metrics = list(metrics)
        isa = isinstance
        for i in metrics:
            if not isa(i, Metric):
                raise ValueError(f"metrics must be '{Metric.__name__}' "
                                 f"objects, not '{type(i).__name__}'")
        Metric.check_indices_unique(metrics)
        self._metrics = metrics
    
    def metric(self, getter=None, name=None, units=None, element=None):
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
        element : object, optional
                  Element being evaluated. Works mainly for bookkeeping. 
                  Defaults to 'Biorefinery'.
        
        Notes
        -----
        This method works as a decorator.
        
        """
        
        if isinstance(getter, Metric):
            if name is None: name = getter.name
            if units is None: units = getter.units
            if element is None: element = getter.element
            getter = getter.getter
        elif isinstance(getter, MockFeature):
            if element is None: element = getter.element
            if name is None: name = getter.name
            if units is None: units = getter.units
        elif not getter: 
            return lambda getter: self.metric(getter, name, units, element)
        metric = Metric(name, getter, units, element)
        Metric.check_index_unique(metric, self._metrics)
        self._metrics.append(metric)
        return metric 
    
    def _sample_hook(self, samples, parameters):
        if any([p.hook for p in parameters]):
            return np.array(
                [[(i if p.hook is None else p.hook(i))
                  for p, i in zip(parameters, row)]
                 for row in samples]
            )
        else:
            return samples
    
    def _load_sample_order(self, samples, parameters, distance):
        """
        Sort simulation order to optimize convergence speed
        by minimizing perturbations to the system between simulations.
        This function may take long depending on the number of parameters 
        because it uses single-point sensitivity to inform the sorting 
        algorithm.
        
        """
        N_samples = samples.shape[0]
        if N_samples < 2: 
            self._index = list(range(N_samples))
            return
        if distance is None: distance = 'cityblock'
        length = samples.shape[0]
        columns = [i for i, parameter in enumerate(self._parameters) if parameter.kind == 'coupled']
        parameters = [parameters[i] for i in columns]
        samples = samples.copy()
        samples = samples[:, columns]
        samples_min = samples.min(axis=0)
        samples_max = samples.max(axis=0)
        samples_diff = samples_max - samples_min
        normalized_samples = (samples - samples_min) / samples_diff
        # Note: Not sure if to deprecate or fix.
        # original_parameters = self._parameters
        # if ss: 
        #     def evaluate(sample, material_data, **kwargs):
        #         try:
        #             self._parameters = parameters
        #             for f, s in zip(self._parameters, sample): 
        #                 f.setter(s if f.scale is None else f.scale * s)
        #             diff = self._system.converge(material_data=material_data, **kwargs)
        #         finally:
        #             self._parameters = original_parameters
        #         return diff
        #     sample = [i.baseline for i in parameters]
        #     N_parameters = len(parameters)
        #     index = range(N_parameters)
        #     material_data = self._system.get_material_data()
        #     evaluate(sample, material_data, update_material_data=True)
        #     N_elements = material_data.material_flows.size
        #     diffs = np.zeros([N_parameters, N_elements])
        #     for i in index:
        #         sample_lb = sample.copy()
        #         sample_ub = sample.copy()
        #         lb = samples_min[i]
        #         ub = samples_max[i]
        #         hook = parameters[i].hook
        #         if hook:
        #             sample_lb[i] = hook(lb)
        #             sample_ub[i] = hook(ub)
        #         else:
        #             sample_lb[i] = lb
        #             sample_ub[i] = ub
        #         lb = evaluate(sample_lb, material_data).flatten()
        #         ub = evaluate(sample_ub, material_data).flatten()
        #         diffs[i] = ub - lb
        #     div = np.abs(diffs).max(axis=0)
        #     div[div == 0.] = 1
        #     diffs /= div
        #     normalized_samples = normalized_samples @ diffs
        nearest_arr = cdist(normalized_samples, normalized_samples, metric=distance)
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
        
    def load_samples(self, samples=None, optimize=None, file=None, 
                     autoload=None, autosave=None, distance=None):
        """
        Load samples for evaluation.
        
        Parameters
        ----------
        samples : numpy.ndarray, dim=2, optional
            All parameter samples to evaluate.
        optimize : bool, optional
            Whether to internally sort the samples to optimize convergence speed
            by minimizing perturbations to the system between simulations.
            Defaults to False.
        file : str, optional
            File to load/save samples and simulation order to/from.
        autosave : bool, optional
            Whether to save samples and simulation order to file (when not loaded from file).
        autoload : bool, optional
            Whether to load samples and simulation order from file (if possible).
        distance : str, optional
            Distance metric used for sorting. Defaults to 'cityblock'.
            See scipy.spatial.distance.cdist for options.
        
        Warning
        -------
        Depending on the number of parameters, optimizing sample simulation order 
        using single-point sensitivity may take long.
        
        """
        parameters = self._parameters
        Parameter.check_indices_unique(parameters)
        if autoload:
            try:
                with open(file, "rb") as f: (self._samples, self._index) = pickle.load(f)
            except FileNotFoundError: pass
            else:
                if (samples is None or (samples.shape == self._samples.shape and (samples == self._samples).all())):
                    metrics = self._metrics
                    empty_metric_data = np.zeros((len(samples), len(metrics)))
                    self.table = pd.DataFrame(np.hstack((samples, empty_metric_data)),
                                              columns=var_columns(parameters + metrics),
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
        metrics = self._metrics
        samples = self._sample_hook(samples, parameters)
        if optimize: 
            self._load_sample_order(samples, parameters, distance)
        else:
            self._index = list(range(samples.shape[0]))
        empty_metric_data = np.zeros((len(samples), len(metrics)))
        self.table = pd.DataFrame(np.hstack((samples, empty_metric_data)),
                                  columns=var_columns(parameters + metrics),
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
            etol=0.01, array=False, parameters=None, metrics=None, evaluate=None, 
            **kwargs
        ):
        if parameters is None: parameters = self.parameters
        bounds = [i.bounds for i in parameters]
        sample = [i.baseline for i in parameters]
        N_parameters = len(parameters)
        index = range(N_parameters)
        if metrics is None: metrics = self.metrics
        N_metrics = len(metrics)
        values_lb = np.zeros([N_parameters, N_metrics])
        values_ub = values_lb.copy()
        if evaluate is None: evaluate = self._evaluate_sample
        baseline_1 = np.array(evaluate(sample, **kwargs))
        sys = self.system
        if not sys.isdynamic: kwargs['material_data'] = sys.get_material_data()
        for i in index:
            sample_lb = sample.copy()
            sample_ub = sample.copy()
            lb, ub = bounds[i]
            hook = parameters[i].hook
            if hook:
                sample_lb[i] = hook(lb)
                sample_ub[i] = hook(ub)
            else:
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
                raise RuntimeError(
                    f"inconsistent model; {metrics[idx]} has a value of "
                    f"{baseline_1[idx]} before evaluating sensitivity and "
                    f"{baseline_2[idx]} after"
                )
        baseline = 0.5 * (baseline_1 + baseline_2)
        if array:
            return baseline, values_lb, values_ub
        else:
            metric_index = var_columns(metrics)
            baseline = pd.Series(baseline, index=metric_index)
            df_lb = pd.DataFrame(values_lb, 
                                index=var_columns(parameters),
                                columns=metric_index)
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
        metrics = self._metrics
        table[var_indices(metrics)] = replace_nones(values, [np.nan] * len(metrics))
    
    def evaluate(self, notify=0, file=None, autosave=0, autoload=False,
                 **kwargs):
        """
        Evaluate metrics over the loaded samples and save values to `table`.
        
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
        if notify:
            timer = TicToc()
            timer.tic()
            count = [0]
            def evaluate(sample, count=count, **kwargs):
                count[0] += 1
                values = evaluate_sample(sample, **kwargs)
                if not count[0] % notify:
                    print(f"{count} Elapsed time: {timer.elapsed_time:.0f} sec")
                return values
        else:
            evaluate = evaluate_sample
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
                values = [None] * len(index)   
            else:
                if notify: count[0] = number
        else:
            number = 0
            index = self._index
            values = [None] * len(index)
        
        export = 'export_state_to' in kwargs
        layout = table.index, table.columns
        try:
            for number, i in enumerate(index, number + 1): 
                if export: kwargs['sample_id'] = i
                values[i] = evaluate(samples[i], **kwargs)
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
            table[var_indices(self._metrics)] = replace_nones(values, [np.nan] * len(self.metrics))
    
    def _evaluate_sample(self, sample, **kwargs):
        state_updated = False
        try:
            self._update_state(sample, **kwargs)
            state_updated = True
            return [i() for i in self.metrics]
        except Exception as exception:
            if self.retry_evaluation and state_updated:
                self._reset_system()
                try:
                    self._update_state(sample, **kwargs)
                    self._specification() if self._specification else self._system.simulate(**kwargs)
                    return [i() for i in self.metrics]
                except Exception as new_exception: 
                    if self._exception_hook: 
                        values = self._exception_hook(new_exception, sample)
                        self._reset_system()
                        if isinstance(values, Sized) and len(values) == len(self.metrics):
                            return values
                        elif values is not None:
                            raise RuntimeError('exception hook must return either None or '
                                               'an array of metric values for the given sample')
            if self._exception_hook: 
                values = self._exception_hook(exception, sample)
                self._reset_system()
                if isinstance(values, Sized) and len(values) == len(self.metrics):
                    return values
                elif values is not None:
                    raise RuntimeError('exception hook must return either None or '
                                       'an array of metric values for the given sample')
            return self._failed_evaluation()
    
    def _reset_system(self):
        self._system.empty_recycles()
        self._system.reset_cache()
    
    def _failed_evaluation(self):
        self._reset_system()
        return [np.nan] * len(self.metrics)
    
    def metrics_at_baseline(self):
        """Return metric values at baseline sample."""
        baseline = self.get_baseline_sample()
        return self(baseline)
    
    def evaluate_across_coordinate(self, name, f_coordinate, coordinate,
                                   *, xlfile=None, notify=0, notify_coordinate=True,
                                   multi_coordinate=False, 
                                   f_evaluate=None):
        """
        Evaluate across coordinate and save sample metrics.
        
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
        notify : bool, optional
            If True, notify elapsed time after each coordinate evaluation. Defaults to True.
        f_evaluate : callable, optional
            Function to evaluate model. Defaults to evaluate method.
        
        """
        N_points = len(coordinate)
        if f_evaluate is None: f_evaluate = self.evaluate
        
        # Initialize timer
        if notify_coordinate:
            from biosteam.utils import TicToc
            timer = TicToc()
            timer.tic()
            def evaluate():
                f_evaluate(notify=notify)
                print(f"[Coordinate {n}] Elapsed time: {timer.elapsed_time:.0f} sec")
        else:
            evaluate = f_evaluate
        
        metric_data = None
        for n, x in enumerate(coordinate):
            f_coordinate(*x) if multi_coordinate else f_coordinate(x)
            evaluate()
            # Initialize data containers dynamically in case samples are loaded during evaluation
            if metric_data is None:
                N_samples, _ = self.table.shape
                metric_indices = var_indices(self.metrics)
                shape = (N_samples, N_points)
                metric_data = {i: np.zeros(shape) for i in metric_indices}
            for metric in metric_data:
                metric_data[metric][:, n] = self.table[metric]
        
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
                for metric in self.metrics:
                    data[:] = metric_data[metric.index]
                    data.to_excel(writer, sheet_name=metric.short_description)
        return metric_data
    
    def spearman(self, parameters=None, metrics=None):
        warn(DeprecationWarning('this method will be deprecated in biosteam 2.25; '
                                'use spearman_r instead'), stacklevel=2)
        return self.spearman_r(parameters, metrics)[0]
    
    def spearman_r(self, parameters=None, metrics=None, filter=None, **kwargs): # pragma: no cover
        """
        Return two DataFrame objects of Spearman's rho and p-values between metrics 
        and parameters.
        
        Parameters
        ----------
        parameters : Iterable[Parameter], optional
            Parameters to be correlated with metrics. Defaults to all parameters.
        metrics : Iterable[Metric], optional 
            Metrics to be correlated with parameters. Defaults to all metrics.
        
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
        return self._correlation(spearmanr, parameters, metrics, filter, kwargs)
    
    def pearson_r(self, parameters=None, metrics=None, filter=None, **kwargs):
        """
        Return two DataFrame objects of Pearson's rho and p-values between metrics 
        and parameters.
        
        Parameters
        ----------
        parameters : Iterable[Parameter], optional
            Parameters to be correlated with metrics. Defaults to all parameters.
        metrics : Iterable[Metric], optional 
            Metrics to be correlated with parameters. Defaults to all metrics.
            
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
        return self._correlation(pearsonr, parameters, metrics, filter, kwargs)
    
    def kendall_tau(self, parameters=None, metrics=None, filter=None, **kwargs):
        """
        Return two DataFrame objects of Kendall's tau and p-values between metrics 
        and parameters.
        
        Parameters
        ----------
        parameters : Iterable[Parameter], optional
            Parameters to be correlated with metrics. Defaults to all parameters.
        metrics : Iterable[Metric], optional 
            Metrics to be correlated with parameters. Defaults to all metrics.
            
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
        return self._correlation(kendalltau, parameters, metrics, filter, kwargs)
    
    def kolmogorov_smirnov_d(self, parameters=None, metrics=None, thresholds=[],
                             filter=None, **kwargs):
        """
        Return two DataFrame objects of Kolmogorov–Smirnov's D and p-values
        with the given thresholds for the metrics.
        
        For one particular parameter, all of the sampled values will be divided into two sets,
        one where the resulting metric value is smaller than or equal to the threshold,
        and the other where the resulting value is larger than the threshold.
        
        Kolmogorov–Smirnov test will then be performed for these two sets of values
        for the particular parameter.
        
        Parameters
        ----------
        parameters : Iterable[Parameter], optional
            Parameters to be correlated with metrics. Defaults to all parameters.
        metrics : Iterable[Metric], optional 
            Metrics to be correlated with parameters. Defaults to all metrics.
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
        metrics = metrics or self.metrics
        if len(thresholds) != len(metrics):
            raise ValueError(f'The number of metrics {len(metrics)} must match '
                             f'the number of thresholds ({len(thresholds)}).')
        kwargs['thresholds'] = thresholds
        return self._correlation(kstest, parameters, metrics, filter, kwargs)
    
    def _correlation(self, f, parameters, metrics, filter, kwargs):
        """
        Return two DataFrame objects of statistics and p-values between metrics 
        and parameters.
        
        Parameters
        ----------
        f : Callable
            Function with signature f(x, y) -> stat, p
        parameters : Iterable[Parameter], optional
            Parameters to be correlated with metrics. Defaults to all parameters.
        metrics : Iterable[Metric], optional 
            Metrics to be correlated with parameters. Defaults to all metrics.
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
        metric_indices = var_indices(metrics or self.metrics)
        metric_data = [values[index(i)] for i in metric_indices]
                
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
            data = np.array([[corr(p.astype(float), m.astype(float)) for m in metric_data] for p in parameter_data])
        else: # KS test
            thresholds = kwargs.pop('thresholds')
            data = np.array(
                [
                    [corr(parameter_data[n_p][metric_data[n_m]<=thresholds[n_m]],
                          parameter_data[n_p][metric_data[n_m]>thresholds[n_m]]) \
                     for n_m in range(len(metric_data))] \
                for n_p in range(len(parameter_data))])
        
        index = indices_to_multiindex(parameter_indices, ('Element', 'Parameter'))
        columns = indices_to_multiindex(metric_indices, ('Element', 'Metric'))        
        return [pd.DataFrame(i, index=index, columns=columns) for i in (data[..., 0], data[..., 1])]
        
    def create_fitted_model(self, parameters, metrics): # pragma: no cover
        from pipeml import FittedModel
        Xdf = self.table[[i.index for i in parameters]]
        ydf_index = metrics.index if isinstance(metrics, Metric) else [i.index for i in metrics]
        ydf = self.table[ydf_index]
        return FittedModel.from_dfs(Xdf, ydf)
    
    def __call__(self, sample):
        """Return pandas Series of metric values at given sample."""
        super().__call__(sample)
        return pd.Series({i.index: i() for i in self._metrics})
    
    def _repr(self, m):
        clsname = type(self).__name__
        newline = "\n" + " "*(len(clsname)+2)
        return f'{clsname}: {newline.join([i.describe() for i in self.metrics])}'
    
    def __repr__(self):
        return f'<{type(self).__name__}: {len(self.parameters)}-parameters, {len(self.metrics)}-metrics>'
    
    def _info(self, p, m):
        info = f'{type(self).__name__}:'
        parameters = self._parameters
        if parameters: 
            if p is None: p = len(parameters)
            ptitle = 'parameters: '
            info += '\n' + ptitle
            newline = "\n" + " "*(len(ptitle))
            info += newline.join([
                parameters[i].describe(
                    distribution=False, bounds=False
                ) for i in range(p)
            ])
        else:
            info += '\n(No parameters)'
        metrics = self._metrics
        if metrics: 
            if m is None: m = len(metrics)
            mtitle = 'metrics: '
            info += '\n' + mtitle
            newline = "\n" + " "*(len(mtitle))
            info += newline.join([
                metrics[i].describe() for i in range(m)
            ])
        else:
            info += '\n(No metrics)'
        return info
    
    def show(self, p=None, m=None):
        """Return information on p-parameters and m-metrics."""
        print(self._info(p, m))
    _ipython_display_ = show