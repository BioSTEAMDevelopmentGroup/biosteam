# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>,
#                          Yalin Li <zoe.yalin.li@gmail.com>
#
# This module implements a filtering feature from the stats module of the QSDsan library:
# QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems
# Copyright (C) 2020-2021, Yalin Li <zoe.yalin.li@gmail.com>
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
from ._parameter import Parameter
from ._utils import var_indices, var_columns, indices_to_multiindex
from biosteam.exceptions import FailedEvaluation
from warnings import warn
from collections.abc import Sized
import pickle

__all__ = ('Model',)


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
        Loads speficications once all parameters are set. Specification should 
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
        if self._N_parameters_cache:
            copy.table = self.table.copy()
        else:
            copy._samples = copy.table = None
        return copy
    
    def _erase(self):
        """Erase cached data."""
        super()._erase()
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
    
    def metric(self, getter=None, name=None, units=None, element='Biorefinery'):
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
        if not getter: return lambda getter: self.metric(getter, name, units, element)
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
    
    def _load_sample_order(self, samples, parameters):
        key = lambda x: samples[x, i]
        index = list(range(samples.shape[0]))
        for i in range(self._N_parameters_cache-1,  -1, -1):
            if not parameters[i].system: continue
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
        samples = self._sample_hook(samples, parameters)
        self._load_sample_order(samples, parameters)
        empty_metric_data = np.zeros((len(samples), len(metrics)))
        self.table = pd.DataFrame(np.hstack((samples, empty_metric_data)),
                                  columns=var_columns(parameters + metrics))
        self._samples = samples
        
    def single_point_sensitivity(self, etol=0.01):
        parameters = self.parameters
        bounds = [i.bounds for i in parameters]
        sample = self.get_baseline_sample()
        N_parameters = len(parameters)
        index = range(N_parameters)
        metrics = self.metrics
        N_metrics = len(metrics)
        values_lb = np.zeros([N_parameters, N_metrics])
        values_ub = values_lb.copy()
        evaluate = self._evaluate_sample
        baseline_1 = np.array(evaluate(sample, True))
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
            values_lb[i, :] = evaluate(sample_lb, True)
            values_ub[i, :] = evaluate(sample_ub, True)
        baseline_2 = np.array(evaluate(sample, True))
        error = np.abs(baseline_2 - baseline_1)
        index, = np.where(error > 1e-6)
        error = error[index]
        relative_error = error / np.maximum.reduce([np.abs(baseline_1[index]), np.abs(baseline_2[index])])
        for i, idx in enumerate(index):
            if relative_error[i] > etol:
                raise RuntimeError(
                    f"inconsistent model; {metrics[idx]} has a value of "
                    f"{baseline_1[idx]} before evaluating sensitivty and "
                    f"{baseline_2[idx]} after"
                )
        metric_index = var_columns(metrics)
        baseline = pd.Series(0.5 * (baseline_1 + baseline_2), index=metric_index)
        df_lb = pd.DataFrame(values_lb, 
                            index=var_columns(parameters),
                            columns=metric_index)
        df_ub = df_lb.copy()
        df_ub[:] = values_ub
        return baseline, df_lb, df_ub
    
    def evaluate(self, thorough=True, notify=0, 
                 file=None, autosave=0, autoload=False):
        """
        Evaluate metrics over the loaded samples and save values to `table`.
        
        Parameters
        ----------
        thorough : bool
            If True, simulate the whole system with each sample.
            If False, simulate only the affected parts of the system and
            skip simulation for repeated states.
        notify=0 : int, optional
            If 1 or greater, notify elapsed time after the given number of sample evaluations. 
        file : str, optional
            Name of file to save/load pickled evaluation results.
        autosave : int, optional
            If 1 or greater, save pickled evaluation results after the given number of sample evaluations.
        autoload : bool, optional
            Whether to load pickled evaluation results from file.
        
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
            from biosteam.utils import TicToc
            timer = TicToc()
            timer.tic()
            count = [0]
            def evaluate(sample, thorough, count=count):
                count[0] += 1
                values = evaluate_sample(sample, thorough)
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
        
        if autosave:
            layout = table.index, table.columns
            for number, i in enumerate(index, number): 
                values[i] = evaluate(samples[i], thorough)
                if not number % autosave: 
                    obj = (number, values, *layout)
                    with open(file, 'wb') as f: pickle.dump(obj, f)
        else:
            for i in index: values[i] = evaluate(samples[i], thorough)
        table[var_indices(self._metrics)] = values
    
    def _evaluate_sample(self, sample, thorough):
        try:
            self._update_state(sample, thorough)
            return [i() for i in self.metrics]
        except Exception as exception:
            self._reset_system()
            if self.retry_evaluation:
                try:
                    self._specification() if self._specification else self._system.simulate()
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
                                   re_evaluate=True, multi_coordinate=False):
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
        table = self.table
        N_samples, N_parameters = table.shape
        N_points = len(coordinate)
        
        # Initialize data containers
        metric_indices = var_indices(self.metrics)
        shape = (N_samples, N_points)
        metric_data = {i: np.zeros(shape) for i in metric_indices}
        f_evaluate = self.evaluate
        
        # Initialize timer
        if re_evaluate: 
            if notify_coordinate:
                from biosteam.utils import TicToc
                timer = TicToc()
                timer.tic()
                def evaluate():
                    f_evaluate(notify=notify)
                    print(f"[Coordinate {n}] Elapsed time: {timer.elapsed_time:.0f} sec")
            else:
                evaluate = f_evaluate
        else:
            evaluate = lambda: None
        
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
                for metric in self.metrics:
                    data[:] = metric_data[metric.index]
                    data.to_excel(writer, sheet_name=metric.short_description)
        return metric_data
    
    def spearman(self, parameters=None, metrics=None):
        warn(DeprecationWarning('this method will be depracated in biosteam 2.25; '
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
            data = np.array([[corr(p, m) for m in metric_data] for p in parameter_data])
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
    
    def _repr(self):
        clsname = type(self).__name__
        newline = "\n" + " "*(len(clsname)+2)
        return f'{clsname}: {newline.join([i.describe() for i in self.metrics])}'
        
    def __repr__(self):
        return f'<{type(self).__name__}: {", ".join([i.name for i in self.metrics])}>'
    
        