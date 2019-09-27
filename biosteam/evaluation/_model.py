# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 23:17:38 2019

@author: yoelr
"""
import numpy as np
import pandas as pd
from ._state import State
from ._metric import Metric
from biosteam.utils import TicToc
from scipy.stats import spearmanr

__all__ = ('Model',)

# %% Functions

varindices = lambda vars: [var.index for var in vars]
varcolumns = lambda vars: pd.MultiIndex.from_tuples(
                                 varindices(vars),
                                 names=('Element', 'Variable'))

    
# %% Grid of simulation blocks

class Model(State):
    """Create a Model object that allows for evaluation over a sample space.
    
    Parameters
    ----------
    system : System
        Should reflect the model state.
    metrics : tuple[Metric]
        Metrics to be evaluated by model.
    skip=False : bool
        If True, skip simulation for repeated states.
    Examples
    --------
    :doc:`../tutorial/Advanced simulation`
    
    """
    __slots__ = ('_table',   # [DataFrame] All arguments and results
                 '_metrics', # tuple[Metric] Metrics to be evaluated by model
                 '_index',   # list[int] Order of sample evaluation for performance
                 '_samples', # [array] Argument sample space
                 '_setters') # list[function] Cached parameter setters
    
    def __init__(self, system, metrics, skip=False):
        super().__init__(system, skip)
        self.metrics = metrics
        self._setters = self._samples = self._table = None
        
    def copy(self):
        """Return copy."""
        copy = super().copy()
        copy._metrics = self._metrics
        if self._update:
            copy._setters = self._setters
            copy._table = self._table.copy()
        else:
            copy._setters = copy._samples = copy._table = None
        return copy
    
    def _erase(self):
        """Erase cached data."""
        self._setters = self._update = self._table = self._samples = None
    
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
        if not self._update: self._loadparams()
        params = self._params
        paramlen = len(params)
        if not isinstance(samples, np.ndarray):
            raise TypeError(f'samples must be an ndarray, not a {type(samples).__name__} object')
        elif samples.ndim != 2:
            raise ValueError('samples must be 2 dimensional')
        elif samples.shape[1] != paramlen:
            raise ValueError(f'number of parameters in samples ({samples.shape[1]}) must be equal to the number of parameters ({len(params)})')
        key = lambda x: samples[x][i]
        N_samples = len(samples)
        index = list(range(N_samples))
        for i in range(paramlen-1,  -1, -1):
            if not params[i].system: break
            index.sort(key=key)
        self._index = index
        empty_metric_data = np.zeros((N_samples, len(self.metrics)))
        self._table = pd.DataFrame(np.hstack((samples, empty_metric_data)),
                                   columns=varcolumns(tuple(params)+self.metrics))
        self._samples = samples
        
    def evaluate(self, thorough=True):
        """Evaluate metric over the argument sample space and save values to ``table``.
        
        Parameters
        ----------
        thorough : bool
            If True, simulate the whole system with each sample.
            If False, simulate only the affected parts of the system.
        
        """
        # Setup before simulation
        funcs = [i.getter for i in self._metrics]
        values = []
        add = values.append
        samples = self._samples
        if samples is None: raise ValueError('must load samples or distribution before evaluating')
        if thorough:
            if self._setters:
                setters = self._setters
            else:
                self._setters = setters = [p.setter for p in self._params]
            simulate = self._system.simulate
            zip_ = zip
            for i in self._index:
                for f, s in zip_(setters, samples[i]): f(s)
                simulate()
                add([i() for i in funcs])
        else:
            update = self._update
            for i in self._index: 
                update(samples[i])
                add([i() for i in funcs])    
        cols = varindices(self._metrics)
        for k, v in zip_(cols, zip_(*values)): self.table[k] = v
    
    def evaluate_across_coordinate(self, f_coordinate, name,
                                   lb, ub, xlfile,
                                   N_samples=1000, N_points=20,  rule='L',
                                   notify=True):
        """Evaluate across coordinate and save sample metrics.
        
        Parameters
        ----------
        f_coordinate : function
            Should change state of system given the coordinate.
        lb : float
            Lower bound of coordinate.
        ub : float
            Upper bound of coordinate.
        xlfile : str
            Name of file to save. File must end with ".xlsx"
        N_samples=1000 : int, optional
            Number of samples at each coordinate.
        N_points=20 : int, optional
            Number of lipid fractions evaluated.
        rule='L' : str
            Sampling rule.
        notify=True : bool, optional
            If True, notify elapsed time after each coordinate evaluation.
        
        """
        # Load samples
        samples = self.sample(N_samples, 'L')
        self.load_samples(samples)
        
        # Points to evaluate
        coordinate = np.linspace(lb, ub, N_points)
        
        # Initialize data containers
        metric_data = {}
        def new_data(key, dct=metric_data):
            data = np.zeros([N_samples, N_points])
            dct[key] = data
            return data
        for i in self.metrics: new_data(i.name)
        
        # Initialize timer
        if notify:
            timer = TicToc()
            timer.tic()
            def evaluate():
                self.evaluate()
                print(f"Elapsed time: {timer.elapsed_time:.0f} sec")
        else:
            evaluate = self.evaluate
        
        # Simulate lipidcane biorefinery    
        table = self.table
        for n, f in enumerate(coordinate):
            f_coordinate(f)
            evaluate()
            for metric in metric_data:
                metric_data[metric][:, n] = table[metric]
        
        # Save data to excel
        data = pd.DataFrame(data=np.zeros([N_samples, N_points]), columns=coordinate)
        data.columns.name = name
        
        with pd.ExcelWriter(xlfile) as writer:
            for metric in metric_data:
                data[:] = metric_data[metric]
                data.to_excel(writer, sheet_name=metric)
    
    def spearman(self, metrics=(), excluded_params=()):
        """Return DataFrame of Spearman's rank correlation for metrics vs parameters.
        
        Parameters
        ----------
        metrics=() : Iterable[Metric], defaults to all metrics
            Names of metrics to be be correlated with parameters.
        excluded_params=() : Iterable[str], optional
            Names of parameters to exclude.
        """
        data = self.table
        param_cols = list(data)
        all_metric_names = varindices(self.metrics)
        if not metrics: 
            metric_names = all_metric_names
        else:
            metric_names = varindices(metrics)
        params = list(self._params)
        for i in excluded_params:
            for j in params:
                if j.name == i:
                    params.remove(j)
                    break
            for j in param_cols: 
                if j[0] == i:
                    param_cols.remove(j)
                    break
        
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
    
    def __call__(self, sample):
        """Return list of metric values."""
        super().__call__(sample)
        return [i.getter() for i in self._metrics]
    
    def _repr(self):
        clsname = type(self).__name__
        newline = "\n" + " "*(len(clsname)+2)
        return f'{clsname}: {newline.join([i.describe() for i in self.metrics])}'
        
    def __repr__(self):
        return f'<{type(self).__name__}: {", ".join([i.name for i in self.metrics])}>'
    
        