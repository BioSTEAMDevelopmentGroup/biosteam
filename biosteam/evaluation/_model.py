# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 23:17:38 2019

@author: yoelr
"""
import numpy as np
import pandas as pd
from ._state import State

__all__ = ('Model',)

# %% Functions

paramindex = lambda params: pd.MultiIndex.from_arrays(
                                ([f.element_name for f in params], 
                                 [f.name for f in params]),
                                names=('Element', 'Parameter'))

    
# %% Grid of simulation blocks

class Model(State):
    """Create a Model object that allows for evaluation over a sample space.
    
    **Parameters**
    
        **ID:** [str] ID of metric        
    
        **system:** [System] Should reflect the model state.
    
        **metric:** [function] Should return metric value.
    
    **Examples**

         :doc:`Advanced simulation`
    
    """
    __slots__ = ('_ID',      # [str] Should be the metric name.
                 '_table',   # [DataFrame] All arguments and results.
                 '_metric',  # [function] Should return metric being evaluated.
                 '_index',   # list[int] Order of sample evaluation for performance.
                 '_samples') # [array] Argument sample space.
    
    def __init__(self, ID, system, metric):
        super().__init__(system)
        self._ID = ID
        self._metric = metric
        self._samples = self._table = None
    
    def _erase(self):
        """Erase cached data."""
        self._model = self._table = self._samples = None
    
    @property
    def table(self):
        """[DataFrame] Table of the sample space with results in the final column."""
        return self._table
    
    def load_samples(self, samples):
        """Load samples for evaluation"""
        if not self._model: self._loadmodel()
        params = self._params
        paramlen = len(params)
        if not isinstance(samples, np.ndarray):
            raise TypeError(f'samples must be an ndarray, not a {type(samples).__name__} object')
        elif samples.ndim != 2:
            raise ValueError('samples must be 2 dimensional')
        elif samples.shape[1] != paramlen:
            raise ValueError(f'number of parameters in samples ({samples.shape[1]}) must be equal to the number of parameters ({len(params)})')
        key = lambda x: samples[x][i]
        index = list(range(len(samples)))
        for i in range(paramlen-1,  -1, -1):
            if not params[i].system: break
            index.sort(key=key)
        self._index = index
        self._table = pd.DataFrame(samples, columns=paramindex(params))
        self._table[self._ID] = None
        self._samples = samples
        
    def evaluate(self, default=None):
        """Evaluate metric over the argument sample space and save values to `table`."""
        if not self._model: self._loadmodel()
        # Setup before simulation
        metric = self._metric
        values = []
        add = values.append
        model = self._model
        samples = self._samples
        if samples is None: raise ValueError('must load samples or distribution before evaluating')
        for i in self._index: 
            try: 
                model(samples[i])
                add(metric())
            except: add(default)
        self.table[self._ID] = values
    
    def __call__(self, sample):
        super().__call__(sample)
        return self._metric()
    
    def _repr(self):
        return f'{type(self).__name__}: {self._ID}'
        
   
    
        