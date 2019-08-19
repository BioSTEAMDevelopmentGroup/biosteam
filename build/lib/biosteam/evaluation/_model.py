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
    
        **system:** [System] Should reflect the model state.
    
        **skip:** [bool] If True, skip simulation for repeated states
    
        ****metrics:** [str: function] ID-metric pairs.
    
    **Examples**

         :doc:`Advanced simulation`
    
    """
    __slots__ = ('_table',   # [DataFrame] All arguments and results
                 '_metrics', # dict[ID: function] Functions should return evaluated metric
                 '_index',   # list[int] Order of sample evaluation for performance
                 '_samples', # [array] Argument sample space
                 '_setters', # list[function] Cached parameter setters
                 '_repeats') # [bool] True if parameter states repeat
    
    def __init__(self, system, skip=False, **metrics):
        super().__init__(system, skip)
        self._metrics = metrics
        self._repeats = self._setters = self._samples = self._table = None
    
    def _erase(self):
        """Erase cached data."""
        self._setters = self._model = self._table = self._samples = None
    
    @property
    def metrics(self):
        """dict[str: function] ID-metric pairs."""
        return self._metrics
    @metrics.setter
    def metrics(self, metrics):
        table = self._table
        if table is not None:
            for i in self._metrics:
                if i in table: del table[i]
        self._metrics = metrics
    
    @property
    def table(self):
        """[DataFrame] Table of the sample space with results in the final column."""
        return self._table
    
    def load_samples(self, samples, has_repeats=None):
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
        self._samples = samples
        if has_repeats is None:
            try:
                j = samples[0]
                for i in samples[1:]:
                    if any(i==j):
                        self._repeats = True
                        break
            except: pass
        else:
            self._repeats = has_repeats
        
    def evaluate(self, default=None):
        """Evaluate metric over the argument sample space and save values to `table`."""
        # Setup before simulation
        funcs = tuple(self._metrics.values())
        values = []
        add = values.append
        samples = self._samples
        if samples is None: raise ValueError('must load samples or distribution before evaluating')
        if self._repeats:
            model = self._model
            for i in self._index: 
                try: 
                    model(samples[i])
                    add([i() for i in funcs])
                except Exception as err:
                    if default: add(default)
                    else: raise err            
        else:
            if self._setters:
                setters = self._setters
            else:
                self._setters = setters = [p.setter for p in self._params]
            add = values.append
            sim = self._system.simulate
            for i in self._index:
                for f, s in zip(setters, samples[i]): f(s)
                sim()
                add([i() for i in funcs])
        for k, v in zip(self._metrics.keys(), zip(*values)): self.table[k] = v
    
    def __call__(self, sample):
        """Return dictionary of metric values."""
        super().__call__(sample)
        return {i:j() for i,j in self._metrics.items()}
    
    def _repr(self):
        clsname = type(self).__name__
        newline = "\n" + " "*(len(clsname)+2)
        return f'{clsname}: {newline.join(self._metrics.keys())}'
        
    def __repr__(self):
        return f'<{type(self).__name__}: {", ".join(self._metrics.keys())}>'
    
        