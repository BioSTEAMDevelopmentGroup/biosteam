# -*- coding: utf-8 -*-
"""
Created on Thu May  9 13:38:57 2019

@author: Guest Group
"""
from ._block import Block
from .. import Unit, Stream
import numpy as np
import pandas as pd
from chaospy import J

__all__ = ('State',)

# %% functions

def param_unit(param):
    element = param.element
    if isinstance(element, Unit): return element
    elif isinstance(element, Stream): return element._sink

def parameter(system, element, setter, kind, name, distribution, units):
    if kind is 'coupled':
        return Block(element, system).parameter(setter, name=name, distribution=distribution, units=units)
    elif kind is 'isolated':
        return Block(element, None).parameter(setter, name=name, distribution=distribution, units=units)
    elif kind is 'design':
        return Block(element, None).parameter(setter, element._summary, name, distribution=distribution, units=units)
    elif kind is 'cost':
        if hasattr(element, '_end'):
            def simulate():
                element._cost()
                element._end()
        else:
            simulate = element._cost
        return Block(element, None).parameter(setter, simulate, name, distribution=distribution, units=units)
    raise ValueError(f"kind must be either 'coupled', 'isolated', 'design', or 'cost' (not {kind}).")

def modelfunction_with_skipping(params):
    cached = None
    zip_= zip
    params = tuple(params)
    def model(sample): 
        nonlocal cached
        try:
            sim = None
            for p, x, same in zip_(params, sample, cached==sample):
                if same: continue
                p.setter(x)
                if sim: continue
                if p.system: sim = p.simulate 
                else: p.simulate()
            if sim: sim()
            cached = sample
        except Exception as Error:
            cached = None
            raise Error
    return model

def modelfunction(params):
    zip_= zip
    params = tuple(params)
    def model(sample): 
        sim = None
        for p, x in zip_(params, sample):
            p.setter(x)
            if sim: continue
            if p.system: sim = p.simulate 
            else: p.simulate()
        if sim: sim()
    return model


# %%
    
class State:
    """Create a State object that can update the `system` state given a sample of parameter states.
    
    **Parameters**
    
        **system:** [System] Reflects the model state.
        
        **skip:** [bool] If True, skip simulation for repeated states
    
    """
    __slots__ = ('_system', # [System] Reflects the model state.
                 '_params', # list[Parameter] All parameters.
                 '_model', # [function] Model function.
                 '_skip') # [bool] If True, skip simulation for repeated states
    
    def __init__(self, system, skip=True):
        self._system = system
        self._params = []
        self._model = None
        self._skip = skip
    
    def __len__(self):
        return len(self._params)
    
    def get_parameters(self):
        if not self._model: self._loadmodel()
        return tuple(self._params)
    
    def get_distribution_summary(self):
        params = self.get_parameters()
        if not params: return None
        index = []
        values = []
        for p in params:
            distribution = p.distribution
            element = p.element_name
            name = p.name
            units = p.units or ''
            if distribution:
                shape = type(distribution).__name__
                index.append((element, name, units, shape))
                values.append(', '.join([f'{i:.4g}' for i in distribution._repr.values()]))
            else:
                index.append((element, name, units, ''))
                values.append('')
        index = pd.MultiIndex.from_tuples(index,
                                          names=('Element',
                                                 'Name',
                                                 'Units',
                                                 'Shape'))
        return pd.DataFrame(values, index, ('Values',))
    
    def parameter(self, setter=None, element=None, kind='isolated',
                  name=None, distribution=None, units=None):
        """Define parameter and add to model.
        
        **Parameters**
            
            **setter:** [function] Should set parameter in the element.
            
            **element:** [Unit or Stream] Element in the system being altered.
            
            **kind:** {'isolated', 'design', 'cost'}
                * 'coupled': parameter is coupled to the system.
                * 'isolated': parameter does not affect the system but does affect the element (if any).
                * 'design': parameter only affects design and cost of the element.
                * 'cost': parameter only affects cost of the element.
                
            **name:** [str] Name of parameter. If None, default to argument name of setter.
            
            **distribution:** [chaospy.Dist] Parameter distribution.
            
            **units:** [str] Parameter units
            
        .. Note::
            
            If kind is 'coupled', account for downstream operations. Otherwise, only account for given element. If kind is 'design' or 'cost', element must be a Unit object.
        
        """
        if not setter: return lambda setter: self.parameter(setter, element, kind, name, distribution, units)
        param = parameter(self._system, element, setter, kind, name, distribution, units)
        self._params.append(param)
        self._erase()
        return param
    
    def sample(self, N, rule): 
        """Sample from parameter distribution.
        
        **Parameters**
        
            **N:** [int] Number of samples.
            
            **rule:** [str] Sampling rule (e.g. "L" for latin hypercube sampling).
        """
        if not self._model: self._loadmodel()
        return J(*[i.distribution for i in self._params]).sample(N, rule).transpose()
    
    def _erase(self):
        """Erase cached data."""
        self._model = None
    
    def _loadmodel(self):
        length = len( self._system._unitnetwork)
        index =  self._system._unitnetwork.index
        self._params.sort(key=lambda x: index(param_unit(x))
                                        if x.system else length)
        if self._skip:
            self._model = modelfunction_with_skipping(self._params)
        else:
            self._model = modelfunction(self._params)
    
    def __call__(self, sample):
        """Update state given sample of parameters."""
        if not self._model: self._loadmodel()
        return self._model(np.asarray(sample))
    
    def _repr(self):
        return f'{type(self).__name__}: {self._system}'
    
    def __repr__(self):
        return '<' + self._repr() + '>'
    
    def _info(self):
        if not self._model: self._loadmodel()
        if not self._params: return f'{self._repr()}\n (No parameters)'
        lines = []
        lenghts_block = []
        lastblk = None
        for i in self._params:
            blk = i.element_name
            element = len(blk)*' ' if blk==lastblk else blk
            lines.append(f" {element}${i.name}\n")
            lastblk = blk
            lenghts_block.append(len(blk))
        maxlen_block = max(lenghts_block)
        out = f'{self._repr()}\n'
        maxlen_block = max(maxlen_block, 7)
        out += ' Element:' + (maxlen_block - 7)*' ' + ' Parameter:\n'
        for newline, len_ in zip(lines, lenghts_block):
            newline = newline.replace('$', ' '*(maxlen_block-len_) + '  ')
            out += newline
        return out.rstrip('\n ')
    
    def show(self):
        """Return information on metric parameters."""
        print(self._info())
    _ipython_display_ = show
    