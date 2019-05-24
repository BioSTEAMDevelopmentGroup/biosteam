# -*- coding: utf-8 -*-
"""
Created on Thu May  9 13:38:57 2019

@author: Guest Group
"""
from ._block import Block
from .. import Unit, Stream
from ._name import elementname
import numpy as np

__all__ = ('State',)

# %% functions

def blockunit(blockf):
    element = blockf._element
    if isinstance(element, Unit): return element
    elif isinstance(element, Stream): return element._sink

def parameter(system, element, setter, kind, name):
    if kind is 'coupled':
        return Block(element, system).parameter(setter, name=name)
    elif kind is 'isolated':
        return Block(element, None).parameter(setter, name=name)
    elif kind is 'design':
        return Block(element, None).parameter(setter, element._summary, name)
    elif kind is 'cost':
        return Block(element, None).parameter(setter, element._finalize, name)
    raise ValueError(f"kind must be either 'coupled', 'isolated', 'design', or 'cost' (not {kind}).")

def modelfunction(params):
    cached = None
    zip_= zip
    params = tuple(params)
    def model(sample): 
        nonlocal cached
        try:
            sim = None
            for p, x, same in zip_(params, sample, cached==sample):
                if same: continue
                p._setter(x)
                if sim: continue
                if p._system: sim = p._simulate 
                else: p._simulate()
            if sim: sim()
            cached = sample
        except Exception as Error:
            cached = None
            raise Error
    return model

# %%
    
class State:
    """Create a State object that can update the `system` state given a set of parameters.
    
    **Parameters**
    
        **system:** [System] Reflects the model state.
        
    **Examples**
    
        :doc:`Model Example`
    
    """
    __slots__ = ('_system', # [System] Reflects the model state.
                 '_params', # list[Parameter] All parameters.
                 '_model')  # [function] Model function.
    
    def __init__(self, system):
        self._system = system
        self._params = []
        self._model = None
    
    def __len__(self):
        return len(self._params)
    
    def get_parameters(self):
        if not self._model: self._loadmodel()
        return tuple(self._params)
    
    def parameter(self, setter=None, element=None, kind='isolated', name=None):
        """Define parameter and add to model.
        
        **Parameters**
            
            **setter:** [function] Should set parameter in the element.
            
            **element:** [Unit or Stream] Element in the system being altered.
            
            **kind:** {'isolated', 'design', 'cost'}
                * 'coupled': parameter is coupled to the system.
                * 'isolated': parameter does not affect the system but does affect the element (if any).
                * 'design': parameter only affects design and cost of the element.
                * 'cost': parameter only affects cost of the element.
                
            **param:** [str] Name of parameter. If None, default to argument name of setter.
            
        .. Note::
            
            If kind is 'coupled', account for downstream operations. Otherwise, only account for given element. If kind is 'design' or 'cost', element must be a Unit object.
        
        """
        if not setter: return lambda setter: self.parameter(setter, element, kind, name)
        param = parameter(self._system, element, setter, kind, name)
        self._params.append(param)
        self._erase()
        return param
    
    def _erase(self):
        """Erase cached data."""
        self._model = None
    
    def _loadmodel(self):
        length = len( self._system._unitnetwork)
        index =  self._system._unitnetwork.index
        self._params.sort(key=lambda x: index(blockunit(x))
                                        if x._system else length)
        self._model = modelfunction(self._params)
    
    def __call__(self, sample):
        if not self._model: self._loadmodel()
        return self._model(np.asarray(sample))
    
    def _repr(self):
        return f'{type(self).__name__}: {self._system}'
    
    def __repr__(self):
        return '<' + self._repr() + '>'
       
    def _info(self):
        if not self._model: self._loadmodel()
        if not self._params:
            return (f'{self._repr()}\n'
                    +' Element:  Parameters:\n'
                    +'  None      None')
        lines = []
        lenghts_block = []
        lastblk = None
        for i in self._params:
            blk = elementname(i._element)
            element = len(blk)*' ' if blk==lastblk else blk
            lines.append(f"  {element}${i._name}\n")
            lastblk = blk
            lenghts_block.append(len(blk))
        maxlen_block = max(lenghts_block)
        out = f'{self._repr()}\n'
        maxlen_block = max(maxlen_block, 8)
        out += ' Element:' + (maxlen_block - 8)*' ' + '  Parameter:\n'
        for newline, len_ in zip(lines, lenghts_block):
            newline = newline.replace('$', ' '*(maxlen_block-len_) + '  ')
            out += newline
        return out.rstrip('\n ')
    
    def show(self):
        """Return information on metric parameters."""
        print(self._info())
    _ipython_display_ = show
    