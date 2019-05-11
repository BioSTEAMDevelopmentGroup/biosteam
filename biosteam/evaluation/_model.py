# -*- coding: utf-8 -*-
"""
Created on Thu May  9 13:38:57 2019

@author: Guest Group
"""
from ._block import Block
from .. import Unit, Stream
import numpy as np
import re

__all__ = ('Model',)

# %% functions

_pretty = lambda param: param.replace('_', ' ').capitalize()

def elementname(element):
    if element:
        if isinstance(element, type):
            return re.sub(r"\B([A-Z])", r" \1", element.__name__.replace('_', ' ')).capitalize()
        elif isinstance(element, str):
            return element.replace('_', ' ')
        else:
            return element.line + '-' + element.ID.replace('_', ' ')
    else:
        return 'None'

def blockunit(blockf):
    element = blockf._element
    if isinstance(element, Unit): return element
    elif isinstance(element, Stream): return element._sink

def blockfunction(system, element, setter, kind, param):
    if kind is 'coupled':
        return Block(element, system)(setter, param=param)
    elif kind is 'isolated':
        return Block(element, None)(setter, param=param)
    elif kind is 'design':
        return Block(element, None)(setter, element._summary, param=param)
    elif kind is 'cost':
        return Block(element, None)(setter, element._finalize, param=param)
    raise ValueError(f"kind must be either 'coupled', 'isolated', 'design', or 'cost' (not {kind}).")

# %%
    
class Model:
    
    __slots__ = ('_system', # [System] Reflects the model state.
                 '_blockf', # list[function] All block functions in the stack.
                 '_reload', # [bool] True if model should be reloaded.
                 '_params') # [list] Cached parameter values of last model run.
    
    def __init__(self, system):
        self._system = system
        self._blockf = []
        self._params = None
    
    def __len__(self):
        return len(self._blockf)
    
    def addparam(self, setter=None, element=None, kind='isolated', param=None):
        """Add parameter to model.
        
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
        if not setter: return lambda setter: self.addparam(setter, element, kind, param)
        f = blockfunction(self._system, element, setter, kind, param)
        self._blockf.append(f)
        self._reload = True
        return f
    
    def _loadmodel(self):
        if self._reload:
            system = self._system
            length = len(system._unitnetwork)
            index = system._unitnetwork.index
            self._blockf.sort(key=lambda x: index(blockunit(x))
                                            if x._system else length,
                              reverse=True)
            self._reload = False
    
    def __call__(self, parameters, default=None):
        self._loadmodel()
        parameters = np.asarray(parameters)
        try: same = self._params==parameters
        except: same = np.zeros_like(parameters)
        sim = True
        for s, f, p in zip(same, self._blockf, parameters):
            if s: continue
            try:
                f._setter(p)
                if sim: f._simulate()
                elif f._system: sim = False
            except: return default
        self._params = parameters
    
    def _repr(self):
        ID = re.sub(r"\B([A-Z])", r" \1", self._system.ID.replace('_', ' ')).capitalize()
        return f'{type(self).__name__}: {ID}'
    
    def __repr__(self):
        sig = ', '.join(i._param for i in self._blockf)
        return f'<{self._system.ID}({sig})>'
       
    def _info(self):
        self._loadmodel()
        if not self._blockf:
            return (f'{self._repr()}\n'
                    +' Element:  Parameters:\n'
                    +'  None      None')
        lines = []
        lenghts_block = []
        lastblk = None
        for i in self._blockf:
            p = i._param
            blk = elementname(i._element)
            element = len(blk)*' ' if blk==lastblk else blk
            lines.append(f"  {element}${_pretty(p)}\n")
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

    