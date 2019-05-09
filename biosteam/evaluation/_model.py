# -*- coding: utf-8 -*-
"""
Created on Thu May  9 13:38:57 2019

@author: Guest Group
"""
from ._block import Block
from .. import Unit, Stream
import numpy as np

__all__ = ('Model',)


# %% functions

_pretty = lambda param: param.replace('_', ' ').capitalize()

def _get_element_name(element):
    if element:
        if isinstance(element, type):
            return element.__name__.replace('_', ' ')
        elif isinstance(element, str):
            return element.replace('_', ' ')
        else:
            return element.line + '-' + element.ID.replace('_', ' ')
    else:
        return 'None'

def _blockunit(blockf):
    element = blockf._element
    if isinstance(element, Unit): return element
    elif isinstance(element, Stream): return element._sink


# %%
    
class Model:
    
    __slots__ = ('_system', '_metric', '_blockf', '_params', '_orderf', '_ID')
    
    def __init__(self, ID, system, metric, *, args=()):
        self._ID = ID
        self._system = system
        self._metric = (lambda: metric(*args)) if args else metric
        self._blockf = []
        self._params = None
        self._orderf =  None
    
    def add(self, setter, kind='isolated', element=None, param=None):
        """Add parameter to model.
        
        **Parameters**
        
            **element:** [Unit or Stream] Element in the system being altered.
            
            **setter:** [function] Should set parameter in the element.
            
            **kind:** {'isolated', 'design', 'cost'}
                * 'coupled': parameter is coupled to the system.
                * 'isolated': parameter does not affect the system but does affect the element (if any).
                * 'design': parameter only affects design and cost of the element.
                * 'cost': parameter only affects cost of the element.
            
        .. Note::
            
            If kind is 'coupled', account for downstream operations. Otherwise, only account for given element. If kind is 'design' or 'cost', element must be a Unit object.
        
        """
        if kind is 'coupled':
            blockfunc = Block(element, self._system)(setter, param=param)
        elif kind is 'isolated':
            blockfunc = Block(element, None)(setter, param=param)
        elif kind is 'design':
            blockfunc = Block(element, None)(setter, element._summary, param=param)
        elif kind is 'cost':
            blockfunc = Block(element, None)(setter, element._finalize, param=param)
        if kind not in ('coupled', 'isolated', 'design', 'cost'):
            raise ValueError(f"kind must be either 'coupled', 'isolated', 'design', or 'cost' (not {kind}).")
        self._blockf.append(blockfunc)
        self._optimized = False
    
    def _optimize_param_order(self):
        if self._orderf != self._blockf:    
            system = self._system
            length = len(system._unitnetwork)
            index = system._unitnetwork.index
            self._orderf = sorted(self._blockf, key=lambda x: index(_blockunit(x))
                                                              if x._system else length,
                                  reverse=True)
    
    def __call__(self, parameters):
        self._optimize_param_order()
        parameters = np.asarray(parameters)
        try:
            reset = self._params!=parameters
        except:
            reset = np.ones_like(parameters)
        sim = True
        for i, j, k in zip(reset, self._orderf, parameters):
            if i:
                j._setter(k)
                if sim: j._simulate()
                else: sim = False
        self._params = parameters
        return self._metric()
    
    def _repr(self):    
        return f'{type(self).__name__}: {self._ID}'
    
    def __repr__(self):
        return f'<{self._repr()}>'
       
    def _info(self):
        self._optimize_param_order()
        if not self._orderf:
            return (f'{self._repr()}\n'
                    +' Element  Parameters\n'
                    +' None     None')
        lines = []
        lenghts_block = []
        lastblk = None
        for i in self._orderf:
            p = i._param
            blk = _get_element_name(i._element)
            element = len(blk)*' ' if blk==lastblk else blk
            lines.append(f" {element}${_pretty(p)}\n")
            lastblk = blk
            lenghts_block.append(len(blk))
        maxlen_block = max(lenghts_block)
        out = f'{self._repr()}\n'
        maxlen_block = max(maxlen_block, 7)
        out += ' Element' + (maxlen_block - 7)*' ' + '  Parameter\n'
        for newline, len_ in zip(lines, lenghts_block):
            newline = newline.replace('$', ' '*(maxlen_block-len_) + '  ')
            out += newline
        return out.rstrip('\n ')
    
    def show(self):
        """Return information on metric parameters."""
        print(self._info())
    
    
    # def create(self, name):
    #     """Return a model function."""
    #     sta
    #     params = [i for i in self._stack]
        
    #     name = getter.__name__
    #     if name[0] == '<': name = 'Lambda'
    #     str2exec =  (f'def {name}({param}):    \n'
    #                + f'    setter({param}) \n'
    #                + f'    simulate()      \n'
    #                + f'    return getter()   ')

    #     globs = {'setter': setter,
    #              'simulate': simulate,
    #              'getter': getter}
    #     locs = {}
    #     exec(str2exec, globs, locs)
    #     blockfunc = locs[name]
    #     blockfunc.__qualname__ = f'{self} {blockfunc.__qualname__}'
    #     blockfunc._param = param
    #     blockfunc._simulate = simulate
    #     blockfunc._element = self._element
    #     blockfunc._system = self._system
    #     blockfunc._setter = setter
    #     blockfunc._getter = getter
    #     return blockfunc