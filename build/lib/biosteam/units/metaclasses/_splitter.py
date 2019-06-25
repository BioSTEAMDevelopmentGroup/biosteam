# -*- coding: utf-8 -*-
"""
Created on Mon May 20 22:04:02 2019

@author: yoelr
"""

from ... import Stream
from ..._unit import metaUnit
from numpy import asarray
from pandas import Series

__all__ = ('splitter', '__init_split__', 'run_split', 'run_split_with_mixing')

class splitter(metaUnit):
    """Create a splitter Unit class which behaves like a splitter regarding mass and energy balances.
    
    **Examples**
    
    Create a sugar cane crushing mill Unit subclass and test:
    
    >>> from biosteam import *
    >>> from biosteam.units.decorators import *
    >>> from biosteam.units.metaclasses import splitter
    >>> # Create subclass
    >>> @cost('Flow rate', units='kg/hr', cost=1.5e6,
    ...       CE=541.7, exp=0.6, S=335e3, kW=2010)
    ... class CrushingMill(Unit, metaclass=splitter): pass
    >>> # Set Stream.species
    >>> species = Species('Water')
    >>> species.Sugar = compounds.Substance('Sugar')
    >>> species.Fiber = compounds.Substance('Fiber')
    >>> Stream.species = species
    >>> # Create CrushingMill object and simulate
    >>> CM = CrushingMill(ins=Stream(Water=700, Fiber=150, Sugar=150, units='kg/min'), 
    ...                   split=(0.86, 0.08, 0.96),
    ...                   order=('Water', 'Fiber', 'Sugar'))
    >>> CM.simulate()
    >>> CM.show()
    CrushingMill: U1
    ins...
    [0] d2
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Water  2.33e+03
                        Sugar  9e+03
                        Fiber  9e+03
    outs...
    [0] d3
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Water  2e+03
                        Sugar  8.64e+03
                        Fiber  720
    [1] d4
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Water  326
                        Sugar  360
                        Fiber  8.28e+03
    
    """
    def __new__(mcl, name, bases, dct):
        if '__init__' in dct:
            raise TypeError(f"cannot use {mcl.__name__} metaclass with an implemented '__init__' method")
        elif 'split' in dct:
            raise TypeError("cannot use {mcl.__name__} metaclass with an implemented 'split' member")
        
        _kwargs = dct.get('_kwargs')
        if _kwargs is None:
            dct['__init__'] = __init_split__
        elif 'split' in _kwargs or 'order' in _kwargs:
            raise TypeError("cannot use {mcl.__name__} metaclass with 'split' or 'order' defined in '{clsname}._kwargs'")
        else:
            # Begin changing __init__ to have kwargs
            inputs = ', '.join([key + '=' + key for key in _kwargs])
            
            # Make a string to execute
            str2exec = f"def __init__(self, ID='', outs=(), ins=None, order=None, {inputs}, *, split):\n"
            str2exec+= f"    self._kwargs = dict({inputs})\n"
            str2exec+= f"    _(self, ID, outs, ins, order, split=split)"
            
            # Execute string and replace __init__
            globs = {'_': __init_split__}
            globs.update(_kwargs)
            exec(str2exec, globs, dct)
        
        if '_run' not in dct:
            dct['_run'] = run_split
        
        if '_N_ins' not in dct: dct['_N_ins'] = 1
        dct['_N_outs'] = 2
        dct['split'] = splitprop
        
        if "__doc__" in dct:
            dct['__doc__'] = dct['__doc__'].replace('**Parameters**', 
       "**Parameters**"
+"\n"    
+"\n        **split:** Should be one of the following"
+"\n            * [float] The fraction of net feed in the 0th output stream"
+"\n            * [array_like] Componentwise split of feed to 0th output stream"
+"\n"
+"\n        **order:** Iterable[str] Species order of split")    
        
        return super().__new__(mcl, name, bases, dct)

@property
def splitprop(self):
    """[array] Componentwise split of feed to 0th output stream."""
    try:
        return self._split_series
    except:
        self._split_series = split = Series(self._split,
                                            index=self._reorder_.__self__._IDs)
        return split

def __init_split__(self, ID='', outs=(), ins=None, order=None, *, split):
    try: self._reorder_ = Stream._cls_species._reorder
    except AttributeError as err:
        if Stream._cls_species: raise err
        else: raise RuntimeError('must set Stream.species first')
    self.ID = ID
    self._split = self._reorder_(split, order) if order else asarray(split)
    self._init_ins(ins)
    self._init_outs(outs)
    self._init_results()
    self._init_heat_utils()
    self._init_power_util()
    self._init()
    self._setup()
    
def run_split(self):
    """Splitter mass and energy balance function based on one input stream."""
    top, bot = self.outs
    feed = self._ins[0]
    net_mol = feed.mol
    top._mol[:] = net_mol * self._split
    bot._mol[:] = net_mol - top._mol
    bot.T = top.T = feed.T
    bot.P = top.P = feed.P
    
def run_split_with_mixing(self):
    """Splitter mass and energy balance function with mixing all input streams."""
    top, bot = self._outs
    ins = self._ins
    if len(ins) > 1: Stream.sum(top, ins)
    else: top.copylike(ins[0])
    bot.copylike(top)
    top._mol[:] *= self._split
    bot._mol[:] -= top._mol