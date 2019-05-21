# -*- coding: utf-8 -*-
"""
Created on Mon May 20 22:04:02 2019

@author: yoelr
"""

from ... import Stream, Unit
from ..._unit import metaUnit

__all__ = ('splitter', '_kwargs_split', '_init_split', '_run_split')

class splitter(metaUnit):
    def __new__(mcl, clsname, superclasses, definitions):
        if len(superclasses) > 1:
            raise RuntimeError('cannot create splitter Unit subclass with more than one super class')
        kwargs = _kwargs_split.copy()
        kwargs.update(definitions.get('_kwargs') or superclasses[0]._kwargs.copy())
        definitions['_kwargs'] = kwargs
        _init = definitions.get('_init') or superclasses[0]._init
        definitions['_init'] = _init_split if _init is Unit._init else _merge(_init, _init_split)
        definitions['reset'] = _reset_split
        definitions['_run'] = definitions.get('_run') or _run_split
        return metaUnit(clsname, superclasses, definitions)
    
def _merge(first, last):
    return lambda self: first(self) or last(self)

_kwargs_split = {'split': None,
                 'order': None}

def _init_split(self):
    self._reorder = Stream._cls_species._reorder
    _kwargs = self._kwargs
    order = _kwargs.pop('order')
    if order:
        _kwargs['split'] = self._reorder(_kwargs['split'], order)

def _reset_split(self, split=None, order=None, **kwargs):
    _kwargs = self._kwargs
    _kwargs['split'] = self._reorder(split, order) if order else split
    if kwargs:
        _kwargs.update(kwargs)
        self._setup()
    
def _run_split(self):
    split = self._kwargs['split']
    top, bot = self.outs
    feed = self._ins[0]
    net_mol = feed.mol
    top._mol[:] = net_mol*split
    bot._mol[:] = net_mol - top._mol
    bot.T = top.T = feed.T
    bot.P = top.P = feed.P
    