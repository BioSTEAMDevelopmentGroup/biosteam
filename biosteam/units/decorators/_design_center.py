# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 21:18:50 2019

@author: yoelr
"""
from ..._stream import mol_flow_dim, mass_flow_dim, vol_flow_dim, _Q, DimensionError

__all__ = ('design_center',)

class Basis:
    __slots__ = ('factory', 'cached')
    def __init__(self, factory, cached):
        self.factory = factory
        self.cached = cached
        
    def __call__(self, args):
        cached = self.cached
        if args in cached:
            func = cached[args]
        else:
            func = self.factory(*args)
            cached[args] = func
        return func
    

class DesignCenter:
    
    def define(self, name, factory, cached=None):
        if name in self:
            raise ValueError(f"basis '{name}' already implemented")
        self.__dict__[name] = Basis(factory, cached or {})
    
    def __call__(self, name, args):
        try:
            basis = getattr(self, name)
        except:
            raise ValueError(f"unknown basis '{name}', basis must be one of the following: {', '.join(self)}")
        return basis(args)

    def __setattr__(self, name, basis):
        if isinstance(basis, Basis):
            self.__dict__[name] = basis
        else:
            raise ValueError(f"can only set 'DesignBasis' objects, not '{type(basis).__name__}'")

    def __contains__(self, basis):
        return basis in self.__dict__
    
    def __iter__(self):
        yield from self.__dict__
    
    def __repr__(self):
        return f"<{type(self).__name__}: {', '.join(self)}>"


design_center = DesignCenter()

def _flow(units, N_ins, N_outs):
    q = _Q(1, units)
    dim = q.dimensionality
    if dim == mol_flow_dim:
        if N_ins == 1:
            func = lambda self: self._ins[0].molnet
        elif N_outs == 1:
            func = lambda self: self._outs[0].molnet
        else: 
            func = lambda self: self._molnet_in
        ubase = 'kmol/hr'
    elif dim == mass_flow_dim:
        if N_ins == 1:
            func = lambda self: self._ins[0].massnet
        elif N_outs == 1:
            func = lambda self: self._outs[0].massnet
        else: 
            func = lambda self: self._massnet_in
        ubase = 'kg/hr'
    elif dim == vol_flow_dim:
        if N_ins == 1:
            func = lambda self: self._ins[0].volnet
        elif N_outs == 1:
            func = lambda self: self._outs[0].volnet
        else: 
            func = lambda self: self._volnet_in
        ubase = 'm3/hr'
    else:
        raise DimensionError(f"dimensions for flow units must be in molar, mass or volumetric flow rates, not '{dim}'")
    factor = 1/q.to(ubase).magnitude
    if factor == 1:
        return func
    else:
        return lambda self: factor*func(self)

design_center.define('Flow rate', _flow)

def _duty(units, N_ins, N_outs):
    if not (N_ins == N_outs == 1):
        raise ValueError(f"number of input and output streams must be 1 for selected basis")
    factor = _Q(1, 'kJ/hr').to(units)
    if factor == 1:
        return lambda self: self._outs[0].H - self._ins[0].H
    else:
        return lambda self: factor*(self._outs[0].H - self._ins[0].H)
    
design_center.define('Duty', _duty)

def _dry_flow(units, N_ins, N_outs):
    if not (N_ins == N_outs == 1):
        raise ValueError(f"number of input and output streams must be 1 for selected basis")
    if units != 'kg/hr':
        raise ValueError(f"units must be in kg/hr for selected basis")
    def dry_flow(self):
        feed = self._ins[0]
        return feed.massnet - feed.mass[feed.indices('7732-18-5', CAS=True)]
    return dry_flow

design_center.define('Dry flow rate', _dry_flow)