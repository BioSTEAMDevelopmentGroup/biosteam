# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import biosteam as bst
import copy
from ._design import design
from math import ceil

__all__ = ('cost', 'copy_algorithm', 'add_cost', 'CostItem')

class CostItem:
    """
    Create a CostItem object which defines exponential scaling factors for an item's purchase cost.
    
    Parameters
    ----------
    basis : str
        Name of size parameter used for scaling.
    units : str
        Units of measure.
    S : float
        Size.
    lb : float, optional
        Lower size bound.
    ub : float
        Upper Size bound.
    CE : float
        Chemical engineering plant cost index.
    cost : float
        Purchase cost of item.
    n : float
        Exponential factor.
    kW : float
        Electricity rate.
    N : str
        Attribute name for number of parallel units.
    f : function, optional
        Should return the cost given the size `S` at the given CE.
    
    """
    __slots__ = ('_basis', '_units', 'S', 'lb', 'ub', 'CE',
                 'cost', 'n', 'kW', 'N', 'f')
    def __init__(self, basis, units, S, lb, ub, CE, cost, n, kW, N, f):
        if f:
            if (cost is not None or n is not None):
                raise ValueError('cannot define parameters `cost` and `n` '
                                 'when the cost function `f` is given')
            self.n = self.cost = None
        else:
            self.cost = 0. if cost is None else float(cost)
            self.n = 1. if n is None else float(n)
        if lb is not None: lb = float(lb)
        if ub is not None: ub = float(ub)
        if N:
            if not isinstance(N, str): # Prevent downstream error for common mistakes
                raise ValueError("N parameter must be a string or None; not a "
                                 "'{type(N).__name__}' object")
        elif ub is not None:
            N = '#'
        else:
            N = None
        self._basis = str(basis)
        self._units = str(units)
        self.S = float(S)
        self.lb = lb
        self.ub = ub
        self.CE = float(CE)
        self.kW = 0. if kW is None else float(kW)
        self.N = N
        self.f = f
    
    __getitem__ = object.__getattribute__
    __setitem__ = object.__setattr__
    
    def copy(self):
        new = CostItem.__new__(CostItem)
        new._basis = self._basis
        new._units = self._units
        new.S = self.S
        new.lb = self.lb
        new.ub = self.ub
        new.CE = self.CE
        new.cost = self.cost
        new.n = self.n
        new.kW = self.kW
        new.N = self.N
        new.f = self.f
        return new
    
    @property
    def basis(self):
        return self._basis
    @property
    def units(self):
        return self._units
    
    def __repr__(self):
        return f"<{type(self).__name__}: {self._basis} ({self._units})>"
    
    def _ipython_display_(self):
        print(f"{type(self).__name__}: {self._basis} ({self._units})\n"
             +f" S     {self.S:.3g}\n"
            +(f" lb    {self.lb:.3g}\n" if self.lb is not None else "")
            +(f" ub    {self.ub:.3g}\n" if self.ub is not None else "")
             +f" CE    {self.CE:.3g}\n"
            +(f" cost  {self.cost:.3g}\n" if self.cost is not None else "")
            +(f" n     {self.n:.3g}\n"  if self.n is not None else "")
             +f" kW    {self.kW:.3g}\n"
            +(f" N     '{self.N}'" if self.N is not None else ""))
    show = _ipython_display_

def _decorated_cost(self):
    D = self.design_results
    C = self.baseline_purchase_costs
    kW = 0
    for i, x in self.cost_items.items():
        S = D[x._basis]
        if x.lb is not None and S < x.lb:
            S = x.lb
            D[x.N or '#' + i] = 1
        elif x.ub is not None:
            D[x.N or '#' + i] = N = ceil(S/x.ub)
            q = S/x.S
            F = q/N
            C[i] = N*bst.CE/x.CE*(x.f(F) if x.f else x.cost*F**x.n)
            kW += x.kW*q
            continue
        if x.N:
            N = getattr(self, x.N, None) or D[x.N]
            F = S/x.S
            C[i] = N*bst.CE/x.CE*(x.f(F) if x.f else x.cost*F**x.n)
            kW += N*x.kW*F
        else:
            F = S/x.S
            C[i] = bst.CE/x.CE*(x.f(F) if x.f else x.cost*abs(F)**x.n)
            kW += x.kW*F
    if kW: self.power_utility(kW)

def copy_algorithm(other, cls=None, run=True, design=True, cost=True):
    if not cls: return lambda cls: copy_algorithm(other, cls, run, design, cost)
    dct = cls.__dict__
    if run:
        if '_run' in dct: raise RuntimeError('run method already implemented')
        cls._run = other._run
    if cost:
        if '_cost' in dct: raise RuntimeError('cost method already implemented')
        cls._cost = other._cost
        try:
            cls._F_BM_default = other._F_BM_default
            cls.cost_items = other.cost_items
        except: pass
    if design:
        if '_design' in dct: raise RuntimeError('design method already implemented')
        cls._design = other._design
        cls._units = other._units
        try:
            cls._decorated_cost = other._decorated_cost
            cls._design_basis_ = other._design_basis_
            cls._decorated_design = other._decorated_design
        except: pass
    return cls

def cost(basis, ID=None, *, CE, cost=None, n=None, S=1., lb=None, ub=None, kW=None, BM=1.,
         units=None, N=None, lifetime=None, f=None):    
    r"""
    Add item (free-on-board) purchase cost based on exponential scale up.
    
    Parameters
    ----------
    basis : str
        Name of size parameter used for scaling.    
    ID : str
        Name of purchase item.
    CE : float
         Chemical engineering plant cost index.
    cost : float
        Purchase cost of item.
    n : float
        Exponential factor.
    S : float, optional
        Size. Defaults to 1.
    ub : float, optional
        Size limit, if any.
    kW : float, optional
        Electricity rate. Defaults to 0.
    BM : float, optional
        Bare module factor (installation factor). Defaults to 1.
    units : str, optional
        Units of measure.
    N : str, optional
        Attribute name for number of parallel units.
    lifetime : int, optional
        Number of operating years until equipment needs to be replaced.
    f : function, optional
        Should return the cost given the size `S` at the given `CE`.
        
    Examples
    --------
    :doc:`../../tutorial/Unit_decorators`
    
    """
    
    return lambda cls: add_cost(cls, ID, basis, units, S, lb, ub, CE, cost, n, kW, BM, N, lifetime, f)

def add_cost(cls, ID, basis, units, S, lb, ub, CE, cost, n, kW, BM, N, lifetime, f):
    # Make sure new _units dictionary is defined
    if '_units' not in cls.__dict__:
        cls._units = cls._units.copy() if hasattr(cls, '_units') else {}
    if basis in cls._units:
        if not units: 
            units = cls._units[basis]
        elif units != cls._units[basis]:
            raise RuntimeError(f"cost basis '{basis}' already defined in class with units '{cls._units[basis]}'")
    elif units:
        design.add_design_basis_to_cls(cls, basis, units)
    else:
        raise RuntimeError(f"units of cost basis '{basis}' not available in '{cls.__name__}._units' dictionary, must pass units or define in class")
    if hasattr(cls, 'cost_items'):
        if 'cost_items' not in cls.__dict__:
            cls.cost_items = copy.deepcopy(cls.cost_items)
        if not ID:
            raise ValueError("must pass an 'ID' for purchase cost item")
        if ID in cls.cost_items:
            raise ValueError(f"ID '{ID}' already in use")
        cls.cost_items[ID] = CostItem(basis, units, S, lb, ub, CE, cost, n, kW, N, f)
        cls._F_BM_default[ID] = BM
        if lifetime: cls._default_equipment_lifetime[ID] = lifetime
    else:
        ID = ID or cls.line
        cls.cost_items = {ID: CostItem(basis, units, S, lb, ub, CE, cost, n, kW, N, f)}
        if '_F_BM_default' not in cls.__dict__:
            cls._F_BM_default = cls._F_BM_default.copy() if hasattr(cls, '_F_BM_default') else {}
        if '_default_equipment_lifetime' not in cls.__dict__:
            cls._default_equipment_lifetime = cls._default_equipment_lifetime.copy() if hasattr(cls, '_default_equipment_lifetime') else {}
        cls._F_BM_default[ID] = BM
        if lifetime: cls._default_equipment_lifetime[ID] = lifetime
        if '_cost' not in cls.__dict__:
            cls._cost = _decorated_cost
        cls._decorated_cost = _decorated_cost
    return cls

