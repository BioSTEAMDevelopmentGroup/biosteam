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
    ub : float
         Size limit.
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
    
    """
    __slots__ = ('_basis', '_units', 'S', 'ub', 'CE',
                 'cost', 'n', 'kW', 'N')
    def __init__(self, basis, units, S, ub, CE, cost, n, kW, N):
        s = str; f = float; b = bool
        if N: # Prevent downstream error for common mistakes
            if isinstance(N, s):
                self.N = N
            else:
                raise ValueError("N parameter must be a string or None; not a "
                                 "'{type(N).__name__}' object")
        elif ub:
            self.N = '#'
        else:
            self.N = None
        self._basis = s(basis)
        self._units = s(units)
        self.S = f(S)
        self.ub = f(ub)
        self.CE = f(CE)
        self.cost = f(cost)
        self.n = f(n)
        self.kW = f(kW)
        self.N = N
    
    __getitem__ = object.__getattribute__
    __setitem__ = object.__setattr__
    
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
             +f" ub    {self.ub:.3g}\n"
             +f" CE    {self.CE:.3g}\n"
             +f" cost  {self.cost:.3g}\n"
             +f" n     {self.n:.3g}\n"
             +f" kW    {self.kW:.3g}\n"
            +(f" N     '{self.N}'" if self.N else ""))
    show = _ipython_display_

def _decorated_cost(self):
    D = self.design_results
    C = self.purchase_costs
    kW = 0
    for i, x in self.cost_items.items():
        S = D[x._basis]
        if x.ub:
            D[x.N or '#' + i] = N = ceil(S/x.ub)
            q = S/x.S
            F = q/N
            C[i] = N*bst.CE/x.CE*x.cost*F**x.n
            kW += x.kW*q
        elif x.N:
            N = getattr(self, x.N, None) or D[x.N]
            F = S/x.S
            C[i] = N*bst.CE/x.CE*x.cost*F**x.n
            kW += N*x.kW*F
        else:
            F = S/x.S
            C[i] = bst.CE/x.CE*x.cost*F**x.n
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
            cls._BM = other._BM
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

def cost(basis, ID=None, *, CE, cost, n, S=1., ub=0., kW=0., BM=1.,
         units=None, N=None, lifetime=None, fsize=None):    
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
    fsize : function, optional
        Accepts a Unit object argument and returns the size parameter. If not 
        given, defaults to function predefined for given name and units of 
        measure.
        
    Examples
    --------
    :doc:`../../tutorial/Unit_decorators`
    
    """
    
    return lambda cls: add_cost(cls, ID, basis, units, S, ub, CE, cost, n, kW, BM, N, lifetime, fsize)

def add_cost(cls, ID, basis, units, S, ub, CE, cost, n, kW, BM, N, lifetime, fsize):
    # Make sure new _units dictionary is defined
    if '_units' not in cls.__dict__:
        cls._units = cls._units.copy() if hasattr(cls, '_units') else {}
    if basis in cls._units:
        if fsize:
            raise RuntimeError(f"cost basis '{basis}' already defined in class, cannot pass 'fsize' argument")
        elif not units: 
            units = cls._units[basis]
        elif units != cls._units[basis]:
            raise RuntimeError(f"cost basis '{basis}' already defined in class with units '{cls._units[basis]}'")
    elif units:
        design.add_design_basis_to_cls(cls, basis, units, fsize)
    else:
        raise RuntimeError(f"units of cost basis '{basis}' not available in '{cls.__name__}._units' dictionary, must pass units or define in class")
    if hasattr(cls, 'cost_items'):
        if 'cost_items' not in cls.__dict__:
            cls.cost_items = copy.deepcopy(cls.cost_items)
        if not ID:
            raise ValueError("must pass an 'ID' for purchase cost item")
        if ID in cls.cost_items:
            raise ValueError(f"ID '{ID}' already in use")
        cls.cost_items[ID] = CostItem(basis, units, S, ub, CE, cost, n, kW, N)
        cls._BM[ID] = BM
        if lifetime: cls._equipment_lifetime[ID] = lifetime
    else:
        ID = ID or cls.line
        cls.cost_items = {ID: CostItem(basis, units, S, ub, CE, cost, n, kW, N)}
        if '_BM' not in cls.__dict__:
            cls._BM = cls._BM.copy() if hasattr(cls, '_BM') else {}
        if '_equipment_lifetime' not in cls.__dict__:
            cls._equipment_lifetime = cls._equipment_lifetime.copy() if hasattr(cls, '_equipment_lifetime') else {}
        cls._BM[ID] = BM
        if lifetime: cls._equipment_lifetime[ID] = lifetime
        if '_cost' not in cls.__dict__:
            cls._cost = _decorated_cost
        cls._decorated_cost = _decorated_cost
    return cls

