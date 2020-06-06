# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
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

__all__ = ('cost', 'add_cost', 'CostItem')

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
        if ub and not N: N='#'
        self._basis = basis
        self._units = units
        self.S = S
        self.ub = ub
        self.CE = CE
        self.cost = cost
        self.n = n
        self.kW = kW
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
            N = getattr(self, x.N)
            F = S/x.S
            C[i] = N*bst.CE/x.CE*x.cost*F**x.n
            kW += N*x.kW*F
        else:
            F = S/x.S
            C[i] = bst.CE/x.CE*x.cost*F**x.n
            kW += x.kW*F
    if kW: self.power_utility(kW)

def cost(basis, ID=None, *, CE, cost, n, S=1, ub=0, kW=0, BM=1,
         units=None, fsize=None, N=None):    
    r"""
    Add item (free-on-board) purchase cost based on exponential scale up:
    
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
    S = 1 : float
        Size.
    ub : float
         Size limit.
    kW : float
        Electricity rate.
    BM = 1: float
        Bare module factor (installation factor).
    units = None : str, optional
        Units of measure.
    fsize = None : function, optional
        Accepts a Unit object argument and returns the size parameter. If None, defaults to function predefined for given name and units of measure.
    N : str
        Attribute name for number of parallel units.
        
    Examples
    --------
    :doc:`../../tutorial/Unit_decorators`
    
    """
    
    return lambda cls: add_cost(cls, ID, basis, units, S, ub, CE, cost, n, kW, BM, fsize, N)

def add_cost(cls, ID, basis, units, S, ub, CE, cost, n, kW, BM, fsize, N):
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
    else:
        ID = ID or cls.line
        cls.cost_items = {ID: CostItem(basis, units, S, ub, CE, cost, n, kW, N)}
        if '_BM' not in cls.__dict__:
            cls._BM = cls._BM.copy() if hasattr(cls, '_BM') else {}
        cls._BM[ID] = BM
        if '_cost' not in cls.__dict__:
            cls._cost = _decorated_cost
        cls._decorated_cost = _decorated_cost
    return cls

