# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import biosteam as bst
from .._unit import Unit
from .._graphics import junction_graphics
from .._power_utility import PowerUtility
from ..utils.piping import Inlets, Outlets

__all__ = ('Terminal',)


# %% Connect between different property packages

class Terminal(Unit):
    """
    Create a Terminal object that serves to denote a stream that is not a
    feed/product.
    
    Parameters
    ----------    
    upstream=None : Stream or str, defaults to missing stream
        Inlet stream.
    downstream="" : Stream or str, defaults to missing stream
        Outlet stream.
    
    Examples
    --------
    Create a Terminal object:
        
    >>> from biosteam import *
    >>> settings.set_thermo(['Water'])
    >>> s1 = Stream('s1', Water=20)
    >>> S1 = Terminal(upstream=s1)
    >>> S1.simulate()
    >>> S1.show()
    Terminal: S1
    ins...
    [0] s1  
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Water  20
        
    """
    _stacklevel = Unit._stacklevel
    _graphics = junction_graphics
    power_utility = PowerUtility()
    design_results = {}
    baseline_purchase_cost = 0.
    baseline_purchase_costs = {}
    purchase_cost = 0.
    purchase_costs = {}
    installed_cost = 0.
    installed_costs = {}
    utility_cost = 0.
    prioritize = False
    
    def __init__(self, ID="", upstream=None, downstream=None, thermo=None):
        if bst.settings.ID_magic and ID == '': ID = bst.utils.infer_variable_assignment(self.__class__)
        self._register(ID)
        thermo = self._load_thermo(thermo)
        self._init_specifications()
        self._isdynamic = False
        self.heat_utilities = []
        if upstream is None:
            N_in = 0
            N_out = 1
            if downstream is None:
                raise ValueError('must pass either upstream or downstream')
        elif downstream is None:
            N_in = 1
            N_out = 0
        else:
            raise ValueError('cannot pass both upstream and downstream')
        self._N_ins = N_in
        self._N_outs = N_out
        self._ins = Inlets(self, N_in, upstream, thermo, False, self._stacklevel)
        self._outs = Outlets(self, N_out, downstream, thermo, False, self._stacklevel)
    
    @property
    def upstream(self):
        if self._N_ins == 0: return None
        return self._ins[0]
    @upstream.setter
    def upstream(self, upstream):
        if self._N_ins == 0: raise AttributeError('Sink-terminal cannot have an upstream')
        self._ins[0] = upstream
        
    @property
    def downstream(self):
        if self._N_outs == 0: return None
        return self._outs[0]
    @downstream.setter
    def downstream(self, downstream):
        if self._N_outs == 0: raise AttributeError('Source-terminal cannot have a downstream')
        self._outs[0] = downstream
        
    def _setup(self): 
        for ps in self._specifications: ps.compile(self)
    
    def _run(self): 
        pass
    
    simulate = Unit.run

    def _get_tooltip_string(self):
        return f"{type(self).__name__}: {self.ID}"

    @property
    def _inlet_utility_indices(self): return {}
    
    @property
    def _outlet_utility_indices(self): return {}

