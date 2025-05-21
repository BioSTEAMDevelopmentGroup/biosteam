# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from .._unit import Unit
from thermosteam._graphics import junction_graphics
from .._power_utility import PowerUtility
from ..exceptions import UndefinedChemical
from ..utils.piping import Inlets, Outlets

__all__ = ('Junction',)


# %% Connect between different property packages

class Junction(Unit):
    """
    Create a Junction object that copies data from `upstream` to `downstream`. 
    This serves to connect streams with different property packages.
    
    Parameters
    ----------    
    upstream=None : Stream or str, defaults to missing stream
        Stream that will be copied to `downstream`.
    downstream="" : Stream or str, defaults to missing stream
        Flow rate, T, P, and phase information
        will be copied from `upstream` to this stream.
        If None, stream will be missing.
    
    Examples
    --------
    Create a Junction object and connect streams with different chemicals:
        
    >>> from biosteam import *
    >>> settings.set_thermo(['Water'])
    >>> s1 = Stream('s1', Water=20)
    >>> settings.set_thermo(['Ethanol', 'Water'], cache=True)
    >>> s2 = Stream('s2') # Note that s2 and s1 have different chemicals defined
    >>> J1 = units.Junction('J1', s1, s2)
    >>> J1.simulate()
    >>> J1.show()
    Junction: J1
    ins...
    [0] s1
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Water  20
    outs...
    [0] s2
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
        self._register(ID)
        thermo = self._load_thermo(thermo)
        self._init_specifications()
        self._isdynamic = False
        self.heat_utilities = []
        self._ins = Inlets(self, 1, upstream, thermo, True, self._stacklevel)
        self._outs = Outlets(self, 1, downstream, thermo, True, self._stacklevel)
    
    @property
    def upstream(self):
        return self._ins[0]
    @upstream.setter
    def upstream(self, upstream):
        self._ins[0] = upstream
        
    @property
    def downstream(self):
        return self._outs[0]
    @downstream.setter
    def downstream(self, downstream):
        self._outs[0] = downstream
        
    def _setup(self): 
        for ps in self._specifications: ps.compile(self)
    
    def _run(self): 
        upstream = self._ins[0]
        downstream = self._outs[0]
        try: downstream.copy_like(upstream)
        except UndefinedChemical:
            self._reset_thermo(self._ins[0]._thermo)
            downstream.copy_like(upstream)
    simulate = Unit.run

    def _get_tooltip_string(self):
        return f"{type(self).__name__}: {self.ID}"

    @property
    def _inlet_utility_indices(self): return {}
    
    @property
    def _outlet_utility_indices(self): return {}

