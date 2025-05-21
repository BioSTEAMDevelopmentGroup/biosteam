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
from thermosteam._graphics import scaler_graphics
from .._power_utility import PowerUtility
from ..utils.piping import Inlets, Outlets

__all__ = ('Scaler',)


# %% Connect between different property packages

class Scaler(Unit):
    """
    Create a Scale object scales the flow rate of the outlet stream.
    
    Parameters
    ----------    
    upstream=None : Stream or str, defaults to missing stream
        Stream that will be copied to `downstream`.
    downstream="" : Stream or str, defaults to missing stream
        Flow rate, T, P, and phase information
        will be copied from `upstream` to this stream.
        If None, stream will be missing.
    scale=None : int, optional
        Downstream flow rate will be the upstream flow rate multiplied by the scale.
    
    Examples
    --------
    Scale outlet by 2:
        
    >>> from biosteam import *
    >>> settings.set_thermo(['Water'])
    >>> s1 = Stream('s1', Water=20)
    >>> settings.set_thermo(['Ethanol', 'Water'], cache=True)
    >>> s2 = Stream('s2') # Note that s2 and s1 have different chemicals defined
    >>> S1 = units.Scaler('S1', s1, s2, scale=2)
    >>> S1.simulate()
    >>> S1.show()
    Scaler: S1
    ins...
    [0] s1
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Water  20
    outs...
    [0] s2
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Water  40
        
    """
    _stacklevel = Unit._stacklevel
    _graphics = scaler_graphics
    power_utility = PowerUtility()
    design_results = {}
    baseline_purchase_cost = 0.
    baseline_purchase_costs = {}
    purchase_cost = 0.
    purchase_costs = {}
    installed_cost = 0.
    installed_costs = {}
    utility_cost = 0.
    
    def __init__(self, ID="", upstream=None, downstream=None, thermo=None, scale=1):
        self._register(ID)
        thermo = self._load_thermo(thermo)
        self._init_specifications()
        self._isdynamic = False
        self.heat_utilities = []
        self._ins = Inlets(self, 1, upstream, thermo, True, self._stacklevel)
        self._outs = Outlets(self, 1, downstream, thermo, True, self._stacklevel)
        self.scale = scale
    
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
        downstream.copy_like(upstream)
        downstream.scale(self.scale)
    simulate = Unit.run

    def _get_tooltip_string(self):
        return f"{type(self).__name__}"

    @property
    def _inlet_utility_indices(self): return {}
    
    @property
    def _outlet_utility_indices(self): return {}

