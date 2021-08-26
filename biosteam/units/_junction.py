# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from .._unit import Unit
from .._graphics import junction_graphics
from .._power_utility import PowerUtility
from thermosteam import MultiStream
from ..utils.piping import Inlets, Outlets

__all__ = ('Junction',)


# %% Connect between different property packages

class Junction(Unit):
    """
    Create a Junction object that copies specifications from `upstream`
    to `downstream`. This serves to connect streams with different
    Species object.
    
    Parameters
    ----------    
    upstream=None : Stream or str, defaults to missing stream
        Stream that will be copied to `downstream`.
    downstream="" : Stream or str, defaults to new stream
        Flow rate, T, P, and phase information
        will be copied from `upstream` to this stream.
        If None, stream will be missing.
    species=None : list[str], defaults to all species in common
        IDs of species to be passed down.
    
    Examples
    --------
    Create a Junction object and connect streams with different Species objects:
        
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
    heat_utilities = ()
    power_utility = PowerUtility()
    baseline_purchase_cost = 0.
    baseline_purchase_costs = {}
    purchase_cost = 0.
    purchase_costs = {}
    installed_cost = 0.
    installed_costs = {}
    utility_cost = 0.
    
    def __init__(self, ID="", upstream=None, downstream=None, thermo=None):
        self._register(ID)
        thermo = self._load_thermo(thermo)
        self._specification = None
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
        pass
    
    def _run(self): 
        self._outs[0].copy_like(self._ins[0])
    simulate = _run
