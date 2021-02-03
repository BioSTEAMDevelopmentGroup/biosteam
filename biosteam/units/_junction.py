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

chemicals_in_common = lambda upstream, downstream: \
    tuple(set(upstream.chemicals.IDs).intersection(downstream.chemicals.IDs))

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
    def __init__(self, ID="", upstream=None, downstream=None, thermo=None):
        thermo = self._load_thermo(thermo)
        self._specification = None
        self._chemicals_in_common = self._past_streams = ()
        self._ins = Inlets(self, 1, upstream, thermo, True, self._stacklevel)
        self._outs = Outlets(self, 1, downstream, thermo, True, self._stacklevel)
        self._register(ID)
    
    def set_spec(self, *args, **kwargs):
        raise TypeError("{type(self).__name__}' does not support design specifications")
        
    def get_spec(self):
        return None
    
    def _get_chemicals_in_common(self, upstream, downstream):
        if (upstream, downstream) == self._past_streams:
            IDs = self._chemicals_in_common
        else:
            self._chemicals_in_common = IDs = chemicals_in_common(upstream, downstream)
        return IDs
    
    def _get_streams(self):
        try:
            upstream, = self._ins
            downstream, = self._outs
        except ValueError as error:
            N_ins = self._ins.size
            N_outs = self._outs.size
            if N_ins != 1:
                raise RuntimeError(f'a Junction object must have 1 input stream, not {N_ins}')
            elif N_outs != 1:
                raise RuntimeError(f'a Junction object must have 1 output stream, not {N_outs}')
            else:
                raise error
        return upstream, downstream
    
    def _run(self):
        upstream, downstream = self._get_streams()
        IDs = self._get_chemicals_in_common(upstream, downstream)
        if isinstance(upstream, MultiStream):
            downstream.phases = upstream.phases
            downstream.imol[..., IDs] = upstream.imol[..., IDs]
        else:
            downstream.phase = upstream.phase
            downstream.imol[IDs] = upstream.imol[IDs]
        downstream.T = upstream.T
        downstream.P = upstream.P
    simulate = _run
