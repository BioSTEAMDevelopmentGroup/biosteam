# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
This module contains unit operations for mixing.

.. contents:: :local:
    
Unit operations
---------------
.. autoclass:: biosteam.units.mixing.Mixer
.. autoclass:: biosteam.units.mixing.SteamMixer

"""
from ..utils import InletPort, OutletPort, ignore_docking_warnings
from .._unit import Unit
from .._graphics import mixer_graphics
import flexsolve as flx
import biosteam as bst

__all__ = ('Mixer', 'SteamMixer')

class Mixer(Unit):
    """
    Create a mixer that mixes any number of streams together.
    
    Parameters
    ----------
    ins : streams
        Inlet fluids to be mixed.
    outs : stream
        Mixed outlet fluid.
    Examples
    --------
    Mix two streams:
    
    >>> from biosteam import units, settings, Stream
    >>> settings.set_thermo(['Ethanol', 'Water'], cache=True)
    >>> s1 = Stream('s1', Water=20, T=350)
    >>> s2 = Stream('s2', Ethanol=30, T=300)
    >>> M1 = units.Mixer('M1', ins=(s1, s2), outs='s3')
    >>> M1.simulate()
    >>> M1.show()
    Mixer: M1
    ins...
    [0] s1
        phase: 'l', T: 350 K, P: 101325 Pa
        flow (kmol/hr): Water  20
    [1] s2
        phase: 'l', T: 300 K, P: 101325 Pa
        flow (kmol/hr): Ethanol  30
    outs...
    [0] s3
        phase: 'l', T: 315.11 K, P: 101325 Pa
        flow (kmol/hr): Ethanol  30
                        Water    20
    
    
    """
    _graphics = mixer_graphics
    _N_outs = 1
    _N_ins = 2
    _ins_size_is_fixed = False
    
    def _assert_compatible_property_package(self): 
        pass # Not necessary for mixing streams
    
    def _run(self):
        s_out, = self.outs
        s_out.mix_from(self.ins)
        
    @ignore_docking_warnings
    def insert(self, stream):
        """
        Insert Mixer object between two units at a given stream connection.
        
        Examples
        --------
        >>> from biosteam import *
        >>> settings.set_thermo(['Water'], cache=True)
        >>> feed = Stream('feed')
        >>> other_feed = Stream('other_feed')
        >>> P1 = Pump('P1', feed, 'pump_outlet')
        >>> H1 = HXutility('H1', P1-0, T=310)
        >>> M1 = Mixer('M1', other_feed, 'mixer_outlet')
        >>> M1.insert(P1-0)
        >>> M1.show()
        Mixer: M1
        ins...
        [0] other_feed
            phase: 'l', T: 298.15 K, P: 101325 Pa
            flow: 0
        [1] pump_outlet  from  Pump-P1
            phase: 'l', T: 298.15 K, P: 101325 Pa
            flow: 0
        outs...
        [0] mixer_outlet  to  HXutility-H1
            phase: 'l', T: 298.15 K, P: 101325 Pa
            flow: 0
        
        """
        sink = stream.sink
        sink.ins.replace(stream, self.outs[0])
        self.ins.append(stream)
        

class SteamMixer(Unit):
    """
    Create a mixer that varies the flow of steam to achieve a specified outlet
    pressure.
    
    Parameters
    ----------
    ins : stream sequence
        [0] Feed    
        [1] Steam
    outs : stream
        Mixed product.
    P : float
        Outlet pressure.
    
    """
    _N_outs = 1
    _N_ins = 2
    _ins_size_is_fixed = False
    _N_heat_utilities = 1
    _graphics = mixer_graphics
    installation_cost = purchase_cost = 0.
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *, P):
        Unit.__init__(self, ID, ins, outs, thermo)
        self.P = P
    
    @property
    def steam(self):
        return self.ins[1]
    
    def reset_cache(self): 
        for utility in bst.HeatUtility.heating_agents:
            if utility.P > self.P: break
        self.steam.copy_like(utility)
    
    def pressure_objective_function(self, steam_mol):
        steam = self.ins[1]
        mixed = self.outs[0]
        steam.imol['7732-18-5'] = steam_mol
        mixed.mix_from(self.ins)
        P_new = mixed.chemicals.Water.Psat(min(mixed.T, mixed.chemicals.Water.Tc - 1))
        return self.P - P_new
    
    def _setup(self):
        super()._setup()
        if self.steam.isempty(): self.reset_cache()
        self.outs[0].P = self.P
    
    def _run(self):
        steam = self.ins[1]
        steam_mol = steam.F_mol
        f = self.pressure_objective_function
        steam_mol = flx.IQ_interpolation(f, *flx.find_bracket(f, 0., steam_mol, None, None), 
                                         xtol=1e-2, ytol=1e-4, 
                                         maxiter=500, checkroot=False)
        
    def _design(self): 
        steam = self.ins[1]
        mixed = self.outs[0]
        self.heat_utilities[0](steam.H, mixed.T)