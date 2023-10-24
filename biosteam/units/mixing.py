# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
.. contents:: :local:

.. autoclass:: biosteam.units.mixing.Mixer
.. autoclass:: biosteam.units.mixing.SteamMixer
.. autoclass:: biosteam.units.mixing.MockMixer

"""
from .._unit import Unit
from .._graphics import mixer_graphics
import flexsolve as flx
import biosteam as bst
from typing import Optional

__all__ = ('Mixer', 'SteamMixer', 'FakeMixer', 'MockMixer')

class Mixer(Unit):
    """
    Create a mixer that mixes any number of streams together.
    
    Parameters
    ----------
    ins : 
        Inlet fluids to be mixed.
    outs : 
        Mixed outlet fluid.
    rigorous :
        Whether to perform vapor-liquid equilibrium.
    
    Notes
    -----
    When streams at different pressures are mixed, BioSTEAM assumes valves 
    reduce the pressure of the streams being mixed to prevent backflow 
    (pressure needs to decrease in the direction of flow according to 
    Bernoulli's principle). The outlet pressure will be the minimum pressure
    of all inlet streams.
    
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
        phase: 'l', T: 315.14 K, P: 101325 Pa
        flow (kmol/hr): Ethanol  30
                        Water    20
    
    
    """
    _graphics = mixer_graphics
    _N_outs = 1
    _N_ins = 2
    _ins_size_is_fixed = False
    
    def _assert_compatible_property_package(self): 
        pass # Not necessary for mixing streams
    
    def _init(self, rigorous: Optional[bool]=False,
              conserve_phases: Optional[bool]=False):
        self.rigorous = rigorous
        self.conserve_phases = conserve_phases
    
    def _run(self):
        s_out, = self.outs
        s_out.mix_from(self.ins, vle=self.rigorous,
                       conserve_phases=self.conserve_phases)
        

class SteamMixer(Unit):
    """
    Create a mixer that varies the flow of steam to achieve a specified outlet
    pressure and varies the flow of process water to achieve a specified 
    solids loading (by wt).
    
    Parameters
    ----------
    ins : 
        * [0] Feed    
        * [1] Steam
        * [2] Process water
    outs : 
        Mixed product.
    P : float
        Outlet pressure.
    T : float
        Outlet temperature.
    solids_loading : float, optional
        Final solids loading after mixing in process water.
    soilds_loading_includes_steam : bool, optional
        Whether to include steam in solids loading calculation.
    
    Examples
    --------
    >>> import biosteam as bst
    >>> bst.settings.set_thermo(['Water', 'Glucose'])
    >>> feed = bst.Stream('feed', Water=10, Glucose=10)
    >>> M1 = bst.SteamMixer(None, ins=[feed, 'steam', 'process_water'], outs='outlet', T=431.15, P=557287.5, solids_loading=0.3)
    >>> M1.simulate()
    >>> M1.show('cwt100') # Note that outlet solids loading is not exactly 0.3 because of the steam requirement.
    SteamMixer
    ins...
    [0] feed
        phase: 'l', T: 298.15 K, P: 101325 Pa
        composition (%): Water    9.09
                         Glucose  90.9
                         -------  1.98e+03 kg/hr
    [1] steam
        phase: 'g', T: 454.77 K, P: 1.041e+06 Pa
        composition (%): Water  100
                         -----  1.28e+03 kg/hr
    [2] process_water
        phase: 'l', T: 298.15 K, P: 101325 Pa
        composition (%): Water  100
                         -----  4.02e+03 kg/hr
    outs...
    [0] outlet
        phase: 'l', T: 431.15 K, P: 557288 Pa
        composition (%): Water    75.3
                         Glucose  24.7
                         -------  7.29e+03 kg/hr
    
    >>> M1.solids_loading_includes_steam = True
    >>> M1.simulate()
    >>> M1.show('cwt100') # Now the outlet solids content is exactly 0.3
    SteamMixer
    ins...
    [0] feed
        phase: 'l', T: 298.15 K, P: 101325 Pa
        composition (%): Water    9.09
                         Glucose  90.9
                         -------  1.98e+03 kg/hr
    [1] steam
        phase: 'g', T: 454.77 K, P: 1.041e+06 Pa
        composition (%): Water  100
                         -----  1.02e+03 kg/hr
    [2] process_water
        phase: 'l', T: 298.15 K, P: 101325 Pa
        composition (%): Water  100
                         -----  3.01e+03 kg/hr
    outs...
    [0] outlet
        phase: 'l', T: 431.15 K, P: 557288 Pa
        composition (%): Water    70
                         Glucose  30
                         -------  6.01e+03 kg/hr
    
    """
    _N_outs = 1
    _N_ins = 3
    _ins_size_is_fixed = False
    _graphics = mixer_graphics
    installation_cost = purchase_cost = 0.
    def _init(self, P, T=None, solids_loading=None, 
             liquid_IDs=['7732-18-5'], solid_IDs=None,
             solids_loading_includes_steam=None):
        self.P = P
        self.T = T
        self.solids_loading = solids_loading
        self.solids_loading_includes_steam = solids_loading_includes_steam
        self.liquid_IDs = tuple(liquid_IDs)
        self.solid_IDs = solid_IDs
    
    @property
    def steam(self):
        return self.ins[1]
    
    def reset_cache(self, isdynamic=None): 
        for utility in bst.HeatUtility.heating_agents:
            if utility.P > self.P: break
        self.steam.copy_like(utility)
    
    def pressure_objective_function(self, steam_mol):
        mixed = self.outs[0]
        self.steam.imol['7732-18-5'] = steam_mol # Only change water
        if self.P: mixed.P = self.P # Assume pumps take care of this
        if self.solids_loading_includes_steam and self.solids_loading:
            solids_loading = self.solids_loading
            feed, steam, process_water, *others = self.ins
            process_water.empty()
            feeds = self.ins
            chemicals = self.chemicals
            index = chemicals.get_index(self.liquid_IDs)
            available_water = 18.01528 * sum([(j.sum() if hasattr((j:=i.mol[index]), 'sum') else j) for i in feeds if i])
            solid_IDs = self.solid_IDs
            if solid_IDs:
                F_mass_solids = sum([i.imass[solid_IDs].sum() for i in feeds if i])
            else:
                F_mass_feed = sum([i.F_mass for i in feeds if i])
                F_mass_solids = F_mass_feed - available_water
            required_water = F_mass_solids * (1. - solids_loading) / solids_loading
            process_water.imol['7732-18-5'] = max(required_water - available_water, 0.) / 18.01528
        
        mixed.mix_from(self.ins, energy_balance=False)
        H = sum([i.H for i in self.ins])
        Tmax = mixed.chemicals.Water.Tc - 1
        mixed.T = Tmax
        Hmax = mixed.H
        if H > Hmax:
            mixed.T = Tmax + (H - Hmax) / mixed.chemicals.Water.Cn('l', Tmax)
        else:
            mixed.H = H
        if self.T:
            return self.T - mixed.T
        else: # If no pressure, assume it is at the boiling point
            P_new = mixed.chemicals.Water.Psat(min(mixed.T, Tmax))
            return self.P - P_new
    
    def _setup(self):
        super()._setup()
        if self.steam.isempty(): self.reset_cache()
    
    def _run(self):
        solids_loading = self.solids_loading
        if solids_loading is not None and not self.solids_loading_includes_steam:
            # Solids loading need to be achieved first before mixing with steam
            # to avoid pumping issues (see Humbird 2011 NREL report).
            feed, steam, process_water, *others = self.ins
            process_water.empty()
            feeds = [feed, *others]
            chemicals = self.chemicals
            index = chemicals.get_index(self.liquid_IDs)
            available_water = 18.01528 * sum([(j.sum() if hasattr((j:=i.mol[index]), 'sum') else j) for i in feeds if i])
            solid_IDs = self.solid_IDs
            if solid_IDs:
                F_mass_solids = sum([i.imass[solid_IDs].sum() for i in feeds if i])
            else:
                F_mass_feed = sum([i.F_mass for i in feeds if i])
                F_mass_solids = F_mass_feed - available_water
            required_water = F_mass_solids * (1. - solids_loading) / solids_loading
            process_water.imol['7732-18-5'] = max(required_water - available_water, 0.) / 18.01528
        else:
            steam = self.steam
        steam_mol = steam.F_mol or 1.
        f = self.pressure_objective_function
        steam_mol = flx.IQ_interpolation(f, *flx.find_bracket(f, 0., steam_mol, None, None), 
                                         xtol=1e-2, ytol=1e-4, 
                                         maxiter=500, checkroot=False)
        self.outs[0].P = self.P
        
    def _design(self): 
        steam = self.ins[1]
        warm_process_water = self.ins[2]
        mixed = self.outs[0]
        self.add_heat_utility(steam.H + warm_process_water.H, mixed.T)
        
class MockMixer(Unit):
    """
    Create a MockMixer object that does nothing when simulated.
    """
    _graphics = Mixer._graphics
    _N_ins = 2
    _N_outs = 1
    _ins_size_is_fixed = False
    
    def _run(self): pass

MockMixer.line = 'Mixer'
FakeMixer = MockMixer    
