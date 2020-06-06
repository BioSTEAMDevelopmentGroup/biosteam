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

__all__ = ('LLEUnit',)

class LLEUnit(bst.Unit, isabstract=True):
    r"""
    Abstract class for simulating liquid-liquid equilibrium.

    Parameters
    ----------
    ins : stream
        Inlet fluid.
    outs : stream sequence
        * [0] 'liquid' phase fluid
        * [1] 'LIQUID' phase fluid
    top_chemical : str, optional
        Identifier of chemical that will be favored in the "liquid" phase.
        If none given, the "liquid" phase will the lightest and the "LIQUID"
        phase will be the heaviest.
    efficiency : float,
        Fraction of feed in liquid-liquid equilibrium.
        The rest of the feed is divided equally between phases.
    
    Examples
    --------
    >>> from biorefineries.lipidcane.chemicals import lipidcane_chemicals
    >>> from biosteam import units, settings, Stream
    >>> settings.set_thermo(lipidcane_chemicals)
    >>> feed = Stream('feed', T=333.15,
    ...               Lipid=0.996, Biodiesel=26.9,
    ...               Methanol=32.9, Glycerol=8.97)
    >>> C1 = units.LLEUnit('C1', ins=feed, outs=('light', 'heavy'))
    >>> C1.simulate()
    >>> C1.show()
    LLEUnit: C1
    ins...
    [0] feed
        phase: 'l', T: 333.15 K, P: 101325 Pa
        flow (kmol/hr): Methanol   32.9
                        Glycerol   8.97
                        Biodiesel  26.9
                        Lipid      0.996
    outs...
    [0] light
        phase: 'l', T: 333.15 K, P: 101325 Pa
        flow (kmol/hr): Methanol   10.2
                        Glycerol   0.0239
                        Biodiesel  26.9
                        Lipid      0.996
    [1] heavy
        phase: 'l', T: 333.15 K, P: 101325 Pa
        flow (kmol/hr): Methanol   22.7
                        Glycerol   8.95
                        Biodiesel  0.0031
    
    """
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 top_chemical=None, efficiency=1.0):
        bst.Unit.__init__(self, ID, ins, outs, thermo)
        #: [str] Identifier of chemical that will be favored in the "liquid" phase.
        #: If none given, the "liquid" phase will the lightest and the "LIQUID"
        #: phase will be the heaviest.
        self.top_chemical = top_chemical
        #: Fraction of feed in liquid-liquid equilibrium.
        #: The rest of the feed is divided equally between phases.
        self.efficiency = efficiency 
        
    def _run(self):
        feed = self.ins[0]
        top, bottom = self.outs
        ms = feed.copy()
        ms.lle(feed.T, top_chemical=self.top_chemical)
        top_chemical = self.top_chemical
        top_phase = 'l'
        bottom_phase = 'L'
        if not top_chemical:
            rho_l = ms['l'].rho
            rho_L = ms['L'].rho
            top_L = rho_L < rho_l
            if top_L:
                top_phase = 'L'
                bottom_phase = 'l'
        top.mol[:] = ms.imol[top_phase]
        bottom.mol[:] = ms.imol[bottom_phase]
        top.T = bottom.T = feed.T
        top.P = bottom.P = feed.P
        efficiency = self.efficiency
        if efficiency < 1.:
            top.mol *= efficiency
            bottom.mol *= efficiency
            mixing = (1 - efficiency) / 2 * feed.mol
            top.mol += mixing
            bottom.mol += mixing