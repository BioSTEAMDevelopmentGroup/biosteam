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
from .design_tools import separations

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
    >>> from biorefineries.lipidcane import chemicals
    >>> from biosteam import units, settings, Stream
    >>> settings.set_thermo(chemicals)
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
        separations.lle(self.ins[0], *self.outs, self.top_chemical, self.efficiency)