# -*- coding: utf-8 -*-
<<<<<<< HEAD
"""
Created on Thu Aug 23 21:23:56 2018

@author: yoelr
=======
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
>>>>>>> cd2c5013aaf9b5bc94bb764b52fd37db183472f1
"""
from .. import Unit
from ._flash import RatioFlash
from ._splitter import Splitter
<<<<<<< HEAD
from .decorators import cost
import thermosteam as tmo
=======
from ._lle_unit import LLEUnit
from .decorators import cost
>>>>>>> cd2c5013aaf9b5bc94bb764b52fd37db183472f1

__all__ = ('LiquidsCentrifuge',
           'LiquidsSplitCentrifuge', 
           'LiquidsRatioCentrifuge',
           'LLECentrifuge')

# electricity kW/(m3/hr) from USDA biosdiesel Super Pro model
# Possibly 1.4  kW/(m3/hr)
# https://www.sciencedirect.com/topics/engineering/disc-stack-centrifuge
# Microalgal fatty acids—From harvesting until extraction H.M. Amaro, et. al.,
# in Microalgae-Based Biofuels and Bioproducts, 2017

@cost('Flow rate', units='m^3/hr', CE=525.4, cost=28100,
<<<<<<< HEAD
      n=0.574, kW=3.66, ub=100, BM=2.03, N='Number of centrifuges')
=======
      n=0.574, kW=1.4, ub=100, BM=2.03, N='Number of centrifuges')
>>>>>>> cd2c5013aaf9b5bc94bb764b52fd37db183472f1
class LiquidsCentrifuge(Unit, isabstract=True):
    r"""
    Abstract class for liquid centrifuges.

    Parameters
    ----------
    ins : stream
        Inlet fluid.
    outs : stream sequence
        * [0] 'liquid' phase fluid
        * [1] 'LIQUID' phase fluid
    
    Notes
    -----
    The f.o.b purchase cost is given by [1]_:

    :math:`C_{f.o.b}^{2007} = 28100 Q^{0.574} (Q < 100 \frac{m^3}{h})` 
    
    References
    ----------
    .. [1] Apostolakou, A. A.; Kookos, I. K.; Marazioti, C.; Angelopoulos, 
        K. C. Techno-Economic Analysis of a Biodiesel Production Process 
        from Vegetable Oils. Fuel Process. Technol. 2009, 90, 1023−1031
    
    """
    _N_outs = 2
    line = 'Liquids centrifuge'


class LiquidsRatioCentrifuge(LiquidsCentrifuge):
    _N_heat_utilities = 0
    line = 'Liquids centrifuge'
    __init__ = RatioFlash.__init__
    _run = RatioFlash._run


class LiquidsSplitCentrifuge(LiquidsCentrifuge):
    r"""
    Create a liquids centrifuge simulated by component splits.

    Parameters
    ----------
    ins : stream
        Inlet fluid.
    outs : stream sequence
        * [0] 'liquid' phase fluid
        * [1] 'LIQUID' phase fluid
    split : Should be one of the following
            * [float] The fraction of net feed in the 0th outlet stream
            * [array_like] Componentwise split of feed to 0th outlet stream
            * [dict] ID-split pairs of feed to 0th outlet stream
    order=None : Iterable[str], defaults to biosteam.settings.chemicals.IDs
        Chemical order of split.
    
    Notes
    -----
    The f.o.b purchase cost is given by [1]_:

    :math:`C_{f.o.b}^{2007} = 28100 Q^{0.574} (Q < 100 \frac{m^3}{h})` 
    
    References
    ----------
    .. [1] Apostolakou, A. A.; Kookos, I. K.; Marazioti, C.; Angelopoulos, 
        K. C. Techno-Economic Analysis of a Biodiesel Production Process 
        from Vegetable Oils. Fuel Process. Technol. 2009, 90, 1023−1031
    
    """
    line = 'Liquids centrifuge'
    __init__ = Splitter.__init__
    _run = Splitter._run
    split = Splitter.split
    
    
<<<<<<< HEAD
class LLECentrifuge(LiquidsCentrifuge):
=======
class LLECentrifuge(LLEUnit, LiquidsCentrifuge):
>>>>>>> cd2c5013aaf9b5bc94bb764b52fd37db183472f1
    r"""
    Create a liquids centrifuge simulated by liquid-liquid equilibrium.

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
<<<<<<< HEAD
=======
    efficiency : float,
        Fraction of feed in liquid-liquid equilibrium.
        The rest of the feed is divided equally between phases.
>>>>>>> cd2c5013aaf9b5bc94bb764b52fd37db183472f1
    
    Notes
    -----
    The f.o.b purchase cost is given by [1]_:

    :math:`C_{f.o.b}^{2007} = 28100 Q^{0.574} (Q < 100 \frac{m^3}{h})` 
    
    References
    ----------
    .. [1] Apostolakou, A. A.; Kookos, I. K.; Marazioti, C.; Angelopoulos, 
        K. C. Techno-Economic Analysis of a Biodiesel Production Process 
        from Vegetable Oils. Fuel Process. Technol. 2009, 90, 1023−1031
    
    Examples
    --------
<<<<<<< HEAD
    >>> from biorefineries.lipidcane.chemicals import lipidcane_chemicals
    >>> from biosteam import units, settings, Stream
    >>> settings.set_thermo(lipidcane_chemicals)
=======
    >>> from biorefineries.lipidcane import chemicals
    >>> from biosteam import units, settings, Stream
    >>> settings.set_thermo(chemicals)
>>>>>>> cd2c5013aaf9b5bc94bb764b52fd37db183472f1
    >>> feed = Stream('feed', T=333.15,
    ...               Lipid=0.996, Biodiesel=26.9,
    ...               Methanol=32.9, Glycerol=8.97)
    >>> C1 = units.LLECentrifuge('C1', ins=feed, outs=('light', 'heavy'))
    >>> C1.simulate()
    >>> C1.show()
    LLECentrifuge: C1
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
    >>> C1.results()
    Liquids centrifuge                          Units       C1
<<<<<<< HEAD
    Power               Rate                       kW     44.7
                        Cost                   USD/hr      3.5
=======
    Power               Rate                       kW     17.1
                        Cost                   USD/hr     1.34
>>>>>>> cd2c5013aaf9b5bc94bb764b52fd37db183472f1
    Design              Flow rate              m^3/hr     12.2
                        Number of centrifuges                1
    Purchase cost       Liquids centrifuge        USD 1.28e+05
    Total purchase cost                           USD 1.28e+05
<<<<<<< HEAD
    Utility cost                               USD/hr      3.5
    
    """
    line = 'Liquids centrifuge'
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, top_chemical=None):
        super().__init__(ID, ins, outs, thermo)
        self._multistream = tmo.MultiStream(None, phases=('L', 'l'),
                                            thermo=thermo)
        self.top_chemical = top_chemical
        
    def _run(self):
        ms = self._multistream
        feed = self.ins[0]
        top, bottom = self.outs
        ms.imol['l'] = feed.mol
        ms.lle(feed.T)
        top_chemical = self.top_chemical
        if top_chemical:
            F_l = ms.imol['l', top_chemical]
            F_L = ms.imol['L', top_chemical]
            top_l = F_l > F_L
        else:
            rho_l = ms['l'].rho
            rho_L = ms['L'].rho
            top_l = rho_l < rho_L
        if top_l:
            top_phase = 'l'
            bottom_phase = 'L'
        else:
            top_phase = 'L'
            bottom_phase = 'l'
        top.mol[:] = ms.imol[top_phase]
        bottom.mol[:] = ms.imol[bottom_phase]
        top.T = bottom.T = feed.T
        top.P = bottom.P = feed.P
=======
    Utility cost                               USD/hr     1.34
    
    """
    line = 'Liquids centrifuge'
>>>>>>> cd2c5013aaf9b5bc94bb764b52fd37db183472f1
