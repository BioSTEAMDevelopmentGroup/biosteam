# -*- coding: utf-8 -*-
<<<<<<< HEAD
"""
Created on Thu Aug 23 14:26:41 2018

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
from .._graphics import mixer_graphics

__all__ = ('Mixer',)

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
    >>> settings.set_thermo(['Ethanol', 'Water'])
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
    
<<<<<<< HEAD
=======
    def _assert_compatible_property_package(self): 
        pass # Not necessary for mixing streams
    
>>>>>>> cd2c5013aaf9b5bc94bb764b52fd37db183472f1
    def _run(self):
        s_out, = self.outs
        s_out.mix_from(self.ins)
