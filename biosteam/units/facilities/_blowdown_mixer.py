# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from . import Facility
from ..._graphics import mixer_graphics

__all__ = ('BlowdownMixer',)

class BlowdownMixer(Facility):
    network_priority = 2
    _graphics = mixer_graphics
    _N_outs = 1
    _N_ins = 2
    _ins_size_is_fixed = False
    
    def _run(self):
        s_out, = self.outs
        s_out.mix_from(self.ins)