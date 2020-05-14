# -*- coding: utf-8 -*-
"""
Created on Thu May  7 14:21:38 2020

@author: yoelr
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