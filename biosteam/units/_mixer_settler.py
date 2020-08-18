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
from .._graphics import mixer_settler_graphics

__all__ = ('MixerSettler',)

class MixerSettler(bst.Unit):
    _N_ins = 2
    _ins_size_is_fixed = False
    _N_outs = 2
    auxiliary_unit_names = ('mixer', 'settler')
    _graphics = mixer_settler_graphics
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, 
                 mixer_data={}, settler_data={}, model="LLE"):
        bst.Unit.__init__(self, ID, ins, outs, thermo)
        self.mixer = mixer = bst.LiquidsMixingTank(None, None, (None,),
                                                   self.thermo, **mixer_data)
        self.mixed_stream = mixed_stream = mixer-0
        mixer._ins = self._ins
        model = model.lower()
        if model == 'lle':
            Settler = bst.LLESettler
        elif model == 'split':
            Settler = bst.LiquidsSplitSettler
        elif model == 'partition coefficients':
            raise NotImplementedError("partition coefficient model not yet implemented in BioSTEAM")
            Settler = bst.LiquidsPartitionSettler
        self.settler = Settler(None, mixed_stream, None, self.thermo, **settler_data)
        self.settler._outs = self._outs
        self.power_utility = mixer.power_utility
        
    def _run(self):
        self.mixer._run()
        self.settler._run()
        
    def _design(self):
        self.mixer._design()
        self.settler._design()
        
    def _cost(self):
        self.mixer._cost()
        self.settler._cost()
        
class MultiStageMixerSettlers:
    _N_ins = 2
    _N_outs = 2
    def __init__(self, ID="", ins=None, outs=(), thermo=None, *,
                 mixer_kwargs, N_stages):
        super().__init__(ID, ins, outs, thermo)
        self.N_stages = N_stages
        self.stages = None
        
    def _setup(self):
        if self.N_stages != self._N_stages_previous:
            self._N_stages_previous = self.N_stages
            self.reset_cache()
            
    def reset_cache(self):
        self._stages = []