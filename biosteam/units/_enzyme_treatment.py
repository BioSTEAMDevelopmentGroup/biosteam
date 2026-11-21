# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from .tank import MixTank
from .heat_exchange import HXutility
from .._unit import design_input, update_design_input

__all__ = ('EnzymeTreatment',)

class EnzymeTreatment(MixTank):
    """Create an EnzymeTreatment unit that is cost as a MixTank with a heat exchanger."""
    _N_outs = 1
    auxiliary_unit_names = ('heat_exchanger',)
    T = design_input('Operating temperature', 'K', lb=0, default=None)
    V_wf = update_design_input(default=0.9)

    def _init(self):
        self.heat_exchanger = HXutility(None, None, None, T=self.T, thermo=self.thermo)
    
    def _run(self):
        out = self.outs[0]
        if self.T is None:
            out.mix_from(self.ins)
        else:
            out.mix_from(self.ins, energy_balance=False)
            out.T = self.T
        
    def _design(self):
        super()._design()
        self.heat_exchanger.simulate_as_auxiliary_exchanger(self.ins, self.outs, self.Hnet)
        self.parallel['heat_exchanger'] = 1