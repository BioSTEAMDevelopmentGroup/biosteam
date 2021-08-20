# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from .tank import MixTank
from .heat_exchange import HXutility

__all__ = ('EnzymeTreatment',)

class EnzymeTreatment(MixTank):
    """Create an EnzymeTreatment unit that is cost as a MixTank with a heat exchanger."""
    _N_outs = 1
    auxiliary_unit_names = ('heat_exchanger',)
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 vessel_type='Conventional', tau=1.0, V_wf=0.9,
                 vessel_material='Stainless steel', T):
        MixTank.__init__(self, ID, ins, outs, thermo, vessel_type=vessel_type,
                         tau=tau, V_wf=V_wf, vessel_material=vessel_material)
        self.T = T #: Operating temperature
        self.heat_exchanger = HXutility(None, None, None, T=T, thermo=thermo)
        self.heat_utilities = self.heat_exchanger.heat_utilities
    
    def _run(self):
        feed = self.ins[0]
        out = self.outs[0]
        out.copy_like(feed)
        out.T = self.T
        
    def _design(self):
        super()._design()
        self.heat_exchanger.simulate_as_auxiliary_exchanger(self.Hnet, self.outs[0])