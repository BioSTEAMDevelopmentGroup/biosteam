# -*- coding: utf-8 -*-
<<<<<<< HEAD
"""
Created on Thu Aug 23 22:07:01 2018

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
from ._tank import MixTank
from ._hx import HXutility

__all__ = ('EnzymeTreatment',)

class EnzymeTreatment(MixTank):
    """Create an EnzymeTreatment unit that is cost as a MixTank with a heat exchanger."""
    _N_outs = 1
<<<<<<< HEAD
=======
    auxiliary_unit_names = ('heat_exchanger',)
>>>>>>> cd2c5013aaf9b5bc94bb764b52fd37db183472f1
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 vessel_type='Conventional', tau=1.0, V_wf=0.9,
                 vessel_material='Stainless steel', T):
        MixTank.__init__(self, ID, ins, outs, thermo, vessel_type=vessel_type,
                         tau=tau, V_wf=V_wf, vessel_material=vessel_material)
        self.T = T #: Operating temperature
<<<<<<< HEAD
        self.heat_exchanger = hx = HXutility(None, None, None, T=T) 
        self.heat_utilities = hx.heat_utilities
=======
        self.heat_exchanger = HXutility(None, None, None, T=T)
        self.heat_utilities = self.heat_exchanger.heat_utilities
>>>>>>> cd2c5013aaf9b5bc94bb764b52fd37db183472f1
    
    def _run(self):
        feed = self.ins[0]
        out = self.outs[0]
        out.mol[:] = self.mol_in
        out.phase = feed.phase
        out.P = feed.P
        out.T = self.T
        
    def _design(self):
        super()._design()
<<<<<<< HEAD
        self.heat_exchanger.simulate_as_auxiliary_exchanger(self.Hnet, self.outs[0])
        
    def _cost(self):
        super()._cost()
        self.purchase_costs['Heat exchanger'] = self.heat_exchanger.purchase_costs['Heat exchanger'] 
    
=======
        self.heat_exchanger.simulate_as_auxiliary_exchanger(self.Hnet, self.outs[0])
>>>>>>> cd2c5013aaf9b5bc94bb764b52fd37db183472f1
