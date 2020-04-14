# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 22:07:01 2018

@author: yoelr
"""
from ._tank import MixTank
from ._hx import HXutility

__all__ = ('EnzymeTreatment',)

class EnzymeTreatment(MixTank):
    """Create an EnzymeTreatment unit that is cost as a MixTank with a heat exchanger."""
    _N_outs = 1
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 vessel_type='Conventional', tau=1.0, V_wf=0.9,
                 vessel_material='Stainless steel', T):
        MixTank.__init__(self, ID, ins, outs, thermo, vessel_type=vessel_type,
                         tau=tau, V_wf=V_wf, vessel_material=vessel_material)
        self.T = T #: Operating temperature
        self.heat_exchanger = hx = HXutility(None, None, None, T=T) 
        self.heat_utilities = hx.heat_utilities
    
    def _run(self):
        feed = self.ins[0]
        out = self.outs[0]
        out.mol[:] = self.mol_in
        out.phase = feed.phase
        out.P = feed.P
        out.T = self.T
        
    def _design(self):
        super()._design()
        self.heat_exchanger.simulate_as_auxiliary_exchanger(self.Hnet, self.outs[0])
        
    def _cost(self):
        super()._cost()
        self.purchase_costs['Heat exchanger'] = self.heat_exchanger.purchase_costs['Heat exchanger'] 
    