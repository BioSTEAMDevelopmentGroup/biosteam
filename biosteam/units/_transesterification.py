# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from .. import Unit
from .decorators import cost
from thermosteam.reaction import Reaction, ParallelReaction

__all__ = ('Transesterification',)

@cost('Volume', 'Reactor',
      CE=525.4, cost=15000, n=0.55, kW=1.5, BM=4.3,)
class Transesterification(Unit):
    """
    Create a transesterification reactor that converts 'TAG', 'DAG', 'MAG' and 'Methanol'
    to 'Biodiesel' and 'Glycerol'. Finds the amount of catalyst 'NaOCH3'
    required and consumes it to produce 'NaOH' and 'Methanol'.
    
    Parameters
    ----------
    ins : stream sequence
        * [0] Lipid feed
        * [1] Methanol feed (includes catalyst)
    outs : stream
        Reactor effluent.
    efficiency : float
        Efficiency of conversion (on a 'Lipid' basis).
    excess_methanol : float
        Excess methanol fed in addition to the required stoichiometric amount.
    x_catalyst : float
        Catalyst molar fraction in methanol feed.
    T : float
        Operating temperature [K].
    
    """
    _bounds = {'Volume': (0.1, 20)}
    _units = {'Volume': 'm^3'}
    _tau = 1
    _N_ins = 2
    _N_outs = 1
    _N_heat_utilities = 1

    def _get_design_info(self):
        return (('Residence time', self.tau, 'hr'),
                ('Conversion efficiency', self.efficiency, ''),
                ('Working volume fraction', 0.8, ''))

    def __init__(self, ID='', ins=None, outs=(), *,
                 efficiency, excess_methanol, T, x_catalyst):
        Unit.__init__(self, ID, ins, outs)
        #: [:class:`~thermosteam.ParallelReaction`] Transesterification and catalyst consumption reactions.
        self.transesterification = ParallelReaction([
          Reaction('MAG + Methanol -> Biodiesel + Glycerol', reactant='MAG',  X=efficiency),
          Reaction('DAG + 2Methanol -> 2Biodiesel + Glycerol', reactant='DAG',  X=efficiency),
          Reaction('TAG + 3Methanol -> 3Biodiesel + Glycerol', reactant='TAG',  X=efficiency),
          Reaction('NaOCH3 -> NaOH + Methanol', reactant='NaOCH3', X=1)
        ])
        self.x_catalyst = x_catalyst #: [float] Catalyst molar fraction in methanol feed.
        self.excess_methanol = excess_methanol #: [float] Excess methanol fed in addition to the required stoichiometric amount.
        self.T = T #: [float] Operating temperature (K).
    
    @property
    def tau(self):
        """Residence time (hr)."""
        return self._tau
    @tau.setter
    def tau(self, tau):
        self._tau = tau
    
    @property
    def efficiency(self):
        """Transesterification efficiency."""
        return self.transesterification.X[:-1].mean()
    @efficiency.setter
    def efficiency(self, efficiency):
        self.transesterification.X[:-1] = efficiency
        
    def _run(self):
        feed, methanol = self.ins
        product, = self.outs
        x_catalyst = self.x_catalyst
        x_methanol = 1. - x_catalyst
        required_methanol = 1. + self.excess_methanol
        methanol.imol['Methanol'] = x_methanol * required_methanol * (
              feed.imol['MAG'] + 2 * feed.imol['DAG'] + 3 * feed.imol['TAG'] 
        )
        methanol.imol['NaOCH3'] =  x_catalyst / x_methanol * methanol.imol['Methanol'] 
        product.mix_from([feed, methanol], energy_balance=False)
        self.transesterification(product)
        product.T = self.T
        
    def _design(self):
        effluent = self._outs[0]
        self.design_results['Volume'] = self._tau * effluent.F_vol / 0.8
        self.heat_utilities[0](self.Hnet, effluent.T)

        
    
        
        
