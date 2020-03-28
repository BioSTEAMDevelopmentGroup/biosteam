# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 22:49:58 2018

@author: yoelr
"""
from .. import Unit
from .decorators import cost
from thermosteam.reaction import Reaction, ParallelReaction

__all__ = ('Transesterification',)

@cost('Volume', 'Reactor',
      CE=525.4, cost=15000, n=0.55, kW=1.5, BM=4.3,)
class Transesterification(Unit):
    """
    Create a transesterification reactor that converts 'Lipid' and 'Methanol'
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
    methanol2lipid : float
        Methanol feed to lipid molar ratio.
    catalyst_molfrac : float
        Catalyst to methanol molar ratio.
    T : float
        Operating temperature [K].
    
    """
    _bounds = {'Volume': (0.1, 20)}
    _units = {'Volume': 'm^3'}
    _tau = 1
    _N_ins = 2
    _N_outs = 1
    _N_heat_utilities = 1

    def _more_design_specs(self):
        return (('Residence time', self.tau, 'hr'),
                ('Conversion efficiency', self.efficiency, ''),
                ('Working volume fraction', 0.8, ''))

    def __init__(self, ID='', ins=None, outs=(), *,
                 efficiency, methanol2lipid, T, catalyst_molfrac):
        Unit.__init__(self, ID, ins, outs)
        #: [:class:`~thermosteam.ParallelReaction`] Transesterification and catalyst consumption reaction
        self.reaction = ParallelReaction([
          Reaction('Lipid + 3Methanol -> 3Biodiesel + Glycerol',
                   reactant='Lipid',  X=efficiency),
          Reaction('NaOCH3 -> NaOH + Methanol',
                   reactant='NaOCH3', X=1)])
        chemicals = self.chemicals
        self._methanol_composition = chemicals.kwarray(
                dict(Methanol=1-catalyst_molfrac,
                     NaOCH3=catalyst_molfrac))
        self._lipid_index, self._methanol_index, self._catalyst_index = \
                chemicals.get_index(('Lipid', 'Methanol', 'NaOCH3'))
        self._methanol2lipid = methanol2lipid
        self.T = T #: Operating temperature (K).
    
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
        return self.reaction.X[0]
    @efficiency.setter
    def efficiency(self, efficiency):
        self.reaction.X[0] = efficiency
    
    @property
    def methanol2lipid(self):
        """Methanol feed to lipid molar ratio."""
        return self._methanol2lipid
    @methanol2lipid.setter
    def methanol2lipid(self, ratio):
        self._methanol2lipid = ratio

    @property
    def catalyst_molfrac(self):
        """Catalyst molar fraction in methanol feed."""
        return self._methanol_composition[self._catalyst_index]
    @catalyst_molfrac.setter
    def catalyst_molfrac(self, molfrac):
        meoh = self._methanol_composition
        meoh[self._catalyst_index] = molfrac
        meoh[self._methanol_index] = 1-molfrac

    def _run(self):
        feed, methanol = self.ins
        product, = self.outs
        methanol.mol[:] = (self._methanol_composition
                         * feed.mol[self._lipid_index]
                         * self._methanol2lipid)
        product.mol[:] = feed.mol + methanol.mol
        self.reaction(product.mol)
        product.T = self.T
        
    def _design(self):
        effluent = self._outs[0]
        self.design_results['Volume'] = self._tau * effluent.F_vol / 0.8
        self.heat_utilities[0](self.Hnet, effluent.T)

        
    
        
        
