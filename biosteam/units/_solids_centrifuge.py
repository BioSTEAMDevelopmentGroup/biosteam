# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import numpy as np
from ..utils.unit_warnings import lb_warning
from .decorators import cost
from ._splitter import Splitter

__all__ = ('SolidsCentrifuge',)

@cost('Solids loading', cost=68040, CE=567, n=0.50, ub=40, BM=2.03,
      N='Number of centrifuges')
class SolidsCentrifuge(Splitter):
    """
    Create a solids centrifuge that separates out solids according to
    user defined split. Assume a continuous scroll solid bowl. 
    
    Parameters
    ----------
    ins : stream
        Inlet fluid with solids.
    outs : stream sequence
        * [0] Liquid-rich stream.
        * [1] Solids-rich stream.
    split: array_like
           Component splits to 0th output stream
    order=None : Iterable[str]
        Species order of split. Defaults to Stream.chemicals.IDs.
    solids : tuple[str]
             IDs of solids.
    
    References
    ----------
    .. [0] Seider, Warren D., et al. (2017). "Cost Accounting and Capital Cost Estimation". In Product and Process Design Principles: Synthesis, Analysis, and Evaluation (pp. 481-485). New York: Wiley.
    
    """
    _units = {'Solids loading': 'ton/hr',
              'Flow rate': 'm3/hr'}
    minimum_solids_loading = 2
    kWhr_per_m3 = 1.40


    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 split, order=None, solids=None):
        super().__init__(ID, ins, outs, thermo, split=split, order=order)
        self.solids = solids
    
    @property
    def solids(self):
        return self._solids
    @solids.setter
    def solids(self, solids):
        self._solids = tuple(solids)
    
    def _design(self):
        mass_solids = 0
        solids = self._solids
        for s in self.ins:
            mass_solids += s.imass[solids]
        ts = np.asarray(mass_solids).sum() # Total solids
        ts *= 0.0011023 # To short tons (2000 lbs/hr)
        self.design_results['Solids loading'] = ts
        lb = self.minimum_solids_loading
        if ts < lb:
            lb_warning(self, 'Solids loading', ts, 'ton/hr', lb)
        self.design_results['Flow rate'] = F_vol_in = self.F_vol_in
        self.power_utility(F_vol_in * self.kWhr_per_m3)
            
    
