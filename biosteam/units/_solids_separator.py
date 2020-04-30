# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 22:14:01 2018

@author: yoelr
"""
from ._splitter import Splitter, run_split_with_mixing
from ..exceptions import InfeasibleRegion

__all__ = ('SolidsSeparator',)

class SolidsSeparator(Splitter):
    """Create SolidsSeparator object.
    
    Parameters
    ----------
    ins : stream 
        Inlet fluid with solids.
    outs : stream sequnece
        * [0] Retentate.
        * [1] Permeate.
    split : array_like
           Component splits to 0th output stream
    moisture_content : float
                       Fraction of water in solids
    
    """
    _N_ins = 2
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 order=None, split, moisture_content):
        Splitter.__init__(self, ID, ins, outs, thermo, order=order, split=split)
        #: Moisture content of retentate
        self.mositure_content = moisture_content
        assert self.isplit['7732-18-5'] == 0, 'cannot define water split, only moisture content'
    
    def _run(self):
        run_split_with_mixing(self)
        retentate, permeate = self.outs
        solids = retentate.F_mass
        mc = self.mositure_content
        retentate.imol['7732-18-5'] = water = (solids * mc/(1-mc))/18.01528
        permeate.imol['7732-18-5'] -= water
        if permeate.imol['7732-18-5'] < 0:
            raise InfeasibleRegion('not enough water; permeate moisture content')
    
