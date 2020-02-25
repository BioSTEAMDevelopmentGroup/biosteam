# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 22:14:01 2018

@author: yoelr
"""
from ._splitter import Splitter, run_split_with_mixing
from ._CAS import H2O_CAS

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
    
    def __init__(self, ID='', ins=None, outs=(), *,
                 order=None, split, moisture_content):
        Splitter.__init__(self, ID, ins, outs, order=order, split=split)
        #: Moisture content of retentate
        self.mositure_content = moisture_content
        assert self.isplit[H2O_CAS] == 0, 'cannot define water split, only moisture content'
    
    def _run(self):
        run_split_with_mixing(self)
        retentate, permeate = self.outs
        solids = retentate.F_mass
        mc = self.mositure_content
        retentate.imol[H2O_CAS] = water = (solids * mc/(1-mc))/18.01528
        permeate.imol[H2O_CAS] -= water
        if permeate.imol[H2O_CAS] < water:
            raise ValueError(f'not enough water for {repr(self)}')
    
