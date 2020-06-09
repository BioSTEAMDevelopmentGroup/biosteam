# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from ._splitter import Splitter
from ..exceptions import InfeasibleRegion
from .design_tools import separations

__all__ = ('SolidsSeparator',)

class SolidsSeparator(Splitter):
    """
    Create SolidsSeparator object.
    
    Parameters
    ----------
    ins : streams
        Inlet fluids with solids.
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
        self.moisture_content = moisture_content
        assert self.isplit['7732-18-5'] == 0, 'cannot define water split, only moisture content'
    
    def _run(self):
        separations.mix_and_split_with_moisture_content(
            self.ins, *self.outs, self.split, self.moisture_content
        )
        