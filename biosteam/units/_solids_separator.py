# -*- coding: utf-8 -*-
<<<<<<< HEAD
"""
Created on Thu Aug 23 22:14:01 2018

@author: yoelr
"""
from ._splitter import Splitter, run_split_with_mixing
from ..exceptions import InfeasibleRegion
=======
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
>>>>>>> cd2c5013aaf9b5bc94bb764b52fd37db183472f1

__all__ = ('SolidsSeparator',)

class SolidsSeparator(Splitter):
<<<<<<< HEAD
    """Create SolidsSeparator object.
    
    Parameters
    ----------
    ins : stream 
        Inlet fluid with solids.
=======
    """
    Create SolidsSeparator object.
    
    Parameters
    ----------
    ins : streams
        Inlet fluids with solids.
>>>>>>> cd2c5013aaf9b5bc94bb764b52fd37db183472f1
    outs : stream sequnece
        * [0] Retentate.
        * [1] Permeate.
    split : array_like
<<<<<<< HEAD
           Component splits to 0th output stream
    moisture_content : float
                       Fraction of water in solids
=======
        Component splits to 0th output stream
    moisture_content : float
        Fraction of water in solids
>>>>>>> cd2c5013aaf9b5bc94bb764b52fd37db183472f1
    
    """
    _N_ins = 2
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 order=None, split, moisture_content):
        Splitter.__init__(self, ID, ins, outs, thermo, order=order, split=split)
        #: Moisture content of retentate
<<<<<<< HEAD
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
    
=======
        self.moisture_content = moisture_content
        assert self.isplit['7732-18-5'] == 0, 'cannot define water split, only moisture content'
    
    def _run(self):
        separations.mix_and_split_with_moisture_content(
            self.ins, *self.outs, self.split, self.moisture_content
        )
        
>>>>>>> cd2c5013aaf9b5bc94bb764b52fd37db183472f1
