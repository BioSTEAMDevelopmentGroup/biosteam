# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from .._unit import Unit
from ..utils import streams_from_units

class SummaryUnit(Unit):
    line = ''
    _N_ins = 0
    _N_outs = 0
    
    def __init__(self, ID='', thermo=None, *, units):
        ins = []
        outs = []
        self.units = set(units)
        for s in streams_from_units(units):
            source = s._source
            sink = s._sink
            if source in units and sink not in units:
                outs.append(s)
            elif sink in units and source not in units:
                ins.append(s)
        super().__init__(ID, ins, outs, thermo)
    
    @property
    def auxiliary_unit_names(self):
        return [i.ID for i in self.units]
    
    @property
    def auxiliary_units(self):
        return self.units
    
    def _run(self): pass