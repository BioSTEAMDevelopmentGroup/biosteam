# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import biosteam as bst
from ._lle_unit import LLEUnit
from ._splitter import Splitter
from .design_tools import PressureVessel

__all__ = ('LiquidsSettler', 'LLESettler', 'LiquidsSplitSettler')

class LiquidsSettler(bst.Unit, PressureVessel, isabstract=True):
    _N_ins = 1
    _N_outs = 2
    _N_heat_utilities = 0
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 area_to_feed=0.1, 
                 length_to_diameter=4,
                 vessel_material='Carbon steel',
                 vessel_type='Horizontal'):
        bst.Unit.__init__(self, ID, ins, outs, thermo)
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type
        self.length_to_diameter = length_to_diameter #: Length to diameter ratio
        self.area_to_feed = area_to_feed #: [ft2/gpm] Diameter * length per gpm of feed
    
    @staticmethod
    def _default_vessel_type():
        return 'Horizontal'
    
    def _design(self):
        feed = self.ins[0]
        F_vol_gpm = feed.get_total_flow('gpm')
        area = self.area_to_feed * F_vol_gpm
        length_to_diameter = self.length_to_diameter
        P = feed.get_property('P', 'psi')
        D = (area / length_to_diameter) ** 0.5
        L = length_to_diameter * D
        self.design_results.update(self._vessel_design(P, D, L))
        
    def _cost(self):
        D = self.design_results
        self.purchase_costs.update(
            self._vessel_purchase_cost(D['Weight'], D['Diameter'], D['Length'])
        )
        
class LLESettler(LLEUnit, LiquidsSettler):
    line = 'Settler'
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 area_to_feed=0.1, 
                 length_to_diameter=4,
                 vessel_material='Carbon steel',
                 vessel_type='Horizontal',
                 top_chemical=None,
                 efficiency=1.0):
        LLEUnit.__init__(self, ID, ins, outs, thermo, top_chemical, efficiency)
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type
        self.length_to_diameter = length_to_diameter
        self.area_to_feed = area_to_feed
    
class LiquidsSplitSettler(LiquidsSettler):
    line = 'Settler'
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 split, order=None,
                 area_to_feed=0.1, 
                 length_to_diameter=4,
                 vessel_material='Carbon steel',
                 vessel_type='Horizontal'):
        bst.Unit.__init__(self, ID, ins, outs, thermo)
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type
        self.length_to_diameter = length_to_diameter
        self.area_to_feed = area_to_feed
        self._isplit = self.chemicals.isplit(split, order)        
    split = Splitter.split
    isplit = Splitter.isplit
    _run = Splitter._run