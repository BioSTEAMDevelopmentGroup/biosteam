# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import biosteam as bst
from biosteam.units.design_tools import PressureVessel
from math import sqrt, pi

class AdsorbtionColumnTSA(bst.Splitter, PressureVessel):
    _N_ins = 1
    _N_outs = 2
    _N_heat_utilities = 1
    
    def __init__(self, 
            ID='', ins=None, outs=(), thermo=None, *,
            mean_velocity=4., # Alan Gabelman (2017) Adsorption basics Part 1. AICHE
            cycle_time=3, # 1-2 hours required for thermal-swing-adsorption (TSA) for silica gels (add 1 hr for conservativeness); Seader, J. D., Separation Process Principles: Chemical and Biochemical Operations,” 3rd ed., Wiley, Hoboken, NJ (2011).
            rho_adsorbent=480, # (in kg/m3) Common for silica gels https://www.daisogelusa.com/technical-notes/approximate-packing-density-for-daisogel-bulk-silica-gel/
            adsorbent_capacity=0.1, # Conservative heuristic from Seider et. al. (2017) Product and Process Design Principles. Wiley
            T_regeneration=120 + 273.15, # For silica gels; Seader, J. D., Separation Process Principles: Chemical and Biochemical Operations,” 3rd ed., Wiley, Hoboken, NJ (2011).
            vessel_material='Stainless steel 316',
            vessel_type='Vertical',
            adsorbent_ID, 
            order=None, 
            split,
        ):
        bst.Splitter.__init__(self, ID, ins, outs, thermo, order=order, split=split)
        self.mean_velocity = mean_velocity
        self.cycle_time = cycle_time
        self.rho_adsorbent = rho_adsorbent
        self.adsorbent_capacity = adsorbent_capacity
        self.adsorbent_ID = adsorbent_ID
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type
        
    @property
    def effluent(self):
        return self.outs[0]
        
    @property
    def regeneration_purge(self):
        return self.outs[1]
    
    def _run(self):
        bst.Splitter._run(self)
        self.outs[1].T = T_regeneration
        
    @property
    def online_time(self):
        return 2 * self.cycle_time # 3 columns (1 quard, 1 active, 1 regenerating)
    
    def _design(self):
        design_results = self.design_results
        feed = self.ins[0]
        rho_adsorbent = self.rho_adsorbent
        mean_velocity = self.mean_velocity
        adsorbent_capacity = self.adsorbent_capacity
        F_vol_feed = feed.F_vol
        F_mass_adsorbate = feed.imass[self.adsorbent_ID].value
        diameter = 2 * sqrt(F_vol_feed / (mean_velocity * pi))
        length = (
            online_time * F_mass_adsorbate / (adsorbent_capacity * rho_adsorbent * diameter)
        ) + 0.61 # equilibrium length plust 2 ft accounting for mass transfer limitations
        design_results['Length'] = length
        design_results['Diameter'] = diameter
    
        P = feed.P * 0.000145038 # Pa to psi
        design_results['Number of reactors'] = 3
        design_results.update(self._vessel_design(float(P), float(diameter), float(length)))    
    
        # TODO: Think about heating duty associated to flowing dry air for regeneration.
        # TODO: Add tests for sizing and costs for the pressure vessels
    
    def _cost(self):
        design_results = self.design_results
        baseline_purchase_costs = self.baseline_purchase_costs
        baseline_purchase_costs.update(self._vessel_purchase_cost(
            design_results['Weight'], design_results['Diameter'], design_results['Length']))
        for i, j in baseline_purchase_costs.items():
            baseline_purchase_costs[i] *= design_results['Number of reactors']