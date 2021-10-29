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

__all__ = ('AdsorbtionColumnTSA',)

class AdsorbtionColumnTSA(PressureVessel, bst.Splitter):
    """
    Examples
    --------
    >>> import biosteam as bst
    >>> bst.settings.set_thermo(['Water', 'O2', 'N2'])
    >>> feed = bst.Stream('feed', O2=0.32, N2=0.78, Water=0.1, units='kg/hr', total_flow=1000.)
    >>> A1 = AdsorbtionColumnTSA('A1', [feed, 'air'], 
    ...     split=dict(Water=0., O2=0.99, N2=0.99),
    ...     adsorbate_ID='Water',
    ...  )
    >>> A1.simulate()
    >>> A1.results()
    Adsorbtion column TSA                           Units                   A1
    Low pressure steam  Duty                        kJ/hr              2.3e+05
                        Flow                      kmol/hr                 5.94
                        Cost                       USD/hr                 1.41
    Design              Vessel diameter                                  0.494
                        Vessel length                                     10.8
                        Number of reactors                                   3
                        Vessel type                                   Vertical
                        Length                         ft                  3.3
                        Diameter                       ft                0.151
                        Weight                         lb                 16.7
                        Wall thickness                 in                 0.25
                        Vessel material                    Stainless steel 316
    Purchase cost       Vertical pressure vessel      USD             1.52e+04
                        Platform and ladders          USD                  675
    Total purchase cost                               USD             1.59e+04
    Utility cost                                   USD/hr                 1.41
    
    """
    # NOTE: Unit ignores cost of compressing air.
    _N_ins = 2
    _N_outs = 2
    _N_heat_utilities = 1
    
    def __init__(self, 
            ID='', ins=None, outs=(), thermo=None, *,
            mean_velocity=14.4, # m / hr; Alan Gabelman (2017) Adsorption basics Part 1. AICHE
            regeneration_air_velocity=540, # Lower bound common velocity range for gasses, m / hr; Alan Gabelman (2017) Adsorption basics Part 1. AICHE
            cycle_time=3, # 1-2 hours required for thermal-swing-adsorption (TSA) for silica gels (add 1 hr for conservativeness); Seader, J. D., Separation Process Principles: Chemical and Biochemical Operations,” 3rd ed., Wiley, Hoboken, NJ (2011).
            rho_adsorbent=480, # (in kg/m3) Common for silica gels https://www.daisogelusa.com/technical-notes/approximate-packing-density-for-daisogel-bulk-silica-gel/
            adsorbent_capacity=0.1, # Conservative heuristic from Seider et. al. (2017) Product and Process Design Principles. Wiley
            T_regeneration=120 + 273.15, # For silica gels; Seader, J. D., Separation Process Principles: Chemical and Biochemical Operations,” 3rd ed., Wiley, Hoboken, NJ (2011).
            vessel_material='Stainless steel 316',
            vessel_type='Vertical',
            adsorbate_ID, 
            order=None, 
            split,
        ):
        bst.Splitter.__init__(self, ID, ins, outs, thermo, order=order, split=split)
        self.mean_velocity = mean_velocity
        self.regeneration_air_velocity = regeneration_air_velocity
        self.cycle_time = cycle_time
        self.rho_adsorbent = rho_adsorbent
        self.adsorbent_capacity = adsorbent_capacity
        self.adsorbate_ID = adsorbate_ID
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type
        self.T_regeneration = T_regeneration
        
    @property
    def effluent(self):
        return self.outs[0]
        
    @property
    def regeneration_purge(self):
        return self.outs[1]
    
    def _run(self):
        feed, air = self.ins
        effluent, purge = self.outs 
        feed.split_to(effluent, purge, self.split)
        F_vol_feed = feed.F_vol
        mean_velocity = self.mean_velocity
        self.design_results['Vessel diameter'] = diameter = 2 * sqrt(F_vol_feed / (mean_velocity * pi))
        regeneration_air_velocity = self.regeneration_air_velocity
        diameter = 2 * sqrt(F_vol_feed / (mean_velocity * pi))
        radius = diameter * 0.5
        F_vol_air = radius * radius * regeneration_air_velocity * pi
        purge.T = air.T = self.T_regeneration
        air.P = 10.
        air.imass['N2', 'O2'] = [0.78, 0.32]
        air.phase = 'g'
        air.F_vol = F_vol_air
        purge.mol += air.mol
        purge.phase = 'g'
        
    @property
    def online_time(self):
        return 2 * self.cycle_time # 3 columns (1 quard, 1 active, 1 regenerating)
    
    def _design(self):
        design_results = self.design_results
        feed = self.ins[0]
        rho_adsorbent = self.rho_adsorbent
        adsorbent_capacity = self.adsorbent_capacity
        F_mass_adsorbate = feed.imass[self.adsorbate_ID].value
        diameter = design_results['Vessel diameter']
        online_length = (
            self.online_time * F_mass_adsorbate / (adsorbent_capacity * rho_adsorbent * diameter)
        ) + 0.61 # equilibrium length plust 2 ft accounting for mass transfer limitations
        design_results['Vessel length'] = length = online_length / 2. # Size of each column
        design_results['Number of reactors'] = 3
        design_results.update(
            self._vessel_design(
                feed.P * 0.000145038, # Pa to psi
                float(diameter) / 3.28084, # m to ft
                float(length) / 3.28084, # m to ft
            )
        )    
        self.heat_utilities[0](self.regeneration_purge.H, 273.15 + 25, self.T_regeneration)
        # TODO: Think about heating duty associated to flowing dry air for regeneration.
        # TODO: Add tests for sizing and costs for the pressure vessels
    
    def _cost(self):
        design_results = self.design_results
        baseline_purchase_costs = self.baseline_purchase_costs
        baseline_purchase_costs.update(self._vessel_purchase_cost(
            design_results['Weight'], design_results['Diameter'], design_results['Length']))
        for i, j in baseline_purchase_costs.items():
            baseline_purchase_costs[i] *= design_results['Number of reactors']