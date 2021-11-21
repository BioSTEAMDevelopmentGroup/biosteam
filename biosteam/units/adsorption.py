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
from .splitting import Splitter
from biosteam.units.design_tools import PressureVessel
from math import sqrt, pi

__all__ = ('AdsorptionColumnTSA',)

class AdsorptionColumnTSA(PressureVessel, Splitter):
    """
    Create a temperature swing adsorption (TSA) column. Default parameters
    are heuristic values for adsorption of water and polar components using
    silica gel. The design and cost algorithm of the adsorption column is 
    based on [1]_, [2]_, and [4]_. Three vessels are used: a lead, a guard, 
    and a standby vessel. The lead is in equilibrium, guard holds the 
    breakthrough curve and mass transfer zone, and the standby is in regeneration.
    Once the standby vessel is regenerated, it becomes the guard, the lead becomes
    the standby, and the guard becomes the lead. Therefore, the cycle time is
    the regeneration time.
    
    Parameters
    ----------
    ins : stream sequence
        [0] Feed
        [1] Air
    outs : stream sequence
        [0] Effluent
        [1] Purge
    mean_velocity : float, optional
        Mean velocity of the feed. The diameter of the receiving vessel adjusts
        accordingly. Defaults to 7.2 [m / hr]. Typical velocities are 4 to 14.4 m / hr for liquids [1]_.
    regeneration_air_velocity : float, optional
        Mean velocity of the air used to regenerate the bed. Defaults to 540 [m / hr]. 
        Common velocity range for gasses is 504-2160 m / hr [1]_.
    cycle_time : float, optional
        Time at which the receiving vessel is switched. Defaults to 3 [hr]. 
        Note that 1-2 hours required for thermal-swing-adsorption 
        (TSA) for silica gels. One hr is added to be conservative [2]_. 
    rho_adsorbent : float, optional
        The density of the adsorbent. Defaults to 480 [kg/m3], which is common 
        for silica gels [3]_.
    adsorbent_capacity : float, optional
        Fraction of absorbate that the absorbent can hold. Defaults to 0.1, a 
        conservative heuristic from [2]_.
    T_regeneration : float, optional
        Temperature during the regeneration phase. Defaults to 393.15 [K], which
        is used for silica gels [2]_.
    vessel_material : float, optional
        Vessel material. Defaults to 'Stainless steel 316',
    vessel_type : float, optional
        Vessel type. Defaults to 'Vertical'.
    adsorbate_ID : string
        Name of adsorbate.
    order : tuple[str], optional
        Order of component splits.
    split : dict[str, float] or list[float], optional
        Component splits towards the effluent (0th outlet).
    
    Examples
    --------
    >>> import biosteam as bst
    >>> bst.settings.set_thermo(['Water', 'O2', 'N2', 'Hexane'])
    >>> feed = bst.Stream('feed', Hexane=0.9, Water=0.1, units='kg/hr', total_flow=1000.)
    >>> A1 = bst.AdsorptionColumnTSA('A1', [feed, 'air'], 
    ...     split=dict(Water=0., Hexane=1.0),
    ...     adsorbate_ID='Water',
    ...  )
    >>> A1.simulate()
    >>> A1.results()
    Adsorption column TSA                           Units                   A1
    Low pressure steam  Duty                        kJ/hr             2.75e+05
                        Flow                      kmol/hr                  7.1
                        Cost                       USD/hr                 1.69
    Design              Vessel diameter                                   0.51
                        Vessel length                                     12.6
                        Number of reactors                                   3
                        Vessel type                                   Vertical
                        Length                         ft                 3.83
                        Diameter                       ft                0.156
                        Weight                         lb                 19.9
                        Wall thickness                 in                 0.25
                        Vessel material                    Stainless steel 316
    Purchase cost       Vertical pressure vessel      USD             1.68e+04
                        Platform and ladders          USD                  803
    Total purchase cost                               USD             1.76e+04
    Utility cost                                   USD/hr                 1.69
    
    References
    ----------
    [1] Adsorption basics Alan Gabelman (2017) Adsorption basics Part 1. AICHE
    [2] Seader, J. D., Separation Process Principles: Chemical and Biochemical Operations,” 3rd ed., Wiley, Hoboken, NJ (2011).
    [3] https://www.daisogelusa.com/technical-notes/approximate-packing-density-for-daisogel-bulk-silica-gel/
    [4] Seider, W. D., Lewin,  D. R., Seader, J. D., Widagdo, S., Gani,
        R., & Ng, M. K. (2017). Product and Process Design Principles. Wiley.
        Cost Accounting and Capital Cost Estimation (Chapter 16)
    
    """
    # NOTE: Unit ignores cost of compressing air.
    _N_ins = 2
    _N_outs = 2
    _N_heat_utilities = 1
    
    def __init__(self, 
            ID='', ins=None, outs=(), thermo=None, *,
            mean_velocity=7.2, # m / hr; typical velocities are 4 to 14.4 m /hr for liquids; Adsorption basics Alan Gabelman (2017) Adsorption basics Part 1. AICHE
            regeneration_air_velocity=1332, # Mid point in velocity range for gasses, m / hr; Alan Gabelman (2017) Adsorption basics Part 1. AICHE
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
        F_mass_adsorbate = float(feed.imass[self.adsorbate_ID])
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
    
    def _cost(self):
        design_results = self.design_results
        baseline_purchase_costs = self.baseline_purchase_costs
        baseline_purchase_costs.update(self._vessel_purchase_cost(
            design_results['Weight'], design_results['Diameter'], design_results['Length']))
        for i, j in baseline_purchase_costs.items():
            baseline_purchase_costs[i] *= design_results['Number of reactors']