# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
@author: yoelcp and sarangbhagwat
"""
import biosteam as bst
from .splitting import Splitter
from biosteam.units.design_tools import PressureVessel
from thermosteam.equilibrium import phase_fraction
from math import sqrt, pi, ceil
import numpy as np
import flexsolve as flx

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
        [1] Regeneration fluid
    outs : stream sequence
        [0] Effluent
        [1] Purge
    mean_velocity : float, optional
        Mean velocity of the feed. The diameter of the receiving vessel adjusts
        accordingly. Defaults to 7.2 [m / hr]. Typical velocities are 4 to 14.4 m / hr for liquids [1]_.
    regeneration_velocity : float, optional
        Mean velocity of the fluid used to regenerate the bed. Defaults to 540 [m / hr]. 
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
    regeneration_fluid : dict[str, float]
        Arguments to initialize fluid used to regenerate the bed.
    split : dict[str, float] or list[float], optional
        Component splits towards the effluent (0th outlet).
    
    Examples
    --------
    >>> import biosteam as bst
    >>> bst.default()
    >>> bst.settings.set_thermo(['Water', 'O2', 'N2', 'Hexane'], cache=True)
    >>> feed = bst.Stream('feed', Hexane=0.9, Water=0.1, units='kg/hr', total_flow=1000.)
    >>> A1 = bst.AdsorptionColumnTSA('A1', [feed, 'air'], 
    ...     split=dict(Water=0., Hexane=1.0),
    ...     adsorbate_ID='Water',
    ...  )
    >>> A1.simulate()
    >>> A1.results()
    Adsorption column TSA                           Units                   A1
    Low pressure steam  Duty                        kJ/hr             2.75e+05
                        Flow                      kmol/hr                 7.07
                        Cost                       USD/hr                 1.68
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
    Utility cost                                   USD/hr                 1.68
    
    References
    ----------
    [1] Adsorption basics Alan Gabelman (2017) Adsorption basics Part 1. AICHE
    [2] Seader, J. D., Separation Process Principles: Chemical and Biochemical Operations,” 3rd ed., Wiley, Hoboken, NJ (2011).
    [3] https://www.daisogelusa.com/technical-notes/approximate-packing-density-for-daisogel-bulk-silica-gel/
    [4] Seider, W. D., Lewin,  D. R., Seader, J. D., Widagdo, S., Gani,
        R., & Ng, M. K. (2017). Product and Process Design Principles. Wiley.
        Cost Accounting and Capital Cost Estimation (Chapter 16)
    
    """
    # In $/ft3
    adsorbent_cost = {
        'Activated alumina': 72,
        'Activated carbon': 41,
        'Silica gel': 210,
        'Molecular sieves': 85,
    }
    
    # TODO: Update later plus ref
    _default_equipment_lifetime = {
        'Activated alumina': 1,
        'Activated carbon': 1,
        'Silica gel': 1,
        'Molecular sieves': 1,
    }
    
    # NOTE: Unit ignores cost of compressing 2.
    _N_ins = 3
    _N_outs = 3
    _N_heat_utilities = 1
    
    def __init__(self, 
            ID='', ins=None, outs=(), thermo=None, *,
            mean_velocity=7.2, # m / hr; typical velocities are 4 to 14.4 m /hr for liquids; Adsorption basics Alan Gabelman (2017) Adsorption basics Part 1. AICHE
            regeneration_velocity=1332, # Mid point in velocity range for gasses, m / hr; Alan Gabelman (2017) Adsorption basics Part 1. AICHE
            cycle_time=3, # 1-2 hours required for thermal-swing-adsorption (TSA) for silica gels (add 1 hr for conservativeness); Seader, J. D., Separation Process Principles: Chemical and Biochemical Operations,” 3rd ed., Wiley, Hoboken, NJ (2011).
            rho_adsorbent=480, # (in kg/m3) Common for silica gels https://www.daisogelusa.com/technical-notes/approximate-packing-density-for-daisogel-bulk-silica-gel/
            adsorbent_capacity=0.1, # Conservative heuristic from Seider et. al. (2017) Product and Process Design Principles. Wiley
            T_regeneration=30 + 273.15, # For silica gels; Seader, J. D., Separation Process Principles: Chemical and Biochemical Operations,” 3rd ed., Wiley, Hoboken, NJ (2011).
            vessel_material='Stainless steel 316',
            vessel_type='Vertical',
            regeneration_fluid=dict(N2=0.78, O2=0.32, phase='g', units='kg/hr'),
            void_fraction=0.35, # Only matters when K given; 0.30 - 0.35 for activated carbon
            length_plus=1.219, # Additional length of a column to account for mass transfer limitations (due to unused bed). Defaults to +2 ft per column.
            drying_time=0., # Time for drying after regeneration
            T_air = 100 + 273.15,
            adsorbent='Activated carbon',
            air_velocity = 1332,
            target_recovery=None,
            K=None,
            converge_adsorption_recovery=False,
            adsorbate_ID, 
            order=None, 
            wet_retention=0.01,
            split,
        ):
        bst.Splitter.__init__(self, ID, ins, outs, thermo, order=order, split=split)
        self.mean_velocity = mean_velocity
        self.cycle_time = cycle_time
        self.rho_adsorbent = rho_adsorbent
        self.adsorbent_capacity = adsorbent_capacity
        self.adsorbate_ID = adsorbate_ID
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type
        self.T_regeneration = T_regeneration
        self.regeneration_velocity = regeneration_velocity
        self.regeneration_fluid = regeneration_fluid
        self.void_fraction = void_fraction
        self.length_plus = length_plus
        self.T_air = T_air
        self.air_velocity = air_velocity
        self.drying_time = drying_time
        self.converge_adsorption_recovery = converge_adsorption_recovery
        self.wet_retention = wet_retention
        self.target_recovery = target_recovery
        self.adsorbent = adsorbent
        self.K = K
        
    @property
    def effluent(self):
        return self.outs[0]
        
    @property
    def regeneration_purge(self):
        return self.outs[1]
    
    def _run(self):
        feed, regen, dry_air = self.ins
        effluent, purge, air_purge = self.outs 
        feed.split_to(effluent, purge, self.split)
        F_vol_feed = feed.F_vol
        mean_velocity = self.mean_velocity
        rho_adsorbent = self.rho_adsorbent
        adsorbent_capacity = self.adsorbent_capacity
        adsorbate_ID = self.adsorbate_ID
        F_mass_adsorbate = purge.imass[adsorbate_ID]
        if self.K:
            target_recovery = self.target_recovery
            K = self.K # (g adsorbate / mL solvent)  /  (g adsorbate / g adsorbent)
            def f_efficiency(mean_velocity):
                self.diameter = diameter = 2 * sqrt(F_vol_feed / (mean_velocity * pi))
                self.area = area = pi * diameter * diameter / 4
                total_length = (
                    self.cycle_time * F_mass_adsorbate / (adsorbent_capacity * rho_adsorbent * area)
                ) + self.length_plus # length of equilibrium section plus unused bed (LES + LUB)
                self.length = length = total_length / 2 # Size of each column
                regeneration_velocity = self.regeneration_velocity
                self.radius = radius = diameter * 0.5
                self._F_vol_regen = F_vol_regen = radius * radius * regeneration_velocity * pi
                self.vessel_volume = vessel_volume = length * 0.25 * diameter * diameter
                self.void_volume = void_volume = self.void_fraction * vessel_volume
                N_washes = ceil(F_vol_regen * (self.cycle_time - self.drying_time) / void_volume)
                solvent = void_volume * 1e6 # m3 -> mL
                adsorbent = 1000 * vessel_volume * rho_adsorbent # g
                total_adsorbate = adsorbate = adsorbent * adsorbent_capacity # g
                self.equilibrium_stages = stages = ceil(length / 0.1) # 0.1 m per stage
                self.N_washes = N_washes
                solvent /= stages
                adsorbate_arr = adsorbate * np.ones(stages) / stages
                adsorbent_stage = adsorbent  / stages
                for i in range(stages * N_washes):
                    # R * y + A * x = total
                    # K = y / x
                    # R * y  + A * y / K = total
                    # y = total / (R + A / K)
                    adsorbate_collected = 0
                    for j in range(stages):
                        y = (adsorbate_arr[j] + adsorbate_collected) / (solvent + adsorbent_stage / K)
                        adsorbate_recovered = y * solvent - adsorbate_collected
                        adsorbate_arr[j] -= adsorbate_recovered
                        adsorbate_collected += adsorbate_recovered
                self.recovery = recovery = 1 - adsorbate_arr.sum() / total_adsorbate
                if target_recovery is None:
                    return mean_velocity
                return target_recovery - recovery
            if self.target_recovery is not None:
                if (y1:=f_efficiency(14.4)) <= 0.:
                    self.mean_velocity = mean_velocity = 14.4 # typical velocities are 4 to 14.4 m /hr for liquids; Adsorption basics Alan Gabelman (2017) Adsorption basics Part 1. AICHE
                elif (y0:=f_efficiency(0.1)) >= 0.:
                    self.mean_velocity = mean_velocity = 0.1 # set minimum velocity
                else:
                    self.mean_velocity = mean_velocity = flx.IQ_interpolation(f_efficiency, 0.1, 14.4, x=self.mean_velocity, xtol=0.01, ytol=0.001, y0=y0, y1=y1)
                
            else:
                self.mean_velocity = f_efficiency(mean_velocity)
            purge.T = regen.T = self.T_regeneration
            regen.reset_flow(**self.regeneration_fluid)
            regen.F_vol = self._F_vol_regen
            TAL = feed.imol[adsorbate_ID]
            split = self.isplit[adsorbate_ID]
            self.actual_split = actual_split = 1 - (self.recovery * (1 - split))
            effluent.imol[adsorbate_ID] = actual_split * TAL
            purge.mol += regen.mol
            purge.imol[adsorbate_ID] = feed.imol[adsorbate_ID] - effluent.imol[adsorbate_ID]
        else:
            purge.T = regen.T = self.T_regeneration
            regen.reset_flow(**self.regeneration_fluid)
            regen.F_vol = self._F_vol_regen
            diameter = 2 * sqrt(F_vol_feed / (self.mean_velocity * pi))
            radius = diameter * 0.5
            purge.mol += regen.mol
        
        if self.drying_time:
            radius = self.radius
            air_purge.empty()
            air_purge.T = self.T_air
            air_purge.P = 101325
            air_purge.imass['N2', 'O2'] = [0.78, 0.32]
            air_purge.phase = dry_air.phase = 'g'
            air_purge.F_vol= radius * radius * self.air_velocity * pi * self.drying_time / self.cycle_time
            dry_air.copy_like(air_purge)
            
            retained_ethanol_mol = self.wet_retention * regen.mol / self.N_washes
            air_purge.mol += retained_ethanol_mol
            H_out = air_purge.H
            regen._property_cache.clear()
            H_in0 = self.wet_retention * regen.H / self.N_washes
            try:
                dry_air.H = H_out - H_in0
            except:
                breakpoint()
            purge.mol -= retained_ethanol_mol
    
    def _design(self):
        feed = self.ins[0]
        design_results = self.design_results
        diameter = self.diameter
        length = self.length
        design_results['Number of reactors'] = 3
        design_results.update(
            self._vessel_design(
                feed.P * 0.000145038, # Pa to psi
                diameter / 3.28084, # m to ft
                length / 3.28084, # m to ft
            )
        )    
        self.heat_utilities[0](self.regeneration_purge.H, 273.15 + 25, self.T_regeneration)
    
    def _cost(self):
        design_results = self.design_results
        baseline_purchase_costs = self.baseline_purchase_costs
        baseline_purchase_costs.update(self._vessel_purchase_cost(
            design_results['Weight'], design_results['Diameter'], design_results['Length']))
        N_reactors = design_results['Number of reactors']
        for i, j in baseline_purchase_costs.items():
            baseline_purchase_costs[i] *= N_reactors
        baseline_purchase_costs[self.adsorbent] = N_reactors * 35.3147 * self.vessel_volume * self.adsorbent_cost[self.adsorbent]