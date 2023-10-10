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
from math import sqrt, pi, ceil
import numpy as np
from numba import njit
import flexsolve as flx

__all__ = ('AdsorptionColumnTSA',)


@njit(cache=True)
def f_recovery_objective(
        superficial_velocity, target_recovery, F_mass_adsorbate, F_vol_feed, cycle_time, drying_time,
        adsorbent_capacity, rho_adsorbent, void_fraction, length_unused, regeneration_velocity,
        K
    ):
    diameter = 2 * sqrt(F_vol_feed / (superficial_velocity * pi))
    area = pi * diameter * diameter / 4
    total_length = (
         cycle_time * F_mass_adsorbate / (adsorbent_capacity * rho_adsorbent * area)
    ) + length_unused # length of equilibrium section plus unused bed (LES + LUB)
    length = total_length / 2 # Size of each column
    F_vol_regen = area * regeneration_velocity
    if K is not None:
        vessel_volume = length * 0.25 * diameter * diameter
        void_volume = void_fraction * vessel_volume
        N_washes = ceil(F_vol_regen * (cycle_time - drying_time) / void_volume)
        solvent = void_volume * 1e6 # m3 -> mL
        adsorbent = 1000 * vessel_volume * rho_adsorbent # g
        total_adsorbate = adsorbate = adsorbent * adsorbent_capacity # g
        stages = ceil(length / 0.1) # 0.1 m per stage
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
        recovery = 1 - adsorbate_arr.sum() / total_adsorbate
    
    if target_recovery is None:
        error = 0
    else:
        error = target_recovery - recovery
    return (error, diameter, area, length, 
            F_vol_regen, vessel_volume, void_volume,
            stages, N_washes, recovery)


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
    ins : 
        * [0] Feed
        * [1] Regeneration fluid
    outs : 
        * [0] Effluent
        * [1] Purge
    superficial_velocity : float, optional
        Superficial velocity of the feed. The diameter of the receiving vessel adjusts
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
    >>> A1.results(with_units=False)
    Low pressure steam   Duty                                                    1.68e+03
                         Flow                                                      0.0436
                         Cost                                                      0.0104
    Design               Number of reactors                                             3
                         Vessel type                                             Vertical
                         Length                                                      52.1
                         Diameter                                                    1.68
                         Weight                                                   2.9e+03
                         Wall thickness                                              0.25
                         Vessel material                              Stainless steel 316
    Purchase cost        Vertical pressure vessel                                1.47e+05
                         Platform and ladders                                    2.95e+04
                         Silica gel                                              7.23e+04
                         Heat exchanger regeneration - Double pipe                    117
    Total purchase cost                                                          2.49e+05
    Utility cost                                                                   0.0104
    Name: A1, dtype: object
    
    References
    ----------
    [1] Adsorption basics Alan Gabelman (2017) Adsorption basics Part 1. AICHE
    [2] Seader, J. D., Separation Process Principles: Chemical and Biochemical Operations,” 3rd ed., Wiley, Hoboken, NJ (2011).
    [3] https://www.daisogelusa.com/technical-notes/approximate-packing-density-for-daisogel-bulk-silica-gel/
    [4] Seider, W. D., Lewin,  D. R., Seader, J. D., Widagdo, S., Gani,
        R., & Ng, M. K. (2017). Product and Process Design Principles. Wiley.
        Cost Accounting and Capital Cost Estimation (Chapter 16)
    
    """
    auxiliary_unit_names = ('heat_exchanger_regeneration',
                            'heat_exchanger_drying')
    
    # In $/ft3
    adsorbent_cost = {
        'Activated alumina': 72,
        'Activated carbon': 41,
        'Silica gel': 210,
        'Molecular sieves': 85,
    }
    
    # TODO: Update later plus ref
    _default_equipment_lifetime = {
        'Activated alumina': 10,
        'Activated carbon': 10,
        'Silica gel': 10,
        'Molecular sieves': 10,
    }
    
    # NOTE: Unit ignores cost of compressing 2.
    _N_ins = 3
    _N_outs = 3
    
    def _init(self, 
            adsorbate_ID, split,
            superficial_velocity=7.2, # m / hr; typical velocities are 4 to 14.4 m /hr for liquids; Adsorption basics Alan Gabelman (2017) Adsorption basics Part 1. AICHE
            regeneration_velocity=1332, # Mid point in velocity range for gasses, m / hr; Alan Gabelman (2017) Adsorption basics Part 1. AICHE
            cycle_time=3, # 1-2 hours required for thermal-swing-adsorption (TSA) for silica gels (add 1 hr for conservativeness); Seader, J. D., Separation Process Principles: Chemical and Biochemical Operations,” 3rd ed., Wiley, Hoboken, NJ (2011).
            rho_adsorbent=480, # Bulk density, including void fraction (in kg/m3) Common for silica gels https://www.daisogelusa.com/technical-notes/approximate-packing-density-for-daisogel-bulk-silica-gel/
            rho_adsorbent_solid = None, # Solid density, excluding void fraction (in kg/m3); used only if rho_adsorbent = None or False or 0
            adsorbent_capacity=0.1, # Conservative heuristic from Seider et. al. (2017) Product and Process Design Principles. Wiley
            T_regeneration=30 + 273.15, # For silica gels; Seader, J. D., Separation Process Principles: Chemical and Biochemical Operations,” 3rd ed., Wiley, Hoboken, NJ (2011).
            vessel_material='Stainless steel 316',
            vessel_type='Vertical',
            regeneration_fluid=dict(N2=0.78, O2=0.32, phase='g', units='kg/hr'),
            void_fraction=0.47, # Only matters when K given; 0.30 - 0.35 for activated carbon
            length_unused=1.219, # Additional length of a column to account for mass transfer limitations (due to unused bed). Defaults to +2 ft per column.
            drying_time=0., # Time for drying after regeneration
            T_air = 100 + 273.15,
            adsorbent='Silica gel',
            air_velocity = 1332,
            target_recovery=None,
            K=None,
            converge_adsorption_recovery=False,
            order=None, 
            wet_retention=0.01,
        ):
        bst.Splitter._init(self, order=order, split=split)
        self.superficial_velocity = superficial_velocity
        self.cycle_time = cycle_time
        self.adsorbent_capacity = adsorbent_capacity
        self.adsorbate_ID = adsorbate_ID
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type
        self.T_regeneration = T_regeneration
        self.regeneration_velocity = regeneration_velocity
        self.regeneration_fluid = regeneration_fluid
        self.void_fraction = void_fraction
        self.length_unused = length_unused
        self.T_air = T_air
        self.air_velocity = air_velocity
        self.drying_time = drying_time
        self.converge_adsorption_recovery = converge_adsorption_recovery
        self.wet_retention = wet_retention
        self.target_recovery = target_recovery
        self.adsorbent = adsorbent
        self.K = K
        self.rho_adsorbent_solid = rho_adsorbent_solid
        self.rho_adsorbent = rho_adsorbent if rho_adsorbent else rho_adsorbent_solid * (1-void_fraction)
        self.heat_exchanger_regeneration = bst.HXutility(None, None, None, thermo=self.thermo)
        self.heat_exchanger_drying = bst.HXutility(None, None, None, thermo=self.thermo)
        
    @property
    def effluent(self):
        return self.outs[0]
        
    @property
    def regeneration_purge(self):
        return self.outs[1]
    
    def _run(self):
        feed, regen, dry_air = self.ins
        effluent, purge, air_purge = self.outs 
        regen.empty()
        dry_air.empty()
        for i in self.outs: i.empty()

        feed.split_to(effluent, purge, self.split)
        F_vol_feed = feed.F_vol
        superficial_velocity = self.superficial_velocity
        adsorbate_ID = self.adsorbate_ID
        F_mass_adsorbate = purge.imass[adsorbate_ID]
        target_recovery = self.target_recovery
        K = self.K # (g adsorbate / mL solvent)  /  (g adsorbate / g adsorbent)
        if K is not None:
            args = (
                target_recovery, F_mass_adsorbate, 
                F_vol_feed, self.cycle_time, self.drying_time,
                self.adsorbent_capacity, self.rho_adsorbent, self.void_fraction, 
                self.length_unused, self.regeneration_velocity,
                K
            )
            f = lambda *args: f_recovery_objective(*args)[0]
            if target_recovery is not None:
                if (y1:=f(14.4, *args)) <= 0.:
                    self.superficial_velocity = superficial_velocity = 14.4 # typical velocities are 4 to 14.4 m /hr for liquids; Adsorption basics Alan Gabelman (2017) Adsorption basics Part 1. AICHE
                elif (y0:=f(0.1, *args)) >= 0.:
                    self.superficial_velocity = superficial_velocity = 0.1 # set minimum velocity
                else:
                    self.superficial_velocity = superficial_velocity = flx.IQ_interpolation(
                        f, 0.1, 14.4, 
                        # x=self.superficial_velocity,
                        # x=7.25,
                        x=9.20,
                        # xtol=1e-6, ytol=1e-9, 
                        xtol=1e-5, ytol=1e-7, 
                        y0=y0, y1=y1,
                        args=args,
                    )
            (_, self.diameter, self.area, self.length,  
            self._F_vol_regen, self.vessel_volume, self.void_volume,
            self.equilibrium_stages, self.N_washes, self.recovery) = f_recovery_objective(superficial_velocity, *args)
            T_original = regen.T
            regen.reset_flow(**self.regeneration_fluid)
            purge.T = regen.T = self.T_regeneration
            regen.F_vol = self._F_vol_regen
            regen.T = T_original
            TAL = feed.imol[adsorbate_ID]
            split = self.isplit[adsorbate_ID]
            self.actual_split = actual_split = 1 - (self.recovery * (1 - split))
            effluent.imol[adsorbate_ID] = actual_split * TAL
            purge.mol += regen.mol
            purge.imol[adsorbate_ID] = feed.imol[adsorbate_ID] - effluent.imol[adsorbate_ID]
        else:
            self.diameter = diameter = 2 * sqrt(F_vol_feed / (superficial_velocity * pi))
            self.area = area = pi * diameter * diameter / 4
            total_length = (
                self.cycle_time * F_mass_adsorbate / (self.adsorbent_capacity * self.rho_adsorbent * area)
            ) + self.length_unused # length of equilibrium section plus unused bed (LES + LUB)
            self.length = length = total_length / 2 # Size of each column
            self.vessel_volume = length * area
            T_original = regen.T
            regen.reset_flow(**self.regeneration_fluid)
            purge.T = regen.T = self.T_regeneration
            regen.F_vol = area * self.regeneration_velocity
            regen.T = T_original
        
        if self.drying_time:
            air_purge.empty()
            air_purge.T = self.T_air
            air_purge.P = 10 * 101325
            air_purge.imass['N2'] = 1
            air_purge.phase = dry_air.phase = 'g'
            air_purge.F_vol = self.area * self.air_velocity * self.drying_time / self.cycle_time
            dry_air.copy_like(air_purge)
            retained_ethanol_vol = self.wet_retention * self.void_volume 
            retained_ethanol_mol = retained_ethanol_vol * (regen.mol / regen.F_vol)
            ethanol_mol = retained_ethanol_mol / self.cycle_time
            air_purge.mol += ethanol_mol
            H_out = air_purge.H
            regen._property_cache.clear()
            H_in0 = (regen.H / regen.F_vol) * retained_ethanol_vol # H_solvent 
            try:
                dry_air.H = H_out - H_in0
            except:
                breakpoint()
            purge.mol -= ethanol_mol
    
    def _design(self):
        feed, regen, dry_air = self.ins
        design_results = self.design_results
        diameter = self.diameter
        length = self.length
        design_results['Number of reactors'] = 3
        design_results.update(
            self._vessel_design(
                feed.P * 0.000145038, # Pa to psi
                diameter * 3.28084, # m to ft
                length * 3.28084, # m to ft
            )
        )
        hxr = self.heat_exchanger_regeneration
        hxr.ins.empty()
        hxr.outs.empty()
        hxr.ins[0] = regen.copy()
        hxr.outs[0] = regen.copy()
        hxr.T = self.T_regeneration
        hxr.simulate()
        if self.drying_time:
            hxd = self.heat_exchanger_drying
            hxd.ins.empty()
            hxd.outs.empty()
            hxd.ins[0] = fresh_air = dry_air.copy()
            fresh_air.T = 298.15
            hxd.outs[0] = dry_air.copy()
            hxd.T = dry_air.T
            hxd.simulate()
        else:
            self.heat_exchanger_drying._setup()
    
    def _cost(self):
        design_results = self.design_results
        baseline_purchase_costs = self.baseline_purchase_costs
        baseline_purchase_costs.update(self._vessel_purchase_cost(
            design_results['Weight'], design_results['Diameter'], design_results['Length']))
        N_reactors = design_results['Number of reactors']
        for i, j in baseline_purchase_costs.items():
            baseline_purchase_costs[i] *= N_reactors
        baseline_purchase_costs[self.adsorbent] = N_reactors * 35.3147 * self.vessel_volume * self.adsorbent_cost[self.adsorbent]
