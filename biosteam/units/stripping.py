# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>, Brenda Cansino <cansinoloeza@wisc.edu>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import numpy as np
from scipy.stats import gmean
import biosteam as bst
import thermosteam as tmo
from .design_tools.specification_factors import  (
    distillation_column_material_factors,
    tray_material_factor_functions,
    distillation_tray_type_factor,
    material_densities_lb_per_in3)
from . import design_tools as design
from .._graphics import vertical_column_graphics
from warnings import warn
from .phase_equilibrium import MultiStageEquilibrium

__all__ = ('AdiabaticMultiStageVLEColumn', 'Stripper', 'Absorber',)

class AdiabaticMultiStageVLEColumn(MultiStageEquilibrium):
    r"""
    Create an adsorption or stripping column without a reboiler/condenser. 
    The diameter is based on tray separation and flooding velocity. 
    Purchase costs are based on correlations compiled by Warren et. al.

    Parameters
    ----------
    ins :
        * [0] Liquid
        * [1] Vapor
    outs :
        * [0] Vapor
        * [1] Liquid
    P : float
        Operating pressure [Pa].
    vessel_material : str, optional
        Vessel construction material. Defaults to 'Carbon steel'.
    tray_material : str, optional
        Tray construction material. Defaults to 'Carbon steel'.
    tray_type='Sieve' : 'Sieve', 'Valve', or 'Bubble cap'
        Tray type.
    tray_spacing=450 : float
        Typically between 152 to 915 mm.
    stage_efficiency=None :
        User enforced stage efficiency. If None, stage efficiency is
        calculated by the O'Connell correlation [2]_.
    velocity_fraction=0.8 : float
        Fraction of actual velocity to maximum velocity allowable before
        flooding.
    foaming_factor=1.0 : float
        Must be between 0 to 1.
    open_tray_area_fraction=0.1 : float
        Fraction of open area to active area of a tray.
    downcomer_area_fraction=None : float
        Enforced fraction of downcomer area to net (total) area of a tray.
        If None, estimate ratio based on Oliver's estimation [1]_.

    Examples
    --------
    >>> import biosteam as bst
    >>> bst.settings.set_thermo(['AceticAcid', 'EthylAcetate', 'Water', 'MTBE'], cache=True)
    >>> feed = bst.Stream('feed', Water=75, AceticAcid=5, MTBE=20, T=320)
    >>> steam = bst.Stream('steam', Water=100, phase='g', T=390)
    >>> absorption = bst.Absorption("U1",
    ...     N_stages=2, ins=[feed, steam], 
    ...     solute="AceticAcid", outs=['vapor', 'liquid']
    ... )
    >>> absorption.simulate()
    >>> absorption.show()
    AdiabaticMultiStageVLEColumn: U1
    ins...
    [0] feed  
        phase: 'l', T: 320 K, P: 101325 Pa
        flow (kmol/hr): AceticAcid  5
                        Water       75
                        MTBE        20
    [1] steam  
        phase: 'g', T: 390 K, P: 101325 Pa
        flow (kmol/hr): Water  100
    outs...
    [0] vapor  
        phase: 'g', T: 366.33 K, P: 101325 Pa
        flow (kmol/hr): AceticAcid  3.72
                        Water       73.8
                        MTBE        20
    [1] liquid  
        phase: 'l', T: 372.87 K, P: 101325 Pa
        flow (kmol/hr): AceticAcid  1.28
                        Water       101
                        MTBE        0.000309
    
    >>> stripper.results()
    Absorption                                 Units       U1
    Design              Theoretical stages                  2
                        Actual stages                       4
                        Height                    ft     19.9
                        Diameter                  ft        3
                        Wall thickness            in    0.312
                        Weight                    lb 2.71e+03
    Purchase cost       Trays                    USD 4.37e+03
                        Tower                    USD 2.91e+04
                        Platform and ladders     USD 7.52e+03
    Total purchase cost                          USD  4.1e+04
    Utility cost                              USD/hr        0
    
    """
    _graphics = vertical_column_graphics
    _ins_size_is_fixed = False
    _N_ins = 2
    _N_outs = 2
    _units = {'Height': 'ft',
              'Diameter': 'ft',
              'Wall thickness': 'in',
              'Weight': 'lb'}
    _max_agile_design = (
        'Actual stages',       
        'Height',
        'Diameter',
        'Wall thickness',
        'Weight',
    )
    _F_BM_default = {
        'Platform and ladders': 1,
        'Tower': 4.3,
        'Trays': 4.3,
    }
   
    # [dict] Bounds for results
    _bounds = {'Diameter': (3., 24.),
               'Height': (27., 170.),
               'Weight': (9000., 2.5e6)}
   
    _side_draw_names = ('vapor_side_draws', 'liquid_side_draws')
    
    def _init(self,
            N_stages, feed_stages, 
            vapor_side_draws, liquid_side_draws,
            solute, # Needed to compute the Murphree stage efficiency 
            P=101325,  
            partition_data=None, 
            vessel_material='Carbon steel',
            tray_material='Carbon steel',
            tray_type='Sieve',
            tray_spacing=450,
            stage_efficiency=None,
            velocity_fraction=0.8,
            foaming_factor=1.0,
            open_tray_area_fraction=0.1,
            downcomer_area_fraction=None,  
            use_cache=None,
        ):
        super()._init(N_stages=N_stages, feed_stages=feed_stages,
                      top_side_draws=vapor_side_draws, 
                      bottom_side_draws=liquid_side_draws,
                      partition_data=partition_data,
                      phases=("g", "l"),
                      P=P, use_cache=use_cache)
       
        # Construction specifications
        self.solute = solute
        self.vessel_material = vessel_material
        self.tray_type = tray_type
        self.tray_material = tray_material
        self.tray_spacing = tray_spacing
        self.stage_efficiency = stage_efficiency
        self.velocity_fraction = velocity_fraction
        self.foaming_factor = foaming_factor
        self.open_tray_area_fraction = open_tray_area_fraction
        self.downcomer_area_fraction = downcomer_area_fraction
        self._last_args = (
            self.N_stages, self.feed_stages, self.vapor_side_draws, 
            self.liquid_side_draws, self.use_cache, *self._ins, 
            self.solvent_ID, self.partition_data, self.P
        )
        
    def _setup(self):
        super()._setup()
        args = (self.N_stages, self.feed_stages, self.vapor_side_draws, 
                self.liquid_side_draws, self.use_cache, *self._ins, 
                self.solvent_ID, self.partition_data, self.P)
        if args != self._last_args:
            MultiStageEquilibrium._init(
                self, N_stages=self.N_stages,
                feed_stages=self.feed_stages,
                phases=('g', 'l'), P=self.P,
                top_side_draws=self.vapor_side_draws, 
                bottom_side_draws=self.liquid_side_draws,
                partition_data=self.partition_data, 
                use_cache=self.use_cache, 
            )
            self._last_args = args
    
    def reset_cache(self, isdynamic=None):
        self._last_args = None
        
    @property
    def vapor(self):
        return self.outs[0]
    @property
    def liquid(self):
        return self.outs[1]   
     
    @property
    def tray_spacing(self):
        return self._TS
    @tray_spacing.setter
    def tray_spacing(self, TS):
        """Tray spacing (225-600 mm)."""
        self._TS = TS
   
    @property
    def stage_efficiency(self):
        """Enforced user defined stage efficiency."""
        return self._E_eff
    @stage_efficiency.setter
    def stage_efficiency(self, E_eff):
        self._E_eff = E_eff
   
    @property
    def velocity_fraction(self):
        """Fraction of actual velocity to maximum velocity allowable before flooding."""
        return self._f
    @velocity_fraction.setter
    def velocity_fraction(self, f):
        self._f = f
   
    @property
    def foaming_factor(self):
        """Foaming factor (0 to 1)."""
        return self._F_F
    @foaming_factor.setter
    def foaming_factor(self, F_F):
        if not 0 <= F_F <= 1:
            raise ValueError(f"foaming_factor must be between 0 and 1, ({F_F} given).")
        self._F_F = F_F
   
    @property
    def open_tray_area_fraction(self):
        """Fraction of open area, A_h, to active area, A_a."""
        return self._A_ha
    @open_tray_area_fraction.setter
    def open_tray_area_fraction(self, A_ha):
        self._A_ha = A_ha
   
    @property
    def downcomer_area_fraction(self):
        """Enforced fraction of downcomer area to net (total) area.
        If None, the fraction is estimated based on heuristics."""
        return self._A_dn
    @downcomer_area_fraction.setter
    def downcomer_area_fraction(self, A_dn):
        self._A_dn = A_dn
   
    @property
    def tray_type(self):
        """Default 'Sieve'"""
        return self._tray_type
   

    @tray_type.setter
    def tray_type(self, tray_type):
        if tray_type in distillation_tray_type_factor:
            self._tray_type = tray_type
            F_D = self.F_D
            F_D['Trays'] = F_D['Stripper trays']  = distillation_tray_type_factor[tray_type]
        else:
            raise ValueError("tray type must be one of the following: "
                            f"{', '.join(distillation_tray_type_factor)}")
         
    @property
    def tray_material(self):
        """Default 'Carbon steel'"""
        return self._tray_material
    @tray_material.setter
    def tray_material(self, tray_material):
        if tray_material in tray_material_factor_functions:
            self._tray_material = tray_material
            self._F_TM_function = tray_material_factor_functions[tray_material]
        else:
            raise ValueError("tray material must be one of the following: "
                            f"{', '.join(tray_material_factor_functions)}")
         
    @property
    def vessel_material(self):
        """Default 'Carbon steel'"""
        return self._vessel_material
    @vessel_material.setter
    def vessel_material(self, vessel_material):
        if vessel_material in distillation_column_material_factors:
            self._vessel_material = vessel_material
            ### revisar
            F_M = self.F_M
            F_M['Stripper tower'] = F_M['Tower'] = distillation_column_material_factors[vessel_material]            
        else:
            raise ValueError("vessel material must be one of the following: "
                            f"{', '.join(distillation_column_material_factors)}")
   
    def _actual_stages(self):
        """Return a tuple with the actual number of stages for the rectifier and the stripper."""
        eff = self.stage_efficiency
        if eff is None:
            # Calculate Murphree Efficiency
            vapor, liquid = self.outs
            mu = liquid.get_property('mu', 'mPa*s')
            alpha = self._get_relative_volatilities()
            L_Rmol = liquid.F_mol
            V_Rmol = vapor.F_mol
            eff = design.compute_murphree_stage_efficiency(
                mu, alpha, L_Rmol, V_Rmol
            )
            
        # Calculate actual number of stages
        return np.ceil(self.N_stages / eff)
       
    def _get_relative_volatilities(self):
        stages = self.stages
        IDs = stages[0].partition.IDs
        solute = self.solute
        index = IDs.index(solute)
        K = gmean([i.partition.K[index] for i in stages])
        if self.liquid.imol[solute] > self.vapor.imol[solute]:
            self.line = "Stripping"
            alpha = 1 / K
        else:
            self.line = "Absorption"
            alpha = K
        return alpha
       
    def _design(self):
        vapor_out, liquid_out = self.outs
        Design = self.design_results
        TS = self._TS
        A_ha = self._A_ha
        F_F = self._F_F
        f = self._f
       
        ### Get diameter of column based on outlets (assuming they are comparable to each stages) ###
        rho_L = liquid_out.rho
        V = vapor_out.F_mass
        V_vol = vapor_out.get_total_flow('m^3/s')
        rho_V = vapor_out.rho
        L = liquid_out.F_mass # To get liquid going down
        F_LV = design.compute_flow_parameter(L, V, rho_V, rho_L)
        C_sbf = design.compute_max_capacity_parameter(TS, F_LV)
        sigma = liquid_out.get_property('sigma', 'dyn/cm')
        U_f = design.compute_max_vapor_velocity(C_sbf, sigma, rho_L, rho_V, F_F, A_ha)
        A_dn = self._A_dn
        if A_dn is None:
            A_dn = design.compute_downcomer_area_fraction(F_LV)
        diameter = design.compute_tower_diameter(V_vol, U_f, f, A_dn) * 3.28
        Po = self.P * 0.000145078 # to psi
        rho_M = material_densities_lb_per_in3[self.vessel_material]
       
        if Po < 14.68:
            warn('vacuum pressure vessel ASME codes not implemented yet; '
                 'wall thickness may be inaccurate and stiffening rings may be '
                 'required', category=RuntimeWarning)
            
        Design['Theoretical stages'] = self.N_stages
        Design['Actual stages'] = actual_stages = self._actual_stages()
        Design['Height'] = H = design.compute_tower_height(TS, actual_stages) * 3.28
        Design['Diameter'] = Di = diameter
        Design['Wall thickness'] = tv = design.compute_tower_wall_thickness(Po, Di, H)
        Design['Weight'] = design.compute_tower_weight(Di, H, tv, rho_M)

    def _cost(self):
        Design = self.design_results
        Cost = self.baseline_purchase_costs
        F_M = self.F_M
     
        # Cost trays
        N_T = Design['Actual stages']
        Di = Design['Diameter']
        F_M['Trays'] = self._F_TM_function(Di)
        Cost['Trays'] = design.compute_purchase_cost_of_trays(N_T, Di)
           
        # Cost vessel assuming T < 800 F
        W = Design['Weight'] # in lb
        H = Design['Height'] # in ft
        Cost['Tower'] = design.compute_empty_tower_cost(W)
           
        Cost['Platform and ladders'] = design.compute_plaform_ladder_cost(Di, H)
       
Stripper = Absorber = AdiabaticMultiStageVLEColumn