# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>, Brenda Cansino <cansinoloeza@wisc.edu>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""

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



__all__ = ('Stripping',)

class Stripping (MultiStageEquilibrium):
    r"""
    Stripping column class.  
    The diameter is based on tray separation
    and flooding velocity. Purchase costs are based on correlations
    compiled by Warren et. al.

    Parameters
    ----------
    ins :
        Inlet fluids
        * [0] Molar flow rate of solute-free gas (gas carrier)
        * [1]  Molar flow rate of solute-free absorbent (liquid)      
    outs :
        * [0] Molar flow rate of gas stream
        * [1] Molar flow rate of liquid stream

    y : float
        Mole ratio of solute to solute-free gas in the vapor.
    x : float
        Mole ratio of solute to solute-ree absorbent in the liquid.

    P : float
        Operating pressure [Pa].
    T : float
        Operating Temperature.
       
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


    """
    line = 'Stripping'
   
    _graphics = vertical_column_graphics
    _ins_size_is_fixed = False
    _N_ins = 2
    _N_outs = 2
    _units = {
              'Stripper height': 'ft',
              'Stripper diameter': 'ft',
              'Stripper wall thickness': 'in',
              'Stripper weight': 'lb',
              'Height': 'ft',
              'Diameter': 'ft',
              'Wall thickness': 'in',
              'Weight': 'lb'}
    _max_agile_design = (
        'Actual stages',      
        'Stripper height',
        'Stripper diameter',
        'Stripper wall thickness',
        'Stripper weight',      
        'Height',
        'Diameter',
        'Wall thickness',
        'Weight',
    )
    _F_BM_default = { 'Stripper trays': 4.3,
                      'Platform and ladders': 1.,                    
                      'Stripper platform and ladders': 1.,
                      'Tower': 4.3,
                      'Trays': 4.3,
                      }
   
    # [dict] Bounds for results
    _bounds = {'Diameter': (3., 24.),
               'Height': (27., 170.),
               'Weight': (9000., 2.5e6)}
   
    def _init(self,
            ins,
            x ,
            y ,        
            P = None,
            T = None,
           
            product_specification_format=None,
            vessel_material='Carbon steel',
            tray_material='Carbon steel',
            tray_type='Sieve',
            tray_spacing=450,
            stage_efficiency=None,
            velocity_fraction=0.8,
            foaming_factor=1.0,
            open_tray_area_fraction=0.1,
            downcomer_area_fraction=None,  
        ):

       
        # Operation specifications
        self.x = x
        self.y = y
        self.L_feed = ins[1]
        self.P = P
        self.T = T


       
        # Construction specifications
        self.vessel_material = vessel_material
        self.tray_type = tray_type
        self.tray_material = tray_material
        self.tray_spacing = tray_spacing
        self.stage_efficiency = stage_efficiency
        self.velocity_fraction = velocity_fraction
        self.foaming_factor = foaming_factor
        self.open_tray_area_fraction = open_tray_area_fraction
        self.downcomer_area_fraction = downcomer_area_fraction
       
       

   
    @property
    def vapor_out(self):
        return self.outs[0]
    @property
    def liquid_out(self):
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
   
   
    @property
    def V_feed(self):
        return self.ins[0]

       
    def _stripping_column_design(self):
        vapor_out, liquid_out = self.outs
        Design = self.design_results
        Sstages = self._compute_N_stages()
        TS = self._TS
        A_ha = self._A_ha
        F_F = self._F_F
        f = self._f
        V_feed = self._V_feed
       
        ### Get diameter of stripping column based on feed plate ###
        rho_L = liquid_out.rho
   
        V = V_feed.F_mass
        V_vol = V_feed.get_total_flow('m^3/s')
        rho_V = V_feed.rho
        L = liquid_out.F_mass # To get liquid going down
        F_LV = design.compute_flow_parameter(L, V, rho_V, rho_L)
        C_sbf = design.compute_max_capacity_parameter(TS, F_LV)
       
       
        sigma = L.get_property('sigma', 'dyn/cm')
        U_f = design.compute_max_vapor_velocity(C_sbf, sigma, rho_L, rho_V, F_F, A_ha)
        A_dn = self._A_dn
        if A_dn is None:
            A_dn = design.compute_downcomer_area_fraction(F_LV)
        S_diameter = design.compute_tower_diameter(V_vol, U_f, f, A_dn) * 3.28
        Po = self.P * 0.000145078 # to psi
        rho_M = material_densities_lb_per_in3[self.vessel_material]
       
        if Po < 14.68:
            warn('vacuum pressure vessel ASME codes not implemented yet; '
                 'wall thickness may be inaccurate and stiffening rings may be '
                 'required', category=RuntimeWarning)

        else:
            Design['Actual stages'] =  Sstages
            Design['Height'] = H = design.compute_tower_height(TS, Sstages-2) * 3.28
            Design['Diameter'] = Di = max(( S_diameter))
            Design['Wall thickness'] = tv = design.compute_tower_wall_thickness(Po, Di, H)
            Design['Weight'] = design.compute_tower_weight(Di, H, tv, rho_M)
        self._simulate_components()
   


    def _cost(self):
        Design = self.design_results
        Cost = self.baseline_purchase_costs
        F_M = self.F_M
     
        # Cost trays assuming a partial condenser
        N_T = Design['Actual stages'] - 1.
        Di = Design['Diameter']
        F_M['Trays'] = self._F_TM_function(Di)
        Cost['Trays'] = design.compute_purchase_cost_of_trays(N_T, Di)
           
        # Cost vessel assuming T < 800 F
        W = Design['Weight'] # in lb
        H = Design['Height'] # in ft
        Cost['Tower'] = design.compute_empty_tower_cost(W)
           
        Cost['Platform and ladders'] = design.compute_plaform_ladder_cost(Di, H)
       
       
    def _design(self):
        self._stripping_column_design