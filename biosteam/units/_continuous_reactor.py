# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from .. import Unit
from math import ceil
from biosteam.units.design_tools import PressureVessel
from biosteam.units.design_tools.geometry import cylinder_diameter_from_volume

__all__ = ('ContinuousReactor',)


class ContinuousReactor(Unit, PressureVessel, isabstract=True):
    '''    
    Abstract class for reactor unit, modeled as a pressure vessel with 
    a given aspect ratio and residence time.

    Parameters
    ----------
    ins : stream
        Inlet.        
    outs : stream
        Outlet.
    tau=0.5 : float
        Residence time [hr].        
    V_wf=0.8 : float
        Fraction of working volume over total volume.   
    V_max=355 : float
        Maximum volume of a reactor in ft3.
    kW_per_m3=0.985: float
        Power usage of agitator
        (0.985 converted from 5 hp/1000 gal as in [1], for liquid–liquid reaction or extraction).
    vessel_material : str, optional
        Vessel material. Default to 'Stainless steel 316'.
    vessel_type : str, optional
        Vessel type. Can only be 'Horizontal' or 'Vertical'.
        
    References
    ----------
    .. [1] Seider, W. D.; Lewin, D. R.; Seader, J. D.; Widagdo, S.; Gani, R.; 
        Ng, M. K. Cost Accounting and Capital Cost Estimation. In Product 
        and Process Design Principles; Wiley, 2017; pp 470.
    '''
    _N_ins = 1
    _N_outs = 1
    _ins_size_is_fixed = False
    _outs_size_is_fixed = True
    
    _units = {**PressureVessel._units,
              'Residence time': 'hr',
              'Total volume': 'm3',
              'Reactor volume': 'm3',
              'Single reactor volume': 'm3'}
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *, 
                  P=101325, tau=0.5, V_wf=0.8, V_max=355,
                  length_to_diameter=2, kW_per_m3=0.985,
                  vessel_material='Stainless steel 316',
                  vessel_type='Vertical'):
        
        Unit.__init__(self, ID, ins, outs, thermo)
        self.P = P
        self.tau = tau
        self.V_wf = V_wf
        self.V_max = V_max or self.V_max_default
        self.length_to_diameter = length_to_diameter
        self.kW_per_m3 = kW_per_m3
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type

    def _design(self):
        Design = self.design_results
        ins_F_vol = self.F_vol_in
        V_total = ins_F_vol * self.tau / self.V_wf
        P = self.P * 0.000145038 # Pa to psi
        length_to_diameter = self.length_to_diameter
        N = ceil(V_total/self.V_max)
        if N == 0:
            V_reactor = 0
            D = 0
            L = 0
        else:
            V_reactor = V_total / N
            D = cylinder_diameter_from_volume(V_reactor, self.length_to_diameter)
            D *= 3.28084 # convert from m to ft
            L = D * length_to_diameter
        Design['Residence time'] = self.tau
        Design['Total volume'] = V_total
        Design['Single reactor volume'] = V_reactor
        Design['Number of reactors'] = N
        Design.update(self._vessel_design(float(P), float(D), float(L)))
            
    def _cost(self):
        Design = self.design_results
        baseline_purchase_costs = self.baseline_purchase_costs
        if Design['Total volume'] == 0:
            for i, j in baseline_purchase_costs.items():
                baseline_purchase_costs[i] = 0
            self.power_utility(0)
        else:
            baseline_purchase_costs.update(self._vessel_purchase_cost(
                Design['Weight'], Design['Diameter'], Design['Length']))
            for i, j in baseline_purchase_costs.items():
                baseline_purchase_costs[i] *= Design['Number of reactors']
            self.power_utility(self.kW_per_m3 * Design['Total volume'])
    