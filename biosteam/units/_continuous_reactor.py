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
import biosteam as bst

__all__ = ('ContinuousReactor', 'CSTR')


class ContinuousReactor(PressureVessel, Unit, isabstract=True):
    '''    
    Abstract class for reactor unit, modeled as a pressure vessel with 
    a given aspect ratio and residence time. A heat exchanger recirculation 
    loop is used to satisfy the duty, if any. A vacuum system is also 
    automatically added if the operating pressure is at a vacuum.

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
    auxiliary_unit_names = ('heat_exchanger', 'vacuum_system')
    
    _units = {**PressureVessel._units,
              'Residence time': 'hr',
              'Total volume': 'm3',
              'Reactor volume': 'm3',
              'Single reactor volume': 'm3'}
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *, 
                 T=None, P=None, tau=0.5, V_wf=0.8, V_max=355,
                 length_to_diameter=2, kW_per_m3=0.985,
                 vessel_material='Stainless steel 316',
                 vessel_type='Vertical'):
        
        Unit.__init__(self, ID, ins, outs, thermo)
        self.T = T
        self.P = P
        self.tau = tau
        self.V_wf = V_wf
        self.V_max = V_max or self.V_max_default
        self.length_to_diameter = length_to_diameter
        self.kW_per_m3 = kW_per_m3
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type
        self.heat_exchanger = bst.HXutility(None, (None,), (None,), thermo=self.thermo) 

    def _design(self):
        Design = self.design_results
        ins_F_vol = self.F_vol_in
        V_total = ins_F_vol * self.tau / self.V_wf
        P_pascal = (self.P if self.P else self.outs[0].P)
        P_psi = P_pascal * 0.000145038 # Pa to psi
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
        Design.update(self._vessel_design(float(P_psi), float(D), float(L)))
        self.vacuum_system = bst.VacuumSystem(self) if P_pascal < 1e5 else None
        duty = self.Hnet
        self.parallel['self'] = N
        self.parallel['vacuum_system'] = 1 # Not in parallel
        if duty: 
            # Note: Flows and duty are multiplied by scale to simulate an individual
            # heat exchanger, then BioSTEAM accounts for number of units in parallel
            # through the `parallel` attribute.
            self.heat_exchanger.simulate_as_auxiliary_exchanger(
                self.ins, self.outs, duty, vle=False, scale=1 / N,
            )
            
    def _cost(self):
        Design = self.design_results
        baseline_purchase_costs = self.baseline_purchase_costs
        volume = Design['Single reactor volume']
        if volume != 0:
            baseline_purchase_costs.update(
                self._vessel_purchase_cost(
                    Design['Weight'], Design['Diameter'], Design['Length']
                )
            )
            self.add_power_utility(self.kW_per_m3 * volume)
    
CSTR = ContinuousReactor