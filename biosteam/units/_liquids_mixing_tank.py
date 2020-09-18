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
from biosteam.units.design_tools import geometry, PressureVessel
from biosteam.units.decorators import cost

__all__ = ('LiquidsMixingTank',)

@cost('Power', 'Turbine agitator', N='Number of agitators',
      ub=60, CE=567, cost=3730, n=0.54, BM=2.25)
class LiquidsMixingTank(bst.Unit, PressureVessel):
    """
    Create a LiquidsMixingTank for mixing two liquid phases.
    
    Parameters
    ----------
    ins : streams
        Inlet fluids to be mixed.
    outs : stream
        Mixed outlet fluid.
    tau=0.022 : float
        Residence time [hr].
    agitator_kW_per_m3=1.0 : float
        Electricity consumption in kW / m3 of volume.
    vessel_material='Carbon steel' : str, optional
        Vessel construction material.
    vessel_type='Horizontal': 'Horizontal' or 'Vertical', optional
        Vessel type.
    length_to_diameter=1 : float
        Length to diameter ratio.
        
    """
    _units = {**PressureVessel._units,
              'Volume': 'm^3',
              'Power': 'hp'}
    _ins_size_is_fixed = False
    _N_ins = 3
    _N_outs = 1
    
    def __init__(self, ID="", ins=None, outs=(), thermo=None, *,
                 tau=0.022, agitator_kW_per_m3=1.0, length_to_diameter=1,
                 vessel_material='Carbon steel',
                 vessel_type='Vertical'):
        bst.Unit.__init__(self, ID, ins, outs, thermo)
        self.length_to_diameter = length_to_diameter
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type
        self.agitator_kW_per_m3 = agitator_kW_per_m3
        self.tau = tau
        
    def _run(self):
        self.outs[0].mix_from(self.ins)
        
    def _design(self):
        results = self.design_results
        results['Volume'] = volume = self.tau * self.outs[0].F_vol
        P = self.ins[0].get_property('P', 'psi')
        length_to_diameter = self.length_to_diameter
        results = self.design_results
        rate = self.agitator_kW_per_m3 * volume
        self.power_utility(rate)
        results['Power'] = 1.341 * rate # in hp
        D = geometry.cylinder_diameter_from_volume(volume, length_to_diameter)
        L = length_to_diameter * D
        results.update(self._vessel_design(P, D, L))
        
    def _cost(self):
        self._decorated_cost()
        D = self.design_results
        self.purchase_costs.update(
            self._vessel_purchase_cost(D['Weight'], D['Diameter'], D['Length'])
        )