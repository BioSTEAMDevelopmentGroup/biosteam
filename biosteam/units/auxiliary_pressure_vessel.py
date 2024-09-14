# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import biosteam as bst
from .auxiliary import Auxiliary
from . import design_tools as design
from .design_tools.specification_factors import (
    pressure_vessel_material_factors,
    material_densities_lb_per_ft3
)
from ..exceptions import bounds_warning, DesignWarning
from warnings import warn

__all__ = ('AuxiliaryPressureVessel',)

class AuxiliaryPressureVessel(Auxiliary):
    __slots__ = ()
    _bounds = design.PressureVessel._bounds
    _units = design.PressureVessel._units
    
    def __init__(self, 
            pressure, diameter, length,
            pressure_units=None, length_units=None, 
            orientation=None, material=None, 
        ):
        super().__init__()
        design_results = self.design_results
        if orientation is None: orientation = 'Vertical'
        if material is None: material = 'Carbon steel'
        design_results['Orientation'] = orientation
        design_results['material'] = material
        if length_units is not None:
            diameter = self.set_design_result('Diameter', length_units, diameter)
            length = self.set_design_result('Length', length_units, length)
        if pressure_units is not None:
            pressure = self.set_design_result('Pressure', pressure_units, pressure)
        rho_M = material_densities_lb_per_ft3[material]
        if pressure < 14.68:
            warn('vacuum pressure vessel ASME codes not implemented yet; '
                 'wall thickness may be inaccurate and stiffening rings may be '
                 'required', category=DesignWarning)
        weight, thickness = design.compute_vessel_weight_and_wall_thickness(
            pressure, diameter, length, rho_M
        )
        self.F_M['Platform and ladders'] = 1.
        design_results['Weight'] = weight  # lb
        design_results['Wall thickness'] = thickness # in
        baseline_purchase_costs = self.baseline_purchase_costs
        if orientation == 'Horizontal':
            key = 'Vertical pressure vessel'
            bounds_warning(self, 'Horizontal vessel weight', weight, 'lb',
                           self._bounds['Horizontal vessel weight'], 'cost')
            bounds_warning(self, 'Horizontal vessel diameter', diameter, 'ft',
                           self._bounds['Horizontal vessel diameter'], 'cost')
            self.F_BM[key] = 3.05
            self.F_M[key] = pressure_vessel_material_factors[material]
            baseline_purchase_costs[key] = design.compute_horizontal_vessel_purchase_cost(weight)
            baseline_purchase_costs['Platform and ladders'] = design.compute_horizontal_vessel_platform_and_ladders_purchase_cost(diameter)
        elif orientation == 'Vertical':
            key = 'Vertical pressure vessel'
            bounds_warning(self, 'Vertical vessel weight', weight, 'lb',
                           self._bounds['Vertical vessel weight'], 'cost')
            bounds_warning(self, 'Vertical vessel length', length, 'ft',
                           self._bounds['Vertical vessel length'], 'cost')
            self.F_BM[key] = 4.16
            self.F_M[key] = pressure_vessel_material_factors[material]
            baseline_purchase_costs[key] = design.compute_vertical_vessel_purchase_cost(weight)
            baseline_purchase_costs['Platform and ladders'] = design.compute_vertical_vessel_platform_and_ladders_purchase_cost(diameter, length)
        else:
            raise ValueError(
                f"invalid orientation {orientation!r}; orientation must be "
                 "either 'Horizontal' or 'Vertical'"
            )
        