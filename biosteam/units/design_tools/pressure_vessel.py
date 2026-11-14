# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from ..._unit import design_input
from .specification_factors import (pressure_vessel_material_factors,
                                    material_densities_lb_per_ft3)
from . import flash_vessel_design as design
from ...exceptions import bounds_warning, DesignWarning
from warnings import warn

__all__ = ('PressureVessel',)



class PressureVessel:
    """Abstract class for pressure vessels."""

    _units = {'Vessel weight': 'lb',
              'Length': 'ft',
              'Diameter': 'ft',
              'Jacketed diameter': 'ft',
              'Weight': 'lb',
              'Wall thickness': 'in',
              'Pressure': 'psi'}    

    _bounds = {'Vertical vessel weight': (4200, 1e6),
               'Horizontal vessel weight': (1e3, 9.2e5),
               'Horizontal vessel diameter': (3, 21),
               'Vertical vessel length': (12, 40)}
    
    # Bare module factors
    _F_BM_default = {
        'Platform and ladders': 1.
    }
    
    vessel_type = design_input(
        purchase_item='Pressure vessel',
        catalog_type='Module',
        catalog={'Horizontal': 3.05,
                 'Vertical': 4.16},
    )

    vessel_material = design_input(
        purchase_item='Pressure vessel',
        catalog_type='Material',
        catalog=pressure_vessel_material_factors,
    )

    def _default_vessel_type(self):
        raise NotImplementedError('must define _default_vessel_type method')
    
    def _design_and_cost_pressure_vessel(
            self, pressure, diameter, length,
            pressure_units=None, length_units=None, 
            n_vessels=None,
        ):
        # Convenience method for design and costing of a unit operation.
        if length_units is not None:
            diameter = self.convert_design_result('Diameter', length_units, diameter)
            length = self.convert_design_result('Length', length_units, length)
        if pressure_units is not None:
            pressure = self.convert_design_result('Pressure', pressure_units, pressure)
        design = self._vessel_design(pressure, diameter, length)
        self.design_results.update(design)
        purchase_costs = self._vessel_purchase_cost(design['Weight'], design['Diameter'], design['Length'])
        if n_vessels is not None:
            purchase_costs = {i + f' x{n_vessels}': j * n_vessels for i, j in purchase_costs.items()}
        self.baseline_purchase_costs.update(purchase_costs)
    
    def _vessel_design(self, pressure, diameter, length, annular_diameter=0, check=True) -> dict:
        vessel_type = self.vessel_type
        if vessel_type is None: 
            vessel_type = self._default_vessel_type()
        if vessel_type == 'Horizontal':
            method = self._horizontal_vessel_design
        elif vessel_type == 'Vertical':
            method = self._vertical_vessel_design
        else:
            raise RuntimeError('unknown vessel type')
        design = method(pressure, diameter, length, annular_diameter, check)
        if self.vessel_type is None: design['Vessel type'] = vessel_type
        return design
    
    def _horizontal_vessel_design(self, pressure, diameter, length, annular_diameter=0, check=True) -> dict:
        # Calculate vessel weight and wall thickness
        rho_M = material_densities_lb_per_ft3[self._vessel_material]
        VW, VWT = design.compute_vessel_weight_and_wall_thickness(
            pressure, diameter, length, rho_M, annular_diameter
        )
        if check:
            if pressure < 14.68:
                warning = DesignWarning.from_source(
                    self, 
                    'vacuum pressure vessel ASME codes not implemented yet; '
                    'wall thickness may be inaccurate and stiffening rings may be '
                    'required'
                )
                warn(warning)
            bounds_warning(self, 'Vessel weight', VW, 'lb',
                           self._bounds['Horizontal vessel weight'])
            bounds_warning(self, 'Vessel diameter', diameter, 'ft',
                           self._bounds['Horizontal vessel diameter'])
        Design = {}
        Design['Length'] = length  # ft
        Design['Diameter'] = diameter  # ft
        Design['Weight'] = VW  # lb
        Design['Wall thickness'] = VWT  # in
        if annular_diameter: Design['Jacketed diameter'] = annular_diameter + diameter + VWT / 12 # ft
        return Design
    
    def _vertical_vessel_design(self, pressure, diameter, length, annular_diameter=0, check=True) -> dict:
        rho_M = material_densities_lb_per_ft3[self._vessel_material]
        VW, VWT = design.compute_vessel_weight_and_wall_thickness(
            pressure, diameter, length, rho_M, annular_diameter
        )
        Design = {}
        if check:
            if pressure < 14.68:
                warning = DesignWarning.from_source(
                    self, 
                    'vacuum pressure vessel ASME codes not implemented yet; '
                    'wall thickness may be inaccurate and stiffening rings may be '
                    'required'
                )
                warn(warning)
            bounds_warning(self, 'Vessel weight', VW, 'lb',
                           self._bounds['Vertical vessel weight'])
            bounds_warning(self, 'Vessel length', length, 'ft',
                           self._bounds['Vertical vessel length'])
        Design['Length'] = length # ft
        Design['Diameter'] = diameter # ft
        Design['Weight'] = VW # lb
        Design['Wall thickness'] = VWT  # in
        if annular_diameter: Design['Jacketed diameter'] = annular_diameter + diameter + VWT / 12 # ft
        return Design

    def _vessel_purchase_cost(self, weight, diameter, length, jacketed=None) -> dict:
        vessel_type = self.vessel_type
        if vessel_type == 'Horizontal':
            method = self._horizontal_vessel_purchase_cost
        elif vessel_type == 'Vertical':
            method = self._vertical_vessel_purchase_cost
        else:
            raise RuntimeError('unknown vessel type')
        return method(weight, diameter, length, jacketed)

    def _horizontal_vessel_purchase_cost(self, weight, diameter, length=None, jacketed=None) -> dict:
        return {
            'Pressure vessel': design.compute_horizontal_vessel_purchase_cost(weight),
            'Platform and ladders': design.compute_horizontal_vessel_platform_and_ladders_purchase_cost(diameter)
        }
  
    def _vertical_vessel_purchase_cost(self, weight, diameter, length, jacketed=None) -> dict:
        return {
            'Pressure vessel': design.compute_vertical_vessel_purchase_cost(weight),
            'Platform and ladders': design.compute_vertical_vessel_platform_and_ladders_purchase_cost(diameter, length)
        }
        