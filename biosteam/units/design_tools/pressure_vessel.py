# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from .specification_factors import (pressure_vessel_material_factors,
                                    material_densities_lb_per_ft3)
from . import flash_vessel_design as design
from ...utils import bounds_warning, DesignWarning
from warnings import warn

__all__ = ('PressureVessel',)

allowed_vessel_types = {'Vertical', 'Horizontal', None}

class PressureVessel:
    """Abstract class for pressure vessels."""

    _units = {'Vertical vessel weight': 'lb',
              'Horizontal vessel weight': 'lb',
              'Length': 'ft',
              'Diameter': 'ft',
              'Weight': 'lb',
              'Wall thickness': 'in'}    

    _bounds = {'Vertical vessel weight': (4200, 1e6),
               'Horizontal vessel weight': (1e3, 9.2e5),
               'Horizontal vessel diameter': (3, 21),
               'Vertical vessel length': (12, 40)}
    
    # Bare module factors
    _F_BM_default = {'Horizontal pressure vessel': 3.05,
                     'Vertical pressure vessel': 4.16,
                     'Platform and ladders': 1.}
    
    @property
    def vessel_type(self):
        vessel_type = self._vessel_type
        if not vessel_type:
            self._vessel_type = vessel_type = self._default_vessel_type()
        return vessel_type
    @vessel_type.setter
    def vessel_type(self, vessel_type):
        if vessel_type not in allowed_vessel_types:
            raise ValueError("vessel type must be either 'Vertical', "
                             "'Horizontal', or None")
        self._vessel_type = vessel_type
    
    @property
    def vessel_material(self):
        """Vessel construction material."""
        return self._vessel_material
    @vessel_material.setter
    def vessel_material(self, material):
        try: 
            F_M = self.F_M
            F_M['Vertical pressure vessel'] = F_M['Horizontal pressure vessel'] = pressure_vessel_material_factors[material]
        except KeyError:
            raise ValueError(f"no material factor available for '{material}'; "
                              "only the following materials are available: "
                             f"{', '.join(pressure_vessel_material_factors)}")
        self._vessel_material = material  
    
    def _get_design_info(self):
        return (('Vessel material', self._vessel_material, ''),)
    
    def _default_vessel_type(self):
        return None
    
    def _vessel_design(self, pressure, diameter, length) -> dict:
        vessel_type = self.vessel_type
        if vessel_type == 'Horizontal':
            method = self._horizontal_vessel_design
        elif vessel_type == 'Vertical':
            method = self._vertical_vessel_design
        else:
            raise RuntimeError('unknown vessel type')
        return method(pressure, diameter, length)
    
    def _horizontal_vessel_design(self, pressure, diameter, length) -> dict:
        pressure = float(pressure)
        diameter = float(diameter)
        length = float(length)
        # Calculate vessel weight and wall thickness
        rho_M = material_densities_lb_per_ft3[self._vessel_material]
        if pressure < 14.68:
            warn('vacuum pressure vessel ASME codes not implemented yet; '
                 'wall thickness may be inaccurate and stiffening rings may be '
                 'required', category=DesignWarning)
        VW, VWT = design.compute_vessel_weight_and_wall_thickness(
            pressure, diameter, length, rho_M)
        bounds_warning(self, 'Horizontal vessel weight', VW, 'lb',
                       self._bounds['Horizontal vessel weight'], 'cost')
        bounds_warning(self, 'Horizontal vessel diameter', diameter, 'ft',
                       self._bounds['Horizontal vessel diameter'], 'cost')
        Design = {}
        Design['Vessel type'] = 'Horizontal'
        Design['Length'] = length  # ft
        Design['Diameter'] = diameter  # ft
        Design['Weight'] = VW  # lb
        Design['Wall thickness'] = VWT  # in
        return Design
    
    def _vertical_vessel_design(self, pressure, diameter, length) -> dict:
        pressure = float(pressure)
        diameter = float(diameter)
        length = float(length)
        rho_M = material_densities_lb_per_ft3[self._vessel_material]
        if pressure < 14.68:
            warn('vacuum pressure vessel ASME codes not implemented yet; '
                 'wall thickness may be inaccurate and stiffening rings may be '
                 'required', category=DesignWarning)
        VW, VWT = design.compute_vessel_weight_and_wall_thickness(
            pressure, diameter, length, rho_M)
        Design = {}
        bounds_warning(self, 'Vertical vessel weight', VW, 'lb',
                       self._bounds['Vertical vessel weight'],
                       'cost')
        bounds_warning(self, 'Vertical vessel length', length, 'ft',
                       self._bounds['Vertical vessel length'],
                       'cost')
        Design['Vessel type'] = 'Vertical'
        Design['Length'] = length # ft
        Design['Diameter'] = diameter # ft
        Design['Weight'] = VW # lb
        Design['Wall thickness'] = VWT  # in
        return Design

    def _vessel_purchase_cost(self, weight, diameter, length) -> dict:
        weight = float(weight)
        diameter = float(diameter)
        length = float(length)
        vessel_type = self.vessel_type
        if vessel_type == 'Horizontal':
            method = self._horizontal_vessel_purchase_cost
        elif vessel_type == 'Vertical':
            method = self._vertical_vessel_purchase_cost
        else:
            raise RuntimeError('unknown vessel type')
        return method(weight, diameter, length)

    def _horizontal_vessel_purchase_cost(self, weight, diameter, length=None) -> dict:
        return {'Horizontal pressure vessel': design.compute_horizontal_vessel_purchase_cost(float(weight)),
                'Platform and ladders': design.compute_horizontal_vessel_platform_and_ladders_purchase_cost(float(diameter))}

    def _vertical_vessel_purchase_cost(self, weight, diameter, length) -> dict:
        return {'Vertical pressure vessel': design.compute_vertical_vessel_purchase_cost(float(weight)),
                'Platform and ladders': design.compute_vertical_vessel_platform_and_ladders_purchase_cost(float(diameter), float(length))}