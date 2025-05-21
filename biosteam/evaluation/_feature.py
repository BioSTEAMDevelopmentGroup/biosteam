# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
__all__ = ('Feature', 'MockFeature', 'Variable', 'MockVariable')
from ._name import element_name
from ..utils import format_title
from thermosteam import units_of_measure

class Feature:
    """
    Abstract class for a feature in BioSTEAM.
    
    Attributes
    ----------
    name : str
        Name of feature.
    units : str
        Units of measure.
    element : object
        Element corresponding to feature.
        
    """
    __slots__ = ('name', 'units', 'element')
    include_units_in_index = True
    
    def __init__(self, name, units, element):
        self.name = format_title(name)
        self.units = units
        self.element = element
    
    def mockup(self):
        return MockFeature(self.name, self.units, self.element)
    
    @classmethod
    def check_index_unique(cls, feature, features):
        key = (feature.element_name, feature.name)
        keys = {(i.element_name, i.name) for i in features}
        if key in keys:
            raise ValueError(
                     "each feature must have a unique element and name; "
                    f"feature with element {repr(feature.element)} "
                    f"and name {repr(feature.name)} already present"
                )
    
    @classmethod
    def check_indices_unique(cls, features):
        keys = set()
        for i in features:
            key = (i.element, i.name)
            if key in keys:
                kind = cls.__name__.lower()
                raise ValueError(
                        f"each {kind} must have a unique element and name; "
                        f"more than one {kind} with element {repr(i.element)} "
                        f"and name {repr(i.name)} are present"
                    )
            keys.add(key)
    
    @property
    def element_name(self):
        return element_name(self.element)
    
    @property
    def name_with_units(self):
        units = self.units
        name = self.name
        if units: name += f" [{units}]"
        return name
    
    @property
    def index(self):
        name = self.name
        if self.include_units_in_index:
            units = self.units
            if units: name += f" [{units}]"
        return (self.element_name, name)
    
    @property
    def short_description(self):
        element, name = self.index
        name, *_ = name.split(' [')
        if element not in name:
            name = ' '.join([element, name])
        if len(name) > 31:
            words = name.split(' ')
            words = [(i[:4]+'.' if len(i) > 5 else i) for i in words]
            name = ' '.join(words)
        name = name.strip(' ')
        if len(name) > 31: name = name[:31]
        return name
    
    def label(self, element=True, element_dlim='\n', format_units=True):
        return self.describe(
            element_dlim, 
            format_units,
            False, 
            False, 
            element=element,
        )
    
    def describe(self, 
            element_dlim=' - ',
            format_units=False,
            distribution=True, 
            bounds=True, 
            number_format='.3g', 
            element=True,
        ) -> str:
        """Return description of feature."""
        name = self.name
        if element and self.element:
            name = self.element_name + element_dlim + name
        units = self.units
        if format_units and units:
            units = units_of_measure.format_units(units)
        description = name
        units_added = False
        if distribution:
            if getattr(self, 'distribution', None):
                dist_name = type(self.distribution).__name__
                distribution_values = self.distribution._repr.values()
                distribution = ', '.join([format(j, number_format)
                                          for j in distribution_values])
                if units:
                    description += f' ({dist_name}; {distribution} {str(units)})'
                else:
                    description += f' ({dist_name}; {distribution})'
                units_added = True
        elif bounds:
            baseline = getattr(self, 'baseline', None)
            bounds = getattr(self, 'bounds', None)
            if bounds and baseline:
                lb, ub = bounds
                values = ', '.join([format(i, number_format)
                                    for i in (lb, baseline, ub)])
                if units:
                    description += f' ({values} {str(units)})'
                else:
                    description += f' ({values})'
                units_added = True
        if description:
            first_letter = description[0]
            if first_letter.islower(): 
                description = first_letter.upper() + description[1:]
        if units and not units_added: 
            description += ' [' + str(units) + ']'
        return description
    
    def __repr__(self):
        units = f" ({self.units})" if self.units else ""
        element = f" [{self.element_name}]" if self.element else ""
        return f'<{type(self).__name__}:{element} {self.name}{units}>'
    
    def show(self):
        print(self._info())
     

class MockFeature(Feature):
    __slots__ = ()
    
    def __init__(self, name, units, element):
        super().__init__(name, units, element_name(element))
    
    def __repr__(self):
        return f"{type(self).__name__}('{self.name}', '{self.units}', '{self.element}')"

Variable = Feature
MockVariable = MockFeature 