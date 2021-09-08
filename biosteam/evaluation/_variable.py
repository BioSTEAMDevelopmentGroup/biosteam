# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
__all__ = ('Variable', 'MockVariable')
from ._name import element_name

class Variable:
    """
    Abstract class for a variable in BioSTEAM.
    
    Attributes
    ----------
    name : str
        Name of variable.
    units : str
        Units of measure.
    element : object
        Element corresponding to variable.
        
    """
    __slots__ = ('name', 'units', 'element')
    include_units_in_index = True
    
    def __init__(self, name, units, element):
        self.name = name
        self.units = units
        self.element = element
    
    def mockup(self):
        return MockVariable(self.name, self.units, self.element)
    
    @classmethod
    def check_index_unique(cls, variable, variables):
        key = (variable.element, variable.name)
        keys = {(i.element, i.name) for i in variables}
        if key in keys:
            kind = cls.__name__.lower()
            raise ValueError(
                    f"each {kind} must have a unique element and name; "
                    f"{kind} with element {repr(variable.element)} "
                    f"and name {repr(variable.name)} already present"
                )
    
    @classmethod
    def check_indices_unique(cls, variables):
        keys = set()
        for i in variables:
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
        name = ' '.join([element, name])
        if len(name) > 31:
            words = name.split(' ')
            words = [(i[:4]+'.' if len(i) > 5 else i) for i in words]
            name = ' '.join(words)
        name = name.strip(' ')
        return name
    
    def describe(self, number_format='.3g', distribution=True) -> str:
        """Return description of variable."""
        name = self.name
        if not name.isupper():
            name = name.casefold()
        if self.element:
            name = self.element_name + ' ' + name
        if self.units:
            units = (' [' + str(self.units) + ']')
        else:
            units = ''
        if distribution:
            if getattr(self, 'distribution', None):
                dist_name = type(self.distribution).__name__
                distribution_values = self.distribution._repr.values()
                distribution = ', '.join([format(j, number_format)
                                          for j in distribution_values])
                distribution = f' ({dist_name}; {distribution})'
            else:
                distribution = ''
            description = name + units + distribution
        else:
            baseline = getattr(self, 'baseline', None)
            bounds = getattr(self, 'bounds', None)
            if bounds and baseline:
                lb, ub = bounds
                values = ', '.join([format(i, number_format)
                                    for i in (lb, baseline, ub)])
                description = name + units + f' ({values})'
            else:
                description = name + units
        if description:
            first_letter = description[0]
            if first_letter.islower(): 
                description = first_letter.upper() + description[1:]
        return description
    
    
    def __repr__(self):
        units = f" ({self.units})" if self.units else ""
        element = f" [{self.element_name}]" if self.element else ""
        return f'<{type(self).__name__}:{element} {self.name}{units}>'
    
    def show(self):
        print(self._info())
     

class MockVariable(Variable):
    __slots__ = ()
    def __init__(self, name, units, element):
        self.name = name
        self.units = units
        self.element = element_name(element)
    
    @property
    def element_name(self):
        return self.element
    
    def __repr__(self):
        return f"MockVariable('{self.name}', '{self.units}', '{self.element}')"
