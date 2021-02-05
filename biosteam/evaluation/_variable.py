# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
__all__ = ('Variable',)
from ._name import element_name


class Variable:
    """
    Abstract class for a variable in BioSTEAM.
    
    Attributes
    ----------
    name : str
        Name of variable
    element : object
        Element corresponding to variable
    units : str
        Units of measure
    distribution : chaospy.Dist, optional
        Distribution of variable
        
    """
    __slots__ = ()
    include_units_in_index = True
    
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
    
    def describe(self, number_format='.3g') -> str:
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
        if self.distribution:
            distribution_values = self.distribution._repr.values()
            distribution = ', '.join([format(j, number_format)
                                      for j in distribution_values])
            distribution = ' (' + distribution + ')'
        else:
            distribution = ''
        description = name + units + distribution
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