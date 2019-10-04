# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 06:15:40 2019

@author: yoelr
"""
__all__ = ('Variable',)
from ._name import element_name


class Variable:
    """Abstract class for a variable in BioSTEAM.
    
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
    
    @property
    def element_name(self):
        return element_name(self.element)
    
    @property
    def index(self):
        return (self.element_name, self.name)
    
    def describe(self, number_format='.3g') -> str:
        """Return description of variable."""
        name = self.name
        if not name.isupper():
            name = name.casefold()
        if self.element:
            name = self.element_name + ' ' + name
        if self.units:
            units = (' (' + str(self.units) + ')')
        else:
            units = ''
        try:
            distribution = ', '.join([format(j, number_format)
                                      for j in self.distribution._repr.values()])
            distribution = ' [' + distribution + ']'
        except:
            distribution = ''
        return name + units + distribution
    
    def __repr__(self):
        units = f" ({self.units})" if self.units else ""
        element = f" [{self.element_name}]" if self.element else ""
        return f'<{type(self).__name__}:{element} {self.name}{units}>'
    
    def show(self):
        print(self._info())