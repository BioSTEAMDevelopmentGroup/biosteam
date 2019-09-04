# -*- coding: utf-8 -*-
"""
Created on Tue May 14 14:20:53 2019

@author: yoelr
"""
__all__ = ('Parameter',)
from ._name import elementname

class Parameter:
    """Create a Parameter object that, when called, runs the setter and the simulate functions.
    
    Parameters
    ----------
    name : str
           Name of parameter.
    setter : function
             Should set the parameter.
    simulate : function
               Should simulate parameter effects.
    element : object
              Element associated to parameter.
    system : System
             System associated to parameter.
    distribution : chaospy.Dist
                   Parameter distribution.
    units : str
            Units of parameter.
    
    """
    __slots__ = ('name', 'setter', 'simulate', 'element',
                 'system', 'distribution', 'units')
    
    def __init__(self, name, setter, simulate,
                 element, system, distribution,
                 units):
        self.name = name.replace('_', ' ').capitalize()
        self.setter = setter
        self.simulate = simulate
        self.element = element
        self.system = system
        self.distribution = distribution
        self.units = units
    
    def __call__(self, value):
        self.setter(value)
        self.simulate()
    
    @property
    def element_name(self):
        return elementname(self.element)
    
    def __repr__(self):
        units = f" ({self.units})" if self.units else ""
        element = f" [{self.element_name}]" if self.element else ""
        return f'<{type(self).__name__}:{element} {self.name}{units}>'
    
    def describe(self, number_format='.3g') -> str:
        """Return description of parameter."""
        if self.element:
            name = self.element_name + ' ' + self.name.casefold()
        else:
            name = self.name
        if self.units:
            units = (' [' + str(self.units) + ']')
        else:
            units = ''
        if self.distribution:
            distribution = ', '.join([format(j, number_format) for j in self.distribution._repr.values()])
            distribution = ' (' + distribution + ')'
        else:
            distribution = ''
        return name + units + distribution
    
    def show(self):
        print(f'{type(self).__name__}: {self.name}\n'
           + (f' element: {self.element_name}' if self.element else '')
           + (f' system: {self.system}\n' if self.system else '')
           + (f' units: {self.units}\n' if self.units else '')
           + (f' distribution: {self.distribution}\n' if self.distribution else ''))
                
    
               