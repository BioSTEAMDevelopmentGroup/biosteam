# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
__all__ = ('design', 'add_design')

from thermosteam.units_of_measure import stream_units_of_measure

__all__ = ('design',)

# %% Design Center class

def _design(self):
    D = self.design_results
    U = self._units
    for i, j in self._design_basis_: D[i] = j(self, U[i])

class DesignCenter:
    """
    Create a DesignCenter object that manages all design basis functions.
    When called, it returns a Unit class decorator that adds a design item to
    the given Unit class.
    """
    __slots__ = ('design_basis_functions',)
    
    def __init__(self):
        self.design_basis_functions = {}
    
    def define(self, design_basis):
        """
        Define a new design basis.
        
        Parameters
        ----------
        design_basis : function
            Should accept two arguments, unit object and the units of measure, and return design basis value.
    
        Notes
        -----
        Design basis is registered with the name of the design basis function
        capitalized with underscores replaced by spaces.
        
        """
        name = design_basis.__name__.replace('_', ' ').capitalize()
        functions = self.design_basis_functions
        if name in functions: raise ValueError(f"design basis '{name}' already implemented")
        functions[name] = design_basis
        return design_basis
    
    def __call__(self, name, units):    
        """
        Return a Unit class decorator that adds a size/design requirement to the class.
        
        Parameters
        ----------
        name : str
            Name of design item            
        units : str
            Units of measure of design item.            
        
        """
        return lambda cls: self.add_design_basis_to_cls(cls, name, units)
    
    def add_design_basis_to_cls(self, cls, name, units):
        """
        Add size/design requirement to class.
        
        Parameters
        ----------
        cls : Unit class.    
        name : str
            Name of design item.        
        units : str
            Units of measure of design item.
            
        Examples
        --------
        
        :doc:`Unit decorators`
        
        """
        f = self.design_basis_functions[name.capitalize()]
        
        # Make sure design basis is not defined
        if name in cls._units:
            raise RuntimeError(f"design basis '{name}' already defined in class")
        else:
            cls._units[name] = units
        
        # Add design basis
        if hasattr(cls, '_decorated_design'):
            cls._design_basis_.append((name, f))
        else:
            cls._design_basis_ = [(name, f)]
            cls._decorated_design = _design
        if '_design' not in cls.__dict__:
            cls._design = _design
        
        return cls

    def __contains__(self, basis):
        return basis in self.design_basis_functions
    
    def __iter__(self):
        yield from self.design_basis_functions
    
    def __repr__(self):
        return f"<{type(self).__name__}: {', '.join(self)}>"

# %% Design factories
  
design = DesignCenter() #: Used to decorate classes with new design item

@design.define
def flow_rate(self, units):
    if self._N_ins == 1:
        return self._ins[0].get_total_flow(units)
    elif self._N_outs == 1:
        return self._outs[0].get_total_flow(units)
    elif self._N_ins < self._N_outs: 
        return sum([i.get_total_flow(units) for i in self._ins])
    else:
        return sum([i.get_total_flow(units) for i in self._outs])

H_units = stream_units_of_measure['H']

@design.define
def duty(self, units):
    duty = self.H_out - self.H_in
    self.heat_utilities[0](duty, self.ins[0].T, self.outs[0].T)
    return abs(H_units.conversion_factor(units) * duty)
    
@design.define
def dry_flow_rate(self, units):
    ins = self._ins
    flow_in = sum([i.get_total_flow(units) for i in ins])
    moisture = sum([i.get_flow(units, '7732-18-5') for i in ins])
    return flow_in - moisture

del flow_rate, duty, dry_flow_rate



