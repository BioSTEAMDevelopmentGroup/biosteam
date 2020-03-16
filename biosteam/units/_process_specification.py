# -*- coding: utf-8 -*-
"""
Created on Sat Jul 13 02:24:35 2019

@author: yoelr
"""
from .._unit import Unit
from .._graphics import process_specification_graphics
from ..utils import format_unit_line

__all__ = ('ProcessSpecification',)

class ProcessSpecification(Unit):
    """
    Create a ProcessSpecification object that runs a function when simulated
    to set a process specification.
    
    Parameters
    ----------
    ins : stream, optional
        Inlet stream.
    outs : stream, optional
        Outlet stream.
    specification : callable
        Called during simulation to set a process specification.
    description : 
        Description of process specification. Defaults to name of the `specification` function.
    
    Examples
    --------
    Use a ProcessSpecifcation object to define the amount of denaturant to
    add according to the flow of bioethanol. The final bioethanol product
    must be 2 wt. % denaturant:
    
    >>> from biosteam import settings, Stream, units, main_flowsheet
    >>> main_flowsheet.set_flowsheet('mix_ethanol_with_denaturant')
    >>> settings.set_thermo(['Water', 'Ethanol', 'Octane'])
    >>> ethanol = Stream('ethanol', T=340, Water=200, Ethanol=22500, units='kg/hr')
    >>> denaturant = Stream('denaturant', Octane=1)
    >>> def adjust_denaturant():
    ...     denaturant_over_ethanol_flow = 0.02 / 0.98 # A mass ratio
    ...     denaturant.imass['Octane'] = denaturant_over_ethanol_flow * ethanol.F_mass
    >>> PS1 = units.ProcessSpecification(adjust_denaturant,
    ...                                  ins=ethanol, outs='ethanol_')
    >>> M1 = units.Mixer('M1', ins=(PS1-0, denaturant), outs='denatured_ethanol')
    >>> system = main_flowsheet.create_system()
    >>> system.show()
    System: mix_ethanol_with_denaturant_sys
     path: (adjust_denaturant, M1)
    >>> system.simulate()
    >>> M1.outs[0].show(composition=True, flow='kg/hr')
    Stream: denatured_ethanol from <Mixer: M1>
     phase: 'l', T: 339.32 K, P: 101325 Pa
     composition: Water    0.00863
                  Ethanol  0.971
                  Octane   0.02
                  -------  2.32e+04 kg/hr
    
    """
    _graphics = process_specification_graphics
    _N_ins = _N_outs = 1
    power_utility = None
    results = None
    heat_utilities = ()
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 specification, description=None):
        self._numerical_specification = None
        self._load_thermo(thermo)
        self._init_ins(ins)
        self._init_outs(outs)
        self._assert_compatible_property_package()
        self._register(ID)
        self.specification = specification
        self.description = format_unit_line(description or specification.__name__)
        
    @property
    def specification(self):
        return self._specification
    @specification.setter
    def specification(self, specification):
        assert callable(specification), "specification must be a function"
        self._run = self._specification = specification
        

