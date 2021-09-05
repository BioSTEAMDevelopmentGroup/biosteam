# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from .._unit import Unit
from .._graphics import process_specification_graphics
from ..utils import format_title, static

__all__ = ('ProcessSpecification',)

@static
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
    >>> main_flowsheet.clear() # Remove old units
    >>> main_flowsheet.set_flowsheet('mix_ethanol_with_denaturant')
    >>> settings.set_thermo(['Water', 'Ethanol', 'Octane'], cache=True)
    >>> ethanol = Stream('ethanol', T=340, Water=200, Ethanol=22500, units='kg/hr')
    >>> denaturant = Stream('denaturant', Octane=1)
    >>> def adjust_denaturant():
    ...     denaturant_over_ethanol_flow = 0.02 / 0.98 # A mass ratio
    ...     denaturant.imass['Octane'] = denaturant_over_ethanol_flow * ethanol.F_mass
    >>> PS1 = units.ProcessSpecification('PS1',
    ...                                  ins=ethanol, outs='ethanol_',
    ...                                  specification=adjust_denaturant)
    >>> M1 = units.Mixer('M1', ins=(PS1-0, denaturant), outs='denatured_ethanol')
    >>> system = main_flowsheet.create_system('mix_ethanol_with_denaturant_sys')
    >>> system.show()
    System: mix_ethanol_with_denaturant_sys
    ins...
    [0] denaturant
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Octane  1
    [1] ethanol
        phase: 'l', T: 340 K, P: 101325 Pa
        flow (kmol/hr): Water    11.1
                        Ethanol  488
    outs...
    [0] denatured_ethanol
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow: 0
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
    _N_heat_utilities = 0
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 specification, description=None):
        Unit.__init__(self, ID, ins, outs, thermo)
        self.specification = specification
        self.description = description or format_title(specification.__name__)
        
    @property
    def specification(self):
        return self.__specification
    @specification.setter
    def specification(self, specification):
        assert callable(specification), "specification must be a function"
        self.__specification = specification
    
    def run(self):
        self.outs[0].copy_flow(self.ins[0])
        self.__specification()
        
    def _run(self):
        self.outs[0].copy_flow(self.ins[0])
        self.__specification()
        

