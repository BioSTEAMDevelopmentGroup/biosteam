# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>, Ben Portner <github.com/BenPortner>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
.. contents:: :local:

.. autoclass:: biosteam.units.valve.Valve
.. autoclass:: biosteam.units.valve.IsenthalpicValve

"""
from .. import Unit
from thermosteam._graphics import valve_graphics
from warnings import warn

__all__ = (
    'Valve',
    'IsenthalpicValve',
)

class Valve(Unit, isabstract=True):
    """
    Abstract class for valves. Child classes should implement the `_run` 
    method for mass and energy balances. No design or costing algorithms 
    are implemented (yet). Any contributions to design and costing
    is certainly welcome. For now, it is assumed that the cost of valves in 
    a production process is negligible in relation to other unit operations. 
    Additionally, valves serve a level of detail that is above the accuracy 
    of cost correlations in BioSTEAM (which serve preliminary techno-economic 
    analysis purposes).
    
    """
    _graphics = valve_graphics

class IsenthalpicValve(Valve, new_graphics=False):
    """
    Create an isenthalpic valve. Reduces the pressure of a fluid while keeping the enthalpy
    constant (adiabatic flash). 

    Parameters
    ----------
    ins : 
        Inlet fluid.
    outs : 
        Outlet fluid.
    P : float
        Outlet pressure [Pa].
    vle : bool
        Whether to perform phase equilibrium calculations on
        the outflow. If False, the outlet will be assumed to be the same
        phase as the inlet.

    Warnings
    --------
    No design or costing algorithms have been implemented (yet). For now, it 
    is assumed that the cost of valves in a production process is negligible 
    in relation to other unit operations. Additionally, valves serve a level of 
    detail that is above the accuracy of cost correlations in BioSTEAM 
    (which serve preliminary techno-economic analysis purposes).

    """

    def _init(self, P, vle=False):
        self.P: float = P  #: Outlet pressure [Pa].

        #: Whether to perform phase equilibrium calculations on the outflow.
        #: If False, the outlet will be assumed to be the same phase as the inlet.
        self.vle: bool = vle

    def _run(self):
        feed = self.ins[0]
        out = self.outs[0]
        out.copy_like(feed)
        if feed.P > self.P:
            if self.vle is True:
                out.vle(H=feed.H, P=self.P)
            else:
                out.P = self.P
                out.H = feed.H
        else:
            warn(f"feed pressure ({feed.P:.5g} Pa) is lower or equal than outlet "
                 f"specification ({self.P:.5g} Pa); valve {self.ID} is ignored", RuntimeWarning)
