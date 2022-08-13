# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2022, Yoel Cortes-Pena <yoelcortes@gmail.com>, Ben Portner <github.com/BenPortner>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
.. contents:: :local:

.. autoclass:: biosteam.units._valve.Valve
.. autoclass:: biosteam.units._valve.IsenthalpicValve

"""
import biosteam as bst
from .. import Unit
from warnings import warn
from math import log, exp, ceil
from typing import NamedTuple, Tuple, Callable, Dict
from .heat_exchange import HX
from ..utils import list_available_names
from ..exceptions import DesignWarning, bounds_warning

__all__ = (
    'Valve',
    'IsenthalpicValve',
)

class ValveCostAlgorithm(NamedTuple):
    #: Defines preliminary correlation algorithm for a valve type
    psig_max: float #: Maximum achievable pressure in psig (to autodermine valve type and/or issue warning)
    hp_bounds: Tuple[float, float] #: Horse power per machine (not a hard limit for costing, but included here for completion)
    acfm_bounds: Tuple[float, float] #: Actual cubic feet per minute (hard limit for parallel units)
    cost: Callable #: function(horse_power) -> Baseline purchase cost
    CE: float #: Chemical engineering price cost index.


class Valve(Unit, isabstract=True):
    """
    Abstract class for valves that includes design and costing. Child classes
    should implement the `_run` method for mass and energy balances.
    
    """
    _N_ins = 1
    _N_outs = 1
    _N_heat_utilities = 1
    _units = {
        'Ideal power': 'kW',
        'Ideal duty': 'kJ/hr',
    }
    material_factors = {
        'Carbon steel': 1.0,
        'Stainless steel': 2.5,
        'Nickel alloy': 5.0,
    }
    _F_BM_default = {
        'Valve(s)': 2.15,
    }
    #: dict[str, ValveCostAlgorithm] Cost algorithms by valve type.
    baseline_cost_algorithms = { 
        'Default': ValveCostAlgorithm(
                psig_max=999999999.,
                acfm_bounds=(0., 999999999),
                hp_bounds=(0., 999999999.),
                cost=lambda Pc: 0,
                CE=567,
            ),
    }

    def __init__(self, ID='', ins=None, outs=(), thermo=None, *, 
                 P, vle=False, valve_type=None,
                 material=None):
        Unit.__init__(self, ID, ins, outs, thermo)
        self.P = P  #: Outlet pressure [Pa].

        #: Whether to perform phase equilibrium calculations on the outflow.
        #: If False, the outlet will be assumed to be the same phase as the inlet.
        self.vle = vle
        self.material = 'Carbon steel' if material is None else material 
        self.valve_type = 'Default' if valve_type is None else valve_type

    @property
    def valve_type(self):
        return self._valve_type
    @valve_type.setter
    def valve_type(self, valve_type):
        """[str] Type of valve. If 'Default', the type will be determined based on the outlet pressure."""
        valve_type = valve_type.capitalize()
        if valve_type not in self.baseline_cost_algorithms and valve_type != 'Default':
            raise ValueError(
                f"valve type {repr(valve_type)} not available; "
                f"only {list_available_names(self.baseline_cost_algorithms)} are available"
            )
        self._valve_type = valve_type

    @property
    def material(self):
        """[str]. Construction material. Defaults to 'Carbon steel'."""
        return self._material
    @material.setter
    def material(self, material):
        try:
            self.F_M['Valve(s)'] = self.material_factors[material]
        except KeyError:
            raise AttributeError("material must be one of the following: "
                                 f"{list_available_names(self.material_factors)}")
        self._material = material

    def _determine_valve_type(self):
        psig = (self.P - 101325.) * 14.6959 / 101325.
        cost_algorithms = self.baseline_cost_algorithms
        for name, alg in cost_algorithms.items():
            if psig < alg.psig_max: return name
        warn('no valve available that is recommended for a pressure of '
            f'{self.P:.5g}; defaulting to {name.lower()} valve', DesignWarning)
        return name
    
    def _design(self):
        design_results = self.design_results
        valve_type = self.valve_type
        if valve_type == 'Default': valve_type = self._determine_valve_type()
        design_results['Type'] = valve_type
        alg = self.baseline_cost_algorithms[valve_type]
        acfm_lb, acfm_ub = alg.acfm_bounds
        acfm = self.ins[0].get_total_flow('cfm')
        design_results['Valves in parallel'] = ceil(acfm / acfm_ub) if acfm > acfm_ub else 1
    
    def _cost(self):
        # Note: Must run `_set_power` before running parent cost algorithm
        design_results = self.design_results
        alg = self.baseline_cost_algorithms[design_results['Type']]
        acfm_lb, acfm_ub = alg.acfm_bounds
        Pc = self.power_utility.get_property('consumption', 'hp')
        N = design_results['Valves in parallel']
        Pc_per_valve = Pc / N
        self.baseline_purchase_costs['Valve(s)'] = N * bst.CE / alg.CE * alg.cost(Pc_per_valve)


class IsenthalpicValve(Valve):
    """
    Create an isenthalpic valve. Reduces the pressure of a fluid while keeping the enthalpy
    constant (adiabatic flash).

    Parameters
    ----------
    ins : stream
        Inlet fluid.
    outs : stream
        Outlet fluid.
    P : float
        Outlet pressure [Pa].
    vle : bool
        Whether to perform phase equilibrium calculations on
        the outflow. If False, the outlet will be assumed to be the same
        phase as the inlet.
    type: str
        Type of valve. If None, the type
        will be determined automatically.

    """
    _N_heat_utilities = 0

    def _run(self):
        feed = self.ins[0]
        out = self.outs[0]
        out.copy_like(feed)
        if self.vle is True:
            out.vle(H=feed.H, P=self.P)
        else:
            out.P = self.P
            out.H = feed.H
        
    def _design(self):
        super()._design()
