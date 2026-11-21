# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>, Ben Portner <github.com/BenPortner>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
.. contents:: :local:

.. autoclass:: biosteam.units._turbine.Turbine
.. autoclass:: biosteam.units._turbine.IsentropicTurbine

"""

import biosteam as bst
from warnings import warn
from math import log, exp, ceil
from typing import NamedTuple, Tuple, Callable, Dict
from .heat_exchange import HX
from .. import Unit
from ..exceptions import DesignWarning, bounds_warning
from .._graphics import turbine_graphics


__all__ = (
    'Turbine',
    'IsentropicTurbine',
)

#: TODO:
#: * Implement estimate of isentropic efficiency when not given (is this possible?).
class TurbineCostAlgorithm(NamedTuple):
    #: Defines preliminary correlation algorithm for a turbine type
    psig_max: float #: Maximum achievable pressure in psig (to autodermine turbine type and/or issue warning)
    hp_bounds: Tuple[float, float] #: Horse power per machine (not a hard limit for costing, but included here for completion)
    acfm_bounds: Tuple[float, float] #: Actual cubic feet per minute (hard limit for parallel units)
    cost: Callable #: function(horse_power) -> Baseline purchase cost
    CE: float #: Chemical engineering price cost index.


class Turbine(Unit, isabstract=True):
    """
    Abstract class for turbines that includes design and costing. Child classes
    should implement the `_run` method for mass and energy balances. Preliminary 
    design and costing is estimated according to [1]_.
    
    """
    _graphics = turbine_graphics
    _N_ins = 1
    _N_outs = 1
    
    _F_BM_default = {
        'Turbine(s)': 2.15,
    }
    #: dict[str, TurbineCostAlgorithm] Cost algorithms by turbine type.
    baseline_cost_algorithms = {
        'Default': TurbineCostAlgorithm(
                psig_max=1e12,
                acfm_bounds=(0., 1e12),
                hp_bounds=(0., 1e12),
                cost=lambda Pc: 0,
                CE=567,
            ),
    }
    material_factors = {
        'Carbon steel': 1.0,
        'Stainless steel': 2.5,
        'Nickel alloy': 5.0
    }

    P = bst.design_input('Outlet pressure', 'Pa', lb=0)
    eta = bst.design_input('Isentropic efficiency', lb=0, ub=1, default=0.3)
    turbine_type = bst.design_input(
        catalog=baseline_cost_algorithms, default='Default',
        doc="Type of turbine. If 'Default', the type will be determined "
            "based on the outlet pressure.",
    )
    efficiency = bst.design_input(
        'Driver efficiency', lb=0, ub=1, default=0.65,
        doc="Efficiency of driver of work to electricity conversion. "
            "Defaults to 0.65.",
    )
    material = bst.design_input(
        catalog=material_factors,
        purchase_item='Turbines(s)',
        catalog_type='Material',
        default='Carbon steel',
    )

    def _init(self, vle=False):
        #: Whether to perform phase equilibrium calculations on the outflow.
        #: If False, the outlet will be assumed to be the same phase as the inlet.
        self.vle = vle

    def _determine_turbine_type(self):
        psig = (self.P - 101325.) * 14.6959 / 101325.
        cost_algorithms = self.baseline_cost_algorithms
        for name, alg in cost_algorithms.items():
            if psig < alg.psig_max: return name
        warn('no turbine available that is recommended for a pressure of '
            f'{self.P:.5g}; defaulting to {name.lower()} turbine', DesignWarning)
        return name

    def _calculate_ideal_power_and_duty(self):
        feed = self.ins[0]
        out = self.outs[0]
        if feed.P > out.P: raise RuntimeError('inlet pressure is above outlet')
        dH = out.H - feed.H
        Q = TdS = feed.T * (out.S - feed.S)  # Duty [kJ/hr]
        power_ideal = (dH - TdS) / 3600.  # Power [kW]
        return power_ideal, Q

    def _set_power(self, power):
        self.design_results['Motor efficiency'] = efficiency = self._efficiency
        self.power_utility(power * efficiency)
    
    def _design(self):
        if self.P > self.feed.P:
            self.power_utility(0.)
            return
        design_results = self.design_results
        turbine_type = self.turbine_type
        if turbine_type == 'Default': turbine_type = self._determine_turbine_type()
        design_results['Type'] = turbine_type
        alg = self.baseline_cost_algorithms[turbine_type]
        acfm_lb, acfm_ub = alg.acfm_bounds
        acfm = self.ins[0].get_total_flow('cfm')
        design_results['Turbines in parallel'] = ceil(acfm / acfm_ub) if acfm > acfm_ub else 1
    
    def _cost(self):
        if self.P > self.feed.P: return
        # Note: Must run `_set_power` before running parent cost algorithm
        design_results = self.design_results
        alg = self.baseline_cost_algorithms[design_results['Type']]
        acfm_lb, acfm_ub = alg.acfm_bounds
        Pc = abs(self.power_utility.get_property('rate', 'hp'))
        N = design_results['Turbines in parallel']
        Pc_per_turbine = Pc / N
        bounds_warning(self, 'power', Pc, 'hp', alg.hp_bounds)
        self.baseline_purchase_costs['Turbines(s)'] = N * bst.CE / alg.CE * alg.cost(Pc_per_turbine)


class IsentropicTurbine(Turbine, new_graphics=False):
    """
    Create an isentropic turbine.

    Parameters
    ----------
    ins : 
        Inlet fluid.
    outs : 
        Outlet fluid.
    P : float
        Outlet pressure [Pa].
    eta : float
        Isentropic efficiency.
    vle : bool
        Whether to perform phase equilibrium calculations on
        the outflow. If False, the outlet will be assumed to be the same
        phase as the inlet.
    type: str
        Type of turbine. If None, the type
        will be determined automatically.

    """

    def _run(self):
        feed = self.ins[0]
        out = self.outs[0]
        out.copy_like(feed)
        if self.P < feed.P:
            if self.vle is True:
                out.vle(S=feed.S, P=self.P)
            else:
                out.P = self.P
                out.S = feed.S
            self.T_isentropic = out.T
            dH_isentropic = out.H - feed.H
            self.design_results['Ideal power'] = dH_isentropic / 3600. # kW
            self.design_results['Ideal duty'] = 0.
            dH_actual = dH_isentropic * self.eta
            H = feed.H + dH_actual
            if self.vle is True:
                out.vle(H=H, P=out.P)
            else:
                out.H = H
        else:
            warn(f"feed pressure ({feed.P:.5g} Pa) is lower or equal than outlet "
                 f"specification ({self.P:.5g} Pa); turbine {self.ID} is ignored", RuntimeWarning)

        
    def _design(self):
        super()._design()
        if 'Ideal power' in self.design_results: # Ideal power is not set if self.P > feed.P
            self._set_power(self.design_results['Ideal power'] * self.eta)
