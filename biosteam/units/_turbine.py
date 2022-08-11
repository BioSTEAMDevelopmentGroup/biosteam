# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2022, Yoel Cortes-Pena <yoelcortes@gmail.com>, Ben Portner <github.com/BenPortner>
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
from .. import Unit
from warnings import warn
from math import log, exp, ceil
from typing import NamedTuple, Tuple, Callable, Dict
from .heat_exchange import HX
from ..utils import list_available_names
from ..exceptions import DesignWarning, bounds_warning

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
    efficiencies: Dict[str, float] #: Heuristic efficiencies at 1,000 hp.
    driver: float #: Default driver (e.g., electric motor, steam turbine or gas turbine).
    CE: float #: Chemical engineering price cost index.


class Turbine(Unit, isabstract=True):
    """
    Abstract class for turbines that includes design and costing. Child classes
    should implement the `_run` method for mass and energy balances. Preliminary 
    design and costing is estimated according to [1]_.
    
    """
    _N_ins = 1
    _N_outs = 1
    _N_heat_utilities = 1
    _units = {
        'Ideal power': 'kW',
        'Ideal duty': 'kJ/hr',
    }
    design_factors = {
        'Electric motor': 1.0,
        'Steam turbine': 1.15,
        'Gas turbine': 1.25,
    }
    material_factors = {
        'Carbon steel': 1.0,
        'Stainless steel': 2.5,
        'Nickel alloy': 5.0,
    }
    _F_BM_default = {
        'Turbine(s)': 2.15,
    }
    #: dict[str, TurbineCostAlgorithm] Cost algorithms by turbine type.
    baseline_cost_algorithms = {
        'Default': TurbineCostAlgorithm(
                psig_max=999999999.,
                acfm_bounds=(0., 999999999.),
                hp_bounds=(0., 999999999.),
                cost=lambda Pc: 0,
                efficiencies={
                    'Electric motor': 0.80,
                    'Steam turbine': 0.65,
                    'Gas turbine': 0.35,
                },
                driver='Gas turbine',
                CE=567,
            ),
    }

    def __init__(self, ID='', ins=None, outs=(), thermo=None, *, 
                 P, eta=0.3, vle=False, turbine_type=None,
                 driver=None, material=None, driver_efficiency=None):
        Unit.__init__(self, ID, ins, outs, thermo)
        self.P = P  #: Outlet pressure [Pa].
        self.eta = eta  #: Isentropic efficiency.

        #: Whether to perform phase equilibrium calculations on the outflow.
        #: If False, the outlet will be assumed to be the same phase as the inlet.
        self.vle = vle
        self.material = 'Carbon steel' if material is None else material 
        self.turbine_type = 'Default' if turbine_type is None else turbine_type
        self.driver = 'Default' if driver is None else driver
        self.driver_efficiency = 'Default' if driver_efficiency is None else driver_efficiency

    @property
    def turbine_type(self):
        return self._turbine_type
    @turbine_type.setter
    def turbine_type(self, turbine_type):
        """[str] Type of turbine. If 'Default', the type will be determined based on the outlet pressure."""
        turbine_type = turbine_type.capitalize()
        if turbine_type not in self.baseline_cost_algorithms and turbine_type != 'Default':
            raise ValueError(
                f"turbine type {repr(turbine_type)} not available; "
                f"only {list_available_names(self.baseline_cost_algorithms)} are available"
            )
        self._turbine_type = turbine_type

    @property
    def driver(self):
        return self._driver
    @driver.setter
    def driver(self, driver):
        """[str] Type of turbine. If 'Default', the type will be determined
        based on type of turbine used."""
        driver = driver.capitalize()
        if driver not in self.design_factors and driver != 'Default':
            raise ValueError(
                f"driver {repr(driver)} not available; "
                f"only {list_available_names(self.design_factors)} are available"
            )
        self._driver = driver

    @property
    def driver_efficiency(self):
        return self._driver_efficiency
    @driver_efficiency.setter
    def driver_efficiency(self, driver_efficiency):
        """[str] Efficiency of driver (e.g., steam turbine or electric motor). 
        If 'Default', a heuristic efficiency will be selected based
        on the turbine type and the driver."""
        if isinstance(driver_efficiency, str):
            if driver_efficiency != 'Default':
                raise ValueError(
                    f"driver efficiency must be a number or 'Default'; not {repr(driver_efficiency)}"
                )
        else:
            driver_efficiency = float(driver_efficiency)
        self._driver_efficiency = driver_efficiency

    @property
    def material(self):
        """[str]. Construction material. Defaults to 'Carbon steel'."""
        return self._material
    @material.setter
    def material(self, material):
        try:
            self.F_M['Turbines(s)'] = self.material_factors[material]
        except KeyError:
            raise AttributeError("material must be one of the following: "
                                 f"{list_available_names(self.material_factors)}")
        self._material = material

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
        driver_efficiency = self._driver_efficiency
        if driver_efficiency == 'Default': 
            turbine_type = self.design_results['Type']
            driver = self.design_results['Driver']
            alg = self.baseline_cost_algorithms[turbine_type]
            driver_efficiency = alg.efficiencies[driver]
        self.design_results['Driver efficiency'] = driver_efficiency
        self.power_utility.consumption = power / driver_efficiency
    
    def _design(self):
        design_results = self.design_results
        turbine_type = self.turbine_type
        if turbine_type == 'Default': turbine_type = self._determine_turbine_type()
        design_results['Type'] = turbine_type
        alg = self.baseline_cost_algorithms[turbine_type]
        acfm_lb, acfm_ub = alg.acfm_bounds
        acfm = self.ins[0].get_total_flow('cfm')
        design_results['Turbines in parallel'] = ceil(acfm / acfm_ub) if acfm > acfm_ub else 1
        design_results['Driver'] = alg.driver if self._driver == 'Default' else self._driver 
    
    def _cost(self):
        # Note: Must run `_set_power` before running parent cost algorithm
        design_results = self.design_results
        alg = self.baseline_cost_algorithms[design_results['Type']]
        acfm_lb, acfm_ub = alg.acfm_bounds
        Pc = abs(self.power_utility.get_property('consumption', 'hp'))
        N = design_results['Turbines in parallel']
        Pc_per_turbine = Pc / N
        bounds_warning(self, 'power', Pc, 'hp', alg.hp_bounds, 'cost')
        self.baseline_purchase_costs['Turbines(s)'] = N * bst.CE / alg.CE * alg.cost(Pc_per_turbine)
        self.F_D['Turbines(s)'] = self.design_factors[design_results['Driver']]


class IsentropicTurbine(Turbine):
    """
    Create an isentropic turbine.

    Parameters
    ----------
    ins : stream
        Inlet fluid.
    outs : stream
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
    _N_heat_utilities = 0

    def _run(self):
        feed = self.ins[0]
        out = self.outs[0]
        out.copy_like(feed)
        out.P = self.P
        out.S = feed.S
        if self.vle is True: out.vle(S=out.S, P=out.P)
        self.T_isentropic = out.T
        dH_isentropic = out.H - feed.H
        self.design_results['Ideal power'] = dH_isentropic / 3600. # kW
        self.design_results['Ideal duty'] = 0.
        dH_actual = dH_isentropic * self.eta
        out.H = feed.H + dH_actual        
        if self.vle is True: out.vle(H=out.H, P=out.P)
        
    def _design(self):
        super()._design()
        self._set_power(self.design_results['Ideal power'] * self.eta)
