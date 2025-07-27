# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2022, Yoel Cortes-Pena <yoelcortes@gmail.com>, Ben Portner <github.com/BenPortner>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
.. contents:: :local:

.. autoclass:: biosteam.units.compressor.Compressor
.. autoclass:: biosteam.units.compressor.IsothermalCompressor
.. autoclass:: biosteam.units.compressor.IsentropicCompressor
.. autoclass:: biosteam.units.compressor.PolytropicCompressor
.. autoclass:: biosteam.units.compressor.MultistageCompressor

References
----------
.. [1] Seider, W. D., Lewin,  D. R., Seader, J. D., Widagdo, S., Gani, R.,
    & Ng, M. K. (2017). Product and Process Design Principles. Wiley.
    Cost Accounting and Capital Cost Estimation (Chapter 16)
.. [2] Sinnott, R. and Towler, G (2019). "Chemical Engineering Design: SI Edition (Chemical Engineering Series)". 6th Edition. Butterworth-Heinemann.
.. [3] Schultz, J. (1962). "The Polytropic Analysis of Centrifugal Compressors". J. Eng. Power., 84(1): 69-82 (14 pages)
.. [4] Hundseid, O., Bakken, L. E. and Helde, T. (2006). “A Revised Compressor Polytropic Performance Analysis,” Proceedings of ASME GT2006, Paper Number 91033, ASME Turbo Expo 2006.

"""
import biosteam as bst
import numpy as np
from warnings import warn
from math import log, exp, ceil
from typing import NamedTuple, Tuple, Callable, Dict
from thermosteam.constants import R
from .heat_exchange import HX
from ..utils import list_available_names
from ..exceptions import DesignWarning, bounds_warning
from .. import Unit
from thermosteam._graphics import compressor_graphics
from thermosteam import VariableNode

__all__ = (
    'Compressor',
    'IsentropicCompressor', 
    'IsothermalCompressor', 
    'PolytropicCompressor',
    'MultistageCompressor'
)

#: TODO:
#: * Implement estimate of isentropic efficiency when not given (is this possible?).

class CompressorCostAlgorithm(NamedTuple): 
    #: Defines preliminary correlation algorithm for a compressor type
    psig_max: float #: Maximum achievable pressure in psig (to autodermine compressor type and/or issue warning)
    hp_bounds: Tuple[float, float] #: Horse power per machine (not a hard limit for costing, but included here for completion)
    acfm_bounds: Tuple[float, float] #: Actual cubic feet per minute (hard limit for parallel units)
    cost: Callable #: function(horse_power) -> Baseline purchase cost
    efficiencies: Dict[str, float] #: Heuristic efficiencies at 1,000 hp.
    driver: float #: Default driver (e.g., electric motor, steam turbine or gas turbine).
    CE: float #: Chemical engineering price cost index.


class Compressor(Unit, isabstract=True):
    """
    Abstract class for compressors that includes design and costing. Child classes
    should implement the `_run` method for mass and energy balances. Preliminary 
    design and costing is estimated according to [1]_.
    
    """
    _graphics = compressor_graphics
    _N_ins = 1
    _N_outs = 1
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
        'Compressor(s)': 2.15,
    }
    #: dict[str, CompressorCostAlgorithm] Cost algorithms by compressor type.
    baseline_cost_algorithms = { 
        'Screw': CompressorCostAlgorithm(
                psig_max=400.,
                acfm_bounds=(800., 2e4),
                hp_bounds=(10., 750.),
                cost=lambda Pc: exp(8.2496 + 0.7243 * log(Pc)),
                efficiencies={
                    'Electric motor': 0.80,
                    'Steam turbine': 0.65,
                    'Gas turbine': 0.35,
                },
                driver='Electric motor',
                CE=567,
            ),
        'Centrifugal': CompressorCostAlgorithm(
                psig_max=5e3,
                acfm_bounds=(1e3, 1.5e5),
                hp_bounds=(200., 3e4),
                cost=lambda Pc: exp(9.1553 + 0.63 * log(Pc)),
                efficiencies={
                    'Electric motor': 0.80,
                    'Steam turbine': 0.65,
                    'Gas turbine': 0.35,
                },
                driver='Steam turbine',
                CE=567,
            ),
        'Reciprocating': CompressorCostAlgorithm(
                psig_max=1e5,
                acfm_bounds=(5., 7000.),
                hp_bounds=(100., 20e3),
                cost=lambda Pc: exp(4.6762 + 1.23 * log(Pc)),
                efficiencies={
                    'Electric motor': 0.85,
                    'Steam turbine': 0.65,
                    'Gas turbine': 0.35,
                },
                driver='Electric motor',
                CE=567,
            ),
    }

    def _init(self, 
            P, eta=0.7, vle=False, compressor_type=None, 
            driver=None, material=None, driver_efficiency=None
        ):
        self.P = P  #: Outlet pressure [Pa].
        self.eta = eta  #: Isentropic efficiency.

        #: Whether to perform phase equilibrium calculations on the outflow.
        #: If False, the outlet will be assumed to be the same phase as the inlet.
        self.vle = vle
        self.material = 'Carbon steel' if material is None else material 
        self.compressor_type = 'Default' if compressor_type is None else compressor_type 
        self.driver = 'Default' if driver is None else driver
        self.driver_efficiency = 'Default' if driver_efficiency is None else driver_efficiency

    @property
    def compressor_type(self):
        return self._compressor_type
    @compressor_type.setter
    def compressor_type(self, compressor_type):
        """[str] Type of compressor. If 'Default', the type will be determined based on the outlet pressure."""
        compressor_type = compressor_type.capitalize()
        if compressor_type not in self.baseline_cost_algorithms and compressor_type != 'Default':
            raise ValueError(
                f"compressor type {repr(compressor_type)} not available; "
                f"only {list_available_names(self.baseline_cost_algorithms)} are available"
            )
        self._compressor_type = compressor_type

    @property
    def driver(self):
        return self._driver
    @driver.setter
    def driver(self, driver):
        """[str] Type of compressor. If 'Default', the type will be determined 
        based on type of compressor used. Centrifugal compressors default to 
        steam turbines while reciprocating and screw compressors default to 
        electric motors."""
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
        on the compressor type and the driver."""
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
            self.F_M['Compressor(s)'] = self.material_factors[material]
        except KeyError:
            raise AttributeError("material must be one of the following: "
                                 f"{list_available_names(self.material_factors)}")
        self._material = material

    def _determine_compressor_type(self):
        psig = (self.P - 101325.) * 14.6959 / 101325.
        cost_algorithms = self.baseline_cost_algorithms
        for name, alg in cost_algorithms.items():
            if psig < alg.psig_max: return name
        warn('no compressor available that is recommended for a pressure of '
            f'{self.P:.5g}; defaulting to {name.lower()} compressor', DesignWarning)
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
        driver = self.design_results['Driver']
        driver_efficiency = self._driver_efficiency
        if driver_efficiency == 'Default': 
            compressor_type = self.design_results['Type']
            alg = self.baseline_cost_algorithms[compressor_type]
            driver_efficiency = alg.efficiencies[driver]
        self.design_results['Driver efficiency'] = driver_efficiency
        if driver == 'Electric motor':
            self.power_utility(power / driver_efficiency)
        else:
            # The turbine produces the power that the compressor consumes.
            # This may not be the most elegant way showing this, but it 
            # makes it easy for design and costing.
            self.power_utility.consumption = self.power_utility.production = power / driver_efficiency
            if driver == 'Steam turbine':
                # Use high pressure steam utility as the driver. 
                # Assume that the recondenser cost is negligible and that 
                # heat integration is used (which are commonly the case).
                hps = bst.settings.get_heating_agent('high_pressure_steam')
                self.add_heat_utility(self.power_utility.consumption, T_in=298.15, agent=hps)
            elif driver == 'Gas turbine':
                # TODO: Possibly have an optional inlet stream that can work
                # as either steam or gas feed to the turbine.
                raise RuntimeError('gas turbine driver is not yet available in BioSTEAM')
            else:
                raise RuntimeError(f"invalid driver '{driver}'")
        
    def _design(self):
        design_results = self.design_results
        compressor_type = self.compressor_type 
        if compressor_type == 'Default': compressor_type = self._determine_compressor_type()
        design_results['Type'] = compressor_type
        alg = self.baseline_cost_algorithms[compressor_type]
        acfm_lb, acfm_ub = alg.acfm_bounds
        acfm = self.ins[0].get_total_flow('cfm')
        design_results['Compressors in parallel'] = ceil(acfm / acfm_ub) if acfm > acfm_ub else 1
        design_results['Driver'] = alg.driver if self._driver == 'Default' else self._driver 
    
    def _cost(self):
        # Note: Must run `_set_power` before running parent cost algorithm
        design_results = self.design_results
        alg = self.baseline_cost_algorithms[design_results['Type']]
        acfm_lb, acfm_ub = alg.acfm_bounds
        Pc = self.power_utility.get_property('consumption', 'hp')
        N = design_results['Compressors in parallel']
        Pc_per_compressor = Pc / N
        bounds_warning(self, 'power', Pc, 'hp', alg.hp_bounds, 'cost')
        if Pc_per_compressor < 1.:
            self.baseline_purchase_costs['Compressor(s)'] = 0.
        else:
            self.baseline_purchase_costs['Compressor(s)'] = N * bst.CE / alg.CE * alg.cost(Pc_per_compressor)
        self.F_D['Compressor(s)'] = self.design_factors[design_results['Driver']]


class IsothermalCompressor(Compressor, new_graphics=False):
    """
    Create an isothermal compressor.

    Parameters
    ----------
    ins : 
        Inlet fluid.
    outs : 
        Outlet fluid.
    P : float
        Outlet pressure [Pa].
    eta : float
        Isothermal efficiency.
    vle : bool
        Whether to perform phase equilibrium calculations on
        the outflow. If False, the outlet will be assumed to be the same
        phase as the inlet.
    type: str
        Type of compressor : blower/centrifugal/reciprocating. If None, the type
        will be determined automatically.

    Notes
    -----
    Default compressor selection, design and cost algorithms are adapted from [2]_.

    Examples
    --------
    Simulate reversible isothermal compression of gaseous hydrogen. Note that we set
    `include_excess_energies=True` to correctly account for the non-ideal behavior of
    hydrogen at high pressures. We further use the Soave-Redlich-Kwong (SRK) equation
    of state instead of the default Peng-Robinson (PR) because it is more accurate in
    this regime.

    >>> import biosteam as bst
    >>> from thermo import SRK
    >>> thermo = bst.Thermo([bst.Chemical('H2', eos=SRK)])
    >>> thermo.mixture.include_excess_energies = True
    >>> bst.settings.set_thermo(thermo)
    >>> feed = bst.Stream('feed', H2=1, T=298.15, P=20e5, phase='g')
    >>> K = bst.units.IsothermalCompressor('K', ins=feed, outs='outlet', P=350e5, eta=1)
    >>> K.simulate()
    >>> K.show()
    IsothermalCompressor: K
    ins...
    [0] feed
        phase: 'g', T: 298.15 K, P: 2e+06 Pa
        flow (kmol/hr): H2  1
    outs...
    [0] outlet
        phase: 'g', T: 298.15 K, P: 3.5e+07 Pa
        flow (kmol/hr): H2  1
    
    >>> K.results()
    Isothermal compressor                          Units               K
    Electricity         Power                         kW            2.47
                        Cost                      USD/hr           0.193
    Chilled water       Duty                       kJ/hr       -7.26e+03
                        Flow                     kmol/hr            7.53
                        Cost                      USD/hr          0.0363
    Design              Type                               Reciprocating
                        Compressors in parallel                        1
                        Driver                            Electric motor
                        Ideal power                   kW             2.1
                        Ideal duty                 kJ/hr       -7.26e+03
                        Driver efficiency                           0.85
    Purchase cost       Compressor(s)                USD             470
    Total purchase cost                              USD             470
    Utility cost                                  USD/hr            0.23

    """
    
    def _run(self):
        feed = self.ins[0]
        out = self.outs[0]
        out.copy_like(feed)
        out.P = self.P
        out.T = feed.T
        if self.vle is True: out.vle(T=out.T, P=out.P)
        self.ideal_power, self.ideal_duty = self._calculate_ideal_power_and_duty()

    def _design(self):
        super()._design()
        feed = self.ins[0]
        outlet = self.outs[0]
        ideal_power, ideal_duty = self._calculate_ideal_power_and_duty()
        Q = ideal_duty / self.eta
        self.add_heat_utility(unit_duty=Q, T_in=feed.T, T_out=outlet.T)
        self.design_results['Ideal power'] = ideal_power # kW
        self.design_results['Ideal duty'] = ideal_duty # kJ / hr
        self._set_power(ideal_power / self.eta)


class IsentropicCompressor(Compressor, new_graphics=False):
    """
    Create an isentropic compressor.

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
        Type of compressor : blower/centrifugal/reciprocating. If None, the type
        will be determined automatically.

    Notes
    -----
    Default compressor selection, design and cost algorithms are adapted from [2]_.

    Examples
    --------
    Simulate isentropic compression of gaseous hydrogen with 70% efficiency:

    >>> import biosteam as bst
    >>> bst.settings.set_thermo(["H2"])
    >>> feed = bst.Stream('feed', H2=1, T=25 + 273.15, P=101325, phase='g')
    >>> K = bst.units.IsentropicCompressor('K1', ins=feed, outs='outlet', P=50e5, eta=0.7)
    >>> K.simulate()
    >>> K.show()
    IsentropicCompressor: K1
    ins...
    [0] feed
        phase: 'g', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): H2  1
    outs...
    [0] outlet
        phase: 'g', T: 1152 K, P: 5e+06 Pa
        flow (kmol/hr): H2  1

    >>> K.results()
    Isentropic compressor                          Units             K1
    Electricity         Power                         kW              0
                        Cost                      USD/hr              0
    High pressure steam Duty                       kJ/hr           12.7
                        Flow                     kmol/hr       0.000396
                        Cost                      USD/hr       0.000126
    Design              Ideal power                   kW           4.92
                        Ideal duty                 kJ/hr              0
                        Type                                Centrifugal
                        Compressors in parallel                       1
                        Driver                            Steam turbine
                        Driver efficiency                          0.65
    Purchase cost       Compressor(s)                USD       5.87e+04
    Total purchase cost                              USD       5.87e+04
    Utility cost                                  USD/hr       0.000126


    Per default, the outlet phase is assumed to be the same as the inlet phase. If phase changes are to be accounted for,
    set `vle=True`:

    >>> import biosteam as bst
    >>> bst.settings.set_thermo(["H2O"])
    >>> feed = bst.MultiStream('feed', T=372.75, P=1e5, l=[('H2O', 0.1)], g=[('H2O', 0.9)])
    >>> K = bst.units.IsentropicCompressor('K2', ins=feed, outs='outlet', P=100e5, eta=1.0, vle=True)
    >>> K.simulate()
    >>> K.show()
    IsentropicCompressor: K2
    ins...
    [0] feed
        phases: ('g', 'l'), T: 372.75 K, P: 100000 Pa
        flow (kmol/hr): (g) H2O  0.9
                        (l) H2O  0.1
    outs...
    [0] outlet
        phases: ('g', 'l'), T: 798.63 K, P: 1e+07 Pa
        flow (kmol/hr): (g) H2O  1

    >>> K.results()
    Isentropic compressor                          Units             K2
    Electricity         Power                         kW              0
                        Cost                      USD/hr              0
    High pressure steam Duty                       kJ/hr            9.8
                        Flow                     kmol/hr       0.000305
                        Cost                      USD/hr       9.66e-05
    Design              Ideal power                   kW           5.42
                        Ideal duty                 kJ/hr              0
                        Type                                Centrifugal
                        Compressors in parallel                       1
                        Driver                            Steam turbine
                        Driver efficiency                          0.65
    Purchase cost       Compressor(s)                USD       4.98e+04
    Total purchase cost                              USD       4.98e+04
    Utility cost                                  USD/hr       9.66e-05

    >>> K.results()
    Isentropic compressor                          Units             K2
    Electricity         Power                         kW              0
                        Cost                      USD/hr              0
    High pressure steam Duty                       kJ/hr            9.8
                        Flow                     kmol/hr       0.000305
                        Cost                      USD/hr       9.66e-05
    Design              Ideal power                   kW           5.42
                        Ideal duty                 kJ/hr              0
                        Type                                Centrifugal
                        Compressors in parallel                       1
                        Driver                            Steam turbine
                        Driver efficiency                          0.65
    Purchase cost       Compressor(s)                USD       4.98e+04
    Total purchase cost                              USD       4.98e+04
    Utility cost                                  USD/hr       9.66e-05

    >>> K.results()
    Isentropic compressor                          Units             K2
    Electricity         Power                         kW              0
                        Cost                      USD/hr              0
    High pressure steam Duty                       kJ/hr            9.8
                        Flow                     kmol/hr       0.000305
                        Cost                      USD/hr       9.66e-05
    Design              Ideal power                   kW           5.42
                        Ideal duty                 kJ/hr              0
                        Type                                Centrifugal
                        Compressors in parallel                       1
                        Driver                            Steam turbine
                        Driver efficiency                          0.65
    Purchase cost       Compressor(s)                USD       4.98e+04
    Total purchase cost                              USD       4.98e+04
    Utility cost                                  USD/hr       9.66e-05


    Per default, the outlet phase is assumed to be the same as the inlet phase. If phase changes are to be accounted for,
    set `vle=True`:

    >>> import biosteam as bst
    >>> bst.settings.set_thermo(["H2O"])
    >>> feed = bst.MultiStream('feed', T=372.75, P=1e5, l=[('H2O', 0.1)], g=[('H2O', 0.9)])
    >>> K = bst.units.IsentropicCompressor('K2', ins=feed, outs='outlet', P=100e5, eta=1.0, vle=True)
    >>> K.simulate()
    >>> K.show()
    IsentropicCompressor: K2
    ins...
    [0] feed
        phases: ('g', 'l'), T: 372.75 K, P: 100000 Pa
        flow (kmol/hr): (g) H2O  0.9
                        (l) H2O  0.1
    outs...
    [0] outlet
        phases: ('g', 'l'), T: 798.63 K, P: 1e+07 Pa
        flow (kmol/hr): (g) H2O  1

    >>> K.results()
    Isentropic compressor                          Units             K2
    Electricity         Power                         kW              0
                        Cost                      USD/hr              0
    High pressure steam Duty                       kJ/hr            9.8
                        Flow                     kmol/hr       0.000305
                        Cost                      USD/hr       9.66e-05
    Design              Ideal power                   kW           5.42
                        Ideal duty                 kJ/hr              0
                        Type                                Centrifugal
                        Compressors in parallel                       1
                        Driver                            Steam turbine
                        Driver efficiency                          0.65
    Purchase cost       Compressor(s)                USD       4.98e+04
    Total purchase cost                              USD       4.98e+04
    Utility cost                                  USD/hr       9.66e-05

    """
    equation_node_names = (
        'overall_material_balance_node',
        'energy_balance_node',
        'isentropic_compression_phenomenode',
    )
    _energy_variable = 'T'

    @property
    def T_node(self):
        if hasattr(self, '_T_node'): return self._T_node
        self._T_node = var = VariableNode(f"{self.node_tag}.T", lambda: self.outs[0].T)
        return var 
    
    def get_E_node(self, stream):
        return self.T_node
    
    @property
    def E_node(self):
        return self.T_node

    @property
    def Q_node(self):
        if hasattr(self, '_Q_node'): return self._Q_node
        self._Q_node = var = VariableNode(f"{self.node_tag}.Q", lambda: self.Q)
        return var

    def initialize_isentropic_compression_phenomenode(self):
        self.isentropic_compression_phenomenode.set_equations(
            inputs=(
                self.T_node, 
                *[i.T_node for i in (*self.ins, *self.outs)],
                *[i.F_node for i in (*self.ins, *self.outs)],
            ),
            outputs=[self.Q_node],
        )

    def initialize_overall_material_balance_node(self):
        self.overall_material_balance_node.set_equations(
            inputs=[j for i in self.ins if (j:=i.F_node)],
            outputs=[i.F_node for i in self.outs],
        )
        
    def initialize_energy_balance_node(self):
        self.energy_balance_node.set_equations(
            inputs=(
                self.Q_node,
                *[i.T_node for i in (*self.ins, *self.outs)],
                *[i.F_node for i in (*self.ins, *self.outs)],
            ),
            outputs=[self.E_node],
        )

    def _run(self):
        feed = self.ins[0]
        out = self.outs[0]
        out.copy_like(feed)
        out.P = self.P
        out.S = feed.S
        self.T_isentropic = out.T
        if self.vle is True: out.vle(S=out.S, P=out.P)
        dH_isentropic = out.H - feed.H
        self.design_results['Ideal power'] = dH_isentropic / 3600. # kW
        self.design_results['Ideal duty'] = 0.
        self.Q = dH_actual = dH_isentropic / self.eta
        out.H = feed.H + dH_actual        
        if self.vle is True: out.vle(H=out.H, P=out.P)
        
    def _mass_and_energy_balance_specifications(self):
        return 'Compressor', [
            ('Isentropic efficiency', 100 * self.eta, '%'),
            ('P', self.P, 'Pa'),
        ]
        
    def _design(self):
        super()._design()
        self._set_power(self.design_results['Ideal power'] / self.eta)
    
    def _get_energy_departure_coefficient(self, stream):
        feed = self.ins[0]
        source = feed.source
        if source is None or source._energy_variable != 'T': return None
        return (stream, -stream.C)
    
    def _update_nonlinearities(self):
        feed = self.ins[0]
        out = self.outs[0]
        if len(out.phases) > 1:
            raise NotImplementedError('energy departure equation with multiple phase not yet implemented for isentropic compressors')
        mixture = self.mixture
        T_isentropic = getattr(self, 'T_isentropic', feed.T)
        self.T_isentropic = mixture.solve_T_at_SP(out.phase, out.mol, feed.S, T_isentropic, out.P)
        dH_isentropic = mixture('H', out, T=self.T_isentropic) - feed.H
        self.design_results['Ideal power'] = dH_isentropic / 3600. # kW
        self.design_results['Ideal duty'] = 0.
        self.Q = dH_isentropic / self.eta
        
    def _create_energy_departure_equations(self):
        feed = self.ins[0]
        out = self.outs[0]
        coeff = {self: out.C}
        feed._update_energy_departure_coefficient(coeff)
        return [(coeff, feed.H - out.H + self.Q)]
    
    def _create_material_balance_equations(self, composition_sensitive):
        fresh_inlets, process_inlets, equations = self._begin_equations(composition_sensitive)
        outlet, = self.outs
        if len(outlet.phases) > 1:
            raise NotImplementedError('energy departure equation with multiple phase not yet implemented for isentropic compressors')
        ones = np.ones(self.chemicals.size)
        minus_ones = -ones
        zeros = np.zeros(self.chemicals.size)
        
        # Overall flows
        eq_overall = {outlet: ones}
        for i in process_inlets: eq_overall[i] = minus_ones
        equations.append(
            (eq_overall, sum([i.mol for i in fresh_inlets], zeros))
        )
        return equations
    
    def _update_energy_variable(self, departure):
        self.outs[0].T += departure


class PolytropicCompressor(Compressor, new_graphics=False):
    """
    Create a polytropic compressor.

    Parameters
    ----------
    ins : 
        Inlet fluid.
    outs : 
        Outlet fluid.
    P : float
        Outlet pressure [Pa].
    eta : float
        Polytropic efficiency.
    vle : bool
        Whether to perform phase equilibrium calculations on
        the outflow. If False, the outlet will be assumed to be the same
        phase as the inlet.
    method: str
        'schultz'/'hundseid'. Calculation method for polytropic work. 'hundseid' is recommend
        for real gases at high pressure ratios.
    n_steps: int
        Number of virtual steps used in numerical integration for hundseid method.
    compressor_type: str
        Type of compressor : blower/centrifugal/reciprocating. If None, the type
        will be determined automatically.

    Notes
    -----
    Default compressor selection, design and cost algorithms are adapted from [2]_.

    Examples
    --------
    Simulate polytropic compression of hydrogen with 70% efficiency using the Schultz method [3]_:

    >>> import biosteam as bst
    >>> thermo = bst.Thermo([bst.Chemical('H2')])
    >>> thermo.mixture.include_excess_energies = True
    >>> bst.settings.set_thermo(thermo)
    >>> feed = bst.Stream('feed', H2=1, T=25 + 273.15, P=20e5, phase='g')
    >>> K = bst.units.PolytropicCompressor('K1', ins=feed, outs='outlet', P=350e5, eta=0.7, method='schultz')
    >>> K.simulate()
    >>> K.show(T='degC:.3g')
    PolytropicCompressor: K1
    ins...
    [0] feed
        phase: 'g', T: 25 degC, P: 2e+06 Pa
        flow (kmol/hr): H2  1
    outs...
    [0] outlet
        phase: 'g', T: 713 degC, P: 3.5e+07 Pa
        flow (kmol/hr): H2  1
    >>> K.results()
    Polytropic compressor                         Units              K1
    Electricity         Power                        kW            6.76
                        Cost                     USD/hr           0.529
    Design              Polytropic work                        2.07e+04
                        Type                              Reciprocating
                        Compressors in parallel                       1
                        Driver                           Electric motor
                        Driver efficiency                          0.85
    Purchase cost       Compressor(s)               USD        1.62e+03
    Total purchase cost                             USD        1.62e+03
    Utility cost                                 USD/hr           0.529


    Repeat using Hundseid method [4]_:

    >>> K = bst.units.PolytropicCompressor('K1', ins=feed, outs='outlet', P=350e5, eta=0.7, method='hundseid', n_steps=200)
    >>> K.simulate()
    >>> K.show()
    PolytropicCompressor: K1
    ins...
    [0] feed
        phase: 'g', T: 298.15 K, P: 2e+06 Pa
        flow (kmol/hr): H2  1
    outs...
    [0] outlet
        phase: 'g', T: 958.07 K, P: 3.5e+07 Pa
        flow (kmol/hr): H2  1
    >>> K.results()
    Polytropic compressor                         Units              K1
    Electricity         Power                        kW            6.48
                        Cost                     USD/hr           0.507
    Design              Polytropic work                        1.98e+04
                        Type                              Reciprocating
                        Compressors in parallel                       1
                        Driver                           Electric motor
                        Driver efficiency                          0.85
    Purchase cost       Compressor(s)               USD        1.54e+03
    Total purchase cost                             USD        1.54e+03
    Utility cost                                 USD/hr           0.507
    >>> K.results()
    Polytropic compressor                         Units              K1
    Electricity         Power                        kW            6.48
                        Cost                     USD/hr           0.507
    Design              Polytropic work                        1.98e+04
                        Type                              Reciprocating
                        Compressors in parallel                       1
                        Driver                           Electric motor
                        Driver efficiency                          0.85
    Purchase cost       Compressor(s)               USD        1.54e+03
    Total purchase cost                             USD        1.54e+03
    Utility cost                                 USD/hr           0.507


    Repeat using Hundseid method [4]_:

    >>> K = bst.units.PolytropicCompressor('K1', ins=feed, outs='outlet', P=350e5, eta=0.7, method='hundseid', n_steps=200)
    >>> K.simulate()
    >>> K.show()
    PolytropicCompressor: K1
    ins...
    [0] feed
        phase: 'g', T: 298.15 K, P: 2e+06 Pa
        flow (kmol/hr): H2  1
    outs...
    [0] outlet
        phase: 'g', T: 958.07 K, P: 3.5e+07 Pa
        flow (kmol/hr): H2  1
    
    >>> K.results()
    Polytropic compressor                         Units              K1
    Electricity         Power                        kW            6.48
                        Cost                     USD/hr           0.507
    Design              Polytropic work                        1.98e+04
                        Type                              Reciprocating
                        Compressors in parallel                       1
                        Driver                           Electric motor
                        Driver efficiency                          0.85
    Purchase cost       Compressor(s)               USD        1.54e+03
    Total purchase cost                             USD        1.54e+03
    Utility cost                                 USD/hr           0.507


    """
    available_methods = {'schultz', 'hundseid'}
    def _init(self, P, eta=0.7, 
                 vle=False, compressor_type=None, method=None, n_steps=None):
        Compressor._init(self, P=P, eta=eta, vle=vle, compressor_type=compressor_type)
        self.method = "schultz" if method is None else method
        self.n_steps = 100 if n_steps is None else n_steps
        
    @property
    def method(self):
        return self._method
    @method.setter
    def method(self, method):
        method = method.lower()
        if method not in self.available_methods:
            raise ValueError(
                f"method {repr(method)} not available; "
                f"only {list_available_names(self.available_methods)} are available"
            )
        self._method = method
        
    def _schultz(self):
        # calculate polytropic work using Schultz method
        feed = self.ins[0]
        out = self.outs[0]

        # calculate polytropic exponent and real gas correction factor
        out.P = self.P
        out.S = feed.S
        k = log(out.P / feed.P) / log(feed.V / out.V)
        n_1_n = (k - 1) / k # n: polytropic exponent
        W_poly = feed.P * feed.V / n_1_n * ((out.P/feed.P)**n_1_n - 1) # kJ/kmol
        W_isen = (out.H - feed.H) # kJ/kmol
        f = W_isen/W_poly # f: correction factor for real gases

        # calculate non-reversible polytropic work (accounting for eta)
        n_1_n = n_1_n / self.eta
        W_actual = f * feed.P * feed.V / n_1_n * ((out.P/feed.P)**n_1_n - 1) / self.eta * out.F_mol # kJ/kmol -> kJ/hr

        # calculate outlet state
        out.H = feed.H + W_actual # kJ/hr

        return W_actual # kJ/hr

    def _hundseid(self):
        # calculate polytropic work using Hundseid method
        feed = self.ins[0].copy()
        out = self.outs[0]
        n_steps = self.n_steps

        pr = (self.P / feed.P) ** (1 / n_steps) # pressure ratio between discrete steps
        W_actual = 0
        for i in range(n_steps):
            # isentropic pressure change
            out.P = feed.P * pr
            out.S = feed.S
            dH_isen_i = out.H - feed.H
            # efficiency correction
            dH_i = dH_isen_i / self.eta
            out.H = feed.H + dH_i
            W_actual += dH_i # kJ/hr
            # next step
            feed.P = out.P
            feed.T = out.T

        return W_actual # kJ/hr

    def _run(self):
        feed = self.ins[0]
        out = self.outs[0]
        out.copy_like(feed)
        name = '_' + self._method
        method = getattr(self, name)
        self.design_results['Polytropic work'] = method() # Polytropic work [kJ/hr]
        if self.vle is True: out.vle(H=out.H, P=out.P)
        
    def _design(self):
        super()._design()
        self._set_power(self.design_results['Polytropic work'] / 3600) # kJ/hr -> kW


class MultistageCompressor(Unit):
    """
    Create a multistage compressor. Models multistage polytropic or isentropic compression with intermittent cooling.

    There are two setup options:

    * Option 1: Define `pr` and `n_stages` (optionally `eta`, `vle`, `type`).
        Creates `n_stages` identical isentropic compressors. Each compressor is followed by
        a heat exchanger, which cools the effluent to inlet temperature.

    * Option 2: Define `compressors` and `hxs`.
        Takes a list of pre-defined compressors and heat exchangers and connects them in series. This option
        allows more flexibility in terms of the type of compressor (isentropic/polytropic) and parameterization
        (e.g. each stage can have different efficiencies and outlet temperatures).

    Parameters
    ----------
    ins : 
        Inlet fluid.
    outs : 
        Outlet fluid.
    pr : float
        (setup option 1) Pressure ratio between isentropic stages.
    n_stages: float
        (setup option 1) Number of isentropic stages.
    eta : float
        (setup option 1) Isentropic efficiency.
    vle : bool
        (setup option 1) Whether to perform phase equilibrium calculations on
        the outflow of each stage. If False, the outlet will be assumed to be the same
        phase as the inlet.
    compressor_type: str
        (setup option 1) Type of compressor : blower/centrifugal/reciprocating. If None, the type
        will be determined automatically.
    compressors: list[_CompressorBase]
        (setup option 2) List of compressors to use for each stage.
    hxs: list[HX]
        (setup option 2) List of heat exchangers to use for each stage.

    Notes
    -----
    Default compressor selection, design and cost algorithms are adapted from [2]_.

    Examples
    --------
    Simulate multistage compression of gaseous hydrogen (simple setup). Hydrogen is compressed
    isentropically (with an efficiency of 70%) from 20 bar to 320 bar in four stages
    (pressure ratio of two in each stage):

    >>> import biosteam as bst
    >>> thermo = bst.Thermo([bst.Chemical('H2')])
    >>> thermo.mixture.include_excess_energies = True
    >>> bst.settings.set_thermo(thermo)
    >>> feed = bst.Stream('feed', H2=1, T=298.15, P=20e5, phase='g')
    >>> K = bst.units.MultistageCompressor('K', ins=feed, outs='outlet', pr=2, n_stages=4, eta=0.7)
    >>> K.simulate()
    >>> K.show()
    MultistageCompressor: K
    ins...
    [0] feed
        phase: 'g', T: 298.15 K, P: 2e+06 Pa
        flow (kmol/hr): H2  1
    outs...
    [0] outlet
        phase: 'g', T: 298.15 K, P: 3.2e+07 Pa
        flow (kmol/hr): H2  1
    
    >>> K.results()
    Multistage compressor                       Units                      K
    Electricity         Power                      kW                      0
                        Cost                   USD/hr                      0
    High pressure steam Duty                    kJ/hr                   5.68
                        Flow                  kmol/hr               0.000177
                        Cost                   USD/hr                5.6e-05
    Chilled water       Duty                    kJ/hr              -1.12e+04
                        Flow                  kmol/hr                   7.42
                        Cost                   USD/hr                 0.0559
    Design              Type                           Multistage compressor
                        Area                     ft^2                   1.54
    Purchase cost       K k1 - Compressor(s)      USD               1.45e+04
                        K h1 - Double pipe        USD                    568
                        K k2 - Compressor(s)      USD               1.46e+04
                        K h2 - Double pipe        USD                    675
                        K k3 - Compressor(s)      USD               1.48e+04
                        K h3 - Double pipe        USD                    956
                        K k4 - Compressor(s)      USD               1.52e+04
                        K h4 - Double pipe        USD               1.78e+03
    Total purchase cost                           USD                6.3e+04
    Utility cost                               USD/hr                  0.056

    Show the fluid state at the outlet of each heat exchanger:
    
    >>> for hx in K.hxs:
    ...  hx.outs[0].show()
    Stream: K_H1__K_K2 from <HXutility: K_H1> to <IsentropicCompressor: K_K2>
    phase: 'g', T: 298.15 K, P: 4e+06 Pa
    flow (kmol/hr): H2  1
    Stream: K_H2__K_K3 from <HXutility: K_H2> to <IsentropicCompressor: K_K3>
    phase: 'g', T: 298.15 K, P: 8e+06 Pa
    flow (kmol/hr): H2  1
    Stream: K_H3__K_K4 from <HXutility: K_H3> to <IsentropicCompressor: K_K4>
    phase: 'g', T: 298.15 K, P: 1.6e+07 Pa
    flow (kmol/hr): H2  1
    Stream: outlet from <MultistageCompressor: K>
    phase: 'g', T: 298.15 K, P: 3.2e+07 Pa
    flow (kmol/hr): H2  1

    If we want to setup more complex multistage compression schemes, we can pre-define the compressors and
    heat exchangers and pass them as a list to `MultistageCompressor`:

    >>> ks = [
    ...     bst.units.IsentropicCompressor(P=30e5, eta=0.6),
    ...     bst.units.PolytropicCompressor(P=50e5, eta=0.65),
    ...     bst.units.IsentropicCompressor(P=90e5, eta=0.70),
    ...     bst.units.PolytropicCompressor(P=170e5, eta=0.75),
    ...     bst.units.IsentropicCompressor(P=320e5, eta=0.80),
    ... ]
    >>> hxs = [bst.units.HXutility(T=T) for T in [310, 350, 400, 350, 298]]
    >>> K = bst.units.MultistageCompressor('K2', ins=feed, outs='outlet', compressors=ks, hxs=hxs)
    >>> K.simulate()
    >>> K.show()
    MultistageCompressor: K2
    ins...
    [0] feed
        phase: 'g', T: 298.15 K, P: 2e+06 Pa
        flow (kmol/hr): H2  1
    outs...
    [0] outlet
        phase: 'g', T: 298 K, P: 3.2e+07 Pa
        flow (kmol/hr): H2  1
    
    >>> K.results()
    Multistage compressor                        Units                     K2
    Electricity         Power                       kW                      0
                        Cost                    USD/hr                      0
    High pressure steam Duty                     kJ/hr                   6.48
                        Flow                   kmol/hr               0.000202
                        Cost                    USD/hr               6.39e-05
    Chilled water       Duty                     kJ/hr              -5.63e+03
                        Flow                   kmol/hr                   3.73
                        Cost                    USD/hr                 0.0282
    Cooling water       Duty                     kJ/hr              -7.15e+03
                        Flow                   kmol/hr                   4.88
                        Cost                    USD/hr                0.00238
    Design              Type                            Multistage compressor
                        Area                      ft^2                   1.14
    Purchase cost       K2 k1 - Compressor(s)      USD               1.11e+04
                        K2 h1 - Double pipe        USD                    297
                        K2 k2 - Compressor(s)      USD                1.3e+04
                        K2 h2 - Double pipe        USD                    199
                        K2 k3 - Compressor(s)      USD               1.44e+04
                        K2 h3 - Double pipe        USD                    129
                        K2 k4 - Compressor(s)      USD               1.64e+04
                        K2 h4 - Double pipe        USD                    753
                        K2 k5 - Compressor(s)      USD               1.45e+04
                        K2 h5 - Double pipe        USD               2.03e+03
    Total purchase cost                            USD               7.27e+04
    Utility cost                                USD/hr                 0.0306

    Show the fluid state at the outlet of each heat exchanger:
    
    >>> for hx in K.hxs:
    ...  hx.outs[0].show()
    Stream: K2_H1__K2_K2 from <HXutility: K2_H1> to <PolytropicCompressor: K2_K2>
    phase: 'g', T: 310 K, P: 3e+06 Pa
    flow (kmol/hr): H2  1
    Stream: K2_H2__K2_K3 from <HXutility: K2_H2> to <IsentropicCompressor: K2_K3>
    phase: 'g', T: 350 K, P: 5e+06 Pa
    flow (kmol/hr): H2  1
    Stream: K2_H3__K2_K4 from <HXutility: K2_H3> to <PolytropicCompressor: K2_K4>
    phase: 'g', T: 400 K, P: 9e+06 Pa
    flow (kmol/hr): H2  1
    Stream: K2_H4__K2_K5 from <HXutility: K2_H4> to <IsentropicCompressor: K2_K5>
    phase: 'g', T: 350 K, P: 1.7e+07 Pa
    flow (kmol/hr): H2  1
    Stream: outlet from <MultistageCompressor: K2>
    phase: 'g', T: 298 K, P: 3.2e+07 Pa
    flow (kmol/hr): H2  1

    Show the fluid state at the outlet of each heat exchanger:
    
    >>> for hx in K.hxs:
    ...  hx.outs[0].show()
    Stream: K2_H1__K2_K2 from <HXutility: K2_H1> to <PolytropicCompressor: K2_K2>
    phase: 'g', T: 310 K, P: 3e+06 Pa
    flow (kmol/hr): H2  1
    Stream: K2_H2__K2_K3 from <HXutility: K2_H2> to <IsentropicCompressor: K2_K3>
    phase: 'g', T: 350 K, P: 5e+06 Pa
    flow (kmol/hr): H2  1
    Stream: K2_H3__K2_K4 from <HXutility: K2_H3> to <PolytropicCompressor: K2_K4>
    phase: 'g', T: 400 K, P: 9e+06 Pa
    flow (kmol/hr): H2  1
    Stream: K2_H4__K2_K5 from <HXutility: K2_H4> to <IsentropicCompressor: K2_K5>
    phase: 'g', T: 350 K, P: 1.7e+07 Pa
    flow (kmol/hr): H2  1
    Stream: outlet from <MultistageCompressor: K2>
    phase: 'g', T: 298 K, P: 3.2e+07 Pa
    flow (kmol/hr): H2  1

    If we want to setup more complex multistage compression schemes, we can pre-define the compressors and
    heat exchangers and pass them as a list to `MultistageCompressor`:

    >>> ks = [
    ...     bst.units.IsentropicCompressor(P=30e5, eta=0.6),
    ...     bst.units.PolytropicCompressor(P=50e5, eta=0.65),
    ...     bst.units.IsentropicCompressor(P=90e5, eta=0.70),
    ...     bst.units.PolytropicCompressor(P=170e5, eta=0.75),
    ...     bst.units.IsentropicCompressor(P=320e5, eta=0.80),
    ... ]
    >>> hxs = [bst.units.HXutility(T=T) for T in [310, 350, 400, 350, 298]]
    >>> K = bst.units.MultistageCompressor('K2', ins=feed, outs='outlet', compressors=ks, hxs=hxs)
    >>> K.simulate()
    >>> K.show()
    MultistageCompressor: K2
    ins...
    [0] feed
        phase: 'g', T: 298.15 K, P: 2e+06 Pa
        flow (kmol/hr): H2  1
    outs...
    [0] outlet
        phase: 'g', T: 298 K, P: 3.2e+07 Pa
        flow (kmol/hr): H2  1
    
    >>> K.results()
    Multistage compressor                        Units                     K2
    Electricity         Power                       kW                      0
                        Cost                    USD/hr                      0
    High pressure steam Duty                     kJ/hr                   6.48
                        Flow                   kmol/hr               0.000202
                        Cost                    USD/hr               6.39e-05
    Chilled water       Duty                     kJ/hr              -5.63e+03
                        Flow                   kmol/hr                   3.73
                        Cost                    USD/hr                 0.0282
    Cooling water       Duty                     kJ/hr              -7.15e+03
                        Flow                   kmol/hr                   4.88
                        Cost                    USD/hr                0.00238
    Design              Type                            Multistage compressor
                        Area                      ft^2                   1.14
    Purchase cost       K2 k1 - Compressor(s)      USD               1.11e+04
                        K2 h1 - Double pipe        USD                    297
                        K2 k2 - Compressor(s)      USD                1.3e+04
                        K2 h2 - Double pipe        USD                    199
                        K2 k3 - Compressor(s)      USD               1.44e+04
                        K2 h3 - Double pipe        USD                    129
                        K2 k4 - Compressor(s)      USD               1.64e+04
                        K2 h4 - Double pipe        USD                    753
                        K2 k5 - Compressor(s)      USD               1.45e+04
                        K2 h5 - Double pipe        USD               2.03e+03
    Total purchase cost                            USD               7.27e+04
    Utility cost                                USD/hr                 0.0306

    """

    _N_ins = 1
    _N_outs = 1
    _units = {
        **Compressor._units,
        **HX._units,
    }

    def _init(
            self, pr=None, n_stages=None, eta=0.7, vle=False, compressor_type=None,
            compressors=None, hxs=None,
        ):
        # setup option 1: list of compressors and list of heat exchangers
        if compressors is not None and hxs is not None:
            if not isinstance(compressors[0], Compressor):
                print(compressors[0].__class__)
                raise RuntimeError(f"invalid parameterization of {self.ID}: `compressors` must "
                                   f"be a list of compressor objects.")
            elif not isinstance(hxs[0], bst.HX):
                raise RuntimeError(f"invalid parameterization of {self.ID}: `hxd` must "
                                   f"be a list of heat exchanger objects.")
            elif len(compressors) != len(hxs):
                raise RuntimeError(f"invalid parameterization of {self.ID}: `compressors` and `hxs` "
                                   f"must have the same length.")
            else:
                self.compressors = tuple(compressors)
                self.hxs = tuple(hxs)
                self.pr = None
                self.n_stages = None

        # setup option 2: fixed pressure ratio and number of stages
        elif pr is not None and n_stages is not None:
            self.pr = pr
            self.n_stages = n_stages
            self.eta = eta
            self.vle=vle
            self.compressor_type=compressor_type
            self.compressors = None
            self.hxs = None
        else:
            raise RuntimeError(f"invalid parameterization of {self.ID}: Must specify `pr` and "
                               f"`n_stages` or `compressors` and `hxs`.")
        self._old_specifications = None

    def _overwrite_subcomponent_id(self, subcomponent, i_stage):
        # overwrite subcomponent id
        ID = f"{self.ID}_{subcomponent.ticket_name}{i_stage}"
        subcomponent.ID = ID

        # overwrite inlet id if not multistage inlet
        if i_stage == 1 and isinstance(subcomponent, Compressor):
            pass
        else:
            subcomponent.ins[0].ID = f"{subcomponent.ins[0].ID}__{ID}"

        # overwrite outlet id if not multistage outlet
        if i_stage == (self.n_stages or len(self.compressors)) and isinstance(subcomponent, bst.HX):
            pass
        else:
            subcomponent.outs[0].ID = f"{ID}"

    def reset_cache(self, **kwargs):
        super().reset_cache(**kwargs)
        self._old_specifications = None

    def _setup(self):
        super()._setup()
        feed = self._ins[0]
        compressors = self.compressors
        hxs = self.hxs
        pr = self.pr
        n_stages = self.n_stages
        specifications = (pr, n_stages, compressors, hxs)
        if (specifications == self._old_specifications):
            return # Skip setup (already done)
        
        # setup option 1: rewire compressors and heat exchangers
        if pr is None and n_stages is None:
            last_hx = None
            for n, (c, hx) in enumerate(zip(compressors, hxs)):
                if last_hx is not None: c.ins[0] = last_hx.outs[0]
                hx.ins[0] = c.outs[0]
                last_hx = hx
                self._overwrite_subcomponent_id(c, n+1)
                self._overwrite_subcomponent_id(hx, n+1)
        # setup option 2: create connected compressor and hx objects
        else:
            T = feed.T
            compressors = []
            hxs = []
            
            # Temporarily register all units/streams in this flowsheet 
            # (instead of the main flowsheet) to prevent system creation problems
            self.flowsheet = bst.Flowsheet('Multistage_compressor_' + self.ID)
            with self.flowsheet.temporary(): 
                hx = None; P = feed.P
                for n in range(self.n_stages):
                    inflow = hx.outs[0] if hx else None
                    P *= pr
                    c = IsentropicCompressor(
                        ins=inflow, P=P, eta=self.eta,
                        vle=self.vle, compressor_type=self.compressor_type
                    )
                    self._overwrite_subcomponent_id(c, n+1)
                    hx = bst.HXutility(
                        ins=c.outs[0], T=T, rigorous=self.vle
                    )
                    self._overwrite_subcomponent_id(hx, n+1)
                    compressors.append(c)
                    hxs.append(hx)
            self.compressors = compressors = tuple(compressors)
            self.hxs =  hxs = tuple(hxs)

        # set inlet and outlet reference
        compressors[0]._ins = self._ins
        hxs[-1]._outs = self._outs

        # set auxillary units
        units = [u for t in zip(compressors, hxs) for u in t]
        self.auxiliary_unit_names = tuple([u.ID for u in units])
        for u in units: self.__setattr__(u.ID, u)
        self._old_specifications = (pr, n_stages, compressors, hxs)

    def _run(self):
        # calculate results

        # helper variables
        units = [u for t in zip(self.compressors, self.hxs) for u in t]

        # simulate all subcomponents
        for u in units: 
            u._setup() 
            u._run()

    def _design(self):
        self.design_results["Type"] = "Multistage compressor"

        # design all subcomponents
        units = [u for t in zip(self.compressors, self.hxs) for u in t]
        for u in units: u._summary()
        
        # sum up design values
        sum_fields = [
            "Power", "Duty",
            "Area", "Tube side pressure drop", "Shell side pressure drop"
        ]
        for u in units:
            for k,v in u.design_results.items():
                if k in sum_fields:
                    if k in self.design_results:
                        self.design_results[k] += v
                    else:
                        self.design_results[k] = v
