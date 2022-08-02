# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
import biosteam as bst
from .. import Unit
from .decorators import cost
from warnings import warn
from math import log, exp, ceil
from typing import NamedTuple, Tuple, Callable, Dict
from ..utils import list_available_names
from ..exceptions import DesignWarning

__all__ = (
    'Compressor',
    'IsentropicCompressor', 
    'IsothermalCompressor', 
    'PolytropicCompressor',
)

#: TODO:
#: * Implement estimate of isentropic efficiency when not given (is this possible?).

class CompressorCostAlgorithm(NamedTuple): 
    #: Defines preliminary correlation algorithm for however many stages
    psig_max: float #: Maximum achievable pressure in psig (to autodermine compressor type and/or issue warning)
    hp_bounds: Tuple[float, float] #: Horse power per machine (not a hard limit for costing, but included here for completion)
    acfm_bounds: Tuple[float, float] #: Actual cubic feet per minute (hard limit for parallel units)
    cost: Callable #: function(horse_power) -> Baseline purchase cost (includes multiple stages/machines)
    efficiencies: Dict[str, float] #: Heuristic efficiencies of compressor types at 1,000 hp.
    driver: float #: Default driver (e.g., electric motor, steam turbine or gas turbine).
    CE: float #: Chemical engineering price cost index.


class Compressor(Unit, isabstract=True):
    """
    Abstract class for compressors that includes design and costing. Child classes
    should implement the `_run` method for mass and energy balances. Preliminary 
    design and costing is estimated according to [1]_.
    
    References
    ----------
    .. [1] Seider, W. D., Lewin,  D. R., Seader, J. D., Widagdo, S., Gani, R.,
    & Ng, M. K. (2017). Product and Process Design Principles. Wiley.
    Cost Accounting and Capital Cost Estimation (Chapter 16)
    
    """
    _N_ins = 1
    _N_outs = 1
    _N_heat_utilities = 1
    _units = {
        'Ideal power': 'kW',
        'Ideal duty': 'kJ/hr',
    }
    _F_D = { # Design factors
        'Electric motor': 1.0,
        'Steam turbine': 1.15,
        'Gas turbine': 1.25,
    }
    _F_M = { # Material factors
        'Carbon steel': 1.0,
        'Stainless steel': 2.5,
        'Nickel alloy': 5.0,
    }
    _F_BM_default = {
        'Compressor(s)': 2.15,
    }
    baseline_cost_algorithms = { #: [str] Compressor type: [CompressorCostAlgorithm]
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
        'Recriprocating': CompressorCostAlgorithm(
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

    def __init__(self, ID='', ins=None, outs=(), thermo=None, *, 
                 P, eta=0.7, vle=False, compressor_type=None, 
                 driver=None, material=None):
        Unit.__init__(self, ID, ins, outs, thermo)
        self.P = P  #: Outlet pressure [Pa].
        self.eta = eta  #: Isentropic efficiency.

        #: Whether to perform phase equilibrium calculations on the outflow.
        #: If False, the outlet will be assumed to be the same phase as the inlet.
        self.vle = vle
        self.material = 'Carbon steel' if material is None else material 
        self.compressor_type = 'Default' if compressor_type is None else compressor_type 
        self.driver = 'Default' if driver is None else driver

    @property
    def compressor_type(self):
        return self._compressor_type
    @compressor_type.setter
    def compressor_type(self, compressor_type):
        """[str] Type of compressor. If 'Default', the type will be determined based on power."""
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
        """[str] Type of compressor. If 'Default', the type will be determined based on power."""
        driver = driver.capitalize()
        if driver not in self._F_D and driver != 'Default':
            raise ValueError(
                f"driver {repr(driver)} not available; "
                f"only {list_available_names(self._F_D)} are available"
            )
        self._driver = driver

    @property
    def material(self):
        """Defaults to 'Carbon steel'"""
        return self._material
    @material.setter
    def material(self, material):
        try:
            self.F_M['Compressor(s)'] = self._F_M[material]
        except KeyError:
            raise AttributeError("material must be one of the following: "
                                 f"{list_available_names(self._F_M)}")
        self._material = material

    def _determine_compressor_type(self):
        psig = (self.P - 101325) * 14.6959
        cost_algorithms = self.baseline_cost_algorithms
        for name, alg in cost_algorithms.items():
            if psig < alg.psig_max: return name
        warn('no compressor available that is recommended for a pressure of '
            f'{self.P: .5g}; defaulting to {repr(name)}', DesignWarning)
        return name

    def _calculate_ideal_power_and_duty(self):
        feed = self.ins[0]
        out = self.outs[0]
        dH = out.H - feed.H
        Q = TdS = feed.T * (out.S - feed.S)  # Duty [kJ/hr]
        power_ideal = (dH - TdS) / 3600.  # Power [kW]
        return power_ideal, Q

    def _design(self):
        # Note: Must set power utility before running parent design algorithm
        design_results = self.design_results
        compressor_type = self.compressor_type 
        if compressor_type == 'Default': compressor_type = self._determine_compressor_type()
        design_results['Type'] = compressor_type
        alg = self.baseline_cost_algorithms[compressor_type]
        acfm_lb, acfm_ub = alg.acfm_bounds
        acfm = self.outs[0].get_total_flow('cfm')
        design_results['Compressors in parallel'] = ceil(acfm / acfm_ub) if acfm > acfm_ub else 1
        design_results['Driver'] = alg.driver if self._driver == 'Default' else self._driver 
    
    def _cost(self):
        design_results = self.design_results
        alg = self.baseline_cost_algorithms[design_results['Type']]
        acfm_lb, acfm_ub = alg.acfm_bounds
        Pc = self.power_utility.get_property('consumption', 'hp')
        N = design_results['Compressors in parallel']
        F = Pc / N
        self.baseline_purchase_costs['Compressor(s)'] = N * bst.CE / alg.CE * alg.cost(F)
        self.F_D['Compressor(s)'] = self._F_D[design_results['Driver']]


class IsothermalCompressor(Compressor):
    """
    Create an isothermal compressor.

    Parameters
    ----------
    ins : stream
        Inlet fluid.
    outs : stream
        Outlet fluid.
    P : float
        Outlet pressure [Pa].
    eta : float
        Isothermal efficiency.
    vle : bool
        Whether to perform phase equilibrium calculations on
        the outflow. If False, the outlet will be assumed to be the same
        phase as the inlet.

    Notes
    -----
    Default compressor selection, design and cost algorithms are adapted from [0]_.

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
    >>> feed = bst.Stream(H2=1, T=298.15, P=20e5, phase='g')
    >>> K = bst.units.IsothermalCompressor('K', ins=feed, P=350e5, eta=1)
    >>> K.simulate()
    >>> K.show()
    IsothermalCompressor: K1
    ins...
    [0] s1
        phase: 'g', T: 298.15 K, P: 2e+06 Pa
        flow (kmol/hr): H2  1
    outs...
    [0] s2
        phase: 'g', T: 298.15 K, P: 3.5e+07 Pa
        flow (kmol/hr): H2  1
    
    >>> K.results()
    Isothermal compressor                          Units              K1
    Power               Rate                          kW             2.1
                        Cost                      USD/hr           0.164
    Chilled water       Duty                       kJ/hr       -7.26e+03
                        Flow                     kmol/hr            7.53
                        Cost                      USD/hr          0.0363
    Design              Type                              Recriprocating
                        Compressors in parallel                        1
                        Driver                            Electric motor
                        Ideal power                   kW             2.1
                        Ideal duty                 kJ/hr       -7.26e+03
    Purchase cost       Compressor(s)                USD             385
    Total purchase cost                              USD             385
    Utility cost                                  USD/hr           0.201

    References
    ----------
    .. [0] Sinnott, R. and Towler, G (2019). "Chemical Engineering Design: SI Edition (Chemical Engineering Series)". 6th Edition. Butterworth-Heinemann.

    """

    def _run(self):
        feed = self.ins[0]
        out = self.outs[0]
        out.copy_like(feed)
        out.P = self.P
        out.T = feed.T
        if self.vle is True: out.vle(T=out.T, P=out.P)

    def _design(self):
        feed = self.ins[0]
        outlet = self.outs[0]
        ideal_power, ideal_duty = self._calculate_ideal_power_and_duty()
        self.power_utility.consumption = ideal_power / self.eta
        Q = ideal_duty / self.eta
        super()._design()
        self.heat_utilities[0](unit_duty=Q, T_in=feed.T, T_out=outlet.T)
        self.design_results['Ideal power'] = ideal_power # kW
        self.design_results['Ideal duty'] = ideal_duty # kJ / hr


class IsentropicCompressor(Compressor):
    """
    Create an isentropic compressor.

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

    Notes
    -----
    Default compressor selection, design and cost algorithms are adapted from [0]_.

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
        phase: 'g', T: 1151.3 K, P: 5e+06 Pa
        flow (kmol/hr): H2  1

    >>> K.results()
    Isentropic compressor                         Units              K1
    Power               Rate                         kW            7.03
                        Cost                     USD/hr            0.55
    Design              Ideal power                  kW            4.92
                        Ideal duty                kJ/hr               0
                        Type                             Recriprocating
                        Compressors in parallel                       1
                        Driver                           Electric motor
    Purchase cost       Compressor(s)               USD         1.7e+03
    Total purchase cost                             USD         1.7e+03
    Utility cost                                 USD/hr            0.55


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
        phases: ('g', 'l'), T: 797.75 K, P: 1e+07 Pa
        flow (kmol/hr): (g) H2O  1

    >>> K.results()
    Isentropic compressor                         Units              K2
    Power               Rate                         kW            5.41
                        Cost                     USD/hr           0.423
    Design              Ideal power                  kW            5.41
                        Ideal duty                kJ/hr               0
                        Type                             Recriprocating
                        Compressors in parallel                       1
                        Driver                           Electric motor
    Purchase cost       Compressor(s)               USD        1.23e+03
    Total purchase cost                             USD        1.23e+03
    Utility cost                                 USD/hr           0.423

    References
    ----------
    .. [0] Sinnott, R. and Towler, G (2019). "Chemical Engineering Design: SI Edition (Chemical Engineering Series)". 6th Edition. Butterworth-Heinemann.

    """
    _N_heat_utilities = 0

    def _run(self):
        feed = self.ins[0]
        out = self.outs[0]
        out.copy_like(feed)
        out.P = self.P
        out.S = feed.S
        if self.vle is True:
            out.vle(S=out.S, P=out.P)
            T_isentropic = out.T
        else:
            T_isentropic = out.T
        self.T_isentropic = T_isentropic
        dH_isentropic = out.H - feed.H
        self.design_results['Ideal power'] = ideal_power = dH_isentropic / 3600. # kW
        self.design_results['Ideal duty'] = 0.
        dH_actual = dH_isentropic / self.eta
        out.H = feed.H + dH_actual        
        if self.vle is True: out.vle(H=out.H, P=out.P)
        self.power_utility.consumption = ideal_power  / self.eta


class PolytropicCompressor(Compressor):
    """
    Create a polytropic compressor.

    Parameters
    ----------
    ins : stream
        Inlet fluid.
    outs : stream
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

    Notes
    -----
    Default compressor selection, design and cost algorithms are adapted from [0]_.

    Examples
    --------
    Simulate polytropic compression of hydrogen with 70% efficiency using the Schultz method [1]_:

    >>> import biosteam as bst
    >>> thermo = bst.Thermo([bst.Chemical('H2')])
    >>> thermo.mixture.include_excess_energies = True
    >>> bst.settings.set_thermo(thermo)
    >>> feed = bst.Stream('feed', H2=1, T=25 + 273.15, P=20e5, phase='g')
    >>> K = bst.units.PolytropicCompressor('K1', ins=feed, outs='outlet', P=350e5, eta=0.7, method='schultz')
    >>> K.simulate()
    >>> K.show()
    PolytropicCompressor: K1
    ins...
    [0] feed
        phase: 'g', T: 298.15 K, P: 2e+06 Pa
        flow (kmol/hr): H2  1
    outs...
    [0] outlet
        phase: 'g', T: 961.98 K, P: 3.5e+07 Pa
        flow (kmol/hr): H2  1
    >>> K.results()
    Polytropic compressor                         Units              K1
    Power               Rate                         kW            5.54
                        Cost                     USD/hr           0.433
    Design              Type                             Recriprocating
                        Compressors in parallel                       1
                        Driver                           Electric motor
    Purchase cost       Compressor(s)               USD        1.27e+03
    Total purchase cost                             USD        1.27e+03
    Utility cost                                 USD/hr           0.433


    Repeat using Hundseid method [2]_:

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
        phase: 'g', T: 958.12 K, P: 3.5e+07 Pa
        flow (kmol/hr): H2  1
    >>> K.results()
    Polytropic compressor                         Units              K1
    Power               Rate                         kW            5.51
                        Cost                     USD/hr           0.431
    Design              Type                             Recriprocating
                        Compressors in parallel                       1
                        Driver                           Electric motor
    Purchase cost       Compressor(s)               USD        1.26e+03
    Total purchase cost                             USD        1.26e+03
    Utility cost                                 USD/hr           0.431


    References
    ----------
    .. [0] Sinnott, R. and Towler, G. (2019). "Chemical Engineering Design: SI Edition (Chemical Engineering Series)". 6th Edition. Butterworth-Heinemann.
    .. [1] Schultz, J. (1962). "The Polytropic Analysis of Centrifugal Compressors". J. Eng. Power., 84(1): 69-82 (14 pages)
    .. [2] Hundseid, O., Bakken, L. E. and Helde, T. (2006). “A Revised Compressor Polytropic Performance Analysis,” Proceedings of ASME GT2006, Paper Number 91033, ASME Turbo Expo 2006.

    """
    _N_heat_utilities = 0
    available_methods = {'schultz', 'hundseid'}
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *, P, eta=0.7, 
                 vle=False, compressor_type=None, method=None, n_steps=None):
        Compressor.__init__(self, ID=ID, ins=ins, outs=outs, thermo=thermo, P=P, eta=eta, vle=vle, compressor_type=compressor_type)
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

        return W_actual / out.F_mol # kJ/hr

    def _run(self):
        feed = self.ins[0]
        out = self.outs[0]
        out.copy_like(feed)
        name = '_' + self._method
        method = getattr(self, name)
        W = method() # Polytropic work [kJ/hr]
        if self.vle is True: out.vle(H=out.H, P=out.P)
        self.power_utility.consumption = W / 3600 # kJ/hr -> kW