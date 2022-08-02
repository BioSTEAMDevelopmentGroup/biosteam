# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
import biosteam as bst
from .. import Unit
from biosteam.units.heat_exchange import HX, HXutility
import warnings
from numpy import log

__all__ = ('IsentropicCompressor', 'IsothermalCompressor', 'PolytropicCompressor', 'MultistageCompressor')


#: TODO:
#: * Implement estimate of isentropic efficiency when not given.
#: * Add option to default efficiencies to heuristic values for each type of compressor.
#: * Implement bare-module, design, and material factors.
#: * Maybe use cost correlations from Warren's Process Development and Design for
#:   consistency with factors.
#: * Move cost coefficients to a dictionary.
#: * Only calculate volumetric flow rate if type is Blower.
#: * Only calculate power if type is not blower.
class _CompressorBase(Unit):
    """
    Abstract base class for all compressor types.
    """
    _N_ins = 1
    _N_outs = 1
    _N_heat_utilities = 1
    _F_BM_default = {'Compressor': 1.0}
    _units = {
        'Type': '-',
        'Power': 'kW',
        'Duty': 'kJ/kmol',
        'Outlet Temperature': 'K',
        'Volumetric Flow Rate': 'm^3/hr',
        'Ideal Power': 'kW',
        'Ideal Duty': 'kJ/kmol',
        'Ideal Outlet Temperature': 'K',
    }

    def __init__(self, ID='', ins=None, outs=(), thermo=None, *, P, eta=0.7, vle=False, type=None):
        super().__init__(ID=ID, ins=ins, outs=outs, thermo=thermo)
        self.P = P  #: Outlet pressure [Pa].
        self.eta = eta  #: Isentropic efficiency.

        #: Whether to perform phase equilibrium calculations on the outflow.
        #: If False, the outlet will be assumed to be the same phase as the inlet.
        self.vle = vle

        #: Type of compressor : blower/centrifugal/reciprocating.
        #: If None, the type will be determined automatically.
        self.type = type

        # make sure user-given types are not overwritten
        if type is None:
            self._overwrite_type = True
        else:
            self._overwrite_type = False

    def _setup(self):
        super()._setup()

    def _determine_compressor_type(self):
        # Determine compressor type based on power specification

        # don't overwrite user input
        if not self._overwrite_type:
            return self.type

        # determine type based on power
        power = self.power
        if 0 <= power < 93:
            self.type = 'Blower'
        elif 93 <= power < 16800:
            self.type = 'Reciprocating'
        elif 16800 <= power <= 30000:
            self.type = 'Centrifugal'
        else:
            raise RuntimeError(
                f"power requirement ({power / 1e3:.3g} MW) is outside cost "
                "correlation range (0, 30 MW). No fallback for this case has "
                "been implemented yet"
            )
        return self.type

    def _calculate_ideal_power(self):
        feed = self.ins[0]
        out = self.outs[0]
        dH = out.H - feed.H
        self.Q = TdS = feed.T * (out.S - feed.S)  # kJ/hr
        self.power = (dH - TdS) / 3600  # kW

    def _run(self):
        super()._run()

    def _design(self):
        feed = self.ins[0]
        out = self.outs[0]

        # set power utility
        power = self.power
        self.power_utility(power)

        # determine compressor type depending on power rating
        type = self._determine_compressor_type()

        # write design parameters
        self.design_results["Type"] = type
        self.design_results['Power'] = power
        self.design_results['Duty'] = self.Q
        self.design_results['Outlet Temperature'] = out.T
        self.design_results['Volumetric Flow Rate'] = feed.F_vol

    def _cost(self):
        # cost calculation adapted from Sinnott & Towler: Chemical Engineering Design, 6th Edition, 2019, p.296-297
        # all costs on U.S. Gulf Coast basis, Jan. 2007 (CEPCI = 509.7)
        cost = self.baseline_purchase_costs
        flow_rate = self.design_results['Volumetric Flow Rate']
        power = self.design_results['Power']
        if self.type == "Blower":
            a = 3800
            b = 49
            n = 0.8
            S = flow_rate
        elif self.type == "Reciprocating":
            a = 220000
            b = 2300
            n = 0.75
            S = power
        elif self.type == "Centrifugal":
            a = 490000
            b = 16800
            n = 0.6
            S = power
        else:
            a = b = n = S = 0
        cost["Compressor"] = bst.CE / 509.7 * (a + b * S ** n)


class IsothermalCompressor(_CompressorBase):
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
    type: str
        Type of compressor : blower/centrifugal/reciprocating. If None, the type
        will be determined automatically.

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
    >>> K = bst.units.IsothermalCompressor(ins=feed, P=350e5, eta=1)
    >>> K.simulate()
    >>> K.results()
    Isothermal compressor                           Units        K3
    Power               Rate                           kW       2.1
                        Cost                       USD/hr     0.164
    Chilled water       Duty                        kJ/hr -7.26e+03
                        Flow                      kmol/hr      7.53
                        Cost                       USD/hr    0.0363
    Design              Type                            -    Blower
                        Power                          kW       2.1
                        Duty                      kJ/kmol -7.26e+03
                        Outlet Temperature              K       298
                        Volumetric Flow Rate       m^3/hr      1.24
                        Ideal Power                    kW       2.1
                        Ideal Duty                kJ/kmol -7.26e+03
                        Ideal Outlet Temperature        K       298
    Purchase cost       Compressor                    USD   4.3e+03
    Total purchase cost                               USD   4.3e+03
    Utility cost                                   USD/hr     0.201

    References
    ----------
    .. [0] Sinnott, R. and Towler, G (2019). "Chemical Engineering Design: SI Edition (Chemical Engineering Series)". 6th Edition. Butterworth-Heinemann.

    """

    def _run(self):
        feed = self.ins[0]
        out = self.outs[0]
        out.copy_like(feed)

        # calculate isothermal state change
        out.P = self.P
        out.T = feed.T

        # check phase equilibirum
        if self.vle is True:
            out.vle(T=out.T, P=out.P)

        # calculate ideal power demand and duty
        self._calculate_ideal_power()
        self.ideal_power = self.power
        self.ideal_duty = self.Q

        # calculate actual power and duty (incl. efficiency)
        self.power = self.power / self.eta
        self.Q = self.Q / self.eta

    def _design(self):
        # set default design parameters
        super()._design()

        # set heat utility
        feed = self.ins[0]
        out = self.outs[0]
        u = bst.HeatUtility(heat_transfer_efficiency=1, heat_exchanger=None)
        u(unit_duty=self.Q, T_in=feed.T, T_out=out.T)
        self.heat_utilities = (u, bst.HeatUtility(), bst.HeatUtility())

        # save other design parameters
        F_mol = self.outs[0].F_mol
        self.design_results['Ideal Power'] = self.ideal_power # kW
        self.design_results['Ideal Duty'] = self.ideal_duty / feed.F_mol # kJ/hr -> kJ/kmol
        self.design_results['Ideal Outlet Temperature'] = feed.T  # K


class IsentropicCompressor(_CompressorBase):
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
    type: str
        Type of compressor : blower/centrifugal/reciprocating. If None, the type
        will be determined automatically.

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
    Isentropic compressor                           Units       K1
    Power               Rate                           kW     7.03
                        Cost                       USD/hr     0.55
    Design              Type                            -   Blower
                        Power                          kW     7.03
                        Duty                      kJ/kmol 9.63e-09
                        Outlet Temperature              K 1.15e+03
                        Volumetric Flow Rate       m^3/hr     24.5
                        Ideal Power                    kW     4.92
                        Ideal Duty                kJ/kmol        0
                        Ideal Outlet Temperature        K      901
    Purchase cost       Compressor                    USD 4.94e+03
    Total purchase cost                               USD 4.94e+03
    Utility cost                                   USD/hr     0.55


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
    Isentropic compressor                           Units       K2
    Power               Rate                           kW     5.41
                        Cost                       USD/hr    0.423
    Design              Type                            -   Blower
                        Power                          kW     5.41
                        Duty                      kJ/kmol 6.67e-07
                        Outlet Temperature              K      798
                        Volumetric Flow Rate       m^3/hr     27.9
                        Ideal Power                    kW     5.41
                        Ideal Duty                kJ/kmol        0
                        Ideal Outlet Temperature        K      798
    Purchase cost       Compressor                    USD 5.01e+03
    Total purchase cost                               USD 5.01e+03
    Utility cost                                   USD/hr    0.423

    References
    ----------
    .. [0] Sinnott, R. and Towler, G (2019). "Chemical Engineering Design: SI Edition (Chemical Engineering Series)". 6th Edition. Butterworth-Heinemann.

    """
    _N_heat_utilities = 0

    def _run(self):
        feed = self.ins[0]
        out = self.outs[0]
        out.copy_like(feed)

        # calculate isentropic outlet state
        out.P = self.P
        out.S = feed.S
        T_isentropic = out.T

        # check phase equilibirum
        if self.vle is True:
            out.vle(S=out.S, P=out.P)
            T_isentropic = out.T

        # calculate ideal power demand
        self._calculate_ideal_power()

        # calculate actual state change (incl. efficiency)
        dh_isentropic = out.h - feed.h
        dh_actual = dh_isentropic / self.eta
        out.h = feed.h + dh_actual

        # check phase equilibirum again
        if self.vle is True:
            out.vle(H=out.H, P=out.P)

        # save values for _design
        self.power = self.power / self.eta
        self.T_isentropic = T_isentropic
        self.dh_isentropic = dh_isentropic

    def _design(self):
        # set default design parameters
        super()._design()

        # set isentropic compressor specific design parameters
        F_mol = self.outs[0].F_mol
        self.design_results['Ideal Power'] = (self.dh_isentropic * F_mol) / 3600  # kJ/kmol * kmol/hr / 3600 s/hr -> kW
        self.design_results['Ideal Duty'] = 0 # kJ/kmol
        self.design_results['Ideal Outlet Temperature'] = self.T_isentropic # K


class PolytropicCompressor(_CompressorBase):
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
        for real gases at high ressure ratios.
    n_steps: int
        Number of virtual steps used in numerical integration for hundseid method.
    type: str
        Type of compressor : blower/centrifugal/reciprocating. If None, the type
        will be determined automatically.

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
    Polytropic compressor                       Units      K1
    Power               Rate                       kW    5.54
                        Cost                   USD/hr   0.433
    Design              Type                        -  Blower
                        Power                      kW    5.54
                        Duty                  kJ/kmol       0
                        Outlet Temperature          K     962
                        Volumetric Flow Rate   m^3/hr    1.24
    Purchase cost       Compressor                USD 4.3e+03
    Total purchase cost                           USD 4.3e+03
    Utility cost                               USD/hr   0.433


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
    Polytropic compressor                       Units      K1
    Power               Rate                       kW    5.51
                        Cost                   USD/hr   0.431
    Design              Type                        -  Blower
                        Power                      kW    5.51
                        Duty                  kJ/kmol       0
                        Outlet Temperature          K     958
                        Volumetric Flow Rate   m^3/hr    1.24
    Purchase cost       Compressor                USD 4.3e+03
    Total purchase cost                           USD 4.3e+03
    Utility cost                               USD/hr   0.431


    References
    ----------
    .. [0] Sinnott, R. and Towler, G. (2019). "Chemical Engineering Design: SI Edition (Chemical Engineering Series)". 6th Edition. Butterworth-Heinemann.
    .. [1] Schultz, J. (1962). "The Polytropic Analysis of Centrifugal Compressors". J. Eng. Power., 84(1): 69-82 (14 pages)
    .. [2] Hundseid, O., Bakken, L. E. and Helde, T. (2006). “A Revised Compressor Polytropic Performance Analysis,” Proceedings of ASME GT2006, Paper Number 91033, ASME Turbo Expo 2006.

    """

    _N_heat_utilities = 0

    def __init__(self, ID='', ins=None, outs=(), thermo=None, *, P, eta=0.7, vle=False, type=None, method="schultz",
                 n_steps=100):
        super().__init__(ID=ID, ins=ins, outs=outs, thermo=thermo, P=P, eta=eta, vle=vle, type=type)

        if method == "schultz":
            self._calculate = self._schultz_method
        elif method == "hundseid":
            self._calculate = self._hundseid_method
            self.n_steps = n_steps
        else:
            raise RuntimeError(f"Unrecognized method f{method} for PolytropicCompressor "
                               f"f{self.ID}. Must be 'schultz' or 'hundseid'.")

    def _schultz_method(self):
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

    def _hundseid_method(self):
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

        # calculate polytropic work and outlet state
        W = self._calculate() # kJ/hr

        # check phase equilibirum
        if self.vle is True:
            out.vle(H=out.H, P=out.P)

        # save values for _design
        self.Q = 0
        self.power = W / 3600 # kJ/hr -> kW


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
    ins : stream
        Inlet fluid.
    outs : stream
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
    type: str
        (setup option 1) Type of compressor : blower/centrifugal/reciprocating. If None, the type
        will be determined automatically.
    compressors: list[_CompressorBase]
        (setup option 2) List of compressors to use for each stage.
    hxs: list[HX]
        (setup option 2) List of heat exchangers to use for each stage.

    Notes
    -----
    Default compressor selection, design and cost algorithms are adapted from [0]_.

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
    >>> K = bst.units.MultistageCompressor(ins=feed, pr=2, n_stages=4, eta=0.7)
    >>> K.simulate()
    >>> K.show()
    MultistageCompressor: K4
    ins...
    [0] feed
        phase: 'g', T: 298.15 K, P: 2e+06 Pa
        flow (kmol/hr): H2  1
    outs...
    [0] s11
        phase: 'g', T: 298.15 K, P: 3.2e+07 Pa
        flow (kmol/hr): H2  1
    >>> K.results()
    Multistage compressor                           Units                     K4
    Power               Rate                           kW                   3.14
                        Cost                       USD/hr                  0.246
    Chilled water       Duty                        kJ/hr              -1.12e+04
                        Flow                      kmol/hr                   7.41
                        Cost                       USD/hr                 0.0559
    Design              Type                            -  Multistage compressor
                        Power                          kW                   3.14
                        Duty                      kJ/kmol              -1.12e+04
                        Area                         ft^2                   1.54
                        Tube side pressure drop       psi                     12
                        Shell side pressure drop      psi                     20
                        Outlet Temperature              K                    298
                        Volumetric Flow Rate       m^3/hr                   1.24
    Purchase cost       K5 - Compressor               USD                4.3e+03
                        H1 - Double pipe              USD                    568
                        K6 - Compressor               USD               4.27e+03
                        H2 - Double pipe              USD                    675
                        K7 - Compressor               USD               4.25e+03
                        H3 - Double pipe              USD                    956
                        K8 - Compressor               USD               4.24e+03
                        H4 - Double pipe              USD               1.78e+03
    Total purchase cost                               USD                2.1e+04
    Utility cost                                   USD/hr                  0.301

    Show the fluid state at the outlet of each heat exchanger:
    >>> for hx in K.hxs:
    ...  hx.outs[0].show()
    ...
    Stream: s5 from <HXutility: H1> to <IsentropicCompressor: K6>
     phase: 'g', T: 298.15 K, P: 4e+06 Pa
     flow (kmol/hr): H2  1
    Stream: s7 from <HXutility: H2> to <IsentropicCompressor: K7>
     phase: 'g', T: 298.15 K, P: 8e+06 Pa
     flow (kmol/hr): H2  1
    Stream: s9 from <HXutility: H3> to <IsentropicCompressor: K8>
     phase: 'g', T: 298.15 K, P: 1.6e+07 Pa
     flow (kmol/hr): H2  1
    Stream: s11 from <HXutility: H4>
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
    >>> K = bst.units.MultistageCompressor(ins=feed, compressors=ks, hxs=hxs)
    >>> K.simulate()
    >>> K.show()
    MultistageCompressor: K14
    ins...
    [0] feed
        phase: 'g', T: 298.15 K, P: 2e+06 Pa
        flow (kmol/hr): H2  1
    outs...
    [0] s21
        phase: 'g', T: 298 K, P: 3.2e+07 Pa
        flow (kmol/hr): H2  1
    >>> K.results()
    Multistage compressor                           Units                    K14
    Power               Rate                           kW                   3.57
                        Cost                       USD/hr                  0.279
    Chilled water       Duty                        kJ/hr              -5.63e+03
                        Flow                      kmol/hr                   3.73
                        Cost                       USD/hr                 0.0282
    Cooling water       Duty                        kJ/hr              -7.12e+03
                        Flow                      kmol/hr                   4.87
                        Cost                       USD/hr                0.00237
    Design              Type                            -  Multistage compressor
                        Power                          kW                   3.57
                        Duty                      kJ/kmol              -1.28e+04
                        Area                         ft^2                   1.14
                        Tube side pressure drop       psi                     15
                        Shell side pressure drop      psi                     25
                        Outlet Temperature              K                    298
                        Volumetric Flow Rate       m^3/hr                   1.24
    Purchase cost       K9 - Compressor               USD                4.3e+03
                        H5 - Double pipe              USD                    298
                        K10 - Compressor              USD               4.28e+03
                        H6 - Double pipe              USD                    198
                        K11 - Compressor              USD               4.27e+03
                        H7 - Double pipe              USD                    129
                        K12 - Compressor              USD               4.26e+03
                        H8 - Double pipe              USD                    752
                        K13 - Compressor              USD               4.24e+03
                        H9 - Double pipe              USD               2.03e+03
    Total purchase cost                               USD               2.47e+04
    Utility cost                                   USD/hr                   0.31


    References
    ----------
    .. [0] Sinnott, R. and Towler, G. (2019). "Chemical Engineering Design: SI Edition (Chemical Engineering Series)". 6th Edition. Butterworth-Heinemann.

    """

    _N_ins = 1
    _N_outs = 1
    _N_heat_utilities = 0
    _units = {
        **_CompressorBase._units,
        **HX._units,
    }

    def __init__(
            self, ID='', ins=None, outs=(), thermo=None, *,
            pr=None, n_stages=None, eta=0.7, vle=False, type=None,
            compressors=None, hxs=None,
    ):
        super().__init__(ID=ID, ins=ins, outs=outs, thermo=thermo)

        # setup option 1: list of compressors and list of heat exchangers
        if compressors is not None and hxs is not None:
            if not isinstance(compressors[0], _CompressorBase):
                raise RuntimeError(f"Invalid parameterization of {self.ID}: `compressors` must "
                                   f"be a list of compressor objects.")
            elif not isinstance(hxs[0], HX):
                raise RuntimeError(f"Invalid parameterization of {self.ID}: `hxd` must "
                                   f"be a list of heat exchanger objects.")
            elif len(compressors) != len(hxs):
                raise RuntimeError(f"Invalid parameterization of {self.ID}: `compressors` and `hxs` "
                                   f"must have the same length.")
            else:
                self.compressors = compressors
                self.hxs = hxs
                self.pr = None
                self.n_stages = None

        # setup option 2: fixed pressure ratio and number of stages
        elif pr is not None and n_stages is not None:
            self.pr = pr
            self.n_stages = n_stages
            self.eta = eta
            self.vle=vle
            self.type=type
            self.compressors = None
            self.hxs = None
        else:
            raise RuntimeError(f"Invalid parameterization of {self.ID}: Must specify `pr` and "
                               f"`n_stages` or `compressors` and `hxs`.")

    def _setup(self):
        super()._setup()

        # helper variables
        feed = self.ins[0]

        # setup option 1: create connections between compressors and hxs
        if self.compressors is not None and self.hxs is not None:
            for n, (c, hx) in enumerate(zip(self.compressors, self.hxs)):
                if n == 0:
                    inflow = feed
                else:
                    inflow = self.hxs[n-1].outs[0]
                c.ins[0] = inflow
                hx.ins[0] = c.outs[0]

        # setup option 2: create connected compressor and hx objects
        elif self.pr is not None and self.n_stages is not None:
            self.compressors = []
            self.hxs = []
            for n in range(self.n_stages):
                if n==0:
                    inflow = feed
                    P = feed.P * self.pr
                else:
                    inflow = hx.outs[0]
                    P = P * self.pr

                c = IsentropicCompressor(
                    ins=inflow, P=P, eta=self.eta, vle=self.vle,
                    type=self.type
                )
                hx = HXutility(
                    ins=c.outs[0], T=inflow.T, rigorous=self.vle
                )
                self.compressors.append(c)
                self.hxs.append(hx)

        # set inlet and outlet reference
        self._ins = self.compressors[0].ins
        self._outs = self.hxs[-1].outs

        # set auxillary units
        units = [u for t in zip(self.compressors,self.hxs) for u in t]
        self.auxiliary_unit_names = tuple([u.ID for u in units])
        for u in units:
            self.__setattr__(u.ID, u)

    def _run(self):
        # calculate results

        # helper variables
        units = [u for t in zip(self.compressors, self.hxs) for u in t]

        # simulate all subcomponents
        for u in units:
            u.simulate()

        # set power utility = sum of power utilities
        self.power_utility = sum([u.power_utility for u in units])

        # set heat utility = sum of heat utilities
        # split by agent
        self.heat_utilities = sum([h for u in units for h in u.heat_utilities])


    def _design(self):
        self.design_results["Type"] = "Multistage compressor"

        # sum up design values
        units = [u for t in zip(self.compressors,self.hxs) for u in t]
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

        # add heat exchanger duties
        self.design_results["Duty"] += sum([-hx.Q for hx in self.hxs])

        # manual additions
        self.design_results["Outlet Temperature"] = self.outs[0].T
        self.design_results["Volumetric Flow Rate"] = self.ins[0].F_vol