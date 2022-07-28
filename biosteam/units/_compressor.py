# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
import biosteam as bst
from .. import Unit
import warnings
from numpy import log

__all__ = ('IsentropicCompressor', 'IsothermalCompressor')


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
        'Outlet Temperature': 'K',
        'Volumetric Flow Rate': 'm^3/hr',
    }

    def __init__(self, ID='', ins=None, outs=(), thermo=None, *, P, vle=False, type=None):
        Unit.__init__(self, ID=ID, ins=ins, outs=outs, thermo=thermo)
        self.P = P  #: Outlet pressure [Pa].

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

        # calculate power demand
        self._calculate_ideal_power()

    def _design(self):
        # set default design parameters
        super()._design()

        # set heat utility
        feed = self.ins[0]
        out = self.outs[0]
        u = bst.HeatUtility(heat_transfer_efficiency=1, heat_exchanger=None)
        u(unit_duty=self.Q, T_in=feed.T, T_out=out.T)
        self.heat_utilities = (u, bst.HeatUtility(), bst.HeatUtility())


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
    Isentropic compressor                               Units       K1
    Power               Rate                               kW     7.03
                        Cost                           USD/hr     0.55
    Design              Power                              kW     7.03
                        Isentropic Power                   kW     4.92
                        Outlet Temperature                  K 1.15e+03
                        Isentropic Outlet Temperature       K      901
                        Volumetric Flow Rate           m^3/hr     24.5
    Purchase cost       Compressor                        USD 4.94e+03
    Total purchase cost                                   USD 4.94e+03
    Utility cost                                       USD/hr     0.55

    Per default, the outlet pahse is assumed to be the same as the inlet phase. If phase changes are to be accounted for,
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
    Isentropic compressor                               Units       K2
    Power               Rate                               kW     5.41
                        Cost                           USD/hr    0.423
    Design              Power                              kW     5.41
                        Isentropic Power                   kW     5.41
                        Outlet Temperature                  K      798
                        Isentropic Outlet Temperature       K      798
                        Volumetric Flow Rate           m^3/hr     27.9
    Purchase cost       Compressor                        USD 5.01e+03
    Total purchase cost                                   USD 5.01e+03
    Utility cost                                       USD/hr    0.423


    References
    ----------
    .. [0] Sinnott, R. and Towler, G (2019). "Chemical Engineering Design: SI Edition (Chemical Engineering Series)". 6th Edition. Butterworth-Heinemann.

    """
    _N_heat_utilities = 0
    _units = {
        'Type': '-',
        'Power': 'kW',
        'Isentropic Power': 'kW',
        'Outlet Temperature': 'K',
        'Isentropic Outlet Temperature': 'K',
        'Volumetric Flow Rate': 'm^3/hr',
    }

    def __init__(self, ID='', ins=None, outs=(), thermo=None, *, P, eta=0.7, vle=False):
        super().__init__(ID=ID, ins=ins, outs=outs, thermo=thermo, P=P, vle=vle)
        self.eta = eta  #: Isentropic efficiency.

    def _run(self):
        feed = self.ins[0]
        out = self.outs[0]
        out.copy_like(feed)

        # calculate isentropic state change
        out.P = self.P
        out.S = feed.S
        T_isentropic = out.T

        # check phase equilibirum
        if self.vle is True:
            out.vle(S=out.S, P=out.P)
            T_isentropic = out.T

        # calculate power demand
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
        self.design_results['Isentropic Power'] = (self.dh_isentropic * F_mol) / 3600  # kJ/kmol * kmol/hr / 3600 s/hr -> kW
        self.design_results['Isentropic Outlet Temperature'] = self.T_isentropic