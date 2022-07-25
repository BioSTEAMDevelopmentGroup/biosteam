# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
from .. import Unit, CE
from flexsolve.open_solvers import aitken_secant
import warnings

__all__ = ('IsentropicCompressor',)

class IsentropicCompressor(Unit):
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
        Isentropic efficiency [-].
    vle : bool
        Whether to perform phase equilibrium calculations on
        the outflow or not. Out phase = in phase if False.  [-].

    Notes
    -----
    Default compressor selection, design and cost algorithms are adapted from [0]_.

    References
    ----------
    .. [0] Sinnott, R. and Towler, G (2019). "Chemical Engineering Design: SI Edition (Chemical Engineering Series)". 6th Edition. Butterworth-Heinemann.
    """
    _N_ins = 1
    _N_outs = 1
    _N_heat_utilities = 0
    _units = {
        'Power': 'kW',
        'Isentropic Power': 'kW',
        'Outlet Temperature': 'K',
        'Isentropic Outlet Temperature': 'K',
        'Volumetric Flow Rate': 'm^3/hr',
    }
    _F_BM_default = {'Compressor': 1.0}

    def __init__(self, ID='', ins=None, outs=(), thermo=None, *, P, eta, vle=False):
        Unit.__init__(self, ID, ins, outs, thermo)
        self.P = P  #: Outlet pressure [Pa].
        self.eta = eta  #: Isentropic efficiency [-].
        self.vle = vle  #: Whether to perform phase equilibrium calculations on the outflow or not. Out phase = in phase if False.  [-].
        self.type = None  #: Which type of compressor (determined during cost calculation): blower/centrifugal/reciprocating

    def _setup(self):
        super()._setup()

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

        # calculate actual state change (incl. efficiency)
        dh_isentropic = out.h - feed.h
        dh_actual = dh_isentropic / self.eta
        out.h = feed.h + dh_actual

        # check phase equilibirum again
        if self.vle is True:
            out.vle(H=out.H, P=out.P)

        # save values for _design
        self.T_isentropic = T_isentropic
        self.dh_isentropic = dh_isentropic

    def _design(self):
        feed = self.ins[0]
        out = self.outs[0]

        # set design parameters
        self.design_results['Power'] = power = (out.H - feed.H) / 3600  # kW
        self.design_results['Isentropic Power'] = (self.dh_isentropic * out.F_mol) / 3600  # kJ/kmol * kmol/hr / 3600 s/hr -> kW
        self.design_results['Outlet Temperature'] = out.T
        self.design_results['Isentropic Outlet Temperature'] = self.T_isentropic
        self.design_results['Volumetric Flow Rate'] = feed.F_vol

        # determine compressor type based on power specification
        if 0 <=  power < 93:
            self.type = 'Blower'
        elif 93 <= power < 16800:
            self.type = 'Reciprocating'
        elif 16800 <= power <= 30000:
            self.type = 'Centrifugal'
        else:
            warnings.warn(f"Compressor power {power} outside cost correlation range (0, 30 MW). "
                          f"Setting cost of compressor {self.ID} to zero.")

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
        cost["Compressor"] = CE / 509.7 * (a + b*S**n)
        pass
