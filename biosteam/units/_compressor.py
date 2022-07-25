# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
from .. import Unit
from flexsolve.open_solvers import aitken_secant

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

    """
    _N_ins = 1
    _N_outs = 1
    _N_heat_utilities = 0
    _units = {
        'Power': 'kW',
        'Isentropic Power': 'kW',
        'Outlet Temperature': 'K',
        'Isentropic Outlet Temperature': 'K',
    }

    def __init__(self, ID='', ins=None, outs=(), thermo=None, *, P, eta, vle=False):
        Unit.__init__(self, ID, ins, outs, thermo)
        self.P = P #: Outlet pressure [Pa].
        self.eta = eta  #: Isentropic efficiency [-].
        self.vle = vle

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
        self.design_results['Power'] = (out.H-feed.H)/3600 # kW
        self.design_results['Isentropic Power'] = (self.dh_isentropic * out.F_mol) / 3600  # kJ/kmol * kmol/hr / 3600 s/hr -> kW
        self.design_results['Outlet Temperature'] = out.T
        self.design_results['Isentropic Outlet Temperature'] = self.T_isentropic

    def _cost(self):
        # Todo
        pass