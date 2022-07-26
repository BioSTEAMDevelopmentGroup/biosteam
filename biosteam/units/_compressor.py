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

__all__ = ('IsentropicCompressor',)

#: TODO: 
#: * Implement estimate of isentropic efficiency when not given.
#: * Add option to default efficiencies to heuristic values for each type of compressor. 
#: * Implement bare-module, design, and material factors.
#: * Maybe use cost correlations from Warren's Process Development and Design for 
#:   consistency with factors.
#: * Move cost coefficients to a dictionary.
#: * Allow user to enforce a compressor type.
#: * Only calculate volumetric flow rate if type is Blower.
#: * Only calculate power if type is not blower.
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
        self.eta = eta  #: Isentropic efficiency.
        
        #: Whether to perform phase equilibrium calculations on the outflow.
        #: If False, the outlet will be assumed to be the same phase as the inlet.
        self.vle = vle  
        
        #: Type of compressor (determined during cost calculation):
        #: blower/centrifugal/reciprocating
        self.type = None  

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
        self.power_utility(power)
        self.design_results['Isentropic Power'] = (self.dh_isentropic * out.F_mol) / 3600  # kJ/kmol * kmol/hr / 3600 s/hr -> kW
        self.design_results['Outlet Temperature'] = out.T
        self.design_results['Isentropic Outlet Temperature'] = self.T_isentropic
        self.design_results['Volumetric Flow Rate'] = feed.F_vol

        # Determine compressor type based on power specification
        if 0 <=  power < 93:
            self.type = 'Blower'
        elif 93 <= power < 16800:
            self.type = 'Reciprocating'
        elif 16800 <= power <= 30000:
            self.type = 'Centrifugal'
        else:
            raise RuntimeError(
                f"power requirement ({power / 1e3:.3g} MW) is above is outside cost "
                 "correlation range (0, 30 MW). No fallback for this case has "
                 "been implemented yet"
            )

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
        cost["Compressor"] = bst.CE / 509.7 * (a + b*S**n)
        
