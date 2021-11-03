# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import numpy as np
from .. import Unit
from .design_tools import size_batch
from .decorators import cost
import flexsolve as flx
from math import ceil

__all__ = ('BatchBioreactor',)

@cost('Recirculation flow rate', 'Recirculation pumps', kW=30, S=77.22216,
      cost=47200, n=0.8, BM=2.3, CE=522, N='Number of reactors')
@cost('Reactor volume', 'Cleaning in place', CE=521.9,
      cost=421e3, S=3785, n=0.6, BM=1.8)
@cost('Reactor volume', 'Agitators', CE=521.9, cost=52500,
      S=3785, n=0.5, kW=22.371, BM=1.5, N='Number of reactors')
@cost('Reactor volume', 'Reactors', CE=521.9, cost=844000,
      S=3785, n=0.5, BM=1.5, N='Number of reactors')
@cost('Reactor duty', 'Heat exchangers', CE=522, cost=23900,
      S=20920000.0, n=0.7, BM=2.2, N='Number of reactors') # Based on a similar heat exchanger
class BatchBioreactor(Unit, isabstract=True):
    """
    Abstract Bioreactor class. Conversion is based on reaction time, `tau`.
    Cleaning and unloading time,`tau_0`, fraction of working volume, `V_wf`,
    and number of reactors, `N_reactors`, are attributes that can be changed.
    The cost of a reactor is based on the NREL batch fermentation tank cost 
    assuming volumetric scaling with a 6/10th exponent [1]_. 
    
    **Abstract methods**
        
    kinetic_model(z, t, *kinetic_constants): 
        A staticmethod that returns effluent concentrations (z_t; kg/m^3)
        given the starting concentration (z; kg/m^3), reaction time (t; hr),
        and kinetic constants.
    
    _run():
        Must set outlet stream flows.
        
    _cost():
        Must set purchase cost results.
    
    Parameters
    ----------
    ins : streams
        Inlet fluids to be mixed into the reactor.
    outs : stream sequence
        * [0] Vent
        * [1] Effluent
    tau : float
        Reaction time [hr].
    N : int, optional
        Number of batch reactors
    V : float, optional
        Target volume of reactors [m^3].
    T=305.15 : float
        Operating temperature of reactor [K].
    P=101325 : float
        Operating pressure of reactor [Pa].
    Nmin=2 : int
        Minimum number of fermentors.
    Nmax=36: int
        Maximum number of fermentors.  
    
    Notes
    -----
    Either N or V must be given.
    
    References
    ----------
    .. [1] D. Humbird, R. Davis, L. Tao, C. Kinchin, D. Hsu, and A.
        Aden National. Renewable Energy Laboratory Golden, Colorado. P.
        Schoen, J. Lukas, B. Olthof, M. Worley, D. Sexton, and D. Dudgeon.
        Harris Group Inc. Seattle, Washington and Atlanta, Georgia. Process 
        Design and Economics for Biochemical Conversion of Lignocellulosic
        Biomass to Ethanol Dilute-Acid Pretreatment and Enzymatic Hydrolysis 
        of Corn Stover. May 2011. Technical Report NREL/TP-5100-47764
    
    """
    _units = {'Reactor volume': 'm3',
              'Cycle time': 'hr',
              'Batch time': 'hr',
              'Loading time': 'hr',
              'Total dead time': 'hr',
              'Reactor duty': 'kJ/hr',
              'Recirculation flow rate': 'm3/hr'}
    _N_ins = _N_outs = 2
    _N_heat_utilities = 1
    
    #: [bool] If True, number of reactors (N) is chosen as to minimize installation cost in every simulation. Otherwise, N remains constant.
    autoselect_N = False
    
    #: [float] Cleaning and unloading time (hr).
    tau_0 = 3
    
    #: [float] Fraction of filled tank to total tank volume.
    V_wf = 0.9
    
    def _get_design_info(self):
        return (('Cleaning and unloading time', self.tau_0, 'hr'),
                ('Working volume fraction', self.V_wf, ''))
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 tau=None, N=None, V=None, T=305.15, P=101325,
                 Nmin=2, Nmax=36):
        Unit.__init__(self, ID, ins, outs, thermo)
        self._N = N; self._V = V
        
        #: [float] Reaction time [hr].
        self.tau = tau
        
        #: [int] Number of batch reactors
        if N: self.N = N
        
        #: [float] Target volume of a fermentor
        if V: self.V = V
        
        #: [float] Operating temperature of reactor [K].
        self.T = T
        
        #: [float] Operating pressure of reactor [Pa].
        self.P = P
        
        #: [int] Minimum number of fermentors
        self.Nmin = Nmin
        
        #: [int] Maximum number of fermentors
        self.Nmax = Nmax
        
    def _setup(self):
        super()._setup()
        vent, effluent = self.outs
        vent.phase = 'g'
        vent.T = effluent.T = self.T
        vent.P = effluent.P = self.P
        
    @property
    def vent(self):
        return self.outs[0]
    @property
    def effluent(self):
        return self.outs[1]
        
    @property
    def N(self):
        """[int] Number of reactor."""
        return self._N
    @N.setter
    def N(self, N):
        if N <= 1:
            raise ValueError(f"number of reactors must be greater than 1, value {N} is infeasible")
        assert not self._V, 'cannot specify both reactor volume and number of reactors'
        self._N = ceil(N)

    @property
    def V(self):
        """[float] Reactor volume."""
        return self._V
    @V.setter
    def V(self, V):
        if V <= 1:
            raise ValueError(f"reactor volume must be greater than 1, value {V} is infeasible")
        assert not self._N, 'cannot specify both reactor volume and number of reactors'
        self._V = V

    @property
    def tau(self):
        return self._tau
    @tau.setter
    def tau(self, tau):
        self._tau = tau
    
    @property
    def N_at_minimum_capital_cost(self):
        cost_old = np.inf
        self.autoselect_N = False
        self._N, N = 2, self._N
        cost_new = self.purchase_cost
        self._summary()
        while cost_new < cost_old:
            self._N += 1
            self._summary()
            cost_old = cost_new
            cost_new = self.purchase_cost
        self._N, N = N, self._N
        self.autoselect_N = True
        return N - 1
        
    def _design(self):
        effluent = self.effluent
        v_0 = effluent.F_vol
        tau = self._tau
        tau_0 = self.tau_0
        V_wf = self.V_wf
        Design = self.design_results
        if self.autoselect_N:
            N = self.N_at_minimum_capital_cost
        elif self.V:
            f = lambda N: v_0 / N / V_wf * (tau + tau_0) / (1 - 1 / N) - self.V
            if f(self.Nmax) > 0.:
                N = self.Nmax
            elif f(self.Nmin) < 0.:
                N = self.Nmin
            else:
                N = flx.IQ_interpolation(f, self.Nmin, self.Nmax,
                                         xtol=0.01, ytol=0.5, checkbounds=False)
                N = ceil(N)
        else:
            N = self._N
        Design.update(size_batch(v_0, tau, tau_0, N, V_wf))
        Design['Number of reactors'] = N
        Design['Recirculation flow rate'] = v_0 / N
        duty = self.Hnet
        Design['Reactor duty'] = abs(duty)
        self.heat_utilities[0](duty, self.T)
        
    