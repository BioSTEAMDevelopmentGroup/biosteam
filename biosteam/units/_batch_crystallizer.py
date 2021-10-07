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

__all__ = ('BatchCrystallizer',)

crystallizer_material_factors = {
    'Carbon steel': 1.0,
    'Rubber lined carbon steel': 1.48,
    'Stainless steel': 3.43,
}

@cost('Crystallizer volume', 'Crystallizer',
      CE=444., S=0.003785411784, # Originally 1 gal
      BM=2.0, N='Number of crystallizers',
      f=lambda S: 222.4 * S**0.71 + 35150)
class BatchCrystallizer(Unit):
    """
    Abstract unit operation for atmospheric batch crystallization. Subclasses
    must implement the mass and energy balance through the `_run` method.
    
    Parameters
    ----------
    ins : stream
        Inlet.
    outs : stream
        Effluent.
    tau : float
        Reaction time [hr].
    N : int, optional
        Number of vessels. Either N or V must be given.
    V : float, optional
        Target volume of vessels [m^3].
    T=305.15 : float
        Operating temperature [K].
    Nmin=2 : int
        Minimum number of vessels.
    Nmax=36: int
        Maximum number of vessels.  
    kW : float, optional
        Electricity usage per volume in kW/gal. Defaults to 0.00746, a 
        heuristic value for suspension of solids.
    
    Notes
    -----
    The cost correlation for the crystallizer taken directly from the 
    Matches engineering firm website, which provides educational cost
    correlation tools for a variety of equipment [1]_. The default
    value of kW was retrieved from  [2]_.
    
    References
    ----------
    .. [1] Matches; certified engineering company at Edmond, Oklahoma.
    http://www.matche.com/equipcost/Crystallizer.html
    
    .. [2] Seider, W. D., Lewin,  D. R., Seader, J. D., Widagdo, S., Gani,
        R., & Ng, M. K. (2017). Product and Process Design Principles. Wiley.
        Cost Accounting and Capital Cost Estimation (Chapter 16)
    
    """
    _units = {'Crystallizer volume': 'm3',
              'Batch time': 'hr',
              'Loading time': 'hr'}
    
    _N_ins = _N_outs = 1
    _N_heat_utilities = 1
    
    #: [float] Cleaning and unloading time (hr).
    tau_0 = 1
    
    #: [float] Fraction of filled tank to total tank volume.
    V_wf = 0.9
    
    def _get_design_info(self):
        return (('Cleaning and unloading time', self.tau_0, 'hr'),
                ('Working volume fraction', self.V_wf, ''))
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 tau=None, N=None, V=None, T=305.15,
                 Nmin=2, Nmax=36, vessel_material='Carbon steel',
                 kW=0.00746):
        Unit.__init__(self, ID, ins, outs, thermo)
        self._N = N; self._V = V
        
        #: [float] Reaction time [hr].
        self.tau = tau
        
        #: [int] Number of vessels
        if N: self.N = N
        
        #: [float] Target volume of a vessels
        if V: self.V = V
        
        #: [float] Operating temperature of reactor [K].
        self.T = T
        
        #: [int] Minimum number of vessels
        self.Nmin = Nmin
        
        #: [int] Maximum number of vessels
        self.Nmax = Nmax
        
        #: [float] Electricity usage per volume in kW/gal
        self.kW = kW
        
    def _setup(self):
        super()._setup()
        effluent, = self.outs
        effluent.T = self.T
        effluent.P = 101325
    
    @property
    def effluent(self):
        return self.outs[0]
        
    @property
    def N(self):
        """[int] Number of reactor."""
        return self._N
    @N.setter
    def N(self, N):
        if N <= 1:
            raise ValueError(f"number of crystallizers must be greater than 1, value {N} is infeasible")
        assert not self._V, 'cannot specify both crystallizer volume and number of crystallizers'
        self._N = ceil(N)

    @property
    def V(self):
        """[float] Crystallizer volume."""
        return self._V
    @V.setter
    def V(self, V):
        if V <= 1:
            raise ValueError(f"crystallizer volume must be greater than 1, value {V} is infeasible")
        assert not self._N, 'cannot specify both crystallizer volume and number of crystallizers'
        self._V = V

    @property
    def tau(self):
        return self._tau
    @tau.setter
    def tau(self, tau):
        self._tau = tau
        
    def _design(self):
        effluent = self.effluent
        v_0 = effluent.F_vol
        tau = self._tau
        tau_0 = self.tau_0
        V_wf = self.V_wf
        Design = self.design_results
        if self.V:
            f = lambda N: v_0 / N / V_wf * (tau + tau_0) / (1 - 1 / N) - self.V
            N = flx.IQ_interpolation(f, self.Nmin, self.Nmax,
                                     xtol=0.01, ytol=0.5, checkbounds=False)
            N = ceil(N)
        else:
            N = self._N
        dct = size_batch(v_0, tau, tau_0, N, V_wf)
        Design['Crystallizer volume'] = volume = dct.pop('Reactor volume')
        Design.update(dct)
        Design['Number of crystallizers'] = N
        self.heat_utilities[0](self.Hnet, self.T)
        self.power_utility.consumption = self.kW * V_wf * volume * N
        
    