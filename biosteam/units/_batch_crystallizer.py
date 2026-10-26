# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
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
    ins : 
        Inlet.
    outs : 
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
    
    Examples
    --------
    Batch cooling crystallization of sucrose massecuite after
    vacuum-pan crystallization:

    >>> import biosteam as bst
    >>> bst.settings.set_thermo(['Water', 'Sucrose'])
    >>> massecuite = bst.Stream(
    ...     'massecuite', Water=5000, Sucrose=18000,
    ...     T=338.15, units='kg/hr', # 23 t/hr at 65 degC
    ... )
    >>> # Saturated mother liquor at 65 degC (kg sucrose/kg water).
    >>> dissolved = 3.029 * massecuite.imass['Water']
    >>> massecuite.phases = ('l', 's')
    >>> massecuite.imass['l', 'Sucrose'] = dissolved
    >>> massecuite.imass['s', 'Sucrose'] = 18000 - dissolved
    >>> C1 = bst.BatchCrystallizer(
    ...     'C1', ins=massecuite, outs='cooled_massecuite',
    ...     tau=20, V=175, T=313.15, # hr, m3, 40 degC
    ... )
    >>> @C1.add_specification
    ... def crystallize_sucrose():
    ...     C1.run()
    ...     feed = C1.ins[0]
    ...     effluent = C1.outs[0]
    ...     # Saturated mother liquor at 40 degC (kg sucrose/kg water).
    ...     dissolved = 2.381 * feed.imass['Water']
    ...     effluent.imass['l', 'Sucrose'] = dissolved
    ...     effluent.imass['s', 'Sucrose'] = (
    ...         feed.imass['Sucrose'] - dissolved
    ...     )
    >>> C1.simulate()
    >>> C1.show()
    >>> C1.results()
    BatchCrystallizer: C1
    ins...
    [0] massecuite  
        phases: ('l', 's'), T: 338.15 K, P: 101325 Pa
        flow (kmol/hr): (l) Water    278
                            Sucrose  44.2
                        (s) Sucrose  8.34
    outs...
    [0] cooled_massecuite  
        phases: ('l', 's'), T: 313.15 K, P: 101325 Pa
        flow (kmol/hr): (l) Water    278
                            Sucrose  34.8
                        (s) Sucrose  17.8
    
    >>> C1.results()
    Batch crystallizer                                      Units        C1
    Electricity              Power                             kW      2.65
                             Cost                          USD/hr     0.207
    Cooling water            Duty                           kJ/hr -1.69e+06
                             Flow                         kmol/hr  8.06e+03
                             Cost                          USD/hr      3.93
    Design                   Crystallizer volume               m3       132
                             Number of crystallizers                      3
                             Batch time                        hr        21
                             Loading time                      hr         0
                             Cleaning and unloading time       hr         1
                             Working volume fraction                    0.9
    Purchase cost            Crystallizer (x3)                USD  1.56e+06
    Total purchase cost                                       USD  1.56e+06
    Installed equipment cost                                  USD  3.13e+06
    Utility cost                                           USD/hr      4.14
    
    References
    ----------
    .. [1] Matches; certified engineering company at Edmond, Oklahoma.
    http://www.matche.com/equipcost/Crystallizer.html
    
    .. [2] Seider, W. D., Lewin,  D. R., Seader, J. D., Widagdo, S., Gani,
        R., & Ng, M. K. (2017). Product and Process Design Principles. Wiley.
        Cost Accounting and Capital Cost Estimation (Chapter 16)
    
    .. [3] Sugar Process Tech. "Vertical Crystalliser Design
        Calculation for Sugar Massecuite Cooling."
        https://www.sugarprocesstech.com/vertical-crystalliser-design/

    .. [4] Sugar Process Tech. "Low Grade B and C Massecuite
        Treatment Process in Sugar Crystallization."
        https://www.sugarprocesstech.com/low-grade-massecuite/

    .. [5] Sugar Process Tech. "Crystallizers Application in
        Sugar Industry: Crystallizer Capacity Calculation."
        https://www.sugarprocesstech.com/crystallizer/

    .. [6] Fraga, G., and Broadfoot, R. (2024).
        "Guidelines for Cost-Effectively Improving Sugar
        Recovery from C Massecuite Cooling Crystallisers."
        Australian Society of Sugar Cane Technologists.
        https://assct.com.au/document/guidelines-for-cost-effectively-improving-sugar-recovery-from-c-massecuite-cooling-crystallisers/
    
    """
    _units = {'Crystallizer volume': 'm3',
              'Batch time': 'hr',
              'Loading time': 'hr'}
    
    _N_ins = _N_outs = 1
    
    #: [float] Cleaning and unloading time (hr).
    tau_0 = 1
    
    #: Optional[float] Loading time (hr).
    tau_loading = None
    
    #: [float] Fraction of filled tank to total tank volume.
    V_wf = 0.9
    
    def _get_design_info(self):
        return (('Cleaning and unloading time', self.tau_0, 'hr'),
                ('Working volume fraction', self.V_wf, ''))
    
    def _init(self, 
             tau=None, N=None, V=None, T=305.15,
             Nmin=2, Nmax=36, vessel_material='Carbon steel',
             kW=0.00746, tau_0=None, tau_loading=None,
        ):
        self._N = N; self._V = V
        
        #: [float] Reaction time [hr].
        self.tau = tau
        
        #: [int] Number of vessels
        if N: self.N = N
        
        #: [float] Target volume of a vessel
        if V: self.V = V
        
        #: [float] Operating temperature of reactor [K].
        self.T = T
        
        #: [int] Minimum number of vessels
        self.Nmin = Nmin
        
        #: [int] Maximum number of vessels
        self.Nmax = Nmax
        
        #: [float] Electricity usage per volume in kW/gal
        self.kW = kW
        
        # [str] Vessel material
        self.vessel_material = vessel_material
        
        if tau_0 is not None: self.tau_0 = tau_0
        if tau_loading is not None: self.tau_loading = tau_loading
        
    @property
    def vessel_material(self):
        return self._vessel_material
    @vessel_material.setter
    def vessel_material(self, vessel_material):
        if vessel_material in crystallizer_material_factors:
            self.F_M['Crystallizer'] = crystallizer_material_factors[vessel_material]
        else:
            raise AttributeError(
                f'invalid vessel material {vessel_material!r}; '
                f'valid materials include {set(crystallizer_material_factors)!r}'
            )
        self._vessel_material = vessel_material
    
    def _setup(self):
        super()._setup()
        effluent, = self.outs
        effluent.T = self.T
        effluent.P = 101325
    
    def _run(self):
        effluent, = self.outs
        effluent.copy_like(self.ins[0])
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
        dct = size_batch(v_0, tau, tau_0, V_wf, N_reactors=self.N, V_max=self.V, loading_time=self.tau_loading)
        Design['Crystallizer volume'] = volume = dct.pop('Reactor volume')
        Design['Number of crystallizers'] = N = dct.pop('Number of reactors')
        Design.update(dct)
        self.add_heat_utility(self.Hnet, self.T)
        self.add_power_utility(self.kW * V_wf * volume * N)
        
    