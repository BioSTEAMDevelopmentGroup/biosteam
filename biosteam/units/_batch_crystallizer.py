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
from .._unit import design_input
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
    Batch crystallizer                                      Units            C1
    Electricity              Power                             kW          3.54
                             Cost                          USD/hr         0.277
    Cooling water            Duty                           kJ/hr     -1.69e+06
                             Flow                         kmol/hr      8.06e+03
                             Cost                          USD/hr          3.93
    Design inputs            Crystallization time              hr            20
                             Cleaning and unloading time       hr             1
                             Working volume fraction                        0.9
                             Operating temperature              K           313
                             Power per volume              kW/gal       0.00746
                             Vessel material                       Carbon steel
                             Crystallizer volume               m3           175
    Design results           Crystallizer volume               m3           132
                             Number of crystallizers                          4
                             Batch time                        hr            28
                             Loading time                      hr             7
    Purchase cost            Crystallizer (x4)                USD      2.09e+06
    Total purchase cost                                       USD      2.09e+06
    Installed equipment cost                                  USD      4.17e+06
    Utility cost                                           USD/hr          4.21
    
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
    
    tau = design_input('Crystallization time', 'hr', lb=0, default=None)
    tau_0 = design_input('Cleaning and unloading time', 'hr', lb=0, default=1)
    tau_loading = design_input(
        'Loading time', 'hr', lb=0,
        doc="Loading time of batch crystallizer [hr]. If not given, it will "
            "assume each vessel is constantly being filled.",
        default=None,
    )
    V_wf = design_input('Working volume fraction', lb=0, ub=1, default=0.9)
    T = design_input('Operating temperature', 'K', lb=0, default=305.15)
    kW = design_input(
        'Power per volume', 'kW/gal', lb=0, default=0.00746,
        doc="Electricity usage per volume in kW/gal. Defaults to 0.00746, "
            "a heuristic value for suspension of solids.",
    )
    vessel_material = design_input(
        catalog=crystallizer_material_factors,
        purchase_item='Crystallizer',
        catalog_type='Material',
        default='Carbon steel',
    )
    
    @design_input('Number of crystallizers', lb=2, default=None)
    def N(self) -> int:
        """Number of crystallizers. Either N or V may be given, not both."""
        return self._N
    @N.setter
    def N(self, N):
        if N is not None and getattr(self, '_V', None):
            raise AttributeError('cannot specify both crystallizer volume and number of crystallizers')
        self._N = N

    @design_input('Crystallizer volume', 'm3', lb=1, default=None)
    def V(self):
        """Target volume of crystallizers [m3]. Either N or V may be given, not both."""
        return self._V
    @V.setter
    def V(self, V):
        if V is not None and getattr(self, '_N', None):
            raise AttributeError('cannot specify both crystallizer volume and number of crystallizers')
        self._V = V
    
    def _init(self, Nmin=2, Nmax=36):
        #: [int] Minimum number of vessels
        self.Nmin = Nmin
        
        #: [int] Maximum number of vessels
        self.Nmax = Nmax
    
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
        
    