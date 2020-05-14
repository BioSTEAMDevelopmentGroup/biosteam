# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 22:45:47 2018

@author: yoelr
"""
import numpy as np
from ._hx import HXutility
from .. import Unit
from .decorators import cost
from .design_tools import size_batch

__all__ = ('Bioreactor',)

@cost('Reactor volume', 'Cleaning in place', CE=521.9,
      cost=421e3, S=3785, n=0.6, BM=1.8, N='N')
@cost('Reactor volume', 'Agitators', CE=521.9, cost=52500,
      S=3785, n=0.5, kW=22.371, BM=1.5, N='N')
@cost('Reactor volume', 'Reactors', CE=521.9, cost=844000,
      S=3785, n=0.5, BM=1.5, N='N')
class Bioreactor(Unit, isabstract=True):
    """
    Abstract Bioreactor class. Conversion is based on reaction time, `tau`.
    Cleaning and unloading time,`tau_0`, fraction of working volume, `V_wf`,
    and number of reactors, `N_reactors`, are attributes that can be changed.
    Cost of a reactor is based on the NREL batch fermentation tank cost
    assuming volumetric scaling with a 6/10th exponent [1]_. 
    
    **Abstract methods**
        
    kinetic_model(z, t, *kinetic_constants) : 
        A staticmethod that returns effluent concentrations (z_t; kg/m^3)
        given the starting concentration (z; kg/m^3), reaction time (t; hr),
        and kinetic constants.
    
    _run():
        Must set outlet stream flows.
    
    Parameters
    ----------
    ins : streams
        Inlet fluids to be mixed into the fermentor.
    outs : stream sequence
        * [0] Vent
        * [1] Effluent
    tau : float
        Reaction time.
    efficiency=0.9 : float, optional
        User enforced efficiency.
    iskinetic=False: bool, optional
        If True, `Fermenation.kinetic_model` will be used.
    N : int
        Number of batch reactors
    T=305.15 : float
        Temperature of reactor [K].
    
    References
    ----------
    .. [1] D. Humbird, R. Davis, L. Tao, C. Kinchin, D. Hsu, and A. Aden National. Renewable Energy Laboratory Golden, Colorado. P. Schoen, J. Lukas, B. Olthof, M. Worley, D. Sexton, and D. Dudgeon. Harris Group Inc. Seattle, Washington and Atlanta, Georgia. Process Design and Economics for Biochemical Conversion of Lignocellulosic Biomass to Ethanol Dilute-Acid Pretreatment and Enzymatic Hydrolysis of Corn Stover. May 2011. Technical Report NREL/TP-5100-47764
    
    """
    _units = {'Reactor volume': 'm3',
              'Cycle time': 'hr',
              'Batch time': 'hr',
              'Loading time': 'hr',
              'Total dead time': 'hr'}
    _N_ins = _N_outs = 2
    _N_heat_utilities = 1
    _has_power_utility = True
    
    #: [bool] If True, number of reactors (N) is chosen as to minimize installation cost in every simulation. Otherwise, N remains constant.
    autoselect_N = False
    
    #: [float] Cleaning and unloading time (hr).
    tau_0 = 3
    
    #: [float] Fraction of filled tank to total tank volume.
    V_wf = 0.9
    
    def _get_design_info(self):
        return (('Cleaning and unloading time', self.tau_0, 'hr'),
                ('Working volume fraction', self.V_wf, ''),
                ('Number of reactors', self.N, ''))
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *, 
                 tau,  N, iskinetic=False, T=305.15):
        Unit.__init__(self, ID, ins, outs, thermo)
        self.tau = tau
        self.N = N
        self.T = T
        self.cooler = HXutility(None)
        
    def _setup(self):
        vent, effluent = self.outs
        vent.phase = 'g'
        self.cooler._outs[0].T = effluent.T = vent.T = self.T
       
    @property
    def N(self):
        """[int] Number of reactors"""
        return self._N
    @N.setter
    def N(self, N):
        if N <= 1:
            raise ValueError(f"number of reactors must be greater than 1, value {N} is infeasible")
        self._N = N

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
        effluent = self.outs[1]
        v_0 = effluent.F_vol
        tau = self._tau
        tau_0 = self.tau_0
        Design = self.design_results
        if self.autoselect_N: self._N = self.N_at_minimum_capital_cost
        N = self._N
        Design.update(size_batch(v_0, tau, tau_0, N, self.V_wf))
        hx_effluent = effluent.copy()
        hx_effluent.phase = 'l'
        hx_effluent.mol[:] /= N
        cooler = self.cooler
        cooler.simulate_as_auxiliary_exchanger(self.Hnet/N, hx_effluent)
        hu_bioreactor, = self.heat_utilities
        hu_cooler, = cooler.heat_utilities
        hu_bioreactor.copy_like(hu_cooler)
        hu_bioreactor.scale(N)
        self.purchase_costs['Coolers'] = self.cooler.purchase_costs['Heat exchanger'] * N
        
    