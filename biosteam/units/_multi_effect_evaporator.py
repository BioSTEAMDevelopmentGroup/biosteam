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
import biosteam as bst
from .. import Unit
from .mixing import Mixer
from .heat_exchange import HXutility
from ._flash import Flash, Evaporator_PQ
from .design_tools import (
    compute_vacuum_system_power_and_cost,
    compute_heat_transfer_area
)    
from thermosteam import MultiStream, settings
import flexsolve as flx
from warnings import warn
from .design_tools import heat_transfer as ht

__all__ = ('MultiEffectEvaporator',)

log = np.log
exp = np.exp

# Table 22.32 Product process and design (pg 592)
# Name: ('Area range (m2)', 'Cost(A) (USD)', 'U (kJ/(hr*m2*K)))', 'Material')
evaporators = {'Horizontal tube': 
                    ((9.29, 743.224),
                     lambda A, CE: CE*2.304*A**0.53,
                     4906.02,
                     'Carbon steel'),
               'Long-tube vertical':
                   ((9.29, 743.224),
                    lambda A, CE: CE*3.086*A**0.55,
                    8176.699,
                    'Carbon steel'),
               'Forced circulation': 
                   ((13.935, 8000),
                    lambda A, CE: CE/500*exp(8.2986 + 0.5329*log(A*0.0929)-0.000196*log(A*0.0929)**2),
                    10731.918,
                    'Carbon steel'),
               'Falling film': 
                   ((13.935, 371.612),
                    lambda A, CE: CE*7.416*A**0.55,
                    10220.874,
                    'Stainless steel tubes/Carbon steel shell')}


class MultiEffectEvaporator(Unit):
    """
    Creates evaporatorators with pressures given by P (a list of pressures). 
    Adjusts first evaporator vapor fraction to satisfy an overall fraction
    evaporated. All evaporators after the first have zero duty. Condenses
    the vapor coming out of the last evaporator. Pumps all liquid streams
    to prevent back flow in later parts. All liquid evaporated is ultimately
    recondensed. Cost is based on required heat transfer area. Vacuum system
    is based on air leakage. Air leakage is based on volume, as given by
    residence time `tau` and flow rate to each evaporator.

    Parameters
    ----------
    ins : stream
        Inlet.
    outs : stream sequence
        * [0] Solid-rich stream.
        * [1] Condensate stream.
    P : tuple[float]
        Pressures describing each evaporator (Pa).
    V : float
        Molar fraction evaporated as specified in `V_definition` 
        (either overall or in the first effect).
    V_definition : str, optional
        * 'Overall' - `V` is the overall molar fraction evaporated.
        * 'First-effect' - `V` is the molar fraction evaporated in the first effect.
    
    Examples
    --------
    Concentrate sugar setting vapor fraction at the first effect:
    
    >>> import biosteam as bst
    >>> from biorefineries.cornstover import chemicals
    >>> bst.settings.set_thermo(chemicals)
    >>> feed = bst.Stream('feed', Water=1000, Glucose=100, 
    ...                   AceticAcid=0.5, HMF=0.1, Furfural=0.1,
    ...                   units='kg/hr')
    >>> E1 = bst.MultiEffectEvaporator('E1', ins=feed, outs=('solids', 'liquid'), 
    ...                                V=0.1, V_definition='First-effect',
    ...                                P=(101325, 73581, 50892, 32777, 20000))
    >>> E1.simulate()
    >>> E1.show()
    MultiEffectEvaporator: E1
    ins...
    [0] feed
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Water       55.5
                        AceticAcid  0.00833
                        Furfural    0.00104
                        HMF         0.000793
                        Glucose     0.555
    outs...
    [0] solids
        phase: 'l', T: 333.24 K, P: 20000 Pa
        flow (kmol/hr): Water       20.6
                        AceticAcid  0.00189
                        Furfural    7.39e-05
                        HMF         0.000793
                        Glucose     0.555
    [1] liquid
        phase: 'l', T: 352.12 K, P: 101325 Pa
        flow (kmol/hr): Water       34.9
                        AceticAcid  0.00643
                        Furfural    0.000967
    
    >>> E1.results()
    Multi-Effect Evaporator                 Units        E1
    Power               Rate                   kW      5.72
                        Cost               USD/hr     0.447
    Low pressure steam  Duty                kJ/hr   5.8e+05
                        Flow              kmol/hr      14.9
                        Cost               USD/hr      3.55
    Cooling water       Duty                kJ/hr -3.49e+05
                        Flow              kmol/hr       239
                        Cost               USD/hr     0.116
    Design              Area                  m^2        11
                        Volume                m^3      1.64
    Purchase cost       Condenser             USD  5.35e+03
                        Evaporators           USD  9.59e+03
                        Liquid-ring pump      USD  1.24e+04
    Total purchase cost                       USD  2.74e+04
    Utility cost                           USD/hr      4.12
    
    Concentrate sugar setting overall vapor fraction:
    
    >>> import biosteam as bst
    >>> from biorefineries.cornstover import chemicals
    >>> bst.settings.set_thermo(chemicals)
    >>> feed = bst.Stream('feed', Water=1000, Glucose=100, 
    ...                   AceticAcid=0.5, HMF=0.1, Furfural=0.1,
    ...                   units='kg/hr')
    >>> E1 = bst.MultiEffectEvaporator('E1', ins=feed, outs=('solids', 'liquid'), 
    ...                                V=0.1, V_definition='Overall',
    ...                                P=(101325, 73581, 50892, 32777, 20000))
    >>> E1.simulate()
    >>> E1.show()
    MultiEffectEvaporator: E1
    ins...
    [0] feed
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Water       55.5
                        AceticAcid  0.00833
                        Furfural    0.00104
                        HMF         0.000793
                        Glucose     0.555
    outs...
    [0] solids
        phase: 'l', T: 354.94 K, P: 50892 Pa
        flow (kmol/hr): Water       50
                        AceticAcid  0.0069
                        Furfural    0.000579
                        HMF         0.000793
                        Glucose     0.555
    [1] liquid
        phase: 'l', T: 361.2 K, P: 101325 Pa
        flow (kmol/hr): Water       5.55
                        AceticAcid  0.00143
                        Furfural    0.000462
    
    >>> E1.results()
    Multi-Effect Evaporator                 Units        E1
    Power               Rate                   kW      5.72
                        Cost               USD/hr     0.447
    Low pressure steam  Duty                kJ/hr  3.82e+05
                        Flow              kmol/hr      9.85
                        Cost               USD/hr      2.34
    Cooling water       Duty                kJ/hr -1.15e+05
                        Flow              kmol/hr      78.5
                        Cost               USD/hr    0.0383
    Design              Area                  m^2      1.64
                        Volume                m^3     0.986
    Purchase cost       Condenser             USD  3.89e+03
                        Evaporators           USD  2.77e+03
                        Liquid-ring pump      USD  1.24e+04
    Total purchase cost                       USD  1.91e+04
    Utility cost                           USD/hr      2.83
    
    """
    line = 'Multi-Effect Evaporator'
    _units = {'Area': 'm^2',
              'Volume': 'm^3'}
    _BM = {'Evaporators': 2.45,
           'Liquid-ring pump': 1.0,
           'Condenser': 3.17}
    _N_outs = 2
    _N_heat_utilities = 2

    #: Residence time (hr)
    tau = 0.30

    # Evaporator type
    _Type = 'Forced circulation'
    
    # Data for simmulation and costing
    _evap_data = evaporators[_Type]

    @property
    def Type(self):
        """Evaporation type."""
        return self._Type
    @Type.setter
    def Type(self, evap_type):
        try:
            self._evap_data = evaporators[evap_type]
        except KeyError:
            dummy = str(evaporators.keys())[11:-2]
            raise ValueError(f"Type must be one of the following: {dummy}")
        self._Type = evap_type

    @property
    def V_definition(self):
        """[str] Must be one of the following:
        * 'Overall' - Defines attribute `V` as the overall molar fraction evaporated.
        * 'First-effect' - Defines attribute `V` as the molar fraction evaporated in the first effect.
        """
        return self._V_definition
    @V_definition.setter
    def V_definition(self, V_definition):
        V_definition = V_definition.capitalize()
        if V_definition in ('Overall', 'First-effect'):
            self._V_definition = V_definition
        else:
            raise ValueError("V_definition must be either 'Overall' or 'First-effect'")

    def __init__(self, ID='', ins=None, outs=(), thermo=None, *, P, V, V_definition='Overall'):
        Unit.__init__(self, ID, ins, outs, thermo)
        self.P = P #: tuple[float] Pressures describing each evaporator (Pa).
        self.V = V #: [float] Molar fraction evaporated.
        self.V_definition = V_definition
        self._V_first_effect = None
        
    def reset_cache(self):
        components = self.components
        evaporators = components['evaporators']
        for i in evaporators: i.reset_cache()
        components['condenser'].reset_cache()
        
    def load_components(self):
        P = self.P
        thermo = self.thermo
        
        # Create components
        self._N_evap = n = len(P) # Number of evaporators
        first_evaporator = Flash(None, outs=(None, None), P=P[0], thermo=thermo)
        
        # Put liquid first, then vapor side stream
        evaporators = [first_evaporator]
        for i in range(1, n):
            evap = Evaporator_PQ(None, outs=(None, None, None), P=P[i], Q=0, thermo=thermo)
            evaporators.append(evap)
        
        condenser = HXutility(None, outs=[None], thermo=thermo, V=0)
        self.heat_utilities = (first_evaporator.heat_utilities[0],
                               condenser.heat_utilities[0])
        mixer = Mixer(None, outs=[None], thermo=thermo)
        
        self.components = {'evaporators': evaporators,
                           'condenser': condenser,
                           'mixer': mixer}
        
        # Set-up components
        other_evaporators = evaporators[1:]
        first_evaporator.ins[:] = [i.copy() for i in self.ins]
        
        # Put liquid first, then vapor side stream
        ins = [first_evaporator.outs[1], first_evaporator.outs[0]]
        for evap in other_evaporators:
            evap.ins[:] = ins
            ins = [evap.outs[1], evap.outs[0]]
        
    def _V_overall(self, V_first_effect):
        first_evaporator, *other_evaporators = self.components['evaporators']
        first_evaporator.V = V_overall = V_first_effect
        first_evaporator._run()
        for evap in other_evaporators:
            evap._run()
            V_overall += (1. - V_overall) * evap.V
        return V_overall
        
    def _V_overall_objective_function(self, V_first_effect):
        return self._V_overall(V_first_effect) - self.V
    
    def _run(self):
        out_wt_solids, liq = self.outs
        ins = self.ins
        self.load_components()

        if self.V == 0:
            out_wt_solids.copy_like(ins[0])
            for i in self.heat_utilities: 
                i.empty(); i.heat_exchanger = None
            liq.empty()
            return
        
        if self.V_definition == 'Overall':
            P = tuple(self.P)
            self.P = list(P)
            for i in range(self._N_evap-1):
                if self._V_overall(0.) > self.V:
                    self.P.pop()
                    self.load_components()
                else:
                    break
            self.P = P
            self._V_first_effect = flx.IQ_interpolation(self._V_overall_objective_function,
                                                        0., 1., None, None, self._V_first_effect, 
                                                        xtol=1e-9, ytol=1e-6,
                                                        checkiter=False)
            V_overall = self.V
        else: 
            V_overall = self._V_overall(self.V)
    
        n = self._N_evap  # Number of evaporators
        components = self.components
        evaporators = components['evaporators']
        condenser = components['condenser']
        mixer = components['mixer']
        last_evaporator = evaporators[-1]
    
        # Condensing vapor from last effector
        outs_vap = last_evaporator.outs[0]
        condenser.ins[:] = [outs_vap]
        condenser._run()
        outs_liq = [condenser.outs[0]]  # list containing all output liquids

        # Unpack other output streams
        out_wt_solids.copy_like(last_evaporator.outs[1])
        for i in range(1, n):
            evap = evaporators[i]
            outs_liq.append(evap.outs[2])

        # Mix liquid streams
        mixer.ins[:] = outs_liq
        mixer._run()
        liq.copy_like(mixer.outs[0])
        
        mixed_stream = MultiStream(None, thermo=self.thermo)
        mixed_stream.copy_flow(self.ins[0])
        mixed_stream.vle(P=last_evaporator.P, V=V_overall)
        out_wt_solids.mol = mixed_stream.imol['l']
        liq.mol = mixed_stream.imol['g']
        
    def _design(self):
        if self.V == 0: return
        
        # This functions also finds the cost
        A_range, C_func, U, _ = self._evap_data
        components = self.components
        evaporators = components['evaporators']
        Design = self.design_results
        Cost = self.baseline_purchase_costs
        CE = bst.CE
        
        first_evaporator = evaporators[0]
        heat_exchanger = first_evaporator.heat_exchanger
        hu = heat_exchanger.heat_utilities[0]
        duty = heat_exchanger.Q = first_evaporator.H_out - first_evaporator.H_in
        Q = abs(duty)
        Tci = first_evaporator.ins[0].T
        Tco = first_evaporator.outs[0].T
        hu(duty, Tci, Tco)
        Th = hu.inlet_utility_stream.T
        LMTD = ht.compute_LMTD(Th, Th, Tci, Tco)
        ft = 1
        A = abs(compute_heat_transfer_area(LMTD, U, Q, ft))
        first_evaporator.baseline_purchase_costs['Evaporator'] = C = C_func(A, CE)
        self._evap_costs = evap_costs = [C]
        
        # Find condenser requirements
        condenser = components['condenser']
        condenser._design()
        condenser._cost()
        Cost['Condenser'] = condenser.purchase_cost
        
        # Find area and cost of evaporators
        As = [A]
        A_min, A_max = A_range
        evap = evaporators[-1]
        for evap in evaporators[1:]:
            Q = evap.design_results['Heat transfer']
            if Q <= 1e-12: 
                As.append(0.)
                evap_costs.append(0.)
            else:
                Tc = evap.outs[0].T
                Th = evap.outs[2].T
                LMTD = Th - Tc
                A = compute_heat_transfer_area(LMTD, U, Q, 1.)
                As.append(A)
                if settings.debug and not A_min < A < A_max:
                    warn(f'area requirement ({A}) is out of range, {A_range}')
                evap_costs.append(C_func(A, CE))
        self._As = As
        Design['Area'] = A = sum(As)
        Design['Volume'] = total_volume = self._N_evap * self.tau * self.ins[0].F_vol
        Cost['Evaporators'] = sum(evap_costs)
        
        # Calculate power
        power, cost = compute_vacuum_system_power_and_cost(
            F_mass=0, F_vol=0, P_suction=evap.outs[0].P,
            vessel_volume=total_volume,
            vacuum_system_preference='Liquid-ring pump')
        Cost['Liquid-ring pump'] = cost
        self.power_utility(power)
        
        
