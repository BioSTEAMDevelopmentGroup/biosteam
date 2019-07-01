# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 21:43:13 2018

@author: yoelr
"""
import numpy as np
import biosteam as bst
from .. import Unit, Stream
from scipy.optimize import brentq
from . import Mixer, HXutility
from ._flash import Evaporator_PV, Evaporator_PQin
from .designtools import vacuum_system
import ht
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
    """Creates evaporatorators with pressures given by P (a list of pressures). Adjusts first evaporator vapor fraction to satisfy an overall fraction evaporated. All evaporators after the first have zero duty. Condenses the vapor coming out of the last evaporator. Pumps all liquid streams to prevent back flow in later parts. All liquid evaporated is ultimately recondensed. Cost is based on required heat transfer area. Vacuum system is based on air leakage. Air leakage is based on volume, as given by residence time `tau` and flow rate to each evaporator.

    **Parameters**

        **component:** *[str]* Component being evaporated
         
        **P:** *[tuple]* Pressures describing each evaporator (Pa)
         
        **V:** *[float]* Overall molar fraction of component evaporated
         
        **P_liq:** *[tuple]* Liquid pressure after pumping (Pa)
    
    """
    _units = {'Area': 'm^2',
              'Volume': 'm^3'}
    _has_power_utility = True
    _N_heat_utilities = 2
    BM = 2.45
    line = 'Multi-Effect Evaporator'
    _kwargs = {'component': 'Water',
               'P': (101325,),  
               'V': 0.5, 
               'P_liq': 101325}

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

    def _init(self):
        # Unpack
        component, P, V, P_liq = (self._kwargs[i]
                                  for i in ('component', 'P', 'V', 'P_liq'))
        out_wt_solids, liq = self.outs
        
        # Create components
        self._N_evap = n = len(P) # Number of evaporators
        Stream.species = liq.species
        evap0 = Evaporator_PV(None, outs=(None, None),
                             component=component, P=P[0])
        
        evaporators = [evap0]
        for i in range(1, n):
            evap = Evaporator_PQin(None,
                                   outs=(None, None, None),
                                   component=component, P=P[i], Qin=0)
            evaporators.append(evap)
        condenser = HXutility(None, outs=Stream(None), V=0)
        evap0._heat_utilities[0], condenser._heat_utilities[0] = self._heat_utilities
        mixer = Mixer(None, outs=Stream(None))
        
        def V_error(v1):
            # Run first evaporator
            v_test = v1
            evap0._kwargs['V'] = v1
            evap0._run()
            # Put liquid first, then vapor side stream
            ins = [evap0.outs[1], evap0.outs[0]]
            for i in range(1, n):
                evap = evaporators[i]
                evap._ins[:] = ins
                evap._run()
                v_test += (1-v_test) * evap._V
                # Put liquid first, then vapor side stream
                ins = [evap.outs[1], evap.outs[0]]
            return V - v_test
        self._V_error = V_error
        self.components = {'evaporators': evaporators,
                           'condenser': condenser,
                           'mixer': mixer}
        

    def _run(self):
        component, P, V, P_liq = (self._kwargs[i]
                                  for i in ('component', 'P', 'V', 'P_liq'))
        out_wt_solids, liq = self.outs
        ins = self.ins

        n = self._N_evap  # Number of evaporators

        # Set-up components
        components = self.components
        evaporators = components['evaporators']
        evaporators[0].ins[:] = [Stream.like(i, None) for i in ins]
        condenser = components['condenser']
        mixer = components['mixer']
        brentq(self._V_error, 0.0001, 0.9909, xtol=0.0001)
        
        # Condensing vapor from last effector
        
        outs_vap = evaporators[-1].outs[0]
        condenser.ins[:] = [outs_vap]
        condenser._run()
        outs_liq = [condenser.outs[0]]  # list containing all output liquids

        # Unpack other output streams
        out_wt_solids.copylike(evaporators[-1].outs[1])
        for i in range(1, n):
            evap = evaporators[i]
            outs_liq.append(evap.outs[2])

        # Mix liquid streams
        mixer.ins[:] = outs_liq
        mixer._run()
        liq.copylike(mixer.outs[0])
        
    def _design(self):
        # This functions also finds the cost
        A_range, C_func, U, _ = self._evap_data
        components = self.components
        evaporators = components['evaporators']
        Design = self._Design
        Cost = self._Cost
        CE = bst.CE
        
        evap0 = evaporators[0]
        hu = evap0._heat_utilities[0]
        duty = evap0._H_out - evap0._H_in
        hu(duty, evap0.ins[0].T, evap0.outs[0].T)
        Q = abs(duty)
        Tci = evap0.ins[0].T
        Tco = evap0.outs[0].T
        Th = evap0._heat_utilities[0]._fresh.T
        LMTD = ht.LMTD(Th, Th, Tci, Tco)
        ft = 1
        A = HXutility._calc_area(LMTD, U, Q, ft)
        self._evap_costs = evap_costs = [C_func(A, CE)]
        
        # Find condenser requirements
        condenser = components['condenser']
        condenser._design()
        condenser._cost()
        Cost['Condenser'] = condenser._Cost['Heat exchanger']
        
        # Find area and cost of evaporators
        As = [A]
        for evap in evaporators[1:]:
            Q = evap._Qin
            Tc = evap.outs[0].T
            Th = evap.outs[2].T
            LMTD = Th - Tc
            A = HXutility._calc_area(LMTD, U, Q, ft)
            As.append(A)
            if not A_range[0] < A < A_range[1]:
                print(f'WARNING, area requirement ({A}) is out of range, {A_range}')
            evap_costs.append(C_func(A, CE))
        self._As = As
        Design['Area'] = A = sum(As)
        Design['Volume'] = vol = self._N_evap * self.tau * self.ins[0].volnet
        Cost['Evaporators'] = sum(evap_costs)
        
        # Calculate power
        power, cost = vacuum_system(massflow=0, volflow=0,
                                    P_suction=evap.outs[0].P, vol=vol,
                                    vacuum_system_preference='Liquid-ring pump')
        Cost['Vacuum liquid-ring pump'] = cost
        self._power_utility(power)
        
        
