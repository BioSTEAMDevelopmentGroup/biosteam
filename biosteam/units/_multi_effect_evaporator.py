# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 21:43:13 2018

@author: yoelr
"""
import numpy as np
import biosteam as bst
from .. import Unit
from ._mixer import Mixer
from ._hx import HXutility
from ._flash import Evaporator_PV, Evaporator_PQ
from .design_tools import (
    compute_vacuum_system_power_and_cost,
    compute_heat_transfer_area
)    
from thermosteam import MultiStream, Stream, settings
from flexsolve import IQ_interpolation
from warnings import warn
import ht

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
    component : str
                Component being evaporated.
    P : tuple[float]
        Pressures describing each evaporator (Pa).
    V : float
        Overall molar fraction of component evaporated.
    P_liq : tuple
            Liquid pressure after pumping (Pa).
    
    """
    _units = {'Area': 'm^2',
              'Volume': 'm^3'}
    _N_outs = 2
    _N_heat_utilities = 2
    BM = 2.45
    line = 'Multi-Effect Evaporator'

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

    def __init__(self, ID='', ins=None, outs=(), thermo=None, *, P, V):
        Unit.__init__(self, ID, ins, outs, thermo)
        # Unpack
        out_wt_solids, liq = self.outs
        self.V = V #: [float] Overall molar fraction of component evaporated.
        self._V1 = V/2.
        
        # Create components
        self._N_evap = n = len(P) # Number of evaporators
        first_evaporator = Evaporator_PV(None, outs=(None, None), P=P[0])
        
        # Put liquid first, then vapor side stream
        evaporators = [first_evaporator]
        for i in range(1, n):
            evap = Evaporator_PQ(None, outs=(None, None, None), P=P[i], Q=0)
            evaporators.append(evap)
        
        condenser = HXutility(None, outs=Stream(None), V=0)
        self.heat_utilities = (first_evaporator.heat_utilities[0],
                               condenser.heat_utilities[0])
        mixer = Mixer(None, outs=Stream(None))
        
        self.components = {'evaporators': evaporators,
                           'condenser': condenser,
                           'mixer': mixer}
        
    def _run(self):
        out_wt_solids, liq = self.outs
        ins = self.ins

        n = self._N_evap  # Number of evaporators

        # Set-up components
        components = self.components
        evaporators = components['evaporators']
        first_evaporator, *other_evaporators = evaporators
        first_evaporator.ins[:] = [i.copy() for i in ins]
        condenser = components['condenser']
        mixer = components['mixer']
        
        # Put liquid first, then vapor side stream
        ins = [first_evaporator.outs[1], first_evaporator.outs[0]]
        for evap in other_evaporators:
            evap.ins[:] = ins
            ins = [evap.outs[1], evap.outs[0]]
        
        def compute_overall_vapor_fraction(v1):
            v_overall = v1
            first_evaporator.V = v1
            first_evaporator._run()
            for evap in other_evaporators:
                evap._run()
                v_overall += (1-v_overall) * evap.V
            return v_overall
        
        x0 = 0.0001
        x1 = 0.9990
        y0 = compute_overall_vapor_fraction(x0)
        y1 = compute_overall_vapor_fraction(x1)
        self._V1 = IQ_interpolation(compute_overall_vapor_fraction,
                                    x0, x1, y0, y1, self._V1, self.V, 
                                    xtol=0.0001, ytol=0.001)
        # Condensing vapor from last effector
        outs_vap = evaporators[-1].outs[0]
        condenser.ins[:] = [outs_vap]
        condenser._run()
        outs_liq = [condenser.outs[0]]  # list containing all output liquids

        # Unpack other output streams
        out_wt_solids.copy_like(evaporators[-1].outs[1])
        for i in range(1, n):
            evap = evaporators[i]
            outs_liq.append(evap.outs[2])

        # Mix liquid streams
        mixer.ins[:] = outs_liq
        mixer._run()
        liq.copy_like(mixer.outs[0])
        
        mixed_stream = MultiStream(thermo=self.thermo)
        mixed_stream.copy_flow(self.ins[0])
        mixed_stream.vle(P=evaporators[-1].P, V=self.V)
        out_wt_solids.mol = mixed_stream.imol['l']
        liq.mol = mixed_stream.imol['g']
        
    def _design(self):
        # This functions also finds the cost
        A_range, C_func, U, _ = self._evap_data
        components = self.components
        evaporators = components['evaporators']
        Design = self.design_results
        Cost = self.purchase_costs
        CE = bst.CE
        
        first_evaporator = evaporators[0]
        hu = first_evaporator.heat_utilities[0]
        duty = first_evaporator.H_out - first_evaporator.H_in
        Q = abs(duty)
        Tci = first_evaporator.ins[0].T
        Tco = first_evaporator.outs[0].T
        hu(duty, Tci, Tco)
        Th = hu.inlet_utility_stream.T
        LMTD = ht.LMTD(Th, Th, Tci, Tco)
        ft = 1
        A = abs(compute_heat_transfer_area(LMTD, U, Q, ft))
        self._evap_costs = evap_costs = [C_func(A, CE)]
        
        # Find condenser requirements
        condenser = components['condenser']
        condenser._design()
        condenser._cost()
        Cost['Condenser'] = condenser.purchase_costs['Heat exchanger']
        
        # Find area and cost of evaporators
        As = [A]
        A_min, A_max = A_range
        for evap in evaporators[1:]:
            Q = evap._Q
            Tc = evap.outs[0].T
            Th = evap.outs[2].T
            LMTD = Th - Tc
            A = compute_heat_transfer_area(LMTD, U, Q, ft)
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
        Cost['Vacuum liquid-ring pump'] = cost
        self.power_utility(power)
        
        
