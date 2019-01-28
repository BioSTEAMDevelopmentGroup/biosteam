# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 21:43:13 2018

@author: yoelr
"""
from biosteam import Unit, Stream, np
from scipy.optimize import brentq
from biosteam.units import Flash_PV, Flash_PQin, Pump, Mixer, HXutility
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
    """Creates evaporatorators with pressures given by P (a list of pressures). Adjusts first evaporator vapor fraction to satisfy an overall fraction evaporated. All evaporators after the first have zero duty. Condenses the vapor coming out of the last evaporator. Pumps all liquid streams to prevent back flow in later parts. All liquid evaporated is ultimately recondensed.

    **Parameters**

         **component:** *[str]* Component being evaporated
         
         **P:** *[tuple]* Pressures describing each evaporator (Pa)
         
         **V:** *[float]* Overall molar fraction of component evaporated
         
         **P_liq:** *[tuple]* Liquid pressure after pumping (Pa)
    
    """
    kwargs = {'component': 'Water',
              'P': (101325,),  
              'V': 0.5, 
              'P_liq': 101325}

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
            raise ValueError(f"Evaporator type must be one of the following: {dummy}")
        self._Type = evap_type

    def setup(self):
        # Unpack
        component, P, V, P_liq = (self.kwargs[i]
                                  for i in ('component', 'P', 'V', 'P_liq'))
        [out_wt_solids, liq] = self.outs
        
        # Create components
        n = len(P) # Number of evaporators
        evap = Flash_PV('Evaporator0', outs=('vap', 'liq'),
                        component=component, P=P[0])
        evaporators = [evap]
        for i in range(1, n):
            evap = Flash_PQin('Evaporator' + str(i),
                              outs=('vap', 'liq', 'utility'),
                              component=component, P=P[i], Qin=0)
            evaporators.append(evap)
        condenser = HXutility('Condenser', outs='condensate', V=0)
        pumps = Pump('Pump', outs=('liq',)*(n+1), P=P_liq)
        mixer = Mixer('Mixer', outs='liq')

        evap0 = evaporators[0]
        def V_error(v1):
            # Run first evaporator
            v_test = v1
            evap0.kwargs['V'] = v1
            evap0.run()
            # Put liquid first, then vapor side stream
            ins = [evap0.outs[1], evap0.outs[0]]
            for i in range(1, n):
                evap = evaporators[i]
                evap.ins = ins
                evap.run()
                v_test += (1-v_test) * evap._V
                # Put liquid first, then vapor side stream
                ins = [evap.outs[1], evap.outs[0]]
            return V - v_test
        self._V_error = V_error
        self.components = {'evaporators': evaporators,
                           'condenser': condenser, 'pumps': pumps, 'mixer': mixer}

    def run(self):
        component, P, V, P_liq = (self.kwargs[i]
                                  for i in ('component', 'P', 'V', 'P_liq'))
        [out_wt_solids, liq] = self.outs
        ins = self.ins

        n = len(P)  # Number of evaporators

        # Set-up components
        components = self.components
        evaporators = components['evaporators']
        evaporators[0].ins = [Stream.like(i, '*') for i in ins]
        condenser = components['condenser']
        pumps = components['pumps']
        mixer = components['mixer']
        brentq(self._V_error, 0.001, 0.999, xtol=0.0001)
        
        # Condensing vapor from last effector
        
        outs_vap = evaporators[-1].outs[0]
        condenser.ins = [outs_vap]
        condenser.run()
        outs_liq = [condenser.outs[0]]  # list containing all output liquids

        # Unpack other output streams
        outwtsolids = (evaporators[-1].outs[1])
        for i in range(1, n):
            evap = evaporators[i]
            outs_liq.append(evap.outs[2])

        # Pump output streams (These could be eductors or multiple pumps)
        pumps.ins = [outwtsolids] + outs_liq
        pumps.run()
        out_wt_solids.copy_like(pumps.outs[0])
        outs_liq = pumps.outs[1:]

        # Mix liquid streams
        mixer.ins = outs_liq
        mixer.run()
        liq.copy_like(mixer.outs[0])

    def _lazy_run(self):
        feed = self.ins[0]
        out_wt_solids, liq = self.outs
        component, V = (self.kwargs[i] for i in ('component', 'V'))

        component_pos = feed.Settings.ID_index_dictionary[component]
        liq.mass[component_pos] = V*feed.mol[component_pos]
        out_wt_solids.mol = feed.mol - liq.mol

    def operation(self):
        components = self.components
        
        # Find utilities
        evap0 = components['evaporators'][0]
        evap0.operation()
        condenser = components['condenser']
        condenser.operation()
        
        # Attach utilities
        self.heat_utilities = evap0.heat_utilities + condenser.heat_utilities
        
    def design(self):
        """
        * 'Area': Total area of all evaporators (m^2)
        """
        # This functions also finds the cost
        A_range, C_func, U, _ = self._evap_data
        components = self.components
        evaporators = components['evaporators']
        r = self.results
        Design = r['Design']
        Cost = r['Cost']
        CE = self.CEPCI
        
        evap0 = evaporators[0]
        Q = abs(evap0.heat_utilities[0].results['Duty'])
        Tci = evap0.ins[0].T
        Tco = evap0.outs[0].T
        Th = evap0.heat_utilities[0]._fresh.T
        LMTD = ht.LMTD(Th, Th, Tci, Tco)
        ft = 1
        A = HXutility._calc_area(LMTD, U, Q, ft)
        self._evap_costs = evap_costs = [C_func(A, CE)]
        
        # Find condenser requirements
        condenser = components['condenser']
        condenser.design()
        Cost['Condenser'] = condenser.cost()['Heat exchanger']
        
        # Find area ands cost of evaporators
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
            
        Design['Area'] = sum(As)
        Cost['Evaporators'] = sum(evap_costs)
        S = condenser.outs[-1].volnet*0.5886 # To ft3/min
        Cost['Vacuum liquid ring pump'] = self.CEPCI/567*8250*S**0.37
        # pumps = components['pumps']
        # pumps_cost = 0
        # for i in range(len(As)+1):
        #     pumps.operation(); pumps.design()
        #     pumps_cost += sum(pumps.cost().values())
        #     del pumps.ins[0]
        #     del pumps.outs[0]
        # Cost['Pumps'] = pumps_cost
        return Design
            
    def cost(self):
        """
        * 'Evaporators': Sum of all evaporator costs (USD)
        * 'Condenser': (USD)
        * 'Vacuum liquid ring pump': (USD)
        """
        return self.results['Cost']
        
        
        
        
