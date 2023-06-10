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
import biosteam as bst
from .. import Unit
from .mixing import Mixer
from .heat_exchange import HXutility
from ._flash import Flash, Evaporator_PQ, Evaporator_PV
from .design_tools import (
    compute_vacuum_system_power_and_cost,
    compute_heat_transfer_area
)    
from thermosteam import MultiStream, settings
import flexsolve as flx
from warnings import warn
import numpy as np
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
    ins : 
        Inlet.
    outs : 
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
    >>> bst.process_tools.default()
    >>> from biorefineries.cellulosic import create_cellulosic_ethanol_chemicals
    >>> bst.settings.set_thermo(create_cellulosic_ethanol_chemicals())
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
        phase: 'l', T: 333.21 K, P: 20000 Pa
        flow (kmol/hr): Water       20.5
                        AceticAcid  0.00183
                        Furfural    5.96e-05
                        HMF         0.000793
                        Glucose     0.555
    [1] liquid
        phase: 'l', T: 352.11 K, P: 20000 Pa
        flow (kmol/hr): Water       35
                        AceticAcid  0.0065
                        Furfural    0.000981
    
    >>> E1.results()
    Multi-effect evaporator                                    Units       E1
    Electricity         Power                                     kW     5.72
                        Cost                                  USD/hr    0.447
    Low pressure steam  Duty                                   kJ/hr 5.83e+05
                        Flow                                 kmol/hr     15.1
                        Cost                                  USD/hr     3.58
    Cooling water       Duty                                   kJ/hr -3.5e+05
                        Flow                                 kmol/hr      239
                        Cost                                  USD/hr    0.117
    Design              Area                                     m^2       11
                        Volume                                   m^3     1.24
    Purchase cost       Evaporators                              USD 9.59e+03
                        Condenser - Double pipe                  USD 5.36e+03
                        Vacuum system - Liquid-ring pump...      USD 1.24e+04
    Total purchase cost                                          USD 2.74e+04
    Utility cost                                              USD/hr     4.15
    
    Concentrate sugar setting overall vapor fraction:
    
    >>> import biosteam as bst
    >>> from biorefineries.cellulosic import create_cellulosic_ethanol_chemicals
    >>> bst.settings.set_thermo(create_cellulosic_ethanol_chemicals())
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
        phase: 'l', T: 354.91 K, P: 50892 Pa
        flow (kmol/hr): Water       50
                        AceticAcid  0.00683
                        Furfural    0.000535
                        HMF         0.000793
                        Glucose     0.555
    [1] liquid
        phase: 'l', T: 361.13 K, P: 50892 Pa
        flow (kmol/hr): Water       5.55
                        AceticAcid  0.0015
                        Furfural    0.000506
    
    >>> E1.results()
    Multi-effect evaporator                                    Units        E1
    Electricity         Power                                     kW      5.72
                        Cost                                  USD/hr     0.447
    Low pressure steam  Duty                                   kJ/hr  3.84e+05
                        Flow                                 kmol/hr      9.94
                        Cost                                  USD/hr      2.36
    Cooling water       Duty                                   kJ/hr -1.15e+05
                        Flow                                 kmol/hr      78.8
                        Cost                                  USD/hr    0.0384
    Design              Area                                     m^2      1.64
                        Volume                                   m^3      6.52
    Purchase cost       Evaporators                              USD  2.77e+03
                        Condenser - Double pipe                  USD   3.9e+03
                        Vacuum system - Liquid-ring pump...      USD  1.24e+04
    Total purchase cost                                          USD  1.91e+04
    Utility cost                                              USD/hr      2.85
    
    >>> E1.results()
    Multi-effect evaporator                                    Units        E1
    Electricity         Power                                     kW      5.72
                        Cost                                  USD/hr     0.447
    Low pressure steam  Duty                                   kJ/hr  3.84e+05
                        Flow                                 kmol/hr      9.94
                        Cost                                  USD/hr      2.36
    Cooling water       Duty                                   kJ/hr -1.15e+05
                        Flow                                 kmol/hr      78.8
                        Cost                                  USD/hr    0.0384
    Design              Area                                     m^2      1.64
                        Volume                                   m^3      6.52
    Purchase cost       Evaporators                              USD  2.77e+03
                        Condenser - Double pipe                  USD   3.9e+03
                        Vacuum system - Liquid-ring pump...      USD  1.24e+04
    Total purchase cost                                          USD  1.91e+04
    Utility cost                                              USD/hr      2.85
    
    >>> E1.results()
    Multi-effect evaporator                                    Units        E1
    Electricity         Power                                     kW      5.72
                        Cost                                  USD/hr     0.447
    Low pressure steam  Duty                                   kJ/hr  3.84e+05
                        Flow                                 kmol/hr      9.94
                        Cost                                  USD/hr      2.36
    Cooling water       Duty                                   kJ/hr -1.15e+05
                        Flow                                 kmol/hr      78.8
                        Cost                                  USD/hr    0.0384
    Design              Area                                     m^2      1.64
                        Volume                                   m^3      6.52
    Purchase cost       Evaporators                              USD  2.77e+03
                        Condenser - Double pipe                  USD   3.9e+03
                        Vacuum system - Liquid-ring pump...      USD  1.24e+04
    Total purchase cost                                          USD  1.91e+04
    Utility cost                                              USD/hr      2.85
    
    Concentrate sugar setting overall vapor fraction:
    
    >>> import biosteam as bst
    >>> from biorefineries.cellulosic import create_cellulosic_ethanol_chemicals
    >>> bst.settings.set_thermo(create_cellulosic_ethanol_chemicals())
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
        phase: 'l', T: 354.91 K, P: 50892 Pa
        flow (kmol/hr): Water       50
                        AceticAcid  0.00683
                        Furfural    0.000535
                        HMF         0.000793
                        Glucose     0.555
    [1] liquid
        phase: 'l', T: 361.13 K, P: 50892 Pa
        flow (kmol/hr): Water       5.55
                        AceticAcid  0.0015
                        Furfural    0.000506
    
    >>> E1.results()
    Multi-effect evaporator                                    Units        E1
    Electricity         Power                                     kW      5.72
                        Cost                                  USD/hr     0.447
    Low pressure steam  Duty                                   kJ/hr  3.84e+05
                        Flow                                 kmol/hr      9.94
                        Cost                                  USD/hr      2.36
    Cooling water       Duty                                   kJ/hr -1.15e+05
                        Flow                                 kmol/hr      78.8
                        Cost                                  USD/hr    0.0384
    Design              Area                                     m^2      1.64
                        Volume                                   m^3      6.52
    Purchase cost       Evaporators                              USD  2.77e+03
                        Condenser - Double pipe                  USD   3.9e+03
                        Vacuum system - Liquid-ring pump...      USD  1.24e+04
    Total purchase cost                                          USD  1.91e+04
    Utility cost                                              USD/hr      2.85
    
    """
    line = 'Multi-effect evaporator'
    auxiliary_unit_names = ('condenser', 'mixer', 'vacuum_system')
    _units = {'Area': 'm^2',
              'Volume': 'm^3'}
    _F_BM_default = {'Evaporators': 2.45,
                     'Vacuum system': 1.0,
                     'Condenser': 3.17}
    _N_outs = 2

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

    def __init__(self, ID='', ins=None, outs=(), thermo=None, *, P, V, V_definition='Overall',
                 flash=True, chemical='7732-18-5'):
        Unit.__init__(self, ID, ins, outs, thermo)
        self.P = P #: tuple[float] Pressures describing each evaporator (Pa).
        self.V = V #: [float] Molar fraction evaporated.
        self.V_definition = V_definition
        self.flash = flash #: [bool] Whether to perform a flash calculation to account for volatile components.
        self._V_first_effect = None
        self._reload_components = True
        self.chemical = chemical
        
    def reset_cache(self, isdynamic=None):
        self._reload_components = True
        
    def _load_components(self):
        P = self.P
        thermo = self.thermo
        
        # Create components
        self._N_evap = n = len(P) # Number of evaporators
        
        if self.flash:
            first_evaporator = Flash(None, outs=(None, None), P=P[0], thermo=thermo)
        else:
            first_evaporator = Evaporator_PV(None, outs=(None, None), P=P[0], thermo=thermo, chemical=self.chemical)
        # Put liquid first, then vapor side stream
        self.evaporators = evaporators = [first_evaporator]
        for i in range(1, n):
            evap = Evaporator_PQ(None, outs=(None, None, None), P=P[i], Q=0, thermo=thermo, chemical=self.chemical)
            evaporators.append(evap)
        
        self.condenser = HXutility(None, outs=[None], thermo=thermo, V=0)
        self.mixer = Mixer(None, outs=[None], thermo=thermo)
        
        # Set-up components
        other_evaporators = evaporators[1:]
        first_evaporator.ins[:] = [i.copy() for i in self.ins]
        
        # Put liquid first, then vapor side stream
        ins = [first_evaporator.outs[1], first_evaporator.outs[0]]
        for evap in other_evaporators:
            evap.ins[:] = ins
            ins = [evap.outs[1], evap.outs[0]]
        
    def _V_overall(self, V_first_effect):
        first_evaporator, *other_evaporators = self.evaporators
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
        if self.V == 0:
            out_wt_solids.copy_like(ins[0])
            liq.empty()
            self._reload_components = True
            return
        
        if self._reload_components:
            self._load_components()
            self._reload_components = False
        else:
            self.evaporators[0].ins[:] = [i.copy() for i in ins]
        
        if self.V_definition == 'Overall':
            P = tuple(self.P)
            self.P = list(P)
            for i in range(self._N_evap - 1):
                if self._V_overall(0.) > self.V:
                    self.P.pop()
                    self._load_components()
                    self._reload_components = True
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
        evaporators = self.evaporators
        condenser = self.condenser
        mixer = self.mixer
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
        mixer.simulate()
        liq.copy_like(mixer.outs[0])
        
        if self.flash:
            mixed_stream = MultiStream(None, thermo=self.thermo)
            mixed_stream.copy_flow(self.ins[0])
            mixed_stream.vle(P=last_evaporator.P, V=V_overall)
            out_wt_solids.mol = mixed_stream.imol['l']
            liq.mol = mixed_stream.imol['g']
        liq.P = out_wt_solids.P
        
    def _design(self):
        if self.V == 0: 
            for i in self.auxiliary_units: i._setup()
            return
        
        # This functions also finds the cost
        A_range, C_func, U, _ = self._evap_data
        evaporators = self.evaporators
        Design = self.design_results
        Cost = self.baseline_purchase_costs
        CE = bst.CE
        
        first_evaporator = evaporators[0]
        first_evaporator.simulate(run=False)
        hx = first_evaporator.heat_exchanger
        self.heat_utilities.append(hx.heat_utilities[0])
        
        # Cost first evaporators
        duty = hx.total_heat_transfer
        Q = abs(duty)
        Tci = first_evaporator.heat_exchanger.ins[0].T
        Tco = first_evaporator.outs[0].T
        hu = first_evaporator.heat_utilities[0]
        Th = hu.inlet_utility_stream.T
        LMTD = ht.compute_LMTD(Th, Th, Tci, Tco)
        ft = 1
        A = abs(compute_heat_transfer_area(LMTD, U, Q, ft))
        self._evap_costs = evap_costs = [C_func(A, CE)]
        
        # Find condenser requirements
        condenser = self.condenser
        condenser.simulate(run=False)
        
        # Find area and cost of evaporators
        As = [A]
        evap = evaporators[-1]
        bounds_warning = bst.exceptions.bounds_warning
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
                bounds_warning(self, 'heat transfer area requirement', A, 'ft2', A_range, 'cost')
                evap_costs.append(C_func(A, CE))
        self._As = As
        Design['Area'] = A = sum(As)
        vapor_sep_design = first_evaporator.design_results
        L = vapor_sep_design['Length']
        D = vapor_sep_design['Diameter']
        R = D / 2.
        volume = 0.0283168466 * np.pi * L * R * R # m3
        Design['Volume'] = total_volume = self._N_evap * volume
        Cost['Evaporators'] = sum(evap_costs)
        
        self.vacuum_system = bst.VacuumSystem(
            self, 'Steam-jet ejector', vessel_volume=total_volume, P_suction=self.outs[0].P,
        )
            
        
        
