# -*- coding: utf-8 -*-
"""
.. contents:: :local:

.. autoclass:: biosteam.units.aerated_bioreactor.AeratedBioreactor
.. autoclass:: biosteam.units.aerated_bioreactor.GasFedBioreactor

References
----------
.. [1] Benz, G. T. Optimize Power Consumption in Aerobic Fermenters. 
    Chem. Eng. Progress 2003, 99 (5), 100–103.

.. [2] Benz, G. T. Bioreactor Design for Chemical Engineers. Chem. Eng.\
    Progress 2011, 21–26.

.. [3] Seider, W. D., Lewin,  D. R., Seader, J. D., Widagdo, S., Gani, R.,
    & Ng, M. K. (2017). Product and Process Design Principles. Wiley.

"""
import biosteam as bst
from .stirred_tank_reactor import AbstractStirredTankReactor
from math import pi
import numpy as np
from scipy.constants import g
import flexsolve as flx
from warnings import filterwarnings, catch_warnings
from scipy.optimize import minimize_scalar, minimize, least_squares, differential_evolution
from biosteam.units.design_tools import aeration

__all__ = (
    'AeratedBioreactor', 'ABR',
    'GasFedBioreactor', 'GFB',
)

class AeratedBioreactor(AbstractStirredTankReactor):
    """
    Same as StirredTankReactor but includes aeration. The agitator power may
    vary to minimize the total power requirement of both the compressor and agitator
    yet achieve the required oxygen transfer rate.
    
    Examples
    --------
    >>> import biosteam as bst
    >>> from biorefineries.sugarcane import chemicals
    >>> bst.settings.set_thermo(chemicals)
    >>> feed = bst.Stream('feed',
    ...                   Water=1.20e+05,
    ...                   Glucose=2.5e+04,
    ...                   units='kg/hr',
    ...                   T=32+273.15)
    >>> # Model oxygen uptake as combustion
    >>> rxn = bst.Rxn('Glucose + O2 -> H2O + CO2', reactant='Glucose', X=0.5, correct_atomic_balance=True) 
    >>> R1 = bst.AeratedBioreactor(
    ...     'R1', ins=[feed, bst.Stream('air', phase='g')], outs=('vent', 'product'), tau=12, V_max=500,
    ...     reactions=rxn,
    ... )
    >>> R1.simulate()
    >>> R1.show()
    AeratedBioreactor: R1
    ins...
    [0] feed  
        phase: 'l', T: 305.15 K, P: 101325 Pa
        flow (kmol/hr): Water    6.66e+03
                        Glucose  139
    [1] air  
        phase: 'g', T: 305.15 K, P: 101325 Pa
        flow (kmol/hr): O2  1.02e+03
                        N2  3.85e+03
    outs...
    [0] vent  
        phase: 'g', T: 305.15 K, P: 101325 Pa
        flow (kmol/hr): Water  135
                        CO2    416
                        O2     606
                        N2     3.85e+03
    [1] product  
        phase: 'l', T: 305.15 K, P: 101325 Pa
        flow (kmol/hr): Water    6.94e+03
                        Glucose  69.4
    
    """
    _N_ins = 2
    _N_outs = 2
    _ins_size_is_fixed = False
    auxiliary_unit_names = (
        'compressor',
        'air_cooler',
        *AbstractStirredTankReactor.auxiliary_unit_names
    )
    T_default = 273.15 + 32 
    P_default = 101325
    kW_per_m3_default = 0.2955 # Reaction in homogeneous liquid; reference [1]
    batch_default = True
    default_methods = {
        'Stirred tank': 'Riet',
        'Bubble column': 'Dewes',
    }
    def _init(
            self, reactions, theta_O2=0.5, Q_O2_consumption=None,
            optimize_power=None, design=None, method=None, kLa_kwargs=None,
            cooler_pressure_drop=None, compressor_isentropic_efficiency=None,
            **kwargs,
        ):
        if compressor_isentropic_efficiency is None: compressor_isentropic_efficiency = 0.85
        #: Isentropic efficiency of the compressor. Defaults to 0.85.
        self.compressor_isentropic_efficiency = compressor_isentropic_efficiency 
        #: Pressure drop at the cooler [Pa]. Defaults to 20684.28 Pa, a heuristic value for a gas.
        self.cooler_pressure_drop = 20684.28 if cooler_pressure_drop is None else cooler_pressure_drop
        AbstractStirredTankReactor._init(self, **kwargs)
        self.reactions = reactions
        self.theta_O2 = theta_O2 # Average concentration of O2 in the liquid as a fraction of saturation.
        self.Q_O2_consumption = Q_O2_consumption # Forced duty per O2 consummed [kJ/kmol].
        self.optimize_power = True if optimize_power is None else optimize_power
        if design is None: 
            design = 'Stirred tank'
        elif design not in aeration.kLa_method_names:
            raise ValueError(
                f"{design!r} is not a valid design; only "
                f"{list(aeration.kLa_method_names)} are valid"
            )
        self.design = design
        if method is None:
            method = self.default_methods[design]
        if (key:=(design, method)) in aeration.kLa_methods:
            self.kLa = aeration.kLa_methods[key]
        elif hasattr(method, '__call__'):
            self.kLa = method
        else:
            raise ValueError(
                f"{method!r} is not a valid kLa method; only "
                f"{aeration.kLa_method_names[design]} are valid"
            )
        self.kLa_kwargs = {} if kLa_kwargs is None else kLa_kwargs
    
    def get_kLa(self):
        if self.kLa is aeration.kLa_stirred_Riet:
            V = self.get_design_result('Reactor volume', 'm3') * self.V_wf
            operating_time = self.tau / self.design_results.get('Batch time', 1.)
            N_reactors = self.parallel['self']
            P = 1000 * self.kW_per_m3 * V # W
            air_in = self.sparged_gas
            D = self.get_design_result('Diameter', 'm')
            F = air_in.get_total_flow('m3/s') / N_reactors / operating_time
            R = 0.5 * D
            A = pi * R * R
            self.superficial_gas_flow = U = F / A # m / s 
            return aeration.kLa_stirred_Riet(P, V, U, **self.kLa_kwargs) # 1 / s 
        elif self.kLa is aeration.kla_bubcol_Dewes:
            V = self.get_design_result('Reactor volume', 'm3') * self.V_wf
            operating_time = self.tau / self.design_results.get('Batch time', 1.)
            N_reactors = self.parallel['self']
            air_in = self.sparged_gas
            D = self.get_design_result('Diameter', 'm')
            F = air_in.get_total_flow('m3/s') / N_reactors / operating_time
            R = 0.5 * D
            A = pi * R * R
            self.superficial_gas_flow = U = F / A # m / s 
            feed = self.ins[0]
            try:
                return aeration.kla_bubcol_Dewes(U, feed.get_property('mu', 'mPa*s'), air_in.get_property('rho', 'kg/m3'))
            except:
                breakpoint()
        else:
            raise NotImplementedError('kLa method has not been implemented in BioSTEAM yet')
    
    def get_agitation_power(self, kLa):
        if self.kLa is aeration.kLa_stirred_Riet:
            air_in = self.sparged_gas
            N_reactors = self.parallel['self']
            operating_time = self.tau / self.design_results.get('Batch time', 1.)
            V = self.get_design_result('Reactor volume', 'm3') * self.V_wf
            D = self.get_design_result('Diameter', 'm')
            F = air_in.get_total_flow('m3/s') / N_reactors / operating_time
            R = 0.5 * D
            A = pi * R * R
            self.superficial_gas_flow = U = F / A # m / s
            return aeration.P_at_kLa_Riet(kLa, V, U, **self.kLa_kwargs)
        else:
            raise NotImplementedError('kLa method has not been implemented in BioSTEAM yet')
    
    def _get_duty(self):
        if self.Q_O2_consumption is None:
            H_in = sum([i.H for i in self.ins if i.phase != 'g'], self.air_cooler.outs[0].H)
            return self.H_out - H_in + self.Hf_out - self.Hf_in
        else:
            return self.Q_O2_consumption * (
                sum([i.imol['O2'] for i in self.ins])
                - sum([i.imol['O2'] for i in self.outs])
            )
    
    @property
    def feed(self):
        return self._ins[0]
    
    @property
    def air(self):
        for i in self._ins:
            if i.phase == 'g': return i
    
    @property
    def sparged_gas(self):
        return self.air_cooler.outs[0]
    
    @property
    def vent(self):
        return self._outs[0]
    
    def load_auxiliaries(self):
        super().load_auxiliaries()
        compressor = self.auxiliary(
            'compressor', bst.IsentropicCompressor, 
            self.air, eta=self.compressor_isentropic_efficiency, P=2 * 101325,
        )
        self.auxiliary(
            'air_cooler', bst.HXutility, compressor-0, T=self.T, 
            dP=False, 
        )
        
    def _run_vent(self, vent, effluent):
        vent.receive_vent(effluent, energy_balance=False, ideal=True)
        
    def _run(self):
        air = self.air
        if air is None:
            air = bst.Stream(phase='g', thermo=self.thermo)
            self.ins.insert(1, air)
            self.compressor.ins[0] = self.auxlet(air)
        feeds = [i for i in self.ins if i.phase != 'g']
        vent, effluent = self.outs
        air.P = vent.P = effluent.P = self.P
        air.T = vent.T = effluent.T = self.T
        vent.empty()
        vent.phase = 'g'
        air.phase = 'g'
        air.empty()
        compressor = self.compressor
        effluent.mix_from(feeds, energy_balance=False)
        self._run_reactions(effluent)
        effluent_no_air_data = effluent.get_data()
        OUR = -effluent.get_flow('mol/s', 'O2') # Oxygen uptake rate
        if OUR <= 1e-2:
            if OUR > 0: effluent.imol['O2'] = 0.
            self._run_vent(vent, effluent)
            return
        air_cc = self.sparged_gas
        air_cc.copy_like(air)
        air_cc.P = compressor.P = self._inlet_air_pressure()
        air_cc.T = self.T
        
        if self.optimize_power:
            def total_power_at_oxygen_flow(O2):
                air.set_flow([O2, O2 * 79. / 21.], 'mol/s', ['O2', 'N2'])
                air_cc.copy_flow(air) # Skip simulation of air cooler
                compressor.simulate()
                effluent.set_data(effluent_no_air_data)
                effluent.mix_from([effluent, air_cc], energy_balance=False)
                vent.empty()
                self._run_vent(vent, effluent)
                total_power = self._solve_total_power(OUR)
                return total_power
            
            f = total_power_at_oxygen_flow
            minimize_scalar(f, 1.2 * OUR, bounds=[OUR, 10 * OUR], tol=OUR * 1e-3)
        else:
            def air_flow_rate_objective(O2):
                air.set_flow([O2, O2 * 79. / 21.], 'mol/s', ['O2', 'N2'])
                air_cc.copy_flow(air) # Skip simulation of air cooler
                compressor.simulate()
                effluent.set_data(effluent_no_air_data)
                effluent.mix_from([effluent, air_cc], energy_balance=False)
                vent.empty()
                self._run_vent(vent, effluent)
                return OUR - self.get_OTR()
            
            f = air_flow_rate_objective
            y0 = air_flow_rate_objective(OUR)
            if y0 <= 0.: # Correlation is not perfect and special cases lead to OTR > OUR
                return
            flx.IQ_interpolation(f, x0=OUR, x1=10 * OUR, 
                                 y0=y0, ytol=1e-3, xtol=1e-3)
        
    def _run_reactions(self, effluent):
        self.reactions.force_reaction(effluent)
    
    def _solve_total_power(self, OUR): # For OTR = OUR [mol / s]
        air_in = self.sparged_gas
        N_reactors = self.parallel['self']
        operating_time = self.tau / self.design_results.get('Batch time', 1.)
        V = self.get_design_result('Reactor volume', 'm3') * self.V_wf
        vent = self.vent
        P_O2_air = air_in.get_property('P', 'bar') * air_in.imol['O2'] / air_in.F_mol
        P_O2_vent = 0. if vent.isempty() else vent.get_property('P', 'bar') * vent.imol['O2'] / vent.F_mol
        C_O2_sat_air = aeration.C_O2_L(self.T, P_O2_air) # mol / kg
        C_O2_sat_vent = aeration.C_O2_L(self.T, P_O2_vent) # mol / kg
        theta_O2 = self.theta_O2
        LMDF = aeration.log_mean_driving_force(C_O2_sat_vent, C_O2_sat_air, theta_O2 * C_O2_sat_vent, theta_O2 * C_O2_sat_air)
        kLa = OUR / (LMDF * V * self.effluent_density * N_reactors * operating_time) 
        P = self.get_agitation_power(kLa)
        agitation_power_kW = P / 1000
        total_power_kW = (agitation_power_kW + self.compressor.power_utility.consumption / N_reactors) / V
        self.kW_per_m3 = agitation_power_kW / V 
        return total_power_kW
    
    def get_OUR(self):
        """Return the oxygen uptake rate in mol/s."""
        feeds = [i for i in self.ins if i.phase != 'g']
        effluent = self.effluent.copy()
        effluent.mix_from(feeds, energy_balance=False)
        self._run_reactions(effluent)
        return -effluent.get_flow('mol/s', 'O2') # Oxygen uptake rate
    
    def get_OTR(self):
        """Return the oxygen transfer rate in mol/s."""
        V = self.get_design_result('Reactor volume', 'm3') * self.V_wf
        operating_time = self.tau / self.design_results.get('Batch time', 1.)
        N_reactors = self.parallel['self']
        kLa = self.get_kLa()
        air_in = self.sparged_gas
        vent = self.vent
        P_O2_air = air_in.get_property('P', 'bar') * air_in.imol['O2'] / air_in.F_mol
        P_O2_vent = vent.get_property('P', 'bar') * vent.imol['O2'] / vent.F_mol
        C_O2_sat_air = aeration.C_O2_L(self.T, P_O2_air) # mol / kg
        C_O2_sat_vent = aeration.C_O2_L(self.T, P_O2_vent) # mol / kg
        theta_O2 = self.theta_O2
        LMDF = aeration.log_mean_driving_force(C_O2_sat_vent, C_O2_sat_air, theta_O2 * C_O2_sat_vent, theta_O2 * C_O2_sat_air)
        OTR = kLa * LMDF * self.effluent_density * V * N_reactors * operating_time # mol / s
        return OTR
        
    def _inlet_air_pressure(self):
        AbstractStirredTankReactor._design(self)
        liquid = bst.Stream(None, thermo=self.thermo)
        liquid.mix_from([i for i in self.ins if i.phase != 'g'], energy_balance=False)
        liquid.copy_thermal_condition(self.outs[0])
        self.effluent_density = rho = liquid.rho
        length = self.get_design_result('Length', 'm') * self.V_wf
        return g * rho * length + 101325 + self.cooler_pressure_drop # Pa
    
    def _design(self):
        AbstractStirredTankReactor._design(self)
        if self.air.isempty(): return
        liquid = bst.Stream(None, thermo=self.thermo)
        liquid.mix_from([i for i in self.ins if i.phase != 'g'], energy_balance=False)
        liquid.copy_thermal_condition(self.outs[0])
        rho = liquid.rho
        length = self.get_design_result('Length', 'm') * self.V_wf
        compressor = self.compressor
        compressor.P = g * rho * length + 101325
        compressor.simulate()
        air_cooler = self.air_cooler
        air_cooler.T = self.T
        air_cooler.simulate()
        self.parallel['compressor'] = 1
        self.parallel['air_cooler'] = 1
        # For robust process control, do not include in HXN
        for unit in self.auxiliary_units:
            for hu in unit.heat_utilities: hu.hxn_ok = False


class GasFedBioreactor(AbstractStirredTankReactor):
    """
    Same as AbstractStirredTankReactor but includes multiple gas feeds. The agitator power may
    vary to minimize the total power requirement of both the compressor and agitator
    yet achieve the required oxygen transfer rate.
    
    Examples
    --------
    >>> import biosteam as bst
    >>> bst.settings.set_thermo(['H2', 'CO2', 'N2', 'O2', 'H2O', 'AceticAcid'])
    >>> media = bst.Stream(ID='media', H2O=10000, units='kg/hr')
    >>> H2 = bst.Stream(ID='H2', H2=100, units='kg/hr', phase='g')
    >>> fluegas = bst.Stream(ID='fluegas', N2=70, CO2=25, H2O=3, O2=2, units='kg/hr', phase='g')
    >>> # Model acetic acid production from H2 and CO2
    >>> rxn = bst.Rxn('H2 + CO2 -> AceticAcid + H2O', reactant='H2', correct_atomic_balance=True) 
    >>> brxn = rxn.backwards(reactant='AceticAcid')
    >>> R1 = bst.GasFedBioreactor(
    ...     'R1', ins=[media, H2, fluegas], outs=('vent', 'product'), tau=68, V_max=500,
    ...     reactions=rxn, backward_reactions=brxn,
    ...     gas_substrates=('H2', 'CO2'),
    ...     titer={'AceticAcid': 5},
    ...     optimize_power=False,
    ...     kW_per_m3=0.,
    ... )
    >>> R1.simulate()
    >>> R1.show()
    GasFedBioreactor: R1
    ins...
    [0] media  
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): H2O  188
    [1] H2  
        phase: 'g', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): H2  49.6
    [2] fluegas  
        phase: 'g', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): CO2  0.568
                        N2   2.5
                        O2   0.0625
                        H2O  0.167
    outs...
    [0] vent  
        phase: 'g', T: 305.15 K, P: 101325 Pa
        flow (kmol/hr): H2          48.5
                        N2          2.5
                        O2          0.0625
                        H2O         1.89
                        AceticAcid  0.00183
    [1] product  
        phase: 'l', T: 305.15 K, P: 101325 Pa
        flow (kmol/hr): H2O         187
                        AceticAcid  0.282
    
    """
    _N_ins = 2
    _N_outs = 2
    _ins_size_is_fixed = False
    auxiliary_unit_names = (
        'sparger',
        'compressors',
        'gas_coolers',
        *AbstractStirredTankReactor.auxiliary_unit_names
    )
    T_default = 273.15 + 32 
    P_default = 101325
    kW_per_m3_default = 0.2955 # Reaction in homogeneous liquid
    batch_default = True
    default_methods = AeratedBioreactor.default_methods
    get_kLa = AeratedBioreactor.get_kLa
    get_agitation_power = AeratedBioreactor.get_agitation_power
    
    def _init(self, 
            reactions, gas_substrates, 
            # Can vary either liquid or gas flows
            titer=None, 
            # Only for variable gas flows
            backward_reactions=None, 
            variable_gas_feeds=(), 
            # General design/performance arguments
            design=None, method=None, kLa_kwargs=None,
            theta=0.5, Q_consumption=None,
            cooler_pressure_drop=None,
            # Only for agitated bioreactors (not bubble column)
            optimize_power=None, 
            **kwargs,
        ):
        self.cooler_pressure_drop = 20684.28 if cooler_pressure_drop is None else cooler_pressure_drop
        self.reactions = reactions
        self.backward_reactions = backward_reactions
        self.theta = theta # Average concentration of gas substrate in the liquid as a fraction of saturation.
        self.Q_consumption = Q_consumption # Forced duty per gas substrate consummed [kJ/kmol].
        self.kLa_kwargs = {} if kLa_kwargs is None else kLa_kwargs
        self.optimize_power = True if optimize_power is None else optimize_power
        self.variable_gas_feeds = variable_gas_feeds # list[int|Stream] Feed index or stream.
        self.gas_substrates = gas_substrates
        self.titer = titer # dict[str, float] g / L
        AbstractStirredTankReactor._init(self, **kwargs)
        if design is None: 
            if self.kW_per_m3 == 0:
                design = 'Bubble column'
            else:
                design = 'Stirred tank'
        elif design not in aeration.kLa_method_names:
            raise ValueError(
                f"{design!r} is not a valid design; only "
                f"{list(aeration.kLa_method_names)} are valid"
            )
        self.design = design
        if method is None:
            method = self.default_methods[design]
        if (key:=(design, method)) in aeration.kLa_methods:
            self.kLa = aeration.kLa_methods[key]
        elif hasattr(method, '__call__'):
            self.kLa = method
        else:
            raise ValueError(
                f"{method!r} is not a valid kLa method; only "
                f"{aeration.kLa_method_names[design]} are valid"
            )
    
    def _get_duty(self):
        if self.Q_consumption is None:
            H_in = sum(
                [i.H for i in self.ins if i.phase != 'g']
                + [i.outs[0].H for i in self.gas_coolers]
            )
            return self.H_out - H_in + self.Hf_out - self.Hf_in
        else:
            return self.Q_consumption * (
                sum([i.imol['O2'] for i in self.ins])
                - sum([i.imol['O2'] for i in self.outs])
            )
    
    @property
    def vent(self):
        return self._outs[0]
    
    @property
    def variable_gas_feeds(self):
        return [(i if isinstance(i, bst.Stream) else self.ins[i]) for i in self._variable_gas_feeds]
    @variable_gas_feeds.setter
    def variable_gas_feeds(self, variable_gas_feeds):
        self._variable_gas_feeds = variable_gas_feeds
    
    @property
    def normal_gas_feeds(self):
        variable = set(self.variable_gas_feeds)
        return [i for i in self.ins if i not in variable and i.phase == 'g']
    
    @property
    def liquid_feeds(self):
        return [i for i in self.ins if i.phase != 'g']
    
    @property
    def sparged_gas(self):
        return self.sparger-0
    
    def load_auxiliaries(self):
        super().load_auxiliaries()
        self.compressors = []
        self.gas_coolers = []
        for i in self.variable_gas_feeds: i.phase = 'g'
        for inlet in self.ins:
            if inlet.phase != 'g': continue
            compressor = self.auxiliary(
                'compressors', bst.IsentropicCompressor, inlet, eta=0.85, P=2 * 101325
            )
            self.auxiliary(
                'gas_coolers', bst.HXutility, compressor-0, T=self.T
            )
        self.auxiliary(
            'sparger', bst.Mixer, [i-0 for i in self.gas_coolers]
        )
        
    def _run_vent(self, vent, effluent):
        vent.receive_vent(effluent, energy_balance=False, ideal=True)
        
    def get_SURs(self, effluent):
        F_vol = effluent.ivol['Water'] # m3 / hr
        produced = bst.Stream(None, thermo=self.thermo)
        for ID, concentration in self.titer.items():
            produced.imass[ID] = F_vol * concentration
        consumed = produced.copy()
        self.backward_reactions.force_reaction(consumed)
        return consumed.get_flow('mol/s', self.gas_substrates), consumed, produced
    
    def _load_gas_feeds(self):
        for i in self.compressors: i.simulate()
        for i in self.gas_coolers: i.simulate()
        self.sparger.simulate()
    
    def _run(self):
        variable_gas_feeds = self.variable_gas_feeds
        vent, effluent = self.outs
        vent.P = effluent.P = self.P
        vent.T = effluent.T = self.T
        vent.empty()
        vent.phase = 'g'
        if not self.titer:
            liquid_feeds = [i for i in self.ins if i.phase != 'g']
            T = self.T
            P = self._inlet_gas_pressure()
            for i in self.compressors: i.P = P
            for i in self.gas_coolers: i.T = T
            self._load_gas_feeds()
            STRs = self.get_STRs()
            effluent.mix_from(liquid_feeds, energy_balance=False)
            effluent.set_flow(STRs, units='mol/s', key=self.gas_substrates)
            self._run_reactions(effluent)
            vent.empty()
            self._run_vent(vent, effluent)
        elif variable_gas_feeds:
            liquid_feeds = [i for i in self.ins if i.phase != 'g']
            effluent.mix_from(liquid_feeds, energy_balance=False)
            T = self.T
            P = self._inlet_gas_pressure()
            for i in self.compressors: i.P = P
            for i in self.gas_coolers: i.T = T
            effluent_liquid_data = effluent.get_data()
            SURs, s_consumed, s_produced = self.get_SURs(effluent) # Gas substrate uptake rate [mol / s]
            if (SURs <= 1e-2).all():
                effluent.imol[self.gas_substrates] = 0.
                self._run_vent(vent, effluent)
                return
            x_substrates = []
            for gas, ID in zip(self.variable_gas_feeds, self.gas_substrates):
                x_substrates.append(gas.get_molar_fraction(ID))
            index = range(len(self.gas_substrates))
            
            def load_flow_rates(F_feeds):
                for i in index:
                    gas = variable_gas_feeds[i]
                    gas.set_total_flow(F_feeds[i], 'mol/s')
                self._load_gas_feeds()
                effluent.set_data(effluent_liquid_data)
                effluent.mix_from([self.sparged_gas, -s_consumed, s_produced, *liquid_feeds], energy_balance=False)
                if (effluent.mol < 0).any(): breakpoint()
                vent.empty()
                self._run_vent(vent, effluent)
            
            baseline_feed = bst.Stream.sum(self.normal_gas_feeds, energy_balance=False)
            baseline_flows = baseline_feed.get_flow('mol/s', self.gas_substrates)
            bounds = np.array([[max(1.01 * SURs[i] - baseline_flows[i], 1e-6), 10 * SURs[i]] for i in index])
            if self.optimize_power:
                def total_power_at_substrate_flow(F_substrates):
                    load_flow_rates(F_substrates)
                    total_power = self._solve_total_power(SURs)
                    return total_power
                
                f = total_power_at_substrate_flow
                with catch_warnings():
                    filterwarnings('ignore')
                    results = minimize(f, 1.2 * SURs, bounds=bounds, tol=SURs.max() * 1e-6)
                    load_flow_rates(results.x / x_substrates)
            else:
                def gas_flow_rate_objective(F_substrates):
                    F_feeds = F_substrates / x_substrates
                    load_flow_rates(F_feeds)
                    STRs = self.get_STRs() # Must meet all substrate demands
                    F_ins = F_substrates + baseline_flows
                    mask = STRs - F_ins > 0
                    STRs[mask] = F_ins[mask]
                    diff = SURs - STRs
                    diff[diff > 0] *= 1e3 # Force transfer rate to meet uptake rate
                    SE = (diff * diff).sum()
                    return SE
                
                f = gas_flow_rate_objective
                with catch_warnings():
                    filterwarnings('ignore')
                    bounds = bounds.T
                    results = least_squares(f, 1.2 * SURs, bounds=bounds, ftol=SURs.min() * 1e-6)
                self._results = results
                load_flow_rates(results.x / x_substrates)
        else:
            try:
                feed, = [i for i in self.ins if i.phase != 'g']
            except:
                raise RuntimeError('gas-fed bioreactor must have exactly on liquid feed')
            T = self.T
            P = self._inlet_gas_pressure()
            for i in self.compressors: i.P = P
            for i in self.gas_coolers: i.T = T
            self._load_gas_feeds()
            substrates = sum([i.get_flow(units='mol/s', key=self.gas_substrates) for i in self.ins])
            F_liquid_max = self._initialize_variable_liquid_guess(feed, effluent, vent)
            product, titer = next(iter(self.titer.items()))
            def liquid_flow_rate_objective(F_feed):
                feed.F_mass = F_feed
                
                def f(vent_flow_rates):
                    self.vent.mol = vent_flow_rates
                    STRs = self.get_STRs()
                    STRs = np.minimum(STRs, substrates)
                    effluent.mix_from(self.ins, energy_balance=False)
                    effluent.set_flow(STRs, units='mol/s', key=self.gas_substrates)
                    remaining = substrates - STRs
                    self._run_reactions(effluent)
                    vent.empty()
                    self.vent.set_flow(remaining, units='mol/s', key=self.gas_substrates)
                    self._run_vent(vent, effluent)
                    return self.vent.mol.to_array()
                    
                flx.aitken(f, self.vent.mol.to_array(), checkiter=False, xtol=1e-3, checkconvergence=False)
                return effluent.imass[product] / effluent.ivol['Water'] - titer
                
            flx.IQ_interpolation(liquid_flow_rate_objective, 0.05 * F_liquid_max, F_liquid_max, ytol=1e-3)
        # self.show()
        # breakpoint()
        
    def _run_reactions(self, effluent):
        data = effluent.get_data()
        rxns = self.reactions
        rxns.force_reaction(effluent)
        for reactant in self.gas_substrates:
            if effluent.imol[reactant] < 0:
                if isinstance(rxns, bst.Rxn):
                    rxns.reactant = reactant
                else:
                    for rxn in rxns:
                        if rxn.istoichiometry[reactant] < 0:
                            rxn.reactant = reactant
            
                effluent.set_data(data)
                rxns.force_reaction(effluent)
                break
        
    def _initialize_variable_liquid_guess(self, feed, effluent, vent):
        effluent.mix_from(self.ins)
        data = vent.get_data()
        vent.empty()
        self._run_reactions(effluent)
        product, titer = next(iter(self.titer.items()))
        F_liquid_max = 1000 * effluent.imass[product] / titer # kg / hr
        vent.set_data(data)
        return F_liquid_max
        
    def _solve_total_power(self, SURs): # For STR = SUR [mol / s]
        gas_in = self.sparged_gas
        N_reactors = self.parallel['self']
        operating_time = self.tau / self.design_results.get('Batch time', 1.)
        V = self.get_design_result('Reactor volume', 'm3') * self.V_wf
        D = self.get_design_result('Diameter', 'm')
        F = gas_in.get_total_flow('m3/s') / N_reactors / operating_time
        R = 0.5 * D
        A = pi * R * R
        self.superficial_gas_flow = U = F / A # m / s 
        vent = self.vent
        Ps = []
        for gas_substrate, SUR in zip(self.gas_substrates, SURs):
            Py_gas = gas_in.get_property('P', 'bar') * gas_in.imol[gas_substrate] / gas_in.F_mol
            Py_vent = 0. if vent.isempty() else vent.get_property('P', 'bar') * vent.imol[gas_substrate] / vent.F_mol
            C_sat_gas = aeration.C_L(self.T, Py_gas, gas_substrate) # mol / kg
            C_sat_vent = aeration.C_L(self.T, Py_vent, gas_substrate) # mol / kg
            theta = self.theta
            LMDF = aeration.log_mean_driving_force(C_sat_vent, C_sat_gas, theta * C_sat_vent, theta * C_sat_gas)
            kLa = SUR / (LMDF * V * self.effluent_density * N_reactors * operating_time)
            Ps.append(aeration.P_at_kLa_Riet(kLa, V, U, **self.kLa_kwargs))
        P = max(Ps)  
        agitation_power_kW = P / 1000
        compressor_power_kW = sum([i.power_utility.consumption for i in self.compressors]) / N_reactors
        total_power_kW = (agitation_power_kW + compressor_power_kW) / V
        self.kW_per_m3 = agitation_power_kW / V 
        return total_power_kW
    
    def get_STRs(self):
        """Return the gas substrate transfer rate in mol/s."""
        V = self.get_design_result('Reactor volume', 'm3') * self.V_wf
        operating_time = self.tau / self.design_results.get('Batch time', 1.)
        N_reactors = self.parallel['self']
        gas_in = self.sparged_gas
        kLa = self.get_kLa() # 1 / s 
        vent = self.vent
        P_gas = gas_in.get_property('P', 'bar')
        P_vent = vent.get_property('P', 'bar')
        STRs = []
        for ID in self.gas_substrates:
            Py_gas = P_gas * gas_in.imol[ID] / gas_in.F_mol
            Py_vent = P_vent * vent.imol[ID] / (vent.F_mol or 1)
            C_sat_gas = aeration.C_L(self.T, Py_gas, ID) # mol / kg
            C_sat_vent = aeration.C_L(self.T, Py_vent, ID) # mol / kg
            theta = self.theta
            LMDF = aeration.log_mean_driving_force(C_sat_vent, C_sat_gas, theta * C_sat_vent, theta * C_sat_gas)
            STRs.append(
                kLa * LMDF * self.effluent_density * V * N_reactors * operating_time # mol / s
            )
        return np.array(STRs)
        
    def _inlet_gas_pressure(self):
        AbstractStirredTankReactor._design(self, size_only=True)
        liquid = bst.Stream(None, thermo=self.thermo)
        liquid.mix_from([i for i in self.ins if i.phase != 'g'], energy_balance=False)
        liquid.copy_thermal_condition(self.outs[0])
        self.effluent_density = rho = liquid.rho
        length = self.get_design_result('Length', 'm') * self.V_wf
        return g * rho * length + 101325 + self.cooler_pressure_drop # Pa
    
    def _design(self):
        AbstractStirredTankReactor._design(self)
        self.parallel['sparger'] = 1
        self.parallel['mixers'] = 1
        self.parallel['compressors'] = 1
        self.parallel['gas_coolers'] = 1
        # For robust process control, do not include in HXN
        for unit in self.auxiliary_units:
            for hu in unit.heat_utilities: hu.hxn_ok = False
    
GFB = GasFedBioreactor
ABR = AeratedBioreactor
