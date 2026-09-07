# -*- coding: utf-8 -*-
"""
.. contents:: :local:

.. autoclass:: biosteam.units.aerated_bioreactor.AeratedBioreactor

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
from .abstract_stirred_tank_reactor import AbstractStirredTankReactor
from math import pi
from scipy.constants import g
import flexsolve as flx
from scipy.optimize import minimize_scalar
from biosteam.units.design_tools import aeration

__all__ = (
    'AeratedBioreactor', 'ABR',
)

class AeratedBioreactor(AbstractStirredTankReactor):
    """
    Create an aerated bioreactor which satisfies the oxygen mass tranfer 
    requirement of the mass balance. The reactor is designed as a pressure vessel with a given aspect ratio and 
    residence time. A pump-heat exchanger recirculation loop can be used to satisfy 
    the duty, if any. By default, a turbine agitator is also included if the 
    power usage, `kW_per_m3`, is positive. A vacuum system is also 
    automatically added if the operating pressure is at a vacuum. 

    Parameters
    ----------
    theta_O2 : 
        Fraction of oxygen saturation in the broth. Defaults to 0.5.
    Q_O2_consumption :
        Forced duty per O2 consummed [kJ/kmol].
    optimize_power :
        If true, the agitator power is solved to minimize the total power 
        requirement of both the compressor and agitator such that the 
        required oxygen transfer rate is met.
    design : 
        Bioreactor design configuration. Valid options include 'Stirred tank'
        and 'Bubble column'. Defaults to the former.
    method :
        Method to calculate the overall mass transfer coefficient, kLa. 
        Can be a name or a function. Valid method names are listed in 
        `biosteam.aeration.kLa_methods`.
        For stirred tanks, defaults to the 'Riet'. 
        For bubble columns, defaults to 'Dewes'.
    kLa_kwargs :
        Additional arguments to pass to the kLa method.
    cooler_pressure_drop :
        Pressure drop at the cooler [Pa]. Defaults to 20684.28 Pa, 
        a heuristic value for a gas.
    compressor_isentropic_efficiency :
        Isentropic efficiency of the compressor. Defaults to 0.85.
    tau :
        Residence time [hr].
    T : 
        Operating temperature [K].
    P : 
        Operating pressure [Pa].
    V_wf : 
        Fraction of working volume over total volume. Defaults to 0.8.
    length_to_diameter :
        Length to diameter ratio of bioreactor.
    V_max :
        Maximum volume of a reactor [m3]. Defaults to 355.
    kW_per_m3 : 
        Power usage of agitator. Defaults to 0.2955 [kW / m3] from [1]_.
    vessel_material : 
        Vessel material. Defaults to 'Stainless steel 316'.
    vessel_type : 
        Vessel type. Valid options are 'Horizontal' or 'Vertical'. Defaults to 'Vertical'
    batch :
        Whether to use batch operation mode. If False, operation mode is continuous.
        Defaults to `continuous`.
    tau_0 : 
        Cleaning and unloading time (if batch mode). Defaults to 3 hr.
    N :
        Number of reactors.
    heat_exchanger_configuration : 
        What kind of heat exchanger to default to (if any). Valid options include 
        'jacketed', 'recirculation loop', and 'internal coil'. Defaults to 'recirculation loop'.
    dT_hx_loop : 
        Maximum change in temperature for the heat exchanger loop. Defaults to 5 K.
    jacket_annular_diameter :
        Annular diameter of heat exchanger jacket to vessel [m]. Defaults to 0.1 m.
    loading_time :
        Loading time of batch reactor. If not given, it will assume each vessel is constantly
        being filled.
        
    Notes
    -----
    The heat exchanger configuration can be one of the following:

    * 'recirculation loop': 
        The recirculation loop takes into account the required flow rate needed to
        reach the maximum temperature change of the heat exchanger, `dT_hx_loop`. 
        Increasing `dT_hx_loop` decreases the required recirculation flow rate and
        therefore decreases pump costs.
        
        When parallel reactors are required, one recirculation loop (each with a
        pump and heat exchanger) is assumed. Although it is possible to use the
        same recirculation loop for all reactors, this conservative assumption allows
        for each reactor to be operated independently from each other.

    * 'jacketed':
        The jacket does not account for the heat transfer area requirement. 
        It simply assumes that a full jacket can provide the necessary heat transfer
        area to meet the duty requirement. A heuristic annular diameter is assumed
        through `jacket_annular_diameter` (which can be adjusted by the user).
        The temperature at the wall is assumed to be the operating temperature.
        The weight of the jacket is added to the weight of the vessel and the
        cost is compounded together as a jacketed vessel.
        
    * 'internal coil':
        The internal coil is costed as an ordinary helical tube heat exchanger
        with the added assumption that the temperature at the wall is the 
        operating temperature. This method is still not implemented in BioSTEAM
        yet.

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
    >>> rxn = bst.Rxn(
    ...     'Glucose + O2 -> H2O + CO2', reactant='Glucose', X=0.5,
    ...     correct_atomic_balance=True
    ... ) 
    >>> R1 = bst.AeratedBioreactor(
    ...     ins=[feed, bst.Stream('air', phase='g')],
    ...     outs=('vent', 'product'), tau=12, V_max=500,
    ...     reactions=rxn,
    ... )
    >>> R1.simulate()
    >>> R1.show()
    AeratedBioreactor: R2
    ins...
    [0] feed  
        phase: 'l', T: 305.15 K, P: 101325 Pa
        flow (kmol/hr): Water    6.66e+03
                        Glucose  139
    [1] air  
        phase: 'g', T: 305.15 K, P: 101325 Pa
        flow (kmol/hr): O2  1.02e+03
                        N2  3.84e+03
    outs...
    [0] vent  
        phase: 'g', T: 305.15 K, P: 101325 Pa
        flow (kmol/hr): Water  240
                        CO2    416
                        O2     605
                        N2     3.84e+03
    [1] product  
        phase: 'l', T: 305.15 K, P: 101325 Pa
        flow (kmol/hr): Water    6.84e+03
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
            self, theta_O2=0.5, Q_O2_consumption=None,
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
        self.theta_O2 = theta_O2 # Average concentration of O2 in the liquid as a fraction of saturation.
        self.Q_O2_consumption = Q_O2_consumption # Forced duty per O2 consummed [kJ/kmol].
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
        self.optimize_power = (
            design == 'Stirred tank' 
            if optimize_power is None else
            optimize_power
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
            raise NotImplementedError(
                'cannot solve for the required agitation power using the '
               f'kLa method {self.kLa!r} in BioSTEAM yet'
            )
    
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
        aeration.vent_broth(vent, effluent)
        
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
        compressor.P = self._inlet_air_pressure()
        air_cc.P = compressor.P - self.cooler_pressure_drop
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
            minimize_scalar(f, 1.2 * OUR, bounds=[OUR, 10 * OUR], options=dict(xtol=OUR * 1e-3))
        else:
            def air_flow_rate_objective(O2):
                air.set_flow([O2, O2 * 79. / 21.], 'mol/s', ['O2', 'N2'])
                air_cc.copy_flow(air) # Skip simulation of air cooler
                effluent.set_data(effluent_no_air_data)
                effluent.mix_from([effluent, air_cc], energy_balance=False)
                vent.empty()
                self._run_vent(vent, effluent)
                return OUR - self.get_OTR()
            
            f = air_flow_rate_objective
            x0 = OUR
            y0 = air_flow_rate_objective(OUR)
            
            # Correlation is not perfect and special cases lead to OTR > OUR
            if y0 <= 0.: return 
            
            f = air_flow_rate_objective
            x1 = 10 * OUR
            y1 = air_flow_rate_objective(x1)
            
            if y1 > 0:
                x1 = 50 * OUR
                y1 = air_flow_rate_objective(x1)
            
            # It is possible an infinite flow rate of air is not enough to 
            # satisfy mass transfer if the titer is super high or gas O2 concentration
            # is too low.
            if y1 > 0: 
                raise RuntimeError(
                    'bioreactor conversion/titer cannot be satisfied; '
                    'even an infinite flow rate of gas is not enough'
                )
            
            # There is a known solution because y1 < 0 < y0
            flx.IQ_interpolation(
                f, x0=x0, x1=x1, 
                y0=y0, y1=y1, 
                ytol=1e-3, xtol=1e-3
            )
        
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
        AbstractStirredTankReactor._design(self, size_only=True)
        liquid = bst.Stream(None, thermo=self.thermo)
        liquid.mix_from([i for i in self.ins if i.phase != 'g'], energy_balance=False)
        liquid.copy_thermal_condition(self.outs[0])
        self.effluent_density = rho = liquid.rho
        length = self.get_design_result('Length', 'm') * self.V_wf
        return g * rho * length + 101325 + self.cooler_pressure_drop # Pa
    
    def _design(self):
        AbstractStirredTankReactor._design(self)
        if self.air.isempty(): return
        self.compressor.simulate()
        air_cooler = self.air_cooler
        air_cooler.T = self.T
        air_cooler.dP = self.cooler_pressure_drop
        air_cooler.simulate()
        self.parallel['compressor'] = 1
        self.parallel['air_cooler'] = 1
        # For robust process control, do not include in HXN
        for unit in self.auxiliary_units:
            for hu in unit.heat_utilities: hu.hxn_ok = False

ABR = AeratedBioreactor
