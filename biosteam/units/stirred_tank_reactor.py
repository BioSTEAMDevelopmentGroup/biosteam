# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
.. contents:: :local:

.. autoclass:: biosteam.units.stirred_tank_reactor.StirredTankReactor

References
----------
.. [1] Seider, W. D.; Lewin, D. R.; Seader, J. D.; Widagdo, S.; Gani, R.; 
    Ng, M. K. Product and Process Design Principles. Wiley 2017.

.. [2] Benz, G. T. Optimize Power Consumption in Aerobic Fermenters. 
    Chem. Eng. Progress 2003, 99 (5), 100–103.

.. [3] Benz, G. T. Bioreactor Design for Chemical Engineers. Chem. Eng.\
    Progress 2011, 21–26.

"""
from .. import Unit
from typing import Optional
from math import ceil
from warnings import catch_warnings
from scipy.optimize import minimize_scalar, minimize, least_squares, fsolve, shgo, differential_evolution
from biosteam.units.design_tools import aeration
from biosteam.units.design_tools import (
    PressureVessel, compute_closed_vessel_turbine_purchase_cost, size_batch
)
from biosteam.units.design_tools.geometry import cylinder_diameter_from_volume
from scipy.constants import g
import flexsolve as flx
import biosteam as bst
from math import pi
import numpy as np

__all__ = (
    'StirredTankReactor', 'STR',
    'ContinuousStirredTankReactor', 'CSTR',
    'AeratedBioreactor', 'ABR',
    'GasFedBioreactor', 'GFB',
)


class StirredTankReactor(PressureVessel, Unit, isabstract=True):
    '''    
    Abstract class for a stirred tank reactor, modeled as a pressure vessel with 
    a given aspect ratio and residence time. A pump-heat exchanger recirculation 
    loop is used to satisfy the duty, if any. By default, a turbine agitator is
    also included if the power usage,`kW_per_m3`, is positive. A vacuum 
    system is also automatically added if the operating pressure is at a vacuum. 

    Parameters
    ----------
    tau :
        Residence time [hr].
    T : 
        Operating temperature [K].
    P : 
        Operating pressure [Pa].
    dT_hx_loop : 
        Maximum change in temperature for the heat exchanger loop. Defaults to 5 K.
    V_wf : 
        Fraction of working volume over total volume. Defaults to 0.8.
    V_max :
        Maximum volume of a reactor [m3]. Defaults to 355.
    kW_per_m3: 
        Power usage of agitator. Defaults to 0.985 [kW / m3] converted from 
        5 hp/1000 gal as in [1]_, for liquid–liquid reaction or extraction.
    vessel_material : 
        Vessel material. Defaults to 'Stainless steel 316'.
    vessel_type : 
        Vessel type. Valid options are 'Horizontal' or 'Vertical'. Defaults to 'Vertical'
    batch :
        Whether to use batch operation mode. If False, operation mode is continuous.
        Defaults to `continuous`.
    tau_0 : 
        Cleaning and unloading time (if batch mode). Defaults to 3 hr.
    
    Notes
    -----
    The recirculation loop takes into account the required flow rate needed to
    reach the maximum temperature change of the heat exchanger, `dT_hx_loop`. 
    Increasing `dT_hx_loop` decreases the required recirculation flow rate and
    therefore decreases pump costs.
    
    When parallel reactors are required, one recirculation loop (each with a
    pump and heat exchanger) is assumed. Although it is possible to use the
    same recirculation loop for all reactors, this conservative assumption allows
    for each reactor to be operated independently from each other.
    
    The capital cost for agitators are not yet included in 
    
    Examples
    --------
    Inherit from ContinuousStirredTankReactor to create a new class that
    simulates the continuous fermentative production of ethanol from sugarcane
    juice:
        
    >>> import biosteam as bst
    >>> class ContinuousFermentation(bst.CSTR):
    ...     _N_ins = 1
    ...     _N_outs = 2
    ...     T_default = 32. + 273.15
    ...     P_default = 101325.
    ...     tau_default = 8.
    ...    
    ...     def _setup(self):
    ...         super()._setup()        
    ...         chemicals = self.chemicals
    ...         self.hydrolysis_reaction = bst.Reaction('Sucrose + Water -> 2Glucose', 'Sucrose', 1.00, chemicals)
    ...         self.fermentation_reaction = bst.Reaction('Glucose -> 2Ethanol + 2CO2',  'Glucose', 0.9, chemicals)
    ...         self.cell_growth_reaction = cell_growth = bst.Reaction('Glucose -> Yeast', 'Glucose', 0.70, chemicals, basis='wt')
    ...     
    ...     def _run(self):
    ...         vent, effluent = self.outs
    ...         effluent.mix_from(self.ins, energy_balance=False)
    ...         self.hydrolysis_reaction(effluent)
    ...         self.fermentation_reaction(effluent)
    ...         self.cell_growth_reaction(effluent)
    ...         effluent.T = vent.T = self.T
    ...         effluent.P = vent.P = self.P
    ...         vent.phase = 'g'
    ...         vent.empty()
    ...         vent.receive_vent(effluent, energy_balance=False)
    ...
    >>> from biorefineries.sugarcane import chemicals
    >>> bst.settings.set_thermo(chemicals)
    >>> feed = bst.Stream('feed',
    ...                   Water=1.20e+05,
    ...                   Glucose=1.89e+03,
    ...                   Sucrose=2.14e+04,
    ...                   DryYeast=1.03e+04,
    ...                   units='kg/hr',
    ...                   T=32+273.15)
    >>> R1 = ContinuousFermentation('R1', ins=feed, outs=('CO2', 'product'))
    >>> R1.simulate()
    >>> R1.show()
    ContinuousFermentation: R1
    ins...
    [0] feed
        phase: 'l', T: 305.15 K, P: 101325 Pa
        flow (kmol/hr): Water    6.66e+03
                        Glucose  10.5
                        Sucrose  62.5
                        Yeast    456
    outs...
    [0] CO2
        phase: 'g', T: 305.15 K, P: 101325 Pa
        flow (kmol/hr): Water    9.95
                        Ethanol  3.71
                        CO2      244
    [1] product
        phase: 'l', T: 305.15 K, P: 101325 Pa
        flow (kmol/hr): Water    6.59e+03
                        Ethanol  240
                        Glucose  4.07
                        Yeast    532
    
    >>> R1.results()
    Continuous fermentation                                    Units                   R1
    Electricity         Power                                     kW             1.04e+03
                        Cost                                  USD/hr                 81.5
    Chilled water       Duty                                   kJ/hr            -1.41e+07
                        Flow                                 kmol/hr             9.42e+03
                        Cost                                  USD/hr                 70.3
    Design              Reactor volume                            m3                  319
                        Residence time                            hr                    8
                        Vessel type                                              Vertical
                        Length                                    ft                 50.5
                        Diameter                                  ft                 16.8
                        Weight                                    lb             5.39e+04
                        Wall thickness                            in                0.363
                        Vessel material                               Stainless steel 316
    Purchase cost       Vertical pressure vessel (x4)            USD             1.18e+06
                        Platform and ladders (x4)                USD             2.12e+05
                        Heat exchanger - Floating head (x4)      USD             1.61e+05
                        Recirculation pump - Pump (x4)           USD              3.9e+04
                        Recirculation pump - Motor (x4)          USD             2.78e+03
                        Agitator - Agitator (x4)                 USD             4.53e+05
    Total purchase cost                                          USD             2.05e+06
    Utility cost                                              USD/hr                  152
    
    '''
    auxiliary_unit_names = (
        'heat_exchanger', 
        'vacuum_system', 
        'recirculation_pump',
        'splitter',
        'agitator',
        'scaler',
    )
    
    _units = {**PressureVessel._units,
              'Batch time': 'hr',
              'Loading time': 'hr',
              'Residence time': 'hr',
              'Total volume': 'm3',
              'Reactor volume': 'm3'}
    
    #: Default operating temperature [K]
    T_default: Optional[float] = None
    
    #: Default operating pressure [K]
    P_default: Optional[float] = None
    
    #: Default residence time [hr]
    tau_default: Optional[float] = None
    
    #: Default maximum change in temperature for the heat exchanger loop.
    dT_hx_loop_default: Optional[float] = 5
    
    #: Default fraction of working volume over total volume.
    V_wf_default: Optional[float] = 0.8
    
    #: Default maximum volume of a reactor in m3.
    V_max_default: Optional[float] = 355
    
    #: Default length to diameter ratio.
    length_to_diameter_default: Optional[float] = 3
    
    #: Default power consumption for agitation [kW/m3].
    kW_per_m3_default: Optional[float] = 0.985
    
    #: Default cleaning and unloading time (hr).
    tau_0_default: Optional[float]  = 3
    
    @property
    def effluent(self):
        return self.outs[-1]
    product = effluent
    
    def _init(
            self, 
            T: Optional[float]=None, 
            P: Optional[float]=None, 
            dT_hx_loop: Optional[float]=None,
            tau: Optional[float]=None,
            V_wf: Optional[float]=None, 
            V_max: Optional[float]=None,
            length_to_diameter: Optional[float]=None, 
            kW_per_m3: Optional[float]=None,
            vessel_material: Optional[str]=None,
            vessel_type: Optional[str]=None,
            batch: Optional[bool]=None,
            tau_0: Optional[float]=None,
            adiabatic: Optional[bool]=None,
        ):
        if adiabatic is None: adiabatic = False
        self.T = self.T_default if (T is None and not adiabatic) else T
        self.adiabatic = adiabatic
        self.P = self.P_default if P is None else P
        self.dT_hx_loop = self.dT_hx_loop_default if dT_hx_loop is None else abs(dT_hx_loop)
        self.tau = self.tau_default if tau is None else tau
        self.V_wf = self.V_wf_default if V_wf is None else V_wf
        self.V_max = self.V_max_default if V_max is None else V_max
        self.length_to_diameter = self.length_to_diameter_default if length_to_diameter is None else length_to_diameter
        self.kW_per_m3 = self.kW_per_m3_default if kW_per_m3 is None else kW_per_m3
        self.vessel_material = 'Stainless steel 316' if vessel_material is None else vessel_material
        self.vessel_type = 'Vertical' if vessel_type is None else vessel_type
        self.tau_0 = self.tau_0_default if tau_0 is None else tau_0
        self.batch = batch
        self.load_auxiliaries()

    def load_auxiliaries(self):
        if self.adiabatic: return
        pump = self.auxiliary('recirculation_pump', bst.Pump)
        if self.batch:
            self.auxiliary('heat_exchanger', bst.HXutility, pump-0) 
        else:
            # Split is updated later
            splitter = self.auxiliary('splitter', bst.Splitter, pump-0, split=0.5)
            self.auxiliary('heat_exchanger', bst.HXutility, splitter-0) 
            self.auxiliary('scaler', bst.Scaler, splitter-1, self.outs[-1])

    def _get_duty(self):
        return self.Hnet

    def _design(self):
        Design = self.design_results
        ins_F_vol = sum([i.F_vol for i in self.ins if i.phase != 'g'])
        P_pascal = (self.P if self.P else self.outs[0].P)
        P_psi = P_pascal * 0.000145038 # Pa to psi
        length_to_diameter = self.length_to_diameter
        
        if self.batch:
            v_0 = ins_F_vol
            tau = self.tau
            tau_0 = self.tau_0
            V_wf = self.V_wf
            Design = self.design_results
            V_max = self.V_max
            N = v_0 / V_max / V_wf * (tau + tau_0) + 1
            if N < 2:
                N = 2
            else:
                N = ceil(N)
            Design.update(size_batch(v_0, tau, tau_0, N, V_wf))
            V_reactor = Design['Reactor volume']
        else:
            V_total = ins_F_vol * self.tau / self.V_wf
            N = ceil(V_total/self.V_max)
            if N == 0:
                V_reactor = 0
            else:
                V_reactor = V_total / N
            Design['Reactor volume'] = V_reactor
            
        D = cylinder_diameter_from_volume(V_reactor, self.length_to_diameter)
        D *= 3.28084 # Convert from m to ft
        L = D * length_to_diameter
        Design['Residence time'] = self.tau
        Design.update(self._vessel_design(float(P_psi), float(D), float(L)))
        self.vacuum_system = bst.VacuumSystem(self) if P_pascal < 1e5 else None
        self.parallel['self'] = N
        self.parallel['vacuum_system'] = 1 # Not in parallel
        if self.adiabatic: return
        duty = self._get_duty()
        if duty:
            # Note: Flow and duty are rescaled to simulate an individual
            # heat exchanger, then BioSTEAM accounts for number of units in parallel
            # through the `parallel` attribute.
            reactor_duty = duty / N
            dT_hx_loop = self.dT_hx_loop
            reactor_product = self.effluent.copy()
            reactor_product.scale(1 / N)
            hx_inlet = reactor_product.copy()
            hx_outlet = hx_inlet.copy()
            hx_outlet.T += (dT_hx_loop if duty > 0. else -dT_hx_loop)
            dH = hx_outlet.H - hx_inlet.H
            recirculation_ratio = reactor_duty / dH # Recirculated flow over net product flow
            hx_inlet.scale(recirculation_ratio)
            hx_outlet.scale(recirculation_ratio)
            if self.batch:
                self.recirculation_pump.ins[0].copy_like(hx_inlet)
                self.recirculation_pump.simulate()
            else:
                self.recirculation_pump.ins[0].mix_from([hx_inlet, reactor_product])
                self.recirculation_pump.simulate()
                self.splitter.split = recirculation_ratio / (1 + recirculation_ratio)
                self.splitter.simulate()
                self.scaler.scale = N
                self.scaler.simulate()
            self.heat_exchanger.T = hx_outlet.T
            self.heat_exchanger.simulate()
            
    def _cost(self):
        Design = self.design_results
        baseline_purchase_costs = self.baseline_purchase_costs
        volume = Design['Reactor volume']
        if volume != 0:
            baseline_purchase_costs.update(
                self._vessel_purchase_cost(
                    Design['Weight'], Design['Diameter'], Design['Length'],
                )
            )
            kW = self.kW_per_m3 * volume * self.V_wf
            if kW > 0: self.agitator = bst.Agitator(kW)
    
ContinuousStirredTankReactor = CSTR = STR = StirredTankReactor


class AeratedBioreactor(StirredTankReactor):
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
        flow (kmol/hr): O2  730
                        N2  2.75e+03
    outs...
    [0] vent
        phase: 'g', T: 305.15 K, P: 101325 Pa
        flow (kmol/hr): Water  109
                        CO2    416
                        O2     314
                        N2     2.75e+03
    [1] product
        phase: 'l', T: 305.15 K, P: 101325 Pa
        flow (kmol/hr): Water    6.97e+03
                        Glucose  69.4
    
    """
    _N_ins = 2
    _N_outs = 2
    _ins_size_is_fixed = False
    auxiliary_unit_names = (
        'compressor',
        'air_cooler',
        *StirredTankReactor.auxiliary_unit_names
    )
    T_default = 273.15 + 32 
    P_default = 101325
    kW_per_m3_default = 0.2955 # Reaction in homogeneous liquid; reference [1]
    
    def _init(
            self, reactions, theta_O2=0.5, Q_O2_consumption=None,
            optimize_power=None, kLa_coefficients=None, **kwargs,
        ):
        StirredTankReactor._init(self, **kwargs)
        self.reactions = reactions
        self.theta_O2 = theta_O2 # Average concentration of O2 in the liquid as a fraction of saturation.
        self.Q_O2_consumption = Q_O2_consumption # Forced duty per O2 consummed [kJ/kmol].
        self.optimize_power = True if optimize_power is None else optimize_power
        self.kLa_coefficients = kLa_coefficients
    
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
    def cooled_compressed_air(self):
        return self.air_cooler.outs[0]
    
    @property
    def vent(self):
        return self._outs[0]
    
    def load_auxiliaries(self):
        super().load_auxiliaries()
        compressor = self.auxiliary(
            'compressor', bst.IsentropicCompressor, self.air, eta=0.85, P=2 * 101325
        )
        self.auxiliary(
            'air_cooler', bst.HXutility, compressor-0, T=self.T
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
        air_cc = self.cooled_compressed_air
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
        air_in = self.cooled_compressed_air
        N_reactors = self.parallel['self']
        operating_time = self.tau / self.design_results.get('Batch time', 1.)
        V = self.get_design_result('Reactor volume', 'm3') * self.V_wf
        D = self.get_design_result('Diameter', 'm')
        F = air_in.get_total_flow('m3/s') / N_reactors / operating_time
        R = 0.5 * D
        A = pi * R * R
        self.superficial_gas_flow = U = F / A # m / s 
        vent = self.vent
        P_O2_air = air_in.get_property('P', 'bar') * air_in.imol['O2'] / air_in.F_mol
        P_O2_vent = 0. if vent.isempty() else vent.get_property('P', 'bar') * vent.imol['O2'] / vent.F_mol
        C_O2_sat_air = aeration.C_O2_L(self.T, P_O2_air) # mol / kg
        C_O2_sat_vent = aeration.C_O2_L(self.T, P_O2_vent) # mol / kg
        theta_O2 = self.theta_O2
        LMDF = aeration.log_mean_driving_force(C_O2_sat_vent, C_O2_sat_air, theta_O2 * C_O2_sat_vent, theta_O2 * C_O2_sat_air)
        kLa = OUR / (LMDF * V * self.effluent_density * N_reactors * operating_time)
        P = aeration.P_at_kLa(kLa, V, U, self.kLa_coefficients)
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
        P = 1000 * self.kW_per_m3 * V # W
        air_in = self.cooled_compressed_air
        D = self.get_design_result('Diameter', 'm')
        F = air_in.get_total_flow('m3/s') / N_reactors / operating_time
        R = 0.5 * D
        A = pi * R * R
        self.superficial_gas_flow = U = F / A # m / s 
        kLa = aeration.kLa(P, V, U, self.kLa_coefficients) # 1 / s 
        air_in = self.cooled_compressed_air
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
        StirredTankReactor._design(self)
        liquid = bst.Stream(None, thermo=self.thermo)
        liquid.mix_from([i for i in self.ins if i.phase != 'g'], energy_balance=False)
        liquid.copy_thermal_condition(self.outs[0])
        self.effluent_density = rho = liquid.rho
        length = self.get_design_result('Length', 'm') * self.V_wf
        return g * rho * length + 101325 # Pa
    
    def _design(self):
        StirredTankReactor._design(self)
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


class GasFedBioreactor(StirredTankReactor):
    """
    Same as StirredTankReactor but includes multiple gas feeds. The agitator power may
    vary to minimize the total power requirement of both the compressor and agitator
    yet achieve the required oxygen transfer rate.
    
    # Examples
    # --------
    # >>> import biosteam as bst
    # >>> bst.settings.set_thermo(['H2', 'CO2', 'N2', 'O2', 'H2O', 'AceticAcid'])
    # >>> media = bst.Stream(ID='media', H2O=10000, units='kg/hr')
    # >>> H2 = bst.Stream(ID='H2', phase='g')
    # >>> fluegas = bst.Stream(ID='fluegas', phase='g')
    # >>> recycle = bst.Stream(ID='recycle', phase='g', N2=70, CO2=23, H2O=3, O2=4, total_flow=10, units='kg/hr')
    # >>> # Model acetic acid production from H2 and CO2
    # >>> rxn = bst.Rxn('H2 + CO2 -> AceticAcid + H2O', reactant='H2', correct_atomic_balance=True) 
    # >>> brxn = rxn.backwards(reactant='AceticAcid')
    # >>> R1 = bst.GasFedBioreactor(
    # ...     'R1', ins=[media, H2, fluegas, recycle], outs=('vent', 'product'), tau=68, V_max=500,
    # ...     reactions=rxn, backward_reactions=brxn,
    # ...     feed_gas_compositions={
    # ...         1: dict(H2=100, units='kg/hr'),
    # ...         2: dict(N2=70, CO2=25, H2O=3, O2=2, units='kg/hr'),
    # ...     },
    # ...     gas_substrates=('H2', 'CO2'),
    # ...     titer={'AceticAcid': 5},
    # ...     mixins={2: [3]}, # Recycle gets mixed with fluegas
    # ...     optimize_power=False,
    # ...     kW_per_m3=0.,
    # ... )
    # >>> R1.simulate()
    # >>> R1.show()
    # GasFedBioreactor: R1
    # ins...
    # [0] media  
    #     phase: 'l', T: 298.15 K, P: 101325 Pa
    #     flow (kmol/hr): H2O  555
    # [1] H2  
    #     phase: 'g', T: 298.15 K, P: 101325 Pa
    #     flow (kmol/hr): H2  49.6
    # [2] fluegas  
    #     phase: 'g', T: 298.15 K, P: 101325 Pa
    #     flow (kmol/hr): CO2  2.01
    #                     N2   8.83
    #                     O2   0.221
    #                     H2O  0.589
    # [3] recycle  
    #     phase: 'g', T: 298.15 K, P: 101325 Pa
    #     flow (kmol/hr): CO2  0.0523
    #                     N2   0.25
    #                     O2   0.0125
    #                     H2O  0.0167
    # outs...
    # [0] vent  
    #     phase: 'g', T: 305.15 K, P: 101325 Pa
    #     flow (kmol/hr): H2          46.3
    #                     CO2         0.387
    #                     N2          9.08
    #                     O2          0.233
    #                     H2O         2.39
    #                     AceticAcid  0.00231
    # [1] product  
    #     phase: 'l', T: 305.15 K, P: 101325 Pa
    #     flow (kmol/hr): H2O         555
    #                     AceticAcid  0.834
    
    """
    _N_ins = 2
    _N_outs = 2
    _ins_size_is_fixed = False
    auxiliary_unit_names = (
        'mixers',
        'sparger',
        'compressors',
        'gas_coolers',
        *StirredTankReactor.auxiliary_unit_names
    )
    T_default = 273.15 + 32 
    P_default = 101325
    kW_per_m3_default = 0.2955 # Reaction in homogeneous liquid
    
    def _init(self, 
            reactions, gas_substrates, titer, backward_reactions, 
            feed_gas_compositions, 
            theta=0.5, Q_consumption=None,
            optimize_power=None, kLa_coefficients=None, 
            mixins=None,
            **kwargs,
        ):
        self.reactions = reactions
        self.backward_reactions = backward_reactions
        self.theta = theta # Average concentration of gas substrate in the liquid as a fraction of saturation.
        self.Q_consumption = Q_consumption # Forced duty per gas substrate consummed [kJ/kmol].
        self.kLa_coefficients = kLa_coefficients
        self.optimize_power = True if optimize_power is None else optimize_power
        self.feed_gas_compositions = feed_gas_compositions # dict[int, dict] Feed index and composition pairs.
        self.gas_substrates = gas_substrates
        self.titer = titer # dict[str, float] g / L
        self.mixins = {} if mixins is None else mixins # dict[int, tuple[int]] Pairs of variable feed gas index and inlets that will be mixed.
        StirredTankReactor._init(self, **kwargs)
    
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
        return [self.ins[i] for i in self.feed_gas_compositions]
    
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
        self.mixers = []
        mixins = self.mixins
        for i in self.feed_gas_compositions:
            if i in mixins:
                other_ins = [self.ins[j] for j in mixins[i]]
                mixer = self.auxiliary(
                    'mixers', bst.Mixer, (self.ins[i], *other_ins),
                )
                inlet = mixer-0
            else:
                inlet = self.ins[i]
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
        F_vol = effluent.F_vol # m3 / hr
        produced = bst.Stream(None, thermo=self.thermo)
        for ID, concentration in self.titer.items():
            produced.imass[ID] = F_vol * concentration
        consumed = produced.copy()
        self.backward_reactions.force_reaction(consumed)
        return consumed.get_flow('mol/s', self.gas_substrates), consumed, produced
    
    def _run(self):
        variable_gas_feeds = self.variable_gas_feeds
        for i in variable_gas_feeds: i.phase = 'g'
        liquid_feeds = [i for i in self.ins if i.phase != 'g']
        vent, effluent = self.outs
        vent.P = effluent.P = self.P
        vent.T = effluent.T = self.T
        vent.empty()
        vent.phase = 'g'
        effluent.mix_from(liquid_feeds, energy_balance=False)
        effluent_liquid_data = effluent.get_data()
        SURs, s_consumed, s_produced = self.get_SURs(effluent) # Gas substrate uptake rate [mol / s]
        if (SURs <= 1e-2).all():
            effluent.imol[self.gas_substrates] = 0.
            self._run_vent(vent, effluent)
            return
        T = self.T
        P = self._inlet_gas_pressure()
        for i in self.compressors: i.P = P
        for i in self.gas_coolers: i.T = T
        x_substrates = []
        for (i, dct), ID in zip(self.feed_gas_compositions.items(), self.gas_substrates):
            gas = self.ins[i]
            gas.reset_flow(**dct)
            x_substrates.append(gas.get_molar_fraction(ID))
        index = range(len(self.gas_substrates))
        
        def run_auxiliaries():
            for i in self.mixers: i.simulate()
            for i in self.compressors: i.simulate()
            for i in self.gas_coolers: i.simulate()
            self.sparger.simulate()
        
        def load_flow_rates(F_feeds):
            for i in index:
                gas = variable_gas_feeds[i]
                gas.set_total_flow(F_feeds[i], 'mol/s')
            run_auxiliaries()
            effluent.set_data(effluent_liquid_data)
            effluent.mix_from([self.sparged_gas, -s_consumed, s_produced, *liquid_feeds], energy_balance=False)
            vent.empty()
            self._run_vent(vent, effluent)
        
        baseline_feed = bst.Stream.sum(self.normal_gas_feeds, energy_balance=False)
        baseline_flows = baseline_feed.get_flow('mol/s', self.gas_substrates)
        bounds = np.array([[max(SURs[i] - baseline_flows[i], 0), 5 * SURs[i]] for i in index])
        if self.optimize_power:
            def total_power_at_substrate_flow(F_substrates):
                load_flow_rates(F_substrates)
                total_power = self._solve_total_power(SURs)
                return total_power
            
            f = total_power_at_substrate_flow
            with catch_warnings(action='ignore'):
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
                return (diff * diff).sum()
            
            f = gas_flow_rate_objective
            bounds = bounds.T
            with catch_warnings(action='ignore'):
                results = least_squares(f, 1.2 * SURs, bounds=bounds, ftol=SURs.min() * 1e-6)
                load_flow_rates(results.x / x_substrates)
        
        # self.show()
        # breakpoint()
        
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
            Ps.append(aeration.P_at_kLa(kLa, V, U, self.kLa_coefficients))
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
        P = 1000 * self.kW_per_m3 * V # W
        gas_in = self.sparged_gas
        D = self.get_design_result('Diameter', 'm')
        F = gas_in.get_total_flow('m3/s') / N_reactors / operating_time
        R = 0.5 * D
        A = pi * R * R
        self.superficial_gas_flow = U = F / A # m / s 
        kLa = aeration.kLa(P, V, U, self.kLa_coefficients) # 1 / s 
        vent = self.vent
        P_gas = gas_in.get_property('P', 'bar')
        P_vent = vent.get_property('P', 'bar')
        STRs = []
        for ID in self.gas_substrates:
            Py_gas = P_gas * gas_in.imol[ID] / gas_in.F_mol
            Py_vent = P_vent * vent.imol[ID] / vent.F_mol
            C_sat_gas = aeration.C_L(self.T, Py_gas, ID) # mol / kg
            C_sat_vent = aeration.C_L(self.T, Py_vent, ID) # mol / kg
            theta = self.theta
            LMDF = aeration.log_mean_driving_force(C_sat_vent, C_sat_gas, theta * C_sat_vent, theta * C_sat_gas)
            STRs.append(
                kLa * LMDF * self.effluent_density * V * N_reactors * operating_time # mol / s
            )
        return np.array(STRs)
        
    def _inlet_gas_pressure(self):
        StirredTankReactor._design(self)
        liquid = bst.Stream(None, thermo=self.thermo)
        liquid.mix_from([i for i in self.ins if i.phase != 'g'], energy_balance=False)
        liquid.copy_thermal_condition(self.outs[0])
        self.effluent_density = rho = liquid.rho
        length = self.get_design_result('Length', 'm') * self.V_wf
        return g * rho * length + 101325 # Pa
    
    def _design(self):
        StirredTankReactor._design(self)
        self.parallel['sparger'] = 1
        self.parallel['mixers'] = 1
        self.parallel['compressors'] = 1
        self.parallel['gas_coolers'] = 1
        # For robust process control, do not include in HXN
        for unit in self.auxiliary_units:
            for hu in unit.heat_utilities: hu.hxn_ok = False
    
GFB = GasFedBioreactor
ABR = AeratedBioreactor
    
    
    
