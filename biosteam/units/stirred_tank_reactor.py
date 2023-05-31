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

"""
from .. import Unit
from typing import Optional
from math import ceil
from biosteam.units.design_tools import aeration
from biosteam.units.design_tools import (
    PressureVessel, compute_closed_vessel_turbine_purchase_cost, size_batch
)
from biosteam.units.design_tools.geometry import cylinder_diameter_from_volume
from scipy.constants import g
import flexsolve as flx
import biosteam as bst
from math import pi

__all__ = (
    'StirredTankReactor', 'STR',
    'ContinuousStirredTankReactor', 'CSTR',
    'AeratedBioreactor', 'ABR',
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
                        Yeast    415
    outs...
    [0] CO2
        phase: 'g', T: 305.15 K, P: 101325 Pa
        flow (kmol/hr): Water    10
                        Ethanol  3.73
                        CO2      244
    [1] product
        phase: 'l', T: 305.15 K, P: 101325 Pa
        flow (kmol/hr): Water    6.59e+03
                        Ethanol  240
                        Glucose  4.07
                        Yeast    484
    
    >>> R1.results()
    Continuous fermentation                                    Units                   R1
    Electricity         Power                                     kW              1.4e+03
                        Cost                                  USD/hr                  110
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
                        Turbine (x4)                             USD             5.15e+05
                        Heat exchanger - Floating head (x4)      USD             1.61e+05
                        Recirculation pump - Pump (x4)           USD             5.23e+04
                        Recirculation pump - Motor (x4)          USD             8.38e+03
    Total purchase cost                                          USD             2.13e+06
    Utility cost                                              USD/hr                  180
    
    References
    ----------
    .. [1] Seider, W. D.; Lewin, D. R.; Seader, J. D.; Widagdo, S.; Gani, R.; 
        Ng, M. K. Cost Accounting and Capital Cost Estimation. In Product 
        and Process Design Principles; Wiley, 2017; pp 470.
    
    '''
    auxiliary_unit_names = ('heat_exchanger', 'vacuum_system', 'recirculation_pump')
    
    _units = {**PressureVessel._units,
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
    def liquid_product(self):
        for i in self.outs:
            if i.phase == 'l': 
                return i
        raise AttributeError('no liquid product available')
    
    def __init__(
            self, ID='', ins=None, outs=(), thermo=None, *, 
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
        ):
        
        Unit.__init__(self, ID, ins, outs, thermo)
        self.T = self.T_default if T is None else T
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
        self.recirculation_pump = pump = bst.Pump(None, (None,), (None,), thermo=self.thermo)
        self.splitter = splitter = bst.Splitter(None, pump-0, split=0.5) # Split is updated later
        self.heat_exchanger = bst.HXutility(None, splitter-0, (None,), thermo=self.thermo) 

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
            tau = self._tau
            tau_0 = self.tau_0
            V_wf = self.V_wf
            Design = self.design_results
            V_max = self.V_max
            N_min = 2
            N_max = v_0 / V_wf * (tau + tau_0) / V_max
            f = lambda N: v_0 / N / V_wf * (tau + tau_0) / (1 - 1 / N) - V_max
            if f(N_min) < 0.:
                N = N_min
            else:
                N = flx.IQ_interpolation(f, N_min, N_max,
                                         xtol=0.01, ytol=0.5, checkbounds=False)
                N = ceil(N)
            Design.update(size_batch(v_0, tau, tau_0, N, V_wf))
        else:
            V_total = ins_F_vol * self.tau / self.V_wf
            N = ceil(V_total/self.V_max)
            if N == 0:
                V_reactor = 0
                D = 0
                L = 0
            else:
                V_reactor = V_total / N
                D = cylinder_diameter_from_volume(V_reactor, self.length_to_diameter)
                D *= 3.28084 # Convert from m to ft
                L = D * length_to_diameter
            Design['Reactor volume'] = V_reactor
        
        Design['Residence time'] = self.tau
        Design.update(self._vessel_design(float(P_psi), float(D), float(L)))
        self.vacuum_system = bst.VacuumSystem(self) if P_pascal < 1e5 else None
        duty = self._get_duty()
        self.parallel['self'] = N
        self.parallel['vacuum_system'] = 1 # Not in parallel
        if duty:
            # Note: Flow and duty are rescaled to simulate an individual
            # heat exchanger, then BioSTEAM accounts for number of units in parallel
            # through the `parallel` attribute.
            reactor_duty = duty / N
            dT_hx_loop = self.dT_hx_loop
            reactor_product = self.liquid_product.copy()
            reactor_product.scale(1 / N)
            hx_inlet = reactor_product.copy()
            hx_outlet = hx_inlet.copy()
            hx_outlet.T += (dT_hx_loop if duty > 0. else -dT_hx_loop)
            dH = hx_outlet.H - hx_inlet.H
            recirculation_ratio = reactor_duty / dH # Recirculated flow over net product flow
            hx_inlet.scale(recirculation_ratio)
            hx_outlet.scale(recirculation_ratio)
            self.recirculation_pump.ins[0].mix_from([hx_inlet, reactor_product])
            self.recirculation_pump.simulate()
            self.splitter.split = recirculation_ratio / (1 + recirculation_ratio)
            self.splitter.simulate()
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
            kW = self.kW_per_m3 * volume
            self.add_power_utility(kW)
            if kW > 0:
                hp = kW * 1.34102
                baseline_purchase_costs['Turbine'] = compute_closed_vessel_turbine_purchase_cost(hp)
    
ContinuousStirredTankReactor = CSTR = STR = StirredTankReactor

class AeratedBioreactor(StirredTankReactor):
    """
    Same as StirredTankReactor but includes aeration.
    
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
    ...     'R1', ins=[feed, 'air'], outs=('vent', 'product'), tau=12, V_max=500,
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
    [1] s40
        phase: 'g', T: 305.15 K, P: 101325 Pa
        flow (kmol/hr): O2  914
                        N2  3.44e+03
    [2] air
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow: 0
    outs...
    [0] vent
        phase: 'g', T: 305.15 K, P: 101325 Pa
        flow (kmol/hr): Water  108
                        O2     498
                        N2     3.44e+03
    [1] product
        phase: 'l', T: 305.15 K, P: 101325 Pa
        flow (kmol/hr): Water    5.52e+03
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
    kW_per_m3_default = 0.2955 # Reaction in homogeneous liquid
    
    def __init__(
            self, ID='', ins=None, outs=(), thermo=None,  
            *, reactions, theta_O2=0.5, Q_O2_consumption=None, **kwargs,
        ):
        StirredTankReactor.__init__(self, ID, ins, outs, thermo, **kwargs)
        self.reactions = reactions
        self.theta_O2 = theta_O2 # Concentration of O2 in the liquid as a fraction of saturation.
        self.Q_O2_consumption = Q_O2_consumption # Forced duty per O2 consummed [kJ/kmol].
    
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
    def vent(self):
        return self._outs[0]
    
    @property
    def effluent(self):
        return self._outs[1]
    
    def load_auxiliaries(self):
        super().load_auxiliaries()
        self.compressor = compressor = bst.IsentropicCompressor(None, (None,), (None,), thermo=self.thermo, eta=0.85, P=2 * 101325) 
        self.air_cooler = bst.HXutility(None, compressor-0, (None,), thermo=self.thermo, T=self.T)
        
    def _run_vent(self, vent, effluent):
        vent.receive_vent(effluent, energy_balance=False)
        
    def _run(self):
        air = self.air
        if air is None:
            air = bst.Stream(phase='g', thermo=self.thermo)
            self.ins.insert(1, air)
        feeds = [i for i in self.ins if i.phase != 'g']
        vent, effluent = self.outs
        air.P = vent.P = effluent.P = self.P
        air.T = vent.T = effluent.T = self.T
        vent.empty()
        vent.phase = 'g'
        air.phase = 'g'
        air.empty()
        effluent.mix_from(feeds, energy_balance=False)
        self.run_reactions(effluent)
        self.OUR = OUR = -effluent.get_flow('mol/s', 'O2') # Oxygen uptake rate
        if OUR < 0.:
            self.OUR = 0.
            return
        # print('OUR', format(OUR, '.2f'))
        
        def air_flow_rate_objective(flow):
            air.set_flow([flow, flow * 79. / 21.], 'mol/s', ['O2', 'N2'])
            effluent.imol['O2', 'N2'] = 0.
            effluent.mix_from([effluent, air], energy_balance=False)
            effluent.set_flow(flow - OUR, 'mol/s', 'O2')
            vent.empty()
            self._run_vent(vent, effluent)
            self._design()
            return OUR - self.get_OTR()
        
        y0 = air_flow_rate_objective(OUR)
        if y0 > 0.: # Correlation is not perfect and special cases lead to OTR > OUR
            flx.IQ_interpolation(air_flow_rate_objective, 
                                 x0=OUR, x1=10 * OUR, 
                                 y0=y0, ytol=1e-3, xtol=1e-3)
        
    def run_reactions(self, effluent):
        self.reactions.force_reaction(effluent)
    
    def get_OTR(self):
        V = self.get_design_result('Reactor volume', 'm3') * self.V_wf
        P = 1000 * self.compressor.power_utility.consumption # W
        F = self.air.get_total_flow('m3/s')
        D = self.get_design_result('Diameter', 'm')
        R = 0.5 * D
        A = pi * R * R
        U = F / A
        # print('U', format(U, '.2f'))
        ka_L = aeration.ka_L(P, V, U)
        air_in = self.compressor.outs[0]
        vent = self.vent
        P_O2_air = air_in.get_property('P', 'bar') * air_in.imol['O2'] / air_in.F_mol
        P_O2_vent = vent.get_property('P', 'bar') * vent.imol['O2'] / vent.F_mol
        P_O2_ave = 0.5 * (P_O2_air + P_O2_vent)
        C_O2 = aeration.C_O2_L(self.T, P_O2_ave) # mol / kg
        # print('ka_L', format(ka_L, '.2f'))
        OTR = ka_L * C_O2 * (1. - self.theta_O2) * V * 1000 # mol / s
        # print('OTR', format(OTR, '.2f'))
        return OTR 
        
    def _design(self):
        StirredTankReactor._design(self)
        if self.OUR == 0.: return
        liquid = bst.Stream(None, thermo=self.thermo)
        liquid.mix_from([i for i in self.ins if i.phase != 'g'], energy_balance=False)
        liquid.copy_thermal_condition(self.outs[0])
        rho = liquid.rho
        length = self.get_design_result('Length', 'm') * self.V_wf
        compressor = self.compressor
        compressor.P = g * rho * length + 101325
        compressor.ins[0].copy_like(self.air)
        compressor.simulate()
        air_cooler = self.air_cooler
        air_cooler.T = self.T
        air_cooler.simulate()
        self.parallel['compressor'] = 1
        self.parallel['air_cooler'] = 1
    
ABR = AeratedBioreactor
    
    
    