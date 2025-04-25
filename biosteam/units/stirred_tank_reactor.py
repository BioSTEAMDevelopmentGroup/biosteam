# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
.. contents:: :local:

.. autoclass:: biosteam.units.stirred_tank_reactor.AbstractStirredTankReactor

References
----------
.. [1] Seider, W. D.; Lewin, D. R.; Seader, J. D.; Widagdo, S.; Gani, R.; 
    Ng, M. K. Product and Process Design Principles. Wiley 2017.

"""
from .. import Unit
from typing import Optional
from math import ceil
from biosteam.units.design_tools.geometry import cylinder_diameter_from_volume
import biosteam as bst
from biosteam.units.design_tools import (
    PressureVessel, compute_closed_vessel_turbine_purchase_cost, size_batch
)

__all__ = (
    'AbstractStirredTankReactor', 
    'StirredTankReactor', 'STR',
    'ContinuousStirredTankReactor', 'CSTR',
)


class AbstractStirredTankReactor(PressureVessel, Unit, isabstract=True):
    '''    
    Abstract class for a stirred tank reactor, modeled as a pressure vessel with 
    a given aspect ratio and residence time. A pump-heat exchanger recirculation 
    loop is used to satisfy the duty, if any. By default, a turbine agitator is
    also included if the power usage,`kW_per_m3` , is positive. A vacuum 
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
    N_reactors :
        Number of reactors.
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
    Inherit from AbstractStirredTankReactor to create a new class that
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
    
    #: Whether to default operation in batch mode or continuous
    batch_default = False
    
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
        self.batch = self.batch_default if batch is None else batch
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

    def _design(self, size_only=False):
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
        self.N_reactors = N
        D = cylinder_diameter_from_volume(V_reactor, self.length_to_diameter)
        D *= 3.28084 # Convert from m to ft
        L = D * length_to_diameter
        Design['Residence time'] = self.tau
        Design.update(self._vessel_design(float(P_psi), float(D), float(L)))
        self.vacuum_system = bst.VacuumSystem(self) if P_pascal < 1e5 else None
        self.parallel['self'] = N
        self.parallel['vacuum_system'] = 1 # Not in parallel
        if self.adiabatic or size_only: return
        duty = self._get_duty()
        if duty:
            # Note: Flow and duty are rescaled to simulate an individual
            # heat exchanger, then BioSTEAM accounts for number of units in parallel
            # through the `parallel` attribute.
            if N == 0:
                raise RuntimeError(
                    'stirred tank reactor has heating/cooling duty but no '
                    'liquid to recirculate to heat exchanger'
                )
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

ContinuousStirredTankReactor = CSTR = STR = StirredTankReactor = AbstractStirredTankReactor
