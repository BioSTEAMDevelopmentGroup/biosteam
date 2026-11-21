# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
.. contents:: :local: 

.. autoclass:: biosteam.units.abstract_stirred_tank_reactor.AbstractStirredTankReactor

References
----------
.. [1] Seider, W. D.; Lewin, D. R.; Seader, J. D.; Widagdo, S.; Gani, R.; 
    Ng, M. K. Product and Process Design Principles. Wiley 2017.

"""
from .. import Unit
from .._unit import design_input, update_design_input
from typing import Optional
from math import ceil
from biosteam.units.design_tools.geometry import cylinder_diameter_from_volume
import biosteam as bst
from biosteam.units.design_tools import (
    PressureVessel, size_batch
)

__all__ = (
    'AbstractStirredTankReactor', 
    'ContinuousStirredTankReactor',
    'CSTR', 'STR',
)

class AbstractStirredTankReactor(PressureVessel, Unit, isabstract=True):
    '''
    The reactor is designed as a pressure vessel with a given aspect ratio and 
    residence time. A pump-heat exchanger recirculation loop can be used to satisfy 
    the duty, if any. By default, a turbine agitator is also included if the 
    power usage, `kW_per_m3`, is positive. A vacuum system is also 
    automatically added if the operating pressure is at a vacuum. 

    Parameters
    ----------
    tau :
        Residence time [hr].
    T : 
        Operating temperature [K].
    P : 
        Operating pressure [Pa].
    V_wf : 
        Fraction of working volume over total volume. Defaults to 0.8.
    V_max :
        Maximum volume of a reactor [m3]. Defaults to 355.
    length_to_diameter :
        Length to diameter ratio of bioreactor.
    kW_per_m3 : 
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
    N :
        Number of reactors.
    heat_exchanger_configuration : 
        What kind of heat exchanger to default to (if any). Valid options include 
        'Jacketed', 'recirculation loop', and 'Internal'. Defaults to 'recirculation loop'.
    dT_hx_loop : 
        Maximum change in temperature for the heat exchanger loop. Defaults to 5 K.
    jacket_annular_diameter :
        Annular diameter of heat exchanger jacket to vessel [m]. Defaults to 0.1 m.
    loading_time :
        Loading time of batch reactor. If not given, it will assume each vessel is constantly
        being filled.
    N_vessels :
        Number of vessels.

    Notes
    -----
    The heat exchanger configuration can be one of the following:

    * 'Recirculation loop': 
        The recirculation loop takes into account the required flow rate needed to
        reach the maximum temperature change of the heat exchanger, `dT_hx_loop`. 
        Increasing `dT_hx_loop` decreases the required recirculation flow rate and
        therefore decreases pump costs.
        
        When parallel reactors are required, one recirculation loop (each with a
        pump and heat exchanger) is assumed. Although it is possible to use the
        same recirculation loop for all reactors, this conservative assumption allows
        for each reactor to be operated independently from each other.

    * 'Jacketed':
        The jacket does not account for the heat transfer area requirement. 
        It simply assumes that a full jacket can provide the necessary heat transfer
        area to meet the duty requirement. A heuristic annular diameter is assumed
        through `jacket_annular_diameter` (which can be adjusted by the user).
        The temperature at the wall is assumed to be the operating temperature.
        The weight of the jacket is added to the weight of the vessel and the
        cost is compounded together as a jacketed vessel.
        
    * 'Internal':
        The internal is costed as an ordinary helical tube heat exchanger
        with the added assumption that the temperature at the wall is the 
        operating temperature. This method is still not implemented in BioSTEAM
        yet.


    Examples
    --------
    Inherit from AbstractStirredTankReactor to create a new class that
    simulates the continuous fermentative production of ethanol from sugarcane
    juice:
        
    >>> import biosteam as bst
    >>> class ContinuousFermentation(bst.AbstractStirredTankReactor):
    ...     _N_ins = 1
    ...     _N_outs = 2
    ...     T = bst.update_design_input(default=32. + 273.15)
    ...     P = bst.update_design_input(default=101325.)
    ...     tau = bst.update_design_input(default=8.)
    ...     
    ...     def _run(self):
    ...         vent, effluent = self.outs
    ...         effluent.mix_from(self.ins, energy_balance=False)
    ...         self.reactions(effluent)
    ...         effluent.T = vent.T = self.T
    ...         effluent.P = vent.P = self.P
    ...         vent.phase = 'g'
    ...         vent.empty()
    ...         vent.receive_vent(effluent, energy_balance=False)
    ...
    >>> from biorefineries.sugarcane import chemicals
    >>> bst.settings.set_thermo(chemicals)
    >>> feed = bst.Stream('feed',
    ...     Water=1.20e+05,
    ...     Glucose=1.89e+03,
    ...     Sucrose=2.14e+04,
    ...     DryYeast=1.03e+04,
    ...     units='kg/hr',
    ...     T=32+273.15
    ... )
    >>> R1 = ContinuousFermentation('R1',
    ...     ins=feed, outs=('CO2', 'product'), 
    ...     reactions = bst.ReactionSystem(
    ...         bst.Reaction('Sucrose + Water -> 2Glucose', 'Sucrose', 1.00),  # Hydrolysis
    ...         bst.Reaction('Glucose -> 2Ethanol + 2CO2',  'Glucose', 0.9),   # Production
    ...         bst.Reaction('Glucose -> Yeast', 'Glucose', 0.70, basis='wt'), # Growth
    ...         basis='mol',
    ...     )
    ... )
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
    Continuous fermentation                                         Units                   R1
    Electricity              Power                                     kW             1.04e+03
                             Cost                                  USD/hr                 81.5
    Chilled water            Duty                                   kJ/hr            -1.41e+07
                             Flow                                 kmol/hr             9.42e+03
                             Cost                                  USD/hr                 70.3
    Design inputs            Vessel type                                              Vertical
                             Vessel material                               Stainless steel 316
                             Operating temperature                      K                  305
                             Operating pressure                        Pa             1.01e+05
                             Residence time                            hr                    8
                             Working volume fraction                                       0.8
                             Maximum reactor volume                    m3                  355
                             Length to diameter ratio                                        3
                             Agitator power                         kW/m3                0.985
                             Cleaning and unloading time               hr                    3
                             Residence time phase                                       Liquid
                             Batch                                                       False
                             Heat exchanger configuration                   Recirculation loop
                             HX-loop temperature change                 K                    5
    Design results           Reactor volume                            m3                  319
                             Number of reactors                                              4
                             Residence time                            hr                    8
                             Length                                    ft                 50.5
                             Diameter                                  ft                 16.8
                             Weight                                    lb             5.39e+04
                             Wall thickness                            in                0.363
    Purchase cost            Pressure vessel (x4)                     USD             1.18e+06
                             Platform and ladders (x4)                USD             2.12e+05
                             Heat exchanger - Floating head (x4)      USD             1.61e+05
                             Recirculation pump - Pump (x4)           USD              3.9e+04
                             Recirculation pump - Motor (x4)          USD             2.78e+03
                             Agitator - Agitator (x4)                 USD             4.53e+05
    Total purchase cost                                               USD             2.05e+06
    Installed equipment cost                                          USD             4.25e+06
    Utility cost                                                   USD/hr                  152
    
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
    
    #: Basis of residence time. Defaults to inlet flow rates
    tau_basis: str = 'ins' # Alternatively, 'outs'
    T = design_input('Operating temperature', 'K', lb=0, default=None)
    P = design_input('Operating pressure', 'Pa', lb=0, default=None)
    tau = design_input('Residence time', 'hr', lb=0, default=None)
    V_wf = design_input('Working volume fraction', lb=0, ub=1, default=0.8)
    V_max = design_input('Maximum reactor volume', 'm3', lb=0, default=355)
    length_to_diameter = design_input('Length to diameter ratio', lb=0, default=3)
    kW_per_m3 = design_input('Agitator power', 'kW/m3', lb=0, default=0.985)
    tau_0 = design_input('Cleaning and unloading time', 'hr', lb=0, default=3)
    N: int = design_input('Number of reactors', default=None)
    vessel_material = update_design_input(default='Stainless steel 316')
    vessel_type = update_design_input(default='Vertical')
    tau_phase = design_input(
        'Residence time phase', default='Liquid', 
        catalog={'Liquid', 'Gas', 'Solid'},
        doc="Phase accounted for in residence time. "
            "If None, residence time applies to all phases. "
            "Defaults to 'Liquid'.",
    )
    loading_time = design_input(
        'Loading time', 'hr', lb=0, default=None,
        doc="Loading time of batch reactor [hr]. If not given, it will assume "
            "each vessel is constantly being filled.",
    )
    batch = design_input(
        default=False, 
        doc="Whether to use batch operation mode. If False, operation mode "
            "is continuous. Defaults to False.",
    )
    
    @design_input(default=False)
    def adiabatic(self):
        return self._adiabatic or None # If not adiabatic it will not appear in design inputs
    @adiabatic.setter
    def adiabatic(self, adiabatic):
        self._adiabatic = adiabatic

    @design_input(
        catalog={'Jacketed', 'Recirculation loop', 'Internal'}, 
        default='Recirculation loop',
    )
    def heat_exchanger_configuration(self):
        # If adiabatic, hx configuration is not active and will not appear in design inputs
        if not self.adiabatic: return self._heat_exchanger_configuration
    @heat_exchanger_configuration.setter
    def heat_exchanger_configuration(self, configuration):
        if configuration == 'Internal':
            raise AttributeError(
                "'Internal' heat exchanger configuration not implemented "
                " in BioSTEAM yet"
            )
        else:
            self._heat_exchanger_configuration = configuration

    @design_input(
        'HX-loop temperature change', 'K', lb=0, default=5,
    )
    def dT_hx_loop(self):
        """Maximum change in temperature for the heat exchanger loop [K].
        Defaults to 5 K."""
        dT = self._dT_hx_loop
        if self.heat_exchanger_configuration == 'Recirculation loop': return dT
    @dT_hx_loop.setter
    def dT_hx_loop(self, dT):
        self._dT_hx_loop = dT
    
    @design_input(
        'Jacket annular diameter', 'm', lb=0, default=0.05,
    )
    def jacket_annular_diameter(self):
        d = self._jacket_annular_diameter
        if self.heat_exchanger_configuration == 'Jacketed': return d
    @jacket_annular_diameter.setter
    def jacket_annular_diameter(self, diameter):
        self._jacket_annular_diameter = diameter

    @property
    def effluent(self):
        return self.outs[-1]
    product = effluent
    
    def _init(
            self,
            reactions: Optional[bst.ReactionSystem]=None,
        ):
        self.reactions = reactions
    
    def _setup(self):
        super()._setup()
        self.load_auxiliaries()
    
    def load_auxiliaries(self):
        if self.adiabatic or not self.heat_exchanger_configuration == 'Recirculation loop': 
            self.recirculation_loop = None
            self.heat_exchanger = None
            self.scaler = None
            self.splitter = None
            return
        pump = self.auxiliary('recirculation_pump', bst.Pump)
        if self.batch:
            self.auxiliary('heat_exchanger', bst.HXutility, pump-0) 
            self.recirculation_loop = None
            self.scaler = None
            self.splitter = None
        else:
            # Split is updated later
            splitter = self.auxiliary('splitter', bst.Splitter, pump-0, split=0.5)
            self.auxiliary('heat_exchanger', bst.HXutility, splitter-0) 
            self.auxiliary('scaler', bst.Scaler, splitter-1, self.outs[-1])

    def _get_duty(self): # User can replace
        return self.Hnet

    def _design(self, size_only=False):
        Design = self.design_results
        tau_basis = self.tau_basis
        streams = getattr(self, tau_basis)
        tau_phase = self.tau_phase
        if tau_phase: 
            if tau_phase == 'Liquid':
                phase = 'l'
            elif tau_phase == 'Gas': 
                phase = 'g'
            elif tau_phase == 'Solid': 
                phase = 's'
            else:
                raise RuntimeError('invalid `tau_phase`')
            F_vol = sum([i.F_vol for i in streams if i.phase == phase]) 
        else:
            F_vol = sum([i.F_vol for i in streams])
        P_pascal = (self.P if self.P else self.outs[0].P)
        P_psi = P_pascal * 0.000145038 # Pa to psi
        length_to_diameter = self.length_to_diameter
        if self.batch:
            v_0 = F_vol
            tau = self.tau
            tau_0 = self.tau_0
            V_wf = self.V_wf
            Design = self.design_results
            Design.update(
                size_batch(
                    v_0, tau, tau_0, V_wf, 
                    self.V_max, self.N, self.loading_time
                )
            )
            V_reactor = Design['Reactor volume']
            N = Design['Number of reactors']
        else:
            if self.N is None:
                V_total = F_vol * self.tau / self.V_wf
                N = ceil(V_total/self.V_max)
                if N == 0:
                    V_reactor = 0
                else:
                    V_reactor = V_total / N
            else:
                V_total = F_vol * self.tau / self.V_wf
                V_reactor = V_total / self.N
            Design['Reactor volume'] = V_reactor
            Design['Number of reactors'] = N
        D = cylinder_diameter_from_volume(V_reactor, self.length_to_diameter)
        D *= 3.28084 # Convert from m to ft
        L = D * length_to_diameter
        Design['Residence time'] = self.tau
        hx_config = self.heat_exchanger_configuration
        if hx_config == 'Jacketed':
            annular_diameter = self.jacket_annular_diameter * 3.28084 # Annular diameter [ft]
            dct = self._vessel_design(P_psi, D, L, annular_diameter, check=not size_only)
        else:
            dct = self._vessel_design(P_psi, D, L, check=not size_only)
        Design.update(dct)
        self.vacuum_system = bst.VacuumSystem(self) if P_pascal < 1e5 else None
        self.parallel['self'] = N
        self.parallel['vacuum_system'] = 1 # Not in parallel
        if self.adiabatic or size_only: return
        duty = self._get_duty()
        if not duty: return
        if hx_config == 'Recirculation loop':
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
        elif hx_config == 'Jacketed':
            self.add_heat_utility(duty / N, self.effluent.T)
        else:
            raise RuntimeError(
                f"heat exchanger configuration {hx_config!r} not yet implemented"
            )
            
    def _cost(self):
        Design = self.design_results
        baseline_purchase_costs = self.baseline_purchase_costs
        volume = Design['Reactor volume']
        if volume != 0:
            if self.heat_exchanger_configuration == 'Jacketed':
                baseline_purchase_costs.update(
                    self._vessel_purchase_cost(
                        Design['Weight'], 
                        Design['Diameter'], 
                        Design['Length'], 
                        jacketed=True,
                    )
                )
            else:
                baseline_purchase_costs.update(
                    self._vessel_purchase_cost(
                        Design['Weight'], Design['Diameter'], Design['Length'], 
                    )
                )
            kW = self.kW_per_m3 * volume * self.V_wf
            if kW > 0: self.agitator = bst.Agitator(kW)

ContinuousStirredTankReactor = CSTR = STR = StirredTankReactor = AbstractStirredTankReactor
