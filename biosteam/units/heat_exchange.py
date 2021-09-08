# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
This module contains heat exchanger unit operations.

.. contents:: :local:
    
Unit operations
---------------
.. autoclass:: biosteam.units.heat_exchange.HX
.. autoclass:: biosteam.units.heat_exchange.HXutility
.. autoclass:: biosteam.units.heat_exchange.HXprocess 

"""
from .. import Unit
from .._graphics import utility_heat_exchanger_graphics, process_heat_exchanger_graphics
from .design_tools.specification_factors import (
    shell_and_tube_material_factor_coefficients,
    compute_shell_and_tube_material_factor)
from .design_tools import heat_transfer as ht
import numpy as np
import biosteam as bst
from math import exp, log as ln

__all__ = ('HX', 'HXutility', 'HXprocess')

# Lenght factor 
x = np.array((8, 13, 16, 20)) 
y = np.array((1.25, 1.12,1.05,1))
p2 = np.polyfit(x, y, 2)

# %% Purchase price

def compute_floating_head_purchase_price(A, CE):
    return exp(12.0310 - 0.8709*ln(A) + 0.09005 * ln(A)**2)*CE/567

def compute_fixed_head_purchase_price(A, CE):
    return exp(11.4185 - 0.9228*ln(A) + 0.09861 * ln(A)**2)*CE/567

def compute_u_tube_purchase_price(A, CE):
    return exp(11.5510 - 0.9186*ln(A) + 0.09790 * ln(A)**2)*CE/567

def compute_kettle_vaporizer_purchase_price(A, CE):
    return exp(12.3310 - 0.8709*ln(A) + 0.09005 * ln(A)**2)*CE/567

def compute_double_pipe_purchase_price(A, CE):
    return exp( 7.2718 + 0.16*ln(A))*CE/567

# Purchase price
Cb_dict = {'Floating head': compute_floating_head_purchase_price,
           'Fixed head': compute_fixed_head_purchase_price,
           'U tube': compute_u_tube_purchase_price,
           'Kettle vaporizer': compute_kettle_vaporizer_purchase_price,
           'Double pipe': compute_double_pipe_purchase_price}

# %% Classes

class HX(Unit, isabstract=True):
    """
    Abstract class for counter current heat exchanger.

    **Abstract methods**
    
    get_streams()
        Should return two inlet streams and two outlet streams that exchange
        heat.

    """
    line = 'Heat exchanger'
    _units = {'Area': 'ft^2',
              'Overall heat transfer coefficient': 'kW/m^2/K',
              'Log-mean temperature difference': 'K',
              'Tube side pressure drop': 'psi',
              'Shell side pressure drop': 'psi',
              'Operating pressure': 'psi',
              'Total tube length': 'ft'}
    _N_ins = 1
    _N_outs = 1
    _N_heat_utilities = 1
    _F_BM_default = {'Double pipe': 1.8,
                     'Floating head': 3.17,
                     'Fixed head': 3.17,
                     'U tube': 3.17,
                     'Kettle vaporizer': 3.17}
    
    @property
    def material(self):
        """Default 'Carbon steel/carbon steel'"""
        return self.material
    @material.setter
    def material(self, material):
        try:
            self._F_Mab = shell_and_tube_material_factor_coefficients[material]
        except KeyError:
            raise AttributeError("material must be one of the following: "
                                 f"{', '.join(shell_and_tube_material_factor_coefficients)}")
        self._material = material  
    
    @property
    def heat_exchanger_type(self):
        """[str] Heat exchanger type. Purchase cost depends on this selection."""
        return self._heat_exchanger_type
    @heat_exchanger_type.setter
    def heat_exchanger_type(self, heat_exchanger_type):
        try:
            self._Cb_func = Cb_dict[heat_exchanger_type]
        except KeyError:
            raise AttributeError("heat exchange type must be one of the following: "
                                 f"{', '.join(Cb_dict)}")
        self._heat_exchanger_type = heat_exchanger_type     

    def reset_cache(self):
        for i in self.outs: i.reset_cache()

    def _assert_compatible_property_package(self):
        assert all([i.chemicals is j.chemicals for i, j in zip(self._ins, self._outs) if (i and j)]), (
            "inlet and outlet stream chemicals are incompatible; "
            "try using the `thermo` keyword argument to initialize the unit operation "
            "with a compatible thermodynamic property package"
        )

    def _design(self):
        # Get duty (kW)
        Q = abs(self.Q) / 3600
        
        if Q <= 1e-12: 
            self.design_results.clear()
            return
        
        ###  Use LMTD correction factor method  ###
        Design = self.design_results
        
        # Get cold and hot inlet and outlet streams
        ci, hi, co, ho = ht.order_streams(*self.get_streams())
       
        # Get log mean temperature difference
        Tci = ci.T
        Thi = hi.T
        Tco = co.T
        Tho = ho.T
        LMTD = ht.compute_LMTD(Thi, Tho, Tci, Tco)
        
        # Get correction factor
        ft = self.ft
        if not ft:
            N_shells = self.N_shells
            ft = ht.compute_Fahkeri_LMTD_correction_factor(Tci, Thi, Tco, Tho, N_shells)
        
        # Get overall heat transfer coefficient
        U = self.U or ht.heuristic_overall_heat_transfer_coefficient(ci, hi, co, ho)
        dP_tube, dP_shell = ht.heuristic_tubeside_and_shellside_pressure_drops(ci, hi, co, ho)
        
        # TODO: Complete design of heat exchanger to find L
        # For now assume lenght is 20 ft
        L = 20
        
        # Design pressure
        P = max((ci.P, hi.P))
        Design['Area'] = 10.763 * ht.compute_heat_transfer_area(abs(LMTD), U, Q, ft)
        Design['Overall heat transfer coefficient'] = U
        Design['Log-mean temperature difference'] = LMTD
        Design['Fouling correction factor'] = ft
        Design['Tube side pressure drop'] = dP_tube
        Design['Shell side pressure drop'] = dP_shell
        Design['Operating pressure'] = P * 14.7/101325 # psi
        Design['Total tube length'] = L

    def _cost(self):
        Design = self.design_results
        if not Design: return
            
        A = Design['Area']
        L = Design['Total tube length']
        P = Design['Operating pressure']
        
        if A < 150: # Double pipe
            P = P/600
            F_p = 0.8510 + 0.1292*P + 0.0198*P**2
            # Assume outer pipe carbon steel, inner pipe stainless steel
            F_m = 2 
            A_min = 2.1718
            if A < A_min:
                F_l = A/A_min
                A = A_min
            else:    
                F_l = 1
            heat_exchanger_type = 'Double pipe'
            C_b = compute_double_pipe_purchase_price(A, bst.CE)
        else: # Shell and tube
            F_m = compute_shell_and_tube_material_factor(A,  *self._F_Mab)
            F_l = 1 if L > 20 else np.polyval(p2, L)
            P = P/100
            F_p = 0.9803 + 0.018*P + 0.0017*P**2
            heat_exchanger_type = self.heat_exchanger_type
            C_b = self._Cb_func(A, bst.CE)
        
        # Free on board purchase prize
        self.F_M[heat_exchanger_type] = F_m
        self.F_P[heat_exchanger_type] = F_p
        self.F_D[heat_exchanger_type] = F_l
        self.baseline_purchase_costs[heat_exchanger_type] = C_b


class HXutility(HX):
    """
    Create a heat exchanger that changes temperature of the
    outlet stream using a heat utility.

    Parameters
    ----------
    ins : stream
        Inlet.
    outs : stream
        Outlet.
    T=None : float
        Temperature of outlet stream [K].
    V=None : float
        Vapor fraction of outlet stream.
    rigorous=False : bool
        If true, calculate vapor liquid equilibrium
    U=None : float, optional
        Enforced overall heat transfer coefficent [kW/m^2/K].
    heat_exchanger_type : str, optional
        Heat exchanger type. Defaults to "Floating head".
    N_shells : int, optional
        Number of shells. Defaults to 2.
    ft : float, optional
        User imposed correction factor.
    heat_only : bool, optional
        If True, heat exchanger can only heat.
    cool_only : bool, optional
        If True, heat exchanger can only cool.
        
    Notes
    -----
    Must specify either `T` or `V` when creating a HXutility object.
    
    Examples
    --------
    Run heat exchanger by temperature:
    
    >>> from biosteam.units import HXutility
    >>> from biosteam import Stream, settings
    >>> settings.set_thermo(['Water', 'Ethanol'], cache=True)
    >>> feed = Stream('feed', Water=200, Ethanol=200)
    >>> hx = HXutility('hx', ins=feed, outs='product', T=50+273.15,
    ...                rigorous=False) # Ignore VLE
    >>> hx.simulate()
    >>> hx.show()
    HXutility: hx
    ins...
    [0] feed
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Water    200
                        Ethanol  200
    outs...
    [0] product
        phase: 'l', T: 323.15 K, P: 101325 Pa
        flow (kmol/hr): Water    200
                        Ethanol  200
    >>> hx.results()
    Heat exchanger                                            Units       hx
    Low pressure steam  Duty                                  kJ/hr 1.01e+06
                        Flow                                kmol/hr     26.1
                        Cost                                 USD/hr     6.21
    Design              Area                                   ft^2       57
                        Overall heat transfer coefficient  kW/m^2/K      0.5
                        Log-mean temperature difference           K      101
                        Fouling correction factor                          1
                        Tube side pressure drop                 psi      1.5
                        Shell side pressure drop                psi        5
                        Operating pressure                      psi       50
                        Total tube length                        ft       20
    Purchase cost       Double pipe                             USD 4.74e+03
    Total purchase cost                                         USD 4.74e+03
    Utility cost                                             USD/hr     6.21
    
    Run heat exchanger by vapor fraction:
    
    >>> feed = Stream('feed', Water=200, Ethanol=200)
    >>> hx = HXutility('hx', ins=feed, outs='product', V=1,
    ...                rigorous=True) # Include VLE
    >>> hx.simulate()
    >>> hx.show()
    HXutility: hx
    ins...
    [0] feed
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Water    200
                        Ethanol  200
    outs...
    [0] product
        phase: 'g', T: 357.45 K, P: 101325 Pa
        flow (kmol/hr): Water    200
                        Ethanol  200
    >>> hx.results()
    Heat exchanger                                            Units       hx
    Low pressure steam  Duty                                  kJ/hr 1.94e+07
                        Flow                                kmol/hr      499
                        Cost                                 USD/hr      119
    Design              Area                                   ft^2      680
                        Overall heat transfer coefficient  kW/m^2/K        1
                        Log-mean temperature difference           K     80.8
                        Fouling correction factor                          1
                        Tube side pressure drop                 psi      1.5
                        Shell side pressure drop                psi      1.5
                        Operating pressure                      psi       50
                        Total tube length                        ft       20
    Purchase cost       Floating head                           USD 2.61e+04
    Total purchase cost                                         USD 2.61e+04
    Utility cost                                             USD/hr      119

    """
    line = 'Heat exchanger'
    _graphics = utility_heat_exchanger_graphics
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 T=None, V=None, rigorous=False, U=None, H=None,
                 heat_exchanger_type="Floating head",
                 material="Carbon steel/carbon steel",
                 N_shells=2,
                 ft=None,
                 heat_only=None,
                 cool_only=None):
        super().__init__(ID, ins, outs, thermo)
        self.T = T #: [float] Temperature of outlet stream (K).
        self.V = V #: [float] Vapor fraction of outlet stream.
        self.H = H #: [float] Enthalpy of outlet stream.
        
        #: [bool] If true, calculate vapor liquid equilibrium
        self.rigorous = rigorous
        
        #: [float] Enforced overall heat transfer coefficent (kW/m^2/K)
        self.U = U
        
        #: [int] Number of shells for LMTD correction factor method.
        self.N_shells = N_shells
        
        #: [float] User imposed correction factor.
        self.ft = ft

        #: [bool] If True, heat exchanger can only heat. 
        self.heat_only = heat_only
        
        #: [bool] If True, heat exchanger can only cool. 
        self.cool_only = cool_only
        
        self.material = material
        self.heat_exchanger_type = heat_exchanger_type
    
    def _init_utils(self):
        # tuple[HeatUtility] All heat utilities associated to unit
        self._heat_utilities = tuple([bst.HeatUtility(heat_exchanger=self)
                                      for i in range(self._N_heat_utilities)])
        
        # [PowerUtility] Electric utility associated to unit
        self.power_utility = bst.PowerUtility()
    
    @property
    def Q(self):
        """[float] Total heat transfered."""
        return abs(self.heat_utilities[0].unit_duty)
    
    @property
    def heat_utilities(self):
        return self._heat_utilities
    @heat_utilities.setter
    def heat_utilities(self, heat_utilities):
        self._heat_utilities = heat_utilities = tuple(heat_utilities)
        for i in heat_utilities: i.heat_exchanger = self
    
    def simulate_as_auxiliary_exchanger(self, duty, stream):
        self.outs[0] = stream.proxy()
        self.ins[0] = stream.proxy()
        hu = self.heat_utilities[0]
        hu.heat_exchanger = None
        hu(duty, stream.T)
        super()._design()
        self._cost()
        self._load_capital_costs()
        
    def _run(self):
        feed = self.ins[0]
        outlet = self.outs[0]
        outlet.copy_flow(feed)
        outlet.P = feed.P
        T = self.T
        V = self.V
        H = self.H
        T_given = T is not None
        V_given = V is not None
        H_given = H is not None
        N_given = T_given + V_given + H_given
        if N_given == 0:
            raise RuntimeError("no specification available; must define at either "
                               "temperature 'T', vapor fraction, 'V', or enthalpy 'H'")
        if self.rigorous:
            if N_given > 1:
                raise RuntimeError("may only specify either temperature, 'T', "
                                   "vapor fraction 'V', or enthalpy 'H', "
                                   "in a rigorous simulation")
            if V_given:
                if V == 0:
                    outlet.phase = 'l'
                    outlet.T = outlet.bubble_point_at_P().T
                elif V == 1:
                    outlet.phase = 'g'
                    outlet.T = outlet.dew_point_at_P().T
                elif 0 < V < 1:
                    outlet.vle(V=V, P=outlet.P)
                else:
                    raise RuntimeError("vapor fraction, 'V', must be a "
                                       "positive fraction")
            elif T_given:
                if outlet.isempty():
                    outlet.T = T
                else:
                    try:
                        outlet.vle(T=T, P=outlet.P)
                    except RuntimeError as e:
                        if len(outlet.phases) > 1: raise e
                        T_bubble = outlet.bubble_point_at_P().T
                        if T <= T_bubble:
                            outlet.phase = 'l'
                        else:
                            T_dew = outlet.dew_point_at_P().T
                            if T_dew >= T:
                                outlet.phase = 'g'
                            else:
                                raise RuntimeError('outlet in vapor-liquid equilibrium, but stream is linked')
                        outlet.T = T
            else:
                outlet.vle(H=H, P=outlet.P)
        else:
            if T_given and H_given:
                raise RuntimeError("cannot specify both temperature, 'T' "
                                   "and enthalpy 'H'")
            if T_given:
                outlet.T = T
            else:
                outlet.T = feed.T
            if V_given:
                if V == 0:
                    outlet.phase = 'l'
                elif V == 1:
                    outlet.phase = 'g'
                else:
                    raise RuntimeError("vapor fraction, 'V', must be either "
                                       "0 or 1 in a non-rigorous simulation")
                if V == 1 and feed.vapor_fraction < 1. and (outlet.T + 1e-6) < feed.T:
                    raise ValueError('outlet cannot be cooler than inlet if boiling')
                if V == 0 and feed.vapor_fraction > 0. and outlet.T > feed.T + 1e-6:
                    raise ValueError('outlet cannot be hotter than inlet if condensing')
            else:
                phase = feed.phase
                if len(phase) == 1: outlet.phase = phase
            if H_given:
                outlet.H = H
        if self.heat_only:
            if outlet.H - feed.H < 0.: outlet.copy_like(feed)
        if self.cool_only:
            if outlet.H - feed.H > 0.: outlet.copy_like(feed)

    def get_streams(self):
        """
        Return inlet and outlet streams.
        
        Returns
        -------
        in_a : Stream
            Inlet a.
        in_b : Stream
            Inlet b.
        out_a : Stream
            Outlet a.
        out_b : Stream
            Outlet b.
        
        """
        in_a = self.ins[0]
        out_a = self.outs[0]
        hu = self.heat_utilities[0]
        in_b = hu.inlet_utility_stream
        out_b = hu.outlet_utility_stream
        return in_a, in_b, out_a, out_b

    def _design(self, duty=None):
        # Set duty and run heat utility
        if duty is None:
            duty = self.H_out - self.H_in
        inlet = self.ins[0]
        outlet = self.outs[0] 
        T_in = inlet.T
        T_out = outlet.T
        iscooling = duty < 0.
        if any(['gl' in i.phase for i in (inlet, outlet)]):
            if iscooling:
                if T_out > T_in: T_in = T_out
            else:
                if T_out < T_in: T_out = T_in
        self.heat_utilities[0](duty, T_in, T_out)
        super()._design()


class HXprocess(HX):
    """
    Counter current heat exchanger for process fluids. Rigorously transfers
    heat until the pinch temperature or a user set temperature limit is reached.
    
    Parameters
    ----------
    ins : stream sequence
        * [0] Inlet process fluid a
        * [1] Inlet process fluid b  
    outs : stream sequence
        * [0] Outlet process fluid a        
        * [1] Outlet process fluid b
    U=None : float, optional
        Enforced overall heat transfer coefficent [kW/m^2/K].
    dT=5. : float
        Pinch temperature difference (i.e. dT = abs(outs[0].T - outs[1].T)).
    T_lim0 : float, optional
        Temperature limit of outlet stream at index 0.
    T_lim1 : float, optional
        Temperature limit of outlet stream at index 1.
    heat_exchanger_type : str, optional
        Heat exchanger type. Defaults to 'Floating head'.
    N_shells=2 : int, optional
        Number of shells.
    ft=None : float, optional
        User enforced correction factor.
    phase0=None : 'l' or 'g', optional
        User enforced phase of outlet stream at index 0.
    phase1=None : 'l' or 'g', optional
        User enforced phase of outlet stream at index 1.
    H_lim0 : float, optional
        Enthalpy limit of outlet stream at index 0.
    H_lim1 : float, optional
        Enthalpy limit of outlet stream at index 1.
    
    Examples
    --------
    Rigorous heat exchange until pinch temperature is reached:
        
    >>> from biosteam.units import HXprocess
    >>> from biosteam import Stream, settings
    >>> settings.set_thermo(['Water', 'Ethanol'])
    >>> in_a = Stream('in_a', Ethanol=50, T=351.43, phase='g')
    >>> in_b = Stream('in_b', Water=200)
    >>> hx = HXprocess('hx', ins=(in_a, in_b), outs=('out_a', 'out_b'))
    >>> hx.simulate()
    >>> hx.show()
    HXprocess: hx
    ins...
    [0] in_a
        phase: 'g', T: 351.43 K, P: 101325 Pa
        flow (kmol/hr): Ethanol  50
    [1] in_b
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Water  200
    outs...
    [0] out_a
        phases: ('g', 'l'), T: 351.39 K, P: 101325 Pa
        flow (kmol/hr): (g) Ethanol  31.3
                        (l) Ethanol  18.7
    [1] out_b
        phase: 'l', T: 346.43 K, P: 101325 Pa
        flow (kmol/hr): Water  200
    
    >>> hx.results()
    Heat exchanger                                            Units       hx
    Design              Area                                   ft^2      213
                        Overall heat transfer coefficient  kW/m^2/K      0.5
                        Log-mean temperature difference           K     20.4
                        Fouling correction factor                          1
                        Tube side pressure drop                 psi      1.5
                        Shell side pressure drop                psi        5
                        Operating pressure                      psi     14.7
                        Total tube length                        ft       20
    Purchase cost       Floating head                           USD 2.06e+04
    Total purchase cost                                         USD 2.06e+04
    Utility cost                                             USD/hr        0
    
    Sensible fluids case with user enfored outlet phases 
    (more computationally efficient):
        
    >>> from biosteam.units import HXprocess
    >>> from biosteam import Stream, settings
    >>> settings.set_thermo(['Water', 'Ethanol'])
    >>> in_a = Stream('in_a', Water=200, T=350)
    >>> in_b = Stream('in_b', Ethanol=200)
    >>> hx = HXprocess('hx', ins=(in_a, in_b), outs=('out_a', 'out_b'),
    ...                phase0='l', phase1='l')
    >>> hx.simulate()
    >>> hx.show()
    HXprocess: hx
    ins...
    [0] in_a
        phase: 'l', T: 350 K, P: 101325 Pa
        flow (kmol/hr): Water  200
    [1] in_b
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Ethanol  200
    outs...
    [0] out_a
        phase: 'l', T: 303.15 K, P: 101325 Pa
        flow (kmol/hr): Water  200
    [1] out_b
        phase: 'l', T: 327.92 K, P: 101325 Pa
        flow (kmol/hr): Ethanol  200
    >>> hx.results()
    Heat exchanger                                            Units       hx
    Design              Area                                   ft^2      367
                        Overall heat transfer coefficient  kW/m^2/K      0.5
                        Log-mean temperature difference           K     11.5
                        Fouling correction factor                          1
                        Tube side pressure drop                 psi        5
                        Shell side pressure drop                psi        5
                        Operating pressure                      psi     14.7
                        Total tube length                        ft       20
    Purchase cost       Floating head                           USD 2.23e+04
    Total purchase cost                                         USD 2.23e+04
    Utility cost                                             USD/hr        0
    
    """
    line = 'Heat exchanger'
    _graphics = process_heat_exchanger_graphics
    _N_heat_utilities = 0
    _N_ins = 2
    _N_outs = 2
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 U=None, dT=5., T_lim0=None, T_lim1=None,
                 material="Carbon steel/carbon steel",
                 heat_exchanger_type="Floating head",
                 N_shells=2, ft=None, 
                 phase0=None,
                 phase1=None,
                 H_lim0=None,
                 H_lim1=None):
        super().__init__(ID, ins, outs, thermo)
        
        #: [float] Enforced overall heat transfer coefficent (kW/m^2/K)
        self.U = U
        
        #: [float] Total heat transfered.
        self.Q = None
        
        #: Number of shells for LMTD correction factor method.
        self.N_shells = N_shells
        
        #: User imposed correction factor.
        self.ft = ft
        
        #: [float] Pinch temperature difference.
        self.dT = dT 
        
        #: [float] Temperature limit of outlet stream at index 0.
        self.T_lim0 = T_lim0
        
        #: [float] Temperature limit of outlet stream at index 1.
        self.T_lim1 = T_lim1
        
        #: [float] Temperature limit of outlet stream at index 0.
        self.H_lim0 = H_lim0
        
        #: [float] Temperature limit of outlet stream at index 1.
        self.H_lim1 = H_lim1
        
        #: Enforced phase of outlet at index 0
        self.phase0 = phase0
        
        #: Enforced phase of outlet at index 1
        self.phase1 = phase1
        
        self.material = material
        self.heat_exchanger_type = heat_exchanger_type
        
    def get_streams(self):
        s_in_a, s_in_b = self.ins
        s_out_a, s_out_b  = self.outs
        return s_in_a, s_in_b, s_out_a, s_out_b
    
    def _setup(self):
        super()._setup()
        for i in self._ins:
            if i.source: i.empty()
            
    def simulate(self):
        self._run()
        self._summary()
    
    def _run(self):
        s1_in, s2_in = self._ins
        s1_out, s2_out = self._outs
        if s1_in.isempty():
            s1_out.empty()
            s2_out.copy_like(s2_in)
            self.Q = 0.
        elif s2_in.isempty():
            s2_out.empty()
            s1_out.copy_like(s1_in)
            self.Q = 0.
        else:
            s1_out.copy_like(s1_in)
            s2_out.copy_like(s2_in)
            self._run_counter_current_heat_exchange()
            for s_out in (s1_out, s2_out):
                if isinstance(s_out, bst.MultiStream):
                    phase = s_out.phase
                    if len(phase) == 1: s_out.phase = phase
    
    def _run_counter_current_heat_exchange(self):
        self.Q = ht.counter_current_heat_exchange(*self._ins, *self._outs,
                                                  self.dT, self.T_lim0, self.T_lim1,
                                                  self.phase0, self.phase1, 
                                                  self.H_lim0, self.H_lim1)
