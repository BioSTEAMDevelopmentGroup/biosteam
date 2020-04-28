# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 14:38:34 2018

@author: yoelr
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

__all__ = ('HXutility', 'HXprocess')

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
    """Abstract class for counter current heat exchanger.

    **Abstract methods**
    
    get_streams()
        Should return two inlet streams and two outlet streams that exchange heat.

    """
    line = 'Heat Exchanger'
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
    BM_douple_pipe = 1.80
    BM_shell_and_tube = 3.17
    
    @property
    def BM(self):
        if self._shell_and_tube_type == "Double pipe":
            return self.BM_douple_pipe
        else:
            return self.BM_shell_and_tube

    @property
    def material(self):
        """Default 'Carbon steel/carbon steel'"""
        return self.material
    @material.setter
    def material(self, material):
        try:
            self._F_Mab = shell_and_tube_material_factor_coefficients[material]
        except KeyError:
            raise ValueError("material must be one of the following: "
                            f"{', '.join(shell_and_tube_material_factor_coefficients)}")
        self._material = material  
    
    @property
    def shell_and_tube_type(self):
        """[str] Shell and tube type. Purchase cost depends on this selection."""
        return self._shell_and_tube_type
    @shell_and_tube_type.setter
    def shell_and_tube_type(self, shell_and_tube_type):
        try:
            self._Cb_func = Cb_dict[shell_and_tube_type]
        except KeyError:
            raise ValueError("heat exchange type must be one of the following: "
                            f"{', '.join(Cb_dict)}")
        self._shell_and_tube_type = shell_and_tube_type     

    def _design(self):
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
        
        # Get duty (kW)
        Q = self.Q / 3600
        
        # Get overall heat transfer coefficient
        U = self.U or ht.heuristic_overall_heat_transfer_coefficient(ci, hi, co, ho)
        dP_tube, dP_shell = ht.heuristic_tubeside_and_shellside_pressure_drops(ci, hi, co, ho)
        
        # TODO: Complete design of heat exchanger to find L
        # For now assume lenght is 20 ft
        L = 20
        
        # Design pressure
        P = max((ci.P, hi.P))
        Design['Area'] = 10.763 * ht.compute_heat_transfer_area(abs(LMTD), U, abs(Q), ft)
        Design['Overall heat transfer coefficient'] = U
        Design['Log-mean temperature difference'] = LMTD
        Design['Fouling correction factor'] = ft
        Design['Tube side pressure drop'] = dP_tube
        Design['Shell side pressure drop'] = dP_shell
        Design['Operating pressure'] = P * 14.7/101325 # psi
        Design['Total tube length'] = L

    def _cost(self):
        Design = self.design_results
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
            self.shell_and_tube_type = 'Double pipe'
        else: # Shell and tube
            a, b = self._F_Mab
            F_m = a +  (A/100)**b
            F_l = 1 if L > 20 else np.polyval(p2, L)
            P = P/100
            F_p = 0.9803 + 0.018*P + 0.0017*P**2
        
        C_b_func = self._Cb_func
        C_b = C_b_func(A, bst.CE)
        
        # Free on board purchase prize 
        self.purchase_costs['Heat exchanger'] = F_p * F_l * F_m * C_b


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
    shell_and_tube_type : str, optional
        Heat exchanger type. Defaults to "Floating head".
    N_shells=2 : int
        Number of shells.
    ft=None : float
        User imposed correction factor.
    
    Notes
    -----
    Must specify either `T` or `V` when creating a HXutility object.
    
    Examples
    --------
    Run heat exchanger by temperature:
    
    >>> from biosteam.units import HXutility
    >>> from biosteam import Stream, settings
    >>> settings.set_thermo(['Water', 'Ethanol'])
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
    Heat Exchanger                                            Units       hx
    Low pressure steam  Duty                                  kJ/hr 1.01e+06
                        Flow                                kmol/hr     26.1
                        Cost                                 USD/hr     6.21
    Design              Area                                   ft^2     57.4
                        Overall heat transfer coefficient  kW/m^2/K      0.5
                        Log-mean temperature difference           K      100
                        Fouling correction factor                          1
                        Tube side pressure drop                 psi      1.5
                        Shell side pressure drop                psi        5
                        Operating pressure                      psi       50
                        Total tube length                        ft       20
    Purchase cost       Heat exchanger                          USD 4.75e+03
    Total purchase cost                                         USD 4.75e+03
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
        phases: ('g', 'l'), T: 357.45 K, P: 101325 Pa
        flow (kmol/hr): (g) Water    200
                            Ethanol  200
    >>> hx.results()
    Heat Exchanger                                            Units       hx
    Low pressure steam  Duty                                  kJ/hr 2.07e+07
                        Flow                                kmol/hr      532
                        Cost                                 USD/hr      126
    Design              Area                                   ft^2      733
                        Overall heat transfer coefficient  kW/m^2/K        1
                        Log-mean temperature difference           K     80.1
                        Fouling correction factor                          1
                        Tube side pressure drop                 psi      1.5
                        Shell side pressure drop                psi      1.5
                        Operating pressure                      psi       50
                        Total tube length                        ft       20
    Purchase cost       Heat exchanger                          USD 2.67e+04
    Total purchase cost                                         USD 2.67e+04
    Utility cost                                             USD/hr      126
    
    """
    line = 'Heat Exchanger'
    _graphics = utility_heat_exchanger_graphics
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 T=None, V=None, rigorous=False, U=None,
                 shell_and_tube_type="Floating head",
                 material="Carbon steel/carbon steel",
                 N_shells=2,
                 ft=None):
        super().__init__(ID, ins, outs, thermo)
        self.T = T #: Temperature of outlet stream (K).
        self.V = V #: Vapor fraction of outlet stream.
        
        #: [bool] If true, calculate vapor liquid equilibrium
        self.rigorous = rigorous
        
        #: [float] Enforced overall heat transfer coefficent (kW/m^2/K)
        self.U = U
        
        #: [float] Total heat transfered.
        self.Q = None
        
        #: Number of shells for LMTD correction factor method.
        self.N_shells = N_shells
        
        #: User imposed correction factor.
        self.ft = ft
        
        self.material = material
        self.shell_and_tube_type = shell_and_tube_type
    
    def simulate_as_auxiliary_exchanger(self, duty, stream):
        self.outs[0] = stream.proxy()
        self.ins[0] = stream.proxy()
        self.heat_utilities[0](duty, stream.T)
        self.Q = duty
        super()._design()
        self._cost()
        
    def _run(self):
        feed = self.ins[0]
        s = self.outs[0]
        s.copy_like(feed)
        T = self.T
        V = self.V
        V_given = V is not None
        if self.rigorous:
            if T and V_given:
                raise ValueError("may only define either temperature, 'T', or vapor fraction 'V', in a rigorous simulation")
            if V_given:
                s.vle(V=V, P=s.P)
            else:
                s.vle(T=T, P=s.P)
        elif (T or V_given):
            if V_given:
                if V == 0:
                    s.phase = 'l'
                elif V == 1:
                    s.phase = 'g'
                else:
                    s.phases = ('g', 'l')
                    mol = s.mol
                    s.imol['g'] = vap_mol = V * mol
                    s.imol['l'] = mol - vap_mol
            if T:
                s.T = T
        else:
            raise ValueError("must define at least one of the following: 'T', 'V'")

    def get_streams(self):
        """Return inlet and outlet streams.
        
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
        self.heat_utilities[0](duty, self.ins[0].T, self.outs[0].T)
        self.Q = duty
        super()._design()


class HXprocess(HX):
    """
    Counter current heat exchanger for process fluids. Condenses/boils latent fluid until sensible fluid reaches pinch temperature.
    
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
    fluid_type : None, 'ss', 'll', 'ls'
        * **None:** Rigorously transfers heat until pinch temperature (not implemented yet).
        * **'ss':** Sensible-sensible fluids. Heat is exchanged until the pinch temperature is reached.
        * **'ll':** Latent-latent fluids. Heat is exchanged until one fluid completely changes phase.
        * **'ls':** Latent-sensible fluids. Heat is exchanged until either the pinch temperature is reached or the latent fluid completely changes phase.
    shell_and_tube_type : str, optional
        Heat exchanger type. Defaults to 'Floating head'.
    N_shells=2 : int
        Number of shells.
    ft=None : float
        User imposed correction factor.
    
    Examples
    --------
    Sensible fluids case:
        
    >>> from biosteam.units import HXprocess
    >>> from biosteam import Stream, settings
    >>> settings.set_thermo(['Water', 'Ethanol'])
    >>> in_a = Stream('in_a', Water=200, T=350)
    >>> in_b = Stream('in_b', Ethanol=200)
    >>> hx = HXprocess('hx', ins=(in_a, in_b), outs=('out_a', 'out_b'),
    ...                fluid_type='ss')
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
        phase: 'l', T: 329.64 K, P: 101325 Pa
        flow (kmol/hr): Ethanol  200
    >>> hx.results()
    Heat Exchanger                                            Units       hx
    Design              Area                                   ft^2      207
                        Overall heat transfer coefficient  kW/m^2/K      0.5
                        Log-mean temperature difference           K     20.4
                        Fouling correction factor                          1
                        Tube side pressure drop                 psi        5
                        Shell side pressure drop                psi        5
                        Operating pressure                      psi     14.7
                        Total tube length                        ft       20
    Purchase cost       Heat exchanger                          USD 2.05e+04
    Total purchase cost                                         USD 2.05e+04
    Utility cost                                             USD/hr        0
    
    One latent fluid and one sensible fluid case:
    
    >>> in_a = Stream('in_a', Ethanol=50, T=351.43, phase='g')
    >>> in_b = Stream('in_b', Water=200)
    >>> hx = HXprocess('hx', ins=(in_a, in_b), outs=('out_a', 'out_b'),
    ...                fluid_type='ls')
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
        flow (kmol/hr): (g) Ethanol  32.94
                        (l) Ethanol  17.06
    [1] out_b
        phase: 'l', T: 346.43 K, P: 101325 Pa
        flow (kmol/hr): Water  200
    >>> hx.results()
    Heat Exchanger                                            Units       hx
    Design              Area                                   ft^2      625
                        Overall heat transfer coefficient  kW/m^2/K      0.5
                        Log-mean temperature difference           K     20.4
                        Fouling correction factor                          1
                        Tube side pressure drop                 psi      1.5
                        Shell side pressure drop                psi        5
                        Operating pressure                      psi     14.7
                        Total tube length                        ft       20
    Purchase cost       Heat exchanger                          USD 2.53e+04
    Total purchase cost                                         USD 2.53e+04
    Utility cost                                             USD/hr        0
    
    Latent fluids case:
    
    >>> in_a = Stream('in_a', Ethanol=50, T=351.43, phase='l')
    >>> in_b = Stream('in_b', Water=200, T=373.15, phase='g')
    >>> hx = HXprocess('hx', ins=(in_a, in_b), outs=('out_a', 'out_b'),
    ...                fluid_type='ll')
    >>> hx.simulate()
    >>> hx.show()
    HXprocess: hx
    ins...
    [0] in_a
        phase: 'l', T: 351.43 K, P: 101325 Pa
        flow (kmol/hr): Ethanol  50
    [1] in_b
        phase: 'g', T: 373.15 K, P: 101325 Pa
        flow (kmol/hr): Water  200
    outs...
    [0] out_a
        phase: 'g', T: 351.43 K, P: 101325 Pa
        flow (kmol/hr): Ethanol  50
    [1] out_b
        phases: ('g', 'l'), T: 373.12 K, P: 101325 Pa
        flow (kmol/hr): (g) Water  150.8
                        (l) Water  49.19
    >>> hx.results()
    Heat Exchanger                                            Units       hx
    Design              Area                                   ft^2      293
                        Overall heat transfer coefficient  kW/m^2/K        1
                        Log-mean temperature difference           K     21.7
                        Fouling correction factor                          1
                        Tube side pressure drop                 psi      1.5
                        Shell side pressure drop                psi      1.5
                        Operating pressure                      psi     14.7
                        Total tube length                        ft       20
    Purchase cost       Heat exchanger                          USD 2.14e+04
    Total purchase cost                                         USD 2.14e+04
    Utility cost                                             USD/hr        0
    
    """
    line = 'Heat Exchanger'
    _graphics = process_heat_exchanger_graphics
    _N_heat_utilities = 0
    _N_ins = 2
    _N_outs = 2
    dT = 5 #: [float] Pinch temperature difference.
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 U=None, fluid_type='ss',
                 material="Carbon steel/carbon steel",
                 shell_and_tube_type="Floating head",
                 N_shells=2,
                 ft=None):
        super().__init__(ID, ins, outs, thermo)
        
        #: [float] Enforced overall heat transfer coefficent (kW/m^2/K)
        self.U = U
        
        #: [float] Total heat transfered.
        self.Q = None
        
        #: Number of shells for LMTD correction factor method.
        self.N_shells = N_shells
        
        #: User imposed correction factor.
        self.ft = ft
        
        self.fluid_type = fluid_type
        self.material = material
        self.shell_and_tube_type = shell_and_tube_type
    
    @property
    def fluid_type(self):
        """
        [None, 'ss', 'll', or 'ls']
            * **None** Rigorously transfers heat until pinch temperature (not implemented yet).
            * **'ss'** Sensible-sensible fluids. Heat is exchanged until the pinch temperature is reached.
            * **'ll'** Latent-latent fluids. Heat is exchanged until one fluid completely changes phase.
            * **'ls'** Latent-sensible fluids. Heat is exchanged until either the pinch temperature is reached or the latent fluid completely changes phase.
        """
        return self._fluid_type
    @fluid_type.setter
    def fluid_type(self, fluid_type):
        if fluid_type not in ('ss', 'ls', 'll'):
            raise ValueError(f"fluid type must be either 'ss', 'ls', or 'll', not {repr(fluid_type)}")
        self._fluid_type = fluid_type
        
    def get_streams(self):
        s_in_a, s_in_b = self.ins
        s_out_a, s_out_b  = self.outs
        return s_in_a, s_in_b, s_out_a, s_out_b
    
    def _run(self):
        so0, so1 = self._outs
        si0, si1 = self._ins
        s0_F_mol = si0.F_mol
        s1_F_mol = si1.F_mol
        if not s0_F_mol:
            so1.copy_like(si1)
        elif not s1_F_mol:
            so0.copy_like(si0)
        else:
            fluid_type = self._fluid_type
            so0.copy_like(si0)
            so1.copy_like(si1)
            if fluid_type == 'ss':
                self._run_ss()
            elif fluid_type == 'ls':
                self._run_ls()
            elif fluid_type == 'll':
                self._run_ll()
    
    def _run_ss(self):
        dT = self.dT
        s1f, s2f = self.outs
        s1, s2 = self.ins
        s1_hot = s1.T > s2.T  # s2 energy will increase
        if s1.C < s2.C:
            s1f.T = (s2.T + dT) if s1_hot else (s2.T - dT)
            Q = s1.H - s1f.H
            s2f.T += Q/s2.C
        else:
            s2f.T = (s1.T - dT) if s1_hot else (s1.T + dT)
            Q = s2.H - s2f.H
            s1f.T += Q/s1.C
        self.Q = Q
    
    def _run_ls(self):
        s1_in, s2_in = self.ins
        s1_out, s2_out = self.outs
        dT = self.dT
        
        if s2_out.T > s1_in.T:
            # Stream s1 is boiling
            boiling = True
            s1_out.phase = 'g'
            dp = s1_out.dew_point_at_P()
            s1_out.T = dp.T
            T_pinch = s1_in.T + dT # Minimum
        else:
            # Stream s1 is condensing
            boiling = False
            s1_out.phase = 'l'
            bp = s1_out.bubble_point_at_P()
            s1_out.T = bp.T
            T_pinch = s1_in.T - dT # Maximum
        
        # Calculate maximum latent heat and new temperature of sensible stream
        Q = s1_in.H - s1_out.H
        T_s2_new = s2_out.T + Q/s2_out.C
        s1_out.copy_like(s1_in)
        
        if boiling and T_s2_new < T_pinch:
            # Partial boiling if temperature falls below pinch
            H0 = s2_in.H
            s2_out.T = T_pinch
            delH1 = H0 - s2_out.H
            s1_out.vle(P=s1_out.P, H=s1_in.H + delH1)
        elif not boiling and T_s2_new > T_pinch:
            # Partial condensation if temperature goes above pinch
            H0 = s2_in.H
            s2_out.T = T_pinch
            delH1 = H0 - s2_out.H
            s1_out.vle(P=s1_out.P, H=s1_in.H + delH1)
        elif boiling:
            s1_out.phase ='g'
            s2_out.T = T_s2_new
        elif not boiling:
            s1_out.phase = 'l'
            s2_out.T = T_s2_new
        self.Q = Q
    
    def _run_ll(self):
        s1_in, s2_in = self.ins
        s1_out, s2_out = self.outs
        
        if s1_in.T > s2_in.T and ('g' in s1_in.phase) and ('l' in s2_in.phase):
            for s in self.outs: s.phases = ('g', 'l')
            # s2 boiling, s1 condensing
            boiling = s2_out
            delH1 = s1_out['g'].Hvap
            delH2 = s2_out['l'].Hvap
        elif s1_in.T < s2_in.T and ('l' in s1_in.phase) and ('g' in s2_in.phase):
            for s in self.outs: s.phases = ('g', 'l')
            # s1 boiling, s2 condensing
            boiling = s1_out
            delH1 = s1_out['l'].Hvap
            delH2 = s2_out['g'].Hvap
        else:
            raise ValueError(f"no latent heat streams available for heat exchange with shell_and_tube_type='ll'")
        
        # sc: Stream that completely changes phase
        # sp: Stream that partialy changes phase
        if delH1 > delH2:
            sc_in = s2_in
            sc_out = s2_out
            sp_in = s1_in
            sp_out = s1_out
        else:
            sc_in = s1_in
            sc_out = s1_out
            sp_in = s2_in
            sp_out = s2_out
        
        if sc_out is boiling:
            sc_out.phase = 'g'
            dp = sc_out.dew_point_at_P()
            sc_out.T = dp.T
        else:
            sc_out.phase = 'l'
            bp = sc_out.bubble_point_at_P()
            sc_out.T = bp.T
        
        # VLE
        Q = (sc_in.H - sc_out.H)
        sp_out.vle(P=sp_out.P, H=sp_in.H + Q)
        self.Q = abs(Q)


