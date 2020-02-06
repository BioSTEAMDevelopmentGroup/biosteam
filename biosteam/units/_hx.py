# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 14:38:34 2018

@author: yoelr
"""
from .. import Unit
from thermosteam import Stream
from fluids import nearest_pipe
import ht
import numpy as np
import biosteam as bst
ln = np.log
exp = np.exp
pi = np.pi

# Lenght factor 
x = np.array((8, 13, 16, 20)) 
y = np.array((1.25, 1.12,1.05,1))
p2 = np.polyfit(x,y,2)

# Materials of Construction Shell/Tube	       a and b in Eq. (16.44)
F_Mdict =  {'Carbon steel/carbon steel':       (0, 0),
            'Carbon steel/brass':	            (1.08, 0.05),
            'Carbon steel/stainles steel':	  (1.75, 0.13),
            'Carbon steel/Monel':	            (2.1, 0.13),
            'Carbon steel/titanium':	       (5.2, 0.16),
            'Carbon steel/Cr-Mo steel':        (1.55, 0.05),
            'Cr-Mo steel/Cr-Mo steel':	       (1.7, 0.07),
            'Stainless steel/stainless steel': (2.7, 0.07),
            'Monel/Monel':	                 (3.3, 0.08),
            'Titanium/titanium':	            (9.6, 0.06)}

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
        
        
class HX(Unit, isabstract=True):
    """Abstract class for counter current heat exchanger.

    **Abstract methods**
    
    _get_streams()
        Should return two input streams and two output streams that exchange heat.

    **Abstract attributes**
    
    _duty : float
        The heat transfer requirement (kJ/hr)
    U : float
        Enforced overall heat transfer coefficent (kW/m^2/K)

    """
    line = 'Heat Exchanger'
    _units = {'Area': 'ft^2',
              'Overall heat transfer coefficient': 'kW/m^2/K',
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
        if self._Type == "Double pipe":
            return self.BM_douple_pipe
        else:
            return self.BM_shell_and_tube
    
    # Heat exchanger type
    _Type = 'Floating head'
    
    # Number of shells
    _N_shells = 2
    
    # Material factor function
    _F_Mstr = 'Carbon steel/carbon steel'
    _Cb_func = staticmethod(Cb_dict['Floating head'])
    _F_Mab = (0, 0)
    
    # Correction factor
    _ft = None

    @property
    def N_shells(self):
        """Number of shells for LMTD correction factor method."""
        return self._N_shells
    @N_shells.setter
    def N_shells(self, N):
        self._N_shells = N

    @property
    def ft(self):
        """User imposed correction factor."""
        return self._ft
    @ft.setter
    def ft(self, ft):
        self._ft = ft

    @property
    def material(self):
        """Default 'Carbon steel/carbon steel'"""
        return self._F_Mstr
    @material.setter
    def material(self, material):
        try:
            self._F_Mab = F_Mdict[material]
        except KeyError:
            dummy = str(F_Mdict.keys())[11:-2]
            raise ValueError(f"material must be one of the following: {dummy}")
        self._F_Mstr = material  
    
    @property
    def Type(self):
        return self._Type
    @Type.setter
    def Type(self, Type):
        try:
            self._Cb_func = Cb_dict[Type]
        except KeyError:
            dummy = str(Cb_dict.keys())[11:-2]
            raise ValueError(f"heat exchange type must be one of the following: {dummy}")
        self._Type = Type
        
    @staticmethod
    def _U_table(ci, hi, co, ho):
        """Return an estimate of the overall heat transfer coefficient (kW/m^2/K)."""
        # TODO: Base U on Table 18.5, Warren D. Seider et. al. Product and Process Design Principles. (2016)
        cip, hip, cop, hop = ci.phase, hi.phase, co.phase, ho.phase
        phases = cip + hip + cop + hop
        if 'g' in phases:
            if ('g' in hip and 'l' in hop) and ('l' in cip and 'g' in cop):
                return 1.0
            else:
                return 0.5
        else:
            return 0.5

    @staticmethod
    def _get_dP(inlet_phase, outlet_phase):
        """Return pressure drop (psi) based on heuristics.
        
        Parameters
        ----------
        inlet_phase: str
        outlet_phase : str
            
        """
        if ('l' in inlet_phase and 'g' in outlet_phase) or ('g' in inlet_phase and 'l' in outlet_phase):
            # Latent fluid (boiling or condensing)
            return 1.5
        elif inlet_phase == 'l':
            # Sensible liquid
            return 5
        elif outlet_phase == 'g':
            # Sensible vapor
            return 3
        
    @classmethod
    def _Dp_table(cls, ci, hi, co, ho, inside_heating):
        """Return an estimate of pressure drops.
        
        Returns
        -------  
        dP_in : float
            Pressure drop inside (psi)
        dP_out : float
            Pressure drop outside (psi)
        
        """
        cip, hip, cop, hop = ci.phase, hi.phase, co.phase, ho.phase
        dP_c = cls._get_dP(cip, cop)
        dP_h = cls._get_dP(hip, hop)
        if inside_heating:
            dP_in = dP_c
            dP_out = dP_h
        else:
            dP_in = dP_h
            dP_out = dP_c
        return dP_in, dP_out

    @staticmethod
    def _inside_isheating(ci, hi, co, ho) -> bool:
        """Return True if cold stream goes in tube side (as opposed to shell side)."""
        return abs(298.15 - ci.T) - abs(hi.T - 298.15) > 0

    @staticmethod
    def _shellntube_streams(ci, hi, co, ho, inside_heating):
        """Return shell and tube streams for calculating non-dimensional numbers."""
        # Mean temperatures
        Tci, Thi, Tco, Tho = ci.T, hi.T, co.T, ho.T
        Tc_ave = (Tci + Tco)/2
        Th_ave = (Thi + Tho)/2
         
        # Choose which fluid goes in tube side to minimize heat losses (configuration with smallest temperature difference between shell and surroundings).
        s_shell = Stream('s_shell')
        s_tube = Stream('s_tube')
        
        if inside_heating:
            # Cold stream goes in tube side
            s_tube.copy_like(ci); s_tube.T = Tc_ave
            s_shell.copy_like(hi); s_shell.T = Th_ave
        else:
            # Hot stream goes in tube side
            s_tube.copy_like(hi); s_tube.T = Th_ave
            s_shell.copy_like(ci); s_shell.T = Tc_ave
        return s_tube, s_shell

    @staticmethod
    def _concentric_tubes(s_tube, s_shell, Re_i, Re_o, inside_heating):
        """Return overall heat transfer coefficient in kW/m2/K for concentric tubes."""
        raise NotImplementedError()
        # Get the h value for calculating U 
        # Need: fluid properties, which are calculated at mean temperature between inlet and outlet both tube and shell side
        
        # Use continuity equation to get Tid
        mass_i = s_tube.massnet/3600 # kg/s
        rho_i = s_tube.rho
        mu_i = s_tube.mu
        Tid = (4/pi * mass_i*mu_i/Re_i)**(0.5)/rho_i
         
        # Get Tube Outer diameter
        tx = 0.036576 # Assumption
        Tod = Tid + 2*tx            
        
        # TODO: Use this for the case with multiple tubes
        #NPS, Tid_new, Tod, tx = nearest_pipe(Di=Tid)    
        
        # # Calculate velocity according to nominal pipe size (NPS)
        # A_in = pi/4*Tid**2
        # v_i = mass_i/(A_in*Tid)
        
        # # Recalculate Re and Pr for tube
        # Re_i = (rho_i * v_i * Tid) / mu_i
        Pr_i = s_tube.Pr
        
        # #  For the outer tube (shell)
        # mass_o = s_shell.massnet/3600 # kg/s
        # rho_o = s_shell.rho
        # mu_o = s_shell.mu
        # Sid = (4/pi * mass_o*mu_o/Re_i)**(0.5)/rho_o
        # v_o = Re_o*rho_o/(mu_o*Sid)
        Pr_o = s_shell.Pr
        
        # Hydraulic diameter for shell side 
        D_eq = Tod-Tid 
        
        # Get nusselt number based on correlation
        
        if Re_i <= 2300:
            Nu_i = ht.conv_internal.laminar_T_const()
        elif Re_i > 2300: 
            # For turbulent flow, the Nusselt correlation change if the fluid is heated or cooled. When using this formula check if the fluid inside the inner tube is heated or cooled
            Nu_i = ht.conv_internal.turbulent_Dittus_Boelter(Re=Re_i, Pr=Pr_i, heating=inside_heating, revised=True)
        elif 10000 < Re_i < 100000 and 0.5 < Pr_i < 3:
            Nu_i = ht.conv_internal.turbulent_Colburn(Re=Re_i, Pr=Pr_i)
            
        # Nussel coefficient shell side
        Nu_o = ht.conv_external.Nu_cylinder_Zukauskas(Re=Re_o, Pr=Pr_o, Prw=None)
        
        # Conductivity
        k_in = s_tube.k
        k_out = s_shell.k
        
        # Calculate h for internal, out and in/out
        hi = Nu_i*k_in/Tid # Tube-side coefficient
        ho = Nu_o*k_out/D_eq # Shell-side coefficient
        hio = hi * (Tid/Tod)
         
        # Fouling resitance 
        # Available excel file with fouling factor for different fluids taken from Perry
        # TODO: Link to excel "FoulingFactor"
        Rif = 0.00009
        Rof = 0.000175
        
        #Calculate U 
        U_clean = (1/hio + 1/ho)**(-1)
        return (1/U_clean + Rif + Rof)**(-1) /1000
    
        
    @staticmethod
    def _order_streams(in1, in2, out1, out2):
        """Return cold and hot inlet and outlet streams.
        
        **Parameters**
        
            **in1:** [Stream] Inlet 1
            
            **in2:** [Stream] Inlet 2
                
            **out1:** [Stream] Outlet 1
            
            **out2:** [Stream] Outlet 2
        
        **Returns**
        
            **ci:** [Stream] Cold inlet
            
            **hi:** [Stream] Hot inlet
            
            **co:** [Stream] Cold outlet
            
            **ho:** [Stream] Hot outlet
        
        """
        if in1.T < in2.T:
            return in1, in2, out1, out2
        else:
            return in2, in1, out2, out1
        
    @staticmethod
    def _calc_ft(Tci, Thi, Tco, Tho, N_shells):
        """Return LMTD correction factor."""
        ft = 1
        if abs(Tco - Tci)/Tco > 0.01 and abs(Thi-Tho)/Tho > 0.01:
            try:
                ft = ht.F_LMTD_Fakheri(Thi, Tho, Tci, Tco,
                                       shells=N_shells)
            except ValueError: pass
            if ft > 1: ft = 0.6 # Accounts for worst case scenario
        return ft
    
    @staticmethod
    def _calc_area(LMTD, U, Q, ft):
        """Return Area by LMTD correction factor method.
        
        Parameters
        ----------
        LMTD : float
            Log mean temperature difference
        U : float
            Heat transfer coefficient
        Q : float
            Duty
        
        """
        return Q/(U*LMTD*ft)        

    def _design(self):
        ###  Use LMTD correction factor method  ###
        Design = self.design_results
        
        # Get cold and hot inlet and outlet streams
        ci, hi, co, ho = self._order_streams(*self._get_streams())
        inside_heating = self._inside_isheating(ci, hi, co, ho)
       
        # Get log mean temperature difference
        Tci, Thi, Tco, Tho = ci.T, hi.T, co.T, ho.T
        dTF1 = Thi-Tco
        dTF2 = Tho-Tci
        dummy = abs(dTF2/dTF1)
        LMTD = (dTF2 - dTF1)/ln(dummy) if dummy > 1.1 else dTF1
        LMTD = abs(LMTD)
        
        # Get correction factor
        ft = self._ft if self._ft else self._calc_ft(Tci, Thi, Tco, Tho, self._N_shells)
        
        # Get duty (kW)
        Q = abs(self._duty) / 3600
        
        # Get overall heat transfer coefficient
        U = self.U or self._U_table(ci, hi, co, ho)
        Dp_s, Dp_t = self._Dp_table(ci, hi, co, ho, inside_heating)
        
        # TODO: Complete design of heat exchanger to find L
        # For now assume 20 ft
        L = 20
        
        # Design pressure
        P = max((ci.P, hi.P))
        Design['Area'] = self._calc_area(LMTD, U, Q, ft) * 10.763
        Design['Overall heat transfer coefficient'] = U
        Design['Fouling correction factor'] = ft
        Design['Tube side pressure drop'] = Dp_t
        Design['Shell side pressure drop'] = Dp_s
        Design['Operating pressure'] = P*14.7/101325 # psi
        Design['Total tube length'] = L

    def _cost(self):
        Design = self.design_results
        A = Design['Area']
        L = Design['Total tube length']
        P = Design['Operating pressure']
        #Dp_t = Design['Dp_t']
        #Dp_s = Design['Dp_s']
        
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
            self.Type = 'Double pipe'
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
    """Create a heat exchanger that changes temperature of the output stream using a heat utility.

    Parameters
    ----------
    ins
        [0] Input stream
    outs
        [0] Output stream
    Define at least one:
        * T: [float] Temperature of output stream (K).
        * V: [float] Vapor fraction of output stream.
    rigorous=False : bool
        If true, calculate vapor liquid equilibrium
    U=None : float, optional
        Enforced overall heat transfer coefficent (kW/m^2/K)

    Examples
    --------
    :doc:`notebooks/HXutility Example`
    
    """
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 T=None, V=None, rigorous=False, U=None):
        super().__init__(ID, ins, outs, thermo)
        self.T = T #: Temperature of output stream (K).
        self.V = V #: Vapor fraction of output stream.
        
        #: [bool] If true, calculate vapor liquid equilibrium
        self.rigorous = rigorous
        
        #: [float] Enforced overall heat transfer coefficent (kW/m^2/K)
        self.U = U
        
    def _run(self):
        feed = self.ins[0]
        s = self.outs[0]
        s.copy_like(feed)
        T = self.T
        V = self.V
        V_given = V is not None
        if not (T or V_given):
            raise ValueError("must define at least one of the following: 'T', 'V'")
        if self.rigorous:
            if T and V_given:
                raise ValueError("may only define either temperature, 'T', or vapor fraction 'V', in a rigorous simulation")
            if V_given:
                s.vle(V=V, P=s.P)
            else:
                s.vle(T=T, P=s.P)
        else:
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

    def _get_streams(self):
        """Return cold and hot inlet and outlet streams.
        
        Returns
        -------
        ci : Stream
            Cold inlet.
        hi : Stream
            Hot inlet.
        co : Stream
            Cold outlet.
        ho : Stream
            Hot outlet.
        
        """
        in1 = self.ins[0]
        out1 = self.outs[0]
        hu = self.heat_utilities[0]
        in2 = hu.inlet_utility_stream
        out2 = hu.outlet_utility_stream
        return in1, in2, out1, out2

    def _design(self, duty=None):
        # Set duty and run heat utility
        if duty is None: duty = self.H_out - self.H_in
        self._duty = duty
        self.heat_utilities[0](duty, self.ins[0].T, self.outs[0].T)
        super()._design()


class HXprocess(HX):
    """Counter current heat exchanger for process fluids. Condenses/boils latent fluid until sensible fluid reaches pinch temperature.
    
    Parameters
    ----------
    ins  
        [0] Input process fluid 1
        
        [1] Input process fluid 2  
    outs
        [0] Output process fluid 1
        
        [1] Output process fluid 2
    U=None : float, optional
        Enforced overall heat transfer coefficent (kW/m^2/K)
    fluid_type : {None, 'ss', 'll', 'ls'}
        * **None:** Rigorously transfers heat until pinch temperature -not implemented yet-
        * **'ss':** Sensible-sensible fluids. Heat is exchanged until the pinch temperature is reached.
        * **'ll':** Latent-latent fluids. Heat is exchanged until one fluid completely changes phase.
        * **'ls':** Latent-sensible fluids. Heat is exchanged until either the pinch temperature is reached or the latent fluid completely changes phase.
    
    Examples
    --------
    :doc:`notebooks/HXprocess Example`
    
    """
    _N_heat_utilities = 0
    _N_ins = 2
    _N_outs = 2
    dT = 5 #: [float] Pinch temperature difference.
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 U=None, fluid_type='ss'):
        super().__init__(ID, ins, outs, thermo)
        self._Hvaps = self._not_zero_arr = self._cached_run_args = None
        
        #: [float] Enforced overall heat transfer coefficent (kW/m^2/K)
        self.U = U
        
        self.fluid_type = fluid_type
    
    @property
    def fluid_type(self):
        """
        [str] Must be one of the following:
            * **None:** Rigorously transfers heat until pinch temperature **-not implemented yet-**
            * **'ss':** Sensible-sensible fluids. Heat is exchanged until the pinch temperature is reached.
            * **'ll':** Latent-latent fluids. Heat is exchanged until one fluid completely changes phase.
            * **'ls':** Latent-sensible fluids. Heat is exchanged until either the pinch temperature is reached or the latent fluid completely changes phase.
        """
        return self._fluid_type
    @fluid_type.setter
    def fluid_type(self, fluid_type):
        if fluid_type not in ('ss', 'ls', 'll'):
            raise ValueError(f"fluid type must be either 'ss', 'ls', or 'll', not {fluid_type}")
        self._fluid_type = fluid_type
        
    def _get_streams(self):
        s_in1, s_in2 = self.ins
        s_out1, s_out2  = self.outs
        return s_in1, s_in2, s_out1, s_out2
    
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
            duty = s1.H - s1f.H
            s2f.T += duty/s2.C
        else:
            s2f.T = (s1.T - dT) if s1_hot else (s1.T + dT)
            duty = s2.H - s2f.H
            s1f.T += duty/s1.C
        self._duty = duty
    
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
        duty = s1_in.H - s1_out.H
        T_s2_new = s2_out.T + duty/s2_out.C
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
        
        self._duty = duty
    
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
            raise ValueError(f"no latent heat streams available for heat exchange with Type='ll'")
        
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
        duty = (sc_in.H - sc_out.H)
        sp_out.vle(P=sp_out.P, H=sp_in.H + duty)
        self._duty = abs(duty)
        
    
# elif U == 'Concentric tubes':
#     raise NotImplementedError("'Concentric tubes' not yet implemented")
#     Re_i = 30000 # Arbitrary, but turbulent ??
#     Re_o = 30000 # Arbitrary, but turbulent ??
#     s_tube, s_shell = self._shellntube_streams(ci, hi, co, ho, inside_heating)
#     U = self._concentric_tubes(s_tube, s_shell, Re_i, Re_o, inside_heating)
# else:
#     raise ValueError("overall heat transfer coefficient, 'U', should be one of the following: value (kW/m^2/K), 'Tabulated', or 'Concentric tubes'")