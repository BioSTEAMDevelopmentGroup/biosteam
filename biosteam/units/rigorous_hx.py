# -*- coding: utf-8 -*-
"""
Created on Sun Nov 18 18:27:37 2018

@author: yoelr
"""
from biosteam import Unit, Stream, MixedStream
from biosteam.exceptions import biosteamError
from biosteam.utils import SmartBook
from fluids import nearest_pipe
import ht
import numpy as np
pi = np.pi

class HX(Unit):
    """Create a heat exchanger that changes temperature of the output stream.

    **Parameters**

         T: [float] Temperature of output stream (K)

    **ins**
    
        [0] Input stream
        
    **outs**
    
        [1] Output stream

    """
    _N_ins = 1
    _N_outs = 1
    _N_heat_util = 1
    kwargs = {'T': None}

    def _setup(self):
        self.outs[0].T = self.kwargs['T']
        self._cached = {'s_tube': Stream('s_tube'),
                        's_shell': Stream('s_shell')}

    def _run(self):        feed = self.ins[0]
        s = self.outs[0]
        s.mol = feed.mol
        s.P = feed.P
        s.phase = feed.phase

    def _operation(self):
        util = self.heat_utilities[0]
        util(self._H_out-self._H_in, self.ins[0].T, self.outs[0].T)
        return self.results['Operation']

    def _design(self):
        """
        * 'Total Area': (m^2)
        * 'Clean total Area' : (m^2)
        * 'Overall heat transfer coefficent - with fouling' : (W/mK)
        * 'Overall heat transfer coefficent - no fouling'   : (W/mK)
        * 'Thermal Load' : (kW)
        * 'Pressure drop tube side': (Pa/m)
        * 'Pressure drop shell side' : (Pa/m)
        * 'Tube diameter in': (m)
        * 'HX Lenght' : (m)
        * 'Shell diameter': (m)
        * 'Reynolds tube' : ()
        * 'Reynolds shell' : ()
        * 'Nusselt tube' : ()
        * 'Nusselt shell': ()
        * 'Hydraulic diameter': (m)
        * 'ft': ()
        * 'Thickness': (m)
        """
        s_in = self.ins[0]
        s_out = self.outs[0]
        cached = self._cached
        util = self.heat_utilities[0]
        cached['u_in'] = u_in = util._fresh
        vap_out = util._vap
        liq_out = util._liq
        cached['u_out'] = u_out = MixedStream('u_out')
        u_out.T = liq_out.T
        u_out.P = liq_out.P
        u_out.liquid_mol = liq_out.mol
        u_out.vapor_mol = vap_out.mol
        Re_i = 5000
        Re_o = 6000
        try:
            return self._shellnTube(s_in, u_in, s_out, u_out, Re_i, Re_o)
        except:
            print(f'Problem in {self.ID} design')

    def _shellnTube(self, s_in, u_in, s_out, u_out, Re_i, Re_o) -> 'Design[dict]':
        """Return all design requirements for a shell and tube.
        
        **Parameters**
        
            s_in: [Stream] Stream s entering
            
            s_out: [Stream] Stream s leaving
                
            u_in: [Stream] Stream u entering
            
            u_out: [Stream] Stream u leaving
            
            Re_i: [float] Reynolds number of fluid in inner tube (m/s)
            
            Re_o: [float] Reynolds number of fluid in shell (m/s)
            
        The overall heat transfer coefficient is directy calculated. From Perry, chapter 11-3 'Thermal Design of Heat-Transfer Equipment', general value can be found for comparison.
        
        **Returns**
        
            * 'Total Area': (m^2)
            * 'Clean total Area' : (m^2)
            * 'Overall heat transfer coefficent - with fouling' : (W/mK)
            * 'Overall heat transfer coefficent - no fouling'   : (W/mK)
            * 'Thermal Load' : (kW)
            * 'Pressure drop tube side': (Pa/m)
            * 'Pressure drop shell side' : (Pa/m)
            * 'Tube diameter in': (m)
            * 'HX Lenght' : (m)
            * 'Shell diameter': (m)
            * 'Reynolds tube' : ()
            * 'Reynolds shell' : ()
            * 'Nusselt tube' : ()
            * 'Nusselt shell': ()
            * 'Hydraulic diameter': (m)
            * 'ft': ()
            * 'Thickness': (m)
        
        """
        # Find cold and hot entrances and exit streams
        sci, shi, sco, sho = (s_in, u_in, s_out, u_out) if s_in.T < u_in.T else (u_in, s_in, u_out, s_out)
        
        # TODO: include code for vapors too
        for s in (sci, shi, sco, sho):
            if s.phase == 'g':
                return
        
        # Get temperatures
        Tci, Tco, Thi, Tho = sci.T, sco.T, shi.T, sho.T
        
        # Choose a TEMA (Tubular Exchanger Manufacturer Association) heat exchanger
        # http://www.thermopedia.com/content/1121/
        # Set number of shells
        number_shells = 1 # TEMA E, most common type has 1 shell
        
        # Get the h value for calculating U 
        # Need: fluid properties, which are calculated at mean temperature between inlet and outlet both tube and                 shell side
        # Mean temperature
        Tc_ave = (Tci + Tco)/2
        Th_ave = (Thi + Tho)/2
         
        # Streams for calculating non-dimensional numbers
        s_shell = self._cached['s_shell']
        s_tube = self._cached['s_shell']
         
        # Choose which fluid goes in tube side to minimize heat losses (configuration with smallest temperature difference between shell and surroundings).
        if abs(298.15 - sci.T) - abs(shi.T - 298.15)> 0:
            # Cold stream goes in tube side
            s_tube_heating= True
            s_tube.copy_like(sci); s_tube.T = Tc_ave
            s_shell.copy_like(shi); s_shell.T = Th_ave
        else:
            # Hot stream goes in tube side
            s_tube_heating= False 
            s_tube.copy_like(shi); s_tube.T = Th_ave
            s_shell.copy_like(sci); s_shell.T = Tc_ave
            
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
        
        
        # Calculate velocity according to nominal pipe size (NPS)
        A_in = pi/4*Tid**2
        v_i = mass_i/(A_in*Tid)
        
        # Recalculate Re and Pr for tube
        Re_i = (rho_i * v_i * Tid) / mu_i
        Pr_i = s_tube.Pr
        
        # For the outer tube (shell)
        mass_o = s_shell.massnet/3600 # kg/s
        rho_o = s_shell.rho
        mu_o = s_shell.mu
        Sid = (4/pi * mass_o*mu_o/Re_i)**(0.5)/rho_o
        v_o = Re_o*rho_o/(mu_o*Sid)
        Pr_o = s_shell.Pr
        
        # Hydraulic diameter for shell side 
        D_eq = Tod-Tid 
        
        # Get nusselt number based on correlation
        
        if Re_i <= 2300:
            Nu_i = ht.conv_internal.laminar_T_const()
        elif Re_i > 2300: 
            # For turbulent flow, the Nusselt correlation change if the fluid is heated or cooled. When using this formula check if the fluid inside the inner tube is heated or cooled
            Nu_i = ht.conv_internal.turbulent_Dittus_Boelter(Re=Re_i, Pr=Pr_i, heating=s_tube_heating, revised=True)
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
        U = (1/U_clean + Rif + Rof) **(-1) 
        
        # Calculate UA by LMTD correction factor method
        LMTD = ht.LMTD(Thi, Tho, Tci, Tco) # log mean temperature difference
        try:
            ft = ht.F_LMTD_Fakheri(Thi, Tho, Tci, Tco, shells= number_shells) # correction factor
        except ValueError:
            ft = 1/1.1# !!!! add 10% or is counted in ft?
        
        Q = abs(s_out.H - s_in.H) 
        Q_kW = Q/1000
        A_tot = Q/(U*LMTD*ft)
        
        # Area without acounting for Fouling 
        A_clean = Q/(U_clean*LMTD*ft)
        
        # Get HX lenght (tube lenght) 
        L = A_tot / (pi*Tod)
         
        # Calculate pressure drops
        # Friction factor tube side
        if Re_i <= 2300 :
            f_t = 64 / Re_i
        elif 2300 < Re_i <= 4000:
            f_t = 0.158 * Re_i **(-0.25)
        elif Re_i > 4000:
            f_t =  0.0792 * Re_i **(-0.25)
            
        #Friction factor shell side
        f_m = 0.4 * (Re_o)**(-0.2)
            
        #density of the fluids    
        rho_in=s_tube.rho
        rho_out=s_shell.rho
        
        #Pressure drops inner tube per unit of lenght
        Dp_t = 4 * f_t * (1/Tid) * 0.5 * rho_in *(v_i**2) 
        
        Dp_s = 4 * f_m * (Sid/D_eq) * 0.5 * rho_out * ((v_o)**2)
        #Dp_s = ht.conv_tube_bank.dP_Zukauskas(Re=Re_o, n= 1, ST=0, SL=0, D=Tod, rho=rho_out, Vmax=v_out)
        
        Design = self.results['Design']
        Design['Total Area'] = A_tot
        Design['Clean total Area'] = A_clean
        Design['Overall heat transfer coefficent'] = U
        Design['Clean overall heat transfer coefficent'] = U_clean
        Design['Thermal Load'] = Q_kW
        Design['Pressure drop tube side'] = Dp_t
        Design['Pressure drop shell side'] = Dp_s
        Design['Tube diameter in'] = Tid
        Design['HX Lenght'] = L
        Design['Shell diameter'] = Sid
        Design['Reynolds tube'] = Re_i
        Design['Reynolds shell'] = Re_o
        Design['Nusselt tube'] = Nu_i,
        Design['Nusselt shell'] = Nu_o
        Design['Hydraulic diameter'] = D_eq
        Design['ft'] = ft
        Design['Thickness'] = tx
        return Design
    
    def _cost(self):
        self.results['Cost']