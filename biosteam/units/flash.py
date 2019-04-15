# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 16:21:56 2018

@author: yoelr
"""
from .. import Unit, MixedStream, PowerUtility
from math import pi, ceil
import copy
import numpy as np
from scipy.optimize import brentq, newton
from thermo import activity
from .designtools import vacuum_system, HNATable, FinalValue, \
                          VesselWeightAndWallThickness, Kvalue
from .hx import HX, HXutility

exp = np.exp
ln = np.log

# USDA Biodiesel: 
#V = np.pi*D**2/4*Ht*0.3048**3
#if not (0.1 < V < 70):
#    raise DesignError(f"Volume is out of bounds for costing")
#lambda V, CE: CE*13*V**0.62 # V (m3)
__all__ = ('Flash', 'SplitFlash', 'PartitionFlash',
           'RatioFlash', 'Evaporator_PQin', 'Evaporator_PV')


# %% Flash solutions

flash_error = activity.Rachford_Rice_flash_error

def analytical_flash_solution(zs, Ks):
    """Solution for 2 or 3 component flash vessel."""
    # Not currently in use
    len_ = len(zs)
    if len_ == 2:
        z1, z2 = zs
        K1, K2 = Ks
        V_over_F = (-K1*z1 - K2*z2 + z1 + z2)/(K1*K2*z1 + K1*K2 *
                                               z2 - K1*z1 - K1*z2 - K2*z1 - K2*z2 + z1 + z2)
    elif len_ == 3:
        z1, z2, z3 = zs
        K1, K2, K3 = Ks
        V_over_F = (-K1*K2*z1/2 - K1*K2*z2/2 - K1*K3*z1/2 - K1*K3*z3/2 + K1*z1 + K1*z2/2 + K1*z3/2 - K2*K3*z2/2 - K2*K3*z3/2 + K2*z1/2 + K2*z2 + K2*z3/2 + K3*z1/2 + K3*z2/2 + K3*z3 - z1 - z2 - z3 - (K1**2*K2**2*z1**2 + 2*K1**2*K2**2*z1*z2 + K1**2*K2**2*z2**2 - 2*K1**2*K2*K3*z1**2 - 2*K1**2*K2*K3*z1*z2 - 2*K1**2*K2*K3*z1*z3 + 2*K1**2*K2*K3*z2*z3 - 2*K1**2*K2*z1*z2 + 2*K1**2*K2*z1*z3 - 2*K1**2*K2*z2**2 - 2*K1**2*K2*z2*z3 + K1**2*K3**2*z1**2 + 2*K1**2*K3**2*z1*z3 + K1**2*K3**2*z3**2 + 2*K1**2*K3*z1*z2 - 2*K1**2*K3*z1*z3 - 2*K1**2*K3*z2*z3 - 2*K1**2*K3*z3**2 + K1**2*z2**2 + 2*K1**2*z2*z3 + K1**2*z3**2 - 2*K1*K2**2*K3*z1*z2 + 2*K1*K2**2*K3*z1*z3 - 2*K1*K2**2*K3*z2**2 - 2*K1*K2**2*K3*z2*z3 - 2*K1*K2**2*z1**2 - 2*K1*K2**2*z1*z2 - 2*K1*K2**2*z1*z3 + 2*K1*K2**2*z2*z3 + 2*K1*K2*K3**2*z1*z2 - 2*K1*K2*K3**2*z1*z3 - 2*K1*K2*K3**2*z2*z3 - 2*K1*K2*K3**2*z3**2 + 4*K1*K2*K3*z1**2 + 4*K1*K2*K3*z1 *
                                                                                                                                                                                                       z2 + 4*K1*K2*K3*z1*z3 + 4*K1*K2*K3*z2**2 + 4*K1*K2*K3*z2*z3 + 4*K1*K2*K3*z3**2 + 2*K1*K2*z1*z2 - 2*K1*K2*z1*z3 - 2*K1*K2*z2*z3 - 2*K1*K2*z3**2 - 2*K1*K3**2*z1**2 - 2*K1*K3**2*z1*z2 - 2*K1*K3**2*z1*z3 + 2*K1*K3**2*z2*z3 - 2*K1*K3*z1*z2 + 2*K1*K3*z1*z3 - 2*K1*K3*z2**2 - 2*K1*K3*z2*z3 + K2**2*K3**2*z2**2 + 2*K2**2*K3**2*z2*z3 + K2**2*K3**2*z3**2 + 2*K2**2*K3*z1*z2 - 2*K2**2*K3*z1*z3 - 2*K2**2*K3*z2*z3 - 2*K2**2*K3*z3**2 + K2**2*z1**2 + 2*K2**2*z1*z3 + K2**2*z3**2 - 2*K2*K3**2*z1*z2 + 2*K2*K3**2*z1*z3 - 2*K2*K3**2*z2**2 - 2*K2*K3**2*z2*z3 - 2*K2*K3*z1**2 - 2*K2*K3*z1*z2 - 2*K2*K3*z1*z3 + 2*K2*K3*z2*z3 + K3**2*z1**2 + 2*K3**2*z1*z2 + K3**2*z2**2)**0.5/2)/(K1*K2*K3*z1 + K1*K2*K3*z2 + K1*K2*K3*z3 - K1*K2*z1 - K1*K2*z2 - K1*K2*z3 - K1*K3*z1 - K1*K3*z2 - K1*K3*z3 + K1*z1 + K1*z2 + K1*z3 - K2*K3*z1 - K2*K3*z2 - K2*K3*z3 + K2*z1 + K2*z2 + K2*z3 + K3*z1 + K3*z2 + K3*z3 - z1 - z2 - z3)
    return V_over_F

# %% Data
    
# Material density (lb∕ft^3)
rho_Mdict = {'Carbon steel': 490,
            'Low-alloy steel': None,
            'Stainless steel 304': 499.4,
            'Stainless steel 316': 499.4,
            'Carpenter 20CB-3': None,
            'Nickel-200': None,
            'Monel-400': None,
            'Inconel-600': None,
            'Incoloy-825': None,
            'Titanium': None}

# Vessel Material
F_Mdict = {'Carbon steel': 1.0,
           'Low-alloy steel': 1.2,
           'Stainless steel 304': 1.7,
           'Stainless steel 316': 2.1,
           'Carpenter 20CB-3': 3.2,
           'Nickel-200': 5.4,
           'Monel-400': 3.6,
           'Inconel-600': 3.9,
           'Incoloy-825': 3.7,
           'Titanium': 7.7}


# %% Flash

class Flash(Unit):
    """Create an equlibrium based flash drum with the option of having light non-keys and heavy non-keys completly separate into their respective phases. Design procedure is based on heuristics by Wayne D. Monnery & William Y. Svrcek [1]. The f.o.b cost is based on classical Guthrie's correlation [2]:

    :math:`C_{fob}^{2007} = 6500 V^{0.62} (0.1 < V < 70 m^3)`

    **Parameters**

        **P:** [float] Operating pressure (Pa)

        **User defines one:**
            * **Qin:** [float] Energy input (kJ/hr)
            * **V:** [float] Overall molar vapor fraction
            * **T:** [float] Operating temperature (K)
            * **x:** [array_like] Molar composition of liquid (for binary mixture)
            * **y:** [array_like] Molar composition of vapor (for binary mixture)     
    
        **Optional**
    
            **species_IDs:** tuple[str] IDs of equilibrium species
            
            **LNK**: tuple[str] Light non-keys
        
            **HNK**: tuple[str] Heavy non-keys

    **ins**
    
        [0] Input stream
        
    **outs**
    
        [0] Vapor product
        
        [1] Liquid product

    **References**

        [1] Design Two-Phase Separators Within the Right Limits, Chemical Engineering Progress Oct, 1993
    
        [2] I.K. Kookos, Preliminary Chemical Process Synthesis and Design, Tziolas Publishing, Thessalonika, Greece, 2008 (book in Greek).
    
    **Examples**
    
        :doc:`Flash Example`
    
    """
    _units = {'Vertical vessel weight': 'lb',
              'Horizontal vessel weight': 'lb',
              'Length': 'ft',
              'Diameter': 'ft',
              'Weight': 'lb',
              'Wall thickness': 'in'}
    _has_power_utility = False
    _N_heat_utilities = 0
    
    # Column material factor
    _F_Mstr = 'Carbon steel'
    _F_M = 1 
    
    #: If a vacuum system is needed, it will choose one according to this preference.
    vacuum_system_preference = 'Liquid-ring pump'
    
    #: [str] 'Horizontal', 'Vertical', or 'Default'
    SepType = 'Default'
    
    #: [float] Time it takes to raise liquid to half full (min)
    HoldupTime = 15  
    
    #: [float] Time it takes to reach from normal to maximum liquied level (min)
    SurgeTime = 7.5
    
    #: [bool] True if using a mist eliminator pad
    Mist = False
    
    _kwargs = {'species_IDs': None,  # Equilibrium species
               'LNK': None,   # light non-keys
               'HNK': None,   # heavy non-keys
               'V': None,   # Vapor fraction
               'T': None,   # Operating temperature
               'Qin': None, # Energy input
               'P': None,   # Operating Pressure
               'y': None,   # Vapor composition (of working species)
               'x': None}   # Liquid composition (of working species)

    _bounds = {'Vertical vessel weight': (4200, 1e6),
               'Horizontal vessel weight': (1e3, 9.2e5),
               'Diameter': (3, 21),
               'Vertical vessel Length': (12, 40)}

    @property
    def vessel_material(self):
        return self._F_Mstr
    @vessel_material.setter
    def vessel_material(self, vessel_material):
        try:
            self._F_M = F_Mdict[vessel_material]
        except KeyError:
            dummy = str(F_Mdict.keys())[11:-2]
            raise ValueError(f"Vessel material must be one of the following: {dummy}")
        self._F_Mstr = vessel_material  

    def _init(self):
        vap, liq = self.outs
        vap.phase = 'g'
        liq.phase = 'l'
        self._heat_exchanger = he = HXutility(None, outs=None) 
        self._heat_utilities = he._heat_utilities
        he._ins = self._ins
        self._mixedstream = ms = MixedStream(None)
        he._outs[0] = ms

    def _setup(self):
        self._has_hx = self._kwargs['Qin'] != 0
        if self._kwargs['P'] < 101325 and not self._power_utility:
            self._power_utility = PowerUtility()

    def _run(self):
        # Unpack
        vap, liq = self.outs
        feed = self.ins[0]

        # Vapor Liquid Equilibrium
        ms = self._mixedstream
        ms.empty()
        ms.liquid_mol = feed.mol
        ms.T = feed.T
        ms.VLE(**self._kwargs)

        # Set Values
        vap.mol = ms.vapor_mol
        liq.mol = ms.liquid_mol
        vap.T = liq.T = ms.T
        vap.P = liq.P = ms.P

    def _design(self):
        # Set horizontal or vertical vessel
        if self.SepType == 'Default':
            if self.outs[0].massnet/self.outs[1].massnet > 0.2:
                self.SepType = 'Vertical'
            else:
                self.SepType = 'Horizontal'

        # Run Vertical or Horizontal:
        if self.SepType == 'Vertical':
            out = self._vertical()
        elif self.SepType == 'Horizontal' or '(USDA) Biodiesel':
            out = self._horizontal()
        else:
            raise ValueError( f"SepType must be either 'Horizontal' or 'Vertical', not '{self.SepType}'")
        if self._has_hx: self._heat_exchanger._design()
        return out

    def _cost(self):
        Cost = self._results['Cost']
        Design = self._results['Design']
        W = Design['Weight']
        D = Design['Diameter']
        L = Design['Length']
        CE = self.CEPCI
        type_ = self.SepType
        
        # C_v: Vessel cost
        # C_pl: Platforms and ladders cost
        if type_ == 'Vertical':
            C_v = exp(7.0132 + 0.18255*ln(W) + 0.02297*ln(W)**2)
            C_pl = 361.8*D**0.7396*L**0.70684
        elif type_ == 'Horizontal':
            C_v = exp(8.9552 - 0.2330*ln(W) + 0.04333*ln(W)**2)
            C_pl = 2005*D**0.20294
        else:
            ValueError(f"SepType ({type_}) must be either 'Vertical', 'Horizontal' or 'Default'.")
            
        Cost['Vessel'] = CE/500*(C_v+C_pl)
        if self._has_hx:
            Cost['Heat exchanger'] = self._heat_exchanger._cost()['Heat exchanger']
        self._cost_vacuum()

    def _cost_vacuum(self):
        P = self._kwargs['P']
        if not P or P > 101320:
            return 
        
        r = self._results
        D = r['Design']
        C = r['Cost']
        vol = 0.02832 * np.pi * D['Length'] * (D['Diameter']/2)**2
        
        # If vapor is condensed, assume vacuum system is after condenser
        vapor = self.outs[0]
        hx = vapor.sink
        if isinstance(hx, HX):
            index = hx._ins.index(vapor)
            stream = hx._outs[index]
            if isinstance(stream, MixedStream):
                massflow = stream.vapor_massnet
                volflow = stream.vapor_volnet
            else:
                if stream.phase == 'g':
                    massflow = stream.massnet
                    volflow = stream.volnet
                else:
                    massflow = 0
                    volflow = 0
        else:
            massflow = vapor.massnet
            volflow = vapor.volnet
        
        power, cost = vacuum_system(massflow, volflow,
                                    P, vol,
                                    self.vacuum_system_preference)
        C['Liquid-ring pump'] = cost
        self._power_utility(power)

    def _design_parameters(self):
        # Retrieve run_args and properties
        vap, liq = self.outs[0:2]
        v = vap.quantity
        l = liq.quantity
        rhov = v('rho').to('lb/ft3').magnitude  # VapDensity
        rhol = l('rho').to('lb/ft3').magnitude  # LLiqDensity
        P = l('P').to('psi').magnitude  # Pressure

        # SepType (str), HoldupTime (min), SurgeTime (min), Mist (bool)
        SepType = self.SepType
        Th = self.HoldupTime
        Ts = self.SurgeTime
        Mist = self.Mist

        # Calculate the volumetric flowrate
        Qv = v('volnet').to('ft3/s').magnitude
        Qll = l('volnet').to('ft3/min').magnitude

        # Calculate Ut and set Uv
        K = Kvalue(P)

        # Adjust K value
        if not Mist and SepType == 'Vertical':
            K = K/2

        xs = liq.massfrac
        for s, x in zip(liq._species, xs):
            if x > 0.1:
                # TODO: Verify unifac groups to determine if glycol and/or amine
                dct = s.UNIFAC_groups
                keys = dct.keys()
                if 14 in keys:
                    OHs = dct[14]
                else:
                    OHs = 0
                if OHs > 1:
                    K = K*0.6
                elif any([i in keys for i in (28, 29, 30, 31, 32, 33)]):
                    K = K*0.8
                else:
                    OHs = 0

        Ut = K*((rhol - rhov) / rhov)**0.5
        Uv = 0.75*Ut
        # Calculate Holdup and Surge volume
        Vh = Th*Qll
        Vs = Ts*Qll
        return rhov, rhol, P, Th, Ts, Mist, Qv, Qll, Ut, Uv, Vh, Vs

    def _vertical(self):
        rhov, rhol, P, Th, Ts, Mist, Qv, Qll, Ut, Uv, Vh, Vs = self._design_parameters()

        # Calculate internal diameter, Dvd
        Dvd = (4.0*Qv/(pi*Uv))**0.5
        if int(Mist) == 1:
            D = FinalValue(Dvd + 0.4)
        else:
            D = FinalValue(Dvd)

        # Obtaining low liquid level height, Hlll
        Hlll = 0.5
        if P < 300:
            Hlll = 1.25

        # Calculate the height from Hlll to  Normal liq level, Hnll
        Hh = Vh/(pi/4.0*Dvd**2)
        if Hh < 1.0:
            Hh = 1.0
            Hh = FinalValue(Hh)

        # Calculate the height from Hnll to  High liq level, Hhll
        Hs = Vs/(pi/4.0*Dvd**2)
        if Hs < 0.5:
            Hs = 0.5
            Hs = FinalValue(Hs)

        # Calculate dN
        Qm = Qll + Qv
        lamda = Qll/Qm
        rhoM = rhol*lamda + rhov*(1-lamda)
        dN = (4*Qm/(pi*60.0/(rhoM**0.5)))**0.5
        dN = FinalValue(dN)

        # Calculate Hlin, assume with inlet diverter
        Hlin = 1.0 + dN

        # Calculate the vapor disengagement height
        Hv = 0.5*Dvd
        if int(Mist) == 1:
            Hv2 = 2.0 + dN/2.0
        else:
            Hv2 = 3.0 + dN/2.0

        if Hv2 < Hv:
            Hv = Hv2
        Hv = ceil(Hv)

        # Calculate total height, Ht
        Hme = 0.0
        if int(Mist) == 1:
            Hme = 1.5
        Ht = Hlll + Hh + Hs + Hlin + Hv + Hme
        Ht = FinalValue(Ht)

        # Check if LD is between 1.5 and 6.0
        converged = False
        LD = Ht/D
        while not converged:
            if (LD < 1.5):
                D = D - 0.5
                converged = False
            elif (LD > 6.0):
                D = D + 0.5
                converged = False
            else:
                converged = True
            LD = Ht/D

        # Calculate Vessel weight and wall thickness
        rho_M = rho_Mdict[self._F_Mstr]
        VW, VWT = VesselWeightAndWallThickness(P, D, Ht, rho_M)

        # Find maximum and normal liquid level
        # Hhll = Hs + Hh + Hlll
        # Hnll = Hh + Hlll

        # Results
        d = self._results['Design']
        self._checkbounds('Vertical vessel weight', VW)
        self._checkbounds('Vertical vessel length', Ht)
        d['SepType'] = 'Vertical'
        d['Length'] = Ht     # ft
        d['Diameter'] = D    # ft
        d['Weight'] = VW     # lb
        d['Wall thickness'] = VWT  # in
        return d
        
    def _horizontal(self):
        rhov, rhol, P, Th, Ts, Mist, Qv, Qll, Ut, Uv, Vh, Vs = self._design_parameters()
        maxIter = 50

        # Initialize LD
        if P > 0 and P <= 264.7:
            LD = 1.5/250.0*(P-14.7)+1.5
        elif P > 264.7 and P <= 514.7:
            LD = 1.0/250.0*(P-14.7)+2.0
        elif P > 514.7:
            LD = 5.0

        D = (4.0*(Vh+Vs)/(0.6*pi*LD))**(1.0/3.0)
        D = round(D)
        if D <= 4.0:
            D = 4.0

        outerIter = 0
        converged = False
        converged1 = False
        while not converged and outerIter < maxIter:
            outerIter += 1
            At = pi*(D**2)/4.0

            # Calculate Lower Liquid Area
            Hlll = round(0.5*D + 7.0)  
            Hlll = Hlll/12.0 # D is in ft but Hlll is in inches
            X = Hlll/D
            Y = HNATable(1, X)
            Alll = Y*At

            # Calculate the Vapor  disengagement area, Av
            Hv = 0.2*D
            if int(Mist) == 1 and Hv <= 2.0:
                Hv = 2.0
            else:
                if Hv <= 1.0:
                    Hv = 1.0
            X = Hv/D
            Y = HNATable(1, X)
            Av = Y*At

            # Calculate minimum length fo surge and holdup
            L = (Vh + Vs)/(At - Av - Alll)
            # Calculate liquid dropout
            Phi = Hv/Uv
            # Calculate actual vapor velocity
            Uva = Qv/Av
            # Calculate minimum length for vapor disengagement
            Lmin = Uva*Phi
            Li = L

            if L < Lmin:
                Li = Lmin

            if L < 0.8*Lmin:
                sign = 1.0
                needToIter = True
            elif L > 1.2*Lmin:
                if int(Mist) == 1 and Hv <= 2.0:
                    Hv = 2.0
                    Li = L
                elif not int(Mist) == 0 and Hv <= 1.0:
                    Hv = 1.0
                    Li = L
                else:
                    sign = -1.0
                    needToIter = True
            else:
                sign = -1.0
                needToIter = False

            if needToIter:
                innerIter = 0
                while not converged1 and innerIter < maxIter:
                    innerIter += 1
                    Hv = Hv + sign*0.5
                    if int(Mist) == 1 and Hv <= 2.0:
                        Hv = 2.0
                    if int(Mist) == 0 and Hv <= 1.0:
                        Hv = 1.0
                    X = Hv/D
                    Y = HNATable(1, X)
                    Av = Y*At

                    X = Hlll/D
                    Y = HNATable(1, X)
                    Alll = Y*At

                    Li = (Vh + Vs)/(At - Av - Alll)
                    Phi = Hv/Uv
                    Uva = Qv/Av
                    Lmin = Uva*Phi
                    if Li < 0.8*Lmin:
                        sign = 1.0
                        converged1 = False
                    elif Li > 1.2*Lmin:
                        sign = -1.0
                        converged1 = False
                    else:
                        Li = Li
                        converged1 = True
            
            L = Li
            LD = L/D
            # Check LD
            if LD < 1.2:
                if D <= 4.0:
                    D = D
                    converged = True
                else:
                    D = D - 0.5
                    converged = False

            if LD > 7.2:
                D = D + 0.5
                converged = False
            else:
                converged = True

        # Recalculate LD so it lies between 1.5 - 6.0
        converged = False
        while not converged:
            LD = L / D
            if (LD < 1.5) and D <= 4.0:
                L = L + 0.5
                converged = False
            elif LD < 1.5:
                D = D - 0.5
                converged = False
            elif (LD > 6.0):
                D = D + 0.5
                converged = False
            else:
                converged = True

        # Calculate vessel weight and wall thickness
        rho_M = rho_Mdict[self._F_Mstr]
        VW, VWT = VesselWeightAndWallThickness(P, D, L, rho_M)

        # # To check minimum Hv value
        # if int(Mist) == 1 and Hv <= 2.0:
        #     Hv = 2.0
        # if int(Mist) == 0 and Hv <= 1.0:
        #     Hv = 1.0

        # Calculate normal liquid level and High liquid level
        # Hhll = D - Hv
        # if (Hhll < 0.0):
        #     Hhll = 0.0
        # Anll = Alll + Vh/L
        # X = Anll/At
        # Y = HNATable(2, X)
        # Hnll = Y*D
        
        # Results
        d = self._results['Design']
        self._checkbounds('Horizontal vessel weight', VW)
        d['SepType'] = 'Horizontal'
        d['Length'] = L  # ft
        d['Diameter'] = D  # ft
        d['Weight'] = VW  # lb
        d['Wall thickness'] = VWT  # in
        return d
    
    

# %% Special


class SplitFlash(Flash):
    line = 'Flash'
    
    _kwargs = {'split': 1,  # component split fractions
               'T': 298.15,  # operating temperature (K)
               'P': 101325,
               'Qin': None}  # operating pressure (Pa)
    
    def _setup(self):
        top, bot = self.outs
        _, Tf, Pf, Qin = self._kwargs.values()
        self._has_hx = Qin != 0
        if self._kwargs['P'] < 101325 and not self._power_utility:
            self._power_utility = PowerUtility()
        bot.T = top.T = Tf
        bot.P = top.P = Pf

    def _run(self):
        split = self._kwargs['split']
        top, bot = self.outs
        net_mol = self._mol_in
        top.mol = net_mol*split
        bot.mol = net_mol - top.mol

    def _design(self):
        if self._has_hx:
            ms = MixedStream.sum(self._mixedstream, self.outs)
            self._heat_exchanger.outs = ms
        super()._design()

class PartitionFlash(Flash):
    """Create a PartitionFlash Unit that approximates outputs based on molar partition coefficients."""
    
    _kwargs = {'species_IDs': None,  # Partition species
               'Ks': None,           # Partition coefficients
               'LNK': None,          # light non-keys
               'HNK': None,          # heavy non-keys
               'P': None,            # Operating Pressure
               'T': None}            # Operating temperature
    
    _has_hx = True
    
    def _setup(self):
        top, bot = self.outs
        species_IDs, Ks, LNK, HNK, P, T = self._kwargs.values()
        index = top._IDs.index
        if P < 101325 and not self._power_utility:
            self._power_utility = PowerUtility()
        self._args = ([index(i) for i in species_IDs],
                      [index(i) for i in LNK],
                      [index(i) for i in HNK],
                       np.asarray(Ks), P, T)
        
    def _set_phases(self):
        top, bot = self.outs
        top.phase = 'g'
        bot.phase = 'l'
    
    def _run(self):
        feed = self.ins[0]
        ph1, ph2 = self.outs
        pos, LNK, HNK, Ks, P, T = self._args
        mol_in = feed.mol[pos]
        ph1.T = ph2.T = (T if T else feed.T)
        ph1.P = ph2.P = (P if P else feed.P)
        ph1.mol[LNK] = feed.mol[LNK]
        ph1.mol[HNK] = feed.mol[HNK]
        
        # Get zs and Ks to solve rashford rice equation
        molnet = sum(mol_in)
        zs = mol_in/molnet
        if hasattr(self, '_V'):
            self._V = V = newton(flash_error, self._V, args=(zs, Ks))
        elif len(zs) > 3:
            self._V = V = brentq(flash_error, 0, 1, (zs, Ks))
        else:
            V = analytical_flash_solution(zs, Ks)
        x = zs/(1 + V*(Ks-1))
        y = x*Ks
        ph1.mol[pos] = molnet*V*y
        ph2.mol[pos] = mol_in - ph1.mol[pos]
    
    _design = SplitFlash._design
    

class RatioFlash(Flash):
    _N_heat_utilities = 1
    _kwargs = {'Kspecies': None,  # list of species that correspond to Ks
               'Ks': None,  # array of molar ratio partition coefficinets,
               # Ks = y/x, where y and x are molar ratios of two different phases
               'top_solvents': [],  # list of species that correspond to top_split
               'top_split': [],  # list of splits for top_solvents
               'bot_solvents': [],  # list of species that correspond to bot_split
               'bot_split': []} # list of splits for bot_solvents  

    def _run(self):
        feed = self.ins[0]
        top, bot = self.outs
        kwargs = self._kwargs
        (Kspecies, Ks, top_solvents, top_split,
         bot_solvents, bot_split) = kwargs.values()
        sp_index = feed._IDs.index
        Kindex = [sp_index(s) for s in Kspecies]
        top_index = [sp_index(s) for s in top_solvents]
        bot_index = [sp_index(s) for s in bot_solvents]
        top_mol = top.mol; bot_mol = bot.mol; feed_mol = feed.mol
        top_mol[top_index] = feed_mol[top_index]*top_split
        bot_mol[top_index] = feed_mol[top_index]-top_mol[top_index]
        bot_mol[bot_index] = feed_mol[bot_index]*bot_split
        top_mol[bot_index] = feed_mol[bot_index]-bot_mol[bot_index]
        topnet = sum(top_mol[top_index])
        botnet = sum(bot_mol[bot_index])
        molnet = topnet+botnet
        top_mol[Kindex] = Ks * topnet * feed_mol[Kindex] / molnet  # solvent * mol ratio
        bot_mol[Kindex] = feed_mol[Kindex] - top_mol[Kindex]
        top.T, top.P = feed.T, feed.P
        bot.T, bot.P = feed.T, feed.P


# %% Single Component

class Evaporator_PQin(Unit):

    _kwargs = {'component': 'Water',  # ID of specie
               'Qin': 0,
               'P': 101325}

    def _setup(self):
        component, P = (self._kwargs[i] for i in ('component', 'P'))
        vapor, liquid = self.outs[:2]
        self._cached = cached = {}
        Tf = getattr(vapor._species, component).Tsat(P)

        # Set-up streams for energy balance
        vapor.empty()
        index = vapor._IDs.index(component)
        vapor.mol[index] = 1
        vapor.__dict__.update(phase='g', T=Tf, P=P)

        liquid.empty()
        liquid.mol[index] = 1
        liquid.__dict__.update(phase='l', T=Tf, P=P)

        no_ph_ch = type(liquid)(species=vapor.species)
        cached['vapor_H'] = vapor.H
        cached['liquid_H'] = liquid.H
        cached['no_ph_ch'] = no_ph_ch

    def _run(self):
        feed = self.ins[0]
        component, Qin = (self._kwargs[i] for i in ('component', 'Qin'))
        vapor, liquid = self.outs[:2]
        cached = self._cached

        # Set-up streams for energy balance
        vapor_H = cached['vapor_H']
        liquid_H = cached['liquid_H']
        no_ph_ch = cached['no_ph_ch']
        no_ph_ch.mol = copy.copy(feed.mol)
        index = no_ph_ch._IDs.index(component)
        no_ph_ch.mol[index] = 0

        no_ph_ch_H = no_ph_ch.H
        feed_H = feed.H
        # Optional if Qin also comes from condensing a side stream
        if len(self.ins) == 2:
            boiler_liq = self.outs[2]
            boiler_vap = self.ins[1]
            boiler_liq.copylike(boiler_vap)
            boiler_liq.phase = 'l'
            Qin = Qin + boiler_vap.H - boiler_liq.H

        # Energy balance to find vapor fraction
        def f(v):
            vapor_Hf = vapor_H*(v*feed.mol[index])
            liquid_Hf = liquid_H*((1-v)*feed.mol[index])
            return (vapor_Hf + liquid_Hf + no_ph_ch_H) - feed_H - Qin

        # Final vapor fraction
        v = newton(f, 0.5, tol=0.001)
        liquid.mol = no_ph_ch.mol
        if v < 0:
            v = 0
        elif v > 1:
            v = 1
        vapor.mol[index] = v*feed.mol[index]
        liquid.mol[index] = (1-v)*feed.mol[index]
        self._Qin = Qin
        self._V = v


class Evaporator_PV(Unit):
    _N_heat_utilities = 1
    _kwargs = {'component': 'Water',
               'V': 0.5,
               'P': 101325}

    def _setup(self):
        vapor, liquid = self.outs
        component, P = (self._kwargs[i] for i in ('component', 'P'))
        Tf = getattr(vapor._species, component).Tsat(P)
        vapor.__dict__.update(phase='g', T=Tf, P=P)
        liquid.__dict__.update(phase='l', T=Tf, P=P)

    def _run(self):
        feed = self.ins[0]
        component, V = (self._kwargs[i] for i in ('component', 'V'))
        vapor, liquid = self.outs
        index = vapor._IDs.index(component)
        vapor.mol[index] = V*feed.mol[index]
        liquid.mol = copy.copy(feed.mol)
        liquid.mol[index] = (1-V)*feed.mol[index]

        