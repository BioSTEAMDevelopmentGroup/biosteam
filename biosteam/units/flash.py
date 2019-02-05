# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 16:21:56 2018

@author: yoelr
"""
from biosteam import Unit, MixedStream
from biosteam.exceptions import DesignError
from math import pi, ceil
import copy
import numpy as np
from scipy.optimize import brentq, newton
from thermo import activity
from biosteam.units.hx import HXutility
from biosteam.units.design_tools import (HNATable, FinalValue,
                          VesselWeightAndWallThickness, Kvalue)

exp = np.exp
ln = np.log

# USDA Biodiesel: 
#V = np.pi*D**2/4*Ht*0.3048**3
#if not (0.1 < V < 70):
#    raise DesignError(f"Volume is out of bounds for costing")
#lambda V, CE: CE*13*V**0.62 # V (m3)


# %% Flash solutions

flash_error = activity.Rachford_Rice_flash_error

def analytical_flash_solution(zs, Ks):
    """Solution for 2 or 3 component flash vessel."""
    # Not currently in used, so no comments
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
            'Stainless steel 304': None,
            'Stainless steel 316': None,
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

        **User defines two:**
            * **Qin:** [float] Energy input (kJ/hr)
            * **V:** [float] Overall molar vapor fraction
            * **T:** [float] Operating temperature (K)
            * **P:** [float] Operating pressure (Pa)
            * **x:** [array_like] Molar composition of liquid (for binary mixture)
            * **y:** [array_like] Molar composition of vapor (for binary mixture)     
    
        **Optional**
    
            **specie_IDs:** tuple[str] IDs of equilibrium species
            
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
    
    """
    _N_heat_util = 0
    
    # Column material factor
    _F_Mstr = 'Carbon steel'
    _F_M = 1 
    
    #: [str] 'Horizontal', 'Vertical', or 'Default'
    SepType = 'Default'
    
    #: [float] Time it takes to raise liquid to half full (min)
    HoldupTime = 15  
    
    #: [float] Time it takes to reach from normal to maximum liquied level (min)
    SurgeTime = 7.5
    
    #: [bool] True if using a mist eliminator pad
    Mist = False
    
    kwargs = {'specie_IDs': None,  # Equilibrium species
              'LNK': None,   # light non-keys
              'HNK': None,   # heavy non-keys
              'V': None,   # Vapor fraction
              'T': None,   # Operating temperature
              'Qin': None, # Energy input
              'P': None,   # Operating Pressure
              'y': None,   # Vapor composition (of working species)
              'x': None}   # Liquid composition (of working species)

    bounds = {'Vertical vessel weight': (4200, 1e6),
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

    def _setup(self):
        vap, liq = self.outs
        vap.phase = 'g'
        liq.phase = 'l'
        self._cached = cached = {}
        cached['mixed stream'] = ms = MixedStream()
        self._heat_exchanger = he = HXutility(self.ID+' hx', None) 
        self.heat_utilities = he.heat_utilities
        he._ins = self._ins
        he-ms

    def _run(self):
        # Unpack
        vap, liq = self.outs
        feed = self.ins[0]
        kwargs = self.kwargs

        # Vapor Liquid Equilibrium
        VLE_kwargs = copy.copy(kwargs)
        ms = self._cached['mixed stream']
        ms.empty()
        ms.liquid_mol = feed.mol
        ms.T = feed.T
        ms.VLE(**VLE_kwargs)

        # Set Values
        vap.mol = ms.vapor_mol
        liq.mol = ms.liquid_mol
        vap.T = liq.T = ms.T
        vap.P = liq.P = ms.P

    def _simple_run(self):
        """Approximate outputs based on molar partition coefficients and overall molar fraction of top phase (0th output)."""
        ph1, ph2 = self.outs
        _mol_in = ph1.mol + ph2.mol

        # Account for light and heave keys
        HNK = ph1.mol == 0
        LNK = ph2.mol == 0
        ph1.mol[LNK] = _mol_in[LNK]
        ph2.mol[HNK] = _mol_in[HNK]

        # Do not include heavy or light non keys
        pos = ~HNK & ~LNK

        ph1_mol = ph1.mol[pos]
        ph2_mol = ph2.mol[pos]
        ph1_molnet = sum(ph1_mol)
        ph2_molnet = sum(ph2_mol)
        molnet = ph1_molnet + ph2_molnet

        # Get zs and Ks to solve rashford rice equation
        zs = _mol_in[pos]/molnet
        if ph1_molnet == 0:
            ph2.mol[pos] = _mol_in[pos]
            return
        if ph2_molnet == 0:
            ph1.mol[pos] = _mol_in[pos]
            return
        ph1_molfrac = ph1_mol/ph1_molnet
        ph2_molfrac = ph2_mol/ph2_molnet
        Ks = ph1_molfrac/ph2_molfrac

        if len(zs) > 3:
            V = brentq(flash_error, 0, 1, (zs, Ks))
        else:
            V = analytical_flash_solution(zs, Ks)

        ph1.mol[pos] = V*molnet*ph1_molfrac
        ph2.mol[pos] = _mol_in[pos] - ph1.mol[pos]

    def _operation(self):
        self._heat_exchanger._operation()

    def _design(self):
        """
        * 'Length': (ft)
        * 'Diameter': (ft)
        * 'Weight': (lb)
        * 'Wall thickness': (in)
        """
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
        
        self._heat_exchanger._design()
        return out

    def _cost(self):
        """
        * 'Vessel': (USD)
        * 'Heat exchanger': (USD)
        """
        Cost = self.results['Cost']
        Design = self.results['Design']
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
        Cost['Heat exchanger'] = self._heat_exchanger._cost()['Heat exchanger']
        return Cost

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

        index, specie_IDs = liq.nonzero_species
        xs = liq.massfrac[index]
        for sp_ID, x in zip(specie_IDs, xs):
            if x > 0.1:
                specie = getattr(liq._species, sp_ID)
                # TODO: Verify unifac groups to determine if glycol and/or amine
                dct = specie.UNIFAC_groups
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

        # Calculate Vessel weight and wall thickness
        rho_M = rho_Mdict[self._F_Mstr]
        VW, VWT = VesselWeightAndWallThickness(P, D, Ht, rho_M)

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

        # Find maximum and normal liquid level
        # Hhll = Hs + Hh + Hlll
        # Hnll = Hh + Hlll

        # Results
        d = self.results['Design']
        d.boundscheck('Vertical vessel weight', VW, 'lb')
        d.boundscheck('Vertical vessel length', Ht, 'lb')
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
                    if int(Mist) == 1 and Hv <= 2.0:   # baru tambah Bang !
                        Hv = 2.0
                    if int(Mist) == 0 and Hv <= 1.0:   # baru tambah Bang !
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
            if LD < (0.8*1.5):
                if D <= 4.0:
                    D = D
                    converged = True
                else:
                    D = D - 0.5
                    converged = False

            if LD > (1.2*6.0):
                D = D + 0.5
                converged = False
            else:
                converged = True

        # Calculate vessel weight and wall thickness
        rho_M = rho_Mdict[self._F_Mstr]
        VW, VWT = VesselWeightAndWallThickness(P, D, L, rho_M)

        # To check minimum Hv value
        if int(Mist) == 1 and Hv <= 2.0:
            Hv = 2.0
        if int(Mist) == 0 and Hv <= 1.0:
            Hv = 1.0

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

        # Calculate normal liquid level and High liquid level
        # Hhll = D - Hv
        # if (Hhll < 0.0):
        #     Hhll = 0.0
        # Anll = Alll + Vh/L
        # X = Anll/At
        # Y = HNATable(2, X)
        # Hnll = Y*D
        
        # Results
        d = self.results['Design']
        d.boundscheck('Horizontal vessel weight', VW, 'lb')
        d['SepType'] = 'Horizontal'
        d['Length'] = L  # ft
        d['Diameter'] = D  # ft
        d['Weight'] = VW  # lb
        d['Wall thickness'] = VWT  # in
        return d

# %% Estimate


class EstimateFlash(Flash):
    _N_heat_util = 1
    kwargs = {'Kspecies': [],  # list of species that correspond to Ks
              'Ks': [],  # list of molar ratio partition coefficinets,
              # Ks = y/x, where y and x are molar ratios of two different phases
              'top_solvents': [],  # list of species that correspond to top_split
              'top_split': [],  # list of splits for top_solvents
              'bot_solvents': [],  # list of species that correspond to bot_split
              'bot_split': []}  # list of splits for bot_solvents

    def _setup(self):
        self._heat_exchanger = he = HXutility(self.ID+' hx', outs_ID=()) 
        self.heat_utilities = he.heat_utilities
        he._ins = self._ins
        he.outs = self._outs

    def _run(self):
        feed = self.ins[0]
        top, bot = self.outs
        Kspecies, Ks, top_solvents, top_split, bot_solvents, bot_split = (
            (self.kwargs[i] for i in ('Kspecies', 'Ks', 'top_solvents', 'top_split', 'bot_solvents', 'bot_split')))
        Kindex = [feed._ID_index[specie] for specie in Kspecies]
        top_index = [feed._ID_index[specie] for specie in top_solvents]
        bot_index = [feed._ID_index[specie] for specie in bot_solvents]
        top.mol[top_index] = feed.mol[top_index]*top_split
        bot.mol[top_index] = feed.mol[top_index]-top.mol[top_index]
        bot.mol[bot_index] = feed.mol[bot_index]*bot_split
        top.mol[bot_index] = feed.mol[bot_index]-bot.mol[bot_index]
        topnet = sum(top.mol[top_index])
        botnet = sum(bot.mol[bot_index])
        for i in range(len(Kindex)):
            index = Kindex[i]
            K = Ks[i]
            top.mol[index] = topnet * feed.mol[index] / \
                (topnet+botnet/K)  # solvent * mol ratio
        bot.mol[Kindex] = feed.mol[Kindex] - top.mol[Kindex]
        top.T, top.P = feed.T, feed.P
        bot.T, bot.P = feed.T, feed.P


# %% Single Component

class Flash_PQin(Unit):

    kwargs = {'component': 'Water',  # ID of specie
              'Qin': 0,
              'P': 101325}

    def _setup(self):
        # Unpack
        component, P = (self.kwargs[i] for i in ('component', 'P'))
        vapor, liquid = self.outs[:2]
        self._cached = cached = {}

        # Find final Temperature
        Tf = getattr(vapor._species, component).Tsat(P)

        # Set-up streams for energy balance
        vapor.empty()
        index = vapor.get_index(component)
        vapor.mol[index] = 1
        vapor.__dict__.update(phase='g', T=Tf, P=P)

        liquid.empty()
        liquid.mol[index] = 1
        liquid.__dict__.update(phase='l', T=Tf, P=P)

        no_ph_ch = type(liquid)()

        cached['vapor_H'] = vapor.H
        cached['liquid_H'] = liquid.H
        cached['no_ph_ch'] = no_ph_ch

    def _run(self):
        # Unpack
        feed = self.ins[0]
        component, Qin = (self.kwargs[i] for i in ('component', 'Qin'))
        vapor, liquid = self.outs[:2]
        cached = self._cached

        # Set-up streams for energy balance
        vapor_H = cached['vapor_H']
        liquid_H = cached['liquid_H']
        no_ph_ch = cached['no_ph_ch']
        no_ph_ch.mol = copy.copy(feed.mol)
        index = no_ph_ch.get_index(component)
        no_ph_ch.mol[index] = 0

        no_ph_ch_H = no_ph_ch.H
        feed_H = feed.H
        # Optional if Qin also comes from condensing a side stream
        if len(self.ins) == 2:
            boiler_liq = self.outs[2]
            boiler_vap = self.ins[1]
            boiler_liq.copy_like(boiler_vap)
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


class Flash_PV(Unit):
    _N_heat_util = 1
    kwargs = {'component': 'ID of specie',
              'V': 'Vapor fraction',
              'P': 'Operating pressure (Pa)'}

    def _setup(self):
        # Unpack
        vapor, liquid = self.outs
        component, P = (self.kwargs[i] for i in ('component', 'P'))

        # Find final Temperature
        Tf = getattr(vapor._species, component).Tsat(P)

        # Outputs
        vapor.__dict__.update(phase='g', T=Tf, P=P)
        liquid.__dict__.update(phase='l', T=Tf, P=P)

    def _run(self):
        # Unpack
        feed = self.ins[0]
        component, V = (self.kwargs[i] for i in ('component', 'V'))

        # Outputs
        vapor, liquid = self.outs
        index = vapor.get_index(component)
        vapor.mol[index] = V*feed.mol[index]
        liquid.mol = copy.copy(feed.mol)
        liquid.mol[index] = (1-V)*feed.mol[index]

    _operation = HXutility._operation

class P69_flash(Flash):
    kwargs = {'split': [],  # component split fractions
              'phase0': 'g',  # phase of 0th outs
              'phase1': 'l',  # phase of 1st outs
              'T': 298.15,  # operating temperature (K)
              'P': 101325}  # operating pressure (Pa)
    
    def _setup(self):
        # Unpack
        top, bot = self.outs
        top.phase, bot.phase, Tf, Pf = (
            self.kwargs[i] for i in ('phase0', 'phase1', 'T', 'P'))
        bot.T = top.T = Tf
        bot.P = top.P = Pf
        
        vap, liq = self.outs
        vap.phase = 'g'
        liq.phase = 'l'
        self._heat_exchanger = he = HXutility(self.ID+' hx', outs=None) 
        self.heat_utilities = he.heat_utilities
        he._ins = self._ins

    def _run(self):
        split = self.kwargs['split']
        top, bot = self.outs
        net_mol = self._mol_in
        top.mol = net_mol*split
        bot.mol = net_mol - top.mol

    def _operation(self):
        ms = MixedStream.sum(MixedStream(species=self.outs[0].species), self.outs)
        self._heat_exchanger.outs = ms
        self._heat_exchanger._operation()
        