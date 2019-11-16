# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 16:21:56 2018

@author: yoelr
"""
from .. import Stream, Unit, MixedStream, PowerUtility
from math import pi, ceil
import numpy as np
from scipy.optimize import brentq, newton
from ..thermo import activity
from .designtools import vacuum_system, HNATable, FinalValue, \
                          VesselWeightAndWallThickness, Kvalue
from ..utils.solvers import IQ_interpolation
from .._equilibrium import V_3N, V_2N, V_error
from ._splitter import Splitter
from ._hx import HX, HXutility
import biosteam as bst

exp = np.exp
ln = np.log

# USDA Biodiesel: 
#V = np.pi*D**2/4*Ht*0.3048**3
#if not (0.1 < V < 70):
#    raise DesignError(f"Volume is out of bounds for costing")
#lambda V, CE: CE*13*V**0.62 # V (m3)
__all__ = ('Flash', 'SplitFlash', 'RatioFlash')


# %% Data
    
# Material density (lb∕ft^3)
material_density = {'Carbon steel': 490,
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
material_factor = {'Carbon steel': 1.0,
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
    """Create an equlibrium based flash drum with the option of having light non-keys and heavy non-keys completly separate into their respective phases. Design procedure is based on heuristics by Wayne D. Monnery & William Y. Svrcek [1]. Purchase costs are based on correlations by Mulet et al. [2, 3] as compiled by Warren et. al. [4].

    Parameters
    ----------
    ins
        [0] Input stream
        
    outs
        [0] Vapor product
        
        [1] Liquid product
    Specify two:
        * **P:** Operating pressure (Pa)
        * **Q:** Energy input (kJ/hr)
        * **T:** Operating temperature (K)
        * **V:** Molar vapor fraction
        * **x:** Molar composition of liquid (for binary mixture)
        * **y:** Molar composition of vapor (for binary mixture)
    species_IDs=None : tuple, optional
        IDs of species in equilibrium.
    LNK=None : tuple[str], optional
        Light non-keys that remain as a vapor (disregards equilibrium).
    LNK=None : tuple[str], optional
        Heavy non-keys that remain as a liquid (disregards equilibrium).

    References
    ----------
    [1] "Design Two-Phase Separators Within the Right Limits", Chemical Engineering Progress Oct, 1993.

    [2] Mulet, A., A. B. Corripio, and L. B. Evans, “Estimate Costs of Pressure Vessels via Correlations,” Chem. Eng., 88(20), 145–150 (1981a).

    [3] Mulet, A., A.B. Corripio, and L.B.Evans, “Estimate Costs of Distillation and Absorption Towers via Correlations,” Chem. Eng., 88(26), 77–82 (1981b).

    [4] Seider, W. D., Lewin,  D. R., Seader, J. D., Widagdo, S., Gani, R., & Ng, M. K. (2017). Product and Process Design Principles. Wiley. Cost Accounting and Capital Cost Estimation (Chapter 16)    
    
    Examples
    --------
    :doc:`notebooks/Flash Example`
    
    """
    _units = {'Vertical vessel weight': 'lb',
              'Horizontal vessel weight': 'lb',
              'Length': 'ft',
              'Diameter': 'ft',
              'Weight': 'lb',
              'Wall thickness': 'in'}
    _has_power_utility = False
    _N_heat_utilities = 0
    
    # Bare module factor
    BM_horizontal = 3.05
    BM_vertical = 4.16
    @property
    def BM(self):
        SepType = self.SetType
        if SepType == 'Vertical':
            return self._BM_vertical
        elif SepType == 'Horizontal':
            return self._BM_horizontal
        else:
            raise AttributeError('SepType not defined')
    
    # Column material factor
    _material = 'Carbon steel'
    _F_material = 1 
    
    #: If a vacuum system is needed, it will choose one according to this preference.
    vacuum_system_preference = 'Liquid-ring pump'
    
    #: True if glycol groups are present in the mixture
    HasGlycolGroups = False
    
    #: True if amine groups are present in the mixture
    HasAmineGroups = False
    
    #: [str] 'Horizontal', 'Vertical', or 'Default'
    SepType = 'Default'
    
    #: [float] Time it takes to raise liquid to half full (min)
    HoldupTime = 15  
    
    #: [float] Time it takes to reach from normal to maximum liquied level (min)
    SurgeTime = 7.5
    
    #: [bool] True if using a mist eliminator pad
    Mist = False
    
    _bounds = {'Vertical vessel weight': (4200, 1e6),
               'Horizontal vessel weight': (1e3, 9.2e5),
               'Diameter': (3, 21),
               'Vertical vessel length': (12, 40)}

    @property
    def material(self):
        """Vessel construction material."""
        return self._material
    @material.setter
    def material(self, material):
        try: self._F_material = material_factor[material]
        except KeyError:
            dummy = str(material_factor.keys())[11:-2]
            raise ValueError(f"material must be one of the following: {dummy}")
        self._material = material  

    def __init__(self, ID='', ins=None, outs=(), *,
                 species_IDs=None, LNK=(), HNK=(),
                 V=None, T=None, Q=None, P=None,
                 y=None, x=None):
        Unit.__init__(self, ID, ins, outs)
        vap, liq = self.outs
        vap._phase = 'g'
        self._mixedstream = MixedStream(None)
        self._heat_exchanger = None
        
        #: tuple[str] IDs of species in thermodynamic equilibrium
        self.species_IDs = species_IDs
        
        #: tuple[str] Light non-keys assumed to remain as a vapor
        self.LNK = LNK
        
        #: tuple[str] Heavy non-keys assumed to remain as a liquid
        self.HNK = HNK
        
        #: Enforced molar vapor fraction
        self.V = V
        
        #: Enforced operating temperature (K)
        self.T = T
        
        #: [array_like] Molar composition of vapor (for binary mixture)
        self.y = y
        
        #: [array_like] Molar composition of liquid (for binary mixture)
        self.x = x
        
        #: Duty (kJ/hr)
        self.Q = Q
        
        #: Operating pressure (Pa)
        self.P = P
        
        
    @property
    def P(self):
        """Operating pressure (Pa)."""
        return self._P
    @P.setter
    def P(self, P):
        if P < 101325 and not self._power_utility:
            self._power_utility = PowerUtility()
        self._P = P
    
    @property
    def Q(self):
        """Enforced duty (kJ/hr)."""
        return self._Q
    @Q.setter
    def Q(self, Q):
        if Q == 0:
            self._heat_exchanger = None
        elif not self._heat_exchanger:
            self._heat_exchanger = he = HXutility(None, outs=None) 
            self._heat_utilities = he._heat_utilities
            he._ins = self._ins
            he._outs[0] = self._mixedstream
        self._Q = Q

    def _run(self):
        # Unpack
        vap, liq = self.outs
        feed = self.ins[0]

        # Vapor Liquid Equilibrium
        ms = self._mixedstream
        ms.empty()
        ms.liquid_mol[:] = feed.mol
        ms.T = feed.T
        ms.VLE(self.species_IDs, self.LNK, self.HNK, self.P,
               self.Q, self.T, self.V, self.x, self.y)

        # Set Values
        vap._mol[:] = ms.vapor_mol
        liq._mol[:] = ms.liquid_mol
        vap.T = liq.T = ms.T
        vap.P = liq.P = ms.P

    def _design(self):
        # Set horizontal or vertical vessel
        SepType = self.SepType
        if SepType == 'Default':
            if self.outs[0].massnet/self.outs[1].massnet > 0.2:
                isVertical = True
            else:
                isVertical = False
        elif SepType == 'Vertical':
            isVertical = True
        elif SepType == 'Horizontal':
            isVertical = False
        else:
            raise ValueError( f"SepType must be either 'Default', 'Horizontal', 'Vertical', not '{self.SepType}'")
        self._isVertical = isVertical

        # Run vertical or horizontal design
        if isVertical: self._vertical()
        else: self._horizontal()
        if self._heat_exchanger: self._heat_exchanger._design()
        self._Design['Material'] = self._material

    def _cost(self):
        Design = self._Design
        W = Design['Weight']
        D = Design['Diameter']
        L = Design['Length']
        CE = bst.CE
        
        # C_v: Vessel cost
        # C_pl: Platforms and ladders cost
        if self._isVertical:
            C_v = exp(7.1390 + 0.18255*ln(W) + 0.02297*ln(W)**2)
            C_pl = 410*D**0.7396*L**0.70684
        else:
            C_v = exp(5.6336 - 0.4599*ln(W) + 0.00582*ln(W)**2)
            C_pl = 2275*D**0.20294
            
        self._Cost['Flash'] = CE/567*(self._F_material*C_v+C_pl)
        if self._heat_exchanger:
            hx = self._heat_exchanger
            hx._cost()
            self._Cost.update(hx._Cost)
        self._cost_vacuum()

    def _cost_vacuum(self):
        P = self.P
        if not P or P > 101320: return 
        
        Design = self._Design
        vol = 0.02832 * np.pi * Design['Length'] * (Design['Diameter']/2)**2
        
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
                if stream._phase == 'g':
                    massflow = stream.massnet
                    volflow = stream.volnet
                else:
                    massflow = 0
                    volflow = 0
        else:
            massflow = vapor.massnet
            volflow = vapor.volnet
        
        power, cost = vacuum_system(massflow, volflow,
                                    P, vol, self.vacuum_system_preference)
        self._Cost['Liquid-ring pump'] = cost
        self._power_utility(power)

    def _design_parameters(self):
        # Retrieve run_args and properties
        vap, liq = self._outs
        rhov = vap.rho*0.06243  # VapDensity (lb/ft3)
        rhol = liq.rho*0.06243  # LLiqDensity (lb/ft3)
        P = liq.P*0.000145  # Pressure (psi)

        # SepType (str), HoldupTime (min), SurgeTime (min), Mist (bool)
        SepType = self.SepType
        Th = self.HoldupTime
        Ts = self.SurgeTime
        Mist = self.Mist

        # Calculate the volumetric flowrate
        Qv = vap.volnet * 0.0098096 # ft3/s
        Qll = liq.volnet * 0.58857 # ft3/min

        # Calculate Ut and set Uv
        K = Kvalue(P)

        # Adjust K value
        if not Mist and SepType == 'Vertical': K /= 2

        # Adjust for amine or glycol groups:
        if self.HasGlycolGroups: K *= 0.6
        elif self.HasAmineGroups: K *= 0.8

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
        if Mist:
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
        else:
            Hh = FinalValue(Hh)

        # Calculate the height from Hnll to  High liq level, Hhll
        Hs = Vs/(pi/4.0*Dvd**2)
        if Hs < 0.5:
            Hs = 0.5
        else:
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
        if Mist:
            Hv2 = 2.0 + dN/2.0
        else:
            Hv2 = 3.0 + dN/2.0

        if Hv2 < Hv: Hv = Hv2
        Hv = ceil(Hv)

        # Calculate total height, Ht
        Hme = 0.0
        if Mist:
            Hme = 1.5
        Ht = Hlll + Hh + Hs + Hlin + Hv + Hme
        Ht = FinalValue(Ht)

        # Check if LD is between 1.5 and 6.0
        while True:
            LD = Ht/D
            if (LD < 1.5): D -= 0.5
            elif (LD > 6.0): D += 0.5
            else: break

        # Calculate Vessel weight and wall thickness
        rho_M = material_density[self._material]
        VW, VWT = VesselWeightAndWallThickness(P, D, Ht, rho_M)

        # Find maximum and normal liquid level
        # Hhll = Hs + Hh + Hlll
        # Hnll = Hh + Hlll

        Design = self._Design
        self._checkbounds('Vertical vessel weight', VW, 'lb', self._bounds['Vertical vessel weight'])
        self._checkbounds('Vertical vessel length', Ht, 'ft', self._bounds['Vertical vessel length'])
        Design['SepType'] = 'Vertical'
        Design['Length'] = Ht     # ft
        Design['Diameter'] = D    # ft
        Design['Weight'] = VW     # lb
        Design['Wall thickness'] = VWT  # in
        
    def _horizontal(self):
        rhov, rhol, P, Th, Ts, Mist, Qv, Qll, Ut, Uv, Vh, Vs = self._design_parameters()

        # Initialize LD
        if P > 0 and P <= 264.7:
            LD = 1.5/250.0*(P-14.7)+1.5
        elif P > 264.7 and P <= 514.7:
            LD = 1.0/250.0*(P-14.7)+2.0
        elif P > 514.7:
            LD = 5.0

        D = (4.0*(Vh+Vs)/(0.6*pi*LD))**(1.0/3.0)
        if D <= 4.0:
            D = 4.0
        else:
            D = FinalValue(D)

        outerIter = 0
        while outerIter < 50:
            outerIter += 1
            At = pi*(D**2)/4.0 # Total area

            # Calculate Lower Liquid Area
            Hlll = round(0.5*D + 7.0)  
            Hlll = Hlll/12.0 # D is in ft but Hlll is in inches
            X = Hlll/D
            Y = HNATable(1, X)
            Alll = Y*At

            # Calculate the Vapor disengagement area, Av
            Hv = 0.2*D
            if Mist and Hv <= 2.0: Hv = 2.0
            elif Hv <= 1.0: Hv = 1.0
            else: Hv = FinalValue(Hv)
            Av = HNATable(1, Hv/D)*At
            
            # Calculate minimum length for surge and holdup
            L = (Vh + Vs)/(At - Av - Alll)
            # Calculate liquid dropout
            Phi = Hv/Uv
            # Calculate actual vapor velocity
            Uva = Qv/Av
            # Calculate minimum length for vapor disengagement
            Lmin = Uva*Phi
            Li = L
            
            innerIter = 0
            while innerIter < 50:
                if L < 0.8*Lmin: Hv += 0.5
                elif L > 1.2*Lmin:
                    if Mist and Hv <= 2.0: Hv = 2.0
                    elif not Mist and Hv <= 1.0: Hv = 1.0
                    else: Hv -= 0.5
                else: break
                Av = HNATable(1, Hv/D)*At
                Alll = HNATable(1, Hlll/D)*At
                Li = (Vh + Vs)/(At - Av - Alll)
                Phi = Hv/Uv
                Uva = Qv/Av
                Lmin = Uva*Phi                
                innerIter += 1
            
            L = Li
            LD = L/D
            # Check LD
            if LD < 1.2:
                if D <= 4.0: break
                else: D -= 0.5

            if LD > 7.2:
                D += 0.5
            else: break

        # Recalculate LD so it lies between 1.5 - 6.0
        while True:
            LD = L / D
            if (LD < 1.5) and D <= 4.0: L += 0.5
            elif LD < 1.5: D -= 0.5
            elif (LD > 6.0): D += 0.5
            else: break

        # Calculate vessel weight and wall thickness
        rho_M = material_density[self._material]
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
        
        Design = self._Design
        self._checkbounds('Horizontal vessel weight', VW, 'lb', self._bounds['Horizontal vessel weight'])
        Design['SepType'] = 'Horizontal'
        Design['Length'] = L  # ft
        Design['Diameter'] = D  # ft
        Design['Weight'] = VW  # lb
        Design['Wall thickness'] = VWT  # in    
        
    def _end_decorated_cost_(self):
        if self._heat_utilities: self._heat_utilities[0](self._Hnet, self.T)
    

# %% Special

class SplitFlash(Flash):
    line = 'Flash' 
    
    def __init__(self, ID='', ins=None, outs=(), *, split,
                 order=None, T=None, P=None, Q=None):
        Splitter.__init__(self, ID, ins, outs, split=split, order=order)
        self._mixedstream = MixedStream(None)
        self._heat_exchanger = None
        self.T = T #: Operating temperature (K)
        self.P = P #: Operating pressure (Pa)
        self.Q = Q #: Duty (kJ/hr)
    
    split = Splitter.split
    V = None
    
    def _run(self):
        top, bot = self.outs
        net_mol = self._ins[0].mol
        top._mol[:] = net_mol * self._split
        bot._mol[:] = net_mol - top._mol
        top.phase = 'g'
        bot.T = top.T = self.T
        bot.P = top.P = self.P

    def _design(self):
        if self._heat_exchanger:
            self._heat_exchanger.outs[0] = ms = self._mixedstream
            ms = MixedStream.sum(ms, self.outs)
        super()._design()


# class PartitionFlash(Flash):
#     """Create a PartitionFlash Unit that approximates outputs based on molar partition coefficients."""
    
#     _kwargs = {'species_IDs': None,  # Partition species
#                'Ks': None,           # Partition coefficients
#                'LNK': None,          # light non-keys
#                'HNK': None,          # heavy non-keys
#                'P': None,            # Operating Pressure
#                'T': None}            # Operating temperature
    
#     def _init(self):
#         top, bot = self.outs
#         species_IDs, Ks, LNK, HNK, P, T = self._kwargs.values()
#         index = top._species._IDs.index
#         if P < 101325 and not self._power_utility:
#             self._power_utility = PowerUtility()
#         self._args = ([index(i) for i in species_IDs],
#                       [index(i) for i in LNK],
#                       [index(i) for i in HNK],
#                        np.asarray(Ks), P, T)
    
#     def _run(self):
#         feed = self.ins[0]
#         ph1, ph2 = self.outs
#         pos, LNK, HNK, Ks, P, T = self._args
#         mol_in = feed.mol[pos]
#         ph1.T = ph2.T = (T if T else feed.T)
#         ph1.P = ph2.P = (P if P else feed.P)
#         ph1.mol[LNK] = feed.mol[LNK]
#         ph1.mol[HNK] = feed.mol[HNK]
        
#         # Get zs and Ks to solve rashford rice equation
#         molnet = sum(mol_in)
#         zs = mol_in/molnet
#         lenzs = len(zs)
#         if hasattr(self, '_V'):
#             self._V = V = newton(V_error, self._V, args=(zs, Ks))
#         elif lenzs > 3:
#             self._V = V = brentq(V_error, 0, 1, (zs, Ks))
#         elif lenzs == 2:
#             self._V = V = V_2N(V, zs, Ks)
#         elif lenzs == 3:
#             self._V = V = V_3N(V, zs, Ks)
            
#         x = zs/(1 + V*(Ks-1))
#         y = x*Ks
#         ph1.mol[pos] = molnet*V*y
#         ph2.mol[pos] = mol_in - ph1.mol[pos]
    
#     _design = SplitFlash._design
    

class RatioFlash(Flash):
    _N_heat_utilities = 1

    def __init__(self, ID='', ins=None, outs=(), *,
                 Kspecies, Ks, top_solvents=(), top_split=(),
                 bot_solvents=(), bot_split=()):
        Unit.__init__(self, ID, ins, outs)
        self.Kspecies = Kspecies
        self.Ks = Ks
        self.top_solvents = top_solvents
        self.top_split = top_split
        self.bot_solvents = bot_solvents
        self.bot_split = bot_split

    def _run(self):
        feed = self.ins[0]
        top, bot = self.outs
        indices = feed.indices
        Kindex = indices(self.Kspecies)
        top_index = indices(self.top_solvents)
        bot_index = indices(self.bot_solvents)
        top_mol = top.mol; bot_mol = bot.mol; feed_mol = feed.mol
        top_mol[top_index] = feed_mol[top_index]*self.top_split
        bot_mol[top_index] = feed_mol[top_index]-top_mol[top_index]
        bot_mol[bot_index] = feed_mol[bot_index]*self.bot_split
        top_mol[bot_index] = feed_mol[bot_index]-bot_mol[bot_index]
        topnet = top_mol[top_index].sum()
        botnet = bot_mol[bot_index].sum()
        molnet = topnet+botnet
        top_mol[Kindex] = self.Ks * topnet * feed_mol[Kindex] / molnet  # solvent * mol ratio
        bot_mol[Kindex] = feed_mol[Kindex] - top_mol[Kindex]
        top.T, top.P = feed.T, feed.P
        bot.T, bot.P = feed.T, feed.P


# %% Single Component

class Evaporator_PQ(Unit):

    @property
    def P(self):
        return self._P
    @P.setter
    def P(self, P):
        vap = self.outs[0]
        liq = self.outs[1]
        liq.T = vap.T = getattr(vap.species, self.component).Tsat(P)
        liq.P = vap.P = P
    
    def __init__(self, ID='', ins=None, outs=(), *,
                 component='Water', Q=0, P=101325):
        super().__init__(ID, ins, outs)
        self.component = component
        self.Q = Q
        self.P = P
        vapor, liquid = self.outs[:2]
        Tf = getattr(vapor._species, component).Tsat(P)

        # Set-up streams for energy balance
        vapor.empty()
        species = vapor._species
        index = species.index(component)
        vapor.mol[index] = 1
        vapor._phase = 'g'
        vapor.T = vapor.T = Tf

        liquid.empty()
        liquid.mol[index] = 1

        no_ph_ch = Stream(None)
        self._vapor_H = vapor.H
        self._liquid_H = liquid.H
        self._no_ph_ch = no_ph_ch
        self._index = species.index(component)
        self._V = 0.5

    def _run(self):
        feed = self.ins[0]
        vapor, liquid = self.outs[:2]

        # Set-up streams for energy balance
        vapor_H = self._vapor_H
        liquid_H = self._liquid_H
        no_ph_ch = self._no_ph_ch
        no_ph_ch._mol[:] = feed.mol
        index = self._index
        no_ph_ch.mol[index] = 0
        Q = self.Q
        no_ph_ch_H = no_ph_ch.H
        feed_H = feed.H
        # Optional if Q also comes from condensing a side stream
        if len(self.ins) == 2:
            boiler_liq = self.outs[2]
            boiler_vap = self.ins[1]
            boiler_liq.copylike(boiler_vap)
            boiler_liq.phase = 'l'
            Q = Q + boiler_vap.H - boiler_liq.H

        # Energy balance to find vapor fraction
        f = feed.mol[index]
        H_actual = lambda v:  vapor_H*(v*f) + liquid_H*((1-v)*f)

        # Final vapor fraction
        V = self._V
        H = feed_H + Q - no_ph_ch_H
        V = IQ_interpolation(H_actual, 0, 1, liquid_H, vapor_H,
                             V, H, xtol=1e-4, ytol=1e-3)
        liquid._mol[:] = no_ph_ch._mol
        if V < 0:
            V = 0
        elif V > 1:
            V = 1
        vapor.mol[index] = V*feed.mol[index]
        liquid.mol[index] = (1-V)*feed.mol[index]
        self._Q = Q
        self._V = V


class Evaporator_PV(Unit):
    _N_heat_utilities = 1
    _kwargs = {'component': 'Water',
               'V': 0.5,
               'P': 101325}

    @property
    def P(self):
        return self._P
    @P.setter
    def P(self, P):
        vap, liq = self.outs
        liq.T = vap.T = getattr(vap.species, self.component).Tsat(P)
        liq.P = vap.P = P
    
    @property
    def component(self):
        return self._component

    def __init__(self, ID='', ins=None, outs=(), *,
                 component='Water', V=0.5, P=101325):
        super().__init__(ID, ins, outs)
        vap, liq = self.outs
        self._index = vap.species.index(component)
        self._component = component
        self.V = V
        self.P = P
        vap._phase = 'g'
        liq.phase = 'l'

    def _run(self):
        feed = self.ins[0]
        vapor, liquid = self.outs
        index = self._index
        vapor._mol[index] = self.V*feed.mol[index]
        liquid._mol[:] = feed.mol
        liquid._mol[index] = (1-self.V)*feed.mol[index]

        