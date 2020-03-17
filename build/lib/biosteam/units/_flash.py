# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 16:21:56 2018

@author: yoelr
"""
from .. import Unit, PowerUtility
from thermosteam import MultiStream
from math import pi, ceil
import numpy as np
from .design_tools import (compute_vacuum_system_power_and_cost,
                           HNATable, ceil_half_step,
                           compute_vertical_vessel_purchase_cost,
                           compute_horizontal_vessel_purchase_cost,
                           compute_vessel_weight_and_wall_thickness,
                           compute_Stokes_law_York_Demister_K_value,
                           pressure_vessel_material_factors,
                           material_densities_lb_per_ft3,)
from ._splitter import Splitter
from ._hx import HX, HXutility
from .._graphics import vertical_vessel_graphics
from ..utils import bounds_warning
from ._CAS import H2O_CAS

exp = np.exp
ln = np.log

# USDA Biodiesel: 
#V = np.pi*D**2/4*Ht*0.3048**3
#if not (0.1 < V < 70):
#    raise DesignError(f"Volume is out of bounds for costing")
#lambda V, CE: CE*13*V**0.62 # V (m3)
__all__ = ('Flash', 'SplitFlash', 'RatioFlash')


# %% Flash

class Flash(Unit):
    """
    Create an equlibrium based flash drum with the option of having light
    non-keys and heavy non-keys completly separate into their respective
    phases. Design procedure is based on heuristics by Wayne D. Monnery & 
    William Y. Svrcek [1]_. Purchase costs are based on correlations by
    Mulet et al. [2]_ [3]_ as compiled by Warren et. al. [4]_.

    Parameters
    ----------
    ins : stream
        Inlet fluid.
    outs : stream sequence
        * [0] Vapor product
        * [1] Liquid product
    P=None : float
        Operating pressure [Pa].
    Q=None : float
        Duty [kJ/hr].
    T=None : float
        Operating temperature [K].
    V=None : float
        Molar vapor fraction.
    x=None : float
        Molar composition of liquid (for binary mixtures).
    y=None : float
        Molar composition of vapor (for binary mixtures).
    vessel_material : str, optional
        Vessel construction material. Defaults to 'Carbon steel'.
    vacuum_system_preference : 'Liquid-ring pump', 'Steam-jet ejector', or 'Dry-vacuum pump'
        If a vacuum system is needed, it will choose one according to this
        preference. Defaults to 'Liquid-ring pump'.
    has_glycol_groups=False : bool
        True if glycol groups are present in the mixture.
    has_amine_groups=False : bool
        True if amine groups are present in the mixture.
    vessel_type='Default' : 'Horizontal', 'Vertical', or 'Default'
        Vessel separation type. If 'Default', the vessel type will be chosen
        according to heuristics.
    holdup_time=15.0 : float
        Time it takes to raise liquid to half full [min].
    surge_time=7.5 : float
        Time it takes to reach from normal to maximum liquied level [min].
    has_mist_eliminator : bool
        True if using a mist eliminator pad.

    Notes
    -----
    You may only specify two of the following parameters: P, Q, T, V, x, and y.
    Additionally, If x or y is specified, the other parameter must be either
    P or T (e.g., x and V is invalid).

    Examples
    --------
    >>> from biosteam.units import Flash
    >>> from biosteam import Stream, settings
    >>> settings.set_thermo(['Water', 'Glycerol'])
    >>> feed = Stream('feed', Glycerol=300, Water=1000)
    >>> bp = feed.bubble_point_at_P() # Feed at bubble point T
    >>> feed.T = bp.T
    >>> F1 = Flash('F1',
    ...            ins=feed,
    ...            outs=('vapor', 'crude_glycerin'),
    ...            P=101325, # Pa
    ...            T=410.15) # K
    >>> F1.simulate()
    >>> F1.show(T='degC', P='atm')
    Flash: F1
    ins...
    [0] feed
        phase: 'l', T: 100.7 degC, P: 1 atm
        flow (kmol/hr): Water     1e+03
                        Glycerol  300
    outs...
    [0] vapor
        phase: 'g', T: 137 degC, P: 1 atm
        flow (kmol/hr): Water     958
                        Glycerol  2.32
    [1] crude_glycerin
        phase: 'l', T: 137 degC, P: 1 atm
        flow (kmol/hr): Water     42.4
                        Glycerol  298
    >>> F1.results()
    Flash                                   Units            F1
    Medium pressure steam Duty              kJ/hr      5.06e+07
                          Flow            kmol/hr       1.4e+03
                          Cost             USD/hr           385
    Design                Vessel type                  Vertical
                          Length               ft          23.5
                          Diameter             ft           6.5
                          Weight               lb      8.44e+03
                          Wall thickness       in         0.375
                          Material                 Carbon steel
    Purchase cost         Flash               USD      5.82e+04
                          Heat exchanger      USD      4.04e+04
    Total purchase cost                       USD      9.86e+04
    Utility cost                           USD/hr           385


    References
    ----------
    .. [1] "Design Two-Phase Separators Within the Right Limits", Chemical
        Engineering Progress Oct, 1993.

    .. [2] Mulet, A., A. B. Corripio, and L. B. Evans, “Estimate Costs of
        Pressure Vessels via Correlations,” Chem. Eng., 88(20), 145–150 (1981a).

    .. [3] Mulet, A., A.B. Corripio, and L.B.Evans, “Estimate Costs of
        Distillation and Absorption Towers via Correlations,” Chem. Eng., 88(26), 77–82 (1981b).

    .. [4] Seider, W. D., Lewin,  D. R., Seader, J. D., Widagdo, S., Gani,
        R., & Ng, M. K. (2017). Product and Process Design Principles. Wiley.
        Cost Accounting and Capital Cost Estimation (Chapter 16)
    
    """
    _units = {'Vertical vessel weight': 'lb',
              'Horizontal vessel weight': 'lb',
              'Length': 'ft',
              'Diameter': 'ft',
              'Weight': 'lb',
              'Wall thickness': 'in'}
    _graphics = vertical_vessel_graphics 
    _N_outs = 2
    _N_heat_utilities = 0
    
    # Bare module factor
    BM_horizontal = 3.05
    BM_vertical = 4.16
    @property
    def BM(self):
        vessel_type = self.vessel_type
        if vessel_type == 'Vertical':
            return self.BM_vertical
        elif vessel_type == 'Horizontal':
            return self.BM_horizontal
        elif vessel_type == 'Default':
            return self.BM_vertical if self._isVertical else self.BM_horizontal 
        else:
            raise AttributeError('vessel_type not defined')
    
    _bounds = {'Vertical vessel weight': (4200, 1e6),
               'Horizontal vessel weight': (1e3, 9.2e5),
               'Diameter': (3, 21),
               'Vertical vessel length': (12, 40)}

    @property
    def vessel_material(self):
        """Vessel construction material."""
        return self._vessel_material
    @vessel_material.setter
    def vessel_material(self, material):
        try: self._F_M = pressure_vessel_material_factors[material]
        except KeyError:
            raise ValueError(f"no material factor available for '{material}'; "
                              "only the following materials are available: "
                             f"{', '.join(pressure_vessel_material_factors)}")
        self._vessel_material = material  

    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 V=None, T=None, Q=None, P=None, y=None, x=None,
                 vessel_material='Carbon steel',
                 vacuum_system_preference='Liquid-ring pump',
                 has_glycol_groups=False,
                 has_amine_groups=False,
                 vessel_type='Default',
                 holdup_time=15,
                 surge_time=7.5,
                 has_mist_eliminator=False):
        Unit.__init__(self, ID, ins, outs, thermo)
        self._multistream = MultiStream(None)
        self.heat_exchanger = None
        
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
        
        #: [str] Vessel construction material
        self.vessel_material = vessel_material

        #: [str] If a vacuum system is needed, it will choose one according to this preference.
        self.vacuum_system_preference = vacuum_system_preference
        
        #: [bool] True if glycol groups are present in the mixture
        self.has_glycol_groups = has_glycol_groups
        
        #: [bool] True if amine groups are present in the mixture
        self.has_amine_groups = has_amine_groups
        
        #: [str] 'Horizontal', 'Vertical', or 'Default'
        self.vessel_type = vessel_type
        
        #: [float] Time it takes to raise liquid to half full (min)
        self.holdup_time = holdup_time
        
        #: [float] Time it takes to reach from normal to maximum liquied level (min)
        self.surge_time = surge_time
        
        #: [bool] True if using a mist eliminator pad
        self.has_mist_eliminator = has_mist_eliminator
        
    @property
    def P(self):
        """Operating pressure (Pa)."""
        return self._P
    @P.setter
    def P(self, P):
        if P < 101325 and not self.power_utility:
            self.power_utility = PowerUtility()
        self._P = P
    
    @property
    def Q(self):
        """Enforced duty (kJ/hr)."""
        return self._Q
    @Q.setter
    def Q(self, Q):
        if Q == 0:
            self.heat_exchanger = None
        elif not self.heat_exchanger:
            self.heat_exchanger = he = HXutility(None, outs=None) 
            self.heat_utilities = he.heat_utilities
            he._ins = self._ins
            he._outs[0] = self._multistream
        self._Q = Q

    def _run(self):
        # Unpack
        vap, liq = self.outs
        feed = self.ins[0]

        # Vapor Liquid Equilibrium
        ms = self._multistream
        ms.empty()
        ms.imol['l'] = feed.mol
        ms.T = feed.T
        Q = self.Q
        H = feed.H + Q if Q is not None else None
        ms.vle(P=self.P, H=H, T=self.T, V=self.V, x=self.x, y=self.y)

        # Set Values
        vap.phase = 'g'
        liq.phase = 'l'
        vap.mol[:] = ms.imol['g']
        liq.mol[:] = ms.imol['l']
        vap.T = liq.T = ms.T
        vap.P = liq.P = ms.P

    def _design(self):
        # Set horizontal or vertical vessel
        vessel_type = self.vessel_type
        if vessel_type == 'Default':
            vap, liq = self.outs
            isVertical = vap.F_mass/liq.F_mass > 0.2
        elif vessel_type == 'Vertical':
            isVertical = True
        elif vessel_type == 'Horizontal':
            isVertical = False
        else:
            raise ValueError( f"vessel_type must be either 'Default', 'Horizontal', 'Vertical', not '{self.vessel_type}'")
        self._isVertical = isVertical

        # Run vertical or horizontal design
        if isVertical: self._design_vertical_vessel()
        else: self._design_horizontal_vessel()
        if self.heat_exchanger: self.heat_exchanger._design()
        self.design_results['Material'] = self._vessel_material

    def _cost(self):
        Design = self.design_results
        W = Design['Weight']
        D = Design['Diameter']
        L = Design['Length']
        F_M = self._F_M
        if self._isVertical:
            Cp = compute_vertical_vessel_purchase_cost(W, D, L, F_M)
        else:
            Cp = compute_horizontal_vessel_purchase_cost(W, D, F_M)
        self.purchase_costs['Flash'] = Cp
        if self.heat_exchanger:
            hx = self.heat_exchanger
            hx._cost()
            self.purchase_costs.update(hx.purchase_costs)
        self._cost_vacuum()

    def _cost_vacuum(self):
        P = self.P
        if not P or P > 101320: return 
        
        Design = self.design_results
        volume = 0.02832 * np.pi * Design['Length'] * (Design['Diameter']/2)**2
        
        # If vapor is condensed, assume vacuum system is after condenser
        vapor = self.outs[0]
        hx = vapor.sink
        if isinstance(hx, HX):
            index = hx.ins.index(vapor)
            stream = hx.outs[index]
            if isinstance(stream, MultiStream):
                vapor = stream['g']
                F_mass = vapor.F_mass
                F_vol = vapor.F_vol
            else:
                if stream.phase == 'g':
                    F_mass = stream.F_mass
                    F_vol = stream.F_vol
                else:
                    F_mass = 0
                    F_vol = 0
        else:
            F_mass = vapor.F_mass
            F_vol = vapor.F_vol
        
        power, cost = compute_vacuum_system_power_and_cost(
                          F_mass, F_vol, P, volume, self.vacuum_system_preference)
        self.purchase_costs['Liquid-ring pump'] = cost
        self.power_utility(power)

    def _design_parameters(self):
        # Retrieve run_args and properties
        vap, liq = self._outs
        rhov = vap.get_property('rho', 'lb/ft3')
        rhol = liq.get_property('rho', 'lb/ft3')
        P = liq.get_property('P', 'psi')  # Pressure (psi)

        vessel_type = self.vessel_type
        Th = self.holdup_time
        Ts = self.surge_time
        has_mist_eliminator = self.has_mist_eliminator

        # Calculate the volumetric flowrate
        Qv = vap.get_total_flow('ft^3 / s')
        Qll = liq.get_total_flow('ft^3 / min')

        # Calculate Ut and set Uv
        K = compute_Stokes_law_York_Demister_K_value(P)

        # Adjust K value
        if not has_mist_eliminator and vessel_type == 'Vertical': K /= 2

        # Adjust for amine or glycol groups:
        if self.has_glycol_groups: K *= 0.6
        elif self.has_amine_groups: K *= 0.8

        Ut = K*((rhol - rhov) / rhov)**0.5
        Uv = 0.75*Ut

        # Calculate Holdup and Surge volume
        Vh = Th*Qll
        Vs = Ts*Qll
        return rhov, rhol, P, Th, Ts, has_mist_eliminator, Qv, Qll, Ut, Uv, Vh, Vs

    def _design_vertical_vessel(self):
        rhov, rhol, P, Th, Ts, has_mist_eliminator, Qv, Qll, Ut, Uv, Vh, Vs = self._design_parameters()

        # Calculate internal diameter, Dvd
        Dvd = (4.0*Qv/(pi*Uv))**0.5
        if has_mist_eliminator:
            D = ceil_half_step(Dvd + 0.4)
        else:
            D = ceil_half_step(Dvd)

        # Obtaining low liquid level height, Hlll
        Hlll = 0.5
        if P < 300:
            Hlll = 1.25

        # Calculate the height from Hlll to Normal liq level, Hnll
        Hh = Vh/(pi/4.0*Dvd**2)
        if Hh < 1.0: Hh = 1.0

        # Calculate the height from Hnll to  High liq level, Hhll
        Hs = Vs/(pi/4.0*Dvd**2)
        if Hs < 0.5: Hs = 0.5

        # Calculate dN
        Qm = Qll + Qv
        lamda = Qll/Qm
        rhoM = rhol*lamda + rhov*(1-lamda)
        dN = (4*Qm/(pi*60.0/(rhoM**0.5)))**0.5
        dN = ceil_half_step(dN)

        # Calculate Hlin, assume with inlet diverter
        Hlin = 1.0 + dN

        # Calculate the vapor disengagement height
        Hv = 0.5*Dvd
        Hv2 = (2.0 if has_mist_eliminator else 3.0) + dN/2.0
        if Hv2 < Hv: Hv = Hv2
        Hv = Hv

        # Calculate total height, Ht
        Hme = 1.5 if has_mist_eliminator else 0.0
        Ht = Hlll + Hh + Hs + Hlin + Hv + Hme
        Ht = ceil_half_step(Ht)

        # Calculate Vessel weight and wall thickness
        rho_M = material_densities_lb_per_ft3[self._vessel_material]
        VW, VWT = compute_vessel_weight_and_wall_thickness(P, D, Ht, rho_M)

        # Find maximum and normal liquid level
        # Hhll = Hs + Hh + Hlll
        # Hnll = Hh + Hlll

        Design = self.design_results
        bounds_warning(self, 'Vertical vessel weight', VW, 'lb',
                       self._bounds['Vertical vessel weight'],
                       'cost')
        bounds_warning(self, 'Vertical vessel length', Ht, 'ft',
                       self._bounds['Vertical vessel length'],
                       'cost')
        Design['Vessel type'] = 'Vertical'
        Design['Length'] = Ht     # ft
        Design['Diameter'] = D    # ft
        Design['Weight'] = VW     # lb
        Design['Wall thickness'] = VWT  # in
        
    def _design_horizontal_vessel(self):
        rhov, rhol, P, Th, Ts, has_mist_eliminator, Qv, Qll, Ut, Uv, Vh, Vs = self._design_parameters()

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
            D = ceil_half_step(D)

        for outerIter in range(50):
            At = pi*(D**2)/4.0 # Total area

            # Calculate Lower Liquid Area
            Hlll = round(0.5*D + 7.0)  
            Hlll = Hlll/12.0 # D is in ft but Hlll is in inches
            X = Hlll/D
            Y = HNATable(1, X)
            Alll = Y*At

            # Calculate the Vapor disengagement area, Av
            Hv = 0.2*D
            if has_mist_eliminator and Hv <= 2.0: Hv = 2.0
            elif Hv <= 1.0: Hv = 1.0
            else: Hv = ceil_half_step(Hv)
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
            
            for innerIter in range(50):
                if L < 0.8*Lmin: Hv += 0.5
                elif L > 1.2*Lmin:
                    if has_mist_eliminator and Hv <= 2.0: Hv = 2.0
                    elif not has_mist_eliminator and Hv <= 1.0: Hv = 1.0
                    else: Hv -= 0.5
                else: break
                Av = HNATable(1, Hv/D)*At
                Alll = HNATable(1, Hlll/D)*At
                Li = (Vh + Vs)/(At - Av - Alll)
                Phi = Hv/Uv
                Uva = Qv/Av
                Lmin = Uva*Phi
            
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
        rho_M = material_densities_lb_per_ft3[self._vessel_material]
        VW, VWT = compute_vessel_weight_and_wall_thickness(P, D, L, rho_M)

        # # To check minimum Hv value
        # if int(has_mist_eliminator) == 1 and Hv <= 2.0:
        #     Hv = 2.0
        # if int(has_mist_eliminator) == 0 and Hv <= 1.0:
        #     Hv = 1.0

        # Calculate normal liquid level and High liquid level
        # Hhll = D - Hv
        # if (Hhll < 0.0):
        #     Hhll = 0.0
        # Anll = Alll + Vh/L
        # X = Anll/At
        # Y = HNATable(2, X)
        # Hnll = Y*D
        
        Design = self.design_results
        bounds_warning(self, 'Horizontal vessel weight', VW, 'lb',
                       self._bounds['Horizontal vessel weight'], 'cost')
        Design['Vessel type'] = 'Horizontal'
        Design['Length'] = L  # ft
        Design['Diameter'] = D  # ft
        Design['Weight'] = VW  # lb
        Design['Wall thickness'] = VWT  # in    
        
    def _end_decorated_cost_(self):
        if self.heat_utilities: self.heat_utilities[0](self.Hnet, self.T)
    

# %% Special

class SplitFlash(Flash):
    line = 'Flash' 
    
    def __init__(self, ID='', ins=None, outs=(), *, split,
                 order=None, T=None, P=None, Q=None,
                 vessel_material='Carbon steel',
                 vacuum_system_preference='Liquid-ring pump',
                 has_glycol_groups=False,
                 has_amine_groups=False,
                 vessel_type='Default',
                 holdup_time=15,
                 surge_time=7.5,
                 has_mist_eliminator=False):
        Splitter.__init__(self, ID, ins, outs, split=split, order=order)
        self._multistream = MultiStream(None)
        
        #: [HXutility] Heat exchanger if needed.
        self.heat_exchanger = None
        self.T = T #: Operating temperature (K)
        self.P = P #: Operating pressure (Pa)
        self.Q = Q #: Duty (kJ/hr)
        
        #: [str] Vessel construction material
        self.vessel_material = vessel_material

        #: [str] If a vacuum system is needed, it will choose one according to this preference.
        self.vacuum_system_preference = vacuum_system_preference
        
        #: [bool] True if glycol groups are present in the mixture
        self.has_glycol_groups = has_glycol_groups
        
        #: [bool] True if amine groups are present in the mixture
        self.has_amine_groups = has_amine_groups
        
        #: [str] 'Horizontal', 'Vertical', or 'Default'
        self.vessel_type = vessel_type
        
        #: [float] Time it takes to raise liquid to half full (min)
        self.holdup_time = holdup_time
        
        #: [float] Time it takes to reach from normal to maximum liquied level (min)
        self.surge_time = surge_time
        
        #: [bool] True if using a mist eliminator pad
        self.has_mist_eliminator = has_mist_eliminator
    
    split = Splitter.split
    V = None
    
    def _run(self):
        top, bot = self.outs
        feed, = self.ins
        feed_mol = feed.mol
        top.mol[:] = top_mol = feed_mol * self.split
        bot.mol[:] = feed_mol - top_mol
        top.phase = 'g'
        bot.phase = 'l'
        bot.T = top.T = self.T
        bot.P = top.P = self.P

    def _design(self):
        if self.heat_exchanger:
            self.heat_exchanger.outs[0] = ms = self._multistream
            ms.mix_from(self.outs)
        super()._design()
    

class RatioFlash(Flash):
    _N_heat_utilities = 1

    def __init__(self, ID='', ins=None, outs=(), *,
                 K_chemicals, Ks, top_solvents=(), top_split=(),
                 bot_solvents=(), bot_split=()):
        Unit.__init__(self, ID, ins, outs)
        self.K_chemicals = K_chemicals
        self.Ks = Ks
        self.top_solvents = top_solvents
        self.top_split = top_split
        self.bot_solvents = bot_solvents
        self.bot_split = bot_split

    def _run(self):
        feed = self.ins[0]
        top, bot = self.outs
        indices = self.chemicals.get_index
        K_index = indices(self.K_chemicals)
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
        top_mol[K_index] = self.Ks * topnet * feed_mol[K_index] / molnet  # solvent * mol ratio
        bot_mol[K_index] = feed_mol[K_index] - top_mol[K_index]
        top.T, top.P = feed.T, feed.P
        bot.T, bot.P = feed.T, feed.P


# %% Single Component

class Evaporator_PQ(Unit):
    _N_ins = 2
    _N_outs = 3
    @property
    def P(self):
        return self._P
    @P.setter
    def P(self, P):
        water = getattr(self.chemicals, H2O_CAS)
        self._T = T = water.Tsat(P)
        self._Hvap = water.Hvap(T)
        self._P = P
    @property
    def T(self):
        return self._T
    @T.setter
    def T(self, T):
        water = getattr(self.chemicals, H2O_CAS)
        self._P = water.Psat(T)
        self._Hvap = water.Hvap(T)
        self._T = T
    @property
    def V(self):
        return self._V
    
    def __init__(self, ID='', ins=None, outs=(), *, Q=0, P=101325):
        super().__init__(ID, ins, outs)
        self.Q = Q
        self.P = P
        self._V = None
    
    def _run(self):
        feed, utility_vapor = self.ins
        vapor, liquid, utility_liquid = self.outs
        
        # Optional if Q also comes from condensing a side stream
        Q = self.Q
        if utility_liquid:
            utility_liquid.copy_like(utility_vapor)
            utility_liquid.phase = 'l'
            Q += utility_vapor.Hvap
        feed_H = feed.H
    
        # Set exit conditions
        vapor.T = liquid.T = self.T
        vapor.P = liquid.P = self.P
        liquid.phase = 'l'
        vapor.phase = 'g'
        liquid.copy_flow(feed, IDs=H2O_CAS, exclude=True)
        
        # Energy balance to find vapor fraction
        f = feed.imol[H2O_CAS]
        H = feed_H + Q - liquid.H
        V = H/(f * self._Hvap)
        if V < 0:
            V = 0
        elif V > 1:
            V = 1
        evaporated = f * V
        vapor.imol[H2O_CAS] = evaporated
        liquid.imol[H2O_CAS] = (1-V)*f
        self._Q = Q
        self._V = V


class Evaporator_PV(Unit):
    _N_heat_utilities = 1

    @property
    def P(self):
        return self._P
    @P.setter
    def P(self, P):
        water = getattr(self.chemicals, H2O_CAS)
        self._T = water.Tsat(P)
        self._P = P
    @property
    def T(self):
        return self._T
    @T.setter
    def T(self, T):
        water = getattr(self.chemicals, H2O_CAS)
        self._P = water.Psat(T)
        self._T = T
    
    def __init__(self, ID='', ins=None, outs=(), *, V=0.5, P=101325):
        super().__init__(ID, ins, outs)
        self.V = V
        self.P = P

    def _run(self):
        feed = self.ins[0]
        vapor, liquid = self.outs
        vapor.T = liquid.T = self.T
        H2O_index = self.chemicals.index(H2O_CAS)
        water_mol = feed.mol[H2O_index]
        vapor.mol[H2O_index] = self.V * water_mol
        liquid_mol = liquid.mol
        liquid_mol[:] = feed.mol
        liquid_mol[H2O_index] = (1-self.V) * water_mol

        