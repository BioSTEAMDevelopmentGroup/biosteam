#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 19:33:20 2018

@author: yoelr
"""
import numpy as np
from thermosteam import MultiStream, Stream
from .. import Unit
from ..utils import approx2step
from scipy.optimize import brentq
from ._hx import HXutility
import matplotlib.pyplot as plt
import biosteam as bst
array = np.array

# %% Equations

# Minimum thickness
x = array((4, 6, 8, 10, 12)) 
y = array((1/4, 5/16, 3/8, 7/16, 1/2))
ts_min_p = np.polyfit(x,y,1)

# %% Dictionary of factors

# Tray Type
F_TTdict = {'Sieve': 1,
            'Valve': 1.18,
            'Bubble cap': 1.87}

# Tray Materials (inner diameter, Di, in ft)

def compute_carbon_steel_material_factor(Di):
    return 1

def compute_stainless_steel_304_material_factor(Di):
    return 1.189 + 0.058*Di

def compute_stainless_steel_316_material_factor(Di):
    return 1.401 + 0.073*Di

def compute_carpenter_20CB3_material_factor(Di):
    return 1.525 + 0.079*Di

def compute_monel_material_factor(Di):
    return 2.306 + 0.112*Di

F_TMdict = {'Carbon steel': compute_carbon_steel_material_factor,
            'Stainless steel 304': compute_stainless_steel_304_material_factor,
            'Stainless steel 316': compute_stainless_steel_316_material_factor,
            'Carpenter 20CB-3': compute_carpenter_20CB3_material_factor,
            'Monel': compute_monel_material_factor}

# Column Material
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

# Material density (lb∕in^3)
rho_Mdict = {'Carbon steel': 0.284 ,
             'Low-alloy steel': None,
             'Stainless steel 304': 0.289,
             'Stainless steel 316': 0.289,
             'Carpenter 20CB-3': None,
             'Nickel-200': None,
             'Monel-400': None,
             'Inconel-600': None,
             'Incoloy-825': None,
             'Titanium': None}

# %% Distillation

class Dist(Unit, isabstract=True):
    """Abstract class for a column."""
    # Bare module factor
    BM = 4.3
    
    # Column material factor
    _F_Mstr = 'Carbon steel'
    _F_M = 1
    
    # Tray type factor
    _F_TTstr = 'Sieve'
    _F_TT = 1
    
    # Tray material factor function
    _F_TMstr = 'Carbon steel'
    _F_TM = staticmethod(F_TMdict['Carbon steel'])
    
    # [float] Tray spacing (225-600 mm)
    _TS = 450 
    
    #: [float] Enforced user defined stage efficiency.        
    _stage_efficiency = None
    
    # [float] Ratio of actual velocity to maximum velocity allowable before flooding.
    _f = 0.8 
    
    # [float] Foaming factor (0 < F_F < 1)
    _F_F = 1
    
    # [float] Ratio of open area, A_h, to active area, A_a
    _A_ha = 0.1
    
    # [float] Enforced ratio of downcomer area to net (total) area. If None, estimate ratio based on Oliver's estimation [1].
    _A_dn = None
    
    # [dict] Bounds for results
    _bounds = {'Diameter': (3., 24.),
               'Height': (27., 170.),
               'Weight': (9000., 2.5e6)}
    
    @property
    def LHK(self):
        return self._LHK
    @LHK.setter
    def LHK(self, LHK):
        # Set light non-key and heavy non-key indices
        self._LHK = LHK = tuple(LHK)
        chemicals = self.chemicals
        LHK_chemicals = LK_chemical, HK_chemical = self.chemicals.retrieve(LHK)
        Tb_light = LK_chemical.Tb
        Tb_heavy = HK_chemical.Tb
        LNK = []
        HNK = []
        if Tb_light > Tb_heavy:
            raise ValueError(f"light key must be lighter than heavy key")
        for chemical in chemicals:
            Tb = chemical.Tb
            if not Tb:
                HNK.append(chemical.ID)
            elif Tb < Tb_light:
                LNK.append(chemical.ID)
            elif Tb > Tb_heavy:
                HNK.append(chemical.ID)
            elif chemical not in LHK_chemicals:
                raise ValueError(f"intermediate volatile specie, '{chemical}', between light and heavy key, ['{LK_chemical}', '{HK_chemical}']")
        self._LNK = tuple(LNK)
        self._HNK = tuple(HNK)
    
    @property
    def y_top(self):
        return self._y_top
    @y_top.setter
    def y_top(self, y_top):
        self._y_top = y_top
        self._y = array([y_top, 1-y_top])
    
    @property
    def x_bot(self):
        return self._x_bot
    @x_bot.setter
    def x_bot(self, x_bot):
        self._x_bot = x_bot
        self._x = array([x_bot, 1-x_bot])
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 P=101325, *, LHK, y_top, x_bot, k):
        super().__init__(ID, ins, outs, thermo)
        self.P = P
        self.LHK = LHK
        self.y_top = y_top
        self.x_bot = x_bot
        self.k = k
    
    @property
    def tray_spacing(self):
        return self._TS
    @tray_spacing.setter
    def tray_spacing(self, TS):
        """Tray spacing (225-600 mm)."""
        self._TS = TS
    
    @property
    def stage_efficiency(self):
        """Enforced user defined stage efficiency."""
        return self._stage_efficiency
    @stage_efficiency.setter
    def stage_efficiency(self, stage_efficiency):
        self._stage_efficiency = stage_efficiency
    
    @property
    def velocity_factor(self):
        """Ratio of actual velocity to maximum velocity allowable before flooding."""
        return self._f
    @velocity_factor.setter
    def velocity_factor(self, f):
        self._f = f
    
    @property
    def foaming_factor(self):
        """Foaming factor (0 to 1)."""
        return self._F_F
    @foaming_factor.setter
    def foaming_factor(self, foaming_factor):
        if foaming_factor > 1 or foaming_factor < 0:
            raise ValueError(f"foaming_factor must be between 0 and 1, ({foaming_factor} given).")
        self._F_F = foaming_factor
    
    @property
    def open_area_factor(self):
        """Ratio of open area, A_h, to active area, A_a."""
        return self._A_ha
    @open_area_factor.setter
    def open_area_factor(self, A_ha):
        self._A_ha = A_ha
    
    @property
    def downcomer_area_factor(self):
        """Enforced ratio of downcomer area to net (total) area. If None, estimate ratio based on Oliver's estimation [1]."""
        return self._A_dn
    @downcomer_area_factor.setter
    def downcomer_area_factor(self, A_dn):
        self._A_dn = A_dn
    
    @property
    def tray_type(self):
        """Default 'Sieve'"""
        return self._F_TTstr
    @tray_type.setter
    def tray_type(self, tray_type):
        try:
            self._F_TT = F_TTdict[tray_type]
        except KeyError:
            dummy = str(F_TTdict.keys())[11:-2]
            raise ValueError(f"tray_type must be one of the following: {dummy}")
        self._F_TTstr = tray_type
        
    @property
    def tray_material(self):
        """Default 'Carbon steel'"""
        return self._F_TMstr
    @tray_material.setter
    def tray_material(self, tray_material):
        try:
            self._F_TM = F_TMdict[tray_material]
        except KeyError:
            dummy = str(F_TMdict.keys())[11:-2]
            raise ValueError(f"tray_material must be one of the following: {dummy}")
        self._F_TMstr = tray_material        

    @property
    def vessel_material(self):
        """Default 'Carbon steel'"""
        return self._F_Mstr
    @vessel_material.setter
    def vessel_material(self, vessel_material):
        try:
            self._F_M = F_Mdict[vessel_material]
        except KeyError:
            dummy = str(F_Mdict.keys())[11:-2]
            raise ValueError(f"vessel_material must be one of the following: {dummy}")
        self._F_Mstr = vessel_material  
    
    def _mass_balance(self):
        # Get all important flow rates (both light and heavy keys and non-keys)
        get_index = self.chemicals.get_index
        LHK_index = get_index(self._LHK)
        LNK_index = get_index(self._LNK)
        HNK_index = get_index(self._HNK)
        mol = self.mol_in
        LHK_mol = mol[LHK_index]
        HNK_mol = mol[HNK_index]
        LNK_mol = mol[LNK_index]

        # Set light and heavy keys by lever rule
        light, heavy = LHK_mol
        LHK_F_mol = light + heavy
        zf = light/LHK_F_mol
        split_frac = (zf-self.x_bot)/(self.y_top-self.x_bot)
        top_net = LHK_F_mol*split_frac

        # Set output streams
        vap, liq = self.outs
        vap.mol[LHK_index] = vap_LHK_mol = top_net * self._y
        liq.mol[LHK_index] = LHK_mol - vap_LHK_mol
        vap.mol[LNK_index] = LNK_mol
        liq.mol[HNK_index] = HNK_mol
    
    def _run(self):
        self._mass_balance()
        vap, liq = self.outs
        vap.phase = 'g'
        liq.phase = 'l'
        vap.P = liq.P = self.P
        self._condensate_dew_point = dp = vap.dew_point_at_P()
        self._boilup_bubble_point = bp = liq.bubble_point_at_P()
        liq.T = bp.T
        vap.T = dp.T

    def _equilibrium_staircase(self, operating_line, x_stairs,
                               y_stairs, T_stairs, x_limit, bubble_T):
        """Find the specifications at every stage of the of the operating line before the maximum liquid molar fraction. Append the light key liquid molar fraction, light key vapor molar fraction, and stage temperatures to x_stairs, y_stairs and T_stairs respectively.
        
        Parameters
        ----------
        operating_line : function
                         Should return the liquid molar fraction of the light
                         key given its vapor molar fraction.
        x_stairs : list
                   Liquid molar compositions at each stage. Last element
                   should be the starting point for the next stage.
        y_stairs : list
                   Vapor molar compositions at each stage. Last element 
                   should be the starting point for the next stage.
        T_stairs : list
                   Temperature at each stage.
            
        """
        P = self.P
        i = 0
        yi = y_stairs[-1]
        xi = x_stairs[-1]
        while xi < x_limit:
            if i > 100:
                raise RuntimeError('cannot meet specifications! stages > 100')
            i += 1
            # Go Up
            T, y = bubble_T(array((xi, 1-xi)), P)
            yi = y[0]
            y_stairs.append(yi)
            T_stairs.append(T)
            # Go Right
            xi = operating_line(yi)
            if xi > 1:
                xi = x_limit
            x_stairs.append(xi)

    def _plot_stages(self):
        """Plot stages, graphical aid line, and equilibrium curve. The plot does not include operating lines nor a legend."""
        vap, liq = self.outs
        if not hasattr(self, 'x_stages'):
            self._design()
        x_stages = self._x_stages
        y_stages = self._y_stages
        LHK = self.LHK
        LK = self.LHK[0]
        P = self.P
        
        # Equilibrium data
        x_eq = np.linspace(0, 1, 100)
        y_eq = np.zeros(100)
        T = np.zeros(100)
        n = 0
        
        bp = vap.get_bubble_point(IDs=LHK)
        solve_Ty = bp.solve_Ty
        for xi in x_eq:
            T[n], y = solve_Ty(array([xi, 1-xi]), P)
            y_eq[n] = y[0]
            n += 1
            
        # Set-up graph
        plt.figure()
        plt.xticks(np.arange(0, 1.1, 0.1), fontsize=12)
        plt.yticks(fontsize=12)
        plt.xlabel('x (' + LK + ')', fontsize=16)
        plt.ylabel('y (' + LK + ')', fontsize=16)
        plt.xlim([0, 1])
        
        # Plot stages
        x_stairs = []
        for x in x_stages:
            x_stairs.append(x)
            x_stairs.append(x)
            
        y_stairs = []
        for y in y_stages:
            y_stairs.append(y)
            y_stairs.append(y)
        x_stairs.pop(-1)
        x_stairs.insert(0, y_stairs[0])
        plt.plot(x_stairs, y_stairs, '--')
        
        # Graphical aid line
        plt.plot([0, 1], [0, 1])
        
        # Vapor equilibrium graph
        plt.plot(x_eq, y_eq, lw=2)

    def _cost_trays(self, N_T, Di):
        """Return total cost of all trays at a CE of 500 given number of trays, `N_T`, and inner diameter, `Di`, in ft."""
        # Note: Can only use this function after running design method.
        C_BT = self._calc_TrayBaseCost(Di)
        F_NT = self._calc_NTrayFactor(N_T)
        return N_T*F_NT*self._F_TT*self._F_TM(Di)*C_BT

    def _cost_tower(self, Di, L, W):
        """Return cost of tower at a CE of 500 given the inner diameter, `Di`, and length, `L` in ft and weight, `W`, in lb."""
        C_V = self._calc_EmptyTowerCost(W)
        C_PL = self._calc_PlaformLadderCost(Di, L)
        return (self._F_M*C_V + C_PL)
    
    @staticmethod
    def _calc_EmptyTowerCost(W):
        """Return C_V the cost (USD) of an empty tower vessel assuming a CE of 500.
        
        **Parameters**
        
            W: Weight (lb)
        
        """
        return np.exp(7.2756 + 0.18255*np.log(W) + 0.02297*np.log(W)**2)
    
    @staticmethod
    def _calc_PlaformLadderCost(Di, L):
        """Return C_PL, the cost (USD) of platforms and ladders assuming a CE of 500.
        
        **Parameters**
        
            Di: Inner diameter (ft)
            L: Legnth (ft)
        
        """
        return 300.9*Di**0.63316*L**0.80161
    
    @staticmethod
    def _calc_Weight(Di, L, tv, rho_M):
        """Return W, the weight (lb) of the tower assuming 2:1 elliptical head.
        
        **Parameters**
        
            Di: Diameter (ft)
            L: Legnth (ft)
            tv: Shell thickness (in)
            rho_M: Density of material (lb/in^3)
        
        """
        Di = Di*12
        L = L*12
        return np.pi*(Di+tv)*(L+0.8*Di)*tv*rho_M
    
    @staticmethod
    def _calc_WallThickness(Po, Di, L, S=15000, E=None, M=29.5):
        """Return tv, the wall thinkness (in) designed to withstand the internal pressure and the wind/earthquake load at the bottom.
        
        Parameters
        ----------
        Po : float
            Operating internal pressure (psi)
        Di : float
            Internal diameter (ft)
        L : float
            Height (ft)
        S : float
            Maximum stress (psi)
        E : float
            Fractional weld efficiency
        M : float
            Elasticity (psi)
            
        """
        # TODO: Incorporate temperature for choosing S and M
        Di = Di*12 # ft to in
        L = L*12
        
        E_check = E is None
        if E_check:
            # Assume carbon steel with thickness more than 1.25 in
            E = 1.0 
        
        # Get design pressure, which should be higher than operating pressure.
        Po_gauge = Po - 14.69
        if Po_gauge < 0:
            # TODO: Double check vacuum calculation
            Pd = -Po_gauge
            tE = 1.3*Di*(Pd*L/M/Di)**0.4
            tEC = L*(0.18*Di - 2.2)*10**-5 - 0.19
            tv = tE + tEC
            return tv
        elif Po_gauge < 5:
            Pd = 10
        elif Po_gauge < 1000:
            Pd = np.exp(0.60608 + 0.91615*np.log(Po)) + 0.0015655*np.log(Po)**2
        else:
            Pd = 1.1*Po_gauge
        
        # Calculate thinkess according to ASME pressure-vessel code.
        ts = Pd*Di/(2*S*E-1.2*Pd)
        
        if E_check:
            # Weld efficiency of 0.85 for low thickness carbon steel
            if ts < 1.25:
                E = 0.85
                ts = Pd*Di/(2*S*E-1.2*Pd)
        
        # Add corrosion allowence
        ts += 1/8
        
        # Minimum thickness for vessel rigidity may be larger
        Di_ft = Di/12
        ts_min = np.polyval(ts_min_p, Di/12) if Di_ft < 4 else 0.25
        if ts < ts_min:
            ts = ts_min
        
        # Calculate thickness to withstand wind/earthquake load
        Do = Di + ts
        tw = 0.22*(Do + 18)*L**2/(S*Do**2)
        tv = max(tw, ts)
        
        # Vessels are fabricated from metal plates with small increments
        if tv < 0.5:
            tv = approx2step(tv, 3/16, 1/16)
        elif tv < 2:
            tv = approx2step(tv, 0.5, 1/8)
        elif tv < 3:
            tv = approx2step(tv, 2, 1/4)
        return tv
    
    @staticmethod
    def _calc_TrayBaseCost(Di):
        """Return C_BT, the base cost of a tray (USD) at a CE of 500.
        
        **Parameters**
        
            Di: Inner diameter (ft)
        
        """
        return 412.6985 * np.exp(0.1482*Di)

    @staticmethod
    def _calc_NTrayFactor(N_T):
        """Return F_NT, the cost factor for number of trays.
        
        **Parameters**
        
            N_T: Number of trays
            
        """
        if N_T < 20:
            F_NT = 2.25/1.0414**N_T
        else:
            F_NT = 1
        return F_NT

    @staticmethod
    def _calc_MurphreeEfficiency(mu, alpha, L, V):
        """Return E_mv, the sectional murphree efficiency.
        
        **Parameters**
            
            mu: Viscosity (mPa*s)
            
            alpha: Relative volatility
            
            L: Liquid flow rate by mol
            
            V: Vapor flow rate by mol
        
        """
        S = alpha*V/L # Stripping factor
        e = 0.503*mu**(-0.226)*(S if S > 1 else 1/S)**(-0.08 )
        if e < 1: return e
        else: return 1
    
    @staticmethod
    def _calc_FlowParameter(L, V, rho_V, rho_L):
        """Return F_LV, the flow parameter.
        
        **Parameters**
        
            L: Liquid flow rate by mass
            V: Vapor flow rate by mass
            rho_V: Vapor density
            rho_L: Liquid density
        
        """
        return L/V*(rho_V/rho_L)**0.5
    
    @staticmethod
    def _calc_MaxCapacityParameter(TS, F_LV):
        """Return C_sbf, the maximum capacity parameter before flooding (m/s).
        
        **Parameters**
        
            TS :
                Tray spacing (mm)
            F_LV :
                Flow parameter
        
        """
        return 0.0105 + 8.127e-4*TS**0.755*np.exp(-1.463*F_LV**0.842)
    
    @staticmethod
    def _calc_MaxVaporVelocity(C_sbf, sigma,
                               rho_L, rho_V, F_F, A_ha):
        """Return U_f, the maximum allowable vapor velocity through the net area of flow before flooding (m/s).
        
        Parameters
        ----------
        
        C_sbf : 
            Maximum Capacity Parameter (m/s)
        sigma : 
            Liquid surface tension (dyn/cm)
        rho_L : 
            Liquid density
        rho_V : 
            Vapor density
        F_F : 
            Foaming factor
        A_ha : 
            Ratio of open area, A_h, to active area, A_a
        
        """
        F_ST = (sigma/20)**0.2 # Surface tension factor
        
        # Working area factor
        if A_ha >= 0.1 and A_ha <= 1:
            F_HA = 1
        elif A_ha >= 0.06:
            F_HA = 5*A_ha + 0.5
        else:
            raise ValueError(f"ratio of open to active area, 'A', must be between 0.06 and 1 ({A_ha} given)") 
        
        return C_sbf * F_HA * F_ST * ((rho_L-rho_V)/rho_V)**0.5
    
    @staticmethod
    def _calc_DowncomerAreaRatio(F_LV):
        """Return the ratio of downcomer area to net (total) area, `A_dn`.
        
        **Parameters**
        
            F_LV : float
                Flow parameter
        
        """
        if F_LV < 0.1:
            A_dn = 0.1
        elif F_LV < 1:
            A_dn = 0.1 + (F_LV-0.1)/9
        else:
            A_dn = 0.2
        return A_dn
    
    @staticmethod
    def _calc_Diameter(V_vol, U_f, f, A_dn):
        """Return D_T, the column diameter (ft).
        
        **Parameters**
        
            V_vol: Vapor volumetric flow rate (m^3/s)
            U_f: Maximum vapor velocity before flooding(m/s)
            f: Ratio of actual velocity to U_f
            A_dn: ratio of downcomer area to net (total) area
        
        """
        Di = (4*V_vol/(f*U_f*np.pi*(1-A_dn)))**0.5
        if Di < 0.914:
            # Make sure diameter is not too small
            Di = 0.914
        Di *= 3.28
        return Di
    
    @staticmethod
    def _calc_Height(TS, Nstages: int, top=True, bot=True):
        """Return H, the height of the column (ft).
        
        **Parameters**
        
            TS: Tray spacing (mm)
            Nstages: Number of stages 
        
        """
        # 3 m bottoms surge capacity, 1.25 m above top tray to remove entrained liquid
        H = TS*Nstages/1000
        if top:
            H += 1.2672
        if bot:
            H += 3
        H *= 3.28
        return H 

    def _calc_condenser(self):
        distillate = self.outs[0]
        condensate = self._condensate # Abstract attribute
        condenser = self._condenser
        s_in = condenser.ins[0]
        s_in.mol[:] = distillate.mol+condensate.mol
        s_in.T = distillate.T
        s_in.P = distillate.P
        ms1 = condenser.outs[0]
        ms1.imol['l'] = condensate.mol
        ms1.T = condensate.T
        ms1.P = condensate.P
        ms1.imol['g'] = distillate.mol
        condenser._design()
        condenser._cost()
        
    def _calc_boiler(self):
        bottoms = self.outs[1]
        boil_up = self._boil_up # Abstract attribute
        boiler = self._boiler
        s_in = boiler.ins[0]
        s_in.copy_like(bottoms)
        s_in.mol += boil_up.mol
        ms1 = boiler.outs[0]
        ms1.T = boil_up.T
        ms1.P = boil_up.P
        ms1.imol['g'] = boil_up.mol
        ms1.imol['l'] = bottoms.mol
        if hasattr(self, '_condenser'):
            boiler._design(self.H_out - self.H_in - self._condenser._duty)
            boiler._cost()
        else:
            boiler._design()
            boiler._cost()
        
    def _cost(self):
        Design = self.design_results
        Cost = self.purchase_costs
        F_CE = bst.CE/500
        
        # Cost trays assuming a partial condenser
        N_T = Design['Actual stages'] - 1
        Di = Design['Diameter']
        Cost['Trays'] = F_CE*self._cost_trays(N_T, Di)
        
        # Cost vessel assuming T < 800 F
        W = Design['Weight'] # in lb
        L = Design['Height']*3.28 # in ft
        Cost['Tower'] = F_CE*self._cost_tower(Di, L, W)
        self._cost_components()
        

class Distillation(Dist):
    """
    Create a Distillation column that assumes all light and heavy non keys
    separate to the top and bottoms product respectively. McCabe-Thiele
    analysis is used to find both the number of stages and the reflux ratio
    given a ratio of actual reflux to minimum reflux [1]_. This assumption
    is good for both binary distillation of highly polar compounds and
    ternary distillation assuming complete separation of light non-keys
    and heavy non-keys with large differences in boiling points. Preliminary
    analysis showed that the theoretical number of stages using this method
    on Methanol/Glycerol/Water systems is off by less than +-1 stage. Other
    methods, such as the Fenske-Underwood-Gilliland method, are more suitable
    for hydrocarbons. The Murphree efficiency is based on the modified
    O'Connell correlation [2]_. The diameter is based on tray separation
    and flooding velocity [1]_. Purchase costs are based on correlations
    by Mulet et al. [3, 4]_ as compiled by Warren et. al. [5]_.

    Parameters
    ----------
    ins
        [:] All input streams
    outs
        [0] Distillate product
        
        [1] Bottoms product
    LHK : tuple[str]
          Light and heavy keys.
    P=101325 : float
        Operating pressure (Pa).
    y_top : float
            Molar fraction of light key in the distillate.
    x_bot : float
            Molar fraction of light key in the bottoms.
    k : float
        Ratio of reflux to minimum reflux.

    References
    ----------
    .. [1] J.D. Seader, E.J. Henley, D.K. Roper. (2011)
        Separation Process Principles 3rd Edition. John Wiley & Sons, Inc. 

    .. [2] M. Duss, R. Taylor. (2018)
        Predict Distillation Tray Efficiency. AICHE 
    
    .. [3] Mulet, A., A. B. Corripio, and L. B. Evans. (1981a).
        Estimate Costs of Pressure Vessels via Correlations.
        Chem. Eng., 88(20), 145–150.

    .. [4] Mulet, A., A.B. Corripio, and L.B.Evans. (1981b).
        Estimate Costs of Distillation and Absorption Towers via Correlations.
        Chem. Eng., 88(26), 77–82.

    .. [5] Seider, W. D., Lewin,  D. R., Seader, J. D., Widagdo, S., Gani, R.,
        & Ng, M. K. (2017). Product and Process Design Principles. Wiley.
        Cost Accounting and Capital Cost Estimation (Chapter 16)    

    Examples
    --------
    Binary distillation assuming 100% separation on non-keys:
    
    >>> from biosteam.units import Distillation
    >>> from thermosteam import Stream, Chemicals, settings
    >>> chemicals = Chemicals(['Water', 'Methanol', 'Glycerol'])
    >>> settings.set_thermo(chemicals)
    >>> feed = Stream('feed', flow=(80, 100, 25))
    >>> bp = feed.bubble_point_at_P()
    >>> feed.T = bp.T # Feed at bubble point T
    >>> D1 = Distillation('D1', ins=feed,
    ...                   LHK=('Methanol', 'Water'),
    ...                   y_top=0.99, x_bot=0.01, k=2)
    >>> D1.is_divided = True
    >>> D1.simulate()
    >>> # See all results
    >>> D1.show(T='degC', P='atm', composition=True)
    Distillation: D1
    ins...
    [0] feed
        phase: 'l', T: 76.129 degC, P: 1 atm
        composition: Water     0.39
                     Methanol  0.488
                     Glycerol  0.122
                     --------  205 kmol/hr
    outs...
    [0] s1
        phase: 'g', T: 64.91 degC, P: 1 atm
        composition: Water     0.01
                     Methanol  0.99
                     --------  100 kmol/hr
    [1] s2
        phase: 'l', T: 100.06 degC, P: 1 atm
        composition: Water     0.754
                     Methanol  0.00761
                     Glycerol  0.239
                     --------  105 kmol/hr
    
    >>> D1.results()
    Distillation                                    Units        D1
    Cooling water       Duty                        kJ/hr -5.11e+06
                        Flow                      kmol/hr  3.49e+03
                        Cost                       USD/hr      1.71
    Low pressure steam  Duty                        kJ/hr  9.49e+06
                        Flow                      kmol/hr       244
                        Cost                       USD/hr      58.1
    Design              Theoretical feed stage                    9
                        Theoretical stages                       13
                        Minimum reflux              Ratio     0.687
                        Reflux                      Ratio      1.37
                        Rectifier stages                         14
                        Stripper stages                          13
                        Rectifier height               ft      33.2
                        Stripper height                ft      31.7
                        Rectifier diameter             ft      3.93
                        Stripper diameter              ft      3.01
                        Rectifier wall thickness       in      0.25
                        Stripper wall thickness        in      0.25
                        Rectifier weight               lb  4.61e+03
                        Stripper weight                lb  3.32e+03
    Purchase cost       Rectifier trays               USD  1.45e+04
                        Stripper trays                USD  1.21e+04
                        Rectifier tower               USD  7.41e+04
                        Stripper tower                USD   6.1e+04
                        Condenser                     USD  3.41e+04
                        Boiler                        USD  2.64e+04
    Total purchase cost                               USD  2.22e+05
    Utility cost                                   USD/hr      59.8
    
    """
    line = 'Distillation'
    _N_heat_utilities = 0
    _graphics = Dist._graphics
    _is_divided = False #: [bool] True if the stripper and rectifier are two separate columns.
    _units_not_divided = {'Minimum reflux': 'Ratio',
                          'Reflux': 'Ratio',
                          'Rectifier height': 'ft',
                          'Rectifier diameter': 'ft',
                          'Rectifier wall thickness': 'in',
                          'Rectifier weight': 'lb',
                          'Stripper height': 'ft',
                          'Stripper diameter': 'ft',
                          'Stripper wall thickness': 'in',
                          'Stripper weight': 'lb',
                          'Height': 'ft',
                          'Diameter': 'ft',
                          'Wall thickness': 'in',
                          'Weight': 'lb'}
    _units = _units_not_divided
    _units_divided = {'Minimum reflux': 'Ratio',
                      'Reflux': 'Ratio',
                      'Rectifier height': 'ft',
                      'Rectifier diameter': 'ft',
                      'Rectifier wall thickness': 'in',
                      'Rectifier weight': 'lb',
                      'Stripper height': 'ft',
                      'Stripper diameter': 'ft',
                      'Stripper wall thickness': 'in',
                      'Stripper weight': 'lb'}
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 P=101325, *, LHK, y_top, x_bot, k):
        super().__init__(ID, ins, outs, thermo, P, LHK=LHK,
                         y_top=y_top, x_bot=x_bot, k=k)
        thermo = self._load_thermo(thermo)
        self._condenser = HXutility(None,
                                    ins=Stream(None, phase='g', thermo=thermo),
                                    outs=MultiStream(None, thermo=thermo),
                                    thermo=thermo)
        self._boiler = HXutility(None,
                                 ins=Stream(None, thermo=thermo),
                                 outs=MultiStream(None, thermo=thermo),
                                 thermo=thermo)
        self.heat_utilities = self._condenser.heat_utilities + self._boiler.heat_utilities
        self._condensate = Stream(None, thermo=thermo)
        self._boil_up = Stream(None, phase='g', thermo=thermo)
        self._vapor_stream =  Stream(None)
        self._McCabeThiele_args = np.zeros(6)
    
    @property
    def is_divided(self):
        """[bool] True if the stripper and rectifier are two separate columns."""
        return self._is_divided
    
    @is_divided.setter
    def is_divided(self, is_divided):
        self._is_divided = is_divided
        self._units = self._units_divided if is_divided else self._units_not_divided
    
    def _calc_Nstages(self):
        """Return a tuple with the actual number of stages for the rectifier and the stripper."""
        vap, liq = self.outs
        Design = self.design_results
        x_stages = self._x_stages
        y_stages = self._y_stages
        R = Design['Reflux']
        N_stages = Design['Theoretical stages']
        feed_stage = Design['Theoretical feed stage']
        liq_mol = self._feed_liqmol
        vap_mol = self._feed_vapmol
        
        stage_efficiency = self.stage_efficiency
        if stage_efficiency:
            return N_stages/stage_efficiency
        else:    
            # Calculate Murphree Efficiency for rectifying section
            vap_F_mol = vap.F_mol
            dp = self._condensate_dew_point
            condensate_x_mol = dp.x
            self._L_Rmol = L_Rmol = R*vap_F_mol
            self._V_Rmol = V_Rmol = (R+1)*vap_F_mol
            condensate = self._condensate
            condensate.imol[dp.IDs] = condensate_x_mol * L_Rmol
            condensate.T = vap.T
            condensate.P = vap.P
            mu = 1000 * condensate.mu # mPa*s
            K_light = y_stages[-1]/x_stages[-1] 
            K_heavy = (1-y_stages[-1])/(1-x_stages[-1])
            alpha = K_light/K_heavy
            self._E_rectifier = E_rectifier = self._calc_MurphreeEfficiency(mu, alpha, L_Rmol, V_Rmol)
            
            # Calculate Murphree Efficiency for stripping section
            mu = 1000 * liq.mu # mPa*s
            self._V_Smol = V_Smol = (R+1)*vap_F_mol - sum(vap_mol)
            self._L_Smol = L_Smol = R*vap_F_mol + sum(liq_mol) 
            K_light = y_stages[0]/x_stages[0] 
            K_heavy = (1-y_stages[0])/(1-x_stages[0] )
            alpha = K_light/K_heavy
            self._E_stripper = E_stripper = self._calc_MurphreeEfficiency(mu, alpha, L_Smol, V_Smol)
            
            # Calculate actual number of stages
            mid_stage = feed_stage - 0.5
            N_rectifier = np.ceil(mid_stage/E_rectifier)
            N_stripper = np.ceil((N_stages-mid_stage)/E_stripper)
        return N_rectifier, N_stripper

    def _McCabeThiele(self):
        distillate, bottoms = self.outs
        chemicals = self.chemicals
        LHK = self._LHK
        LHK_index = chemicals.get_index(LHK)

        # Feed light key mol fraction
        liq_mol = np.zeros(chemicals.size)
        vap_mol = liq_mol.copy()
        for s in self.ins:
            if isinstance(s, MultiStream):
                liq_mol += s.imol['l']
                vap_mol += s.imol['g']
            elif s.phase == 'g':
                vap_mol += s.mol
            elif s.phase.lower() == 'l':
                liq_mol += s.mol
            elif s.phase == 's': pass
            else:
                raise RuntimeError(f'invalid phase encountered in {repr(s)}')
        self._feed_liqmol = liq_mol
        self._feed_vapmol = vap_mol
        LHK_mol = liq_mol[LHK_index] + vap_mol[LHK_index]
        LHK_F_mol = LHK_mol.sum()
        zf = LHK_mol[0]/LHK_F_mol
        
        # Get feed quality
        q = liq_mol[LHK_index].sum()/LHK_F_mol
        
        # Main arguments
        P = self.P
        k = self.k
        y_top = self.y_top
        x_bot = self.x_bot
        
        # Cache
        old_args = self._McCabeThiele_args
        args = np.array([P, k, y_top, x_bot, q, zf])
        if (abs(old_args - args) < np.array([50, 1e-5, 1e-6, 1e-6, 1e-6, 1e-6], float)).all(): return
        self._McCabeThiele_args = args
        
        # Get R_min and the q_line 
        if q == 1:
            q = 1 - 1e-5
        q_line = lambda x: q*x/(q-1) - zf/(q-1)
        self._q_line_args = dict(q=q, zf=zf)
        
        solve_Ty = bottoms.get_bubble_point(LHK).solve_Ty
        Rmin_intersection = lambda x: q_line(x) - solve_Ty(array((x, 1-x)), P)[1][0]
        x_Rmin = brentq(Rmin_intersection, 0, 1)
        y_Rmin = q_line(x_Rmin)
        m = (y_Rmin-y_top)/(x_Rmin-y_top)
        Rmin = m/(1-m)
        if Rmin <= 0:
            R = 0.0001*k
        else:
            R = Rmin*k

        # Rectifying section: Inntersects q_line with slope given by R/(R+1)
        m1 = R/(R+1)
        b1 = y_top-m1*y_top
        rs = lambda y: (y - b1)/m1 # -> x
        
        # y_m is the solution to lambda y: y - q_line(rs(y))
        self._y_m = y_m = (q*b1 + m1*zf)/(q - m1*(q-1))
        self._x_m = x_m = rs(y_m)
        
        # Stripping section: Intersects Rectifying section and q_line and beggins at bottoms liquid composition
        m2 = (x_bot-y_m)/(x_bot-x_m)
        b2 = y_m-m2*x_m
        ss = lambda y: (y-b2)/m2 # -> x        
        
        # Data for staircase
        self._x_stages = x_stages = [x_bot]
        self._y_stages = y_stages = [x_bot]
        self._T_stages = T_stages = []
        self._equilibrium_staircase(ss, x_stages, y_stages, T_stages, x_m, solve_Ty)
        yi = y_stages[-1]
        xi = rs(yi)
        x_stages[-1] = xi if xi < 1 else 0.99999
        self._equilibrium_staircase(rs, x_stages, y_stages, T_stages, y_top, solve_Ty)
        
        # Find feed stage
        for i in range(len(y_stages)-1):
            j = i+1
            if y_stages[i] < y_m and y_stages[j] > y_m:
                feed_stage = i+1
        stages = len(x_stages)
        
        # Results
        Design = self.design_results
        Design['Theoretical feed stage'] = stages - feed_stage
        Design['Theoretical stages'] = stages
        Design['Minimum reflux'] = Rmin
        Design['Reflux'] = R
        

    def _design(self):
        self._McCabeThiele()
        distillate, bottoms = self._outs
        Design = self.design_results
        R = Design['Reflux']
        Rstages, Sstages = self._calc_Nstages()
        calc_Height = self._calc_Height
        is_divided = self.is_divided
        TS = self._TS
        
        ### Get diameter of rectifying section based on top plate ###
        
        condensate = self._condensate
        rho_L = condensate.rho
        sigma = 1000 * condensate.sigma # dyn/cm
        L = condensate.F_mass
        V = L*(R+1)/R
        vapor_stream = self._vapor_stream
        vapor_stream.copy_like(distillate)
        vapor_stream.mol *= R+1
        V_vol = 0.0002778 * vapor_stream.F_vol # m^3/s
        rho_V = distillate.rho
        F_LV = self._calc_FlowParameter(L, V, rho_V, rho_L)
        C_sbf = self._calc_MaxCapacityParameter(TS, F_LV)
        F_F = self._F_F
        A_ha = self._A_ha
        U_f = self._calc_MaxVaporVelocity(C_sbf, sigma, rho_L, rho_V, F_F, A_ha)
        A_dn = self._A_dn
        if A_dn is None:
           self._A_dn = A_dn = self._calc_DowncomerAreaRatio(F_LV)
        f = self._f
        R_diameter = self._calc_Diameter(V_vol, U_f, f, A_dn)
        
        ### Get diameter of stripping section based on feed plate ###
        
        V_mol = self._V_Smol
        rho_L = bottoms.rho
        bp = self._boilup_bubble_point
        boil_up_flow = bp.y * V_mol
        boil_up = self._boil_up
        boil_up.T = bottoms.T
        boil_up.P = bottoms.P
        boil_up.imol[bp.IDs] = boil_up_flow
        V = boil_up.F_mass
        V_vol = 0.0002778 * boil_up.F_vol # m^3/s
        rho_V = boil_up.rho
        L = bottoms.F_mass # To get liquid going down
        F_LV = self._calc_FlowParameter(L, V, rho_V, rho_L)
        C_sbf = self._calc_MaxCapacityParameter(TS, F_LV)
        sigma = 1000 * bottoms.sigma
        U_f = self._calc_MaxVaporVelocity(C_sbf, sigma, rho_L, rho_V, F_F, A_ha)
        A_dn = self._A_dn
        if A_dn is None:
            self.A_dn = self._calc_DowncomerAreaRatio(F_LV)
        S_diameter = self._calc_Diameter(V_vol, U_f, f, A_dn)
        Po = self.P*0.000145078 # to psi
        rho_M = rho_Mdict[self._F_Mstr]
        
        if is_divided:
            Design['Rectifier stages'] = Rstages
            Design['Stripper stages'] =  Sstages
            Design['Rectifier height'] = H_R = calc_Height(TS, Rstages-1)
            Design['Stripper height'] = H_S = calc_Height(TS, Sstages-1)
            Design['Rectifier diameter'] = R_diameter
            Design['Stripper diameter'] = S_diameter
            Design['Rectifier wall thickness'] = tv = self._calc_WallThickness(Po, R_diameter, H_R)
            Design['Stripper wall thickness'] = tv = self._calc_WallThickness(Po, S_diameter, H_S)
            Design['Rectifier weight'] = self._calc_Weight(R_diameter, H_R, tv, rho_M)
            Design['Stripper weight'] = self._calc_Weight(S_diameter, H_S, tv, rho_M)
        else:
            Design['Actual stages'] = Rstages + Sstages
            Design['Height'] = H = calc_Height(TS, Rstages+Sstages-2)
            Design['Diameter'] = Di = max((R_diameter, S_diameter))
            Design['Wall thickness'] = tv = self._calc_WallThickness(Po, Di, H)
            Design['Weight'] = self._calc_Weight(Di, H, tv, rho_M)
        
    def _cost(self):
        if not self.is_divided: return super()._cost()
        Design = self.design_results
        Cost = self.purchase_costs
        F_CE = bst.CE/500
        
        # Number of trays assuming a partial condenser
        N_RT = Design['Rectifier stages'] - 1
        Di_R = Design['Rectifier diameter']
        Cost['Rectifier trays'] = F_CE*self._cost_trays(N_RT, Di_R)
        N_ST = Design['Stripper stages'] - 1
        Di_S = Design['Stripper diameter']
        Cost['Stripper trays'] = F_CE*self._cost_trays(N_ST, Di_S)
        
        # Cost vessel assuming T < 800 F
        W_R = Design['Rectifier weight'] # in lb
        H_R = Design['Rectifier height']*3.28 # in ft
        Cost['Rectifier tower'] = F_CE*self._cost_tower(Di_R, H_R, W_R)
        W_S = Design['Stripper weight'] # in lb
        H_S = Design['Stripper height']*3.28 # in ft
        Cost['Stripper tower'] = F_CE*self._cost_tower(Di_S, H_S, W_S)
        self._cost_components()
    
    def _cost_components(self): 
        # Cost condenser
        self._calc_condenser()
        self.purchase_costs['Condenser'] = self._condenser.purchase_costs['Heat exchanger']
        
        # Cost boiler
        self._calc_boiler()
        self.purchase_costs['Boiler'] = self._boiler.purchase_costs['Heat exchanger']
        
    
    def plot_stages(self):
        """Plot the McCabe Thiele Diagram."""
        # Plot stages, graphical aid and equilibrium curve
        self._plot_stages()
        vap, liq = self.outs
        Design = self.design_results
        if not hasattr(self, '_x_stages'): self._design()
        q_args = self._q_line_args
        zf = q_args['zf']
        q = q_args['q']
        q_line = lambda x: q*x/(q-1) - zf/(q-1)
        y_top = self.y_top
        x_bot = self.x_bot
        stages = Design['Theoretical stages']
        Rmin = Design['Minimum reflux']
        R = Design['Reflux']
        feed_stage = Design['Theoretical feed stage']
        
        # q_line
        intersect2 = lambda x: x - q_line(x)
        x_m2 = brentq(intersect2, 0, 1)
        
        # Graph q-line, Rectifying and Stripping section
        plt.plot([self._x_m, x_m2], [self._y_m, x_m2])
        plt.plot([self._x_m, y_top], [self._y_m, y_top])
        plt.plot([x_bot, self._x_m], [x_bot, self._y_m])
        plt.legend([f'Stages: {stages}, Feed: {feed_stage}', 'Graphical aid', 'eq-line', 'q-line', 'ROL', 'SOL'], fontsize=12)
        plt.title(f'McCabe Thiele Diagram (Rmin = {Rmin:.2f}, R = {R:.2f})')
        plt.show()
        return plt


# class Stripper(Dist):
#     line = 'Stripper'
#     __doc__ = column_doc.replace('{Column Type}', 'Stripper')
#     _N_heat_utilities = 0
#     _graphics = Dist._graphics
#     _units ={'Minimum boil up': 'Ratio',
#              'Boil up': 'Ratio',
#              'Height': 'ft',
#              'Diameter': 'ft',
#              'Wall thickness': 'in',
#              'Weight': 'lb'}
    
#     def _init(self):
#         self._boiler = HXutility(None,
#                                  ins=Stream(None),
#                                  outs=MixedStream(None))
#         self.heat_utilities = self._boiler.heat_utilities
#         self._boil_up = Stream(None)
#         self._McCabeThiele_args = np.array([0, 0, 0, 0])
    
#     def plot_stages(self):
#         # Plot stages, graphical aid and equilibrium curve
#         self._plot_stages()
        
#         # Cached Data
#         vap, liq = self.outs
#         Design = self.design_results
#         cached = self._cached
#         x_stages = cached.get('x_stages')
#         if not x_stages:
#             self._design()
#             x_stages = cached.get('x_stages')
#         y_top = self.y_top
#         x_bot = self.x_bot
#         stages = Design['Theoretical stages']
#         Bmin = Design['Minimum boil up']
#         B = Design['Boil up']
        
#         # Graph Stripping section
#         ss = self._ss
#         plt.plot([x_bot, ss(y_top)], [x_bot, y_top])
#         plt.legend([f'Stages: {stages}', 'Graphical aid', 'eq-line', 'SOL'], fontsize=12)
            
#         # Title
#         plt.title(f'McCabe Thiele Diagram (Bmin = {Bmin:.2f}, B = {B:.2f})')
#         plt.show()
#         return plt
    
#     def _calc_Nstages(self):
#         """Return the actunal number of stages"""
#         vap, liq = self.outs
#         cached = self._cached
#         Design = self.design_results
#         x_stages = cached['x_stages']
#         y_stages = cached['y_stages']
#         LHK_index = cached['LHK_index']
#         B = Design['Boil up']
#         N_stages = Design['Theoretical stages']
        
#         stage_efficiency =self.stage_efficiency
#         if stage_efficiency:
#             return N_stages/stage_efficiency
#         else:
#             # Calculate Murphree Efficiency for stripping section
#             mu = 1000 * liq.mu # mPa*s
#             cached['V_mol'] = V_mol = B*sum(liq.mol[LHK_index])
#             cached['L_mol'] = L_mol = liq.F_mol + V_mol
#             K_light = y_stages[0]/x_stages[0] 
#             K_heavy = (1-y_stages[0])/(1-x_stages[0] )
#             alpha = K_light/K_heavy
#             cached['Stripping Section Efficiency'] = E_stripper = self._calc_MurphreeEfficiency(mu, alpha, L_mol, V_mol)
#         return np.ceil(N_stages/E_stripper)
    
#     def _McCabeThiele(self):
#         distillate, bottoms = self.outs
#         Design = self.design_results
#         cached = self._cached
        
#         # Some important info
#         species = self._LHK_chemicals
#         P = self.P
#         k = self.k
#         y_top = self.y_top
#         x_bot = self.x_bot
        
#         old_args = self._McCabeThiele_args
#         args = np.array([y_top, x_bot, k, P])
#         if (abs(old_args - args) < np.array([0.00001, 0.00001, 0.01, 50], float)).all(): return
#         self._McCabeThiele_args = args
        
#         # Get B_min (Boil up ratio)
#         y_Rmin = y_top
#         T_guess, x = distillate._dew_T(species, [y_Rmin, 1-y_Rmin], P)
#         x_Rmin = x[0]
#         m = (y_Rmin-x_bot)/(x_Rmin-x_bot)
#         Bmin = 1/(m-1)
#         B = 0.1*k if Bmin <= 0 else Bmin*k

#         # Stripping section: Inntersects liquid composition with slope given by (B+1)/B
#         m = (B+1)/B
#         b = x_bot-m*x_bot
#         def ss(y): return (y - b)/m
#         self._ss = ss
                
#         # Data for staircase
#         cached['x_stages'] = x_stages = [x_bot]
#         cached['y_stages'] = y_stages = [x_bot]
#         cached['T_stages'] = T_stages = []
#         self._equilibrium_staircase(ss, x_stages, y_stages, T_stages, ss(y_top))
#         stages = len(x_stages)
        
#         # Results
#         Design['Theoretical stages'] = stages
#         Design['Minimum boil up'] = Bmin
#         Design['Boil up'] = B
    
#     def _design(self):
#         distillate, bottoms = self.outs
#         Design = self.design_results
#         cached = self._cached
        
#         ### Get number of stages and height ###
        
#         self._McCabeThiele()
#         TS = self._TS
#         Design['Actual stages'] = Nstages = self._calc_Nstages()
#         Design['Height'] = H = self._calc_Height(TS, Nstages-1)
        
#         ### Get diameter of stripping section based on feed plate ###
        
#         V_mol = cached['V_mol']
#         L_mol = cached['L_mol']
#         rho_L = bottoms.rho
#         boil_up_flow = cached['boilup_z_mol'] * V_mol
#         boil_up = cached['boil_up']
#         boil_up.T = bottoms.T; boil_up.P = bottoms.P
#         lookup = bottoms._species._compounds.index
#         index_ = [lookup(i) for i in cached['vle_bot']]
#         boil_up.mol[index_] = boil_up_flow
#         V = boil_up.F_mass
#         V_vol = 0.0002778 * boil_up.F_vol # m^3/s
#         rho_V = boil_up.rho
#         L = sum(bottoms.MW*bottoms.z_mol*L_mol) # To get liquid going down
#         F_LV = self._calc_FlowParameter(L, V, rho_V, rho_L)
#         C_sbf = self._calc_MaxCapacityParameter(TS, F_LV)
#         sigma = 1000 * bottoms.sigma # dyn/cm
#         U_f = self._calc_MaxVaporVelocity(C_sbf, sigma, rho_L, rho_V, self._F_F, self._A_ha)
#         A_dn = self._A_dn or self._calc_DowncomerAreaRatio(F_LV)
#         Design['Diameter'] = Di = self._calc_Diameter(V_vol, U_f, self._f, A_dn)
#         Po = self.P/101325*14.7
#         Design['Wall thickness'] = tv = self._calc_WallThickness(Po, Di, H)
#         rho_M = rho_Mdict[self._F_Mstr]
#         Design['Weight'] = self._calc_Weight(Di, H, tv, rho_M)
    
#     def _cost_components(self):
#         # Cost boiler
#         self._calc_boiler()
#         self.purchase_costs['Boiler'] = self._boiler.purchase_costs['Heat exchanger']

    