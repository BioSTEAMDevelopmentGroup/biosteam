#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 19:33:20 2018

@author: yoelr
"""
from biosteam import Unit, np, MixedStream, Stream
from biosteam.utils import get_vals, approx2step
from biosteam.exceptions import biosteam_Error, DesignError
from scipy.optimize import brentq
from biosteam.units.hx import HXutility
import matplotlib.pyplot as plt

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
F_TMdict = {'Carbon steel': lambda Di: 1,
            '303 Stainless steel': lambda Di: 1.189 + 0.058*Di,
            '316 Stainless steel': lambda Di: 1.401 + 0.073*Di,
            'Carpenter 20CB-3': lambda Di: 1.525 + 0.079*Di,
            'Monel': lambda Di: 2.306 + 0.112*Di}

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
            'Stainless steel 304': None,
            'Stainless steel 316': None,
            'Carpenter 20CB-3': None,
            'Nickel-200': None,
            'Monel-400': None,
            'Inconel-600': None,
            'Incoloy-825': None,
            'Titanium': None}

# %% Distillation

# Abstract doc string for columns
column_doc = """Create a {Column Type} column that assumes all light and heavy non keys separate to the top and bottoms product respectively. McCabe-Thiele analysis is used to find both the number of stages and the reflux ratio given a ratio of actual reflux to minimum reflux [1]. This assumption is good for both binary distillation of highly polar compounds and ternary distillation assuming complete separation of light non-keys and heavy non-keys with large differences in boiling points. Preliminary analysis showed that the theoretical number of stages using this method on Methanol/Glycerol/Water systems is off by less than +-1 stage. Other methods, such as the Fenske-Underwood-Gilliland method, are more suitable for hydrocarbons. The Murphree efficiency is based on the modified O'Connell correlation [2]. The diameter is based on tray separation and flooding velocity [1].

    **Parameters**

        **LHK:** *tuple[str]* Light and heavy keys

        **P:** *[float]* Operating pressure (Pa)

        **y_top:** *[float]* Molar fraction of light key in the distillate

        **x_bot:** *[float]* Molar fraction of light key in the bottoms

        **k:** *[float]* Ratio of reflux to minimum reflux 

    **ins**
        
        [:] All input streams

    **outs**
        
        [0] Distillate product
        
        [1] Bottoms product

    **References**

        [1] J.D. Seader, E.J. Henley, D.K. Roper. Separation Process Principles 3rd Edition. John Wiley & Sons, Inc. (2011)
    
        [2] M. Duss, R. Taylor. Predict Distillation Tray Efficiency. AICHE (2018)
    
    """

class Dist(Unit):
    """Abstract class for a column."""
    # Column material factor
    _F_Mstr = 'Carbon steel'
    _F_M = 1
    
    # Tray type factor
    _F_TTstr = 'Sieve'
    _F_TT = 1
    
    # Tray material factor function
    _F_TMstr = 'Carbon steel'
    _F_TM = staticmethod(lambda Di: 1)
    
    # [float] Tray spacing (225-600 mm)
    _TS = 450 
    
    #: [float] Enforced user defined stage efficiency.        
    _StageEfficiency = None
    
    # [float] Ratio of actual velocity to maximum velocity allowable before flooding.
    _f = 0.8 
    
    # [float] Foaming factor (0 < F_F < 1)
    _F_F = 1
    
    # [float] Ratio of open area, A_h, to active area, A_a
    _A_ha = 0.1
    
    # [float] Enforced ratio of downcomer area to net (total) area. If None, estimate ratio based on Oliver's estimation [1].
    _A_dn = None
    
    # [dict] Bounds for results
    bounds = {'Diameter': (3., 24.),
              'Height': (27., 170.),
              'Weight': (9000., 2.5e6),
              'Stripper diameter': (3., 24.),
              'Stripper height': (27., 170.),
              'Stripper weight': (9000., 2.5e6),
              'Rectifier diameter': (3., 24.),
              'Rectifier height': (27., 170.),
              'Rectifier weight': (9000., 2.5e6)}
    
    kwargs = {'P': 101325,
              'LHK': None,
              'y_top': None,
              'x_bot': None,
              'k': 1.25}
    
    @property
    def TS(self):
        return self._TS
    @TS.setter
    def TS(self, TS):
        """Tray spacing (225-600 mm)."""
        self._TS = TS
    
    @property
    def StageEfficiency(self):
        """Enforced user defined stage efficiency."""
        return self._StageEfficiency
    @StageEfficiency.setter
    def StageEfficiency(self, StageEfficiency):
        self._StageEfficiency = StageEfficiency
    
    @property
    def f(self):
        """Ratio of actual velocity to maximum velocity allowable before flooding."""
        return self._f
    @f.setter
    def f(self, f):
        self._f = f
    
    @property
    def F_F(self):
        """Foaming factor (0 < F_F < 1)."""
        return self._F_F
    @F_F.setter
    def F_F(self, F_F):
        if F_F > 1 or F_F < 0:
            raise ValueError(f"Foaming factor, 'F_F', must be between 0 and 1, ({F_F} given).")
        self._F_F = F_F
    
    @property
    def A_ha(self):
        """Ratio of open area, A_h, to active area, A_a."""
        return self._A_ha
    @A_ha.setter
    def A_ha(self, A_ha):
        self._A_ha = A_ha
    
    @property
    def A_dn(self):
        """Enforced ratio of downcomer area to net (total) area. If None, estimate ratio based on Oliver's estimation [1]."""
        return self._A_dn
    @A_dn.setter
    def A_dn(self, A_dn):
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
            raise ValueError(f"Tray type must be one of the following: {dummy}")
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
            raise ValueError(f"Tray material must be one of the following: {dummy}")
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
            raise ValueError(f"Vessel material must be one of the following: {dummy}")
        self._F_Mstr = vessel_material  


    def setup(self):
        vap, liq = self.outs
        kwargs = self.kwargs
        self._cached = cached = {}
        species = vap.species

        # Set stream phase and pressure
        vap.phase = 'g'
        liq.phase = 'l'
        vap.P = liq.P = kwargs['P']
        if not kwargs['LHK']:
            raise ValueError("Must specify light and heavy key, 'LHK'.")

        # Set light non-key and heavy non-key indeces
        LK, HK = kwargs['LHK']
        sp_index = vap._ID_index
        LK_index, HK_index = cached['LHK_index'] = [sp_index[LK], sp_index[HK]]
        IDs = list(vap._specie_IDs)
        indeces = list(range(len(IDs)))
        IDs.pop(LK_index)
        indeces.pop(LK_index)
        if HK_index > LK_index:
            HK_index -= 1
        IDs.pop(HK_index)
        indeces.pop(HK_index)
        Tbs = species._get_Tb(IDs)
        Tb_light = getattr(species, LK).Tb
        Tb_heavy = getattr(species, HK).Tb
        cached['LNK_index'] = LNK_index = []
        cached['HNK_index'] = HNK_index = []
        for Tb, i in zip(Tbs, indeces):
            if Tb < Tb_light:
                LNK_index.append(i)
            elif Tb > Tb_heavy:
                HNK_index.append(i)
            else:
                raise ValueError(f"Intermediate volatile specie, '{vap._specie_IDs[i]}', between light and heavy key, ['{LK}', '{HK}'].")
    
        self._update_composition_requirement(kwargs['y_top'], kwargs['x_bot'])
        self._setup_components()
        
    def _setup_components(self):
        self._condenser = HXutility(self.ID+' condenser', None)
        self._boiler = HXutility(self.ID+' boiler', None)
    
    def _update_composition_requirement(self, y_top, x_bot):
        # Set light and heavy key compositions at top and bottoms product
        cached = self._cached
        cached['y'] = array([y_top, 1-y_top])
        cached['x'] = array([x_bot, 1-x_bot])
    
    def _mass_balance(self):
        vap, liq = self.outs
        kwargs = self.kwargs
        cached = self._cached

        # Get all important flow rates (both light and heavy keys and non-keys)
        LHK_index, LNK_index, HNK_index, y = get_vals(cached, 'LHK_index',
                                                   'LNK_index', 'HNK_index', 'y')
        mol = self._mol_in
        LHK_mol = mol[LHK_index]
        HNK_mol = mol[HNK_index]
        LNK_mol = mol[LNK_index]

        # Set light and heavy keys by lever rule
        light, heavy = LHK_mol
        LHK_molnet = light + heavy
        zf = light/LHK_molnet
        y_top, x_bot = kwargs['y_top'], kwargs['x_bot']
        split_frac = (zf-x_bot)/(y_top-x_bot)
        top_net = LHK_molnet*split_frac

        # Set output streams
        vap.mol[LHK_index] = top_net * y
        liq.mol[LHK_index] = LHK_mol - vap.mol[LHK_index]
        vap.mol[LNK_index] = LNK_mol
        liq.mol[HNK_index] = HNK_mol
    
    def _run(self)
        # Run mass balance
        self._mass_balance()
        
        # Unpack arguments
        vap, liq = self.outs
        cached = self._cached
        
        cached['vle_top'], cached['top_index'] = vle_top, top_index = vap._equilibrium_species()
        cached['vle_bot'], cached['bot_index'] = vle_bot, bot_index = liq._equilibrium_species()

        # Get top and bottom compositions
        vap_mol = vap.mol[top_index]
        y = vap_mol/sum(vap_mol)

        liq_mol = liq.mol[bot_index]
        x = liq_mol/sum(liq_mol)

        # Run top equilibrium to find temperature and composition of condensate
        T, y = vap._dew_point(specie_IDs=vle_top, y=y, P=vap.P)
        cached['condensate_molfrac'] = y
        vap.T = T

        # Run bottoms equilibrium to find temperature
        T, y = liq._bubble_point(specie_IDs=vle_bot, x=x, P=liq.P)
        cached['boilup_molfrac'] = y
        liq.T = T

    def _equilibrium_staircase(self, operating_line, x_stairs,
                              y_stairs, T_stairs, x_limit):
        """Find the specifications at every stage of the of the operating line before the maximum liquid molar fraction. Append the light key liquid molar fraction, light key vapor molar fraction, and stage temperatures to x_stairs, y_stairs and T_stairs respectively.
        
        **Parameters**
        
            operating_line: [function] Should return the liquid molar fraction of the light key given its vapor molar fraction.
            
            x_stairs: [list] Liquid molar compositions at each stage. Last element should be the starting point for the next stage
            
            y_stairs: [list] Vapor molar compositions at each stage. Last element should be the starting point for the next stage
            
            T_stairs: [list] Temperature at each stage.
            
        """
        specie_IDs = self.kwargs['LHK']
        P = self.kwargs['P']
        bubble_point = self.outs[1]._bubble_point
        i = 0
        yi = y_stairs[-1]
        xi = x_stairs[-1]
        while xi < x_limit:
            if i > 100:
                print('Cannot meet specifications! Stages > 100')
                return
            i += 1
            # Go Up
            T_guess, y = bubble_point(specie_IDs, array((xi, 1-xi)), P)
            yi = y[0]
            y_stairs.append(yi)
            T_stairs.append(T_guess)
            # Go Right
            xi = operating_line(yi)
            if xi > 1:
                xi = 0.9999999
            x_stairs.append(xi)

    def _plot_stages(self):
        """Plot stages, graphical aid line, and equilibrium curve. The plot does not include operating lines nor a legend."""
        # Cached Data
        vap, liq = self.outs
        cached = self._cached
        x_stages = cached.get('x_stages')
        if not x_stages:
            self.operation()
            x_stages = cached.get('x_stages')
        y_stages = cached['y_stages']
        kwargs = self.kwargs
        specie_IDs = kwargs['LHK']
        light = specie_IDs[0]
        P = kwargs['P']
        
        # Equilibrium data
        x_eq = np.linspace(0, 1, 100)
        y_eq = np.zeros(100)
        T = np.zeros(100)
        n = 0
        
        bubble_point = vap._bubble_point
        for xi in x_eq:
            T[n], y = bubble_point(specie_IDs, array([xi, 1-xi]), P)
            y_eq[n] = y[0]
            n += 1
            
        # Set-up graph
        plt.figure()
        plt.xticks(np.arange(0, 1.1, 0.1), fontsize=12)
        plt.yticks(fontsize=12)
        plt.xlabel('x (' + light + ')', fontsize=16)
        plt.ylabel('y (' + light + ')', fontsize=16)
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

    def _cost_trays(self, N_T, Di:'ft'):
        """Return total cost of all trays."""
        # Note: Can only use this function after running design method.
        C_BT = self._calc_TrayBaseCost(Di, self.CEPCI)
        F_NT = self._calc_NTrayFactor(N_T)
        return N_T*F_NT*self._F_TT*self._F_TM(Di)*C_BT

    def _cost_tower(self, Di:'ft', L:'ft', W:'lb'):
        C_V = self._calc_EmptyTowerCost(W)
        C_PL = self._calc_PlaformLadderCost(Di, L)
        return (self._F_M*C_V + C_PL)*self.CEPCI/500
    
    @staticmethod
    def _calc_EmptyTowerCost(W:'lb')->'C_V (USD)':
        """Return the cost of an empty tower vessel assuming a CE of 500.
        
        **Parameters**
        
            W: Weight (lb)
        
        """
        return np.exp(7.2756 + 0.18255*np.log(W) + 0.02297*np.log(W)**2)
    
    @staticmethod
    def _calc_PlaformLadderCost(Di:'ft', L:'ft')-> 'C_PL (USD)':
        """Return the cost of platforms and ladders assuming a CE of 500.
        
        **Parameters**
        
            Di: Inner diameter (ft)
            L: Legnth (ft)
        
        """
        return 300.9*Di**0.63316*L**0.80161
    
    @staticmethod
    def _calc_Weight(Di:'ft', L:'ft', tv:'in', rho_M:'lb/in3') -> 'W (lb)':
        """Return the weight of the tower assuming 2:1 elliptical head.
        
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
    def _calc_WallThickness(Po:'psi', Di:'ft', L:'ft', S:'psi'=15000, E=None) -> 'tv (in)':
        """Return the wall thinkness designed to withstand the internal pressure and the wind/earthquake load at the bottom.
        
        **Parameters**
        
            Po: Operating internal pressure (psi)
            Di: Internal diameter (ft)
            L: Height (ft)
            S: Maximum stress (psi)
            E: Fractional weld efficiency
            
        """
        Di = Di*12 # ft to in
        L = L*12
        
        # Get design pressure, which should be higher than operating pressure.
        Po_gauge = Po - 14.69
        if Po_gauge < 0:
            raise DesignError('Vacuum vessels are not implemented yet.')
        elif Po_gauge < 5:
            Pd = 10
        elif Po_gauge < 1000:
            Pd = np.exp(0.60608 + 0.91615*np.log(Po)) + 0.0015655*np.log(Po)**2
        else:
            Pd = 1.1*Po_gauge
        
        # Calculate thinkess according to ASME pressure-vessel code.
        E_check = E is None
        if E_check:
            # Assume carbon steel with thickness more than 1.25 in
            E = 1.0 
        ts = Pd*Di/(2*S*E-1.2*Pd)
        
        if E_check:
            # Weld efficiency of 0.85 for low thickness carbon steel
            if ts < 1.25:
                E = 0.85
                ts = Pd*Di/(2*S*E-1.2*Pd)
        
        # Minimum thickness for vessel rigidity may be larger
        ts_min = np.polyval(ts_min_p, Di/12)
        if ts < ts_min:
            ts = ts_min
        
        # Calculate thickness to withstand wind/earthquake load
        Do = Di + ts
        tw = 0.22*(Do + 18)*L**2/(S*Do**2)
        tv = max(tw, ts)
        
        # Add corrosion allowence
        tv += 1/8
        # Vessels are fabricated from metal plates with small increments
        if tv < 0.5:
            tv = approx2step(tv, 3/16, 1/16)
        elif tv < 2:
            tv = approx2step(tv, 0.5, 1/8)
        elif tv < 3:
            tv = approx2step(tv, 2, 1/4)
        return tv
    
    @staticmethod
    def _calc_TrayBaseCost(Di:'ft', CE) -> 'C_BT':
        """Return base cost of a tray.
        
        **Parameters**
        
            Di: Inner diameter (ft)
            CE: Chemical Engineering Plant Cost Index
        """
        return CE * 0.825397 * np.exp(0.1482*Di)

    @staticmethod
    def _calc_NTrayFactor(N_T) -> 'F_NT':
        """Return cost factor for number of trays.
        
        **Parameters**
        
            N_T: Number of trays
            
        """
        if N_T < 20:
            F_NT = 2.25/1.0414**N_T
        else:
            F_NT = 1
        return F_NT

    @staticmethod
    def _calc_MurphreeEfficiency(mu: 'mPa*s', alpha, L, V) -> 'E_mv':
        """Return the sectional murphree efficiency.
        
        **Parameters**
            
            mu: Viscosity (mPa*s)
            
            alpha: Relative volatility
            
            L: Liquid flow rate by mol
            
            V: Vapor flow rate by mol
        
        """
        S = alpha*V/L # Stripping factor
        if S > 1:
            sigma = S
        else:
            sigma = 1/S
        return 0.503*mu**(-0.226)*sigma**(-0.08 )
    
    @staticmethod
    def _calc_FlowParameter(L, V, rho_V, rho_L) -> 'F_LV':
        """Return the flow parameter.
        
        **Parameters**
        
            L: Liquid flow rate by mass
            V: Vapor flow rate by mass
            rho_V: Vapor density
            rho_L: Liquid density
        
        """
        return L/V*(rho_V/rho_L)**0.5
    
    @staticmethod
    def _calc_MaxCapacityParameter(TS: 'mm', F_LV) -> ' C_sbf':
        """Return the maximum capacity parameter before flooding (m/s).
        
        **Parameters**
        
            TS: Tray spacing (mm)
            F_LV: Flow parameter
        
        """
        return 0.0105 + 8.127e-4*TS**0.755*np.exp(-1.463*F_LV**0.842)
    
    @staticmethod
    def _calc_MaxVaporVelocity(C_sbf: 'm/s', sigma: 'dyn/cm',
                               rho_L, rho_V, F_F, A_ha) -> 'U_f':
        """Return the maximum allowable vapor velocity through the net area of flow before flooding (m/s).
        
        **Parameters**
        
            C_sbf: Maximum Capacity Parameter (m/s)
            sigma: Liquid surface tension (dyn/cm)
            rho_L: Liquid density
            rho_V: Vapor density
            F_F: Foaming factor
            A_ha: Ratio of open area, A_h, to active area, A_a
        
        """
        F_ST = (sigma/20)**0.2 # Surface tension factor
        
        # Working area factor
        if A_ha >= 0.1 and A_ha <= 1:
            F_HA = 1
        elif A_ha >= 0.06:
            F_HA = 5*A_ha + 0.5
        else:
            raise ValueError(f"Ratio of open to active area, 'A', must be between 0.06 and 1 ({A_ha} given).") 
        
        return C_sbf * F_HA * F_ST * ((rho_L-rho_V)/rho_V)**0.5
    
    @staticmethod
    def _calc_DowncomerAreaRatio(F_LV) -> 'A_dn':
        """Ratio of downcomer area to net (total) area.
        
        **Parameters**
        
            F_LV: Flow parameter
        
        """
        if F_LV < 0.1:
            A_dn = 0.1
        elif F_LV < 1:
            A_dn = 0.1 + (F_LV-0.1)/9
        else:
            A_dn = 0.2
        return A_dn
    
    @staticmethod
    def _calc_Diameter(V_vol: 'm3/s', U_f: 'm/s', f, A_dn) -> 'D_T (ft)':
        """Return column diameter.
        
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
    def _calc_Height(TS: 'mm', Nstages: int, top=True, bot=True) -> 'H (ft)':
        """Return the height of the column (ft).
        
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
        condensate = self._cached['condensate'] # Abstract instance
        condenser = self._condenser
        s_in = Stream(ID='*', flow=distillate.mol+condensate.mol, phase='g', T=distillate.T, P=distillate.P)
        ms1 = MixedStream(ID='*', liquid_flow=condensate.mol, 
                          T=condensate.T, P=condensate.P)
        ms1.vapor_mol = distillate.mol
        s_in-condenser-ms1
        condenser.operation()
        condenser.design()
        condenser.cost()
        
    def _calc_boiler(self):
        bottoms = self.outs[1]
        boil_up = self._cached['boil_up'] # Abstract instance
        boiler = self._boiler
        s_in = Stream(ID='*', flow=bottoms.mol+boil_up.mol, phase='g',
                      T=bottoms.T, P=bottoms.P)
        ms1 = MixedStream(ID='*', vapor_flow=boil_up.mol,
                          T=boil_up.T, P=boil_up.P)
        ms1.liquid_mol = bottoms.mol
        s_in-boiler-ms1
        boiler.operation()
        boiler.design()
        boiler.cost()
        
    def cost(self):
        """
        * 'Trays': (USD)
        * 'Tower': (USD)
        * 'Boiler': (USD)
        """
        results = self.results
        Design = results['Design']
        Cost = results['Cost']
        
        # Cost trays assuming a partial condenser
        N_T = Design['Actual stages'] - 1
        Di = Design['Diameter']
        Cost['Trays'] = self._cost_trays(N_T, Di)
        
        # Cost vessel assuming T < 800 F
        W = Design['Weight'] # in lb
        L = Design['Height']*3.28 # in ft
        Cost['Tower'] = self._cost_tower(Di, L, W)
        
        self._cost_components(Cost)
        return Cost

class Distillation(Dist):
    __doc__ = column_doc.replace('{Column Type}', 'Distillation')
    _N_heat_util = 2
    is_divided = False #: [bool] True if the stripper and rectifier are two separate columns.    
    
    def operation(self):
        """
        * 'Theoretical feed stage': (#)
        * 'Theoretical stages': (#)
        * 'Minimum reflux': (Ratio)
        * 'Reflux': (Ratio)
        * 'T_condenser': (K) 
        * 'T_boiler': (K)
        * 'Cooling': (kJ/hr)
        * 'Heating': (kJ/hr)
        """
        vap, liq = self.outs
        kwargs = self.kwargs
        cached = self._cached
        bubble_point = liq._bubble_point

        # Some important info
        LHK_index = cached['LHK_index']
        specie_IDs = kwargs['LHK']
        P = kwargs['P']
        k = kwargs['k']
        y_top, x_bot = kwargs['y_top'], kwargs['x_bot']

        # Feed light key mol fraction
        Nspecies = liq._Nspecies
        liq_mol = np.zeros(Nspecies)
        vap_mol = np.zeros(Nspecies)
        for s in self.ins:
            if s.phase == 'g':
                vap_mol += s.mol
            elif s.phase.lower() == 'l':
                liq_mol += s.mol
            elif s.phase == 's':
                pass
            elif isinstance(s, MixedStream):
                liq_mol += s.liquid_mol
                vap_mol += s.vapor_mol
            else:
                raise biosteam_Error(f'Invalid phase encountered in stream {s.ID}')
        cached['feed_liqmol'] = liq_mol
        cached['feed_vapmol'] = vap_mol
        LHK_mol = liq_mol[LHK_index] + vap_mol[LHK_index]
        LHK_molnet = sum(LHK_mol)
        zf = LHK_mol[0]/LHK_molnet
        
        # Get feed quality
        q = sum(liq_mol[LHK_index])/LHK_molnet
        
        # Get R_min and the q_line 
        if q == 1:
            q = 1 - 1e-5
        self._q_line = q_line = lambda x: q*x/(q-1) - zf/(q-1)
        
        y = None
        def Rmin_intersection(x):
            nonlocal specie_IDs, P, y
            T_guess, y = bubble_point(specie_IDs, array((x, 1-x)), P)
            return q_line(x) - y[0]
        
        x_Rmin = brentq(Rmin_intersection, 0, 1)
        y_Rmin = y[0]
        m = (y_Rmin-y_top)/(x_Rmin-y_top)
        Rmin = m/(1-m)
        if Rmin <= 0:
            R = 0.1*k
        else:
            R = Rmin*k

        # Rectifying section: Inntersects q_line with slope given by R/(R+1)
        m1 = R/(R+1)
        b1 = y_top-m1*y_top
        def rs(y) -> 'x': return (y - b1)/m1
        def intersect1(y) -> 'error': return y - q_line(rs(y))
        self._y_m = y_m = brentq(intersect1, 0, 1)
        self._x_m = x_m = rs(y_m)
        
        # Stripping section: Intersects Rectifying section and q_line and beggins at bottoms liquid composition
        m2 = (x_bot-y_m)/(x_bot-x_m)
        b2 = y_m-m2*x_m
        def ss(y) -> 'x': return (y - b2)/m2
                
        # Data for staircase
        cached['x_stages'] = x_stages = [x_bot]
        cached['y_stages'] = y_stages = [x_bot]
        cached['T_stages'] = T_stages = []
        self._equilibrium_staircase(ss, x_stages, y_stages, T_stages, x_m)
        yi = y_stages[-1]
        x_stages[-1] = rs(yi)
        self._equilibrium_staircase(rs, x_stages, y_stages, T_stages, y_top)
        
        # Find feed stage
        for i in range(len(y_stages)-1):
            j = i+1
            if y_stages[i] < y_m and y_stages[j] > y_m:
                feed_stage = i+1
        stages = len(x_stages)
        
        # Set results
        Operation = self.results['Operation']
        Lr = cached['condensate_molfrac']*(vap.molnet*R)
        Hvapm = vap._species.get_props(cached['vle_top'], 'Hvapm')
        
        # Condenser Duty
        cooling = -sum(Hvapm*Lr)
        
        # Reboiler Duty
        heating = (vap.H + liq.H) - sum(s.H for s in self.ins) - cooling
        
        # Calculate Utilities
        self.heat_utilities[0](cooling, vap.T, vap.T)
        self.heat_utilities[1](heating, liq.T, liq.T)
        
        # Results
        Operation['Theoretical feed stage'] = feed_stage
        Operation['Theoretical stages'] = stages
        Operation['Minimum reflux'] = Rmin
        Operation['Reflux'] = R
        Operation['T_condenser'] = vap.T# K
        Operation['T_boiler'] = liq.T   # K
        Operation['Cooling'] = cooling  # kJ/hr
        Operation['Heating'] = heating  # kJ/hr
        return Operation  
    
    def _calc_Nstages(self) -> 'Nstages':
        """Return a tuple with the actual number of stages for the rectifier and the stripper."""
        vap, liq = self.outs
        kwargs = self.kwargs
        cached = self._cached
        results = self.results
        Operation = results['Operation']
        x_stages = cached['x_stages']
        y_stages = cached['y_stages']
        R = Operation['Reflux']
        N_stages = Operation['Theoretical stages']
        feed_stage = Operation['Theoretical feed stage']
        liq_mol = cached['feed_liqmol']
        vap_mol = cached['feed_vapmol']
        
        stage_efficiency = kwargs.get('StageEfficiency')
        if stage_efficiency:
            return N_stages/stage_efficiency
        else:    
            # Calculate Murphree Efficiency for rectifying section
            vap_molnet = vap.molnet
            condensate_molfrac = cached['condensate_molfrac']
            vle_top = cached['vle_top']
            cached['L_Rmol'] = L_Rmol = R*vap_molnet
            cached['V_Rmol'] = V_Rmol = (R+1)*vap_molnet
            cached['condensate'] = condensate = type(liq)(flow=condensate_molfrac,
                                                          species=vle_top, T=vap.T, P=vap.P)
            condensate.mol *= L_Rmol
            mu = condensate.prop_quantity('mu').to('mPa*s').magnitude
            K_light = y_stages[-1]/x_stages[-1] 
            K_heavy = (1-y_stages[-1])/(1-x_stages[-1] )
            alpha = K_light/K_heavy
            cached['Rectifying Section Efficiency'] = E_rectifier = self._calc_MurphreeEfficiency(mu, alpha, L_Rmol, V_Rmol)
            
            # Calculate Murphree Efficiency for stripping section
            mu = liq.prop_quantity('mu').to('mPa*s').magnitude
            cached['V_Smol'] = V_Smol = (R+1)*vap_molnet - sum(vap_mol)
            cached['L_Smol'] = L_Smol = R*vap_molnet + sum(liq_mol) 
            K_light = y_stages[0]/x_stages[0] 
            K_heavy = (1-y_stages[0])/(1-x_stages[0] )
            alpha = K_light/K_heavy
            cached['Stripping Section Efficiency'] = E_stripper = self._calc_MurphreeEfficiency(mu, alpha, L_Smol, V_Smol)
            
            # Calculate actual number of stages
            mid_stage = feed_stage - 0.5
            N_rectifier = np.ceil(mid_stage/E_rectifier)
            N_stripper = np.ceil( (N_stages-mid_stage)/E_stripper )
        return N_rectifier, N_stripper

    def design(self):
        """
        Divided column:
            * 'Rectifier stages': (#)
            * 'Rectifier height': (ft)
            * 'Rectifier diameter': (ft)
            * 'Rectifier wall thickness': (in)
            * 'Rectifier weight': (lb)
            * 'Stripper stages': (#)
            * 'Stripper height': (ft) 
            * 'Stripper diameter': (ft)
            * 'Stripper wall thickness': (in)
            * 'Stripper weight': (lb)
        One column:
            * 'Rectifier stages': (#)
            * 'Stripper stages': (#)
            * 'Actual stages': (#)
            * 'Height': (ft)
            * 'Diameter': (ft)
            * 'Wall thickness': (in)
            * 'Weight': (lb)
        """
        distillate, bottoms = self.outs
        kwargs = self.kwargs
        results = self.results
        cached = self._cached
        Design = results['Design']
        Operation = results['Operation']
        Rstages, Sstages = self._calc_Nstages()
        calc_Height = self._calc_Height
        is_divided = self.is_divided
        TS = self.TS
        Design['Rectifier stages'], Design['Stripper stages'] = Rstages, Sstages
        Design['Actual stages'] = Rstages + Sstages
        
        ### Get diameter of rectifying section based on top plate ###
        
        R = Operation['Reflux']
        condensate = cached['condensate']
        rho_L = condensate.rho
        sigma = condensate.prop_quantity('sigma').to('dyn/cm').magnitude
        L = condensate.massnet
        V = L*(R+1)/R
        vapor_stream = type(distillate).like(distillate)
        vapor_stream.mol *= R+1
        V_vol = vapor_stream.prop_quantity('volnet').to('m^3/s').magnitude
        rho_V = distillate.rho
        F_LV = self._calc_FlowParameter(L, V, rho_V, rho_L)
        C_sbf = self._calc_MaxCapacityParameter(TS, F_LV)
        F_F = self.F_F
        A_ha = self.A_ha
        U_f = self._calc_MaxVaporVelocity(C_sbf, sigma, rho_L, rho_V, F_F, A_ha)
        A_dn = self.A_dn
        if A_dn is None:
           self.A_dn = A_dn = self._calc_DowncomerAreaRatio(F_LV)
        f = self.f
        R_diameter = self._calc_Diameter(V_vol, U_f, f, A_dn)
        
        ### Get diameter of stripping section based on feed plate ###
        
        V_mol = cached['V_Smol']
        L_mol = cached['L_Smol']
        rho_L = bottoms.rho
        boil_up_flow = cached['boilup_molfrac'] * V_mol
        cached['boil_up'] = boil_up = type(distillate)(flow=boil_up_flow,
              species=cached['vle_bot'], T=bottoms.T, P=bottoms.P, phase='g')
        V = boil_up.massnet
        V_vol = boil_up.prop_quantity('volnet').to('m^3/s').magnitude
        rho_V = boil_up.rho
        L = sum(bottoms._MW*bottoms.molfrac*L_mol) # To get liquid going down
        F_LV = self._calc_FlowParameter(L, V, rho_V, rho_L)
        C_sbf = self._calc_MaxCapacityParameter(TS, F_LV)
        sigma = bottoms.prop_quantity('sigma').to('dyn/cm').magnitude
        F_F = self.F_F
        A_ha = self.A_ha
        U_f = self._calc_MaxVaporVelocity(C_sbf, sigma, rho_L, rho_V, F_F, A_ha)
        A_dn = self.A_dn
        if A_dn is None:
            self.A_dn = self._calc_DowncomerAreaRatio(F_LV)
        f = self.f
        S_diameter = self._calc_Diameter(V_vol, U_f, f, A_dn)
        Po = kwargs['P']/101325*14.7
        rho_M = rho_Mdict[self._F_Mstr]
        
        if is_divided:
            Design['Rectifier height'] = H_R = calc_Height(TS, Rstages-1)
            Design['Stripper height'] = H_S = calc_Height(TS, Sstages-1)
            Design['Rectifier diameter'] = R_diameter
            Design['Stripper diameter'] = S_diameter
            Design['Rectifier wall thickness'] = tv = self._calc_WallThickness(Po, R_diameter, H_R)
            Design['Stripper wall thickness'] = tv = self._calc_WallThickness(Po, S_diameter, H_S)
            Design['Rectifier weight'] = self._calc_Weight(R_diameter, H_R, tv, rho_M)
            Design['Stripper weight'] = self._calc_Weight(S_diameter, H_S, tv, rho_M)
        else:
            Design['Height'] = H = calc_Height(TS, Rstages+Sstages-2)
            Design['Diameter'] = Di = max((R_diameter, S_diameter))
            Design['Wall thickness'] = tv = self._calc_WallThickness(Po, Di, H)
            Design['Weight'] = self._calc_Weight(Di, H, tv, rho_M)
        return Design
        
    def cost(self):
        """
        Divided column:
            * 'Stripper trays': (USD)
            * 'Stripper tower': (USD)
            * 'Rectifier trays': (USD)
            * 'Rectifier tower': (USD)
            * 'Condenser': (USD)
            * 'Boiler': (USD)
        One column:
            * 'Trays': (USD)
            * 'Tower': (USD)
            * 'Condenser': (USD)
            * 'Boiler': (USD)
        """
        if not self.is_divided:
            return super().cost()
        
        Design = self.results['Design']
        Cost = self.results['Cost']
        
        # Number of trays assuming a partial condenser
        N_RT = Design['Rectifier stages'] - 1
        Di_R = Design['Rectifier diameter']
        Cost['Rectifier trays'] = self._cost_trays(N_RT, Di_R)
        N_ST = Design['Stripper stages'] - 1
        Di_S = Design['Stripper diameter']
        Cost['Stripper trays'] = self._cost_trays(N_ST, Di_S)
        
        # Cost vessel assuming T < 800 F
        W_R = Design['Rectifier weight'] # in lb
        H_R = Design['Rectifier height']*3.28 # in ft
        Cost['Rectifier tower'] = self._cost_tower(Di_R, H_R, W_R)
        W_S = Design['Stripper weight'] # in lb
        H_S = Design['Stripper height']*3.28 # in ft
        Cost['Stripper tower'] = self._cost_tower(Di_S, H_S, W_S)
        self._cost_components(Cost)
        return Cost
    
    def _cost_components(self, Cost): 
        # Cost condenser
        self._calc_condenser()
        Cost['Condenser'] = self._condenser.results['Cost']['Heat exchanger']
        
        # Cost boiler
        self._calc_boiler()
        Cost['Boiler'] = self._boiler.results['Cost']['Heat exchanger']
        
    
    def plot_stages(self):
        """Plot the McCabe Thiele Diagram."""
        # Plot stages, graphical aid and equilibrium curve
        self._plot_stages()
        
        # Cached Data
        vap, liq = self.outs
        results = self.results['Operation']
        cached = self._cached
        x_stages = cached.get('x_stages')
        if not x_stages:
            self.operation()
            x_stages = cached.get('x_stages')
        kwargs = self.kwargs
        q_line = self._q_line
        y_top = kwargs['y_top']
        x_bot = kwargs['x_bot']
        stages = results['Theoretical stages']
        Rmin = results['Minimum reflux']
        R = results['Reflux']
        feed_stage = results['Theoretical feed stage']
        
        # q_line
        def intersect2(x): return x - q_line(x)
        x_m2 = brentq(intersect2, 0, 1)
        
        # Graph q-line, Rectifying and Stripping section
        plt.plot([self._x_m, x_m2], [self._y_m, x_m2])
        plt.plot([self._x_m, y_top], [self._y_m, y_top])
        plt.plot([x_bot, self._x_m], [x_bot, self._y_m])
        plt.legend([f'Stages: {stages}, Feed: {feed_stage}', 'Graphical aid', 'eq-line', 'q-line', 'ROL', 'SOL'], fontsize=12)
        plt.title(f'McCabe Thiele Diagram (Rmin = {Rmin:.2f}, R = {R:.2f})')
        plt.show()
        return plt


class Stripper(Dist):
    __doc__ = column_doc.replace('{Column Type}', 'Stripper')
    _N_heat_util = 1    
    
    def _setup_components(self):
        self._boiler = HXutility(self.ID+' boiler', None)
    
    def operation(self):
        """
        * 'Theoretical feed stage': (#)
        * 'Theoretical stages': (#)
        * 'Minimum boil up': (Ratio)
        * 'Boil up': (Ratio)
        * 'T_boiler': (K) 
        * 'Heating': (kJ/hr)
        """
        vap, liq = self.outs
        kwargs = self.kwargs
        cached = self._cached
        
        # Some important info
        specie_IDs = kwargs['LHK']
        P = kwargs['P']
        k = kwargs['k']
        y_top, x_bot = kwargs['y_top'], kwargs['x_bot']
        # Get B_min (Boil up ratio)
        y_Rmin = y_top
        T_guess, x = vap._dew_point(specie_IDs, [y_Rmin, 1-y_Rmin], P)
        x_Rmin = x[0]
        m = (y_Rmin-x_bot)/(x_Rmin-x_bot)
        Bmin = 1/(m-1)
        if Bmin <= 0:
            B = 0.1*k
        else:
            B = Bmin*k

        # Stripping section: Inntersects liquid composition with slope given by (B+1)/B
        m = (B+1)/B
        b = x_bot-m*x_bot
        def ss(y): return (y - b)/m
        self._ss = ss
                
        # Data for staircase
        cached['x_stages'] = x_stages = [x_bot]
        cached['y_stages'] = y_stages = [x_bot]
        cached['T_stages'] = T_stages = []
        self._equilibrium_staircase(ss, x_stages, y_stages, T_stages, ss(y_top))
        stages = len(x_stages)
        
        # Set results
        Operation = self.results['Operation']
        Vs = cached['boilup_molfrac']*(liq.molnet*B)
        Hvapm = vap._species.get_props(cached['vle_bot'], 'Hvapm')
        
        # Reboiler Duty
        heating = sum(Hvapm*Vs)
        
        # Calculate Utilities
        self.heat_utilities[0](heating, liq.T, liq.T)
        
        # Results
        Operation['Theoretical stages'] = stages
        Operation['Minimum boil up'] = Bmin
        Operation['Boil up'] = B
        Operation['T_boiler'] = liq.T # K
        Operation['Heating'] = heating   # kJ/hr
        return Operation  
    
    def plot_stages(self):
        # Plot stages, graphical aid and equilibrium curve
        self._plot_stages()
        
        # Cached Data
        vap, liq = self.outs
        Operation = self.results['Operation']
        cached = self._cached
        x_stages = cached.get('x_stages')
        if not x_stages:
            self.operation()
            x_stages = cached.get('x_stages')
        kwargs = self.kwargs
        y_top = kwargs['y_top']
        x_bot = kwargs['x_bot']
        stages = Operation['Theoretical stages']
        Bmin = Operation['Minimum boil up']
        B = Operation['Boil up']
        
        # Graph Stripping section
        ss = self._ss
        plt.plot([x_bot, ss(y_top)], [x_bot, y_top])
        plt.legend([f'Stages: {stages}', 'Graphical aid', 'eq-line', 'SOL'], fontsize=12)
            
        # Title
        plt.title(f'McCabe Thiele Diagram (Bmin = {Bmin:.2f}, B = {B:.2f})')
        plt.show()
        return plt
    
    def _calc_Nstages(self):
        """Return the actunal number of stages"""
        vap, liq = self.outs
        kwargs = self.kwargs
        cached = self._cached
        results = self.results
        Operation = results['Operation']
        x_stages = cached['x_stages']
        y_stages = cached['y_stages']
        LHK_index = cached['LHK_index']
        B = Operation['Boil up']
        N_stages = Operation['Theoretical stages']
        
        stage_efficiency = kwargs.get('StageEfficiency')
        if stage_efficiency:
            return N_stages/stage_efficiency
        else:
            # Calculate Murphree Efficiency for stripping section
            mu = liq.prop_quantity('mu').to('mPa*s').magnitude
            cached['V_mol'] = V_mol = B*sum(liq.mol[LHK_index])
            cached['L_mol'] = L_mol = liq.molnet + V_mol
            K_light = y_stages[0]/x_stages[0] 
            K_heavy = (1-y_stages[0])/(1-x_stages[0] )
            alpha = K_light/K_heavy
            cached['Stripping Section Efficiency'] = E_stripper = self._calc_MurphreeEfficiency(mu, alpha, L_mol, V_mol)
        return np.ceil(N_stages/E_stripper)
    
    def design(self):
        """
            * 'Actual stages': (#)
            * 'Height': (ft)
            * 'Diameter': (ft)
            * 'Wall thickness': (in)
            * 'Weight': (lb)
        """
        distillate, bottoms = self.outs
        kwargs = self.kwargs
        results = self.results
        cached = self._cached
        Design = results['Design']
        
        ### Get number of stages and height ###
        
        TS = self.TS
        Design['Actual stages'] = Nstages = self._calc_Nstages()
        Design['Height'] = H = self._calc_Height(TS, Nstages-1)
        
        ### Get diameter of stripping section based on feed plate ###
        
        V_mol = cached['V_mol']
        L_mol = cached['L_mol']
        rho_L = bottoms.rho
        boil_up_flow = cached['boilup_molfrac'] * V_mol
        cached['boil_up'] = boil_up = type(distillate)(flow=boil_up_flow,
              species=cached['vle_bot'], T=bottoms.T, P=bottoms.P, phase='g')
        V = boil_up.massnet
        V_vol = boil_up.prop_quantity('volnet').to('m^3/s').magnitude
        rho_V = boil_up.rho
        L = sum(bottoms._MW*bottoms.molfrac*L_mol) # To get liquid going down
        F_LV = self._calc_FlowParameter(L, V, rho_V, rho_L)
        C_sbf = self._calc_MaxCapacityParameter(TS, F_LV)
        sigma = bottoms.prop_quantity('sigma').to('dyn/cm').magnitude
        F_F = self.F_F
        A_ha = self.A_ha
        U_f = self._calc_MaxVaporVelocity(C_sbf, sigma, rho_L, rho_V, F_F, A_ha)
        A_dn = self.A_dn
        if A_dn is None:
            A_dn = self._calc_DowncomerAreaRatio(F_LV)
        f = self.f
        Design['Diameter'] = Di = self._calc_Diameter(V_vol, U_f, f, A_dn)
        Po = kwargs['P']/101325*14.7
        Design['Wall thickness'] = tv = self._calc_WallThickness(Po, Di, H)
        rho_M = rho_Mdict[self._F_Mstr]
        Design['Weight'] = self._calc_Weight(Di, H, tv, rho_M)
        return Design
    
    def _cost_components(self, Cost):
        # Cost boiler
        self._calc_boiler()
        Cost['Boiler'] = self._boiler.results['Cost']['Heat exchanger']

    