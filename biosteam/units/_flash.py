# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from .. import Unit, PowerUtility
from thermosteam import MultiStream, separations
from math import pi
import numpy as np
from . import design_tools as design
from .splitting import Splitter
from .heat_exchange import HX, HXutility
from .._graphics import vertical_vessel_graphics
import biosteam as bst

exp = np.exp
ln = np.log

# USDA Biodiesel: 
#V = np.pi*D**2/4*Ht*0.3048**3
#if not (0.1 < V < 70):
#    raise DesignError(f"Volume is out of bounds for costing")
#lambda V, CE: CE*13*V**0.62 # V (m3)
__all__ = ('Flash', 'SplitFlash', 'RatioFlash')


# %% Flash

class Flash(design.PressureVessel, Unit):
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
    vessel_type=None : 'Horizontal' or 'Vertical', optional
        Vessel separation type. If not specified, the vessel type will be chosen
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
    >>> settings.set_thermo(['Water', 'Glycerol'], cache=True)
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
    Flash                                                   Units            F1
    Medium pressure steam Duty                              kJ/hr      4.82e+07
                          Flow                            kmol/hr      1.33e+03
                          Cost                             USD/hr           367
    Design                Vessel type                                  Vertical
                          Length                               ft          15.5
                          Diameter                             ft           8.5
                          Weight                               lb      9.57e+03
                          Wall thickness                       in         0.438
                          Vessel material                          Carbon steel
    Purchase cost         Heat_exchanger - Floating head      USD      4.24e+04
                          Vertical pressure vessel            USD      4.63e+04
                          Platform and ladders                USD      1.39e+04
    Total purchase cost                                       USD      1.03e+05
    Utility cost                                           USD/hr           367


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
    auxiliary_unit_names = ('heat_exchanger',)
    _units = {'Length': 'ft',
              'Diameter': 'ft',
              'Weight': 'lb',
              'Wall thickness': 'in'}
    _max_agile_design = (
        'Length',
        'Diameter',
        'Weight',
        'Wall thickness',
    )
    _F_BM_default = {'Liquid-ring pump': 1.0,
                     **design.PressureVessel._F_BM_default}
    _graphics = vertical_vessel_graphics 
    _N_outs = 2
    _N_heat_utilities = 0

    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 V=None, T=None, Q=None, P=None, y=None, x=None,
                 vessel_material='Carbon steel',
                 vacuum_system_preference='Liquid-ring pump',
                 has_glycol_groups=False,
                 has_amine_groups=False,
                 vessel_type=None,
                 holdup_time=15,
                 surge_time=7.5,
                 has_mist_eliminator=False):
        Unit.__init__(self, ID, ins, outs, thermo)
        self._load_components()
        
        #: Enforced molar vapor fraction
        self.V = V
        
        #: Enforced operating temperature (K)
        self.T = T
        
        #: [array_like] Molar composition of vapor (for binary mixture)
        self.y = y
        
        #: [array_like] Molar composition of liquid (for binary mixture)
        self.x = x
        
        #: Enforced duty (kJ/hr)
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
        
    def _load_components(self):
        self._multi_stream = ms = MultiStream(None, thermo=self.thermo)
        self.heat_exchanger = hx = HXutility(None, (None,), ms, thermo=self.thermo) 
        hx.owner = self.owner
        self.heat_utilities = (*hx.heat_utilities, bst.HeatUtility(), bst.HeatUtility())
        
    def reset_cache(self):
        self._multi_stream.reset_cache()
        self.heat_exchanger.reset_cache()
        
    @property
    def P(self):
        """Operating pressure (Pa)."""
        return self._P
    @P.setter
    def P(self, P):
        if P and P < 101325 and not self.power_utility:
            self.power_utility = PowerUtility()
        self._P = P

    def _default_vessel_type(self):
        vap, liq = self.outs
        F_mass_vap = vap.F_mass
        F_mass_liq = liq.F_mass 
        assert F_mass_liq, "no liquid effluent; flash vessel must have liquid to default vessel type"
        return 'Vertical'if F_mass_vap / F_mass_liq > 0.1 else 'Horizontal'

    def _run(self):
        separations.vle(self.ins[0], *self.outs, self.T, self.P, self.V, 
                        self.Q, self.x, self.y, self._multi_stream)
            
    def _design(self):
        vessel_type = self.vessel_type
        if vessel_type == 'Vertical': 
            args = self._vertical_vessel_pressure_diameter_and_length()
        elif vessel_type == 'Horizontal': 
            args = self._horizontal_vessel_pressure_diameter_and_length()
        else: raise RuntimeError('unknown vessel type') # pragma: no cover
        feed = self.heat_exchanger.ins[0]
        feed.mix_from(self.ins)
        feed.vle(P=self.outs[0].P, H=feed.H)
        if self.Q == 0:
            self.heat_exchanger.baseline_purchase_costs.clear()
            self.heat_exchanger.purchase_costs.clear()
            self.heat_exchanger.installed_costs.clear()
        else:
            self.heat_exchanger._summary()
        self.design_results.update(
            self._vessel_design(*args)
        )

    def _cost(self):
        D = self.design_results
        self.baseline_purchase_costs.update(
            self._vessel_purchase_cost(D['Weight'], D['Diameter'], D['Length'])
        )
        self._cost_vacuum()

    def _cost_vacuum(self):
        P = self.P
        if not P or P > 101320: return 
        
        Design = self.design_results
        R = Design['Diameter'] * 0.5
        volume = 0.02832 * np.pi * Design['Length'] * R * R
        
        # If vapor is condensed, assume vacuum system is after condenser
        vapor = self.outs[0]
        hx = vapor.sink
        if isinstance(hx, HX):
            index = hx.ins.index(vapor)
            stream = hx.outs[index]
            if isinstance(stream, MultiStream): # pragma: no cover
                vapor = stream['g']
                F_mass = vapor.F_mass
                F_vol = vapor.F_vol
            else:
                if stream.phase == 'g': # pragma: no cover
                    F_mass = stream.F_mass
                    F_vol = stream.F_vol
                else:
                    F_mass = 0
                    F_vol = 0
        else:
            F_mass = vapor.F_mass
            F_vol = vapor.F_vol
        
        vacuum_results = design.compute_vacuum_system_power_and_cost(
            F_mass, F_vol, P, volume, self.vacuum_system_preference
        )
        self.baseline_purchase_costs['Vacuum system'] = vacuum_results['Cost']
        self.design_results['Vacuum system'] = vacuum_results['Name']
        hx_hu, vacuum_steam, vacuum_cooling_water = self.heat_utilities
        vacuum_steam.set_utility_by_flow_rate(vacuum_results['Heating agent'], vacuum_results['Steam flow rate'])
        if vacuum_results['Condenser']: 
            vacuum_cooling_water(-vacuum_steam.unit_duty, 373.15)
        else:
            vacuum_cooling_water.empty()
        self.power_utility(vacuum_results['Work'])

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
        K = design.compute_Stokes_law_York_Demister_K_value(P)

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

    def _vertical_vessel_pressure_diameter_and_length(self):
        rhov, rhol, P, Th, Ts, has_mist_eliminator, Qv, Qll, Ut, Uv, Vh, Vs = self._design_parameters()

        # Calculate internal diameter, Dvd
        Dvd = (4.0*Qv/(pi*Uv))**0.5
        if has_mist_eliminator:
            D = design.ceil_half_step(Dvd + 0.4)
        else:
            D = design.ceil_half_step(Dvd)

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
        dN = design.ceil_half_step(dN)

        # Calculate Hlin, assume with inlet diverter
        Hlin = 1.0 + dN

        # Calculate the vapor disengagement height
        Hv = 0.5*Dvd
        Hv2 = (2.0 if has_mist_eliminator else 3.0) + dN/2.0 # pragma: no cover
        if Hv2 < Hv: Hv = Hv2
        Hv = Hv

        # Calculate total height, Ht
        Hme = 1.5 if has_mist_eliminator else 0.0
        Ht = Hlll + Hh + Hs + Hlin + Hv + Hme
        Ht = design.ceil_half_step(Ht)

        # Find maximum and normal liquid level
        # Hhll = Hs + Hh + Hlll
        # Hnll = Hh + Hlll

        return P, D, Ht
        
    def _horizontal_vessel_pressure_diameter_and_length(self):
        rhov, rhol, P, Th, Ts, has_mist_eliminator, Qv, Qll, Ut, Uv, Vh, Vs = self._design_parameters()

        # Initialize LD
        if P > 0 and P <= 264.7:
            LD = 1.5/250.0*(P-14.7)+1.5
        elif P > 264.7 and P <= 514.7: # pragma: no cover
            LD = 1.0/250.0*(P-14.7)+2.0
        elif P > 514.7: # pragma: no cover
            LD = 5.0

        D = (4.0*(Vh+Vs)/(0.6*pi*LD))**(1.0/3.0)
        if D <= 4.0:
            D = 4.0
        else:
            D = design.ceil_half_step(D)

        for outerIter in range(50):
            At = pi*(D**2)/4.0 # Total area

            # Calculate Lower Liquid Area
            Hlll = round(0.5*D + 7.0)  
            Hlll = Hlll/12.0 # D is in ft but Hlll is in inches
            X = Hlll/D
            Y = design.HNATable(1, X)
            Alll = Y*At

            # Calculate the Vapor disengagement area, Av
            Hv = 0.2*D
            if has_mist_eliminator and Hv <= 2.0: Hv = 2.0
            elif Hv <= 1.0: Hv = 1.0
            else: Hv = design.ceil_half_step(Hv)
            Av = design.HNATable(1, Hv/D)*At
            
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
                Av = design.HNATable(1, Hv/D)*At
                Alll = design.HNATable(1, Hlll/D)*At
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

            if LD > 7.2: # pragma: no cover
                D += 0.5
            else: break

        # Recalculate LD so it lies between 1.5 - 6.0
        while True:
            LD = L / D
            if (LD < 1.5) and D <= 4.0: L += 0.5
            elif LD < 1.5: D -= 0.5
            elif (LD > 6.0): D += 0.5
            else: break

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
        
        return P, D, L


# %% Special

class SplitFlash(Flash):
    line = 'Flash' 
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *, split,
                 order=None, T=None, P=None, Q=None,
                 vessel_material='Carbon steel',
                 vacuum_system_preference='Liquid-ring pump',
                 has_glycol_groups=False,
                 has_amine_groups=False,
                 vessel_type=None,
                 holdup_time=15,
                 surge_time=7.5,
                 has_mist_eliminator=False):
        Splitter.__init__(self, ID, ins, outs, thermo, split=split, order=order)
        self._load_components()
        
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
        self.heat_exchanger.outs[0] = ms = self._multi_stream
        ms.mix_from(self.outs)
        super()._design()
    
# TODO: Remove this in favor of partition coefficients
class RatioFlash(Flash):
    _N_heat_utilities = 1

    def __init__(self, ID='', ins=None, outs=(), *,
                 K_chemicals, Ks, top_solvents=(), top_split=(),
                 bot_solvents=(), bot_split=()):
        Unit.__init__(self, ID, ins, outs)
        self._load_components()
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

    def _design(self): # pragma: no cover
        self.heat_exchanger.outs[0] = ms = self._multi_stream
        ms.mix_from(self.outs)
        super()._design()

# %% Single Component

class Evaporator_PQ(Unit):
    _N_ins = 2
    _N_outs = 3
    @property
    def P(self):
        return self._P
    @P.setter
    def P(self, P):
        water = getattr(self.chemicals, '7732-18-5')
        self._T = T = water.Tsat(P)
        self._Hvap = water.Hvap(T)
        self._P = P
    @property
    def T(self):
        return self._T
    @T.setter
    def T(self, T): # pragma: no cover
        water = getattr(self.chemicals, '7732-18-5')
        self._P = water.Psat(T)
        self._Hvap = water.Hvap(T)
        self._T = T
    @property
    def V(self):
        return self._V
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *, Q=0, P=101325):
        super().__init__(ID, ins, outs, thermo)
        self.Q = Q
        self.P = P
        self._V = None
    
    def _run(self):
        feed, utility_vapor = self.ins
        vapor, liquid, utility_liquid = self.outs
        
        # Optional if Q also comes from condensing a side stream
        Q = self.Q
        if not utility_vapor.isempty():
            utility_liquid.copy_like(utility_vapor)
            utility_liquid.phase = 'l'
            Q += utility_vapor.Hvap
        feed_H = feed.H
    
        # Set exit conditions
        vapor.T = liquid.T = self.T
        vapor.P = liquid.P = self.P
        liquid.phase = 'l'
        vapor.phase = 'g'
        liquid.copy_flow(feed)
        
        # Energy balance to find vapor fraction
        f = feed.imol['7732-18-5']
        H = feed_H + Q - liquid.H
        if f:
            V = H/(f * self._Hvap)
            if V < 0: 
                V = 0
            elif V > 1:
                V = 1
        else:
            V = 0
        evaporated = f * V
        vapor.imol['7732-18-5'] = evaporated
        liquid.imol['7732-18-5'] = (1-V)*f
        self.design_results['Heat transfer'] = Q
        self._V = V


class Evaporator_PV(Unit):
    _N_heat_utilities = 0
    _N_outs = 2
    
    @property
    def P(self):
        return self._P
    @P.setter
    def P(self, P):
        water = getattr(self.chemicals, '7732-18-5')
        self._T = water.Tsat(P)
        self._P = P
    @property
    def T(self):
        return self._T
    @T.setter
    def T(self, T): # pragma: no cover
        water = getattr(self.chemicals, '7732-18-5')
        self._P = water.Psat(T)
        self._T = T
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *, V=0.5, P=101325):
        super().__init__(ID, ins, outs, thermo)
        self._multi_stream = ms = MultiStream(None, thermo=self.thermo)
        self.heat_exchanger = hx = HXutility(None, None, ms, thermo=self.thermo) 
        self.heat_utilities = hx.heat_utilities
        self.V = V
        self.P = P
        hx._ins = self._ins

    def _run(self):
        feed = self.ins[0]
        vapor, liquid = self.outs
        vapor.phase = 'g'
        vapor.T = liquid.T = self.T
        vapor.P = liquid.P = self.P
        H2O_index = self.chemicals.index('7732-18-5')
        water_mol = feed.mol[H2O_index]
        vapor.mol[H2O_index] = self.V * water_mol
        liquid_mol = liquid.mol
        liquid_mol[:] = feed.mol
        liquid_mol[H2O_index] = (1-self.V) * water_mol
        ms = self._multi_stream
        self.design_results['Heat transfer'] = self.heat_exchanger.Q = (vapor.H + liquid.H) - feed.H
        ms['g'].copy_like(vapor)
        ms['l'].copy_like(liquid)