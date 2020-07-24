# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import numpy as np
from math import ceil
import thermosteam as tmo
from .design_tools.specification_factors import  (
    distillation_column_material_factors,
    tray_material_factor_functions,
    distillation_tray_type_factor,
    material_densities_lb_per_in3)
from .design_tools import column_design as design
from .. import Unit
from .._graphics import vertical_column_graphics
from scipy.optimize import brentq
from ._hx import HXutility
from ..exceptions import InfeasibleRegion
import matplotlib.pyplot as plt

__all__ = ('BinaryDistillation',)

def compute_stages_McCabeThiele(P, operating_line,
                                x_stages, y_stages, T_stages,
                                x_limit, solve_Ty):
    """
    Use the McCabe-Thiele method to find the specifications at every stage of
    the operating line before the maximum liquid molar fraction, `x_limit`. 
    Append the light key liquid molar fraction, light key vapor molar
    fraction, and stage temperatures to `x_stages`, `y_stages` and `T_stages`
    respectively.
    
    Parameters
    ----------
    P : float
        Pressure [Pa].
    operating_line : function
                     Should return the liquid molar fraction of the light
                     key given its vapor molar fraction.
    x_stages : list
               Liquid molar compositions at each stage. Last element
               should be the starting point for the next stage.
    y_stages : list
               Vapor molar compositions at each stage. Last element 
               should be the starting point for the next stage.
    x_limit : float
              Maximum value of liquid composition before algorithm stops.
    T_stages : list
               Temperature at each stage.
    solve_Ty : function
               Should return T and y given x.
        
    """
    i = 0
    yi = y_stages[-1]
    xi = x_stages[-1]
    while xi < x_limit:
        if i > 100:
            raise RuntimeError('cannot meet specifications! stages > 100')
        i += 1
        # Go Up
        x = np.array((xi, 1-xi))
        T, y = solve_Ty(x, P)
        yi = y[0]
        y_stages.append(yi)
        T_stages.append(T)
        # Go Right
        xi = operating_line(yi)
        if xi > x_limit:
            xi = x_limit
        x_stages.append(xi)
    

# %% Distillation

class BinaryDistillation(Unit):
    r"""
    Create a binary distillation column that assumes all light and heavy non keys
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
    and flooding velocity [1]_ [3]_. Purchase costs are based on correlations
    compiled by Warren et. al. [4]_.

    Parameters
    ----------
    ins : streams
        Inlet fluids to be mixed into the feed stage.
    outs : stream sequence
        * [0] Distillate
        * [1] Bottoms product
    LHK : tuple[str]
        Light and heavy keys.
    y_top : float
        Molar fraction of light key to the light and heavy keys in the
        distillate.
    x_bot : float
        Molar fraction of light key to the light and heavy keys in the bottoms
        product.
    Lr : float
        Recovery of the light key in the distillate.
    Hr : float
        Recovery of the heavy key in the bottoms product.
    k : float
        Ratio of reflux to minimum reflux.
    Rmin : float, optional
        User enforced minimum reflux ratio. If the actual minimum reflux ratio
        is more than `Rmin`, this enforced value is ignored. Defaults to 0.3.
    product_specification_format=None : "Composition" or "Recovery"
        If composition is used, `y_top` and `x_bot` must be specified.
        If recovery is used, `Lr` and `Hr` must be specified.
    P=101325 : float
        Operating pressure [Pa].
    vessel_material : str, optional
        Vessel construction material. Defaults to 'Carbon steel'.
    tray_material : str, optional
        Tray construction material. Defaults to 'Carbon steel'.
    tray_type='Sieve' : 'Sieve', 'Valve', or 'Bubble cap'
        Tray type.
    tray_spacing=450 : float
        Typically between 152 to 915 mm.
    stage_efficiency=None : 
        User enforced stage efficiency. If None, stage efficiency is
        calculated by the O'Connell correlation [2]_.
    velocity_fraction=0.8 : float
        Fraction of actual velocity to maximum velocity allowable before
        flooding.
    foaming_factor=1.0 : float
        Must be between 0 to 1.
    open_tray_area_fraction=0.1 : float
        Fraction of open area to active area of a tray.
    downcomer_area_fraction=None : float
        Enforced fraction of downcomer area to net (total) area of a tray.
        If None, estimate ratio based on Oliver's estimation [1]_.
    is_divided=False : bool
        True if the stripper and rectifier are two separate columns.

    References
    ----------
    .. [1] J.D. Seader, E.J. Henley, D.K. Roper. (2011)
        Separation Process Principles 3rd Edition. John Wiley & Sons, Inc. 

    .. [2] M. Duss, R. Taylor. (2018)
        Predict Distillation Tray Efficiency. AICHE 
    
    .. [3] Green, D. W. Distillation. In Perry’s Chemical Engineers’
        Handbook, 9 ed.; McGraw-Hill Education, 2018.

    .. [4] Seider, W. D., Lewin,  D. R., Seader, J. D., Widagdo, S., Gani, R.,
        & Ng, M. K. (2017). Product and Process Design Principles. Wiley.
        Cost Accounting and Capital Cost Estimation (Chapter 16)

    Examples
    --------
    Binary distillation assuming 100% separation on non-keys:
    
    >>> from biosteam.units import BinaryDistillation
    >>> from biosteam import Stream, settings
    >>> settings.set_thermo(['Water', 'Methanol', 'Glycerol'])
    >>> feed = Stream('feed', flow=(80, 100, 25))
    >>> bp = feed.bubble_point_at_P()
    >>> feed.T = bp.T # Feed at bubble point T
    >>> D1 = BinaryDistillation('D1', ins=feed,
    ...                         outs=('distillate', 'bottoms_product'),
    ...                         LHK=('Methanol', 'Water'),
    ...                         y_top=0.99, x_bot=0.01, k=2,
    ...                         is_divided=True)
    >>> D1.simulate()
    >>> # See all results
    >>> D1.show(T='degC', P='atm', composition=True)
    BinaryDistillation: D1
    ins...
    [0] feed
        phase: 'l', T: 76.129 degC, P: 1 atm
        composition: Water     0.39
                     Methanol  0.488
                     Glycerol  0.122
                     --------  205 kmol/hr
    outs...
    [0] distillate
        phase: 'g', T: 64.91 degC, P: 1 atm
        composition: Water     0.01
                     Methanol  0.99
                     --------  100 kmol/hr
    [1] bottoms_product
        phase: 'l', T: 100.06 degC, P: 1 atm
        composition: Water     0.754
                     Methanol  0.00761
                     Glycerol  0.239
                     --------  105 kmol/hr
    >>> D1.results()
    Distillation                                    Units        D1
    Cooling water       Duty                        kJ/hr -5.12e+06
                        Flow                      kmol/hr   3.5e+03
                        Cost                       USD/hr      1.71
    Low pressure steam  Duty                        kJ/hr  1.09e+07
                        Flow                      kmol/hr       280
                        Cost                       USD/hr      66.5
    Design              Theoretical feed stage                    9
                        Theoretical stages                       13
                        Minimum reflux              Ratio     0.687
                        Reflux                      Ratio      1.37
                        Rectifier stages                         15
                        Stripper stages                          13
                        Rectifier height               ft      34.7
                        Stripper height                ft      31.7
                        Rectifier diameter             ft      3.93
                        Stripper diameter              ft      3.21
                        Rectifier wall thickness       in     0.312
                        Stripper wall thickness        in     0.312
                        Rectifier weight               lb  6.01e+03
                        Stripper weight                lb  4.46e+03
    Purchase cost       Rectifier trays               USD   1.5e+04
                        Stripper trays                USD  1.25e+04
                        Rectifier tower               USD  8.18e+04
                        Stripper tower                USD  6.81e+04
                        Condenser                     USD  3.39e+04
                        Boiler                        USD  2.72e+04
    Total purchase cost                               USD  2.39e+05
    Utility cost                                   USD/hr      68.2
    
    """
    line = 'Distillation'
    auxiliary_unit_names = ('condenser', 'boiler')
    _graphics = vertical_column_graphics
    _N_heat_utilities = 0
    _ins_size_is_fixed = False
    _N_ins = 1
    _N_outs = 2
    _units = {'Minimum reflux': 'Ratio',
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
    _BM = {'Rectifier tower': 4.3,
           'Stripper tower': 4.3,
           'Rectifier trays': 4.3,
           'Stripper trays': 4.3,
           'Tower': 4.3,
           'Trays': 4.3}
    
    # [dict] Bounds for results
    _bounds = {'Diameter': (3., 24.),
               'Height': (27., 170.),
               'Weight': (9000., 2.5e6)}
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                P=101325, *, LHK, k,
                Rmin=0.3,
                Lr=None,
                Hr=None,
                y_top=None,
                x_bot=None, 
                product_specification_format=None,
                vessel_material='Carbon steel',
                tray_material='Carbon steel',
                tray_type='Sieve',
                tray_spacing=450,
                stage_efficiency=None,
                velocity_fraction=0.8,
                foaming_factor=1.0,
                open_tray_area_fraction=0.1,
                downcomer_area_fraction=None,
                is_divided=False):
        Unit.__init__(self, ID, ins, outs, thermo)
        
        # Operation specifications
        self.k = k
        self.P = P
        self.Rmin = Rmin
        self.LHK = LHK
        self._set_distillation_product_specifications(product_specification_format,
                                                      x_bot, y_top, Lr, Hr)
        
        # Construction specifications
        self.vessel_material = vessel_material
        self.tray_type = tray_type
        self.tray_material = tray_material
        self.tray_spacing = tray_spacing
        self.stage_efficiency = stage_efficiency
        self.velocity_fraction = velocity_fraction
        self.foaming_factor = foaming_factor
        self.open_tray_area_fraction = open_tray_area_fraction
        self.downcomer_area_fraction = downcomer_area_fraction
        self.is_divided = is_divided
        
        # Setup components
        thermo = self.thermo
        
        #: [MultiStream] Overall feed to the distillation column
        self.feed = tmo.MultiStream(None, thermo=thermo)
        
        #: [HXutility] Condenser.
        self.condenser = HXutility(None,
                                   ins=tmo.Stream(None, phase='g', thermo=thermo),
                                   outs=tmo.MultiStream(None, thermo=thermo),
                                   thermo=thermo)
        #: [HXutility] Boiler.
        self.boiler = HXutility(None,
                                ins=tmo.Stream(None, thermo=thermo),
                                outs=tmo.MultiStream(None, thermo=thermo),
                                thermo=thermo)
        self.heat_utilities = self.condenser.heat_utilities + self.boiler.heat_utilities
        self.reset_cache()
        
    def reset_cache(self):
        self._McCabeThiele_args = np.zeros(6)
    
    @property
    def product_specification_format(self):
        return self._product_specification_format
    @product_specification_format.setter
    def product_specification_format(self, spec):
        if spec == 'Composition':
            self._Lr = self._Hr = None
        elif spec == 'Recovery':
            self._y_top = self._x_bot = None
        else:
            raise AttributeError("product specification format must be either "
                                 "'Composition' or 'Recovery'")
        self._product_specification_format = spec
    
    @property
    def condensate(self):
        return self.condenser.outs[0]['l']
    @property
    def boilup(self):
        return self.boiler.outs[0]['g']    
    
    @property
    def LHK(self):
        """tuple[str, str] Light and heavy keys."""
        return self._LHK
    @LHK.setter
    def LHK(self, LHK):
        # Set light non-key and heavy non-key indices
        self._LHK = LHK = tuple(LHK)
        intermediate_volatile_chemicals = []
        chemicals = self.chemicals
        LHK_chemicals = LK_chemical, HK_chemical = self.chemicals[LHK]
        Tb_light = LK_chemical.Tb
        Tb_heavy = HK_chemical.Tb
        LNK = []
        HNK = []
        gases = []
        solids = []
        for chemical in chemicals:
            ID = chemical.ID
            Tb = chemical.Tb
            if not Tb or chemical.locked_state in ('l', 's'):
                solids.append(ID)
            elif chemical.locked_state == 'g':
                gases.append(ID)
            elif Tb < Tb_light:
                LNK.append(ID)
            elif Tb > Tb_heavy:
                HNK.append(ID)
            elif chemical not in LHK_chemicals:
                intermediate_volatile_chemicals.append(chemical.ID)
        self._LNK = LNK = tuple(LNK)
        self._HNK = HNK = tuple(HNK)
        self._gases = gases = tuple(gases)
        self._solids = solids = tuple(solids)
        get_index = self.chemicals.get_index
        self._LHK_index = get_index(LHK)
        self._LNK_index = get_index(LNK)
        self._HNK_index = get_index(HNK)
        self._gases_index = get_index(gases)
        self._solids_index = get_index(solids)
        self._intermediate_volatile_chemicals = tuple(intermediate_volatile_chemicals)
    
    @property
    def Rmin(self):
        """User enforced minimum reflux ratio. If the actual minimum reflux ratio is less than `Rmin`. This enforced value is ignored."""
        return self._Rmin
    @Rmin.setter
    def Rmin(self, Rmin):
        self._Rmin = Rmin
    
    @property
    def y_top(self):
        """Light key composition of at the distillate."""
        return self._y_top
    @y_top.setter
    def y_top(self, y_top):
        assert self.product_specification_format == "Composition", (
            "product specification format must be 'Composition' "
            "to set distillate composition")
        assert 0 < y_top < 1, "light key composition in the distillate must be a fraction" 
        self._y_top = y_top
        self._y = np.array([y_top, 1-y_top])
    
    @property
    def x_bot(self):
        """Light key composition at the bottoms product."""
        return self._x_bot
    @x_bot.setter
    def x_bot(self, x_bot):
        assert self.product_specification_format == "Composition", (
            "product specification format must be 'Composition' to set bottoms "
            "product composition")
        assert 0 < x_bot < 1, "heavy key composition in the bottoms product must be a fraction" 
        self._x_bot = x_bot
        self._x = np.array([x_bot, 1-x_bot])
    
    @property
    def Lr(self):
        """Light key recovery in the distillate."""
        return self._Lr
    @Lr.setter
    def Lr(self, Lr):
        assert self.product_specification_format == "Recovery", (
            "product specification format must be 'Recovery' "
            "to set light key recovery")
        assert 0 < Lr < 1, "light key recovery in the distillate must be a fraction" 
        self._Lr = Lr
    
    @property
    def Hr(self):
        """Heavy key recovery in the bottoms product."""
        return self._Hr
    @Hr.setter
    def Hr(self, Hr):
        assert self.product_specification_format == "Recovery", (
            "product specification format must be 'Recovery' "
            "to set heavy key recovery")
        assert 0 < Hr < 1, "heavy key recovery in the bottoms product must be a fraction" 
        self._Hr = Hr
    
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
        return self._E_eff
    @stage_efficiency.setter
    def stage_efficiency(self, E_eff):
        self._E_eff = E_eff
    
    @property
    def velocity_fraction(self):
        """Fraction of actual velocity to maximum velocity allowable before flooding."""
        return self._f
    @velocity_fraction.setter
    def velocity_fraction(self, f):
        self._f = f
    
    @property
    def foaming_factor(self):
        """Foaming factor (0 to 1)."""
        return self._F_F
    @foaming_factor.setter
    def foaming_factor(self, F_F):
        if not 0 <= F_F <= 1:
            raise ValueError(f"foaming_factor must be between 0 and 1, ({F_F} given).")
        self._F_F = F_F
    
    @property
    def open_tray_area_fraction(self):
        """Fraction of open area, A_h, to active area, A_a."""
        return self._A_ha
    @open_tray_area_fraction.setter
    def open_tray_area_fraction(self, A_ha):
        self._A_ha = A_ha
    
    @property
    def downcomer_area_fraction(self):
        """Enforced fraction of downcomer area to net (total) area.
        If None, the fraction is estimated based on heuristics."""
        return self._A_dn
    @downcomer_area_fraction.setter
    def downcomer_area_fraction(self, A_dn):
        self._A_dn = A_dn
    
    @property
    def tray_type(self):
        """Default 'Sieve'"""
        return self._tray_type
    @tray_type.setter
    def tray_type(self, tray_type):
        if tray_type in distillation_tray_type_factor:
            self._tray_type = tray_type
            self._F_TT = distillation_tray_type_factor[tray_type]
        else:
            raise ValueError("tray type must be one of the following: "
                            f"{', '.join(distillation_tray_type_factor)}")
        
    @property
    def tray_material(self):
        """Default 'Carbon steel'"""
        return self._tray_material
    @tray_material.setter
    def tray_material(self, tray_material):
        if tray_material in tray_material_factor_functions:
            self._tray_material = tray_material
            self._F_TM_function = tray_material_factor_functions[tray_material]
        else:
            raise ValueError("tray material must be one of the following: "
                            f"{', '.join(tray_material_factor_functions)}")
        
    @property
    def vessel_material(self):
        """Default 'Carbon steel'"""
        return self._vessel_material
    @vessel_material.setter
    def vessel_material(self, vessel_material):
        if vessel_material in distillation_column_material_factors:
            self._vessel_material = vessel_material
            self._F_VM = distillation_column_material_factors[vessel_material]
        else:
            raise ValueError("vessel material must be one of the following: "
                            f"{', '.join(distillation_column_material_factors)}")
    
    @property
    def is_divided(self):
        """[bool] True if the stripper and rectifier are two separate columns."""
        return self._is_divided
    @is_divided.setter
    def is_divided(self, is_divided):
        self._is_divided = is_divided
    
    def _set_distillation_product_specifications(self,
                                                 product_specification_format,
                                                 x_bot, y_top, Lr, Hr):
        if not product_specification_format:
            if (x_bot and y_top) and not (Lr or Hr):
                product_specification_format = 'Composition'
            elif (Lr and Hr) and not (x_bot or y_top):
                product_specification_format = 'Recovery'
            else:
                raise ValueError("must specify either x_top and y_top, or Lr and Hr")
        self._product_specification_format = product_specification_format
        if product_specification_format == 'Composition':
            self.y_top = y_top
            self.x_bot = x_bot
        elif product_specification_format == 'Recovery':
            self.Lr = Lr
            self.Hr = Hr
        else:
            raise ValueError("product specification format must be either 'Composition' or 'Recovery'")
    
    def _get_y_top_and_x_bot(self):
        if self.product_specification_format == 'Composition':
            y_top = self.y_top
            x_bot = self.x_bot
        else:
            distillate, bottoms_product = self.outs
            LHK = self._LHK
            y_top, _ = distillate.get_normalized_mol(LHK)
            x_bot, _ = bottoms_product.get_normalized_mol(LHK)
        return y_top, x_bot
    
    def _check_mass_balance(self):
        distillate, bottoms_product = self.outs
        LHK = self._LHK
        LHK_index = self._LHK_index
        LK_distillate, HK_distillate = distillate.mol[LHK_index]
        LK_bottoms, HK_bottoms = bottoms_product.mol[LHK_index]
        if self.product_specification_format == 'Composition':
            if LK_distillate < 0. or LK_bottoms < 0.:
                raise InfeasibleRegion("light key composition")
            if HK_distillate < 0. or HK_bottoms < 0.:
                raise InfeasibleRegion("heavy key composition")
        if tmo.settings.debug:
            intermediate_chemicals = self._intermediate_volatile_chemicals
            intemediates_index = self.chemicals.get_index(intermediate_chemicals)
            intermediate_flows = self.mol_in[intemediates_index]
            minflow = min(LK_distillate, HK_bottoms)
            for flow, chemical in zip(intermediate_flows, intermediate_chemicals):
                assert flow > minflow, ("significant intermediate volatile chemical,"
                                        f"'{chemical}', between light and heavy "
                                        f"key, {', '.join(LHK)}")
    
    def _run_binary_distillation_mass_balance(self):
        # Get all important flow rates (both light and heavy keys and non-keys)
        feed = self.feed
        feed.mix_from(self.ins)
        mol = feed.mol
        LHK_index = self._LHK_index
        LNK_index = self._LNK_index
        HNK_index = self._HNK_index
        gases_index = self._gases_index
        solids_index = self._solids_index
        intermediate_chemicals = self._intermediate_volatile_chemicals
        intemediates_index = self.chemicals.get_index(intermediate_chemicals)
        LHK_mol = mol[LHK_index]
        LNK_mol = mol[LNK_index]
        HNK_mol = mol[HNK_index]
        gases_mol = mol[gases_index]
        solids_mol = mol[solids_index]
        
        # Mass balance for non-keys
        distillate, bottoms_product = self.outs
        distillate.mol[LNK_index] = LNK_mol
        distillate.mol[gases_index] = gases_mol
        bottoms_product.mol[HNK_index] = HNK_mol
        bottoms_product.mol[solids_index] = solids_mol
        
        # Mass balance for keys
        spec = self.product_specification_format
        if  spec == 'Composition':
            # Use lever rule
            light, heavy = LHK_mol
            F_mol_LHK = light + heavy
            zf = light/F_mol_LHK
            distillate_fraction = (zf-self.x_bot)/(self.y_top-self.x_bot)
            F_mol_LHK_distillate = F_mol_LHK * distillate_fraction
            distillate_LHK_mol = F_mol_LHK_distillate * self._y
        elif spec == 'Recovery':
            distillate_LHK_mol = LHK_mol * [self.Lr, (1 - self.Hr)]
        else:
            raise ValueError("invalid specification '{spec}'")
        distillate.mol[LHK_index] = distillate_LHK_mol
        bottoms_product.mol[LHK_index] = LHK_mol - distillate_LHK_mol
        distillate.mol[intemediates_index] = \
        bottoms_product.mol[intemediates_index] = mol[intemediates_index] / 2
        self._check_mass_balance()
    
    def _update_distillate_and_bottoms_temperature(self):
        distillate, bottoms_product = self.outs
        self._condensate_dew_point = dp = distillate.dew_point_at_P()
        self._boilup_bubble_point = bp = bottoms_product.bubble_point_at_P()
        bottoms_product.T = bp.T
        distillate.T = dp.T
    
    def _setup(self):
        distillate, bottoms_product = self.outs
        distillate.P = bottoms_product.P = self.P
        distillate.phase = 'g'
        bottoms_product.phase = 'l'
    
    def _run(self):
        self._run_binary_distillation_mass_balance()
        self._update_distillate_and_bottoms_temperature()

    def _get_feed_quality(self):
        feed = self.feed
        feed = feed.copy()
        H_feed = feed.H
        try: dp = feed.dew_point_at_P()
        except: pass
        else: feed.T = dp.T
        feed.phase = 'g'
        H_vap = feed.H
        try: bp = feed.bubble_point_at_P()
        except: pass
        else: feed.T = bp.T
        feed.phase = 'l'
        H_liq = feed.H
        q = (H_vap - H_feed) / (H_vap - H_liq)
        return q

    def _run_McCabeThiele(self):
        distillate, bottoms = self.outs
        chemicals = self.chemicals
        LHK = self._LHK
        LHK_index = chemicals.get_index(LHK)

        # Feed light key mol fraction
        feed = self.feed
        liq_mol = feed.imol['l']
        vap_mol = feed.imol['g']
        LHK_mol = liq_mol[LHK_index] + vap_mol[LHK_index]
        F_mol_LHK = LHK_mol.sum()
        zf = LHK_mol[0]/F_mol_LHK
        q = self._get_feed_quality()
        
        # Main arguments
        P = self.P
        k = self.k
        y_top, x_bot = self._get_y_top_and_x_bot()
        
        # Cache
        old_args = self._McCabeThiele_args
        args = np.array([P, k, y_top, x_bot, q, zf])
        tol = np.array([50, 1e-5, 1e-6, 1e-6, 1e-2, 1e-6], float)
        if (abs(old_args - args) < tol).all(): return
        self._McCabeThiele_args = args
        
        # Get R_min and the q_line 
        if abs(q - 1) < 1e-4:
            q = 1 - 1e-4
        q_line = lambda x: q*x/(q-1) - zf/(q-1)
        self._q_line_args = dict(q=q, zf=zf)
        
        solve_Ty = bottoms.get_bubble_point(LHK).solve_Ty
        Rmin_intersection = lambda x: q_line(x) - solve_Ty(np.array((x, 1-x)), P)[1][0]
        x_Rmin = brentq(Rmin_intersection, 0, 1)
        y_Rmin = q_line(x_Rmin)
        m = (y_Rmin-y_top)/(x_Rmin-y_top)
        Rmin = m/(1-m)
        if Rmin < self._Rmin:
            Rmin = self._Rmin
        R = k * Rmin

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
        compute_stages_McCabeThiele(P, ss, x_stages, y_stages, T_stages, x_m, solve_Ty)
        yi = y_stages[-1]
        xi = rs(yi)
        x_stages[-1] = xi if xi < 1 else 0.99999
        compute_stages_McCabeThiele(P, rs, x_stages, y_stages, T_stages, y_top, solve_Ty)
        
        # Find feed stage
        N_stages = len(x_stages)
        feed_stage = ceil(N_stages/2)
        for i in range(len(y_stages)-1):
            if y_stages[i] < y_m < y_stages[i+1]:
                feed_stage = i+1
        
        # Results
        Design = self.design_results
        Design['Theoretical feed stage'] = N_stages - feed_stage
        Design['Theoretical stages'] = N_stages
        Design['Minimum reflux'] = Rmin
        Design['Reflux'] = R
        
    def _run_condenser_and_boiler(self):
        feed = self.feed
        distillate, bottoms_product = self.outs
        condenser = self.condenser
        boiler = self.boiler
        R = self.design_results['Reflux']
        
        # Set condenser conditions
        condenser.outs[0].imol['g'] = distillate.mol
        self._F_mol_distillate = F_mol_distillate = distillate.F_mol
        self._F_mol_condensate = F_mol_condensate = R * F_mol_distillate
        dp = self._condensate_dew_point
        condensate_x_mol = dp.x
        condensate = self.condensate
        condensate.empty()
        condensate.imol[dp.IDs] = condensate_x_mol * F_mol_condensate
        condensate.T = dp.T
        condensate.P = dp.P
        vap = condenser.ins[0]
        vap.mol = distillate.mol + condensate.mol
        vap_T = vap.dew_point_at_P().T
        if vap_T < dp.T: vap_T = dp.T + 0.1
        vap.T = vap_T
        vap.P = distillate.P
        
        # Set boiler conditions
        boiler.outs[0].imol['l'] = bottoms_product.mol
        F_vap_feed = feed.imol['g'].sum()
        self._F_mol_boilup = F_mol_boilup = (R+1)*F_mol_distillate - F_vap_feed
        bp = self._boilup_bubble_point
        boilup_flow = bp.y * F_mol_boilup
        boilup = self.boilup
        boilup.T = bp.T
        boilup.P = bp.P
        boilup.imol[bp.IDs] = boilup_flow
        liq = boiler.ins[0]
        liq.phase = 'l'
        liq.mol = bottoms_product.mol + boilup.mol
        liq_T = liq.bubble_point_at_P().T
        if liq_T > bp.T: liq_T = bp.T - 0.1
        liq.T = liq_T
    
    def _simulate_components(self): 
        boiler = self.boiler
        condenser = self.condenser
        Q_condenser = condenser.outs[0].H - condenser.ins[0].H
        H_out = self.H_out
        H_in = self.H_in
        Q_overall_boiler =  H_out - H_in - Q_condenser
        Q_boiler = boiler.outs[0].H - boiler.ins[0].H
        if Q_boiler < Q_overall_boiler:
            liquid = boiler.ins[0]
            H_out_boiler = boiler.outs[0].H
            liquid.H = H_out_boiler - Q_overall_boiler
            boiler._design(Q_overall_boiler)
            condenser._design(Q_condenser)
        else:
            boiler._design(Q_boiler)
            condenser._design(Q_condenser)
        boiler._cost()
        condenser._cost()
    
    def _get_relative_volatilities_LHK(self):
        x_stages = self._x_stages
        y_stages = self._y_stages
        
        K_light = y_stages[-1]/x_stages[-1] 
        K_heavy = (1-y_stages[-1])/(1-x_stages[-1])
        alpha_LHK_distillate = K_light/K_heavy
        
        K_light = y_stages[0]/x_stages[0] 
        K_heavy = (1-y_stages[0])/(1-x_stages[0] )
        alpha_LHK_bottoms = K_light/K_heavy
        
        return alpha_LHK_distillate, alpha_LHK_bottoms
    
    def _compute_N_stages(self):
        """Return a tuple with the actual number of stages for the rectifier and the stripper."""
        feed = self.feed
        vap, liq = self.outs
        Design = self.design_results
        R = Design['Reflux']
        N_stages = Design['Theoretical stages']
        feed_stage = Design['Theoretical feed stage']
        E_eff = self.stage_efficiency
        if E_eff:
            E_rectifier = E_stripper = E_eff
        else:    
            # Calculate Murphree Efficiency for rectifying section
            condensate = self.condensate
            mu = condensate.get_property('mu', 'mPa*s')
            alpha_LHK_distillate, alpha_LHK_bottoms = self._get_relative_volatilities_LHK()
            F_mol_distillate = self._F_mol_distillate
            L_Rmol = self._F_mol_condensate
            V_Rmol = (R+1) * F_mol_distillate
            E_rectifier = design.compute_murphree_stage_efficiency(mu,
                                                            alpha_LHK_distillate,
                                                            L_Rmol, V_Rmol)
            
            # Calculate Murphree Efficiency for stripping section
            mu = liq.get_property('mu', 'mPa*s')
            V_Smol = self._F_mol_boilup
            L_Smol = R*F_mol_distillate + feed.imol['g'].sum()
            E_stripper = design.compute_murphree_stage_efficiency(mu,
                                                           alpha_LHK_bottoms,
                                                           L_Smol, V_Smol)
            
        # Calculate actual number of stages
        mid_stage = feed_stage - 0.5
        N_rectifier = np.ceil(mid_stage/E_rectifier)
        N_stripper = np.ceil((N_stages-mid_stage)/E_stripper)
        return N_rectifier, N_stripper
        
    def _design(self):
        self._run_McCabeThiele()
        self._run_condenser_and_boiler()
        self._complete_distillation_column_design()
        
    def _complete_distillation_column_design(self):
        distillate, bottoms_product = self.outs
        Design = self.design_results
        R = Design['Reflux']
        Rstages, Sstages = self._compute_N_stages()
        is_divided = self.is_divided
        TS = self._TS
        
        ### Get diameter of rectifying section based on top plate ###
        
        condensate = self.condensate
        rho_L = condensate.rho
        sigma = condensate.get_property('sigma', 'dyn/cm')
        L = condensate.F_mass
        V = L*(R+1)/R
        vap = self.condenser.ins[0]
        V_vol = vap.get_total_flow('m^3/s')
        rho_V = distillate.rho
        F_LV = design.compute_flow_parameter(L, V, rho_V, rho_L)
        C_sbf = design.compute_max_capacity_parameter(TS, F_LV)
        F_F = self._F_F
        A_ha = self._A_ha
        U_f = design.compute_max_vapor_velocity(C_sbf, sigma, rho_L, rho_V, F_F, A_ha)
        A_dn = self._A_dn
        if A_dn is None:
           self._A_dn = A_dn = design.compute_downcomer_area_fraction(F_LV)
        f = self._f
        R_diameter = design.compute_tower_diameter(V_vol, U_f, f, A_dn) * 3.28
        
        ### Get diameter of stripping section based on feed plate ###
        rho_L = bottoms_product.rho
        boilup = self.boilup
        V = boilup.F_mass
        V_vol = boilup.get_total_flow('m^3/s')
        rho_V = boilup.rho
        L = bottoms_product.F_mass # To get liquid going down
        F_LV = design.compute_flow_parameter(L, V, rho_V, rho_L)
        C_sbf = design.compute_max_capacity_parameter(TS, F_LV)
        sigma = condensate.get_property('sigma', 'dyn/cm')
        U_f = design.compute_max_vapor_velocity(C_sbf, sigma, rho_L, rho_V, F_F, A_ha)
        A_dn = self._A_dn
        if A_dn is None:
            A_dn = design.compute_downcomer_area_fraction(F_LV)
        S_diameter = design.compute_tower_diameter(V_vol, U_f, f, A_dn) * 3.28
        Po = self.P * 0.000145078 # to psi
        rho_M = material_densities_lb_per_in3[self.vessel_material]
        
        if is_divided:
            Design['Rectifier stages'] = Rstages
            Design['Stripper stages'] =  Sstages
            Design['Rectifier height'] = H_R = design.compute_tower_height(TS, Rstages-1) * 3.28
            Design['Stripper height'] = H_S = design.compute_tower_height(TS, Sstages-1) * 3.28
            Design['Rectifier diameter'] = R_diameter
            Design['Stripper diameter'] = S_diameter
            Design['Rectifier wall thickness'] = tv = design.compute_tower_wall_thickness(Po, R_diameter, H_R)
            Design['Stripper wall thickness'] = tv = design.compute_tower_wall_thickness(Po, S_diameter, H_S)
            Design['Rectifier weight'] = design.compute_tower_weight(R_diameter, H_R, tv, rho_M)
            Design['Stripper weight'] = design.compute_tower_weight(S_diameter, H_S, tv, rho_M)
        else:
            Design['Actual stages'] = Rstages + Sstages
            Design['Height'] = H = design.compute_tower_height(TS, Rstages+Sstages-2) * 3.28
            Design['Diameter'] = Di = max((R_diameter, S_diameter))
            Design['Wall thickness'] = tv = design.compute_tower_wall_thickness(Po, Di, H)
            Design['Weight'] = design.compute_tower_weight(Di, H, tv, rho_M)
        
    def _cost(self):
        Design = self.design_results
        Cost = self.purchase_costs
        Cost.clear() # Prevent having previous results if `is_divided` changed
        F_TT = self._F_TT
        F_VM = self._F_VM
        if self.is_divided:
            # Number of trays assuming a partial condenser
            N_RT = Design['Rectifier stages'] - 1
            Di_R = Design['Rectifier diameter']
            F_TM = self._F_TM_function(Di_R)
            Cost['Rectifier trays'] = design.compute_purchase_cost_of_trays(N_RT, Di_R, F_TT, F_TM)
            N_ST = Design['Stripper stages'] - 1
            Di_S = Design['Stripper diameter']
            F_TM = self._F_TM_function(Di_R)
            Cost['Stripper trays'] = design.compute_purchase_cost_of_trays(N_ST, Di_S, F_TT, F_TM)
            
            # Cost vessel assuming T < 800 F
            W_R = Design['Rectifier weight'] # in lb
            H_R = Design['Rectifier height']*3.28 # in ft
            Cost['Rectifier tower'] = design.compute_purchase_cost_of_tower(Di_R, H_R, W_R, F_VM)
            W_S = Design['Stripper weight'] # in lb
            H_S = Design['Stripper height']*3.28 # in ft
            Cost['Stripper tower'] = design.compute_purchase_cost_of_tower(Di_S, H_S, W_S, F_VM)
        else:
            # Cost trays assuming a partial condenser
            N_T = Design['Actual stages'] - 1
            Di = Design['Diameter']
            F_TM = self._F_TM_function(Di)
            Cost['Trays'] = design.compute_purchase_cost_of_trays(N_T, Di, F_TT, F_TM)
            
            # Cost vessel assuming T < 800 F
            W = Design['Weight'] # in lb
            L = Design['Height']*3.28 # in ft
            Cost['Tower'] = design.compute_purchase_cost_of_tower(Di, L, W, F_VM)
        self._simulate_components()
    
    def _plot_stages(self):
        """Plot stages, graphical aid line, and equilibrium curve. The plot does not include operating lines nor a legend."""
        vap, liq = self.outs
        if not hasattr(self, '_x_stages'):
            raise RuntimeError('cannot plot stages without running McCabe Thiele binary distillation')
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
            T[n], y = solve_Ty(np.array([xi, 1-xi]), P)
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
        y_top, x_bot = self._get_y_top_and_x_bot()
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
