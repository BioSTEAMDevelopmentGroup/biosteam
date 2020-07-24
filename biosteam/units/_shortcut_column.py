# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from ._binary_distillation import BinaryDistillation
import flexsolve as flx
from thermosteam.exceptions import InfeasibleRegion
from thermosteam.equilibrium import DewPoint, BubblePoint
import numpy as np

__all__ = ('ShortcutColumn',)

# %% Functions

@flx.njitable(cache=True)
def geometric_mean(a, b):
    return (a * b) ** 0.5

@flx.njitable(cache=True)
def compute_mean_volatilities_relative_to_heavy_key(K_distillate, K_bottoms, HK_index):
    alpha_distillate = K_distillate / K_distillate[HK_index]
    alpha_bottoms = K_bottoms / K_bottoms[HK_index]
    alpha_mean = geometric_mean(alpha_distillate,
                                alpha_bottoms)
    return alpha_mean

@flx.njitable(cache=True)
def compute_partition_coefficients(y, x):
    x[x <= 1e-16] = 1e-16
    return y / x

@flx.njitable(cache=True)
def compute_distillate_recoveries_Hengsteback_and_Gaddes(d_Lr, b_Hr,
                                                         alpha_mean,
                                                         LHK_index):
    LK_index = LHK_index[0]
    alpha_LK = alpha_mean[LK_index]
    A_dummy = (1. - b_Hr) / b_Hr
    A = np.log10(A_dummy)
    B = np.log10(d_Lr / (1. - d_Lr) / A_dummy) / np.log10(alpha_LK)
    dummy = 10.**A * alpha_mean**B
    distillate_recoveries = dummy / (1. + dummy)
    distillate_recoveries[LHK_index] = [d_Lr, 1. - b_Hr]
    distillate_recoveries[distillate_recoveries < 1e-12] = 0.
    return distillate_recoveries

@flx.njitable(cache=True)
def compute_minimum_theoretical_stages_Fenske(LHK_distillate, LHK_bottoms, alpha_LK):
    LK, HK = LHK_distillate
    LHK_ratio_distillate = LK / HK
    LK, HK = LHK_bottoms
    HLK_ratio_bottoms = HK / LK
    N = np.log10(LHK_ratio_distillate * HLK_ratio_bottoms) / np.log10(alpha_LK)
    return N

@flx.njitable(cache=True)
def objective_function_Underwood_constant(theta, q, z_f, alpha_mean):
    return (alpha_mean * z_f / (alpha_mean - theta)).sum() - 1.0 + q

@flx.njitable(cache=True)
def compute_minimum_reflux_ratio_Underwood(alpha_mean, z_d, theta):
    Rm = (alpha_mean * z_d / (alpha_mean - theta)).sum() - 1.0
    return Rm

@flx.njitable(cache=True)
def compute_theoretical_stages_Gilliland(Nm, Rm, R):
    X = (R - Rm) / (R + 1.)
    Y = 1. - np.exp((1. + 54.4*X) / (11. + 117.2*X) * (X - 1.) / X**0.5)
    N = (Y + Nm) / (1. - Y)
    return np.ceil(N)

@flx.njitable(cache=True)
def compute_feed_stage_Kirkbride(N, B, D,
                                 feed_HK_over_LK,
                                 z_LK_bottoms,
                                 z_HK_distillate):
    m_over_p = (B/D * feed_HK_over_LK * (z_LK_bottoms / z_HK_distillate)**2.) ** 0.206
    return np.floor(N / (m_over_p + 1.))


# %%

class ShortcutColumn(BinaryDistillation,
                     new_graphics=False):
    r"""
    Create a multicomponent distillation column that relies on the
    Fenske-Underwood-Gilliland method to solve for the theoretical design
    of the distillation column and the separation of non-keys [1]_.The Murphree
    efficiency (i.e. column efficiency) is based on the modified O'Connell
    correlation [2]_. The diameter is based on tray separation and flooding 
    velocity [1]_ [3]_. Purchase costs are based on correlations compiled by
    Warren et. al. [4]_.

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
        User enforced minimum reflux ratio. If the actual minimum reflux ratio is less than `Rmin`, this enforced value is ignored. Defaults to 0.6.
    specification="Composition" : "Composition" or "Recovery"
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
    >>> from biosteam.units import ShortcutColumn
    >>> from biosteam import Stream, settings
    >>> settings.set_thermo(['Water', 'Methanol', 'Glycerol'])
    >>> feed = Stream('feed', flow=(80, 100, 25))
    >>> bp = feed.bubble_point_at_P()
    >>> feed.T = bp.T # Feed at bubble point T
    >>> D1 = ShortcutColumn('D1', ins=feed,
    ...                     outs=('distillate', 'bottoms_product'),
    ...                     LHK=('Methanol', 'Water'),
    ...                     y_top=0.99, x_bot=0.01, k=2,
    ...                     is_divided=True)
    >>> D1.simulate()
    >>> # See all results
    >>> D1.show(T='degC', P='atm', composition=True)
    ShortcutColumn: D1
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
    Distillation                                    Units       D1
    Cooling water       Duty                        kJ/hr -7.9e+06
                        Flow                      kmol/hr  5.4e+03
                        Cost                       USD/hr     2.64
    Low pressure steam  Duty                        kJ/hr 1.43e+07
                        Flow                      kmol/hr      368
                        Cost                       USD/hr     87.5
    Design              Theoretical feed stage                   8
                        Theoretical stages                      16
                        Minimum reflux              Ratio     1.06
                        Reflux                      Ratio     2.12
                        Rectifier stages                        13
                        Stripper stages                         26
                        Rectifier height               ft     31.7
                        Stripper height                ft     50.9
                        Rectifier diameter             ft     4.53
                        Stripper diameter              ft     3.67
                        Rectifier wall thickness       in    0.312
                        Stripper wall thickness        in    0.312
                        Rectifier weight               lb 6.46e+03
                        Stripper weight                lb 7.98e+03
    Purchase cost       Rectifier trays               USD 1.52e+04
                        Stripper trays                USD 2.02e+04
                        Rectifier tower               USD 8.44e+04
                        Stripper tower                USD 1.01e+05
                        Condenser                     USD 4.17e+04
                        Boiler                        USD 2.99e+04
    Total purchase cost                               USD 2.92e+05
    Utility cost                                   USD/hr     90.1
    """
    line = 'Distillation'
    _ins_size_is_fixed = False
    _N_ins = 1
    _N_outs = 2     
     
    def _run(self):
        # Initial mass balance
        self._run_binary_distillation_mass_balance()
        
        # Initialize objects to calculate bubble and dew points
        vle_chemicals = self.feed.vle_chemicals
        reset_cache = self._vle_chemicals != vle_chemicals
        if reset_cache:
            self._dew_point = DewPoint(vle_chemicals, self.thermo)
            self._bubble_point = BubblePoint(vle_chemicals, self.thermo)
            self._IDs_vle = self._dew_point.IDs
            self._vle_chemicals = vle_chemicals
            
        # Setup light and heavy keys
        LHK = [i.ID for i in self.chemicals[self.LHK]]
        IDs = self._IDs_vle
        self._LHK_vle_index = np.array([IDs.index(i) for i in LHK], dtype=int)
        
        # Add temporary specification
        composition_spec = self.product_specification_format == 'Composition'
        if composition_spec:
            feed = self.feed
            distillate, bottoms = self.outs
            LK_index, HK_index = LHK_index = self._LHK_index
            LK_feed, HK_feed = feed.mol[LHK_index]
            self._Lr = distillate.mol[LK_index] / LK_feed
            self._Hr = bottoms.mol[HK_index] / HK_feed
            
        # Set starting point for solving column
        if reset_cache:
            self._add_trace_heavy_and_light_non_keys_in_products()
            distillate_recoveries = self._estimate_distillate_recoveries()
            self._distillate_recoveries = distillate_recoveries
            self._update_distillate_recoveries(distillate_recoveries)
        else:
            distillate_recoveries = self._distillate_recoveries
            lb = 1e-6; ub = 1 - 1e-6
            distillate_recoveries[distillate_recoveries < lb] = lb
            distillate_recoveries[distillate_recoveries > ub] = ub
            self._update_distillate_recoveries(distillate_recoveries)
        
        # Solve for new recoveries
        self._solve_distillate_recoveries()
        self._update_distillate_and_bottoms_temperature()
        
        # Remove temporary data
        if composition_spec: self._Lr = self._Hr = None
        
    def reset_cache(self):
        self._vle_chemicals = None

    def plot_stages(self):
        raise TypeError('cannot plot stages for shortcut column')
        
    def _design(self):
        self._run_FenskeUnderwoodGilliland()
        self._run_condenser_and_boiler()
        self._complete_distillation_column_design()
        
    def _run_FenskeUnderwoodGilliland(self):
        LHK_index = self._LHK_index
        alpha_mean = self._estimate_mean_volatilities_relative_to_heavy_key()
        LK_index = self._LHK_vle_index[0]
        alpha_LK = alpha_mean[LK_index]
        feed, = self.ins
        distillate, bottoms = self.outs
        Nm = compute_minimum_theoretical_stages_Fenske(distillate.mol[LHK_index],
                                                       bottoms.mol[LHK_index],
                                                       alpha_LK)
        theta = self._solve_Underwood_constant(alpha_mean, alpha_LK)
        IDs = self._IDs_vle
        z_d = distillate.get_normalized_mol(IDs)
        Rm = compute_minimum_reflux_ratio_Underwood(alpha_mean, z_d, theta)
        if Rm < self.Rmin: Rm = self.Rmin
        R = self.k * Rm
        N = compute_theoretical_stages_Gilliland(Nm, Rm, R)
        feed_HK, feed_LK = feed.mol[LHK_index]
        feed_HK_over_LK = feed_HK / feed_LK
        Bs = bottoms.imol[IDs]
        Ds = distillate.imol[IDs]
        B = Bs.sum()
        D = Ds.sum()
        LK_index, HK_index = LHK_index
        z_LK_bottoms = bottoms.mol[LK_index] / B
        z_HK_distillate = distillate.mol[HK_index] / D
        feed_stage = compute_feed_stage_Kirkbride(N, B, D, 
                                                  feed_HK_over_LK,
                                                  z_LK_bottoms,
                                                  z_HK_distillate)
        design = self.design_results
        design['Theoretical feed stage'] = N - feed_stage
        design['Theoretical stages'] = N
        design['Minimum reflux'] = Rm
        design['Reflux'] = R
        
    def _get_relative_volatilities_LHK(self):
        distillate, bottoms = self.outs
        LHK = self.LHK
        condensate = self.condensate
        K_light, K_heavy = distillate.get_molar_composition(LHK) / condensate.get_molar_composition(LHK)
        alpha_LHK_distillate = K_light/K_heavy
        
        boilup = self.boilup
        K_light, K_heavy = boilup.get_molar_composition(LHK) / bottoms.get_molar_composition(LHK)
        alpha_LHK_distillate = K_light/K_heavy
        alpha_LHK_bottoms = K_light/K_heavy
        
        return alpha_LHK_distillate, alpha_LHK_bottoms
        
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
    
    def _solve_Underwood_constant(self, alpha_mean, alpha_LK):
        q = self._get_feed_quality()
        z_f = self.ins[0].get_normalized_mol(self._IDs_vle)
        args = (q, z_f, alpha_mean)
        ub = np.inf
        lb = -np.inf
        bracket = flx.find_bracket(objective_function_Underwood_constant,
                                   1.0, alpha_LK, lb, ub, args)
        theta = flx.IQ_interpolation(objective_function_Underwood_constant,
                                     *bracket, args=args, checkiter=False,
                                     checkbounds=False)
        return theta
        
    def _add_trace_heavy_and_light_non_keys_in_products(self):
        distillate, bottoms = self.outs
        LNK_index = self._LNK_index
        HNK_index = self._HNK_index
        feed_mol = self.feed.mol
        LNK_mol = feed_mol[LNK_index]
        HNK_mol = feed_mol[HNK_index]
        bottoms.mol[LNK_index] = LNK_trace = 0.0001 * LNK_mol
        distillate.mol[LNK_index] = LNK_mol - LNK_trace
        distillate.mol[HNK_index] = HNK_trace = 0.0001 * HNK_mol
        bottoms.mol[HNK_index] = HNK_mol - HNK_trace
        
    def _estimate_mean_volatilities_relative_to_heavy_key(self):
        # Mean volatilities taken at distillate and bottoms product
        distillate, bottoms = self.outs
        dew_point = self._dew_point
        bubble_point = self._bubble_point
        IDs = self._IDs_vle
        z_distillate = distillate.get_normalized_mol(IDs)
        z_bottoms = bottoms.get_normalized_mol(IDs)
        dp = dew_point(z_distillate, P=self.P)
        bp = bubble_point(z_bottoms, P=self.P)
        K_distillate = compute_partition_coefficients(dp.z, dp.x)
        K_bottoms = compute_partition_coefficients(bp.y, bp.z)
        HK_index = self._LHK_vle_index[1]
        alpha_mean = compute_mean_volatilities_relative_to_heavy_key(K_distillate,
                                                                     K_bottoms,
                                                                     HK_index)
        return alpha_mean
        
    def _estimate_distillate_recoveries(self):
        # Use Hengsteback and Geddes equations
        alpha_mean = self._estimate_mean_volatilities_relative_to_heavy_key()
        return compute_distillate_recoveries_Hengsteback_and_Gaddes(self.Lr, self.Hr,
                                                                    alpha_mean,
                                                                    self._LHK_vle_index)
        
    def _update_distillate_recoveries(self, distillate_recoveries):
        feed = self.feed
        distillate, bottoms = self.outs
        IDs = self._IDs_vle
        feed_mol = feed.imol[IDs]
        distillate.imol[IDs] = distillate_mol = distillate_recoveries * feed_mol
        bottoms.imol[IDs] = feed_mol - distillate_mol
        
    def _solve_distillate_recoveries(self):
        distillate_recoveries = self._distillate_recoveries
        flx.aitken(self._recompute_distillate_recoveries,
                   distillate_recoveries, 1e-8, checkiter=False)
        
    def _recompute_distillate_recoveries(self, distillate_recoveries):
        if np.logical_or(distillate_recoveries > 1., distillate_recoveries < 0.).any():
            raise InfeasibleRegion('distillate composition')
        self._update_distillate_recoveries(distillate_recoveries)
        distillate_recoveries = self._estimate_distillate_recoveries()
        if hasattr(self, '_distillate_recoveries_hook'):
            self._distillate_recoveries_hook(self._IDs_vle, distillate_recoveries)
        self._distillate_recoveries = distillate_recoveries
        return distillate_recoveries
        
    
    