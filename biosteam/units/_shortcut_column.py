# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 09:22:08 2020

@author: yoelr
"""
from math import log10, ceil, exp
from ._binary_distillation import BinaryDistillation
from flexsolve import wegstein, bisection
from thermosteam.equilibrium import DewPoint, BubblePoint

__all__ = ('ShortcutColumn',)

# %% Functions

def geometric_mean(a, b):
    return (a * b) ** (0.5)

def compute_mean_volatilities_relative_to_heavy_key(K_distillate, K_bottoms, HK_index):
    alpha_distillate = K_distillate / K_distillate[HK_index]
    alpha_bottoms = K_bottoms / K_bottoms[HK_index]
    alpha_mean = geometric_mean(alpha_distillate,
                                alpha_bottoms)
    return alpha_mean
    
def compute_partition_coefficients(y, x):
    return y / x
    
def compute_distillate_recoveries_Hengsteback_and_Gaddes(d_Lr, b_Hr,
                                                         alpha_mean,
                                                         LHK_index):
    LK_index = LHK_index[0]
    alpha_LK = alpha_mean[LK_index]
    A_dummy = (1 - b_Hr) / b_Hr
    A = log10(A_dummy)
    B = log10(d_Lr / (1 - d_Lr) / A_dummy) / log10(alpha_LK)
    dummy = 10**A * alpha_mean**B
    distillate_recoveries = dummy / (1 + dummy)
    distillate_recoveries[LHK_index] = [d_Lr, 1 - b_Hr]
    return distillate_recoveries

def compute_minimum_theoretical_stages_Fenske(LHK_distillate, LHK_bottoms, alpha_LK):
    LK, HK = LHK_distillate
    LHK_ratio_distillate = LK / HK
    LK, HK = LHK_bottoms
    HLK_ratio_bottoms = HK / LK
    N = log10(LHK_ratio_distillate * HLK_ratio_bottoms) / log10(alpha_LK)
    return N

def objective_function_Underwood_constant(theta, q, z_f,
                                          alpha_mean):
    return (alpha_mean * z_f / (alpha_mean - theta)).sum() - 1 + q

def compute_minimum_reflux_ratio_Underwood(alpha_mean, z_d, theta):
    Rm = (alpha_mean * z_d / (alpha_mean - theta)).sum() - 1
    return Rm

def compute_theoretical_stages_Gilliland(Nm, Rm, R):
    X = (R - Rm) / (R + 1)
    Y = 1 - exp((1 + 54.4*X) / (11 + 117.2*X) * (X - 1) / X**0.5)
    N = (Y + Nm) / (1 - Y)
    return ceil(N)

def compute_feed_stage_Kirkbride(N, B, D,
                                 feed_HK_over_LK,
                                 z_LK_bottoms,
                                 z_HK_distillate):
    m_over_p = B/D * feed_HK_over_LK * (z_LK_bottoms / z_HK_distillate)**2
    return ceil(N / (m_over_p + 1))
    
    

# %%

class ShortcutColumn(BinaryDistillation,
                     new_graphics=False):
    _N_ins = 1
    _N_outs = 2     
     
    def _run(self):
        # Initialize objects to calculate bubble and dew points
        feed = self.ins[0]
        equilibrium_chemicals = feed.equilibrium_chemicals
        self._dew_point = DewPoint(equilibrium_chemicals, self.thermo)
        self._bubble_point = BubblePoint(equilibrium_chemicals, self.thermo)
        self._IDs = IDs = self._dew_point.IDs
        LK, HK = self.LHK
        LK_index = IDs.index(LK)
        HK_index = IDs.index(HK)
        self._LHK_equilibrium_index = [LK_index, HK_index]
        
        # Set starting point for solving column
        self._run_binary_distillation_mass_balance()
        self._add_trace_heavy_and_light_non_keys_in_products()
        
        # Solve for new recoveries
        distillate_recoveries = self._solve_distillate_recoveries()
        self._update_distillate_recoveries(distillate_recoveries)
        self._update_distillate_and_bottoms_temperature()
        
    def _design(self):
        self._load_feed_flows()
        self._run_FenskeUnderwoodGilliland()
        self._run_condenser_and_boiler()
        self._complete_distillation_column_design()
        
    def _run_FenskeUnderwoodGilliland(self):
        LHK_index = self._LHK_index
        alpha_mean = self._estimate_mean_volatilities_relative_to_heavy_key()
        LK_index = self._LHK_equilibrium_index[0]
        alpha_LK = alpha_mean[LK_index]
        feed, = self.ins
        distillate, bottoms = self.outs
        Nm = compute_minimum_theoretical_stages_Fenske(distillate.mol[LHK_index],
                                                       bottoms.mol[LHK_index],
                                                       alpha_LK)
        theta = self._solve_Underwood_constant(alpha_mean, alpha_LK)
        IDs = self._IDs
        z_d = distillate.get_normalized_mol(IDs)
        Rm = compute_minimum_reflux_ratio_Underwood(alpha_mean, z_d, theta)
        if Rm < 0.1: Rm = 0.1
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
        feed, = self.ins
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
        z_f = self.ins[0].get_normalized_mol(self._IDs)
        theta = bisection(objective_function_Underwood_constant,
                          1, alpha_LK, args=(q, z_f, alpha_mean))
        return theta
        
    def _add_trace_heavy_and_light_non_keys_in_products(self):
        distillate, bottoms = self.outs
        LNK_index = self._LNK_index
        HNK_index = self._HNK_index
        bottoms.mol[LNK_index] = LNK_trace = 0.001 * distillate.mol[LNK_index]
        distillate.mol[LNK_index] -= LNK_trace
        distillate.mol[HNK_index] = HNK_trace = 0.001 * bottoms.mol[HNK_index]
        bottoms.mol[HNK_index] -= HNK_trace
        
    def _estimate_mean_volatilities_relative_to_heavy_key(self):
        # Mean volatilities taken at distillate and bottoms product
        distillate, bottoms = self.outs
        dew_point = self._dew_point
        bubble_point = self._bubble_point
        IDs = self._IDs
        z_distillate = distillate.get_normalized_mol(IDs)
        z_bottoms = bottoms.get_normalized_mol(IDs)
        # print(IDs, 'z_distillate', z_distillate, 'P', self.P, 'z_bottoms', z_bottoms)
        dp = dew_point(z_distillate, P=self.P)
        bp = bubble_point(z_bottoms, P=self.P)
        K_distillate = compute_partition_coefficients(dp.z, dp.x)
        K_bottoms = compute_partition_coefficients(bp.y, bp.z)
        HK_index = self._LHK_equilibrium_index[1]
        alpha_mean = compute_mean_volatilities_relative_to_heavy_key(K_distillate,
                                                                     K_bottoms,
                                                                     HK_index)
        return alpha_mean
        
    def _estimate_distillate_recoveries(self):
        # Use Hengsteback and Geddes equations
        alpha_mean = self._estimate_mean_volatilities_relative_to_heavy_key()
        feed, = self.ins
        distillate, bottoms = self.outs
        LHK_index = self._LHK_index
        LK_index, HK_index = LHK_index
        d_Lr = distillate.mol[LK_index] / feed.mol[LK_index]
        b_Hr = bottoms.mol[HK_index] / feed.mol[HK_index]
        LHK_equilibrium_index = self._LHK_equilibrium_index
        LK_index = LHK_equilibrium_index[0]
        return compute_distillate_recoveries_Hengsteback_and_Gaddes(d_Lr, b_Hr,
                                                                    alpha_mean,
                                                                    LHK_equilibrium_index)
        
    def _update_distillate_recoveries(self, distillate_recoveries):
        feed, = self.ins
        distillate, bottoms = self.outs
        IDs = self._IDs
        feed_mol = feed.imol[IDs]
        distillate.imol[IDs] = distillate_mol = distillate_recoveries * feed_mol
        bottoms.imol[IDs] = feed_mol - distillate_mol
        
    def _solve_distillate_recoveries(self):
        distillate_recoveries = self._estimate_distillate_recoveries()
        return wegstein(self._recompute_distillate_recovery,
                        distillate_recoveries, 1e-4)
    
    def _recompute_distillate_recovery(self, distillate_recoveries):
        self._update_distillate_recoveries(distillate_recoveries)
        return self._estimate_distillate_recoveries()
        

# if __name__ == '__main__':
#     from biosteam import Stream, settings, ShortcutColumn, BinaryDistillation
#     settings.set_thermo(['Water', 'Ethanol', 'Methanol', 'Glycerol'])
#     feed = Stream('feed', flow=(1, 80, 80, 1))
#     bp = feed.bubble_point_at_P()
#     feed.T = bp.T
    
#     D1 = ShortcutColumn('D1', ins=feed,
#                         outs=('distillate', 'bottoms_product'),
#                         LHK=('Methanol', 'Ethanol'),
#                         Lr=0.80, Hr=0.80, k=1.2,
#                         product_specification_format='Recovery',                        
#                         is_divided=True)
#     D1.simulate()
#     D1.show(T='degC', P='atm', composition=True)
#     print(D1.results())
#     print('\n')
#     D2 = BinaryDistillation('D2', ins=feed,
#                             outs=('distillate2', 'bottoms_product2'),
#                             LHK=('Methanol', 'Ethanol'),
#                             Lr=0.80, Hr=0.80, k=1.2,
#                             product_specification_format='Recovery',                        
#                             is_divided=True)
#     D2.simulate()
#     D2.show(T='degC', P='atm', composition=True)
#     print(D2.results())
    
    