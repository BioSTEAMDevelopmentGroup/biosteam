# -*- coding: utf-8 -*-
"""
General functional algorithms based on MESH equations to solve multi-stage.
"""
import numpy as np
from numba import njit
from math import inf
from scipy.optimize import minimize

# %% Diagonal matrix solution

@njit(cache=True)
def solve_TDMA(a, b, c, d): # Tridiagonal matrix solver
    """
    Solve a tridiagonal matrix using Thomas' algorithm.
    
    http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
    http://www.cfd-online.com/Wiki/Tridiagonal_matrix_algorithm_-_TDMA_(Thomas_algorithm)
    
    Notes
    -----
    `a` array starts from a1 (not a0).
    
    """
    n = d.shape[0] - 1 # number of equations minus 1
    for i in range(n):
        inext = i + 1
        m = a[i] / b[i]
        b[inext] = b[inext] - m * c[i] 
        d[inext] = d[inext] - m * d[i]
        
    b[n] = d[n] / b[n]
    for i in range(n-1, -1, -1):
        b[i] = (d[i] - c[i] * b[i+1]) / b[i]
    return b

@njit(cache=True)
def solve_LBDMA(a, b, d): # Left bidiagonal matrix solver
    """
    Solve a left bidiagonal matrix using a reformulation of Thomas' algorithm.
    """
    n = d.shape[0] - 1 # number of equations minus 1
    for i in range(n):
        inext = i + 1
        m = a[i] / b[i]
        d[inext] = d[inext] - m * d[i]
    
    b[n] = d[n] / b[n]

    for i in range(n-1, -1, -1):
        b[i] = d[i] / b[i]
    return b

@njit(cache=True)
def solve_RBDMA_1D_careful(b, c, d):
    n = d.shape[0] - 1 # number of equations minus 1
    bn = b[n]
    dn = d[n]
    if bn == 0:
        if dn == 0:
            b[n] = 0
        else:
            b[n] = inf
    else:
        b[n] = d[n] / b[n]

    for i in range(n-1, -1, -1):
        bi = b[i]
        num = d[i] - c[i] * b[i+1]
        if bi == 0:
            if num == 0:
                b[i] = 0
            else:
                b[i] = inf
        else:
            b[i] = num / bi
    return b

@njit(cache=True)
def solve_RBDMA(b, c, d): # Right bidiagonal matrix solver
    """
    Solve a right bidiagonal matrix using a reformulation of Thomas' algorithm.
    """
    n = d.shape[0] - 1 # number of equations minus 1
    b[n] = d[n] / b[n]

    for i in range(n-1, -1, -1):
        b[i] = (d[i] - c[i] * b[i+1]) / b[i]
    return b

# %% Bulk flow solution

@njit(cache=True)
def hot_start_top_flow_rates(
        bottom_flows, phase_ratios, stage_index, top_feed_flows,
        bottom_feed_flows, asplit_1, bsplit_1, N_stages,
    ):
    """
    Solve a-phase flow rates for a single component across 
    equilibrium stages with side draws. 

    Parameters
    ----------
    bottom_flows : Iterable[1d array]
        Bottom flow rates by stages.
    phase_ratios : 1d array
        Phase ratios by stage. The phase ratio for a given stage is 
        defined as F_a / F_b; where F_a and F_b are the flow rates 
        of phase a (extract or vapor) and b (raffinate or liquid) leaving the stage 
        respectively.
    stage_index : 1d array
        Stage index for phase ratios.
    top_feed_flows : Iterable[1d array]
        Top flow rates of all components fed across stages. Shape should be 
        (N_stages, N_chemicals).
    bottom_feed_flows : Iterable [1d array]
        Bottom flow rates of all components fed across stages. Shape should be 
        (N_stages, N_chemicals).
    asplit_1 : 1d array
        Side draw split from phase a minus 1 by stage.
    bsplit_1 : 1d array
        Side draw split from phase b minus 1 by stage.

    Returns
    -------
    flow_rates_a: 2d array
        Flow rates of phase a with stages by row and components by column.

    """
    d = top_feed_flows.copy()
    b = d.copy()
    c = d.copy()
    for i in range(N_stages): 
        c[i] = asplit_1[i]
        b[i] = 1
    for n in range(stage_index.size):
        i = stage_index[n]
        B = phase_ratios[n]
        if B <= 1e-32:
            b[i] = inf
        else:
            b[i] += 1 / B 
        if i == 0:
            d[i] += bottom_feed_flows[i]
        else:
            d[i] += bottom_feed_flows[i] - bottom_flows[i - 1] * bsplit_1[i - 1]
    return solve_RBDMA(b, c, d)

@njit(cache=True)
def hot_start_bottom_flow_rates(
        top_flows, phase_ratios, stage_index, top_feed_flows,
        bottom_feed_flows, asplit_1, bsplit_1, N_stages
    ):
    """
    Solve a-phase flow rates for a single component across 
    equilibrium stages with side draws. 

    Parameters
    ----------
    bottom_flows : Iterable[1d array]
        Bottom flow rates by stages.
    phase_ratios : 1d array
        Phase ratios by stage. The phase ratio for a given stage is 
        defined as F_a / F_b; where F_a and F_b are the flow rates 
        of phase a (extract or vapor) and b (raffinate or liquid) leaving the stage 
        respectively.
    stage_index : 1d array
        Stage index for phase ratios.
    top_feed_flows : Iterable[1d array]
        Top flow rates of all components fed across stages. Shape should be 
        (N_stages, N_chemicals).
    bottom_feed_flows : Iterable [1d array]
        Bottom flow rates of all components fed across stages. Shape should be 
        (N_stages, N_chemicals).
    asplit_1 : 1d array
        Side draw split from phase a minus 1 by stage.
    bsplit_1 : 1d array
        Side draw split from phase b minus 1 by stage.

    Returns
    -------
    flow_rates_a: 2d array
        Flow rates of phase a with stages by row and components by column.

    """
    d = bottom_feed_flows.copy()
    b = d.copy()
    a = d.copy()
    for i in range(N_stages): 
        a[i] = bsplit_1[i]
        b[i] = 1
    last_stage = N_stages - 1
    for n in range(stage_index.size):
        i = stage_index[n]
        b[i] += phase_ratios[n]
        if i == last_stage:
            d[i] += top_feed_flows[i]
        else:
            d[i] += top_feed_flows[i] - top_flows[i + 1] * asplit_1[i + 1]
    return solve_LBDMA(a, b, d)

# %% Component flow rate solution

@njit(cache=True)
def bottom_flow_rates(
        stripping_factors, 
        feed_flows,
        asplit_1,
        bsplit_1,
        N_stages,
    ):
    """
    Solve a-phase flow rates for a single component across equilibrium stages with side draws. 

    Parameters
    ----------
    stripping_factors : Iterable[1d array]
        The ratio of component flow rates in phase a (extract or vapor) over
        the component flow rates in phase a (raffinate or liquid). 
    feed_flows : Iterable[1d array]
        Flow rates of all components fed across stages. Shape should be 
        (N_stages, N_chemicals).
    asplit_1 : 1d array
        Side draw split from phase a minus 1 by stage.
    bsplit_1 : 1d array
        Side draw split from phase b minus 1 by stage.

    Returns
    -------
    flow_rates_a : 2d array
        Flow rates of phase a with stages by row and components by column.

    """
    b = 1. + stripping_factors
    c = np.expand_dims(asplit_1[1:], -1) * stripping_factors[1:]
    d = feed_flows.copy()
    a = bsplit_1
    return solve_TDMA(a, b, c, d) 

@njit(cache=True)
def liquid_compositions(
        bulk_vapor_flow_rates,
        bulk_liquid_flow_rates,
        partition_coefficients, 
        feed_flows,
        asplit_1,
        bsplit_1,
        N_stages,
    ):
    """
    Solve a-phase flow rates for a single component across equilibrium stages with side draws. 

    Parameters
    ----------
    stripping_factors : Iterable[1d array]
        The ratio of component flow rates in phase a (extract or vapor) over
        the component flow rates in phase a (raffinate or liquid). 
    feed_flows : Iterable[1d array]
        Flow rates of all components fed across stages. Shape should be 
        (N_stages, N_chemicals).
    asplit_1 : 1d array
        Side draw split from phase a minus 1 by stage.
    bsplit_1 : 1d array
        Side draw split from phase b minus 1 by stage.

    Returns
    -------
    flow_rates_a : 2d array
        Flow rates of phase a with stages by row and components by column.

    """
    KV = bulk_vapor_flow_rates * partition_coefficients
    b = bulk_liquid_flow_rates + KV
    c = np.expand_dims(asplit_1[1:], -1) * KV[1:]
    d = feed_flows.copy()
    a = bsplit_1 * bulk_liquid_flow_rates
    return solve_TDMA(a, b, c, d)

@njit(cache=True)
def bottom_flows_mass_balance(
        top_flows, feed_flows, asplit_left, bsplit_left,
        N_stages
    ):
    bottom_flows = 0 * top_flows
    bottom_flows[0] = feed_flows[0] + top_flows[1] * asplit_left[1] - top_flows[0]
    for i in range(1, N_stages-1):
        bottom_flows[i] = (
            feed_flows[i] + bsplit_left[i-1] * bottom_flows[i-1] + 
            top_flows[i+1] * asplit_left[i+1] - top_flows[i]
        )
    bottom_flows[-1] = feed_flows[-1] + bsplit_left[-2] * bottom_flows[-2] - top_flows[-1]
    return bottom_flows


@njit(cache=True)
def top_flows_mass_balance(
        bottom_flows, feed_flows, asplit_left, bsplit_left, N_stages
    ):
    top_flows = 0 * bottom_flows
    top_flows[-1] = feed_flows[-1] + bsplit_left[-2] * bottom_flows[-2] - bottom_flows[-1]
    for i in range(N_stages-2, 0, -1):
        top_flows[i] = (
            feed_flows[i] + bsplit_left[i-1] * bottom_flows[i-1] + 
            top_flows[i+1] * asplit_left[i+1] - bottom_flows[i]
        )
    top_flows[0] = feed_flows[0] + top_flows[1] * asplit_left[1] - bottom_flows[0]
    return top_flows

# %% Energy balance solution

@njit(cache=True)
def vapor_flow_rates(
        L, V, hl, hv, asplit_1, asplit_left, bsplit_left, 
        N_stages, specification_index, H_feeds
    ):
    # hV1*L1*dB1 - hv2*L2*dB2 = Q1 + (hl0 * L0 + hv2 * V2) - (hl1 * L1 + hv1 * V1)
    # hV1*L1*dB1 - hv2*L2*dB2 = Q1 + (hl0 * L0 + hv2 * L2 * B2) - (hl1 * L1 + hv1 * L1 * B1)
    # hV1*L1*dB1 + hv1 * L1 * B1 - hv2*L2*dB2 - hv2 * L2 * B2 = Q1 + (hl0 * L0) - (hl1 * L1)
    # hV1*L1*B1 - hv2*L2*B2 = Q1 + hl0 * L0 - hl1 * L1
    # hV1*V1 - hv2*V2 = Q1 + hl0 * L0 - hl1 * L1
    b = hv
    c = b[1:] * asplit_1[1:]
    Hl_out = hl * L
    Hl_in = (Hl_out * bsplit_left)[:-1]
    d = H_feeds - Hl_out
    d[1:] += Hl_in
    raise NotImplementedError('solution not implemented yet')
    for i, j in enumerate(specification_index):
        b[j] = 1
        d[j] = 0
        jlast = j - 1
        if jlast > 0: c[jlast] = 0
        try: c[j] = 0
        except: pass
    return solve_RBDMA(b, c, d)

@njit(cache=True)
def phase_ratio_departures(
        L, V, hl, hv, asplit_1, asplit_left, bsplit_left, 
        N_stages,
        specification_index,
        H_feeds
    ):
    # hV1*L1*dB1 - hv2*L2*dB2 = Q1 + H_in - H_out
    b = hv * L
    c = b[1:] * asplit_1[1:]
    Hl_out = hl * L
    Hv_out = hv * V
    d = H_feeds - Hl_out - Hv_out
    Hl_in = (Hl_out * bsplit_left)[:-1]
    Hv_in = (Hv_out * asplit_left)[1:]
    d[1:] += Hl_in
    d[:-1] += Hv_in
    for i, j in enumerate(specification_index):
        b[j] = 1
        d[j] = 0
        jlast = j - 1
        if jlast > 0: c[jlast] = 0
        try: c[j] = 0
        except: pass
    return solve_RBDMA(b, c, d)

@njit(cache=True)
def temperature_departures(Cv, Cl, Hv, Hl, asplit_left, bsplit_left,
                           N_stages, H_feeds):
    # ENERGY BALANCE
    # C1dT1 - Cv2*dT2 - Cl0*dT0 = Q1 - H_out + H_in
    b = (Cv + Cl)
    a = -(Cl * bsplit_left)
    c = -(Cv * asplit_left)[1:]
    d = H_feeds - Hl - Hv
    d[1:] += (Hl * bsplit_left)[:-1]
    d[:-1] += (Hv * asplit_left)[1:]
    return solve_TDMA(a, b, c, d)

def get_neighbors(index=None, all_index=None, missing=None, size=None):
    if size is not None:
        all_index = set(range(size))
    elif all_index is not None:
        all_index = set(all_index)
    if sum([i is None for i in (index, all_index, missing)]) > 1:
        raise ValueError('at least two arguments must be given')
    if missing is None: 
        missing = all_index.difference(index)
    else:
        missing = set(missing)
    if index is None:
        index_set = all_index.difference(missing)
    else:
        index_set = set(index)
    if all_index is None:
        all_index = index | missing
    size = len(all_index)
    neighbors = []
    for i in missing:
        lb = i
        while lb > -1:
            lb -= 1
            if lb in index_set: break
        ub = i
        while ub < size:
            ub += 1
            if ub in index_set: break
        if ub == size:
            neighbors.append(
                (i, (lb,))
            )
        elif lb == -1:
            neighbors.append(
                (i, (ub,))
            )
        else:
            neighbors.append(
                (i, (lb, ub))
            )
    return neighbors

def expand(values, index, size):
    new_values = np.zeros(size)
    new_values[index] = values
    return new_values

def fillmissing(all_neighbors, values):
    for i, neighbors in all_neighbors:
        if len(neighbors) == 2:
            lb, ub = neighbors
            lb_distance = i - lb
            ub_distance = ub - i
            sum_distance = lb_distance + ub_distance
            wlb = ub_distance / sum_distance
            wub = lb_distance / sum_distance
            x = wlb * values[lb] + wub * values[ub]
            values[i] = x
        else:
            values[i] = values[neighbors[0]]
    return values

# %% Wang-Henke 

#  Energy balance
def solve_vapor_flow_rates():
    pass

# Material balance
def solve_liquid_compositions():
    pass


# %% Russel's inside-out algorithm (not ready)

@njit(cache=True)
def omega_approx(y, K):
    y_over_K = (y / K)
    return y_over_K / y_over_K.sum()

@njit(cache=True)
def Kb_init(y, K):
    omega = omega_approx(y, K)
    return np.exp((omega * np.log(K)).sum(axis=1))

@njit(cache=True)
def Kb_iter(alpha, x):
    return 1 / (alpha * x).sum(axis=1)

@njit(cache=True)
def alpha_approx(K, Kb):
    return K / Kb

@njit(cache=True)
def fit(x, y):
    xmean = x.mean()
    ymean = y.mean()
    xxmean = x - xmean
    m = (xxmean * (y - ymean)).sum() / (xxmean * xxmean).sum()
    b = ymean - m * xmean
    return m, b

@njit(cache=True)
def fit_partition_model(T, Kb):
    # TODO: update to use adjacent stages
    x = 1 / T
    y = np.log(Kb)
    m, b = fit(x, y)
    M = y.copy()
    M[:] = m
    B = y - M * x
    return M, B

@njit(cache=True)
def h_approx(T, m, b):
    return m * T + b

@njit(cache=True)
def T_approx(Kb, m, b):
    return m / (np.log(Kb) - b)

def solve_inside_loop(
        lnSb, Sb, Sb_index, top_flows, K, B, T, hv, hl, feed_flows,
        asplit_1, bsplit_1, asplit_left, bsplit_left,
        N_stages, specification_index, N_chemicals, H_feeds
    ):
    correct_stages = np.ones(N_stages, bool)
    args = (
        *inside_loop_args(top_flows, K, T, hv, hl),
        feed_flows, asplit_1, bsplit_1,
        asplit_left, bsplit_left, N_stages, 
        specification_index, N_chemicals, H_feeds,
        correct_stages, Sb, Sb_index,
    )
    # result = root(inside_loop, KB.flatten(), 
    #               options=dict(ftol=1e-6), args=args)
    # print(result.x)
    # print(KB_new)
    return minimize(inside_loop_enthalpy_error, lnSb, method='CG', args=args)
   
@njit(cache=True)
def inside_loop_args(
        top_flows, K, T, hv, hl
    ):
    dummy = top_flows.sum(axis=1)
    dummy[dummy == 0] = 1
    y = top_flows / np.expand_dims(dummy, -1)
    Kb = Kb_init(y, K)
    Kb_coef = fit_partition_model(T, Kb)
    hv_coef = fit(T, hv)
    hl_coef = fit(T, hl)
    alpha = alpha_approx(K, np.expand_dims(Kb, -1))
    return alpha, Kb_coef, hv_coef, hl_coef
    
# @njit(cache=True)
def inside_loop_enthalpy_error(
        lnSb, alpha, Kb_coef, hv_coef, hl_coef, 
        feed_flows, asplit_1, bsplit_1,
        asplit_left, bsplit_left, N_stages, 
        specification_index, N_chemicals, H_feeds,
        correct_stages, Sb, Sb_index,
    ):
    safe = Sb_index is None
    if safe: 
        Sb = np.exp(lnSb).reshape([N_stages, N_chemicals])
    else:
        Sb[Sb_index] = np.exp(lnSb).reshape([len(Sb_index), N_chemicals])
    top_flows = top_flow_rates(
        Sb, 
        feed_flows,
        asplit_1,
        bsplit_1,
        N_stages,
        safe,
    )
    bottom_flows = mass_balance(
        top_flows, feed_flows, asplit_left, bsplit_left, 
        correct_stages.copy(), N_stages, N_chemicals
    )
    top_flows_net = top_flows.sum(axis=1)
    bottom_flows_net = bottom_flows.sum(axis=1)
    dummy = bottom_flows_net.copy()
    dummy[dummy == 0] = 1e-12
    x = bottom_flows / dummy[:, np.newaxis]
    dummy = top_flows_net.copy()
    dummy[dummy == 0] = 1e-12
    Kb = Kb_iter(alpha, x)
    T = T_approx(Kb, *Kb_coef)
    print(T)
    hv = h_approx(T, *hv_coef)
    hl = h_approx(T, *hl_coef)
    print(hl)
    Hl_out = hl * bottom_flows_net
    Hv_out = hv * top_flows_net
    d = H_feeds - Hl_out - Hv_out
    Hl_in = (Hl_out * bsplit_left)[:-1]
    Hv_in = (Hv_out * asplit_left)[1:]
    d[1:] += Hl_in
    d[:-1] += Hv_in
    for i in specification_index: d[i] = 0
    print(np.abs(d).sum())
    return np.abs(d).sum() 
    
    