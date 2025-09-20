# -*- coding: utf-8 -*-
"""
General functional algorithms based on MESH equations to solve multi-stage.
"""
import numpy as np
from numba import njit
from math import inf
from scipy.optimize import minimize, lsq_linear

# %% Diagonal matrix solution

@njit(cache=True)
def solve_tridiagonal_matrix(a, b, c, d): # Tridiagonal matrix solver
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
        b[inext] -= m * c[i] 
        d[inext] -= m * d[i]
    
    b[n] = d[n] / b[n]
    for i in range(n-1, -1, -1):
        b[i] = (d[i] - c[i] * b[i+1]) / b[i]
    return b

@njit(cache=True)
def solve_block_tridiagonal_matrix(A, B, C, d):
    """
    Solve a block tridiagonal matrix using Thomas' algorithm.
    
    http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
    http://www.cfd-online.com/Wiki/Tridiagonal_matrix_algorithm_-_TDMA_(Thomas_algorithm)
    
    Notes
    -----
    `a` array starts from a1 (not a0).
    
    """
    d = d.copy()
    n = d.shape[0] - 1 # number of equations minus 1
    for i in range(n):
        inext = i + 1
        M = np.linalg.solve(B[i], A[i])
        B[inext] -= M @ C[i] 
        d[inext] -= M @ d[i]
        
    d[n] = np.linalg.solve(B[n], d[n])
    for i in range(n-1, -1, -1):
        d[i] = np.linalg.solve(B[i], d[i] - C[i] @ d[i+1])
    return d
    
@njit(cache=True)
def solve_left_bidiagonal_matrix(a, b, d): # Left bidiagonal matrix solver
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
def solve_right_bidiagonal_matrix(b, c, d): # Right bidiagonal matrix solver
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
        bottom_feed_flows, neg_asplit, neg_bsplit, N_stages,
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
    neg_asplit : 1d array
        Negative of the fraction of phase a being sent to the next stage.
    neg_bsplit : 1d array
        Negative of the fraction of phase b being sent to the next stage.

    Returns
    -------
    flow_rates_a: 2d array
        Flow rates of phase a with stages by row and components by column.

    """
    d = top_feed_flows.copy()
    b = d.copy()
    c = d.copy()
    for i in range(N_stages): 
        c[i] = neg_asplit[i]
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
            d[i] += bottom_feed_flows[i] - bottom_flows[i - 1] * neg_bsplit[i - 1]
    return solve_right_bidiagonal_matrix(b, c, d)

@njit(cache=True)
def hot_start_bottom_flow_rates(
        top_flows, phase_ratios, stage_index, top_feed_flows,
        bottom_feed_flows, neg_asplit, neg_bsplit, N_stages
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
    neg_asplit : 1d array
        Negative of the fraction of phase a being sent to the next stage.
    neg_bsplit : 1d array
        Negative of the fraction of phase b being sent to the next stage.

    Returns
    -------
    flow_rates_a: 2d array
        Flow rates of phase a with stages by row and components by column.

    """
    d = bottom_feed_flows.copy()
    b = d.copy()
    a = d.copy()
    for i in range(N_stages): 
        a[i] = neg_bsplit[i]
        b[i] = 1
    last_stage = N_stages - 1
    for n in range(stage_index.size):
        i = stage_index[n]
        b[i] += phase_ratios[n]
        if i == last_stage:
            d[i] += top_feed_flows[i]
        else:
            d[i] += top_feed_flows[i] - top_flows[i + 1] * neg_asplit[i + 1]
    return solve_left_bidiagonal_matrix(a, b, d)

@njit(cache=True)
def create_block_tridiagonal_matrix(As, Bs, Cs):
    """
    Create a dense block tridiagonal matrix.

    Parameters
    ----------
    As : list of (m x m) arrays
        Lower-diagonal blocks, length n-1
    Bs : list of (m x m) arrays
        Diagonal blocks, length n
    Cs : list of (m x m) arrays
        Upper-diagonal blocks, length n-1

    Returns
    -------
    M : (n*m, n*m) ndarray
        Dense block tridiagonal matrix
    """
    n = len(Bs)
    m = Bs[0].shape[0]
    n_vars = n * m
    M = np.zeros((n_vars, n_vars))

    # initialize first row slice
    row_start = 0
    row_end = row_start + m
    row_end_next = row_end + m
    row_slice = slice(row_start, row_end)
    last = n-1
    for i in range(n):
        M[row_slice, row_slice] = Bs[i]  # diagonal
        if i == last: break
        next_slice = slice(row_end, row_end_next)
        M[row_slice, next_slice] = Cs[i]  # upper
        M[next_slice, row_slice] = As[i]  # lower
        # prepare for next iteration
        row_slice = next_slice
        row_start = row_end
        row_end = row_end_next
        row_end_next += m
    return M

# %% Component flow rate solution

@njit(cache=True)
def bottom_flow_rates(
        stripping_factors, 
        feed_flows,
        neg_asplit,
        neg_bsplit,
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
    neg_asplit : 1d array
        Negative of the fraction of phase a being sent to the next stage.
    neg_bsplit : 1d array
        Negative of the fraction of phase b being sent to the next stage.

    Returns
    -------
    flow_rates_a : 2d array
        Flow rates of phase a with stages by row and components by column.

    """
    b = 1. + stripping_factors
    c = np.expand_dims(neg_asplit[1:], -1) * stripping_factors[1:]
    d = feed_flows.copy()
    a = neg_bsplit
    return solve_tridiagonal_matrix(a, b, c, d) 

@njit(cache=True)
def liquid_compositions(
        bulk_vapor_flow_rates,
        bulk_liquid_flow_rates,
        partition_coefficients, 
        feed_flows,
        neg_asplit,
        neg_bsplit,
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
    neg_asplit : 1d array
        Negative of the fraction of phase a being sent to the next stage.
    neg_bsplit : 1d array
        Negative of the fraction of phase b being sent to the next stage.

    Returns
    -------
    flow_rates_a : 2d array
        Flow rates of phase a with stages by row and components by column.

    """
    KV = bulk_vapor_flow_rates * partition_coefficients
    b = bulk_liquid_flow_rates + KV
    c = np.expand_dims(neg_asplit[1:], -1) * KV[1:]
    d = feed_flows.copy()
    a = np.expand_dims(neg_bsplit, -1) * bulk_liquid_flow_rates
    return solve_tridiagonal_matrix(a, b, c, d)

@njit(cache=True)
def bottom_flows_mass_balance(
        top_flows, feed_flows, asplit, bsplit,
        N_stages
    ):
    bottom_flows = 0 * top_flows
    bottom_flows[0] = feed_flows[0] + top_flows[1] * asplit[1] - top_flows[0]
    for i in range(1, N_stages-1):
        bottom_flows[i] = (
            feed_flows[i] + bsplit[i-1] * bottom_flows[i-1] + 
            top_flows[i+1] * asplit[i+1] - top_flows[i]
        )
    bottom_flows[-1] = feed_flows[-1] + bsplit[-2] * bottom_flows[-2] - top_flows[-1]
    return bottom_flows


@njit(cache=True)
def top_flows_mass_balance(
        bottom_flows, feed_flows, asplit, bsplit, N_stages
    ):
    top_flows = 0 * bottom_flows
    top_flows[-1] = feed_flows[-1] + bsplit[-2] * bottom_flows[-2] - bottom_flows[-1]
    for i in range(N_stages-2, 0, -1):
        top_flows[i] = (
            feed_flows[i] + bsplit[i-1] * bottom_flows[i-1] + 
            top_flows[i+1] * asplit[i+1] - bottom_flows[i]
        )
    top_flows[0] = feed_flows[0] + top_flows[1] * asplit[1] - bottom_flows[0]
    return top_flows

# %% Energy balance solution

@njit(cache=True)
def bulk_vapor_and_liquid_flow_rates(
        hl, hv, 
        neg_asplit, neg_bsplit, 
        top_split, bottom_split, 
        N_stages, H_feeds, F_feeds,
        variables,
        values,
        bulk_feed,
    ):
    # Equations:
    # - hl0*L0 + hv1*V1 + hl1*L1 - hv2*V2 = Q1
    # - L0 + V1 + L1 - V2 = 0
    N_bulk = 2 * N_stages
    Q = np.zeros(N_bulk)
    L0 = -1 # L0
    V1 = 0 # V1
    L1 = 1 # L1
    V2 = 2 # V2
    A = np.zeros((N_bulk, N_bulk))
    last = N_stages - 1
    for i in range(N_stages):
        variable = variables[i]
        value = values[i] 
        n = 2 * i
        m = n + 1
        ilast = i - 1
        inext = i + 1
        
        # Bulk material balance
        if L0 != -1: A[m, L0] = neg_bsplit[ilast]
        A[m, V1] = 1
        A[m, L1] = 1
        if i != last: A[m, V2] = neg_asplit[inext]
        Q[m] = F_feeds[i]
        
        # Energy balance
        if variable == 'Q':
            if L0 != -1: A[n, L0] = hl[ilast] * neg_bsplit[ilast]
            A[n, V1] = hv[i]
            A[n, L1] = hl[i]
            Q[n] = H_feeds[i] + value
            if i != last: A[n, V2] = hv[inext] * neg_asplit[inext]
        elif variable == 'B':
            A[n, V1] = -1
            A[n, L1] = value
        elif variable == 'F':
            A[n, L1] = 1
            Q[n] = bulk_feed * value
        else:
            raise ValueError('unexpected variable specification')
        L0 += 2
        V1 += 2
        L1 += 2
        V2 += 2
    VLs = np.linalg.solve(A, Q)
    VLs[VLs < 1e-16] = 1e-16
    return VLs[::2], VLs[1::2] # V, L

@njit(cache=True)
def phase_ratio_departures(
        L, V, hl, hv, neg_asplit, asplit, bsplit, 
        N_stages,
        specification_index,
        H_feeds,
        RF=False, # Reflux rate and bottoms flow rate specification
    ):
    # hV1*L1*dB1 - hv2*L2*dB2 = Q1 + H_in - H_out
    b = hv * L
    c = b[1:] * neg_asplit[1:]
    Hl_out = hl * L
    Hv_out = hv * V
    d = H_feeds - Hl_out - Hv_out
    Hl_in = (Hl_out * bsplit)[:-1]
    Hv_in = (Hv_out * asplit)[1:]
    d[1:] += Hl_in
    d[:-1] += Hv_in
    if RF:
        dB = d * 0
        dB[2:] = solve_left_bidiagonal_matrix(b[2:], c[1:], d[1:-1])
        return dB
    else:
        n_last = d.size - 1
        for j in specification_index:
            b[j] = 1
            d[j] = 0
            jlast = j - 1
            if jlast > 0: c[jlast] = 0
            if j < n_last: c[j] = 0
        return solve_right_bidiagonal_matrix(b, c, d)

@njit(cache=True)
def temperature_departures(Cv, Cl, Hv, Hl, asplit, bsplit,
                           N_stages, H_feeds):
    # ENERGY BALANCE
    # C1dT1 - Cv2*dT2 - Cl0*dT0 = Q1 - H_out + H_in
    b = (Cv + Cl)
    a = -(Cl * bsplit)
    c = -(Cv * asplit)[1:]
    d = H_feeds - Hl - Hv
    d[1:] += (Hl * bsplit)[:-1]
    d[:-1] += (Hv * asplit)[1:]
    return solve_tridiagonal_matrix(a, b, c, d)

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

    
    