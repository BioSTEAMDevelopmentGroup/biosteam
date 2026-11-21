# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import pytest
import biosteam as bst
import numpy as np
from numpy.testing import assert_allclose
from biosteam.units.design_tools import MESH

def test_solve_block_tridiagonal_matrix():
    N, M = 3, 2
    rng = np.random.default_rng(0)
    A = np.array([0.5 * np.eye(M), 0.3 * np.eye(M)]) # Scalar sub-diagonal blocks
    B = rng.random((N, M, M)) + 4 * np.eye(M)
    C = rng.random((N - 1, M, M))
    d = rng.random((N, M))
    x = MESH.solve_block_tridiagonal_matrix(A.copy(), B.copy(), C.copy(), d)
    dense = MESH.create_block_tridiagonal_matrix(A, B, C)
    assert_allclose(dense @ x.flatten(), d.flatten(), rtol=1e-6, atol=1e-9)

def test_hot_start_flow_rates():
    N_stages, N_chemicals = 3, 2
    top_feed_flows = np.zeros((N_stages, N_chemicals))
    bottom_feed_flows = np.zeros((N_stages, N_chemicals))
    bottom_feed_flows[0] = [10., 5.]
    top_feed_flows[-1] = [2., 8.]
    neg_asplit = -np.ones(N_stages)
    neg_bsplit = -np.ones(N_stages)
    stage_index = np.arange(N_stages)
    phase_ratios = np.array([1., 1e-40, 2.])
    bottom_flows = np.ones((N_stages, N_chemicals))
    top_flows = MESH.hot_start_top_flow_rates(
        bottom_flows, phase_ratios, stage_index, top_feed_flows,
        bottom_feed_flows, neg_asplit, neg_bsplit, N_stages,
    )
    assert top_flows.shape == (N_stages, N_chemicals)
    assert np.isfinite(top_flows).all()
    bottom_flows = MESH.hot_start_bottom_flow_rates(
        top_flows, phase_ratios, stage_index, top_feed_flows,
        bottom_feed_flows, neg_asplit, neg_bsplit, N_stages,
    )
    assert bottom_flows.shape == (N_stages, N_chemicals)
    assert np.isfinite(bottom_flows).all()

def test_liquid_compositions_and_mass_balance():
    N_stages, N_chemicals = 3, 2
    V = np.array([[5.], [5.], [5.]])
    L = np.array([[5.], [5.], [5.]])
    K = np.array([[2., 0.5]] * N_stages)
    feed_flows = np.zeros((N_stages, N_chemicals))
    feed_flows[1] = [5., 5.]
    neg_split = -np.ones(N_stages)
    x = MESH.liquid_compositions(V, L, K, feed_flows, neg_split, neg_split, N_stages)
    assert x.shape == (N_stages, N_chemicals)
    assert np.isfinite(x).all()
    top_flows = np.array([[3., 1.], [2., 1.], [1., 1.]])
    feed_flows = np.array([[0., 0.], [10., 10.], [0., 0.]])
    split = np.ones(N_stages)
    bottom_flows = MESH.bottom_flows_mass_balance(
        top_flows, feed_flows, split, split, N_stages
    )
    # Overall mass balance: feed = top out of stage 0 + bottom out of last stage
    assert_allclose(
        top_flows[0] + bottom_flows[-1], feed_flows.sum(axis=0), rtol=1e-6
    )

def test_phase_ratio_departures():
    N_stages = 4
    L = np.array([10., 10., 10., 10.])
    V = np.array([10., 10., 10., 10.])
    hl = np.array([1., 1., 1., 1.])
    hv = np.array([2., 2., 2., 2.])
    neg_asplit = -np.ones(N_stages)
    asplit = np.ones(N_stages)
    bsplit = np.ones(N_stages)
    H_feeds = np.array([0., 5., 0., 0.])
    dB = MESH.phase_ratio_departures(
        L, V, hl, hv, neg_asplit, asplit, bsplit, N_stages,
        np.array([0]), H_feeds, False,
    )
    assert dB.shape == (N_stages,) and np.isfinite(dB).all()
    assert dB[0] == 0
    dB = MESH.phase_ratio_departures(
        L, V, hl, hv, neg_asplit, asplit, bsplit, N_stages,
        np.array([0]), H_feeds, True,
    )
    assert dB.shape == (N_stages,) and np.isfinite(dB).all()

def test_bulk_flow_rates_invalid_specification():
    N_stages = 2
    ones = np.ones(N_stages)
    zeros = np.zeros(N_stages)
    with pytest.raises(ValueError):
        MESH.bulk_vapor_and_liquid_flow_rates(
            ones, 2 * ones, -ones, -ones, zeros, zeros, N_stages,
            zeros, ones, 'XX', ones, 2.
        )

def test_interpolation_utilities():
    assert MESH.get_neighbors(index=[0, 3], size=4) == [(1, (0, 3)), (2, (0, 3))]
    neighbors = MESH.get_neighbors(index=[1], all_index=[0, 1, 2])
    assert sorted(neighbors) == [(0, (1,)), (2, (1,))]
    neighbors = MESH.get_neighbors(missing=[1], all_index=[0, 1, 2])
    assert neighbors == [(1, (0, 2))]
    with pytest.raises(ValueError):
        MESH.get_neighbors(index=[0])
    T, = MESH.interpolate_missing_variables([0, -1], 3, np.array([300., 400.]))
    assert_allclose(T, [300., 350., 400.])
    assert_allclose(MESH.expand(np.array([1., 2.]), [0, 2], 3), [1., 0., 2.])

if __name__ == '__main__':
    test_solve_block_tridiagonal_matrix()
    test_hot_start_flow_rates()
    test_liquid_compositions_and_mass_balance()
    test_phase_ratio_departures()
    test_bulk_flow_rates_invalid_specification()
    test_interpolation_utilities()
