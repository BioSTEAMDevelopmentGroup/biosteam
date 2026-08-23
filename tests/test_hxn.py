# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-, Yoel Cortes-Pena <yoelcortes@gmail.com>
# Copyright (C) 2026-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Tests for the heat exchanger network facility.
"""
import warnings
import pytest
import biosteam as bst
import numpy as np
from numpy.testing import assert_allclose

def build_system(N_columns=1):
    """Doctest system of HeatExchangerNetwork; `N_columns > 1` adds more
    ShortcutColumns so auxiliary heat exchangers have duplicate IDs."""
    bst.settings.set_thermo(['Water', 'Methanol', 'Glycerol'], cache=True)
    bst.main_flowsheet.set_flowsheet('test_hxn')
    units = []
    feeds = []
    for i in range(N_columns):
        feed = bst.Stream(f'feed{i}', flow=(8000, 100 * (i + 1), 25))
        D = bst.ShortcutColumn(f'D{i}', ins=feed, LHK=('Methanol', 'Water'),
                               y_top=0.99, x_bot=0.01, k=2, is_divided=True)
        H1 = bst.HXutility(f'D{i}_H1', ins=D.outs[1], T=300)
        H2 = bst.HXutility(f'D{i}_H2', ins=D.outs[0], T=300)
        units.extend([D, H1, H2])
        feeds.append(feed)
    feed2 = bst.Stream('feed_flash', flow=(10000, 1000, 10))
    F1 = bst.Flash('F1', ins=feed2, V=0.9, P=101325)
    HXN = bst.HeatExchangerNetwork('HXN', T_min_app=5.)
    sys = bst.System.from_units('sys', units=[*units, F1, HXN])
    return sys, HXN, feeds[0]

def network_results(HXN):
    return dict(
        heat=HXN.actual_heat_util_load,
        cool=HXN.actual_cool_util_load,
        Q=np.array([hx.Q for hx in HXN.new_HXs]),
        installed=HXN.installed_costs['Heat exchangers'],
    )

def assert_same_results(a, b, rtol=1e-6):
    for key in a:
        assert_allclose(a[key], b[key], rtol=rtol, err_msg=key)

def simulate_cached(sys, HXN):
    HXN.cache_network = True
    HXN_sys = HXN.HXN_sys
    with warnings.catch_warnings():
        warnings.simplefilter('error', RuntimeWarning)
        sys.simulate()
    assert HXN.HXN_sys is HXN_sys, 'cached network was not used'

def test_cache_network_matches_fresh_synthesis():
    sys, HXN, feed = build_system()
    sys.simulate()
    fresh = network_results(HXN)
    assert HXN.actual_heat_util_load < 0.9 * HXN.original_heat_util_load
    simulate_cached(sys, HXN)
    assert_same_results(network_results(HXN), fresh)

def test_cache_network_perturbed_feed():
    sys, HXN, feed = build_system()
    sys.simulate()
    feed.F_mass *= 1.01
    simulate_cached(sys, HXN)
    cached = network_results(HXN)
    HXN.cache_network = False
    sys.simulate()
    assert_same_results(cached, network_results(HXN), rtol=1e-3)

def test_cache_network_duplicate_IDs():
    sys, HXN, feed = build_system(N_columns=2)
    sys.simulate()
    IDs = [hx.ID for hx in HXN.original_heat_exchangers]
    assert len(IDs) != len(set(IDs)), 'test needs duplicate auxiliary IDs'
    fresh = network_results(HXN)
    simulate_cached(sys, HXN)
    assert_same_results(network_results(HXN), fresh)

def test_energy_balance_error_contributions_ignored_none():
    sys, HXN, feed = build_system()
    sys.simulate()
    N = len(HXN.original_heat_utils)
    errors = HXN._energy_balance_error_contributions()
    assert len(errors) == N
    assert HXN.ignored is None

# ---------------------------------------------------------------------------
# Problem-table (pinch) analysis
# ---------------------------------------------------------------------------

from biosteam.facilities.hxn.hxn_synthesis import (
    temperature_interval_pinch_analysis, problem_table,
)

def utility_hx(ID, T, P, phase, T_out, **flow):
    """A simulated HXutility acting as one process stream (kmol/hr flows)."""
    s = bst.Stream(ID + '_in', T=T, P=P, phase=phase, units='kmol/hr', **flow)
    hx = bst.HXutility(ID, ins=s, T=T_out,
                       rigorous=(phase == 'g'))  # liquid streams stay liquid (report's case)
    hx.simulate()
    return hx

def synthetic_units():
    """Report's 4-stream case: two cold, one condensing hot, one hot liquid."""
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    bst.main_flowsheet.set_flowsheet('test_hxn_synthetic')
    return [utility_hx('C1', 300., 101325., 'l', 390., Water=2000.),
            utility_hx('C2', 310., 101325., 'l', 345., Water=500., Ethanol=500.),
            utility_hx('H1', 352., 101325., 'g', 340., Ethanol=300.),
            utility_hx('H2', 420., 5e5, 'l', 320., Water=800.)]

def heat_utilities(units):
    hus = [hx.heat_utilities[0] for hx in units]
    hus.sort(key=lambda hu: hu.duty)
    return hus

def pinch_streams(hus):
    """Inlet/quenched-outlet stream copies exactly as the synthesizer prepares them."""
    streams_inlet = [hu.unit.ins[0].copy() for hu in hus]
    streams_quenched = [hu.unit.outs[0].copy() for hu in hus]
    for s in streams_quenched: s.vle(H=s.H, P=s.P)
    is_hot = [hu.duty < 0 for hu in hus]
    return streams_inlet, streams_quenched, is_hot

def assert_energy_consistent(hus, T_min_app):
    """Invariants of a correct problem table, independent of the synthesizer."""
    unit_duties = np.array([hu.unit_duty for hu in hus])
    table = problem_table(*pinch_streams(hus), T_min_app)
    # (a) each stream's grid contributions telescope to its real duty
    #     (hot: +|dH|, cold: -dH)
    per_stream = table.interval_H.sum(axis=1) + table.point_H.sum(axis=1)
    assert_allclose(per_stream, -unit_duties, rtol=1e-9)
    # (b) targets are non-negative and hot - cold is the net demand
    assert table.hot_util_load >= 0. and table.cold_util_load >= 0.
    assert_allclose(table.hot_util_load - table.cold_util_load,
                    unit_duties.sum(), rtol=1e-9)
    # (c) a target can never exceed the un-integrated load
    assert table.hot_util_load <= unit_duties[unit_duties > 0].sum() * (1 + 1e-9)
    assert table.cold_util_load <= -unit_duties[unit_duties < 0].sum() * (1 + 1e-9)
    # (d) the public wrapper reports the same targets
    pinch_T_arr, hot, cold, *_ = temperature_interval_pinch_analysis(hus, T_min_app)
    assert_allclose([hot, cold], [table.hot_util_load, table.cold_util_load], rtol=1e-12)
    return table

def test_problem_table_energy_consistency_doctest_system():
    sys, HXN, feed = build_system()
    sys.simulate()
    hus = HXN._get_original_heat_utilties()
    hus.sort(key=lambda hu: hu.duty)
    assert_energy_consistent(hus, 5.)

@pytest.mark.parametrize('T_min_app', [5., 10., 20.])
def test_problem_table_energy_consistency_synthetic(T_min_app):
    hus = heat_utilities(synthetic_units())
    table = assert_energy_consistent(hus, T_min_app)
    # the condensing ethanol stream (352 K in, ~351.4 K dew point) must keep
    # its latent heat: hot target well below the un-integrated heating load
    heating = sum(hu.unit_duty for hu in hus if hu.unit_duty > 0)
    assert table.hot_util_load < 0.5 * heating

def test_problem_table_two_streams_closed_form():
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    bst.main_flowsheet.set_flowsheet('test_hxn_two_streams')
    hot = utility_hx('Hw', 400., 5e5, 'l', 300., Water=1000.)
    # Case 1 (threshold problem): the hot stream covers every interval of the
    # smaller cold stream -> no hot utility, cold utility = net surplus.
    cold = utility_hx('Cs', 300., 5e5, 'l', 390., Water=900.)
    hus = heat_utilities([hot, cold])
    table = problem_table(*pinch_streams(hus), 5.)
    Q_hot = -hot.heat_utilities[0].unit_duty
    Q_cold = cold.heat_utilities[0].unit_duty
    assert table.hot_util_load == 0.
    assert_allclose(table.cold_util_load, Q_hot - Q_cold, rtol=1e-9)
    assert table.pinch_T == table.Ts[0]   # no pinch: everything sits below it
    # Case 2: the larger cold stream is short of heat everywhere; the cascade
    # minimum is at the cold inlet (300 K on the shifted scale) and the hot
    # utility is the cold duty minus what the hot stream gives down to 305 K.
    cold = utility_hx('Cb', 300., 5e5, 'l', 390., Water=1200.)
    hus = heat_utilities([hot, cold])
    table = problem_table(*pinch_streams(hus), 5.)
    s = hot.ins[0].copy(); s.vle(T=305., P=s.P)
    Q_hot_above_pinch = hot.ins[0].H - s.H
    Q_cold = cold.heat_utilities[0].unit_duty
    assert_allclose(table.hot_util_load, Q_cold - Q_hot_above_pinch, rtol=1e-9)
    assert table.pinch_T == 300.
    # remaining hot-stream heat below the pinch leaves as cold utility
    assert_allclose(table.cold_util_load, s.H - hot.outs[0].H, rtol=1e-9)

def test_problem_table_non_monotone_stream_is_point_load():
    # A heated stream whose outlet is colder than its inlet (e.g. a column
    # reboiler outlet at VLE): treated as an isothermal load at T_out.
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    a = bst.Stream('a', Water=100., T=372., P=101325., phase='l', units='kmol/hr')
    b = bst.Stream('b', Water=100., T=371., P=101325., phase='g', units='kmol/hr')
    dH = b.H - a.H
    assert dH > 0
    table = problem_table([a], [b], [False], 5.)
    assert_allclose(table.point_H, [[-dH]])
    assert table.interval_H.shape == (1, 0)
    assert_allclose([table.hot_util_load, table.cold_util_load, table.pinch_T],
                    [dH, 0., 371.])
    table = problem_table([b], [a], [True], 5.)
    assert_allclose([table.hot_util_load, table.cold_util_load, table.pinch_T],
                    [0., dH, 372. - 5.])  # outlet T 372 K, shifted by T_min_app

def test_problem_table_point_load_cannot_heat_above_itself():
    # A source at shifted temperature T cannot serve sinks above T: an
    # isothermal condensing hot stream at 400 K (shifted to 395 K) must not
    # cover the 395-398 K segment of a cold stream that runs 392 -> 398 K.
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    hot_in = bst.Stream('hv', Water=100., T=400., P=101325., phase='g',
                        units='kmol/hr')
    hot_out = bst.Stream('hl', Water=100., T=400., P=101325., phase='l',
                         units='kmol/hr')
    P_hot = hot_in.H - hot_out.H
    assert P_hot > 0
    cold_in = bst.Stream('cl', Water=1000., T=392., P=5e5, phase='l',
                         units='kmol/hr')
    cold_out = cold_in.copy('cl_out'); cold_out.vle(T=398., P=5e5)
    s = cold_in.copy(); s.vle(T=395., P=5e5)
    H_392, H_395, H_398 = cold_in.H, s.H, cold_out.H
    assert P_hot > H_398 - H_392   # more than enough heat overall
    table = problem_table([hot_in, cold_in], [hot_out, cold_out],
                          [True, False], 5.)
    assert_allclose(table.Ts, [398., 395., 392.])
    # heat arriving at 395 K, before the point load there, is short by the
    # 395-398 K segment of the cold stream
    assert_allclose(table.hot_util_load, H_398 - H_395, rtol=1e-9)
    assert table.pinch_T == 395.
    assert_allclose(table.cold_util_load, P_hot - (H_395 - H_392), rtol=1e-9)

def test_synthetic_network_reaches_MER():
    units = synthetic_units()
    HXN = bst.HeatExchangerNetwork('HXN', T_min_app=5.)
    sys = bst.System.from_units('sys_synthetic', units=[*units, HXN])
    sys.simulate()
    hus = heat_utilities(units)
    table = problem_table(*pinch_streams(hus), 5.)
    actual_heat = sum(hu.unit_duty for hx in HXN.new_HX_utils
                      for hu in hx.heat_utilities if hu.unit_duty > 0)
    # the greedy heuristic reaches the (corrected) MER on this case; it can
    # never legitimately beat it
    assert_allclose(actual_heat, table.hot_util_load, rtol=1e-2)
    # the lower bound is exact here only because every synthetic inlet is an
    # equilibrium state (so the clipped table is exact) and the synthesizer
    # respects T_min_app; with non-equilibrium inlets the table is conservative
    assert actual_heat >= table.hot_util_load * (1 - 1e-3)

# --- pinch diagram -----------------------------------------------------------

class _FakeStage:
    def __init__(self, unit): self.unit = unit

class _FakeLifeCycle:
    def __init__(self, units, cold=True):
        self.cold = cold
        self.life_cycle = [_FakeStage(u) for u in units]

def test_pinch_diagram_column_order_follows_stream_direction():
    from biosteam.facilities.hxn.hxn_synthesis import _order_exchanger_columns
    h1, h2, h3 = 'h1', 'h2', 'h3'
    # stream A visits h2 then h1; stream B visits h1 then h3 -> h2, h1, h3
    cycles = [_FakeLifeCycle([h2, h1]), _FakeLifeCycle([h1, h3])]
    assert _order_exchanger_columns([h1, h2, h3], cycles) == [h2, h1, h3]
    # exchangers not in the requested subset are ignored, order is stable
    assert _order_exchanger_columns([h3, h1], cycles) == [h1, h3]
    # a hot stream flows right to left, so its stage order is reversed:
    # hot stream visits h1 then h3 -> h3 left of h1
    cycles = [_FakeLifeCycle([h2, h1]), _FakeLifeCycle([h1, h3], cold=False)]
    assert _order_exchanger_columns([h1, h2, h3], cycles) == [h2, h3, h1]
    # contradictory constraints (a cycle) fall back to the given order
    cycles = [_FakeLifeCycle([h1, h2]), _FakeLifeCycle([h2, h1])]
    assert _order_exchanger_columns([h2, h1], cycles) == [h2, h1]

def _gid_artists(ax, prefix):
    return [a for a in ax.findobj() if (a.get_gid() or '').startswith(prefix)]

def test_pinch_diagram_doctest_system():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    sys, HXN, feed = build_system()
    sys.simulate()
    assert HXN.new_HXs_hot_side + HXN.new_HXs_cold_side == HXN.new_HXs
    fig, ax = HXN.plot_pinch_diagram()
    try:
        # one connector per process exchanger
        connectors = _gid_artists(ax, 'HX:')
        assert {a.get_gid() for a in connectors} == {'HX:' + hx.ID for hx in HXN.new_HXs}
        # one utility marker per utility exchanger with a duty above Qmin
        utils = _gid_artists(ax, 'Util:')
        expected = {'Util:' + hx.ID for hx in HXN.new_HX_utils
                    if abs(hx.outs[0].H - hx.ins[0].H) > HXN.Qmin}
        assert {a.get_gid() for a in utils} == expected
        # one row per stream, with inlet temperatures in degC
        texts = {t.get_text() for t in ax.texts}
        for T in HXN.inlet_Ts: assert f'{T - 273.15:.1f}' in texts
        for T in HXN.outlet_Ts: assert f'{T - 273.15:.1f}' in texts
        for i in range(len(HXN.inlet_Ts)): assert str(i) in texts
    finally:
        plt.close(fig)

if __name__ == '__main__':
    test_cache_network_matches_fresh_synthesis()
    test_cache_network_perturbed_feed()
    test_cache_network_duplicate_IDs()
    test_energy_balance_error_contributions_ignored_none()
    test_problem_table_energy_consistency_doctest_system()
    for T_min_app in (5., 10., 20.):
        test_problem_table_energy_consistency_synthetic(T_min_app)
    test_problem_table_two_streams_closed_form()
    test_problem_table_non_monotone_stream_is_point_load()
    test_problem_table_point_load_cannot_heat_above_itself()
    test_synthetic_network_reaches_MER()
    test_pinch_diagram_column_order_follows_stream_direction()
    test_pinch_diagram_doctest_system()
