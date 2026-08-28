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
    temperature_interval_pinch_analysis, problem_table, load_duties, pinch_state,
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

def test_load_duties_non_equilibrium_inlet_conserves_energy():
    # A non-rigorous HXutility can carry a superheated liquid (ethanol at
    # 370 K, 1 atm; bp 351.4 K). Flashing it at the 355 K pinch gives vapor
    # with far more enthalpy than the inlet has; the pinch split must clip
    # to the stream's real enthalpy range so that the hot-side and cold-side
    # loads sum to the stream's duty (and are not phantom latent heat).
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    s_in = bst.Stream('ne_in', Ethanol=700., T=370., P=101325., phase='l',
                      units='kmol/hr')
    s_out = bst.Stream('ne_out', Ethanol=700., T=310., P=101325., phase='l',
                       units='kmol/hr')
    duty = s_in.H - s_out.H
    flashed = s_in.copy(); flashed.vle(T=355., P=101325.)
    assert flashed.H > s_in.H   # the trap: equilibrium enthalpy above H_in
    Q_hot_side, Q_cold_side = {}, {}
    load_duties([s_in], [s_out], np.array([355.]), np.array([310.]), [0],
                lambda i: False, Q_hot_side, Q_cold_side)
    assert Q_hot_side[0][0] == Q_cold_side[0][0] == 'cool'
    assert_allclose(Q_hot_side[0][1] + Q_cold_side[0][1], duty, rtol=1e-9)
    # the stream cannot give any heat above the pinch: its whole enthalpy
    # range lies below the equilibrium state at 355 K
    assert_allclose(Q_hot_side[0][1], 0., atol=1e-6)
    assert_allclose(Q_cold_side[0][1], duty, rtol=1e-9)
    # and the transient state used for matching carries exactly H_in, at
    # equilibrium (two-phase at the bubble point), not at the fictitious
    # 370 K of the superheated liquid, so matching is consistent with the
    # problem table
    s = pinch_state(s_in, s_out, 355.)
    assert_allclose(s.H, s_in.H, rtol=1e-9)
    assert s.T < 355. and 'g' in s.phase and 'l' in s.phase

def test_pinch_state_at_endpoints_uses_real_states():
    # pinch_T == T_in must put the whole duty on one side even when the
    # inlet is not at equilibrium: flashing a superheated liquid at its own
    # T does not reproduce H_in, and the resulting enthalpy can lie inside
    # the stream's range so that clipping alone would not catch it.
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    s_in = bst.Stream('sh_in', Water=100., T=380., P=101325., phase='l',
                      units='kmol/hr')
    s_out = bst.Stream('sh_out', Water=100., T=400., P=101325., phase='g',
                       units='kmol/hr')
    duty = s_out.H - s_in.H
    flashed = s_in.copy(); flashed.vle(T=380., P=101325.)
    assert s_in.H < flashed.H < s_out.H   # the trap: in range, but not H_in
    Q_hot_side, Q_cold_side = {}, {}
    load_duties([s_in], [s_out], np.array([380.]), np.array([400.]), [0],
                lambda i: True, Q_hot_side, Q_cold_side)
    assert_allclose(Q_hot_side[0][1], duty, rtol=1e-9)
    assert_allclose(Q_cold_side[0][1], 0., atol=1e-6)
    assert_allclose(pinch_state(s_in, s_out, 380.).H, s_in.H, rtol=1e-9)
    assert_allclose(pinch_state(s_in, s_out, 400.).H, s_out.H, rtol=1e-9)
    # a phase-mislabelled non-condensable: clipping must return the real
    # end state, never an equilibrium re-flash (N2 "liquid" at 400 K would
    # otherwise come back as gas at thousands of K)
    bst.settings.set_thermo(['Water', 'N2'], cache=True)
    h_in = bst.Stream('n2_in', N2=100., T=400., P=101325., phase='l',
                      units='kmol/hr')
    h_out = bst.Stream('n2_out', N2=100., T=320., P=101325., phase='l',
                       units='kmol/hr')
    s = pinch_state(h_in, h_out, 355.)
    assert h_out.H <= s.H <= h_in.H
    assert 320. <= s.T <= 400.

def test_unordered_network_path_warns_with_context(monkeypatch):
    # thermosteam's Network.sort warns 'network path could not be determined'
    # when its ordering heuristic does not settle; HXN re-raises that as its
    # own warning so the user knows which facility it concerns.
    import thermosteam as tmo
    original = tmo.Network.from_units
    def from_units_with_warning(*args, **kwargs):
        network = original(*args, **kwargs)
        warnings.warn('network path could not be determined', RuntimeWarning)
        return network
    monkeypatch.setattr(tmo.Network, 'from_units', from_units_with_warning)
    units = synthetic_units()
    HXN = bst.HeatExchangerNetwork('HXN', T_min_app=5.)
    sys = bst.System.from_units('sys_unordered', units=[*units, HXN])
    with pytest.warns(RuntimeWarning, match='heat exchanger network path could not be fully ordered'):
        sys.simulate()
    assert abs(HXN.energy_balance_percent_error) < 1e-6

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

def test_pinch_diagram_stream_labels():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from biosteam.facilities.hxn.hxn_synthesis import _auxiliary_name, _stream_label
    sys, HXN, feed = build_system()
    sys.simulate()
    D1 = bst.main_flowsheet.unit.D0
    D1_H1 = bst.main_flowsheet.unit.D0_H1
    assert _auxiliary_name(D1.condenser) == 'condenser'
    assert _auxiliary_name(D1_H1) is None
    # label composition
    assert _stream_label(D1.condenser, True, True, False) == 'D0 - condenser'
    assert _stream_label(D1.condenser, True, False, False) == 'D0'
    assert _stream_label(D1.condenser, False, True, False) == 'condenser'
    assert _stream_label(D1_H1, True, True, True) == 'D0_H1 (' + D1_H1.ins[0].ID + ')'
    assert _stream_label(D1_H1, False, False, True) == D1_H1.ins[0].ID
    assert _stream_label(D1_H1, False, False, False) == ''
    # an unnamed inlet adds nothing
    assert _stream_label(D1.reboiler, True, True, True) == 'D0 - reboiler'
    # every stream gets a label on the figure; toggles remove them
    fig, ax = HXN.plot_pinch_diagram()
    try:
        labels = {a.get_gid(): a.get_text() for a in _gid_artists(ax, 'Label:')}
        hxs = HXN.original_heat_exchangers
        assert labels == {
            f'Label:{i}': _stream_label(hx, True, True, True)
            for i, hx in enumerate(hxs)
        }
        assert any(text.startswith('D0 - condenser') for text in labels.values())
        assert any(text == 'F1 - heat_exchanger (feed_flash)' for text in labels.values())
    finally:
        plt.close(fig)
    fig, ax = HXN.plot_pinch_diagram(show_units=False, show_auxiliary_units=False,
                                     show_stream_IDs=False)
    try:
        assert not _gid_artists(ax, 'Label:')
    finally:
        plt.close(fig)

def test_pinch_diagram_requires_simulation():
    sys, HXN, feed = build_system()
    with pytest.raises(RuntimeError, match='simulate'):
        HXN.plot_pinch_diagram()

def test_pinch_diagram_legend():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    sys, HXN, feed = build_system()
    sys.simulate()
    fig, ax = HXN.plot_pinch_diagram()
    try:
        legend = ax.get_legend()
        assert legend is not None
        labels = [t.get_text() for t in legend.get_texts()]
        assert labels == ['Cold stream', 'Hot stream', 'Process heat exchange',
                          'Hot utility', 'Cold utility', 'Pinch']
        # utility markers are colored by utility type, consistent with the legend
        handles = dict(zip(labels, legend.legend_handles))
        hot_util_color = handles['Hot utility'].get_markeredgecolor()
        cold_util_color = handles['Cold utility'].get_markeredgecolor()
        assert hot_util_color != cold_util_color
        heaters = {hx.ID for hx in HXN.new_HX_utils if hx.outs[0].H > hx.ins[0].H}
        for artist in _gid_artists(ax, 'Util:'):
            heater = artist.get_gid()[len('Util:'):] in heaters
            expected = hot_util_color if heater else cold_util_color
            assert artist.get_markeredgecolor() == expected, artist.get_gid()
    finally:
        plt.close(fig)
    fig, ax = HXN.plot_pinch_diagram(show_legend=False)
    try:
        assert ax.get_legend() is None
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
    test_pinch_diagram_stream_labels()
    test_pinch_diagram_legend()
    test_pinch_diagram_requires_simulation()
