# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-, Yoel Cortes-Pena <yoelcortes@gmail.com>
# Copyright (C) 2026-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Regression tests for heat exchanger network synthesis on synthetic systems.

Ten synthetic systems of increasing complexity (all with phase-changing
streams from case 3 on). For each, the synthesized network must

(i)   close its energy balance (|error| < 0.1 %) without RuntimeWarnings,
(ii)  never beat the minimum-energy-requirement (MER) targets of the problem
      table computed on the same streams, and
(iii) recover at least as much heat as documented in ``CASES`` below, so that
      no future change to ``hxn`` silently makes the synthesizer perform worse.

The documented utility loads were recorded by running this file directly
(``python tests/test_hxn_regression.py`` prints them): cases 1-9 at commit
``1ab689ff`` (branch ``hxn-pinch-diagram``), case 10 after the synthesizer
fixes on ``hxn-regression-tests`` (non-equilibrium inlets clipped to the
stream's enthalpy range; network path ordered by its connections; H_lim
honored at the bubble point). Improvements leave slack; a
maintainer lowers the numbers deliberately when a better network is
intended. Never raise them to make a failing test pass.
"""
import warnings
import pytest
import biosteam as bst
from numpy.testing import assert_allclose
from biosteam.facilities.hxn.hxn_synthesis import problem_table

EB_TOLERANCE = 0.1   # percent
MER_RTOL = 1e-3      # network may not beat the MER target by more than this
DOC_RTOL = 1e-3      # network may not be worse than documented by more than this

def utility_hx(ID, T, P, phase, T_out, rigorous=None, **flow):
    """A simulated HXutility acting as one process stream (kmol/hr flows)."""
    s = bst.Stream(ID + '_in', T=T, P=P, phase=phase, units='kmol/hr', **flow)
    if rigorous is None: rigorous = phase == 'g'
    hx = bst.HXutility(ID, ins=s, T=T_out, rigorous=rigorous)
    hx.simulate()
    return hx

def boiling_hx(ID, T, P, T_out, **flow):
    """A cold liquid stream heated past its bubble point (rigorous VLE)."""
    return utility_hx(ID, T, P, 'l', T_out, rigorous=True, **flow)

def setup(name, chemicals=('Water', 'Ethanol')):
    bst.settings.set_thermo(list(chemicals), cache=True)
    bst.main_flowsheet.set_flowsheet('test_hxn_regression_' + name)

# ---------------------------------------------------------------------------
# Cases
# ---------------------------------------------------------------------------

def case_01_two_liquids():
    """Hot liquid, cold liquid; trivial counter-current match."""
    setup('01')
    return [utility_hx('H1', 400., 5e5, 'l', 320., Water=1000.),
            utility_hx('C1', 300., 101325., 'l', 360., Water=1000.)], 5.

def case_02_pinch_limited():
    """Cold target above the hot inlet: part of the heating must be utility."""
    setup('02')
    return [utility_hx('H1', 360., 5e5, 'l', 320., Water=1000.),
            utility_hx('C1', 300., 5e5, 'l', 380., Water=800.)], 5.

def case_03_condenser_two_colds():
    """Condensing ethanol vapor against two cold liquids."""
    setup('03')
    return [utility_hx('H1', 355., 101325., 'g', 340., Ethanol=400.),
            utility_hx('C1', 300., 101325., 'l', 345., Water=1500.),
            utility_hx('C2', 310., 101325., 'l', 340., Water=300., Ethanol=300.)], 5.

def case_04_report_case():
    """The 4-stream report case: 2 colds, condensing hot, hot liquid."""
    setup('04')
    return [utility_hx('C1', 300., 101325., 'l', 390., Water=2000.),
            utility_hx('C2', 310., 101325., 'l', 345., Water=500., Ethanol=500.),
            utility_hx('H1', 352., 101325., 'g', 340., Ethanol=300.),
            utility_hx('H2', 420., 5e5, 'l', 320., Water=800.)], 5.

def case_05_boiling_cold():
    """A cold stream that boils (water -> steam) and a condensing hot stream."""
    setup('05')
    return [boiling_hx('C1', 330., 101325., 380., Water=300.),
            utility_hx('C2', 300., 101325., 'l', 350., Water=1000.),
            utility_hx('H1', 420., 5e5, 'g', 330., Water=250.),
            utility_hx('H2', 400., 5e5, 'l', 310., Water=1500.)], 5.

def case_06_mixed_pressures():
    """5-bar condensing hot against 1-atm boiling cold; T_min_app = 10."""
    setup('06')
    return [utility_hx('H1', 430., 5e5, 'g', 400., Water=200.),
            utility_hx('H2', 380., 5e5, 'l', 320., Water=1200.),
            boiling_hx('C1', 340., 101325., 375., Water=150.),
            utility_hx('C2', 300., 101325., 'l', 360., Ethanol=800.),
            utility_hx('C3', 320., 101325., 'l', 390., Water=600.)], 10.

def case_07_threshold():
    """Threshold problem: heating dominates so the cold target is ~0;
    one hot stream condenses and subcools."""
    setup('07')
    return [utility_hx('H1', 375., 101325., 'g', 310., Water=80.),
            utility_hx('H2', 360., 101325., 'l', 330., Water=200.),
            utility_hx('C1', 300., 101325., 'l', 370., Water=1500.),
            utility_hx('C2', 305., 101325., 'l', 350., Ethanol=800.),
            boiling_hx('C3', 340., 101325., 355., Ethanol=200.),
            utility_hx('C4', 320., 101325., 'l', 360., Water=700.)], 5.

def case_08_two_condensers():
    """Two condensers at different temperatures (ethanol 1 atm, water 2 bar)
    against three colds, one of which boils."""
    setup('08')
    return [utility_hx('H1', 355., 101325., 'g', 335., Ethanol=300.),
            utility_hx('H2', 400., 2e5, 'g', 360., Water=150.),
            utility_hx('H3', 390., 5e5, 'l', 330., Water=900.),
            utility_hx('C1', 300., 101325., 'l', 345., Water=1200.),
            boiling_hx('C2', 330., 101325., 370., Ethanol=250.),
            utility_hx('C3', 310., 101325., 'l', 380., Water=700.)], 5.

def case_09_near_degenerate():
    """Eight streams with two near-equal pinch candidates, a partially
    condensing hot stream (wet outlet) and a partially boiling cold stream."""
    setup('09')
    return [utility_hx('H1', 380., 101325., 'g', 372., Water=120.),   # partial condensation
            utility_hx('H2', 365., 101325., 'g', 330., Ethanol=250.),
            utility_hx('H3', 410., 5e5, 'l', 340., Water=700.),
            utility_hx('H4', 345., 101325., 'l', 305., Ethanol=900.),
            boiling_hx('C1', 350., 101325., 373.5, Water=200.),        # partial boiling
            utility_hx('C2', 300., 101325., 'l', 340., Water=1500.),
            boiling_hx('C3', 320., 101325., 352., Ethanol=300.),
            utility_hx('C4', 335., 101325., 'l', 395., Water=500.)], 5.

def case_10_ten_streams():
    """Ten streams mixing liquids, condensers, boilers, and pressures."""
    setup('10')
    return [utility_hx('H1', 355., 101325., 'g', 320., Ethanol=300.),
            utility_hx('H2', 420., 5e5, 'g', 340., Water=150.),
            utility_hx('H3', 395., 5e5, 'l', 330., Water=1000.),
            utility_hx('H4', 370., 101325., 'l', 310., Ethanol=700.),
            utility_hx('H5', 380., 101325., 'g', 372.5, Water=100.),   # partial condensation
            utility_hx('C1', 300., 101325., 'l', 360., Water=2000.),
            boiling_hx('C2', 330., 101325., 380., Water=200.),
            boiling_hx('C3', 320., 101325., 352., Ethanol=400.),
            utility_hx('C4', 310., 101325., 'l', 345., Water=400., Ethanol=400.),
            utility_hx('C5', 340., 101325., 'l', 390., Water=600.)], 5.

# name -> (builder, documented hot utility load [kJ/hr], documented cold utility load [kJ/hr])
# Documented values: see module docstring for provenance.
CASES = {
    'case_01_two_liquids':         (case_01_two_liquids,         0, 1.53912e+06),
    'case_02_pinch_limited':       (case_02_pinch_limited,       1.81522e+06, 0),
    'case_03_condenser_two_colds': (case_03_condenser_two_colds, 0, 9.49905e+06),
    'case_04_report_case':         (case_04_report_case,         2.37319e+06, 3.56871e+06),
    'case_05_boiling_cold':        (case_05_boiling_cold,        5.8433e+06, 1.04356e+07),
    'case_06_mixed_pressures':     (case_06_mixed_pressures,     7.05541e+06, 4.93079e+06),
    'case_07_threshold':           (case_07_threshold,           1.85977e+07, 0),
    'case_08_two_condensers':      (case_08_two_condensers,      3.02237e+06, 7.36431e+06),
    'case_09_near_degenerate':     (case_09_near_degenerate,     1.40965e+07, 9.66427e+06),
    'case_10_ten_streams':         (case_10_ten_streams,         1.40742e+07, 8.06488e+06),
}

# ---------------------------------------------------------------------------
# Harness
# ---------------------------------------------------------------------------

def synthesize(builder):
    units, T_min_app = builder()
    HXN = bst.HeatExchangerNetwork('HXN', T_min_app=T_min_app)
    sys = bst.System.from_units('sys', units=[*units, HXN])
    with warnings.catch_warnings():
        warnings.simplefilter('error', RuntimeWarning)
        # thermosteam registry bookkeeping on temporary stream copies; not numerical
        warnings.filterwarnings('ignore', message='.*has been replaced in registry',
                                category=RuntimeWarning)
        sys.simulate()
    return units, HXN, T_min_app

def mer_targets(units, T_min_app):
    hus = [hx.heat_utilities[0] for hx in units]
    hus.sort(key=lambda hu: hu.duty)
    streams_inlet = [hu.unit.ins[0].copy() for hu in hus]
    streams_quenched = [hu.unit.outs[0].copy() for hu in hus]
    for s in streams_quenched: s.vle(H=s.H, P=s.P)
    is_hot = [hu.duty < 0 for hu in hus]
    table = problem_table(streams_inlet, streams_quenched, is_hot, T_min_app)
    return table.hot_util_load, table.cold_util_load

def actual_loads(HXN):
    hus = [hu for hx in HXN.new_HX_utils for hu in hx.heat_utilities]
    heat = sum(hu.unit_duty for hu in hus if hu.unit_duty > 0)
    cool = -sum(hu.unit_duty for hu in hus if hu.unit_duty < 0)
    return heat, cool

@pytest.mark.parametrize('name', list(CASES))
def test_hxn_regression(name):
    builder, doc_heat, doc_cool = CASES[name]
    units, HXN, T_min_app = synthesize(builder)
    # (i) energy balance
    assert abs(HXN.energy_balance_percent_error) < EB_TOLERANCE, name
    # (ii) MER targets are a lower bound; energy identity holds
    heat, cool = actual_loads(HXN)
    hot_target, cold_target = mer_targets(units, T_min_app)
    net_duty = sum(hx.heat_utilities[0].unit_duty for hx in units)
    assert heat >= hot_target * (1 - MER_RTOL), (name, heat, hot_target)
    assert cool >= cold_target * (1 - MER_RTOL), (name, cool, cold_target)
    assert_allclose(heat - cool, net_duty, rtol=1e-3, err_msg=name)
    # (iii) never worse than documented
    assert doc_heat is not None and doc_cool is not None, f'{name}: baseline not recorded'
    assert heat <= doc_heat * (1 + DOC_RTOL) + 1e-9, (name, heat, doc_heat)
    assert cool <= doc_cool * (1 + DOC_RTOL) + 1e-9, (name, cool, doc_cool)

if __name__ == '__main__':
    for name, (builder, *_) in CASES.items():
        units, HXN, T_min_app = synthesize(builder)
        heat, cool = actual_loads(HXN)
        hot_target, cold_target = mer_targets(units, T_min_app)
        print(f"{name}: heat={heat:.6g} cool={cool:.6g} "
              f"(MER hot={hot_target:.6g} cold={cold_target:.6g}; "
              f"EB error={HXN.energy_balance_percent_error:.4f}%)")
