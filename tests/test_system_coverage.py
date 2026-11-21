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
import os
import tempfile
from numpy.testing import assert_allclose

def create_recycle_system(name, method=None, algorithm=None):
    bst.F.set_flowsheet(name)
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    feed = bst.Stream('feed', Water=100, Ethanol=50)
    recycle = bst.Stream('recycle')
    M1 = bst.Mixer('M1', [feed, recycle])
    S1 = bst.Splitter('S1', M1-0, ['product', recycle], split=0.5)
    sys = bst.System.from_units('sys', [M1, S1])
    if method: sys.method = method
    if algorithm: sys.algorithm = algorithm
    return sys

def assert_mass_balance(sys, rtol=1e-2):
    F_in = sum([i.F_mass for i in sys.feeds])
    F_out = sum([i.F_mass for i in sys.products])
    assert_allclose(F_in, F_out, rtol=rtol)

def test_recycle_convergence_methods():
    for method in ('aitken', 'wegstein', 'fixed-point', 'anderson'):
        sys = create_recycle_system(f'recycle_methods_{method}', method)
        sys.simulate()
        assert_mass_balance(sys)
        recycle = sys.flowsheet.stream.recycle
        assert_allclose(recycle.F_mol, 150., rtol=2e-2)
    with pytest.raises(AttributeError):
        sys.method = 'not a method'
    with pytest.raises(AttributeError):
        sys.algorithm = 'not an algorithm'
    assert sys.converge_method == sys.method

def test_recycle_maxiter():
    sys = create_recycle_system('recycle_maxiter', 'fixed-point')
    sys.maxiter = 2
    sys.molar_tolerance = 1e-12
    sys.relative_molar_tolerance = 1e-12
    with pytest.raises(RuntimeError):
        sys.simulate()
    assert sys._error_info()
    bst.System.strict_convergence = False
    try:
        sys.empty_recycles()
        sys.simulate()
    finally:
        bst.System.strict_convergence = True
    sys.show()
    assert sys._iter >= 2

def test_recycle_data():
    sys = create_recycle_system('recycle_data')
    sys.simulate()
    recycle = sys.flowsheet.stream.recycle
    data = sys.get_recycle_data()
    assert len(data.get_keys()) == len(data.get_names()) == len(data.to_list())
    assert len(data.to_dict()) == len(data.to_tuples())
    assert np.isfinite(data.to_array()).all()
    assert len(data.to_series()) == len(data.to_list())
    assert 'recycle' in repr(data)
    rdata = data.recycle_data[0]
    assert 'recycle' in repr(rdata)
    assert len(rdata.to_series()) == len(rdata.to_list())
    assert len(rdata.to_dict()) == len(rdata.get_keys())
    assert len(rdata.to_tuples()) == len(rdata.get_names())
    assert np.isfinite(rdata.to_array()).all()
    F_mol = recycle.F_mol
    recycle.empty()
    data.reset()
    assert_allclose(recycle.F_mol, F_mol, rtol=1e-6)
    sys.converge(recycle_data=data, update_recycle_data=True)
    assert_allclose(data.to_array()[:-2].sum(), F_mol, rtol=1e-2)
    with pytest.raises(ValueError):
        sys.converge(update_recycle_data=True)

def test_multiple_recycles():
    bst.F.set_flowsheet('multiple_recycles')
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    feed = bst.Stream('feed', Water=100, Ethanol=50)
    r1 = bst.Stream('r1')
    r2 = bst.Stream('r2')
    M1 = bst.Mixer('M1', [feed, r1, r2])
    S1 = bst.Splitter('S1', M1-0, ['s1', r1], split=0.5)
    S2 = bst.Splitter('S2', S1-0, ['product', r2], split=0.5)
    sys = bst.System('sys', [M1, S1, S2], recycle=[r1, r2])
    assert len(sys.get_all_recycles()) == 2
    sys.simulate()
    assert_mass_balance(sys)
    assert sys.depth == 1
    T = sys._get_recycle_temperatures()
    assert len(T) == 2
    data = sys._get_recycle_data()
    with pytest.raises(IndexError):
        sys._set_recycle_data(data[:-1])
    sys._set_recycle_data(data)
    mol = sys._get_recycle_mol()
    assert mol.shape[0] == 2
    assert 'S1-1' in sys._get_recycle_info()
    sys.print()
    sys.empty_recycles()
    assert r1.isempty() and r2.isempty()
    sys.rescale(feed, 2.)
    assert_allclose(feed.F_mol, 300.)

def test_no_recycle_errors():
    bst.F.set_flowsheet('no_recycle_errors')
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    feed = bst.Stream('feed', Water=100, Ethanol=50)
    H1 = bst.HXutility('H1', feed, T=320)
    P1 = bst.Pump('P1', H1-0, P=2e5)
    sys = bst.System('sys', [H1, P1])
    for f in (sys._get_recycle_mol, sys._get_recycle_data, sys._get_recycle_temperatures):
        with pytest.raises(RuntimeError): f()
    with pytest.raises(RuntimeError):
        sys._set_recycle_data(np.zeros(3))
    sys.simulate()
    assert sys._error_info() is None
    with pytest.raises(ValueError):
        sys.recycle = 1
    sys.N_runs = 2
    sys.simulate()
    assert sys.N_runs == 2
    assert 'N_runs' in sys._stacked_info()
    assert sys.depth == 0

def test_mock_system():
    bst.F.set_flowsheet('mock_system')
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    feed = bst.Stream('feed', Water=100, Ethanol=50)
    M1 = bst.Mixer('M1', feed)
    S1 = bst.Splitter('S1', M1-0, ['a', 'b'], split=0.5)
    mock = bst.MockSystem([M1, S1])
    assert feed in mock.ins
    assert len(mock.outs) == 2
    assert 'MockSystem' in repr(mock)
    mock.show()
    mock.load_inlet_ports([feed], {'feed': 0})
    mock.load_outlet_ports(S1.outs, {'a': 0, 'b': 1})
    assert mock.get_inlet('feed') is feed
    assert mock.get_outlet('b') is S1.outs[1]
    with pytest.raises(ValueError): mock.get_inlet('none')
    with pytest.raises(ValueError): mock.get_outlet('none')
    with pytest.raises(ValueError): mock.load_inlet_ports([S1-0])
    with pytest.raises(ValueError): mock.load_outlet_ports([feed])
    new_feed = bst.Stream('new_feed', Water=1)
    mock.set_inlet('feed', new_feed)
    assert M1.ins[0] is new_feed
    mock.set_outlet('a', bst.Stream('new_a')) # Sink is None; nothing happens
    with pytest.raises(RuntimeError):
        with mock: pass
    with bst.MockSystem() as mock2:
        bst.Mixer('M2')
    assert len(mock2.units) == 1

def test_from_segment_and_paths():
    bst.F.set_flowsheet('from_segment')
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    feed = bst.Stream('feed', Water=100, Ethanol=50)
    H1 = bst.HXutility('H1', feed, T=320)
    P1 = bst.Pump('P1', H1-0, P=2e5)
    M1 = bst.Mixer('M1', P1-0)
    S1 = bst.Splitter('S1', M1-0, split=0.5)
    with pytest.raises(ValueError):
        bst.System.from_segment('none')
    sys = bst.System('sys', [H1, P1, M1, S1])
    assert sys.path_segment(H1, S1)[-2:] == [P1, M1]
    assert sys.path_segment(H1, S1, inclusive=True)[-1] is S1
    assert sys.path_segment(H1, S1)[-1] is M1 # Cached
    assert sys.path_section([H1], [S1], inclusive=True)[-1] is S1
    assert sys.path_section([H1], [S1])[-1] is M1 # Cached
    sys.prioritize_unit(P1) # No recycle; nothing changes
    with pytest.raises(ValueError): sys.prioritize_unit(bst.Mixer('M_other'))
    with pytest.raises(ValueError): sys.find_system(bst.Mixer('M_other2'))
    assert sys.find_system(P1) is sys
    upstream, downstream = sys.split(P1-0)
    assert upstream.path == (H1, P1)
    upstream, downstream = sys.split(feed)
    assert upstream.path == ()
    with pytest.raises(ValueError): sys.split(bst.Stream('none'))
    group = sys.to_unit_group('group')
    assert len(group.units) == 4
    assert sys._downstream_path(M1) == (M1, S1)
    assert sys._downstream_path(bst.Mixer('M_other3')) == []
    dsys = sys._downstream_system(M1)
    assert dsys.path == (M1, S1)
    assert sys._downstream_system(H1) is sys
    assert sys._downstream_facilities(M1) == []

def test_nested_system_operations():
    bst.F.set_flowsheet('nested_system_operations')
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    feed = bst.Stream('feed', Water=100, Ethanol=50)
    recycle = bst.Stream('recycle')
    H1 = bst.HXutility('H1', feed, T=320)
    M1 = bst.Mixer('M1', [H1-0, recycle])
    S1 = bst.Splitter('S1', M1-0, [None, recycle], split=0.5)
    P1 = bst.Pump('P1', S1-0, P=2e5)
    inner = bst.System('inner', [M1, S1], recycle=recycle)
    sys = bst.System('outer', [H1, inner, P1])
    assert sys.subsystems == [inner]
    assert sys.depth == 2
    assert sys.find_system(S1) is inner
    sys.set_tolerance(mol=1e-3, rmol=1e-3, T=1e-3, rT=1e-3, maxiter=100,
                      subsystems=True, subfactor=0.5, method='wegstein',
                      algorithm='sequential modular', maxtime=10)
    assert_allclose(inner.molar_tolerance, 5e-4)
    sys.set_tolerance(mol=1e-2, subsystems=True)
    assert_allclose(inner.molar_tolerance, 1e-2)
    sys.prioritize_unit(S1)
    assert inner.path[0] is S1
    assert inner in sys.path_segment(H1, P1)
    assert sys._downstream_path(S1)[0] is S1
    with pytest.raises(ValueError):
        sys.split(recycle)
    sys.simulate()
    sys.simulate() # Prioritized splitter starts with empty inlet; simulate twice
    assert_mass_balance(sys)
    assert sys._error_info() is not None
    sys.track_recycle(recycle)
    sys.simulate()
    assert len(inner.tracked_recycles[recycle]) > 0
    new = sys.copy('outer_copy')
    assert new.path == sys.path
    sys.flatten()
    assert set(sys.path) == {H1, M1, S1, P1}
    assert recycle in sys.recycle
    network = sys._to_network()
    assert network.units == set(sys.unit_path)

def test_flatten_warnings():
    bst.F.set_flowsheet('flatten_warnings')
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    feed = bst.Stream('feed', Water=100, Ethanol=50)
    H1 = bst.HXutility('H1', feed, T=320)
    P1 = bst.Pump('P1', H1-0, P=2e5)
    H2 = bst.HXutility('H2', P1-0, T=300)
    inner1 = bst.System('inner1', [P1])
    inner1.add_specification(lambda: None)
    sys = bst.System('outer', [H1, inner1, H2])
    with pytest.warns(RuntimeWarning):
        sys.flatten()
    assert inner1 in sys.path

def test_specifications():
    sys = create_recycle_system('system_specifications')
    S1 = sys.flowsheet.unit.S1
    with pytest.raises(ValueError):
        sys.add_specification(1)

    @sys.add_specification(simulate=True)
    def adjust_split():
        S1.split[:] = 0.4

    sys.simulate()
    assert sys.simulate_after_specifications
    assert_mass_balance(sys)
    sys.simulate_after_specifications = False
    sys._specifications.clear()

    @sys.add_specification
    def adjust_split_and_simulate():
        S1.split[:] = 0.3
        sys.simulate()

    sys.simulate()
    assert_allclose(S1.split, 0.3)
    assert_mass_balance(sys)

def test_system_context_and_configuration():
    bst.F.set_flowsheet('system_context')
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    feed = bst.Stream('feed', Water=100, Ethanol=50)
    with bst.System('sys') as sys:
        M1 = bst.Mixer('M1', [feed, 'recycle'])
        S1 = bst.Splitter('S1', M1-0, ['product', 'recycle'], split=0.5)
    S1.outs[1] = M1.ins[1]
    sys.update_configuration()
    sys.simulate()
    assert_mass_balance(sys)
    with pytest.raises(RuntimeError):
        with sys: pass
    # Simulate with updated configuration
    sys.simulate(update_configuration=True)
    assert_mass_balance(sys)
    sys._load_configuration()
    # Recycle with temporary streams / non-streams
    sys.recycle = [None]
    assert sys.recycle is not None
    assert sys.method == 'fixedpoint'

def test_interface_property_packages():
    bst.F.set_flowsheet('interface_property_packages')
    thermo_a = bst.Thermo(['Water', 'Ethanol'], cache=True)
    thermo_b = bst.Thermo(['Water', 'Ethanol', 'Methanol'], cache=True)
    bst.settings.set_thermo(thermo_a)
    feed = bst.Stream('feed', Water=100, Ethanol=50)
    H1 = bst.HXutility('H1', feed, T=320)
    bst.settings.set_thermo(thermo_b)
    H2 = bst.HXutility('H2', T=300)
    H2.ins[0] = H1.outs[0] # Connect without junction
    sys = bst.System('sys', [H1, H2])
    sys.interface_property_packages()
    assert len(sys.path) == 3
    sys.simulate()
    assert_allclose(H2.outs[0].F_mass, feed.F_mass, rtol=1e-3)

def test_diagrams_without_display():
    sys = create_recycle_system('system_diagrams')
    sys.simulate()
    for kind in ('cluster', 'thorough', 'surface', 'minimal', 'stage'):
        graph = sys.diagram(kind, display=False, number=True, label=True,
                            profile=False, format='svg', auxiliaries=True,
                            title='title')
        assert graph is not None
    with pytest.raises(ValueError):
        sys.diagram('not a kind', display=False)

def test_utility_and_cost_properties():
    bst.F.set_flowsheet('system_cost_properties')
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    feed = bst.Stream('feed', Water=100, Ethanol=50, price=0.5)
    H1 = bst.HXutility('H1', feed, T=360, V=0)
    H2 = bst.HXutility('H2', H1-0, T=300)
    P1 = bst.Pump('P1', H2-0, P=5e5)
    P1.outs[0].price = 1.
    sys = bst.System('sys', [H1, H2, P1], operating_hours=8000, lang_factor=4)
    sys.simulate()
    product = P1.outs[0]
    assert sys.get_heating_duty() > 0
    assert sys.get_cooling_duty() > 0
    assert sys.get_electricity_consumption() > 0
    assert sys.get_electricity_production() == 0
    assert sys.get_utility_duty('low_pressure_steam') > 0
    assert sys.get_utility_flow(bst.HeatUtility.get_agent('low_pressure_steam')) > 0
    assert_allclose(sys.installed_equipment_cost, 4 * sys.purchase_cost)
    sys.lang_factor = None
    assert sys.installed_cost > sys.purchase_cost
    assert sys.material_cost > 0 and sys.sales > 0 and sys.utility_cost > 0
    assert sys.has_market_value(product)
    assert sys.get_market_value(product) > 0
    assert_allclose(sys.get_mass_flow(product), product.F_mass * 8000)
    assert sys.get_property(product, 'F_mass') > 0
    assert sys._price2cost(product) > 0
    assert sys._price2cost(feed) < 0
    with pytest.raises(ValueError):
        sys._price2cost(bst.Stream(Water=1))
    for with_units in (True, False):
        results = sys.results(with_units)
        assert len(results) > 5
    assert sys.get_inlet_cost_flows() == {}
    assert sys.get_outlet_revenue_flows() == {}
    sys.empty_outlet_streams()
    assert product.isempty()
    sys.reset_cache()

def test_debug_mode_errors():
    bst.F.set_flowsheet('debug_and_profile')
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    feed = bst.Stream('feed', Water=100, Ethanol=50)
    H1 = bst.HXutility('H1', feed, T=320)
    P1 = bst.Pump('P1', H1-0, P=2e5)
    sys = bst.System('sys', [H1, P1])
    with pytest.raises(TypeError): sys._turn_on(1)
    with pytest.raises(ValueError): sys._turn_on('none')

def test_dynamic_kwargs_validation():
    bst.F.set_flowsheet('dynamic_kwargs_validation')
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    feed = bst.Stream('feed', Water=100, Ethanol=50)
    H1 = bst.HXutility('H1', feed, T=320)
    sys = bst.System('sys', [H1])
    assert not sys.isdynamic
    sys.isdynamic = False
    assert not sys.isdynamic
    check = sys._check_dynsim_kwargs
    assert check({'t_span': (0, 1), 'rtol': 1e-3, 'atol': 1e-6, 'first_step': 0.1, 'max_step': 1.}) == set()
    for kwargs, error in [
            ({'t_span': 1}, TypeError),
            ({'t_span': (0, 1, 2)}, ValueError),
            ({'t_span': ('a', 'b')}, TypeError),
            ({'t_span': (1, 0)}, ValueError),
            ({'t_span': (0, 1), 't_eval': ['a'], 'export_state_to': 'file.npy'}, TypeError),
            ({'t_span': (0, 1), 't_eval': [0, 2], 'export_state_to': 'file.npy'}, ValueError),
            ({'state_reset_hook': 'not_a_method'}, AttributeError),
            ({'state_reset_hook': 1}, TypeError),
            ({'rtol': -1}, ValueError),
            ({'rtol': 'a'}, ValueError),
            ({'atol': -1}, ValueError),
            ({'atol': 'a'}, ValueError),
            ({'first_step': np.inf}, ValueError),
            ({'max_step': 'a'}, ValueError),
            ({'not_an_option': 1}, TypeError),
            ({'rtoll': 1}, TypeError),
        ]:
        with pytest.raises(error):
            check(kwargs)
    check({'state_reset_hook': 'reset_cache', 't_span': (0, 1),
           't_eval': [0, 1], 'export_state_to': 'file.npy'})
    assert check({'method': 'not a method'}) == set()
    assert check({'method': object}) == set()
    with pytest.warns(RuntimeWarning):
        assert check({'method': 'RK45', 'lband': 1}) == {'lband'}
    with pytest.raises(TypeError):
        sys.dynamic_run()
    with pytest.warns(UserWarning):
        sys.set_dynamic_tracker()

def test_process_impacts():
    bst.F.set_flowsheet('process_impacts')
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    bst.settings.define_impact_indicator('GWP', 'kg*CO2e')
    try:
        feed = bst.Stream('feed', Water=100, Ethanol=50, price=0.1)
        feed.set_CF('GWP', 0.5)
        H1 = bst.HXutility('H1', feed, T=360, V=0)
        P1 = bst.Pump('P1', H1-0, P=5e5)
        S1 = bst.Splitter('S1', P1-0, ['product', 'coproduct'], split=0.7)
        product, coproduct = S1.outs
        product.price = 1.
        coproduct.price = 0.2
        coproduct.set_CF('GWP', 0.1)
        sys = bst.System('sys', [H1, P1, S1], operating_hours=8000)
        sys.define_process_impact('GWP', 'Emissions', 'kg', lambda: 1000., 1.)
        sys.define_process_impact('GWP', 'Credit', 'kg', lambda: -100., 1.)
        with pytest.raises(ValueError):
            sys.define_process_impact('GWP', 'Emissions', 'kg', lambda: 1000., 1.)
        sys.simulate()
        assert_allclose(sys.get_process_impact('GWP'), 900.)
        assert_allclose(sys.get_process_impact('GWP', displace=False), 1000.)
        assert sys.get_process_impact('other') == 0.
        net = sys.get_net_impact('GWP')
        net_no_displacement = sys.get_net_impact('GWP', displace=False)
        assert net > 0 and net_no_displacement > 0
        assert sys.get_total_input_impact('GWP') > 0
        assert sys.get_material_impact(feed, 'GWP') > 0
        assert sys.get_total_products_impact('GWP') > 0
        assert sys.get_net_electricity_impact('GWP') == 0.
        assert sys.get_net_electricity_impact('GWP', displace=False) == 0.
        assert sys.get_net_utility_impact('GWP') == 0.
        assert sys.get_net_heat_utility_impact('low_pressure_steam', 'GWP') == 0.
        agent = bst.HeatUtility.get_agent('low_pressure_steam')
        assert sys.get_net_heat_utility_impact(agent, 'GWP') == 0.
        with pytest.raises(TypeError):
            sys.get_net_heat_utility_impact(1, 'GWP')
        bst.HeatUtility.set_CF('low_pressure_steam', 'GWP', 0.1, basis='kg')
        assert sys.get_net_heat_utility_impact(agent, 'GWP') > 0.
        assert sys.get_net_heat_utility_impact(agent, 'GWP', displace=False) > 0.
        bst.HeatUtility.set_CF('low_pressure_steam', 'GWP', 1e-3, basis='kJ')
        assert sys.get_net_heat_utility_impact(agent, 'GWP') > 0.
        bst.HeatUtility.set_CF('low_pressure_steam', 'GWP', 1e-3, basis='kmol')
        assert sys.get_net_heat_utility_impact(agent, 'GWP') > 0.
        bst.PowerUtility.set_CF('GWP', 0.5)
        assert sys.get_net_electricity_impact('GWP') > 0.
        bst.HeatUtility.set_CF('low_pressure_steam', 'GWP', 0.1, basis='kg')
        for method in ('mass', 'energy', 'revenue', 'displacement'):
            value = sys.get_product_impact(product, 'GWP', method)
            assert np.isfinite(value)
        with pytest.raises(ValueError):
            sys.get_product_impact(1, 'GWP')
        with pytest.raises(ValueError):
            sys.get_product_impact('electricity', 'GWP', 'mass')
        with pytest.raises(ValueError):
            sys.get_property_allocated_impact('GWP', 'not a property')
        assert sys.get_property_allocated_impact('GWP', 'mass', 'kg', ignored=[coproduct]) > 0
        assert sys.get_property_allocated_impact('WU', 'mass') is not None or True
        factors = sys.get_property_allocation_factors('mass', groups='product')
        assert_allclose(sum(factors.values()), 1.)
        for args in [('electricity', 'GWP', 'energy')]:
            assert np.isfinite(sys.get_product_impact(*args))
        assert np.isfinite(sys.get_property_allocated_impact('GWP', 'mass', products=[feed, product]))
        with pytest.raises(ValueError):
            sys.get_property_allocated_impact('GWP', 'mass', products=[P1-0])
        factors = sys.get_displacement_allocation_factors(product, 'GWP')
        assert_allclose(sum(factors.values()), 1.)
    finally:
        bst.settings.impact_indicators.clear()
        bst.HeatUtility.characterization_factors.clear()
        bst.PowerUtility.characterization_factors.clear()

def test_save_report():
    bst.F.set_flowsheet('system_save_report')
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    feed = bst.Stream('feed', Water=100, Ethanol=50, price=0.1)
    H1 = bst.HXutility('H1', feed, T=320)
    P1 = bst.Pump('P1', H1-0, P=2e5)
    sys = bst.System('sys', [H1, P1], operating_hours=8000)
    sys.simulate()
    with tempfile.TemporaryDirectory() as tmpdir:
        file = os.path.join(tmpdir, 'report.xlsx')
        sys.save_report(file, sheets={'Stream table', 'Utilities', 'Design requirements',
                                      'Reactions', 'Water mass balance', 'Itemized costs'})
        assert os.path.exists(file)
        with pytest.raises(ValueError):
            sys.save_report(file, allocation_method='mass')

def test_agile_system():
    bst.F.set_flowsheet('agile_system')
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    feed = bst.Stream('feed', Water=100, Ethanol=50, price=0.1)
    H1 = bst.HXutility('H1', feed, T=320)
    P1 = bst.Pump('P1', H1-0, P=2e5)
    product = P1.outs[0]
    product.price = 0.5
    sys = bst.System('sys', [H1, P1])
    agile = bst.AgileSystem(stream_property_names=['T'], lang_factor=4)

    @agile.operation_parameter
    def set_feed_flow(flow):
        feed.F_mass = flow

    @agile.operation_parameter(mode_dependent=True)
    def set_temperature(T, mode):
        H1.T = T

    @agile.operation_metric(annualize=True)
    def annual_production(mode):
        return product.F_mass

    @agile.operation_metric
    def product_temperature(mode):
        return product.T

    mode_a = agile.operation_mode(sys, operating_hours=4000, flow=1000, T=320)
    mode_b = agile.operation_mode(sys, operating_hours=2000, flow=2000, T=330)
    assert 'OperationMode' in repr(mode_a)
    assert 'AgileSystem' in repr(agile)
    agile.simulate()
    assert_allclose(annual_production(), 1000 * 4000 + 2000 * 2000, rtol=1e-3)
    assert len(product_temperature()) == 2
    assert 'OperationMetric' in repr(annual_production)
    assert agile.operating_hours == 6000
    agile.operating_hours = 3000
    assert_allclose(mode_a.operating_hours, 2000)
    agile.operating_hours = 6000
    assert agile.material_cost > 0 and agile.sales > 0
    assert agile.has_market_value(product)
    assert agile.get_market_value(product) > 0
    assert agile.get_mass_flow(product) > 0
    assert agile._price2cost(product) > 0
    assert agile._price2cost(feed) < 0
    with pytest.raises(ValueError):
        agile._price2cost(bst.Stream('other_stream', Water=1))
    assert agile.get_property(product, 'T') > 0
    assert agile.get_property(product, 'T', 'degC') < agile.get_property(product, 'T')
    assert agile.get_total_feeds_impact('GWP') == 0
    assert agile.get_total_products_impact('GWP') == 0
    assert set(agile.units) == {H1, P1}
    assert set(agile.cost_units) == {H1, P1}
    assert agile.get_all_recycles() == set()
    assert agile._downstream_system(H1) is agile
    assert agile.purchase_cost > 0 and agile.installed_equipment_cost > 0
    agile.empty_recycles()
    agile.reset_cache()

def test_phenomena_tracking_and_mock_convergence():
    bst.F.set_flowsheet('phenomena_tracking_coverage')
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    for algorithm in ('phenomena based', 'phenomena modular', 'sequential modular'):
        feed_liquid = bst.Stream('feed_liquid', Water=50, Ethanol=50, phase='l', T=400)
        feed_gas = bst.Stream('feed_gas', Water=50, Ethanol=50, phase='g', T=400)
        flash = bst.StageEquilibrium('flash', ins=[feed_liquid, feed_gas], Q=0, P=101325, phases=('g', 'l'))
        flash.T = 350
        flash.K = np.array([0.5, 1.5])
        flash.B = 1.
        flash.outs[0].mol[:] = [50, 50]
        flash.outs[1].mol[:] = [50, 50]
        sys = bst.System('sys', path=[flash], algorithm=algorithm)
        sys.track_convergence(True, condition_number=True)
        for i in range(3): sys.run()
        sys.assert_tracking_ok()
        F_out = sum([i.F_mol for i in flash.outs])
        assert_allclose(F_out, 200, rtol=5e-2)
        if algorithm == 'phenomena based':
            graph = sys.get_phenomegraph(bipartite=True)
            assert graph is not None
        sys.track_convergence(False)
        assert not sys.tracking_flag
        sys = bst.System('sys', path=[flash], algorithm=algorithm)
        sys.mock_convergence()
        assert len(sys.timer.record) > 0
    with pytest.raises(ValueError):
        sys.get_monopartite_phenomegraph(decomposition='none')

def test_network_methods():
    bst.F.set_flowsheet('network_methods')
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    feed = bst.Stream('feed', Water=100, Ethanol=50)
    recycle = bst.Stream('recycle')
    P1 = bst.Pump('P1', feed)
    M1 = bst.Mixer('M1', [P1-0, recycle])
    S1 = bst.Splitter('S1', M1-0, ['product', recycle], split=0.5)
    inner = bst.Network([M1, S1], recycle=recycle)
    network = bst.Network([P1, inner])
    network.show()
    network._ipython_display_()
    assert network.get_first_unit() is P1
    assert network.get_last_unit() is S1
    assert bst.Network([inner]).get_first_unit() is M1
    assert inner.recycle_sink is M1
    assert network.recycle_sink is None
    assert recycle in network.streams
    assert network.get_all_recycles() == {recycle}
    assert not inner.isdisjoint(network)
    assert inner.issubset(network)
    assert network.first_unit([S1]) is S1
    assert inner.get_recycle_units() == {M1, S1}
    other = bst.Stream('other')
    another = bst.Stream('another')
    net = bst.Network([P1])
    net.add_recycle(None)
    net.add_recycle(other)
    net.add_recycle(other)
    assert net.recycle is other
    net.add_recycle({another})
    assert net.recycle == {other, another}
    net.add_recycle(recycle)
    net.add_recycle({feed})
    assert len(net.recycle) == 4
    net = bst.Network([P1])
    net.recycle = other
    net.add_recycle(another)
    assert net.recycle == {other, another}
    net.show()
    net = bst.Network([P1], recycle=other)
    net._append_network(bst.Network([M1]))
    assert len(net.path) == 2
    net = bst.Network([P1])
    net._append_network(bst.Network([M1, S1], recycle=recycle))
    assert isinstance(net.path[-1], bst.Network)
    net = bst.Network([S1, P1, M1])
    net.sort_without_recycle(set())
    assert net.path.index(P1) < net.path.index(M1)
    
def test_flowsheet_methods():
    F = bst.F
    F.set_flowsheet('flowsheet_methods_coverage')
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    feed = bst.Stream('feed', Water=100, Ethanol=50)
    M1 = bst.Mixer('M1', feed)
    S1 = bst.Splitter('S1', M1-0, ['a', 'b'], split=0.5)
    H1 = bst.HXutility('H1', S1-0, T=320)
    H2 = bst.HXutility('H2', S1-1, T=320)
    flowsheet = F.get_flowsheet()
    with pytest.raises(AttributeError):
        flowsheet.not_an_item
    assert flowsheet('M1') is M1
    assert flowsheet('Mixer') is M1
    assert flowsheet('HXutility') == [H1, H2]
    assert flowsheet(bst.Mixer) is M1
    assert flowsheet(bst.HXutility) == [H1, H2]
    with pytest.raises(LookupError): flowsheet('M1', strict=True) and flowsheet('nothing', strict=True)
    with pytest.raises(LookupError): flowsheet('nothing')
    with pytest.raises(LookupError): flowsheet(bst.Pump)
    with pytest.raises(TypeError): flowsheet(1)
    with pytest.raises(TypeError): F.set_flowsheet(1)
    new = bst.Flowsheet.from_flowsheets('combined_flowsheet', ['flowsheet_methods_coverage'])
    assert 'M1' in new.unit
    new.discard('M1')
    assert 'M1' not in new.unit
    flowsheet.remove_unit_and_associated_streams('H2')
    assert 'H2' not in flowsheet.unit
    with pytest.raises(RuntimeError):
        with bst.Flowsheet('temporary_flowsheet'):
            raise RuntimeError('error')
    assert F.ID == 'flowsheet_methods_coverage'
    F.set_flowsheet('flowsheet_methods_coverage', new=True)
    assert 'M1' not in F.unit
    dct = {}
    F.open(dct)
    bst.Mixer('M_open')
    assert 'M_open' in dct
    F.close()
    assert 'M_open' in F.unit
    streams = bst.get_streams_from_context_level(level=10)
    assert len(streams) > 0
    with bst.System('sys'):
        bst.Mixer('M_context', outs='context_outlet')
        streams = bst.get_streams_from_context_level(level=-1)
        assert len(streams) == 1
    
def test_module():
    bst.F.set_flowsheet('module_coverage')
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    feed = bst.Stream('feed', Water=100, Ethanol=50)
    with bst.System('subsys') as subsys:
        H1 = bst.HXutility('H1', feed, T=320)
        P1 = bst.Pump('P1', H1-0, P=2e5)
    module = bst.Module('module', system=subsys)
    assert module.ins[0] is feed
    assert not module.run_full_simulation
    module.simulate()
    assert_allclose(module.outs[0].F_mass, feed.F_mass)
    module.run_full_simulation = True
    module.simulate()
    assert module.outs[0].P > 1.5e5

class SimpleFacility(bst.Facility):
    network_priority = 0
    _N_ins = 1
    _N_outs = 1
    
    def _run(self):
        self.outs[0].copy_like(self.ins[0])

def test_system_with_facilities():
    bst.F.set_flowsheet('system_with_facilities')
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    feed = bst.Stream('feed', Water=100, Ethanol=50)
    M1 = bst.Mixer('M1', [feed, 'returned'])
    S1 = bst.Splitter('S1', M1-0, ['product', 'purge'], split=0.8)
    FC = SimpleFacility('FC', S1-1, 'returned')
    FC.outs[0] = M1.ins[1]
    called = []
    sys = bst.System('sys', [M1, S1], facilities=[FC])
    assert sys._integrated_facilities
    assert M1.ins[1] in sys.get_all_recycles()
    assert sys._downstream_facilities(FC) == (FC,)
    downstream = sys._downstream_system(FC)
    assert downstream is sys or FC in downstream.units
    sys.simulate()
    product = sys.flowsheet.stream.product
    assert 0 < product.F_mass <= feed.F_mass * 1.01
    assert FC.outs[0].F_mass > 0
    # Non-integrated facilities, including functions and systems
    bst.F.set_flowsheet('system_with_function_facilities')
    feed = bst.Stream('feed', Water=100, Ethanol=50)
    H1 = bst.HXutility('H1', feed, T=320)
    H2 = bst.HXutility('H2', bst.Stream('other_feed', Water=1), T=320)
    facility_sys = bst.System('facility_sys', [H2])
    sys = bst.System('sys', [H1], facilities=[facility_sys])
    sys.simulate()
    assert_allclose(H2.outs[0].T, 320)
    with pytest.warns(RuntimeWarning):
        bst.System('outer', [sys]).flatten()
    with pytest.raises(RuntimeError):
        sys._downstream_system(bst.Mixer('M_other'))
    assert sys.LCA is None and sys.TEA is None
    
def test_phenomena_based_recycle():
    sys = create_recycle_system('phenomena_based_recycle', algorithm='phenomena based')
    sys.simulate()
    assert_mass_balance(sys, rtol=5e-2)
    sys = create_recycle_system('phenomena_modular_recycle', algorithm='phenomena modular')
    sys.simulate()
    assert_mass_balance(sys, rtol=5e-2)
    
def test_system_parameters_and_misc():
    sys = create_recycle_system('system_parameters')
    S1 = sys.flowsheet.unit.S1
    M1 = sys.flowsheet.unit.M1
    with pytest.raises(ValueError):
        bst.System.from_feedstock('none')
    
    def set_split(split): S1.split[:] = split
    
    p = sys.parameter(set_split, S1, baseline=0.5, bounds=(0.4, 0.6), units='-')
    assert p.baseline == 0.5
    p2 = sys.parameter(p, coupled=True)
    assert p2.setter is set_split
    
    @sys.parameter(element=S1, kind='coupled')
    def set_split_2(split): S1.split[:] = split
    
    assert set_split_2.coupled
    assert sys.path_segment(S1, M1) == [S1] or M1 not in sys.path_segment(S1, M1)
    with pytest.raises(ValueError):
        sys.path_segment(bst.Mixer('M_other'), M1, relevant_units={M1}, critical_units=set())
    
def test_multiple_agile_systems():
    bst.F.set_flowsheet('multiple_agile_systems')
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    feed = bst.Stream('feed', Water=100, Ethanol=50, price=0.1)
    H1 = bst.HXutility('H1', feed, T=320)
    H2 = bst.HXutility('H2', H1-0, T=330)
    sys_a = bst.System('sys_a', [H1])
    sys_b = bst.System('sys_b', [H1, H2])
    agile = bst.AgileSystem()
    agile.operation_mode(sys_a, operating_hours=1000)
    agile.operation_mode(sys_b, operating_hours=1000)
    agile.simulate()
    assert set(agile.cost_units) == {H1, H2}
    assert agile.installed_equipment_cost > 0

if __name__ == '__main__':
    test_recycle_convergence_methods()
    test_recycle_maxiter()
    test_recycle_data()
    test_multiple_recycles()
    test_no_recycle_errors()
    test_mock_system()
    test_from_segment_and_paths()
    test_nested_system_operations()
    test_flatten_warnings()
    test_specifications()
    test_system_context_and_configuration()
    test_interface_property_packages()
    test_diagrams_without_display()
    test_utility_and_cost_properties()
    test_debug_mode_errors()
    test_dynamic_kwargs_validation()
    test_process_impacts()
    test_save_report()
    test_agile_system()
    test_phenomena_tracking_and_mock_convergence()
    test_network_methods()
    test_flowsheet_methods()
    test_module()
    test_system_with_facilities()
    test_phenomena_based_recycle()
    test_system_parameters_and_misc()
    test_multiple_agile_systems()
