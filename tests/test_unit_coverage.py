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
from biosteam._unit import (
    SuperpositionInlet, SuperpositionOutlet, DesignInput, NumericalDesignInput, CategoricalDesignInput,
    LogicalDesignInput, UnitInheritanceError,
)

def test_design_input_edge_cases():
    bst.settings.set_thermo(['Water'], cache=True)
    # Callable name with extra arguments is invalid
    with pytest.raises(ValueError):
        bst.design_input(lambda self: 1, units='m')
    # Copy of a design input keeps name and return type
    class A(bst.Unit):
        x = bst.design_input('X', 'm', lb=0, ub=10, default=1.)
        mode = bst.design_input(catalog={'Only': 1.0}, default='Only')
        empty_catalog = bst.design_input(catalog={})
    copy = A.x.copy(default=2.)
    assert copy.default == 2. and copy.__name__ == 'x'
    # Numerical validation errors with values that cannot be formatted
    class Weird(float):
        def __format__(self, spec): raise ValueError()
    x = A.x
    with pytest.raises(AttributeError):
        x.validate(Weird(-1.))
    with pytest.raises(AttributeError):
        x.validate(Weird(20.))
    with pytest.raises(AttributeError):
        x.validate('not a number')
    # Categorical docs and validation with single/empty catalog
    assert "only valid choice" in A.mode._make_doc()
    assert "No valid choice" in A.empty_catalog._make_doc()
    assert A.mode.validate('only') == 'Only'
    with pytest.raises(AttributeError):
        A.mode.validate('other')
    with pytest.raises(AttributeError):
        A.empty_catalog.validate('other')
    assert 'catalog' not in repr(CategoricalDesignInput('Z', None, None, None, None, None, None, None, None))
    # Generic design input validation
    di = DesignInput('Y')
    di._return_type = int
    assert di.validate(2.0) == 2
    with pytest.raises(AttributeError):
        di.validate(2.5)
    di._return_type = None
    assert di.validate('anything') == 'anything'
    # Only getter is invalid
    with pytest.raises(RuntimeError):
        class B(bst.Unit):
            y = bst.design_input(lambda self: 1)
            y.fget = lambda self: 1
            y.fset = None
    # Logical design input validation
    class C(bst.Unit):
        flag = bst.design_input(default=False)
    with pytest.raises(AttributeError):
        C.flag.validate(1)
    with pytest.raises(AttributeError):
        del C().flag

def test_categorical_design_input_cost_factor():
    bst.settings.set_thermo(['Water'], cache=True)
    class Vessel(bst.Unit):
        _F_BM_default = {'Vessel': 2.}
        material = bst.design_input(
            purchase_item='Vessel', catalog_type='Material',
            catalog={'Carbon steel': 1.0, 'Stainless steel': 2.5},
        )
        def _cost(self):
            self.baseline_purchase_costs['Vessel'] = 1000.
    V = Vessel(None, bst.Stream(None, Water=1))
    V.material = 'stainless steel'
    assert V.F_M['Vessel'] == 2.5
    V.material = None
    assert 'Vessel' not in V.F_M
    V.simulate()
    assert_allclose(V.purchase_cost, 1000., rtol=1e-2)

def test_unit_inheritance_errors():
    bst.settings.set_thermo(['Water'], cache=True)
    with pytest.raises(UnitInheritanceError):
        class BadRun(bst.Unit):
            def run(self): pass
    with pytest.raises(UnitInheritanceError):
        class NoRun(bst.Unit):
            _N_ins = 2
            _N_outs = 1
    # Ticket names inferred from line
    names = {
        'MyCentrifuge': 'C', 'MyDistillation': 'D', 'MyEvaporator': 'E',
        'MyFlash': 'F', 'MyCooler': 'H', 'MyCompressor': 'K', 'MyTurbine': 'Ʞ',
        'MyMixer': 'M', 'MyPump': 'P', 'MyReactor': 'R', 'MySplitter': 'S',
        'MyTank': 'T', 'MyValve': 'V', 'MySpecification': 'PS', 'MyThing': 'U',
    }
    for name, ticket in names.items():
        cls = type(name, (bst.Unit,), {})
        assert cls.ticket_name == ticket
    # _N_heat_utilities is deprecated
    with pytest.warns(DeprecationWarning):
        class OldUnit(bst.Unit):
            _N_heat_utilities = 1

def test_unit_feed_and_product_properties():
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    bst.F.set_flowsheet('unit_feed_product')
    M = bst.Mixer('M1', ins=[bst.Stream('a', Water=1), bst.Stream('b', Ethanol=1)])
    with pytest.raises(AttributeError): M.feed
    with pytest.raises(AttributeError): M.feed = bst.Stream()
    S = bst.Splitter('S1', ins=bst.Stream('c', Water=1), split=0.5)
    with pytest.raises(AttributeError): S.product
    with pytest.raises(AttributeError): S.product = bst.Stream()
    new = bst.Stream('d', Water=2)
    S.feed = new
    assert S.feed is new
    P = bst.Pump('P1', ins=new)
    out = bst.Stream('e')
    P.product = out
    assert P.product is out
    class Empty(bst.Unit):
        _N_ins = _N_outs = 0
        _ins_size_is_fixed = _outs_size_is_fixed = False
        def _run(self): pass
    U = Empty('U1')
    with pytest.raises(AttributeError): U.feed
    with pytest.raises(AttributeError): U.feed = bst.Stream()
    with pytest.raises(AttributeError): U.product
    with pytest.raises(AttributeError): U.product = bst.Stream()

def test_auxiliary_unit_methods():
    bst.settings.set_thermo(['Water'], cache=True)
    bst.F.set_flowsheet('unit_auxiliary')
    class Parent(bst.Unit):
        _N_ins = 1
        _N_outs = 1
        auxiliary_unit_names = ('pump', 'mixers')
        def _init(self):
            self.pump = self.auxiliary('pump', bst.Pump, ins=self.ins[0])
            self.mixers = [bst.Mixer(None), bst.Mixer(None)]
        def _run(self):
            self.outs[0].copy_like(self.ins[0])
        def _design(self):
            self.pump.simulate()
    feed = bst.Stream('feed', Water=100)
    P = Parent('P', ins=feed)
    nested = P.nested_auxiliary_units
    assert P.pump in nested
    names = [i for i, j in P.get_nested_auxiliary_units_with_names()]
    assert 'pump' in names and 'mixers[0]' in names
    assert P.get_nested_auxiliary_units_with_names(depth=0) == []
    assert [i for i, j in P._diagram_auxiliary_units_with_names()][-1] == 'mixers[1]'
    # Register auxiliaries
    m = bst.Mixer(None)
    P.register_auxiliary(m, 'extra')
    assert 'extra' in P.auxiliary_unit_names
    P.register_auxiliary([bst.Mixer(None)], 'extras')
    assert 'extras' in P.auxiliary_unit_names
    with pytest.raises(ValueError):
        P.register_auxiliary([bst.Mixer(None)])
    with pytest.raises(ValueError):
        P.register_auxiliary(1)
    sys = bst.System(None, [bst.Mixer('MA')])
    P.register_auxiliary(sys, 'subsys')
    assert 'MA' in P.auxiliary_unit_names
    # Auxlets
    s = P.auxlet(None)
    assert s.source is P and s.sink is P
    s = P.auxlet('aux_stream')
    assert s.source is P
    s = P.auxlet(P.outs[0])
    assert isinstance(s, SuperpositionOutlet)
    s = P.auxlet(P.ins[0])
    assert isinstance(s, SuperpositionInlet)
    other = bst.Stream(None)
    assert P.auxlet(other).source is P
    assert P.auxin(None).sink is P
    # Auxiliary stack
    P.stack = []
    P.auxiliary('stack', bst.Mixer)
    assert len(P.stack) == 1
    P.simulate()
    assert P.purchase_cost > 0

def test_unit_node_and_paths():
    bst.settings.set_thermo(['Water'], cache=True)
    bst.F.set_flowsheet('unit_paths')
    feed = bst.Stream('feed', Water=10)
    P1 = bst.Pump('P1', feed)
    H1 = bst.HXutility('H1', P1-0, T=320)
    P2 = bst.Pump('P2', H1-0)
    assert H1 in P2.path_from(P1)
    assert P1 in P2.path_from(P1, inclusive=True)
    assert H1 in P1.path_until(P2)
    assert P2 in P1.path_until(P2, inclusive=True)
    P1.run_until(P2, inclusive=True)
    assert_allclose(P2.outs[0].T, 320, rtol=1e-3)
    try:
        bst.preferences.minimal_nodes = True
        assert P1.get_node()
    finally:
        bst.preferences.minimal_nodes = False
    P1._get_tooltip_string = None # Tooltip fails silently
    assert P1.get_node()
    del P1._get_tooltip_string

def test_unit_disconnect_insert_and_replace():
    bst.settings.set_thermo(['Water'], cache=True)
    bst.F.set_flowsheet('unit_connections_coverage')
    feed = bst.Stream('feed', Water=10)
    P1 = bst.Pump('P1', feed)
    H1 = bst.HXutility('H1', P1-0, T=320)
    P2 = bst.Pump('P2', H1-0)
    # Disconnect specific inlets/outlets and join ends
    H1.disconnect(inlets=[H1.ins[0]], outlets=[H1.outs[0]], join_ends=True)
    assert P2.ins[0] is P1.outs[0]
    with pytest.raises(ValueError):
        M = bst.Mixer('M', ins=[bst.Stream(), bst.Stream()], outs='m_out')
        M.disconnect(join_ends=True)
    # Insert unit with explicit outlet and inlet
    H2 = bst.HXutility('H2', T=330)
    H2.insert(P1-0, outlet=0)
    assert P2.ins[0] is H2.outs[0]
    with pytest.raises(ValueError):
        bst.HXutility('H3', T=330).insert(P1-0, outlet=bst.Stream())
    with pytest.raises(ValueError):
        bst.HXutility('H4', T=330).insert(H2-0, inlet=bst.Stream())
    S = bst.Splitter('S', split=0.5)
    with pytest.raises(ValueError):
        S.insert(H2-0)
    # Replace with nothing
    H5 = bst.HXutility('H5', ins=P2-0, outs='h5_out', T=300)
    P3 = bst.Pump('P3', H5-0)
    H5.replace_with(discard=True)
    assert P3.ins[0] is P2.outs[0]
    H6 = bst.HXutility('H6', ins=bst.Stream('lonely'), T=300)
    P4 = bst.Pump('P4', H6-0)
    H6.replace_with()
    assert P4.ins[0].ID == 'lonely'
    H7 = bst.HXutility('H7', T=300)
    H7.replace_with(bst.HXutility('H8', T=300), discard=True)

def test_unit_piping_operators():
    bst.settings.set_thermo(['Water'], cache=True)
    bst.F.set_flowsheet('unit_piping_operators')
    M = bst.Mixer('M', ins=[bst.Stream('s1'), bst.Stream('s2')])
    S = bst.Splitter('S', outs=('o1', 'o2'), split=0.5)
    s = bst.Stream('s3')
    assert (S - s) is S and S.outs[0] is s
    outs = [bst.Stream('o3'), bst.Stream('o4')]
    assert (S - outs) is S
    assert (S - [0, 1]) == outs
    assert ([0, 1] - M) == [M.ins[0], M.ins[1]]
    ins = [bst.Stream('i1'), bst.Stream('i2')]
    assert (ins - M) is M
    with pytest.raises(ValueError):
        M.__rsub__(1.5)
    assert (S ** 0) is S.outs[0]

def test_unit_utilities_and_credits():
    bst.settings.set_thermo(['Water'], cache=True)
    bst.F.set_flowsheet('unit_utilities_coverage')
    feed = bst.Stream('feed', Water=100)
    H = bst.HXutility('H', feed, T=350)
    with pytest.raises(ValueError):
        H.define_utility('not a utility', feed)
    with pytest.raises(ValueError):
        H.define_credit('not a credit', feed)
    bst.stream_prices['Dummy'] = 0.1
    try:
        with pytest.raises(ValueError):
            H.define_utility('Dummy', bst.Stream())
        with pytest.raises(ValueError):
            H.define_credit('Dummy', bst.Stream())
        H.define_utility('Dummy', feed)
        H.define_credit('Dummy', H.outs[0], fraction=lambda: 0.5)
        H.simulate()
        assert H.get_inlet_cost_flows()['Dummy'] > 0
        assert H.get_outlet_revenue_flows()['Dummy'] > 0
        results = H.results(include_zeros=True)
        assert results is not None
        assert H.results(with_units=False, include_zeros=True) is not None
        assert H.results(basis='SI') is not None
    finally:
        del bst.stream_prices['Dummy']
    assert H.net_heating_duty > 0 and H.net_cooling_duty == 0
    H.add_power_utility(10)
    H.add_power_utility(-5)
    assert_allclose(H.net_power, 5)
    hu = H.create_heat_utility()
    assert hu in H.heat_utilities
    H.empty()
    assert not H.heat_utilities and H.outs[0].isempty()

def test_unit_results_with_design_inputs_and_units():
    bst.settings.set_thermo(['Water'], cache=True)
    bst.F.set_flowsheet('unit_results_basis')
    class Vessel(bst.Unit):
        _units = {'Pressure': 'psig', 'Volume': 'm3'}
        tau = bst.design_input('Residence time', 'hr', lb=0, default=1)
        def _design(self):
            self.design_results['Pressure'] = 10.
            self.design_results['Volume'] = self.tau * self.ins[0].F_vol
        def _cost(self):
            self.baseline_purchase_costs['Vessel'] = 1000. * self.design_results['Volume']
            self.F_BM['Vessel'] = 2.
    V = Vessel('V', bst.Stream('feed', Water=100))
    V.tau = 2
    V.simulate()
    df = V.results(basis='SI')
    assert ('Design inputs', 'Residence time') in df.index
    assert V.results(with_units=False)[('Design results', 'Volume')] > 0
    assert V.installed_cost > V.purchase_cost > 0
    assert V.baseline_purchase_cost > 0
    V._reevaluate()
    assert V.purchase_cost > 0

def test_unit_parallel_and_agile_design():
    bst.settings.set_thermo(['Water'], cache=True)
    bst.F.set_flowsheet('unit_parallel_agile')
    class Vessel(bst.Unit):
        _F_BM_default = {'Vessel': 2.}
        _materials_and_maintenance = frozenset(['Vessel'])
        _max_agile_design = ('Volume',)
        auxiliary_unit_names = ('pump',)
        def _init(self):
            self.pump = self.auxiliary('pump', bst.Pump, ins=self.ins[0])
        def _design(self):
            self.design_results['Volume'] = self.ins[0].F_vol
            self.pump.simulate()
        def _cost(self):
            self.baseline_purchase_costs['Vessel'] = 1000. * self.design_results['Volume']
    feed = bst.Stream('feed', Water=100)
    V = Vessel('V', feed)
    V.simulate()
    dc1 = V.get_design_and_capital()
    assert dc1.equipment_lifetime == V.equipment_lifetime
    feed.imol['Water'] = 200
    V.simulate()
    dc2 = V.get_design_and_capital()
    agile = V.get_agile_design_and_capital([dc1, dc2])
    assert agile.installed_costs['Vessel'] > 0
    # Parallel units multiply costs
    V.parallel['self'] = 2
    V._costs_loaded = False
    V.purchase_costs['Vessel'] = V.installed_costs['Vessel'] = 1.
    V.baseline_purchase_costs.pop('Pump - Pump')
    V.baseline_purchase_costs.pop('Pump - Motor')
    V._load_costs()
    assert V.purchase_costs['Vessel'] == 2.
    V.parallel['pump'] = 10**6
    V._costs_loaded = False
    V.baseline_purchase_costs.pop('Pump - Pump')
    V.baseline_purchase_costs.pop('Pump - Motor')
    with pytest.raises(RuntimeError):
        V._load_costs()

def test_unit_checks_and_balances():
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    bst.F.set_flowsheet('unit_checks')
    class Bad(bst.Unit):
        def _run(self):
            self.outs[0].copy_like(self.ins[0])
            self.add_power_utility(1)
    U = Bad('U', bst.Stream('feed', Water=1, Ethanol=1))
    U._run()
    with pytest.warns(RuntimeWarning):
        U._check_run()
    with pytest.raises(UnitInheritanceError):
        U._check_setup()
    U.power_utility.empty()
    # Auxiliary heat utilities manually added
    class Bad2(bst.Unit):
        auxiliary_unit_names = ('hx',)
        def _init(self):
            self.hx = self.auxiliary('hx', bst.HXutility, ins=self.ins[0], T=350)
        def _design(self):
            self.hx.simulate()
            self.heat_utilities.extend(self.hx.heat_utilities)
    U2 = Bad2('U2', bst.Stream('feed2', Water=1))
    with pytest.raises(UnitInheritanceError):
        U2.simulate()
    U.run()
    assert abs(U.mass_balance_error()) < 1e-6
    assert all([abs(i) < 1e-6 for i in U.atomic_balance_error().values()])
    assert_allclose(U.z_mass_in.sum(), 1)
    assert_allclose(U.z_mass_out.sum(), 1)
    assert_allclose(U.z_vol_in.sum(), 1)
    assert_allclose(U.z_vol_out.sum(), 1)
    assert U.F_vol_out > 0
    assert set(U.get_available_chemicals()) == set(U.chemicals)
    U.specifications = None
    assert U.specifications == []
    U.specifications = []
    with pytest.raises(AttributeError):
        U.specification
    with pytest.raises(AttributeError):
        U.specification = None

def test_unit_reset_thermo():
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    bst.F.set_flowsheet('unit_reset_thermo')
    # Note: streams do not implement `_reset_thermo`, so only units with
    # missing streams can reset thermo.
    R = bst.MixTank('R', ins=None, outs=None)
    R.rxn = bst.Reaction('Water -> Ethanol', 'Water', 0.1)
    R.others = {'a': [bst.Stream(None)]}
    thermo = bst.Thermo(['Water', 'Ethanol', 'Methanol'], cache=True)
    R._reset_thermo(thermo)
    assert R.thermo is thermo and R.rxn.chemicals is thermo.chemicals

def test_phenomena_oriented_methods():
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    bst.F.set_flowsheet('unit_phenomena')
    feed = bst.Stream('feed', Water=1, Ethanol=1)
    H = bst.HXutility('H', feed, T=320)
    H.run()
    H._recycle_system = None
    H._dmol = np.zeros(2)
    mol_error, T_error = H._simulation_error()
    assert mol_error < 1e-6 and T_error < 1e-6
    H.specified_variable = 'Q'
    coefficients = {}
    assert H._update_energy_coefficient(H.outs[0], coefficients) < 0
    H.specified_variable = 'T'
    assert H._update_energy_coefficient(H.outs[0], coefficients) == 0
    H._update_variable('T', 330)
    assert H.outs[0].T == 330
    with pytest.raises(ValueError):
        H._update_variable('Q', 0)
    eqs = H._create_material_balance_equations(False)
    assert len(eqs) == 1
    eqs = H._create_material_balance_equations(True)
    assert len(eqs) == 1
    eqs = H._create_bulk_balance_equations()
    assert len(eqs) == 1
    H.specified_variable = 'Q'
    H._duty = H.Hnet
    eqs = H._create_energy_balance_equations()
    assert len(eqs) == 1
    H._dmol = np.ones(2)
    eqs = H._create_energy_balance_equations()
    assert len(eqs) == 1
    H.specified_variable = 'T'
    assert H._create_energy_balance_equations() == []
    assert H.node_tag == 'H'
    H.material_balance_specifications_nodes = []
    H.specification_equations = []
    assert H._collect_specification_errors() == []
    assert H.get_E_node(feed) is None

if __name__ == '__main__':
    test_design_input_edge_cases()
    test_categorical_design_input_cost_factor()
    test_unit_inheritance_errors()
    test_unit_feed_and_product_properties()
    test_auxiliary_unit_methods()
    test_unit_node_and_paths()
    test_unit_disconnect_insert_and_replace()
    test_unit_piping_operators()
    test_unit_utilities_and_credits()
    test_unit_results_with_design_inputs_and_units()
    test_unit_parallel_and_agile_design()
    test_unit_checks_and_balances()
    test_unit_reset_thermo()
    test_phenomena_oriented_methods()
