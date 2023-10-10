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
from biosteam._network import Network
from numpy.testing import assert_allclose

def test_auxiliary_unit_owners():
    with bst.System() as sys:
        class NewUnit(bst.Unit):
            auxiliary_unit_names = ('mixer',)
            def __init__(self, *args, **kwargs):
                super().__init__(*args, **kwargs)
                self.mixer = bst.Mixer()
                
        bst.settings.set_thermo(['Water'], cache=True)
        unit = NewUnit()
    # Once unit is in a system, auxiliary units must have an owner
    assert unit.mixer.owner is unit

def test_unit_inheritance_setup_method():
    class NewUnit(bst.Unit):
        def _setup(self):
            pass
    
    bst.settings.set_thermo(['Water'], cache=True)
    U1 = NewUnit()
    U1.add_specification(lambda: None)
    U1.simulate()

def test_process_specifications_linear():
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    with bst.System() as sys:
        T1 = bst.StorageTank('T1', ins=bst.Stream(Water=1000))
        H1 = bst.HXutility('H1', ins=T1-0, T=320)
        T2 = bst.StorageTank('T2', ins=bst.Stream(Ethanol=1000))
        H2 = bst.HXutility('H2', ins=T2-0, T=320)
        T3 = bst.StorageTank('T3', ins=bst.Stream(Water=10))
        H3 = bst.HXutility('H3', ins=T3-0, T=320)
        M1 = bst.Mixer('M1', ins=[H1-0, H2-0, H3-0])
        H4 = bst.HXutility('H4', ins=M1-0, T=350)
    
    # Specification impacting upstream units
    @M1.add_specification(run=True, impacted_units=[T1])
    def adjust_flow_rate():
        water = T1.ins[0]
        ethanol = T2.ins[0]
        water.F_mass = ethanol.F_mass
    sys.simulate()
    assert M1.specifications[0].path == ()
    assert sys._to_network() == Network(
        [T2,
         H2,
         T3,
         H3,
         Network(
            [M1,
             T1,
             H1],
            recycle=H1-0),
         H4]
    )
    assert_allclose(H4.outs[0].mol, sum([i.ins[0].mol for i in (T1, T2, T3)]))
    
    # Specification impacting units in parallel (neither upstream nor downstream units).
    # System simulation order must switch
    old_path_length = len(sys.unit_path)
    M1.specifications.clear() # Remove specification
    tanks = [i for i in sys.units if isinstance(i, bst.Tank)]
    first_tank, mid_tank, last_tank = tanks
    @last_tank.add_specification(run=True, impacted_units=[first_tank, mid_tank])
    def adjust_flow_rate():
        mid_tank.ins[0].F_vol = last_tank.ins[0].F_vol
        first_tank.ins[0].F_mol = 1.0
    
    @mid_tank.add_specification(run=True, impacted_units=[first_tank])
    def adjust_flow_rate():
        first_tank.ins[0].F_vol = mid_tank.ins[0].F_vol
    
    sys.simulate()
    # No upstream unit simulations
    assert mid_tank.specifications[0].path == []
    assert last_tank.specifications[0].path == []
    # Order of tank simulation is now reversed
    assert [last_tank, mid_tank, first_tank] == [i for i in sys.units if isinstance(i, bst.Tank)] 
    # 2 hidden connection source and 2 sinks
    assert len(sys.unit_path) == old_path_length
    assert (H4.outs[0].mol == sum([i.ins[0].mol for i in (T1, T2, T3)])).all()
    assert sys.path[0] is last_tank
    
    # Specification impacting downstream units (it doesn't matter).
    last_tank.specifications.clear() # Remove specification
    mid_tank.specifications.clear()
    
    @T1.add_specification(run=True, impacted_units=[H3])
    def adjust_flow_rate():
        H3.T = 320
    
    sys.simulate()
    assert H3.outs[0].T == H3.T
    assert T1.specifications[0].path == []
    
    # Two overlapping parallel specifications should make
    # just one parallel specification and one upstream specification
    T1.specifications.clear()
    H1.add_specification(lambda: None, run=True, impacted_units=[T2])
    H2.add_specification(lambda: None, run=True, impacted_units=[T1])
    sys.simulate()
    assert sys.unit_path.index(H1) < sys.unit_path.index(T2)
    assert H1.specifications[0].path == ()
    assert H2.specifications[0].path == ()
    # Net simulation order is ..., T2, H2, T1, H1, T2, H2, ...
    
def test_process_specifications_with_recycles():
    bst.F.set_flowsheet('bifurcated_recycle_loops')
    bst.settings.set_thermo(['Water'], cache=True)
    feed_a = bst.Stream('feed_a', Water=10)
    water_a = bst.Stream('water_a', Water=10)
    recycle_a = bst.Stream('recycle_a')
    P1_a = bst.Pump('P1_a', feed_a)
    P2_a = bst.Pump('P2_a', water_a)
    S1_a = bst.Splitter('S1_a', P2_a-0, split=0.5)
    M1_a = bst.Mixer('M1_a', [P1_a-0, S1_a-1, recycle_a])
    S2_a = bst.Splitter('S2_a', M1_a-0, split=0.5)
    M2_a = bst.Mixer('M2_a', [S1_a-0, S2_a-1])
    S3_a = bst.Splitter('S3_a', M2_a-0, ['', recycle_a], split=0.5)
    feedstock_b = bst.Stream('feedstock_b', Water=1000)
    recycle_b = bst.Stream('recycle_b')
    product_b = bst.Stream('product_b')
    side_product_b = bst.Stream('side_product_b')
    P1_b = bst.Pump('P1_b', feedstock_b)
    P2_b = bst.Pump('P2_b', S2_a-0)
    S1_b = bst.Splitter('S1_b', P2_b-0, split=0.5)
    M1_b = bst.Mixer('M1_b', [P1_b-0, S1_b-1, recycle_b])
    S2_b = bst.Splitter('S2_b', M1_b-0, [side_product_b, ''], split=0.5)
    M2_b = bst.Mixer('M2_b', [S1_b-0, S2_b-1, S3_a-0])
    S3_b = bst.Splitter('S3_b', M2_b-0, [product_b, recycle_b], split=0.5)
    recycle_loop_sys = bst.F.create_system('recycle_loop_sys')
    assert recycle_loop_sys.unit_path == [
        P1_b, P1_a, P2_a, S1_a, M1_a, S2_a, M2_a, 
        S3_a, P2_b, S1_b, M1_b, S2_b, M2_b, S3_b,
    ]
    recycle_loop_sys.simulate()
    x_nested_solution = np.vstack([recycle_a.mol, recycle_b.mol])
    # Test parallel process specification
    P2_b.add_specification(f=lambda: None, run=True, impacted_units=[P1_b])
    recycle_loop_sys.simulate()
    assert recycle_loop_sys.unit_path == [ 
        # Reordered path
        P1_a, P2_a, S1_a, M1_a, S2_a, M2_a, S3_a, 
        P2_b, S1_b, P1_b, M1_b, S2_b, M2_b, S3_b,
    ]
    x_reconfigured_solution = np.vstack([recycle_a.mol, recycle_b.mol])
    assert_allclose(x_nested_solution, x_reconfigured_solution, rtol=2e-2)
    # Test upstream process specification
    M1_b.add_specification(lambda: None, run=True, impacted_units=[S2_a, S3_a])
    recycle_loop_sys.simulate()
    # Make sure path includes upstream subsystem and other units in recycle loop
    assert recycle_loop_sys._to_network() == Network(
        [P1_a,
         P2_a,
         S1_a,
         Network(
            [Network(
                [P1_b,
                 Network(
                    [M1_b,
                     S2_a,
                     M2_a,
                     S2_b,
                     S3_a,
                     M2_b,
                     S3_b,
                     M1_a],
                    recycle=M1_a-0),
                 P2_b,
                 S1_b],
                recycle=P1_b-0)],
            recycle=M1_b-0)])
    assert M1_b.specifications[0].path == ()
    
def test_unit_connections():
    from biorefineries import sugarcane as sc
    sc.load()
    f = sc.flowsheet
    u = f.unit
    all_units = set(sc.sys.units).difference(sc.sys.facilities)
    upstream_units = u.R301.get_upstream_units()
    downstream_units = u.R301.get_downstream_units()
    assert u.R301.neighborhood(1) == {u.T301, u.D301, u.S302, u.H301}
    assert u.R301.neighborhood(2) == {u.S302, u.C301, u.H301, u.M302, u.D301, u.R301,
                                    u.T301, u.M301}
    assert u.R301.neighborhood(200) == u.R301.neighborhood(1000) == all_units
    
    recycle_units = set(sc.sys.find_system(u.R301).units)
    assert recycle_units == upstream_units.intersection(downstream_units)
    ins = tuple(u.R301.ins)
    outs = tuple(u.R301.outs)
    u.R301.disconnect()
    assert not any(u.R301.ins + u.R301.outs)
    ins - u.R301 - outs
    
    class DummyUnit(bst.Unit, isabstract=True, new_graphics=False):
        _N_ins = 2; _N_outs = 2
        
    unit = DummyUnit(None)
    unit.take_place_of(u.R301)
    assert tuple(unit.ins + unit.outs) == (ins + outs)
    assert not any(u.R301.ins + u.R301.outs)
    u.R301.take_place_of(unit)

def test_unit_graphics():
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    M = bst.Mixer(None, outs=None)
    assert M._graphics.get_inlet_options(M, 2) == {'headport': 'c'}
    assert M._graphics.get_inlet_options(M, 100) == {'headport': 'c'}
    assert M._graphics.get_outlet_options(M, 0) == {'tailport': 'e'}
    
    S = bst.Splitter(None, outs=None, split=0.5)
    
    GraphicsWarning = bst.exceptions.GraphicsWarning
    with pytest.warns(GraphicsWarning):
        assert S._graphics.get_inlet_options(S, 1) == {'headport': 'c'}
    
    with pytest.warns(GraphicsWarning):
        assert M._graphics.get_outlet_options(M, 1) == {'tailport': 'c'}

def test_cost_decorator():
    from biosteam.units.decorators import cost
    bst.settings.set_thermo(['Water'], cache=True)
    @cost('Flow rate', CE=bst.settings.CEPCI, cost=1, n=0.6, lb=2, ub=10, units='kg/hr', BM=2.)
    class A(bst.Unit): pass
    
    feed = bst.Stream('feed', Water=1, units='kg/hr')
    
    # Test when size factor is under lower bound
    A1 = A('A1', ins=feed)
    A1.simulate()
    assert_allclose(A1.purchase_cost, 2 ** 0.6)
    assert_allclose(A1.installed_cost, 2 ** 1.6)
    
    # Test when size factor is between lower and upper bound
    feed.imass['Water'] = 5
    A1.simulate()
    assert_allclose(A1.purchase_cost, 5 ** 0.6)
    assert_allclose(A1.installed_cost, 2 * 5 ** 0.6)
    
    # Test when size factor is above upper bound
    feed.imass['Water'] = 20
    A1.simulate()
    assert_allclose(A1.purchase_cost, 2 * 10 ** 0.6)
    assert_allclose(A1.installed_cost, 4 * 10 ** 0.6)
    
def test_equipment_lifetimes():
    from biorefineries.sugarcane import create_tea
    bst.settings.set_thermo(['Water'], cache=True)
    
    class A(bst.Unit):
        _F_BM_default = {
            'Equipment A': 2,
            'Equipment B': 3,
        }
        _default_equipment_lifetime = {
            'Equipment A': 10,
            'Equipment B': 5,
        }
        def _cost(self):
            purchase_costs = self.baseline_purchase_costs
            purchase_costs['Equipment A'] = 1e6
            purchase_costs['Equipment B'] = 1e5
            
    class B(bst.Unit):
        _F_BM_default = {
            'Equipment A': 2,
        }
        _default_equipment_lifetime = 15
        def _cost(self):
            self.baseline_purchase_costs['Equipment A'] = 1e6

    class C(bst.Unit):
        _F_BM_default = {
            'Equipment A': 4,
        }
        def _cost(self):
            self.baseline_purchase_costs['Equipment A'] = 1e6
            
    @bst.decorators.cost('Flow rate', units='kmol/hr', S=1, BM=3,
                         cost=1e6, n=0.6, CE=bst.CE, lifetime=8)
    class D(bst.Unit): pass
    
    @bst.decorators.cost('Flow rate', units='kmol/hr', S=1, BM=4,
                         cost=1e6, n=0.6, CE=bst.CE, lifetime=20)
    class E(bst.Unit): pass
    
    D_feed = bst.Stream('D_feed', Water=1)
    E_feed = D_feed.copy('E_feed')
    units = [A(None, 'A_feed'), B(None, 'B_feed'), C(None, 'C_feed'), D(None, D_feed), E(None, E_feed)]
    test_sys = bst.System('test_sys', units)
    test_sys.simulate()
    tea = create_tea(test_sys)
    table = tea.get_cashflow_table()
    C_FCI = table['Fixed capital investment [MM$]']
    # Test with lang factor = 3 (default for sugarcane biorefinery)
    cashflows_FCI = [6.12, 9.18, 0.  , 0.  , 0.  , 0.  , 0.  , 0.3 , 0.  , 0.  , 3.  ,
                     0.  , 3.3 , 0.  , 0.  , 0.  , 0.  , 3.3 , 3.  , 0.  , 0.  , 0.  ]
    assert_allclose(C_FCI, cashflows_FCI)
    # Cashflows include maintainance and others, so all entries are not zero
    cashflows = [
       -6120000., -9945000., -4321300., -4321300., -4321300., -4321300.,
       -4321300., -4621300., -4321300., -4321300., -7321300., -4321300.,
       -7621300., -4321300., -4321300., -4321300., -4321300., -7621300.,
       -7321300., -4321300., -4321300., -3556300.
    ]
    assert_allclose(tea.cashflow_array, cashflows)
    
    # Test with bare module costs
    tea.lang_factor = None
    table = tea.get_cashflow_table()
    C_FCI = table['Fixed capital investment [MM$]']
    cashflows_FCI = [
        6.12, 9.18, 0.  , 0.  , 0.  , 0.  , 0.  , 0.3 , 0.  , 0.  , 3.  ,
        0.  , 2.3 , 0.  , 0.  , 0.  , 0.  , 2.3 , 3.  , 0.  , 0.  , 0.  
    ]
    assert_allclose(C_FCI, cashflows_FCI)
    cashflows = [
       -6120000., -9945000., -4321300., -4321300., -4321300., -4321300.,
       -4321300., -4621300., -4321300., -4321300., -7321300., -4321300.,
       -6621300., -4321300., -4321300., -4321300., -4321300., -6621300.,
       -7321300., -4321300., -4321300., -3556300.
    ]
    assert_allclose(tea.cashflow_array, cashflows)

def test_skipping_unit_simulation_with_empty_inlet_streams():
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    s1 = bst.Stream('s1', Water=10, Ethanol=50)
    s2 = bst.Stream('s2', Water=10)
    M1 = bst.MixTank('M1', ins=(s1, s2))
    
    # Unit must simulate when skip flag is off
    bst.settings.skip_simulation_of_units_with_empty_inlets = False
    M1.simulate()
    assert M1.installed_cost > 0.
    assert M1.utility_cost > 0.
    assert not all([i.isempty() for i in M1.outs])
       
    # Unit must simulate when skip flag is on and inlets are not empty
    bst.settings.skip_simulation_of_units_with_empty_inlets = True
    M1.simulate()
    assert M1.installed_cost > 0.
    assert M1.utility_cost > 0.
    assert not all([i.isempty() for i in M1.outs])
    
    # Unit must have no costs nor outlet flows when skip flag is on and inlets are empty
    bst.settings.skip_simulation_of_units_with_empty_inlets = True
    s1.empty()
    s2.empty()
    M1.simulate()
    assert M1.installed_cost == 0.
    assert M1.utility_cost == 0.
    assert all([i.isempty() for i in M1.outs])
    
    # Unit should not simulate when the skip flag is on 
    # and inlets are empty 
    bst.settings.skip_simulation_of_units_with_empty_inlets = True
    @M1.add_specification
    def f(): raise RuntimeError('unit simulated')
    M1._design = f
    M1._cost = f
    M1.simulate()
    
    # Unit should simulate when the skip flag is on 
    # and inlets are not empty 
    bst.settings.skip_simulation_of_units_with_empty_inlets = False
    with pytest.raises(RuntimeError):
        M1.simulate()
        
if __name__ == '__main__':
    test_auxiliary_unit_owners()
    test_unit_inheritance_setup_method()
    test_process_specifications_linear()
    test_process_specifications_with_recycles()
    test_unit_connections()
    test_unit_graphics()
    test_cost_decorator()
    test_equipment_lifetimes()
    test_skipping_unit_simulation_with_empty_inlet_streams()
