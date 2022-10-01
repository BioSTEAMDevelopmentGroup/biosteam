# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import pytest
import biosteam as bst
from numpy.testing import assert_allclose

def test_unit_inheritance_setup_method():
    class NewUnit(bst.Unit):
        def _setup(self):
            pass
    
    bst.settings.set_thermo(['Water'], cache=True)
    U1 = NewUnit()
    U1.add_specification(lambda: None)
    with pytest.raises(bst.exceptions.UnitInheritanceError):
        U1.simulate()

def test_process_specifications():
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    with bst.System() as sys:
        T1 = bst.StorageTank(ins=bst.Stream(Water=1000))
        H1 = bst.HXutility(ins=T1-0, T=320)
        T2 = bst.StorageTank(ins=bst.Stream(Ethanol=1000))
        H2 = bst.HXutility(ins=T2-0, T=320)
        T3 = bst.StorageTank(ins=bst.Stream(Water=10))
        H3 = bst.HXutility(ins=T3-0, T=320)
        M1 = bst.Mixer(ins=[H1-0, H2-0, H3-0])
        H4 = bst.HXutility(ins=M1-0, T=350)
    
    # Specification impacting upstream units
    @M1.add_specification(run=True, impacted_units=[T1, T2])
    def adjust_flow_rate():
        water = T1.ins[0]
        ethanol = T2.ins[0]
        water.F_mass = ethanol.F_mass
    
    sys.simulate()
    assert (H4.outs[0].mol == sum([i.ins[0].mol for i in (T1, T2, T3)])).all()
    
    # Specification impacting units in parallel (neither upstream nor downstream units).
    # System simulation order must switch
    old_path_length = len(sys.unit_path)
    M1.specifications.pop() # Remove specification
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
    # Order of tank simulation is now reversed
    assert [last_tank, mid_tank, first_tank] == [i for i in sys.units if isinstance(i, bst.Tank)] 
    # 2 hidden connection source and 2 sinks
    assert len(sys.unit_path) == old_path_length
    assert (H4.outs[0].mol == sum([i.ins[0].mol for i in (T1, T2, T3)])).all()
    assert sys.path[0] is last_tank
    
    # Specification impacting downstream units (it doesn't matter).
    last_tank.specifications.pop() # Remove specification
    
    @T1.add_specification(run=True, impacted_units=[H3])
    def adjust_flow_rate():
        H3.T = 320
    
    sys.simulate()
    assert H3.outs[0].T == H3.T
    assert T1.specifications[0].units == []
    

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
    assert u.R301.neighborhood(100) == u.R301.neighborhood(1000) == all_units
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
    
if __name__ == '__main__':
    test_unit_inheritance_setup_method()
    test_process_specifications()
    test_unit_connections()
    test_unit_graphics()
    test_equipment_lifetimes()