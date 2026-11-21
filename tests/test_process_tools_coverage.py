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
import thermosteam as tmo
import numpy as np
from numpy.testing import assert_allclose

def create_units(name):
    bst.F.set_flowsheet(name)
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    feed = bst.Stream('feed', Water=100, Ethanol=50, price=0.1)
    H1 = bst.HXutility('H101', feed, T=360, V=0)
    P1 = bst.Pump('P201', H1-0, P=5e5)
    H2 = bst.HXutility('H202', P1-0, T=320)
    T1 = bst.StorageTank('T301', H2-0, tau=1)
    return feed, H1, P1, H2, T1

def test_segment():
    bst.F.set_flowsheet('segment_coverage')
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    feed = bst.Stream('feed', Water=100)
    M1 = bst.Mixer('M1', feed)
    H1 = bst.HXutility('H1', M1-0, T=320)
    P1 = bst.Pump('P1', H1-0)
    M2 = bst.Mixer('M2', P1-0)
    segment = bst.process_tools.Segment(M1-0, P1-0)
    assert segment.sink is H1 and segment.source is P1
    segment.disconnect(join_ends=True)
    assert M2.ins[0] is M1.outs[0]
    segment.insert(M1-0)
    assert H1.ins[0] is M1.outs[0] and M2.ins[0] is P1.outs[0]
    segment.disconnect()
    assert H1.ins[0].isempty() and not M2.ins[0].source

def test_utils():
    feed, H1, P1, H2, T1 = units = create_units('process_tools_utils')
    units = [H1, P1, H2, T1]
    sys = bst.System('sys', units)
    sys.simulate()
    groups = bst.process_tools.group_units_by_available_chemicals(units)
    assert sum([len(i) for i in groups.values()]) == 4
    assert bst.process_tools.is_storage_unit(T1)
    assert not bst.process_tools.is_storage_unit(P1)
    assert len(bst.process_tools.units_with_costs(units)) == 4
    areas = bst.process_tools.group_by_area(units)
    assert set(areas) == {100, 200, 300}
    assert len(bst.process_tools.group_by_lines(units)) == 3
    assert len(bst.process_tools.group_by_types(units)) == 3
    assert bst.process_tools.filter_by_types(units, bst.Pump) == [P1]
    assert bst.process_tools.filter_by_lines(units, ['Pump']) == [P1]
    heat_utilities = bst.process_tools.get_heat_utilities(units)
    assert len(bst.process_tools.filter_out_heat_utility_savings(heat_utilities)) >= 1
    assert bst.process_tools.get_utility_flow(heat_utilities, 'low_pressure_steam') > 0
    assert bst.process_tools.get_utility_duty(heat_utilities, 'low_pressure_steam') > 0
    assert bst.process_tools.get_heating_duty(heat_utilities) > 0
    assert bst.process_tools.get_cooling_duty(heat_utilities) > 0
    power_utilities = bst.process_tools.get_power_utilities(units)
    assert bst.process_tools.get_electricity_consumption(power_utilities) > 0
    assert bst.process_tools.get_electricity_production(power_utilities) == 0
    assert 0 < bst.process_tools.get_purchase_cost(units) < bst.process_tools.get_installed_cost(units)
    assert len(bst.process_tools.heat_exchanger_utilities_from_units(units)) >= 1
    assert bst.process_tools.volume_of_chemical_in_units(units, 'Water') > 0
    assert T1 in bst.process_tools.get_OSBL(units)
    bst.process_tools.set_construction_material_to_stainless_steel(units, '316')
    assert P1.material == 'Stainless steel'
    with pytest.raises(ValueError):
        bst.process_tools.set_construction_material_to_stainless_steel(units, '1')
    bst.process_tools.set_construction_material_to_carbon_steel(units)
    assert P1.material == 'Cast iron'
    assert H1.material == 'Carbon steel/carbon steel'
    bst.process_tools.default(flowsheet=False)

def test_connect_by_ID():
    bst.F.set_flowsheet('connect_by_ID')
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    M1 = bst.Mixer('M1', bst.Stream('feed', Water=1), outs='link')
    M2 = bst.Mixer('M2')
    M2.ins[0] = bst.Stream('link')
    bst.process_tools.connect_by_ID([M2, M1])
    assert M2.ins[0] is M1.outs[0]

def test_unit_group():
    feed, H1, P1, H2, T1 = create_units('unit_group_coverage')
    units = [H1, P1, H2, T1]
    sys = bst.System('sys', units)
    sys.simulate()
    group = bst.UnitGroup('group', units, filter_savings=False)
    assert len(list(group)) == 4
    assert group.to_system('group_sys').units == units
    upstream, downstream = group.split(P1-0, 'up', 'down')
    assert upstream.units == [H1, P1]
    group.register_utility_agent('low_pressure_steam')
    group.register_utility_agent('cooling_water', 'flow')
    with pytest.raises(ValueError):
        group.register_utility_agent('chilled_water', 'none')
    assert len(group.heat_utilities) >= 2
    assert group.get_utility_duty('low_pressure_steam') > 0
    assert group.get_utility_flow('cooling_water') > 0
    assert group.get_purchase_cost() > 0
    assert group.get_net_electricity_production() < 0
    assert group.get_inlet_flow('kg/hr') > 0
    assert group.get_inlet_flow('kg/hr', 'Water') > 0
    dct = group.to_dict(with_units=False)
    assert all([np.isfinite(i) for i in dct.values()])
    assert len(group.to_series()) == len(dct)
    assert 'UnitGroup' in repr(group)
    group.show()
    bst.UnitGroup('empty').show()
    assert len(bst.UnitGroup.filter_by_types('pumps', units, bst.Pump).units) == 1
    assert len(bst.UnitGroup.filter_by_lines('pumps', units, ['Pump']).units) == 1
    assert len(bst.UnitGroup.group_by_types(units)) == 3
    assert len(bst.UnitGroup.group_by_types(units, {'pumps': bst.Pump})) == 1
    assert len(bst.UnitGroup.group_by_lines(units)) == 3
    assert len(bst.UnitGroup.group_by_lines(units, {'pumps': ['Pump']})) == 1
    groups = bst.UnitGroup.group_by_area(units)
    for i in groups: i.autofill_metrics()
    df = bst.UnitGroup.df_from_groups(groups, fraction=True, scale_fractions_to_positive_values=False)
    assert df.shape[0] == 3
    def set_feed(x): feed.F_mass = x
    def simulate(x):
        set_feed(x)
        sys.simulate()
    df = bst.UnitGroup.df_from_groups_across_coordinate(groups, simulate, [1000, 2000], 'Feed')
    assert df.shape[0] == 3
    assert group.diagram(display=False) is not None
    grouped_units = bst.process_tools.unit_group.create_connected_grouped_units(
        [upstream, downstream], copy_streams=True
    )
    assert len(grouped_units) == 2

def test_system_mesh():
    bst.F.set_flowsheet('system_mesh_coverage')
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)

    @bst.SystemFactory(ins=['feed'], outs=['heated'])
    def create_heating_system(ins, outs, T=350):
        H1 = bst.HXutility(ins=ins[0], outs=outs[0], T=T)

    @bst.SystemFactory(ins=['heated'], outs=[dict(ID='product')])
    def create_pumping_system(ins, outs):
        P1 = bst.Pump(ins=ins[0], outs=outs[0], P=3e5)

    called = []
    def add_note(): called.append(True)

    mesh = bst.SystemMesh()
    mesh.add('heating', create_heating_system, 100, T=340)
    mesh.add('pumping', create_pumping_system)
    mesh.add('note', add_note)
    mesh.show()
    copy = mesh.copy()
    copy.remove('note')
    copy.disconnect('heated')
    with pytest.raises(ValueError): copy.disconnect('none')
    with pytest.raises(ValueError): copy.connect('none', 'heated')
    with pytest.raises(ValueError): copy.connect('heated', 'none')
    copy.connect('heated', 'heated')
    mesh.add('pumping', create_pumping_system) # Replace object
    sys = mesh('mesh_sys', feed=bst.Stream('feed_stream', Water=10, Ethanol=1))
    assert called
    with pytest.raises(ValueError):
        mesh('mesh_sys2', none=bst.Stream())
    sys.simulate()
    print(sys.path, [(i.ID, i.T, i.P) for i in sys.products])
    product, = sys.products
    assert_allclose(product.T, 340, rtol=1e-2)
    assert product.P > 2.5e5
    with pytest.raises(TypeError):
        bst.process_tools.system_mesh.meshable(1)

def test_system_mesh_with_units_and_systems():
    bst.F.set_flowsheet('system_mesh_units_coverage')
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    feed = bst.Stream('mesh_feed', Water=10)
    H1 = bst.HXutility('H1', feed, 'hot', T=340)
    P1 = bst.Pump('P1', 'hot', 'pumped', P=3e5)
    M1 = bst.Mixer('M1', 'pumped')
    other = bst.HXutility('H2', bst.Stream('other_feed', Water=1), T=310)
    subsys = bst.System('subsys', [other])
    mesh = bst.SystemMesh()
    mesh.add('heater', H1)
    mesh.add('pump', P1, 200)
    mesh.add('mixer', [M1])
    mesh.add('subsys', subsys, 300)
    mesh.show()
    sys = mesh('mesh_sys')
    assert P1.ID == 'P201'
    sys.simulate()
    assert P1.ins[0] is H1.outs[0]

def test_system_factory():
    bst.F.set_flowsheet('system_factory_coverage')
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    kwargs = bst.process_tools.stream_kwargs(
        'feed', flow=None, phase='l', T=300, P=101325, units='kg/hr',
        price=0.1, total_flow=10, thermo=None, characterization_factors={},
        Water=1,
    )
    assert kwargs['T'] == 300 and kwargs['phase'] == 'l'
    with pytest.raises(ValueError):
        bst.SystemFactory(lambda a, b: None)
    with pytest.raises(ValueError):
        bst.SystemFactory(lambda ins, outs, area: None)

    def default_outlet():
        return bst.Stream('generated', Water=1)

    @bst.SystemFactory(
        ID='heat_sys', ins=[dict(ID='feed', Water=10, NotAChemical=1)],
        outs=['product', default_outlet],
        fthermo=lambda: bst.Thermo(['Water', 'Ethanol'], cache=True),
    )
    def create_heat_system(ins, outs):
        H1 = bst.HXutility('H1', ins[0], outs[0], T=340)
        H2 = bst.HXutility('H2', bst.Stream(Water=1), outs[1], T=340)

    bst.SystemFactory(lambda ins, outs: None, ins=['feed'], outs=['product']).show()
    bst.SystemFactory(lambda ins, outs: None, ID=None).show()
    assert 'create_heat_system' in repr(create_heat_system)
    sys = create_heat_system()
    sys.simulate()
    assert_allclose(sys.get_outlet('product').T, 340)
    sys, udct = create_heat_system('heat_sys_2', ins=['new_feed'], outs=[None], udct=True, area=200)
    assert 'H1' in udct and len(udct['Heat exchanger']) == 2
    with pytest.raises(ValueError):
        create_heat_system(ins=[bst.Stream(), bst.Stream()])
    with pytest.raises(TypeError):
        create_heat_system(ins=[1])
    mock = create_heat_system(mockup=True, autorename=None)
    assert isinstance(mock, bst.MockSystem)
    with pytest.raises(ValueError):
        create_heat_system(box=True, mockup=True)

class SimpleReactor(bst.Unit):
    _N_ins = 1
    _N_outs = 2
    
    def _init(self, tau=1.):
        self.tau = tau
        self.reaction = tmo.Reaction('2Methanol -> Ethanol + Water', 'Methanol', 0.9)
        
    def _run(self):
        vent, effluent = self.outs
        effluent.copy_like(self.ins[0])
        self.reaction(effluent)
        vent.empty()
        
    def _design(self):
        self.design_results['Volume'] = self.tau * self.outs[1].F_vol

def test_reactor_specification():
    bst.F.set_flowsheet('reactor_specification_coverage')
    bst.settings.set_thermo(['Water', 'Ethanol', 'Methanol'], cache=True)
    feed = bst.Stream('feed', Water=100, Methanol=10)
    R1 = SimpleReactor('R1', feed)
    sys = bst.System('sys', [R1])
    spec = bst.process_tools.ReactorSpecification(
        R1, 'reaction', substrates=('Methanol',), products=('Ethanol',),
        yield_=0.8, titer=20., productivity=1., production=50.,
    )
    assert spec.feed is feed
    assert spec.vent is R1.outs[0]
    assert spec.effluent is R1.outs[1]
    assert 0 < spec.products_over_substrates < 1
    spec.load_specifications()
    R1.simulate()
    effluent = R1.outs[1]
    assert_allclose(effluent.imass['Ethanol'] / effluent.F_vol, 20., rtol=1e-2)
    assert_allclose(effluent.imass['Ethanol'], 50., rtol=1e-2)
    assert_allclose(R1.tau, 20.)
    def volume(): return R1.design_results['Volume']
    data = spec.evaluate_across_productivity([volume], [1., 2.])
    assert data.shape == (1, 2)
    assert data[0, 0] > data[0, 1] > 0
    data = spec.evaluate_across_TRY(sys, np.array([15., 20.]), np.array([0.7, 0.8]), [volume], [1., 2.])
    assert data.shape == (2, 1, 2)
    assert np.isfinite(data).all() and (data > 0).all()
    
def test_process_model():
    
    class HeaterModel(bst.ProcessModel):
        
        class Scenario:
            T: float = (320., '# Heater temperature [K]')
            label: str = 'baseline'
            simple: bool = True
        
        def create_thermo(self):
            return bst.Thermo(['Water', 'Ethanol'], cache=True)
        
        def create_system(self):
            feed = bst.Stream('feed', Water=10, Ethanol=1)
            H1 = bst.HXutility('H1', feed, T=self.scenario.T)
        
        def create_model(self):
            model = bst.Model(self.system)
            
            @model.parameter(baseline=self.scenario.T, bounds=(300, 340))
            def set_T(T): self.H1.T = T
            
            @model.indicator
            def outlet_T(): return self.H1.outs[0].T
            
            return model
    
    Scenario = HeaterModel.Scenario
    assert Scenario.metadata == {'T': 'Heater temperature [K]'}
    scenario = Scenario()
    scenario.show()
    other = scenario.copy(T=330.)
    assert other.T == 330.
    assert dict(other)['T'] == 330.
    assert (scenario - other).right is other
    pm = HeaterModel()
    assert_allclose(pm.H1.outs[0].T, 320.)
    assert HeaterModel() is pm # Cached
    pm2 = HeaterModel(T=330., simulate=False, save=False)
    assert pm2 is not pm
    assert 'HeaterModel' in repr(pm2)
    pm3 = HeaterModel(scenario=other, load=False, save=False)
    assert_allclose(pm3.H1.outs[0].T, 330., rtol=1e-3)
    sample, results = pm3.baseline()
    assert_allclose(results.iloc[0], 330., rtol=1e-3)
    assert len(pm3.parameters) == 1 and len(pm3.indicators) == 1
    with pytest.raises(NotImplementedError):
        HeaterModel(scenario='not a scenario')
    
    class DefaultModel(HeaterModel):
        Scenario = HeaterModel.Scenario
        
        def default_scenario(self=None, **kwargs):
            return HeaterModel.Scenario(T=310.)
        
        @classmethod
        def as_scenario(cls, scenario):
            return HeaterModel.Scenario(T=float(scenario))
    
    pm4 = DefaultModel(simulate=False)
    assert pm4.scenario.T == 310.
    pm5 = DefaultModel(325, simulate=False, label='other')
    assert pm5.scenario.T == 325. and pm5.scenario.label == 'other'
    
    class MissingThermo(bst.ProcessModel): pass
    
    with pytest.raises(NotImplementedError):
        MissingThermo()
    
    class MissingSystem(bst.ProcessModel):
        def create_thermo(self): return bst.Thermo(['Water'], cache=True)
    
    with pytest.raises(NotImplementedError):
        MissingSystem()
    
    class MissingModel(MissingSystem):
        def create_system(self): bst.Mixer('M1', bst.Stream(Water=1))
    
    with pytest.raises(NotImplementedError):
        MissingModel()
    
    class NoneModel(MissingModel):
        def create_model(self): return None
    
    with pytest.raises(RuntimeError):
        NoneModel()

def test_system_mesh_duplicate_names():
    bst.F.set_flowsheet('system_mesh_duplicates')
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    
    def water(): return bst.Stream(Water=1)
    
    @bst.SystemFactory(ins=['feed', water], outs=['product'])
    def create_mixing_system(ins, outs):
        bst.Mixer(ins=ins, outs=outs[0])
    
    mesh = bst.SystemMesh()
    mesh.add('mixing_a', create_mixing_system, autoconnect=False)
    mesh.add('mixing_b', create_mixing_system, autoconnect=False)
    mesh.add('mock', [bst.Mixer('M1'), bst.Mixer('M2'), bst.Mixer('M3'), bst.Mixer('M4')])
    mesh.show()
    assert len(mesh._inlets) == 4 and len(mesh._outlets) == 2
    outlet = [i for i in mesh._outlets if mesh._outlets[i][0] == 'mixing_a'][0]
    inlet = [i for i in mesh._inlets if mesh._inlets[i] == ('mixing_b', 0)][0]
    mesh.connect(outlet, inlet)
    sys = mesh('mesh_sys')
    sys.simulate()
    assert_allclose(sum([i.F_mol for i in sys.products]), sum([i.F_mol for i in sys.feeds]), rtol=1e-2)

if __name__ == '__main__':
    test_segment()
    test_utils()
    test_connect_by_ID()
    test_unit_group()
    test_system_mesh()
    test_system_mesh_with_units_and_systems()
    test_system_factory()
    test_reactor_specification()
    test_process_model()
    test_system_mesh_duplicate_names()
