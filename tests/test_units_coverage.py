# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import matplotlib
matplotlib.use('Agg')
import os
import tempfile
import pytest
import biosteam as bst
import numpy as np
import pandas as pd
from numpy.testing import assert_allclose

def assert_mass_balance(unit, rtol=1e-3):
    F_in = sum([i.F_mass for i in unit.ins])
    F_out = sum([i.F_mass for i in unit.outs])
    assert_allclose(F_in, F_out, rtol=rtol)

def assert_positive_costs(unit):
    assert unit.purchase_cost > 0
    assert np.isfinite(unit.purchase_cost)
    assert unit.installed_cost >= unit.purchase_cost

def test_mixer_phenomena_equations():
    bst.F.set_flowsheet('mixer_phenomena')
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    liquid = bst.Stream(Water=10, T=300)
    vapor = bst.Stream(Ethanol=10, T=400, phase='g')
    M = bst.Mixer(ins=[liquid, vapor], rigorous=True)
    M.simulate()
    assert_mass_balance(M)
    assert M.outs[0].phases == ('g', 'l')
    assert 0 < M.B < np.inf
    assert M.T_node is M.T_node
    assert M.B_node is M.B_node
    assert M.get_E_node(liquid) is M.B_node
    assert len(M._create_material_balance_equations(False)) == 2
    assert len(M._create_energy_balance_equations()) == 1
    assert len(M._create_bulk_balance_equations()) == 2
    M._update_variable('T', 350)
    assert M.outs[0].T == 350
    with pytest.raises(ValueError):
        M._update_variable('B', 1)
    M._update_nonlinearities()
    # Single phase outlet
    M = bst.Mixer(ins=[liquid.copy(), bst.Stream(Ethanol=1)])
    M.simulate()
    assert M.B == 0
    assert M.get_E_node(liquid) is M.T_node
    assert len(M._create_material_balance_equations(False)) == 1
    assert len(M._create_energy_balance_equations()) == 1
    assert len(M._create_bulk_balance_equations()) == 1
    assert M._update_energy_coefficient(liquid, {}) > 0
    # Vapor outlet
    M = bst.Mixer(ins=[bst.Stream(Water=1, T=400, phase='g')], rigorous=True)
    M.simulate()
    assert M.B == np.inf

def test_steam_mixer_solids_loading():
    bst.F.set_flowsheet('steam_mixer')
    bst.settings.set_thermo(['Water', 'Glucose'], cache=True)
    feed = bst.Stream('feed', Water=10, Glucose=10)
    M1 = bst.SteamMixer(ins=[feed, 'steam', 'process_water'], T=431.15, P=557287.5,
                        solids_loading=0.3, solid_IDs=['Glucose'])
    M1.simulate()
    out = M1.outs[0]
    assert_allclose(out.T, 431.15, rtol=1e-3)
    assert 0.2 < out.imass['Glucose'] / out.F_mass < 0.3
    M1.solids_loading_includes_steam = True
    M1.simulate()
    assert_allclose(out.imass['Glucose'] / out.F_mass, 0.3, rtol=1e-2)
    M1.solid_IDs = None
    M1.simulate()
    assert_allclose(out.imass['Glucose'] / out.F_mass, 0.3, rtol=1e-2)

def test_reversed_splitter_and_separator():
    bst.F.set_flowsheet('reversed_splitter')
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    feed = bst.Stream(T=320)
    S = bst.ReversedSplitter(ins=feed, outs=[bst.Stream(Water=1), bst.Stream(Ethanol=2)])
    S.simulate()
    assert_allclose(feed.mol, [1, 2])
    assert S.outs[0].T == 320
    inlet = bst.Stream(T=330)
    outlets = [bst.Stream(Water=3), bst.Stream(Ethanol=1)]
    bst.units.splitting.reversed_split(inlet, outlets)
    assert_allclose(inlet.mol, [3, 1])
    assert outlets[1].T == 330
    # Separator
    feed = bst.Stream(Water=10, Ethanol=10, T=300)
    S = bst.units.splitting.Separator(ins=feed, split={'Ethanol': 0.9}, T=350, P=2e5, phases=('g', 'l'))
    S.simulate()
    top, bottom = S.outs
    assert_allclose(top.imol['Ethanol'], 9)
    assert top.T == bottom.T == 350
    assert top.P == 2e5
    assert top.phase == 'g'
    assert S.isplit['Ethanol'] == 0.9
    S.split = {'Water': 0.1}
    assert S.isplit['Ethanol'] == 0 and S.isplit['Water'] == 0.1
    assert S.equation_node_names == ('overall_material_balance_node', 'separation_material_balance_node')
    assert S._get_energy_departure_coefficient(feed) is None
    assert S._create_energy_departure_equations() == []
    S.T = None
    assert len(S.equation_node_names) == 3
    assert S.T_node is S.E_node
    assert S.get_E_node(feed) is S.T_node
    S.simulate()
    assert S._get_energy_departure_coefficient(feed)[1] < 0
    T = top.T
    S._update_energy_variable(1.)
    assert_allclose(top.T, T + 1)

def test_unit_factories():
    from biosteam.units import factories
    bst.F.set_flowsheet('factories')
    columns = pd.MultiIndex.from_tuples([('Mixer-Blender', 'Blender'), ('Static tank', 'Tank')])
    index = ['Basis', 'Units', 'Size', 'Lower bound', 'Upper bound', 'CEPCI',
             'Cost (USD)', 'Exponent', 'Electricity (kW)', 'Installation factor',
             'Number', 'Lifetime (yr)']
    data = [['Flow rate', 'Flow rate'], ['kg/hr', 'kg/hr'], [1000, 1000],
            [False, 10], [False, 1e5], [500, 500], [1e4, 2e4], [0.6, 0.7],
            [1, False], [2, 1.5], [False, False], [False, 10]]
    df = pd.DataFrame(data, index=index, columns=columns, dtype=object)
    txt = factories.df2txt(df)
    assert 'class Blender(Mixer)' in txt
    assert 'class StaticTank(Unit)' in txt
    assert 'kW=1' in txt
    assert 'lifetime=10' in txt
    with tempfile.TemporaryDirectory() as folder:
        file = os.path.join(folder, 'units.xlsx')
        df.to_excel(file)
        assert 'class Blender(Mixer)' in factories.xl2txt(file)
        bad = df.copy()
        bad.index = [*index[:-1], 'Wrong']
        bad.to_excel(file)
        with pytest.raises(RuntimeError):
            factories.xl2txt(file)
    df.columns = pd.MultiIndex.from_tuples([('NotAUnit-Blender', 'Blender'), ('Static tank', 'Tank')])
    with pytest.raises(ValueError):
        factories.df2txt(df)

def create_fcc_thermo():
    catalyst = bst.Chemical('Catalyst', search_db=False, default=True, phase='s', MW=60, rho=1200, Cp=1.0)
    octene = bst.Chemical('Octene', search_ID='1-octene')
    bst.settings.set_thermo(['Water', 'O2', 'N2', 'CO2', 'Hexadecane', 'Octane', octene, catalyst])

def test_fluidized_catalytic_cracking():
    bst.F.set_flowsheet('fcc')
    create_fcc_thermo()
    feed = bst.Stream(Hexadecane=10, T=300)
    rxn = bst.Reaction('Hexadecane -> Octane + Octene', 'Hexadecane', 0.8)
    U = bst.FluidizedCatalyticCracking(ins=[feed, 'catalyst', 'air', 'steam'], reaction=rxn)
    U.simulate()
    assert_mass_balance(U)
    assert U.feed is feed
    assert U.makeup_catalyst is U.ins[1]
    assert U.air.F_mass > 0
    assert U.steam.F_mass > 0
    assert U.product.T > feed.T
    assert U.flue_gas.T > feed.T
    assert U.spent_catalyst.imass['Catalyst'] > 0
    assert_positive_costs(U)
    assert U.power_utility.rate > 0
    # Bulk reactant
    rxn = bst.Reaction('Hexadecane -> Octane + Octene', 'Hexadecane', 1)
    U = bst.FluidizedCatalyticCracking(
        ins=[feed, 'catalyst', 'air', 'steam'], reaction=rxn,
        bulk_reactant=('Hexadecane', ['Hexadecane'])
    )
    U.simulate()
    assert_mass_balance(U)
    # NOTE: Simulating the same FCC unit twice raises UnitInheritanceError
    # (FluidizedCatalyticCracking._setup does not call super()._setup()),
    # so a new unit is created.
    rxn = bst.Reaction('Hexadecane -> Octane + Octene', 'Hexadecane', 0.5)
    U = bst.FluidizedCatalyticCracking(
        ins=[feed, 'catalyst', 'air', 'steam'], reaction=rxn,
        bulk_reactant=('Hexadecane', ['Hexadecane'])
    )
    with pytest.raises(RuntimeError):
        U.simulate()

def test_auxiliary_pressure_vessel():
    bst.settings.set_thermo(['Water'], cache=True)
    V = bst.AuxiliaryPressureVessel(100, 5, 20, orientation='Horizontal', material='Stainless steel 316')
    assert V.design_results['Weight'] > 0
    assert V.design_results['Wall thickness'] > 0
    assert sum(V.baseline_purchase_costs.values()) > 0
    V = bst.AuxiliaryPressureVessel(2e5, 2, 10, pressure_units='Pa', length_units='m')
    assert V.design_results['Weight'] > 0
    assert sum(V.baseline_purchase_costs.values()) > 0
    with pytest.warns(bst.exceptions.DesignWarning):
        bst.AuxiliaryPressureVessel(10, 5, 20)
    with pytest.raises(ValueError):
        bst.AuxiliaryPressureVessel(100, 5, 20, orientation='Diagonal')

def test_spray_dryer():
    bst.F.set_flowsheet('spray_dryer')
    bst.settings.set_thermo(['Water', 'Glucose'], cache=True)
    feed = bst.Stream(Water=1000, Glucose=100, units='kg/hr')
    D = bst.SprayDryer(ins=feed, moisture_content=0.1)
    D.simulate()
    assert_mass_balance(D)
    water, solids = D.outs
    assert_allclose(solids.imass['Water'] / solids.F_mass, 0.1, rtol=1e-3)
    assert water.phase == 'g'
    assert_positive_costs(D)

def test_drum_dryer_options():
    bst.F.set_flowsheet('drum_dryer')
    bst.settings.set_thermo(['Water', 'Glucose', 'O2', 'N2', 'CO2', 'CH4'], cache=True)
    feed = bst.Stream(Water=1000, Glucose=200, units='kg/hr')
    D = bst.DrumDryer(ins=feed, split={'Glucose': 0.0}, utility_agent='Steam')
    D.simulate()
    assert D.isplit['Glucose'] == 0
    assert D.natural_gas.isempty()
    assert D.outs[1].imass['Water'] > 0
    assert_positive_costs(D)
    assert_allclose(D._convert_air_mol_to_mass(1., D.gas_composition), 1 / (0.79/28.0134 + 0.21/31.9988), rtol=1e-2)
    D = bst.DrumDryer(ins=feed, split={'Glucose': 0.0}, P=101325, T=390)
    with pytest.warns(bst.exceptions.DesignWarning):
        D.simulate()
    assert D.natural_gas.F_mass > 0

def test_thermal_oxidizer():
    bst.F.set_flowsheet('thermal_oxidizer')
    bst.settings.set_thermo(['Water', 'Ethanol', 'O2', 'N2', 'CO2', 'CH4'], cache=True)
    feed = bst.Stream(Water=100, Ethanol=1, N2=100, phase='g', T=350)
    U = bst.ThermalOxidizer(ins=[feed, 'air', 'natural_gas'])
    U.simulate()
    assert_mass_balance(U)
    emissions, = U.outs
    assert emissions.T > feed.T
    assert emissions.imol['Ethanol'] == 0
    assert U.natural_gas.F_mass > 0
    assert U.design_results['Number of vessels'] >= 1
    assert_positive_costs(U)

def test_terminal():
    bst.F.set_flowsheet('terminal')
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    s = bst.Stream(Water=1)
    T = bst.units.Terminal(upstream=s)
    assert T.upstream is s
    assert T.downstream is None
    with pytest.raises(AttributeError):
        T.downstream = bst.Stream()
    T.upstream = s2 = bst.Stream(Ethanol=1)
    assert T.upstream is s2
    T.simulate()
    assert T._get_tooltip_string().startswith('Terminal')
    assert T._inlet_utility_indices == {}
    T = bst.units.Terminal(downstream=bst.Stream())
    assert T.upstream is None
    with pytest.raises(AttributeError):
        T.upstream = bst.Stream()
    T.downstream = s3 = bst.Stream()
    assert T.downstream is s3
    with pytest.raises(ValueError):
        bst.units.Terminal()
    with pytest.raises(ValueError):
        bst.units.Terminal(upstream=bst.Stream(), downstream=bst.Stream())

def test_junction_and_scaler():
    bst.F.set_flowsheet('junction')
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    s1 = bst.Stream(Water=2)
    s2 = bst.Stream()
    J = bst.units.Junction(upstream=s1, downstream=s2)
    assert J.upstream is s1 and J.downstream is s2
    J.upstream = s3 = bst.Stream(Ethanol=1)
    J.downstream = s4 = bst.Stream()
    J.simulate()
    assert_allclose(s4.mol, s3.mol)
    assert J._get_tooltip_string().startswith('Junction')
    assert J._inlet_utility_indices == {}
    S = bst.units.Scaler(upstream=s3, downstream=bst.Stream(), scale=2)
    assert S.upstream is s3
    S.upstream = s1
    S.downstream = s5 = bst.Stream()
    assert S.downstream is s5
    S.simulate()
    assert_allclose(s5.mol, 2 * s1.mol)
    assert S._get_tooltip_string() == 'Scaler'

def test_enzyme_treatment():
    bst.F.set_flowsheet('enzyme_treatment')
    bst.settings.set_thermo(['Water', 'Glucose'], cache=True)
    feed = bst.Stream(Water=1000, Glucose=10, T=300)
    U = bst.units.EnzymeTreatment(ins=feed, T=320)
    U.simulate()
    assert_mass_balance(U)
    assert U.outs[0].T == 320
    assert U.heat_exchanger.heat_utilities[0].duty > 0
    assert_positive_costs(U)
    U = bst.units.EnzymeTreatment(ins=feed)
    U.simulate()
    assert_allclose(U.outs[0].T, 300)

def test_vent_scrubber():
    bst.F.set_flowsheet('vent_scrubber')
    bst.settings.set_thermo(['Water', 'Ethanol', 'CO2'], cache=True)
    vent = bst.Stream(CO2=100, Ethanol=1, Water=1, phase='g')
    U = bst.units.VentScrubber(ins=['water', vent], titer=('Ethanol', 10))
    U.simulate()
    assert_mass_balance(U)
    assert U.ins[0].F_mass > 0
    assert_positive_costs(U)
    U = bst.units.VentScrubber(ins=['water', vent], titer=10)
    U.simulate()
    assert_mass_balance(U)
    assert U.ins[0].F_mass > 0

def test_screw_feeder_and_pellet_mills():
    bst.F.set_flowsheet('screw_feeder')
    bst.settings.set_thermo(['Water', 'Glucose'], cache=True)
    feed = bst.Stream(Glucose=1000, units='kg/hr')
    U = bst.units.ScrewFeeder(ins=feed)
    U.simulate()
    assert U.power_utility.rate > 0
    assert_positive_costs(U)
    U = bst.units.PelletMill(ins=feed)
    U.simulate()
    assert_allclose(U.power_utility.rate, 0.01 * 1000, rtol=1e-2)
    U = bst.units.BagassePelletMill(ins=feed)
    U.simulate()
    assert U.power_utility.rate > 0
    U.ins[0].imass['Glucose'] = 4000
    U.simulate()
    assert U.power_utility.rate > 114.5

def test_single_phase_reactor_options():
    bst.F.set_flowsheet('single_phase_reactor')
    bst.settings.set_thermo(['Water', 'Ethanol', 'Glucose', 'CO2'], cache=True)
    feed = bst.Stream(Water=1000, Glucose=10, T=300, P=2e5)
    rxn = bst.Reaction('Glucose -> 2Ethanol + 2CO2', 'Glucose', 0.5)
    R = bst.SinglePhaseReactor(ins=feed, reactions=rxn, tau=1, P=101325, adiabatic=True)
    R.simulate()
    assert_mass_balance(R)
    assert R.outs[0].T > feed.T # Exothermic
    # NOTE: Outlet pressure specification is overwritten by `mix_from`
    R = bst.SinglePhaseReactor(ins=feed, reactions=rxn, tau=1)
    R.simulate()
    assert R.outs[0].P == feed.P
    assert_allclose(R.outs[0].T, 300, rtol=1e-3)

def test_isenthalpic_valve_low_pressure():
    bst.F.set_flowsheet('valve')
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    feed = bst.Stream(Water=10, P=101325)
    V = bst.IsenthalpicValve(ins=feed, P=2e5)
    with pytest.warns(RuntimeWarning):
        V.simulate()
    assert V.outs[0].P == feed.P

def test_pump_types():
    bst.F.set_flowsheet('pump_types')
    bst.settings.set_thermo(['Water', 'Glycerol'], cache=True)
    feed = bst.Stream(Water=100, T=300)
    P = bst.Pump(ins=feed, P=5e5, ignore_NPSH=False)
    P.simulate()
    assert P.design_results['NPSH'] > 0
    assert_positive_costs(P)
    for kind in ('Gear', 'CentrifugalDouble'):
        P = bst.Pump(ins=feed, P=5e5, pump_type=kind)
        P.simulate()
        assert P.design_results['Type'] == kind
        assert_positive_costs(P)
    viscous = bst.Stream(Glycerol=100, T=300)
    P = bst.Pump(ins=viscous, P=5e5)
    P.simulate()
    assert P.design_results['Type'] != 'Centrifugal'
    assert_positive_costs(P)

def test_vacuum_system_deduced_from_unit():
    bst.F.set_flowsheet('vacuum_system')
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    with pytest.raises(ValueError):
        bst.VacuumSystem(F_mass=1, F_vol=1)
    vacuum = bst.VacuumSystem(F_mass=100, F_vol=1000, P_suction=5000, vessel_volume=10)
    assert sum(vacuum.baseline_purchase_costs.values()) > 0

    class VacuumVessel(bst.Unit):
        _N_outs = 1
        _units = {'Volume': 'm3'}
        P = 5000
        def _run(self):
            out = self.outs[0]
            out.copy_like(self.ins[0])
            out.P = self.P
            out.phase = 'g'
        def _design(self):
            self.design_results['Volume'] = 10
            self.vacuum_system = bst.VacuumSystem(self)

    V = VacuumVessel(ins=bst.Stream(Water=10, T=320))
    V.simulate()
    assert sum(V.vacuum_system.baseline_purchase_costs.values()) > 0
    V._design = lambda: setattr(V, 'vacuum_system', bst.VacuumSystem(V))
    V.design_results.clear()
    with pytest.raises(ValueError):
        V._design()

def create_adsorption_thermo():
    bst.settings.set_thermo([
       'Water', 
       bst.Chemical('Adsorbate', search_db=False, default=True, phase='l'),
       bst.Chemical('ActivatedCarbon', search_db=False, default=True, phase='s')
    ])

def test_adsorption_column_heuristic_design():
    bst.F.set_flowsheet('adsorption_heuristic')
    create_adsorption_thermo()
    feed = bst.Stream(phase='l', T=298, P=1.01e+06, Water=1000, Adsorbate=0.001, units='kg/hr')
    A = bst.AdsorptionColumn(
        ins=feed, cycle_time=1000, superficial_velocity=9.2, 
        adsorbate='Adsorbate', N_columns=2,
    )
    A.simulate()
    assert A.LUB == A.LUB_forced > 0
    assert_allclose(A.column_length, A.LES + A.LUB)
    assert A.outs[0].imol['Adsorbate'] == 0
    assert_positive_costs(A)
    assert A.parallel['self'] == 2
    # No adsorbate
    feed.imol['Adsorbate'] = 0
    A.simulate()
    assert A.design_results['Weight'] == 0
    assert A.purchase_cost == 0
    # Invalid arguments
    with pytest.raises(ValueError):
        bst.AdsorptionColumn(ins=feed, adsorbate='Adsorbate', isotherm_args=(1, 1))
    with pytest.raises(ValueError):
        bst.AdsorptionColumn(ins=feed, adsorbate='Adsorbate', isotherm_model='Unknown')
    with pytest.raises(ValueError):
        bst.AdsorptionColumn(ins=feed, adsorbate='Adsorbate', regeneration_isotherm_args=(1, 1))
    with pytest.raises(ValueError):
        bst.AdsorptionColumn(ins=feed, adsorbate='Adsorbate', k=0.3, LUB_forced=1)

def test_adsorption_isotherm_fits():
    import matplotlib.pyplot as plt
    Ce = np.array([0.06013, 0.31769, 0.63124, 2.03102, 4.16988, 5.00974]) / 1000
    qe = np.array([1.24, 1.87, 2.85, 4.48, 5.5, 6.53])
    t = np.array([0.0, 31.0, 60.0, 90.0, 120.0, 150.0]) / 60
    C = np.array([7.0, 2.45, 1.03, 0.43, 0.26, 0.13])
    AC = bst.AdsorptionColumn
    dct = AC.fit_Freundlich_isotherm(Ce, qe)
    assert dct['K'] > 0 and dct['n'] > 0 and dct['R2'] > 0.9
    dct = AC.fit_Langmuir_isotherm(Ce, qe)
    assert dct['K'] > 0 and dct['qmax'] > 0 and dct['R2'] > 0.9
    for method in ('Langmuir', 'Freundlich'):
        AC.plot_isotherm_fit(Ce, qe, method)
        plt.close('all')
    with pytest.raises(ValueError):
        AC.plot_isotherm_fit(Ce, qe, 'Other')
    fig, axes, params = AC.plot_isotherm_and_mass_transfer_coefficient_fit(
        Ce, qe, t, C, 0.075, 1, C_units='g/L', q_units='g/kg', t_units='hr',
    )
    assert params['k'] > 0
    fig, axes, params = AC.plot_isotherm_and_mass_transfer_coefficient_fit(
        Ce, qe, t, C, 0.075, 1, method='Langmuir',
    )
    assert params['isotherm_model'] == 'Langmuir'
    plt.close('all')

def test_cost_and_design_decorators():
    from biosteam.units.decorators import cost, copy_algorithm, add_cost, CostItem, design
    bst.F.set_flowsheet('decorators')
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    with pytest.raises(ValueError):
        CostItem('Flow rate', 'kg/hr', 1, None, None, 500, 1, 0.6, None, None, lambda x: x, None, None)
    with pytest.raises(ValueError):
        CostItem('Flow rate', 'kg/hr', 1, None, None, 500, 1, 0.6, None, 2, None, None, None)
    item = CostItem('Flow rate', 'kg/hr', 1, 1, 10, 500, 1, 0.6, 1, 'N', None, None, True)
    new = item.copy()
    assert (new.basis, new.units, new.N, new.lb) == ('Flow rate', 'kg/hr', 'N', 1)
    assert 'Flow rate' in repr(item)
    item.show()
    
    @cost('Flow rate', 'Pump', CE=500, cost=1e3, n=0.6, ub=10, N='N_pumps')
    @cost('Flow rate', 'Vessel', units='kg/hr', CE=500, cost=1e4, n=0.6, lb=1e4, kW=1, magnitude=True)
    @cost('Annual flow rate', 'Storage', units='kg/hr', CE=500, cost=1e3, n=0.6)
    class Dummy(bst.Unit):
        _N_ins = 2
        _N_outs = 1
        N_pumps = 2
        def _run(self):
            self.outs[0].mix_from(self.ins)
    
    U = Dummy(ins=[bst.Stream(Water=10), bst.Stream(Ethanol=10)])
    sys = bst.System.from_units(units=[U], operating_hours=8000)
    sys.simulate()
    assert_allclose(U.design_results['Flow rate'], U.outs[0].F_mass)
    assert U.design_results['Annual flow rate'] > 0
    assert U.parallel['Pump'] >= 2
    assert U.power_utility.rate > 0
    assert U.baseline_purchase_costs['Vessel'] > 0
    with pytest.raises(RuntimeError):
        add_cost(Dummy, 'Other', 'Flow rate', 'm3/hr', 1, None, None, 500, 1, 1, None, 1, None, None, None, None, None, None)
    with pytest.raises(RuntimeError):
        add_cost(Dummy, 'Other', 'Not a basis', None, 1, None, None, 500, 1, 1, None, 1, None, None, None, None, None, None)
    with pytest.raises(ValueError):
        add_cost(Dummy, None, 'Flow rate', None, 1, None, None, 500, 1, 1, None, 1, None, None, None, None, None, None)
    with pytest.raises(ValueError):
        add_cost(Dummy, 'Pump', 'Flow rate', None, 1, None, None, 500, 1, 1, None, 1, None, None, None, None, None, None)
    # Flow rate basis by outlets
    @cost('Flow rate', units='kg/hr', CE=500, cost=1e3, n=0.6)
    class Splitter(bst.Unit):
        _N_ins = 1
        _N_outs = 2
        def _run(self): pass
    @cost('Flow rate', units='kg/hr', CE=500, cost=1e3, n=0.6)
    class Junction(bst.Unit):
        _N_ins = 2
        _N_outs = 3
        def _run(self): pass
    @cost('Flow rate', units='kg/hr', CE=500, cost=1e3, n=0.6)
    class Merger(bst.Unit):
        _N_ins = 3
        _N_outs = 2
        def _run(self): pass
    for cls in (Splitter, Junction, Merger):
        U = cls()
        U.design_results.clear()
        U._design()
        assert U.design_results['Flow rate'] == 0
    with pytest.raises(RuntimeError):
        design.add_design_basis_to_cls(Dummy, 'Flow rate', 'kg/hr', None)
    def flow_rate(self, units): return 0
    with pytest.raises(ValueError):
        design.define(flow_rate)
    assert 'Flow rate' in design
    assert 'Flow rate' in list(design)
    assert 'Flow rate' in repr(design)
    # Copy algorithm
    class Copied(bst.Unit, isabstract=True): pass
    copy_algorithm(Dummy, Copied)
    assert Copied._run is Dummy._run
    with pytest.raises(RuntimeError):
        copy_algorithm(Dummy, Copied, design=False, cost=False)
    class Copied2(bst.Unit, isabstract=True):
        def _cost(self): pass
    with pytest.raises(RuntimeError):
        copy_algorithm(Dummy)(Copied2)
    class Copied3(bst.Unit, isabstract=True):
        def _design(self): pass
    with pytest.raises(RuntimeError):
        copy_algorithm(Dummy, Copied3, cost=False)

if __name__ == '__main__':
    test_mixer_phenomena_equations()
    test_steam_mixer_solids_loading()
    test_reversed_splitter_and_separator()
    test_unit_factories()
    test_fluidized_catalytic_cracking()
    test_auxiliary_pressure_vessel()
    test_spray_dryer()
    test_drum_dryer_options()
    test_thermal_oxidizer()
    test_terminal()
    test_junction_and_scaler()
    test_enzyme_treatment()
    test_vent_scrubber()
    test_screw_feeder_and_pellet_mills()
    test_single_phase_reactor_options()
    test_isenthalpic_valve_low_pressure()
    test_pump_types()
    test_vacuum_system_deduced_from_unit()
    test_adsorption_column_heuristic_design()
    test_adsorption_isotherm_fits()
    test_cost_and_design_decorators()
