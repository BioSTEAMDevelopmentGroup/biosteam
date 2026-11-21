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
import pandas as pd
import os
import tempfile
from numpy.testing import assert_allclose

def create_system(ID):
    bst.F.set_flowsheet(ID)
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    feed = bst.Stream('feed', Water=100, Ethanol=10, units='kg/hr', price=0.1)
    with bst.System(ID) as system:
        P1 = bst.Pump('P1', feed, P=5 * 101325)
        H1 = bst.HXutility('H1', P1-0, T=360, V=0)
        S1 = bst.Splitter('S1', H1-0, ['product', 'coproduct'], split=0.8)
    S1.outs[0].price = 0.5
    S1.outs[1].price = 0.2
    system.operating_hours = 8000
    system.simulate()
    return system

class LCAContext:
    # Set and reset global LCA settings.

    def __enter__(self):
        self.impact_indicators = bst.settings.impact_indicators.copy()
        self.power_CFs = bst.PowerUtility.characterization_factors.copy()
        bst.settings.define_impact_indicator('GWP', 'kg*CO2e')
        bst.PowerUtility.set_CF('GWP', 0.4, 0.3)
        return self

    def __exit__(self, *args):
        bst.settings.impact_indicators.clear()
        bst.settings.impact_indicators.update(self.impact_indicators)
        bst.PowerUtility.characterization_factors.clear()
        bst.PowerUtility.characterization_factors.update(self.power_CFs)

def set_stream_CFs(system):
    F = system.flowsheet
    F.feed.set_CF('GWP', 0.5)
    F.coproduct.set_CF('GWP', 0.3)
    system.define_process_impact(
        key='GWP', name='Process emissions', basis='kg',
        inventory=lambda: 10., CF=1.,
    )

def test_FOC_table_builder():
    builder = bst.report.FOCTableBuilder()
    builder.entry('Labor', np.array([1., 2.]), 'note')
    builder.entry('Maintenance', np.array([0.5, 0.25]))
    df = builder.table(['A', 'B'])
    assert df.shape == (3, 3)
    builder = bst.report.FOCTableBuilder()
    builder.entry('Labor', np.array([1.]))
    data, index, columns = builder.table(['A'], dataframe=False)
    assert index[-1] == 'Fixed operating cost (FOC)'
    assert len(columns) == 2
    assert_allclose(data[-1][1], 1., rtol=1e-2)

def test_voc_table():
    system = create_system('voc_table')
    df = bst.report.voc_table(system, 'product')
    assert df.shape[1] == 2
    assert ('Variable operating cost', '') in df.index
    df = bst.report.voc_table([system], [system.flowsheet.product],
                              system_names=['A'], with_products=True)
    assert df.shape[1] == 2
    assert ('Products', 'Product') in df.index
    data, index, columns = bst.report.voc_table(system, 'product', dataframe=False)
    assert len(index) == data.shape[0]
    system_2 = create_system('voc_table_2')
    df = bst.report.voc_table([system, system_2], 'product')
    assert df.shape[1] == 3
    assert np.isfinite(df.iloc[-1, 1:].values.astype(float)).all()

def test_lca_tables():
    system = create_system('lca_tables')
    product = system.flowsheet.product
    coproduct = system.flowsheet.coproduct
    with LCAContext():
        set_stream_CFs(system)
        df = bst.report.lca_displacement_allocation_table(system, 'GWP', [product])
        assert ('Total', '') in df.index
        df = bst.report.lca_displacement_allocation_table(
            [system], 'GWP', [[product]], system_names=['A'],
        )
        assert df.shape[1] == 2
        df = bst.report.lca_displacement_allocation_table(
            [system], 'GWP', [[product]], item_name='product',
        )
        assert df.shape[1] == 2
        # Note: process impacts with a basis other than kg are only tested 
        # in other tables because lca_displacement_allocation_table fails.
        system.define_process_impact(
            key='GWP', name='Catalyst', basis='MJ',
            inventory=lambda: 5., CF=2.,
        )
        df = bst.report.lca_characterization_factor_table(system)
        assert df.shape[1] == len(bst.settings.impact_indicators)
        assert ('Inputs', 'Feed') in df.index
        df = bst.report.lca_inventory_table(system, 'GWP', items=[[product]])
        assert df.shape[1] == 1
        df = bst.report.lca_inventory_table([system], items=[(product,)],
                                            system_names=['A'])
        assert list(df.columns) == ['A inventory [kg/yr]']
        df = bst.report.lca_property_allocation_factor_table([system], 'mass')
        assert_allclose(df.values.astype(float).sum(), 1., rtol=1e-2)
        df = bst.report.lca_property_allocation_factor_table(
            [system], 'mass', system_names=['A'],
        )
        assert list(df.columns) == ['A mass allocation factors']
        df = bst.report.lca_displacement_allocation_factor_table(
            [system], [[product]], 'GWP',
        )
        assert df.shape[1] == 1
        df = bst.report.lca_displacement_allocation_factor_table(
            [system], [[product]], 'GWP', system_names=['A'],
        )
        assert df.shape[1] == 1
        df = bst.report.environmental_impacts_table(
            system, [product], ['kg'], 'mass',
        )
        assert df.shape == (len(bst.settings.impact_indicators), 1)

def test_lca_electricity_production():
    bst.F.set_flowsheet('lca_electricity_production')
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    feed = bst.Stream('feed', Water=100, units='kg/hr', price=0.1)
    with bst.System('lca_electricity_production') as system:
        M1 = bst.Pump('M1', feed, 'product')
    system.operating_hours = 8000
    system.simulate()
    M1.power_utility.consumption = 0.
    M1.power_utility.production = 10.
    with LCAContext():
        feed.set_CF('GWP', 0.5)
        system.define_process_impact(
            key='GWP', name='Process emissions', basis='kg',
            inventory=lambda: 10., CF=1.,
        )
        df = bst.report.lca_characterization_factor_table(system)
        assert ('Outputs', 'Electricity production') in df.index
        df = bst.report.lca_inventory_table(system, 'GWP', items=[])
        assert ('Outputs', 'Electricity [kWhr/yr]') in df.index
        df = bst.report.lca_displacement_allocation_table(system, 'GWP', [M1.outs[0]])
        assert ('Outputs displaced', 'Electricity') in df.index
        df = bst.report.voc_table(system, 'product')
        assert ('Co-products & credits', 'Electricity production') in df.index

def test_unit_tables():
    system = create_system('unit_tables')
    units = system.units
    P1 = system.flowsheet.P1
    P1.reaction = tmo.Reaction('Ethanol -> Water', 'Ethanol', 0.5, basis='wt')
    tables = bst.report.unit_reaction_tables(units)
    assert len(tables) == 1
    del P1.reaction
    tables = bst.report.unit_result_tables(units, include_utilities=True,
                                           include_total_cost=True,
                                           include_installed_cost=True)
    assert len(tables) >= 2
    tables = bst.report.heat_utility_tables(units)
    assert len(tables) == 1
    assert tables[0].shape[1] == 4
    df = bst.report.power_utility_table(units)
    assert df.shape == (1, 3)

def test_other_utilities_table():
    system = create_system('other_utilities_table')
    F = system.flowsheet
    feed = F.feed
    feed.price = 0
    bst.stream_prices['Mock fee'] = 0.01
    bst.stream_prices['Mock credit'] = 0.02
    try:
        F.P1.define_utility('Mock fee', F.P1.ins[0])
        F.S1.define_utility('Mock credit', F.S1.outs[1])
        F.S1.outs[1].price = 0
        tables = bst.report.other_utilities_table(system.units)
        assert len(tables) == 2
        assert all([i.shape == (1, 3) for i in tables])
    finally:
        del bst.stream_prices['Mock fee']
        del bst.stream_prices['Mock credit']

def test_stream_tables():
    system = create_system('stream_tables')
    tables = bst.report.stream_tables(system.streams, P='Pa')
    assert len(tables) == 1
    df = tables[0]
    assert 'Water' in df.index and 'Ethanol' in df.index
    vapor = bst.Stream('vapor', Water=1, phase='g')
    solid = bst.Stream('solid', Water=1, phase='s')
    mixed = bst.Stream('mixed', Water=2)
    mixed.phases = ('L', 'S', 'l')
    empty = bst.Stream('empty')
    tmo.settings.set_thermo(['Water'], cache=True)
    other = bst.Stream('other', Water=1)
    df = bst.report.stream_table([vapor, solid, mixed, empty, other], percent=False)
    assert df.shape[1] == 5
    assert 'Composition:' in df.index

def test_mass_balance_tables():
    system = create_system('mass_balance_tables')
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    df = bst.report.water_mass_balance_table(system)
    assert df.shape[1] == 1
    assert_allclose(df.values[-1, 0], 0, atol=1)
    df = bst.report.mass_balance_table(system, 'Ethanol')
    assert df.shape[1] == 1

def test_tables_to_excel():
    system = create_system('tables_to_excel')
    tables = bst.report.unit_result_tables(system.units)
    with tempfile.TemporaryDirectory() as tmp:
        file = os.path.join(tmp, 'tables.xlsx')
        with pd.ExcelWriter(file) as writer:
            n_row = bst.report.tables_to_excel(tables, writer)
        assert n_row > 1
        assert os.path.exists(file)

if __name__ == '__main__':
    test_FOC_table_builder()
    test_voc_table()
    test_lca_tables()
    test_lca_electricity_production()
    test_unit_tables()
    test_other_utilities_table()
    test_stream_tables()
    test_mass_balance_tables()
    test_tables_to_excel()
