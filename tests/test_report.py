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
import numpy as np
import os
from numpy.testing import assert_allclose

@pytest.fixture
def system():
    bst.main_flowsheet.clear()
    MockChemical = bst.Chemical('MockChemical', phase='l', search_db=False, default=True)
    bst.settings.set_thermo([MockChemical])
    with bst.System('MockSystem') as system:
        feed = bst.Stream('feed', MockChemical=1e6, units='kg/hr', price=0.150)
        product = bst.Stream('product', price=0.151)
        M1 = bst.MixTank('M1', feed, 's1')
        H1 = bst.HXutility('H1', M1-0, product, T=325, V=0.)
    system.simulate()
    return system
    
@pytest.fixture
def tea(system):
    
    class MockTEA(bst.TEA):
        def _DPI(self, installed_equipment_cost): return installed_equipment_cost
        def _TDC(self, DPI): return DPI
        def _FCI(self, TDC): return TDC
        def _FOC(self, FCI): return FCI
        
    tea = MockTEA(system,
        startup_months=0, startup_FOCfrac=0, startup_VOCfrac=0,
        startup_salesfrac=0, finance_interest=0, finance_years=0, 
        finance_fraction=0, IRR=0.15, duration=(2018, 2038),
        depreciation='MACRS7', income_tax=0.35, lang_factor=None,
        operating_days=200, construction_schedule=(0.4, 0.6), WC_over_FCI=0.05,
    )
    return tea

def test_unit_result_tables(system):
    results = bst.report.unit_result_tables(system.units)
    columns = ['Units', 'H1', 'Units', 'M1']
    assert [i for df in results for i in df] == columns

def test_cost_table(tea):
    df = bst.report.cost_table(tea)
    columns = ['Unit operation', 
               'Purchase cost (10^6 USD)', 
               'Utility cost (10^6 USD/yr)',
               'Installed cost (10^6 USD)']
    assert list(df) == columns
    assert list(df.index) == ['H1', 'M1']

def test_heat_utilities_table(system):
    df, = bst.report.heat_utility_tables(system.units)
    columns = ['Unit operation', 'Duty (kJ/hr)', 'Flow (kmol/hr)', 'Cost (USD/hr)']
    assert list(df) == columns
    assert list(df.index) == ['H1']
    
def test_power_utilities_table(system):
    df = bst.report.power_utility_table(system.units)
    columns = ['Unit Operation', 'Rate (kW)', 'Cost (USD/hr)']
    assert list(df) == columns
    assert list(df.index) == ['M1']
    
def test_stream_table(system):
    df = bst.report.stream_table(system.streams)
    values = np.array(
        [['-',      'H1',       'M1'],
         ['M1',     '-',        'H1'],
         ['liquid', 'liquid',   'liquid'],
         [1000000.0, 1000000.0, 1000000.0],
         ['',        '',        ''],
         [100.0,     100.0,     100.0]], 
        dtype=object
    )
    index = ['Source', 'Sink', 'Phase', 'flow (kg/hr)', 'Composition [%]:', 'MockChemical']
    IDs = ['feed', 'product', 's1']
    assert (df.values == values).all()
    assert list(df) == IDs
    assert list(df.index) == index

def test_save_report(system, tea):
    assert system.TEA is tea
    system.save_report("report.xlsx") # Make sure it runs
    os.remove("report.xlsx")
    # TODO: More robust test for this method
