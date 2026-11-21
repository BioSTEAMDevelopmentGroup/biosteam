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

class SimpleTEA(bst.TEA):
    def _FOC(self, FCI):
        return 0.05 * FCI

def create_tea(**kwargs):
    bst.settings.set_thermo(['Water'], cache=True)
    bst.F.set_flowsheet('tea_coverage')
    feed = bst.Stream('feed', Water=1000, price=0.001)
    product = bst.Stream('product', price=0.01)
    P1 = bst.Pump('P1', feed, product)
    sys = bst.System('sys', path=[P1])
    sys.simulate()
    options = dict(
        IRR=0.10, duration=(2020, 2030), depreciation='MACRS7',
        income_tax=0.21, operating_days=330, lang_factor=None,
        construction_schedule=(0.4, 0.6), startup_months=3,
        startup_FOCfrac=1, startup_VOCfrac=0.75, startup_salesfrac=0.5,
        WC_over_FCI=0.05, finance_interest=0.08, finance_years=5,
        finance_fraction=0.4,
    )
    options.update(kwargs)
    return SimpleTEA(sys, **options)

def test_tea_properties():
    tea = create_tea()
    assert tea.units and tea.feeds and tea.products
    assert tea.start_year == 2
    tea.start_year = None
    assert tea.start_year == 2
    tea.startup_months = 6
    assert_allclose(tea.startup_months, 6)
    assert tea.purchase_cost > 0
    assert tea.FOC > 0 and tea.AOC > tea.FOC
    assert_allclose(tea.working_capital, 0.05 * tea.FCI, rtol=1e-6)
    assert tea.annual_depreciation > 0
    assert np.isfinite(tea.ROI) and np.isfinite(tea.PBP)
    assert np.isfinite(tea.net_earnings)
    assert len(tea.net_earnings_array) == 12
    assert tea.CAPEX_table() is NotImplemented
    assert tea.FOC_table() is NotImplemented
    assert 'SimpleTEA' in repr(tea)
    assert callable(tea.save_report.__call__)
    new = tea.copy(bst.System('sys2', path=[tea.system.units[0]]))
    assert new.system._TEA is new

def test_tea_negative_earnings():
    tea = create_tea()
    tea.system.flowsheet.stream.product.price = 0.
    assert tea.net_earnings < 0

def test_tea_with_interest_during_construction_and_inflation():
    tea = create_tea(accumulate_interest_during_construction=True, inflation_rate=0.02)
    table = tea.get_cashflow_table()
    assert table.shape[0] == 12
    assert np.isfinite(tea.NPV)
    IRR = tea.solve_IRR()
    assert 0 < IRR < 100
    IRR_bounded = tea.solve_IRR(bounds=(0.01, 50))
    assert_allclose(IRR, IRR_bounded, rtol=5e-2)
    IRR_unfinanced = tea.solve_IRR(financing=False)
    assert 0 < IRR_unfinanced < 100
    IRR_unfinanced_bounded = tea.solve_IRR(financing=False, bounds=(0.01, 50))
    assert_allclose(IRR_unfinanced, IRR_unfinanced_bounded, rtol=5e-2)

def test_tea_solve_price_and_production_costs():
    tea = create_tea()
    product = tea.system.flowsheet.stream.product
    price = tea.solve_price(product)
    assert price > 0
    with pytest.raises(ValueError):
        tea.solve_price(bst.Stream('empty'))
    costs = tea.production_costs([product], with_annual_depreciation=False)
    assert costs > 0
    df = tea.VOC_table([product])
    assert df is not None
    # Depreciation schedule longer than plant lifetime
    tea.depreciation = 'MACRS20'
    with pytest.raises(RuntimeError):
        tea.get_cashflow_table()
    tea.depreciation = 'SL'
    assert len(tea._get_depreciation_array()) == 10

if __name__ == '__main__':
    test_tea_properties()
    test_tea_negative_earnings()
    test_tea_with_interest_during_construction_and_inflation()
    test_tea_solve_price_and_production_costs()
