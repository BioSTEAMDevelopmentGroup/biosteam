# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2021-2022, Yalin Li <zoe.yalin.li@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""

import pytest, numpy as np
from numpy.testing import assert_allclose
from biosteam import TEA, System, Stream, Unit, Chemical, settings
from biorefineries.tea import create_cellulosic_ethanol_tea

def test_depreciation_schedule():    
    correct_arrs = {
		# https://www.irs.gov/pub/irs-pdf/p946.pdf, table A-1
		'MACRS5': np.array([0.2000, 0.3200, 0.1920, 0.1152, 0.1152, 0.0576]),
		# Below from Excel function (SLN, SYD, DB, DDB)
		'SL5'   : np.array([0.2000, 0.2000, 0.2000, 0.2000, 0.2000]),
		'DDB5'  : np.array([0.4000, 0.2400, 0.1440, 0.0864, 0.0518]),
		'SYD5'  : np.array([0.3333, 0.2667, 0.2000, 0.1333, 0.0667]),
    }
    
    empty_sys = System('Empty')
    tea = TEA(
        system=empty_sys, IRR=0.15,
        duration=(2021, 2031),
        depreciation='MACRS7', income_tax=0.35,
        operating_days=200, lang_factor=3,
        construction_schedule=(0.4, 0.6), 
        startup_months=0, startup_FOCfrac=0, startup_VOCfrac=0,
        startup_salesfrac=0, finance_interest=0, finance_years=0, 
        finance_fraction=0, WC_over_FCI=0.05
    )

    for k, v in correct_arrs.items():
        tea.depreciation = k
        assert_allclose(tea._get_depreciation_array(), v, atol=1e-4)
    
    defined = np.array([0.4000, 0.2400, 0.1400, 0.0900, 0.1300])
    tea.depreciation = defined
    assert_allclose(tea._get_depreciation_array(), defined)
    
    with pytest.raises(TypeError):
        tea.depreciation = np.array([0.5, 0.1, 'bad'])
        
    with pytest.raises(ValueError):
        tea.depreciation = 'MACRS12312'
        
    with pytest.raises(ValueError):
        tea.depreciation = 'bad'

def test_cashflow():
    settings.set_thermo([
        Chemical('Dummy', default=True, phase='s', MW=1, search_db=False)
    ])
    
    cost = Stream(Dummy=9793.983511363867 - 1280.916880939845, price=1) # Includes co-product electricity credits
    ethanol = Stream(flow=[21978.374283953395], price=0.7198608114634679)
    
    class MockCellulosicEthanolBiorefinery(Unit):
        _N_ins = _N_outs = 1
        
        def _run(self): pass
        
        def _cost(self):
            self.baseline_purchase_costs['Biorefinery'] = 85338080.48935215
            
    class OSBL(Unit):
        N_outs = _N_ins = 0
        
        def _run(self): pass
        
        def _cost(self):
            self.baseline_purchase_costs['Biorefinery'] = 122287135.13152598
            
    unit = MockCellulosicEthanolBiorefinery(ins=cost, outs=ethanol)
    osbl = OSBL()
    sys = System.from_units(units=[unit, osbl])
    sys.simulate()
    tea = create_cellulosic_ethanol_tea(sys, OSBL_units=[osbl])
    table = tea.get_cashflow_table()
    assert_allclose(tea.NPV, 32131936.781448975)
    assert_allclose(tea.NPV, table['Cumulative NPV [MM$]'].iloc[-1]*1e6)
    tea.IRR = tea.solve_IRR()
    assert_allclose(tea.NPV, 0, atol=100)
    assert_allclose(tea.IRR, 0.11246761316144724)
    tea.IRR = 0.10
    ethanol.price = tea.solve_price(ethanol)
    assert_allclose(ethanol.price, 0.6952016482242149)
    assert_allclose(tea.NPV, 0, atol=100)

if __name__ == '__main__':
    test_depreciation_schedule()
    test_cashflow()