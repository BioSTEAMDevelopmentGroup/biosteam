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
from biosteam import TEA, System

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
    
    with pytest.raises(ValueError):
        tea.depreciation = np.array([0.5, 0.1, 'bad'])
        
    with pytest.raises(ValueError):
        tea.depreciation = 'MACRS12312'
        
    with pytest.raises(ValueError):
        tea.depreciation = 'bad'


if __name__ == '__main__':
    test_depreciation_schedule()