# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2021-2024, Yalin Li <mailto.yalin.li@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""

import pytest, numpy as np, biosteam as bst
from numpy.testing import assert_allclose

def test_tea():   
    cost = bst.decorators.cost
    # Total installed equipment cost to be $1 MM
    bst.CE = CE = bst.design_tools.CEPCI_by_year[2013]
    @cost('Fake scaler', 'Lumped cost', CE=CE, cost=1e6, S=1, n=1, BM=1)
    class LumpedCost(bst.Unit):
        '''Does nothing but adding given costs.'''
        _units = {'Fake scaler': ''}
        
        def _design(self):
            self.design_results['Fake scaler'] = 1
        
    class TEA(bst.TEA):
        def __init__(
                self,
                system,
                FOC_over_installed=0.5, # annual O&M
                DPI_over_installed=(1+1),
                TDC_over_DPI=(1+0.2)*(1+0.4), # 20% non-installed & 40% indirect
                FCI_over_TDC=1,
                **kwargs,
                ):
            self.FOC_over_installed = FOC_over_installed
            self.DPI_over_installed = DPI_over_installed
            self.TDC_over_DPI = TDC_over_DPI
            self.FCI_over_TDC = FCI_over_TDC
            bst.TEA.__init__(self, system, **kwargs)
    
        def _FOC(self, installed_equipment_cost): # fixed operating cost
            return installed_equipment_cost*self.FOC_over_installed
        
        def _DPI(self, installed_equipment_cost): # direct permanent investment
            return installed_equipment_cost*self.DPI_over_installed
    
        def _TDC(self, DPI): # total depreciable cost
            return DPI*self.TDC_over_DPI
    
        def _FCI(self, TDC): # fixed capital investment
            return TDC*self.FCI_over_TDC
        
    bst.settings.set_thermo([bst.Chemical('Water')])
    reactant = bst.Stream('reactant', Water=1, units='kg/hr')
    # Total annual sales to be $2.5 MM
    product = bst.Stream('product', Water=1, price=2.5e6/365/24, units='kg/hr')
    
    U101 = LumpedCost('U101', ins=reactant, outs=product)
    sys = bst.System('sys', path=(U101,))
    sys.simulate()

    tea = TEA(
        system=sys,
        IRR=0.1,
        duration=(2013, 2013+15), # 15 years
        income_tax=0.21+0.1,
        construction_schedule=(1,),
        depreciation='MACRS7',
        operating_days=365,
        startup_months=0,
        startup_FOCfrac=1,
        startup_VOCfrac=1,
        startup_salesfrac=1,
        lang_factor=None,
        WC_over_FCI=0.05,
        finance_interest=0.08,
        finance_years=10,
        finance_fraction =0.6,
        accumulate_interest_during_construction=False,
        )

    # Below test NPV and discounted cashflow calculation
    table = tea.get_cashflow_table()
    data = np.array([
        # Depreciable capital [MM$]
        # Fixed capital investment [MM$]
        # Working capital [MM$]
        # Depreciation [MM$]
        # Loan [MM$]
       [ 3.36      ,  3.36      ,  0.168     ,  0.        ,  2.016     ,
        # Loan interest payment [MM$]
        # Loan payment [MM$]
        # Loan principal [MM$]
        # Annual operating cost (excluding depreciation) [MM$]
        # Sales [MM$]
         0.16128   ,  0.        ,  2.016     ,  0.        ,  0.        ,
        # Tax [MM$]
        # Incentives [MM$]
        # Taxed earnings [MM$]
        # Forwarded losses [MM$]
        # Net earnings [MM$]
         0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        # Cash flow [MM$]     
        # Discount factor
        # Net present value (NPV) [MM$]
        # Cumulative NPV [MM$]
        -1.67328   ,  1.        , -1.67328   , -1.67328   ],
       [ 0.        ,  0.        ,  0.        ,  0.480144  ,  0.        ,
         0.16128   ,  0.30044345,  1.87683655,  1.68      ,  2.5       ,
         0.01221789,  0.        ,  0.03941255,  0.        ,  0.02719466,
         0.50733866,  0.90909091,  0.46121696, -1.21206304],
       [ 0.        ,  0.        ,  0.        ,  0.822864  ,  0.        ,
         0.15014692,  0.30044345,  1.72654003,  1.68      ,  2.5       ,
         0.        ,  0.        ,  0.        ,  0.        , -0.30330745,
         0.51955655,  0.82644628,  0.42938558, -0.78267746],
       [ 0.        ,  0.        ,  0.        ,  0.587664  ,  0.        ,
         0.1381232 ,  0.30044345,  1.56421978,  1.68      ,  2.5       ,
         0.        ,  0.        ,  0.        , -0.30330745, -0.06810745,
         0.51955655,  0.7513148 ,  0.39035053, -0.39232693],
       [ 0.        ,  0.        ,  0.        ,  0.419664  ,  0.        ,
         0.12513758,  0.30044345,  1.38891391,  1.68      ,  2.5       ,
         0.        ,  0.        ,  0.        , -0.3714149 ,  0.09989255,
         0.51955655,  0.68301346,  0.35486412, -0.03746282],
       [ 0.        ,  0.        ,  0.        ,  0.300048  ,  0.        ,
         0.11111311,  0.30044345,  1.19958358,  1.68      ,  2.5       ,
         0.        ,  0.        ,  0.        , -0.27152235,  0.21950855,
         0.51955655,  0.62092132,  0.32260374,  0.28514093],
       [ 0.        ,  0.        ,  0.        ,  0.299712  ,  0.        ,
         0.09596669,  0.30044345,  0.99510681,  1.68      ,  2.5       ,
         0.05202753,  0.        ,  0.16783075, -0.0520138 ,  0.16781702,
         0.46752902,  0.56447393,  0.26390794,  0.54904887],
       [ 0.        ,  0.        ,  0.        ,  0.300048  ,  0.        ,
         0.07960854,  0.30044345,  0.77427191,  1.68      ,  2.5       ,
         0.06804765,  0.        ,  0.21950855,  0.        ,  0.1514609 ,
         0.4515089 ,  0.51315812,  0.23169546,  0.78074432],
       [ 0.        ,  0.        ,  0.        ,  0.149856  ,  0.        ,
         0.06194175,  0.30044345,  0.53577021,  1.68      ,  2.5       ,
         0.11460717,  0.        ,  0.36970055,  0.        ,  0.25509338,
         0.40494938,  0.46650738,  0.18891187,  0.9696562 ],
       [ 0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
         0.04286162,  0.30044345,  0.27818838,  1.68      ,  2.5       ,
         0.16106253,  0.        ,  0.51955655,  0.        ,  0.35849402,
         0.35849402,  0.42409762,  0.15203646,  1.12169266],
       [ 0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
         0.02225507,  0.30044345,  0.        ,  1.68      ,  2.5       ,
         0.16106253,  0.        ,  0.51955655,  0.        ,  0.35849402,
         0.35849402,  0.38554329,  0.13821496,  1.25990762],
       [ 0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
         0.        ,  0.        ,  0.        ,  1.68      ,  2.5       ,
         0.2542    ,  0.        ,  0.82      ,  0.        ,  0.5658    ,
         0.5658    ,  0.3504939 ,  0.19830945,  1.45821707],
       [ 0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
         0.        ,  0.        ,  0.        ,  1.68      ,  2.5       ,
         0.2542    ,  0.        ,  0.82      ,  0.        ,  0.5658    ,
         0.5658    ,  0.31863082,  0.18028132,  1.63849839],
       [ 0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
         0.        ,  0.        ,  0.        ,  1.68      ,  2.5       ,
         0.2542    ,  0.        ,  0.82      ,  0.        ,  0.5658    ,
         0.5658    ,  0.28966438,  0.16389211,  1.80239049],
       [ 0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
         0.        ,  0.        ,  0.        ,  1.68      ,  2.5       ,
         0.2542    ,  0.        ,  0.82      ,  0.        ,  0.5658    ,
         0.5658    ,  0.26333125,  0.14899282,  1.95138332],
       [ 0.        ,  0.        , -0.168     ,  0.        ,  0.        ,
         0.        ,  0.        ,  0.        ,  1.68      ,  2.5       ,
         0.2542    ,  0.        ,  0.82      ,  0.        ,  0.5658    ,
         0.7338    ,  0.23939205,  0.17566589,  2.1270492 ]])
    
    assert_allclose(table.values, data, atol=1e-4)
    assert_allclose(tea.NPV, table.iloc[-1,-1]*1e6, atol=1e-4)
    total_interest_payment1 = table['Loan payment [MM$]'].sum()
    total_interest_payment2 = (
        table['Loan interest payment [MM$]'].iloc[1:].sum() + # payment during construction (year 0) is equity/cash
            table['Loan principal [MM$]'].iloc[0])
    assert_allclose(total_interest_payment1, total_interest_payment2, atol=1e-4)

    # Below test the depreciation schedule
    correct_arrs = {
		# https://www.irs.gov/pub/irs-pdf/p946.pdf, table A-1
		'MACRS5': np.array([0.2000, 0.3200, 0.1920, 0.1152, 0.1152, 0.0576]),
		# Below from Excel function (SLN, SYD, DB, DDB)
		'SL5'   : np.array([0.2000, 0.2000, 0.2000, 0.2000, 0.2000]),
		'DDB5'  : np.array([0.4000, 0.2400, 0.1440, 0.0864, 0.0518]),
		'SYD5'  : np.array([0.3333, 0.2667, 0.2000, 0.1333, 0.0667]),
    }
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


if __name__ == '__main__':
    test_tea()