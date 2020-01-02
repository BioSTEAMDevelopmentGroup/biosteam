# -*- coding: utf-8 -*-
"""
Created on Sun Dec  8 02:01:41 2019

@author: yoelr
"""

__all__ = ('test_biorefineries', 'test_lipidcane', 'test_cornstover')

def test_biorefineries():
    test_lipidcane()
    test_cornstover()

def test_lipidcane():
    from biorefineries.lipidcane import system
    IRR = system.lipidcane_tea.IRR
    assert 0.194 < IRR < 0.198, ("IRR of the lipid-cane biorefinery changed significantly "
                                f"from 0.1945 to {IRR}")
    
def test_cornstover():
    from biorefineries.cornstover import system
    MESP = system.cornstover_tea.solve_price(system.ethanol)
    assert 0.70 < MESP < 0.73, ("MESP of the corn stover biorefinery changed significantly "
                               f"from 0.747 to {MESP} USD/kg")
    