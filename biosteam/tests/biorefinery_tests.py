# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
__all__ = ('test_biorefineries', 'test_lipidcane', 'test_cornstover')

def test_biorefineries():
    test_lipidcane()
    test_cornstover()

def test_lipidcane():
    from biorefineries.lipidcane import system
    IRR = system.lipidcane_tea.IRR
    assert 0.165 < IRR < 0.180, ("IRR of the lipid-cane biorefinery changed significantly "
                                f"to {IRR}")
    
def test_cornstover():
    from biorefineries.cornstover import system
    MESP = system.cornstover_tea.solve_price(system.ethanol)
    assert 0.65 < MESP < 0.80, ("MESP of the corn stover biorefinery changed significantly "
                               f"to {MESP} USD/kg")
    
if __name__ == '__main__':
    test_biorefineries()
    