# -*- coding: utf-8 -*-
"""
Created on Sun Apr 28 16:26:41 2019

@author: yoelr
"""

import lipidcane as lc
import biosteam as bs
solve_IRR = lc.lipidcane_tea.solve_IRR
grid = bs.evaluation.Grid(system=lc.lipidcane_sys,
                          metric=solve_IRR,
                          ID='IRR')

Lipid_cane = lc.Lipid_cane # The feedstock stream
def set_feed_price(feedstock_price): Lipid_cane.price = feedstock_price
grid.addparam(element=Lipid_cane,
              setter=set_feed_price,
              values=[0.030, 0.035, 0.040],
              isolated=True)

grid.addparam(Lipid_cane, lc.set_lipid_fraction, [0.02, 0.05, 0.10])

U34 = bs.find('U34')
def set_efficiency(fermentation_efficiency): U34.reset(efficiency=fermentation_efficiency)
grid.addparam(U34, set_efficiency, [0.85, 0.90, 0.95])
grid.show()