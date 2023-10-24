# -*- coding: utf-8 -*-
"""
"""
import biosteam as bst
from numpy.testing import assert_allclose

def test_flash_energy_balance_with_pressure_drop():
    settings = bst.settings
    settings.set_thermo(['Water', 'Ethanol'], cache=True)
    settings.mixture.include_excess_energies = True
    feed = bst.Stream(None, Water=900, Ethanol=28000,
                      T=450, P=2000000)
    F1 = bst.Flash(None, ins=feed, P=101325, T=373.15) 
    F1.simulate()
    cooling, = F1.heat_utilities
    assert_allclose(F1.Hnet, cooling.unit_duty)
    
if __name__ == '__main__':
    test_flash_energy_balance_with_pressure_drop()