# -*- coding: utf-8 -*-
"""
"""
import biosteam as bst
from numpy.testing import assert_allclose

def test_tire():
    bst.settings.set_thermo(
        ['Rubber', 'Ash', 'Water', 'Char'],
        db='BioSTEAM'
    )
    bst.settings.chemicals.define_group(
        'Tire',
        ['Rubber', 'Ash', 'Water'], 
        [94.88, 4.2, 0.92],
    )
    tire = bst.Stream(Tire=100, T=350)
    Hnet = tire.Hnet
    assert_allclose(Hnet, -94972.77484287221)
    assert_allclose(tire.F_mass, 115.65184145619916)
    assert_allclose(tire.F_mol, 100)
    
    char = bst.Stream(Char=100, T=350)
    Hnet = char.Hnet
    F_mass = char.F_mass
    assert_allclose(Hnet, 6281.35)
    assert_allclose(F_mass, 1201.07)
    
if __name__ == '__main__':
    test_tire()