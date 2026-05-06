# -*- coding: utf-8 -*-
"""
"""
import numpy as np
import biosteam as bst
from thermosteam.constants import R
from numpy.testing import assert_allclose

def test_equilibrium_reactor():
    bst.settings.set_thermo(['N2O4', 'NO2'], pkg='ideal gas', cache=True)
    scenarios = [
        (273.15, 101325),
    ]
    for T, P in scenarios:
        feed = bst.Stream(N2O4=100, phase='g', T=T, P=P)
        EqR = bst.EquilibriumReactor(
            ins=feed, outs=('gas',), T=T, P=P 
        )
        EqR.simulate()
        actual_product = EqR.outs[0]
        
        theoretical_product = feed.copy()
        rxn = bst.Reaction('N2O4 -> 2NO2', X=1)
        rxn(theoretical_product)
        deltaG = (theoretical_product.G - feed.G) / feed.F_mol
        
        # Calculate Equilibrium Constant (Kp)
        # Formula: Delta G = -RT ln(Kp)
        Kp = np.exp(-deltaG / (R * T))
        
        # Solve for Degree of Dissociation (X)
        # Derived from Kp = (P_NO2^2) / P_N2O4
        # Using mole fractions: Kp = (4 * X^2 * P) / (1 - X^2)
        rxn.X = np.sqrt(Kp / (4 * (P / 101325) + Kp))
        
        # Calculate Equilibrium Moles
        desired_product = feed.copy()
        rxn(desired_product)
        
        assert_allclose(actual_product.mol, desired_product.mol, rtol=1e-2, atol=1e-2)
    
if __name__ == '__main__':
    test_equilibrium_reactor()