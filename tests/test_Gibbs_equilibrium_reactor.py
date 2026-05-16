# -*- coding: utf-8 -*-
"""
"""
import numpy as np
import biosteam as bst
from thermosteam.constants import R
from numpy.testing import assert_allclose
from scipy.optimize import differential_evolution, Bounds

def test_equilibrium_reactor_NO2_N2O4():
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
    
def test_equilibrium_reactor_gas_reforming_and_water_gas_shift():
    bst.settings.set_thermo(['CH4', 'H2O', 'CO', 'H2', 'CO2'], pkg='ideal gas', cache=True)
    scenarios = [
        (273.15 + 1000, 101325),
    ]
    for T, P in scenarios:
        feed = bst.Stream(CH4=100, H2O=300, phase='g', T=T, P=P)
        EqR = bst.EquilibriumReactor(
            ins=feed, outs=('gas',), T=T, P=P,
            method='differential evolution'
        )
        EqR.simulate()
        actual_product = EqR.outs[0]
        
        # Gas reforming
        # theoretical_feed = bst.Stream(CH4=1, H2O=1, T=T, P=P)
        # theoretical_product = theoretical_feed.copy()
        rxn_gas_reforming = bst.Reaction('CH4 + H2O -> CO + 3H2', X=1, reactant='CH4')
        # rxn_gas_reforming(theoretical_product)
        # deltaG = (theoretical_product.G - theoretical_feed.G) / theoretical_feed.F_mol
        
        # Calculate Equilibrium Constant (Kp)
        # Formula: Delta G = -RT ln(Kp)
        # Kp_gas_reforming = np.exp(-deltaG / (R * T))
        Kp_gas_reforming = 9571.032284891147
        
        # Water-gas shift 
        # theoretical_feed = bst.Stream(H2O=1, CO=1, T=T, P=P)
        # theoretical_product = theoretical_feed.copy()
        rxn_water_gas_shift = bst.Reaction('H2O + CO -> CO2 + H2', X=1, reactant='CO')
        # rxn_water_gas_shift(theoretical_product)
        # deltaG = (theoretical_product.G - theoretical_feed.G) / theoretical_feed.F_mol
        
        # Calculate Equilibrium Constant (Kp)
        # Formula: Delta G = -RT ln(Kp)
        # Kp_water_gas_shift = np.exp(-deltaG / (R * T))
        Kp_water_gas_shift = 0.6343839656680074
        
        rxs = bst.SeriesReaction([
            rxn_gas_reforming, rxn_water_gas_shift
        ])
        desired_product = feed.copy()
        desired_product
        def objective(X):
            X[X < 1e-9] = 1e-9
            X[X > 1 - 1e-9] = 1 - 1e-9
            rxs.X[:] = X
            desired_product.copy_like(feed)
            rxs(desired_product)
            imol = desired_product.imol
            F_mol = desired_product.F_mol
            Kp_gas_reforming_actual = (imol['CO'] * imol['H2']**3) / (imol['CH4'] * imol['H2O']) * (P / 1e5 / F_mol)**2
            Kp_water_gas_shift_actual = (imol['CO2'] * imol['H2']) / (imol['CO'] * imol['H2O'])
            return (Kp_gas_reforming_actual - Kp_gas_reforming) ** 2 + (Kp_water_gas_shift_actual - Kp_water_gas_shift) ** 2
        
        lb = np.zeros(2) + 1e-9
        ub = np.ones(2) - 1e-9
        bounds = Bounds(lb, ub)
        solution = differential_evolution(
            objective,
            bounds=bounds,
            tol=1e-9,
        )
        rxs.X[:] = solution.x
        desired_product.copy_like(feed)
        rxs(desired_product)
        assert_allclose(actual_product.mol, desired_product.mol, rtol=1e-2, atol=1e-2)
    
if __name__ == '__main__':
    test_equilibrium_reactor_NO2_N2O4()
    test_equilibrium_reactor_gas_reforming_and_water_gas_shift()