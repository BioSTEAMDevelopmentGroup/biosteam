# -*- coding: utf-8 -*-
"""
"""
import numpy as np
import biosteam as bst
from scipy.optimize import (
    differential_evolution, minimize,
    LinearConstraint, Bounds
)

def Gibbs_equilibrium_objective(
        mol, # Molar flow rate of product chemicals
        IDs, # IDs of chemicals
        product, # Product stream (constant T and P)
        main_phase, # Initial guess phase
        phase_equilibrium, # Name of phase equilibrium. Valid inputs include None, 'vle', 'lle', and 'vlle'
    ):
    mask = mol < 0
    penalty = mask.any()
    if penalty: mol[mask] = 0
    product.empty()
    product[main_phase].imol[IDs] = mol
    if phase_equilibrium is not None:
        eq = getattr(product, phase_equilibrium)
        eq(T=eq.T, P=eq.P)
    return product.G

def minimize_Gibbs_free_energy(
        product, 
        IDs=None, # Names of potential products
    ):
    F_mass = product.F_mass
    z_mol = product.mol / product.F_mol
    chemicals = product.chemicals
    atoms = chemicals.formula_array @ z_mol # Normalize
    index, = np.where(atoms)
    atoms = atoms[index]
    formula_array = chemicals.formula_array[index]
    if IDs is None: 
        IDs = chemicals.IDs
    else:
        formula_array = chemicals.formula_array[:, chemicals.get_index(IDs)]
    N = len(IDs)
    lb = np.zeros(N)
    ub = np.zeros(N)
    for i in range(N):
        formula = formula_array[:, i]
        index, = np.where(formula)
        ub[i] = (atoms[index] / formula[index]).min() 
    bounds = Bounds(lb, ub)
    linear_constraint = LinearConstraint(
        formula_array, atoms, atoms
    )
    match product.phases:
        case ('g', 'l'):
            phase_equilibrium = 'vle' 
            main_phase = 'g'
        case ('L', 'l'):
            phase_equilibrium = 'lle' 
            main_phase = 'l'
        case ('L', 'g', 'l'):
            phase_equilibrium = 'vlle'
            main_phase = 'l'
        case ('l', 's'):
            phase_equilibrium = 'sle'
            main_phase = 'l'
        case [main_phase]:
            phase_equilibrium = None
        case _:
            raise RuntimeError(f'phase equilibrium for {product.phases!r} not supported')
    f_args = (IDs, product, main_phase, phase_equilibrium)
    polish = lambda *args, **kwargs: minimize(
        *args, **kwargs, 
        args=f_args,
        method='trust-constr', 
        options=dict(factorization_method='SVDFactorization'),
    )
    solution = differential_evolution(
        Gibbs_equilibrium_objective,
        args=f_args,
        bounds=bounds,
        constraints=[linear_constraint],
        polish=polish,
        tol=1e-6,
    )
    product.empty()
    product[main_phase].imol[IDs] = solution.x
    if phase_equilibrium is not None:
        eq = getattr(product, phase_equilibrium)
        eq(T=eq.T, P=eq.P)
    product.F_mass = F_mass
    solution.x[:] = product.mol
    return solution
    
def test_Gibbs_equilibrium_reaction_surface():
    bst.settings.set_thermo(['N2O4', 'NO2'], pkg='ideal gas')
    feed = bst.Stream(N2O4=100, T=300, P=101325, phase='g')
    product = feed.copy()
    rxn = bst.Reaction('N2O4 -> 2NO2', X=1)
    
    def solve_opt(T, P):
        product.T = T
        product.P = P
        product.copy_flow(feed)
        solution = minimize_Gibbs_free_energy(product)
        return solution.x
    
    def solve_exact(T, P):
        product.T = T
        product.P = P
        
        def g(X):
            rxn.X = X[0]
            product.copy_flow(feed)
            rxn(product)
            return product.G / product.F_mass
        
        lb = np.zeros(1)
        ub = np.ones(1)
        bounds = Bounds(lb, ub)
        solution = differential_evolution(g, bounds=bounds, tol=1e-6)
        X_actual = solution.x
        g(X_actual)
        return product.mol.to_array()
    
    Ts = np.linspace(273.15, 500, 5)
    Ps = np.linspace(101325, 50 * 101325, 5)
    actual = np.array([[solve_opt(T, P) for T in Ts] for P in Ps])
    desired = np.array([[solve_exact(T, P) for T in Ts] for P in Ps])
    np.testing.assert_allclose(desired, actual, rtol=1e-2, atol=1e-2)
    
    