# -*- coding: utf-8 -*-
"""
"""
import numpy as np
import biosteam as bst
import matplotlib.pyplot as plt
from scipy.optimize import (
    differential_evolution, minimize,
    LinearConstraint, Bounds
)

def Gibbs_equilibrium_objective(
        mass, # Molar flow rate of product chemicals
        MWs, # Molecular weights
        IDs, # IDs of chemicals
        product, # Product stream (constant T and P)
        main_phase, # Initial guess phase
        phase_equilibrium, # Name of phase equilibrium. Valid inputs include None, 'vle', 'lle', and 'vlle'
    ):
    mask = mass < 0
    penalty = mask.any()
    if penalty: mass[mask] = 0
    product.empty()
    product[main_phase].imol[IDs] = mass / MWs
    if phase_equilibrium is not None:
        eq = getattr(product, phase_equilibrium)
        eq(T=eq.T, P=eq.P)
    return product.G

def minimize_Gibbs_free_energy(
        product, 
        IDs=None, # Names of potential products
        method=None,
    ):
    if method is None: method = 'differential evolution'
    F_mass = product.F_mass
    mol_norm = product.mol / F_mass # Normalize
    chemicals = product.chemicals
    atoms = chemicals.formula_array @ mol_norm 
    index, = np.where(atoms)
    atoms = atoms[index]
    formula_array = chemicals.formula_array[index]
    if IDs is None: 
        IDs = chemicals.IDs
        MWs = chemicals.MW
    else:
        index = chemicals.get_index(IDs)
        formula_array = formula_array[:, index]
        MWs = chemicals.MW[index]
    N = len(IDs)
    lb = np.zeros(N)
    ub = np.ones(N)
    bounds = Bounds(lb, ub)
    for i in range(N):
        formula = formula_array[:, i]
        index, = np.where(formula)
        ub_by_atom = (atoms[index] / formula[index]).min() * MWs[i]
        if ub_by_atom < ub[i]: ub[i] = ub_by_atom
    atomic_balance = LinearConstraint(
        formula_array / MWs, atoms, atoms
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
    f_args = (MWs, IDs, product, main_phase, phase_equilibrium)
    if method == 'differential evolution':
        polish = lambda *args, **kwargs: minimize(
            *args, **kwargs, 
            args=f_args,
            method='COBYLA', 
        )
        solution = differential_evolution(
            Gibbs_equilibrium_objective,
            args=f_args,
            bounds=bounds,
            constraints=[atomic_balance],
            polish=polish,
            tol=1e-5,
            seed=0,
        )
    elif method == 'COBYLA':
        solution = minimize(
            Gibbs_equilibrium_objective,
            x0=np.array(mol_norm * MWs),
            args=f_args,
            bounds=bounds,
            method=method, 
            constraints=[atomic_balance],
        )
    else:
        raise ValueError(
           f"invalid method {method!r}; "
            "only 'SHGO' and 'differential evolution' "
            "are valid methods"
        )
    product.empty()
    product[main_phase].imol[IDs] = solution.x / MWs
    if phase_equilibrium is not None:
        eq = getattr(product, phase_equilibrium)
        eq(T=eq.T, P=eq.P)
    product.F_mass = F_mass
    return solution
    
def plot_Gibbs_equilibrium_reaction_surface():
    colors = bst.utils.GG_colors
    bst.settings.set_thermo(['N2O4', 'NO2'], pkg='ideal gas')
    feed = bst.Stream(N2O4=100, T=298.15, P=101325, phase='g')
    product = feed.copy()
    rxn = bst.Reaction('N2O4 -> 2NO2', X=1)
    
    def g(X):
        rxn.X = X
        product.copy_flow(feed)
        rxn(product)
        return product.G / product.F_mass
    
    Xs = np.linspace(0, 1)    
    gs = [g(X) for X in Xs]
    plt.plot(Xs, gs, c=colors.orange.RGBn)
    plt.xlim([0, 1])
    plt.xlabel('X')
    plt.ylabel('G')

def test_Gibbs_equilibrium_reaction_surface():
    bst.settings.set_thermo(['N2O4', 'NO2'], pkg='ideal gas')
    feed = bst.Stream(N2O4=100, T=300, P=101325, phase='g')
    product = feed.copy()
    rxn = bst.Reaction('N2O4 -> 2NO2', X=1)
    
    def solve_opt(T, P):
        product.T = T
        product.P = P
        product.copy_flow(feed)
        minimize_Gibbs_free_energy(product)
        return product.mol.to_array()
    
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
        solution = differential_evolution(g, bounds=bounds, tol=1e-9, seed=0)
        X_actual = solution.x
        g(X_actual)
        return product.mol.to_array()
    
    Ts = np.linspace(273.15, 500, 5)
    Ps = np.linspace(101325, 50 * 101325, 5)
    actual = np.array([[solve_opt(T, P) for T in Ts] for P in Ps])
    desired = np.array([[solve_exact(T, P) for T in Ts] for P in Ps])
    np.testing.assert_allclose(actual, desired, rtol=1e-2, atol=1e-2)
    
    