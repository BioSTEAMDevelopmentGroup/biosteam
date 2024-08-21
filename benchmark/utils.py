# -*- coding: utf-8 -*-
"""
"""
import biosteam as bst
import thermosteam as tmo

def high_precision(f):
    def g(*args, **kwargs):
        m = (
            bst.MultiStageEquilibrium.default_maxiter,
            bst.MultiStageEquilibrium.default_max_attempts,
            bst.MultiStageEquilibrium.default_fallback_maxiter,
            bst.MultiStageEquilibrium.default_molar_tolerance,
            bst.MultiStageEquilibrium.default_relative_molar_tolerance,
            bst.MultiStageEquilibrium.default_molar_tolerance,
            bst.MultiStageEquilibrium.default_relative_molar_tolerance,
            bst.MultiStageMixerSettlers.default_maxiter,
            tmo.VLE.maxiter,
            tmo.VLE.T_tol,
            tmo.VLE.P_tol,
            tmo.VLE.H_hat_tol,
            tmo.VLE.S_hat_tol ,
            tmo.VLE.V_tol,
            tmo.VLE.x_tol,
            tmo.VLE.y_tol,
            tmo.LLE.shgo_options,
            tmo.LLE.differential_evolution_options,
            tmo.LLE.pseudo_equilibrium_outer_loop_options,
            tmo.LLE.pseudo_equilibrium_inner_loop_options,
            tmo.LLE.default_composition_cache_tolerance,
            tmo.LLE.default_temperature_cache_tolerance,
        )
        bst.MultiStageEquilibrium.default_maxiter = 10
        bst.MultiStageEquilibrium.default_max_attempts = 20
        bst.MultiStageEquilibrium.default_fallback_maxiter = 1
        bst.MultiStageEquilibrium.default_molar_tolerance = 1e-6
        bst.MultiStageEquilibrium.default_relative_molar_tolerance = 1e-9
        bst.MultiStageMixerSettlers.default_maxiter = 50
        tmo.VLE.maxiter = 50 # -> 20 [-]
        tmo.VLE.T_tol = 1e-9 # -> 5e-8 [K]
        tmo.VLE.P_tol = 1e-3 # -> 1. [Pa]
        tmo.VLE.H_hat_tol = 1e-9 # -> 1e-6 [J/g]
        tmo.VLE.S_hat_tol = 1e-9 # -> 1e-6 [J/g/K]
        tmo.VLE.V_tol = 1e-9 # -> 1e-6 [mol %]
        tmo.VLE.x_tol = 1e-12 # -> 1e-9 [mol %]
        tmo.VLE.y_tol = 1e-12 # -> 1e-9 [mol %]
        tmo.LLE.shgo_options = dict(f_tol=1e-9, minimizer_kwargs=dict(f_tol=1e-9))
        tmo.LLE.differential_evolution_options = {'seed': 0, 'popsize': 12, 'tol': 1e-9}
        tmo.LLE.pseudo_equilibrium_outer_loop_options = dict(
            xtol=1e-12, maxiter=200, checkiter=False, 
            checkconvergence=False, convergenceiter=20,
        )
        tmo.LLE.pseudo_equilibrium_inner_loop_options = dict(
            xtol=1e-16, maxiter=200, checkiter=False,
            checkconvergence=False, convergenceiter=20,
        )
        tmo.LLE.default_composition_cache_tolerance = 1e-12
        tmo.LLE.default_temperature_cache_tolerance = 1e-9
        try:
            return f(*args, **kwargs)
        finally:
            (bst.MultiStageEquilibrium.default_maxiter,
             bst.MultiStageEquilibrium.default_max_attempts,
             bst.MultiStageEquilibrium.default_fallback_maxiter,
             bst.MultiStageEquilibrium.default_molar_tolerance,
             bst.MultiStageEquilibrium.default_relative_molar_tolerance,
             bst.MultiStageEquilibrium.default_molar_tolerance,
             bst.MultiStageEquilibrium.default_relative_molar_tolerance,
             bst.MultiStageMixerSettlers.default_maxiter,
             tmo.VLE.maxiter,
             tmo.VLE.T_tol,
             tmo.VLE.P_tol,
             tmo.VLE.H_hat_tol,
             tmo.VLE.S_hat_tol ,
             tmo.VLE.V_tol,
             tmo.VLE.x_tol,
             tmo.VLE.y_tol,
             tmo.LLE.shgo_options,
             tmo.LLE.differential_evolution_options,
             tmo.LLE.pseudo_equilibrium_outer_loop_options,
             tmo.LLE.pseudo_equilibrium_inner_loop_options,
             tmo.LLE.default_composition_cache_tolerance,
             tmo.LLE.default_temperature_cache_tolerance) = m
    return g