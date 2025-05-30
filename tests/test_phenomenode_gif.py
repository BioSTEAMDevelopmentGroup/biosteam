# -*- coding: utf-8 -*-
"""
Created on Fri May 16 12:17:20 2025

@author: yoelr
"""

# %% Make simulation rigorous to achieve lower tolerances
import biosteam as bst
import thermosteam as tmo
bst.System.default_maxiter = 200
bst.System.default_molar_tolerance = 1e-6
bst.System.default_relative_molar_tolerance = 1e-6
bst.MultiStageEquilibrium.default_maxiter = 20
bst.MultiStageEquilibrium.default_max_attempts = 10
bst.MultiStageEquilibrium.default_molar_tolerance = 1e-9
bst.MultiStageEquilibrium.default_relative_molar_tolerance = 1e-12
bst.MultiStageMixerSettlers.default_maxiter = 50
tmo.BubblePoint.maxiter = 100 # -> 50 [-]
tmo.DewPoint.maxiter = 100 # -> 50 [-]
tmo.BubblePoint.T_tol = 1e-16 # -> 1e-9 [K]
tmo.DewPoint.T_tol = 1e-16 # -> 1e-9 [K]
tmo.BubblePoint.P_tol = 1e-9 # -> 1e-3 [Pa]
tmo.DewPoint.P_tol = 1e-9 # -> 1e-3 [Pa]
tmo.VLE.T_tol = 1e-12 # -> 5e-8 [K]
tmo.VLE.P_tol = 1e-12 # -> 1. [Pa]
tmo.VLE.maxiter = 50 # -> 20 [-]
tmo.VLE.T_tol = 1e-12 # -> 5e-8 [K]
tmo.VLE.P_tol = 1e-9 # -> 1. [Pa]
tmo.VLE.H_hat_tol = 1e-12 # -> 1e-6 [J/g]
tmo.VLE.S_hat_tol = 1e-12 # -> 1e-6 [J/g/K]
tmo.VLE.V_tol = 1e-12 # -> 1e-6 [mol %]
tmo.VLE.x_tol = 1e-12 # -> 1e-9 [mol %]
tmo.VLE.y_tol = 1e-12 # -> 1e-9 [mol %]
tmo.LLE.shgo_options = dict(f_tol=1e-9, minimizer_kwargs=dict(f_tol=1e-9))
tmo.LLE.differential_evolution_options = {'seed': 0, 'popsize': 12, 'tol': 1e-9}
tmo.LLE.pseudo_equilibrium_outer_loop_options = dict(
    xtol=1e-16, maxiter=200, checkiter=False, 
    checkconvergence=False, convergenceiter=20,
)
tmo.LLE.pseudo_equilibrium_inner_loop_options = dict(
    xtol=1e-12, maxiter=200, checkiter=False,
    checkconvergence=False, convergenceiter=20,
)
tmo.LLE.default_composition_cache_tolerance = 1e-16
tmo.LLE.default_temperature_cache_tolerance = 1e-16
bst.ShortcutColumn.iter_solver_kwargs = dict(
    xtol=1e-16,
    checkiter=False,
    checkconvergence=False, 
    convergenceiter=20,
)

# %%

import numpy as np
from numpy.testing import assert_allclose
solvent_feed_ratio = 1
thermo = bst.Thermo(['Water', 'AceticAcid', 'EthylAcetate'], cache=True)
bst.settings.set_thermo(thermo)
with bst.System(algorithm='phenomena oriented',
                molar_tolerance=1e-6,
                relative_molar_tolerance=1e-6,
                maxiter=100,
                method='fixed-point') as sys:
    feed = bst.Stream('feed', AceticAcid=6660, Water=43600)
    solvent = bst.Stream('solvent', EthylAcetate=65000)
    recycle = bst.Stream('recycle')
    LE = bst.MultiStageEquilibrium(
        N_stages=6, ins=[feed, solvent, recycle],
        feed_stages=(0, -1, -1),
        phases=('L', 'l'),
        maxiter=200,
        use_cache=True,
        method='fixed-point',
    )
    # DAA = bst.MultiStageEquilibrium(N_stages=6, ins=[LE-0], feed_stages=[3],
    #     outs=['vapor', 'liquid'],
    #     stage_specifications={0: ('Reflux', 0.673), -1: ('Boilup', 2.57)},
    #     maxiter=200,
    #     phases=('g', 'l'),
    #     method='fixed-point',
    #     use_cache=True,
    # )
    DEA = bst.MultiStageEquilibrium(N_stages=6, ins=[LE-1], feed_stages=[3],
        outs=['vapor', 'liquid'],
        stage_specifications={0: ('Reflux', 0.673), -1: ('Boilup', 2.57)},
        phases=('g', 'l'),
        maxiter=200,
        method='fixed-point',
        use_cache=True,
    )
    HX = bst.SinglePhaseStage(ins=DEA-0, outs=recycle, T=320, phase='l')
    chemicals = bst.settings.chemicals

    @LE.add_specification(run=True)
    def fresh_solvent_flow_rate():
        broth = feed.F_mass
        EtAc_recycle = recycle.imass['EthylAcetate']
        EtAc_required = broth * solvent_feed_ratio
        if EtAc_required < EtAc_recycle:
            recycle.F_mass *= EtAc_required / EtAc_recycle
            EtAc_recycle = recycle.imass['EthylAcetate']
        EtAc_fresh = EtAc_required - EtAc_recycle
        solvent.imass['EthylAcetate'] = max(
            0, EtAc_fresh
        )

    @solvent.material_balance
    def fresh_solvent_flow_rate():
        s = np.ones(chemicals.size)
        r = np.zeros(chemicals.size)
        v = r.copy()
        index = chemicals.index('EthylAcetate')
        r[index] = 1
        v[index] = solvent_feed_ratio * feed.F_mass / chemicals.EthylAcetate.MW
        return (
            {solvent: s,
             recycle: r},
             v
        )
    
sys.track_convergence(True)
sys.simulate()
phenomena_graph = sys.get_phenomena_graph()
phenomena_graph.convergence_diagram()