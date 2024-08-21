# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 00:05:26 2024

@author: yoelr
"""

# %% 
import biosteam as bst
import thermosteam as tmo
bst.MultiStageEquilibrium.default_maxiter = 10
bst.MultiStageEquilibrium.default_max_attempts = 20
bst.MultiStageEquilibrium.default_fallback_maxiter = 1
bst.MultiStageEquilibrium.default_molar_tolerance = 1e-6
bst.MultiStageEquilibrium.default_relative_molar_tolerance = 1e-9
bst.MultiStageMixerSettlers.default_maxiter = 50
tmo.BubblePoint.maxiter = 200 # -> 50 [-]
tmo.DewPoint.maxiter = 200 # -> 50 [-]
tmo.BubblePoint.T_tol = 1e-12 # -> 1e-9 [K]
tmo.DewPoint.T_tol = 1e-12 # -> 1e-9 [K]
tmo.BubblePoint.P_tol = 1e-6 # -> 1e-3 [Pa]
tmo.DewPoint.P_tol = 1e-6 # -> 1e-3 [Pa]
tmo.VLE.T_tol = 1e-9 # -> 5e-8 [K]
tmo.VLE.P_tol = 1e-3 # -> 1. [Pa]
tmo.VLE.maxiter = 100 # -> 20 [-]
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

# %%

runfile('C:/Users/yoelr/code/biosteam/benchmark/systems/acetic_acid_system.py', wdir='C:/Users/yoelr/code/biosteam/benchmark/systems')
from scipy.optimize import minimize
po = create_acetic_acid_complex_system('phenomena-oriented')
po.flatten()
po.set_tolerance(mol=1e-5, rmol=1e-5, maxiter=50)
po.simulate()
bst.units.stage.PhasePartition.S_relaxation_factor = 0
bst.units.stage.PhasePartition.B_relaxation_factor = 0
bst.units.stage.PhasePartition.K_relaxation_factor = 0
bst.units.stage.PhasePartition.T_relaxation_factor = 0

# %% Optimization

stages = []
shortcuts = []
for i in po.units:
    if hasattr(i, 'stages'):
        stages.extend(i.stages)
    elif isinstance(i, bst.StageEquilibrium):
        stages.append(i)
    elif isinstance(i, bst.Distillation):
        shortcuts.append(i)
        i._update_equilibrium_variables()
all_stages = stages + shortcuts   
N_chemicals = po.units[0].chemicals.size
Sb, safe = bst.units.stage.bottoms_stripping_factors_safe(
    np.array([i.B for i in all_stages]), 
    np.array([i.K for i in all_stages]), 
)
Sb_full = Sb
Sb_index = Sb_index = [
    i for i, j in enumerate(all_stages)
    if not (getattr(j, 'B_specification', None) == 0 or getattr(j, 'B_specification', None) == np.inf)
]
Sb = Sb[Sb_index]
N_Sb_index = len(Sb_index)
lnSb = np.log(Sb).flatten()
options = {'tol': 1e-3, 'method': 'CG', 'options':{'maxiter': 500}}
count = [0]
def lnSb_objective(lnSb):
    Sb_original = np.exp(lnSb)
    Sb = Sb_full
    Sb[Sb_index] = Sb_original.reshape([N_Sb_index, N_chemicals])
    for i in Sb_index:
        all_stages[i].B = 1 # Work around Sb = 1 / K * B
        all_stages[i].K = 1 / Sb[i]
    with po.stage_configuration(aggregated=False) as conf:
        conf._solve_material_flows(composition_sensitive=False)
    for i in stages:
        mixer = i.mixer
        partition = i.partition
        mixer.outs[0].mix_from(
            mixer.ins, energy_balance=False,
        )
        if partition.phases == ('g', 'l'):
            partition._run_decoupled_KTvle()
            partition._run_decoupled_B()
            T = partition.T
            for i in (partition.outs + i.outs): i.T = T
        else:
            partition._run_lle(update=False)
    for i in shortcuts:
        for s in i.outs:
            if s.phase == 'l': 
                bp = s.bubble_point_at_P()
                s.T = bp.T
            elif s.phase == 'g': 
                dp = s.dew_point_at_P()
                s.T = dp.T
    energy_error = lambda stage: abs(
        (sum([i.H for i in stage.outs]) - sum([i.H for i in stage.ins])) / sum([i.C for i in stage.outs])
    )
    total_energy_error = sum([energy_error(stage) for stage in stages if stage.B_specification is None and stage.T_specification is None])
    for i in shortcuts: i.run()
    Sb_new = np.array([1 / (all_stages[i].K * all_stages[i].B) for i in Sb_index]).flatten()
    splits_new = 1 / (Sb_new + 1)
    splits = 1 / (Sb_original + 1)
    diff = splits_new - splits
    total_split_error = np.abs(diff).sum()
    total_error = total_split_error + total_energy_error
    err = np.sqrt(total_error)
    count[0] += 1
    if not count[0] % 100:
        print(err)
        po.show()
    return err

result = minimize(
    lnSb_objective, 
    lnSb,
    **options,
)
# %% Saving results

def dT_error(stage):
    if all([i.isempty() for i in stage.outs]): 
        return 0
    else:
        return (
            sum([i.H for i in stage.outs]) - sum([i.H for i in stage.ins])
        ) / sum([i.C for i in stage.outs])

def streams_and_stages(sys):
    all_stages = []
    adiabatic_stages = []
    streams = []
    past_streams = set()
    for i in sys.units:
        if hasattr(i, 'stages'):
            all_stages.extend(i.stages)
            for j in i.stages:
                new_streams = [i for i in (j.outs + j.ins) if i not in past_streams]
                streams.extend(new_streams)
                past_streams.update(new_streams)
                if j.B_specification is None and j.T_specification is None:
                    adiabatic_stages.append(j)
        else:
            try:
                if i.B_specification is None and i.T_specification is None:
                    adiabatic_stages.append(j)
            except:
                pass
            all_stages.append(i)
            new_streams = [j for j in (i.outs + i.ins) if j not in past_streams]
            streams.extend(new_streams)
            past_streams.update(new_streams)
    return (streams, adiabatic_stages, all_stages)

streams, adiabatic_stages, all_stages = streams_and_stages(po)
total_energy_error = sum([dT_error(i) for i in adiabatic_stages])

# %% 
import numpy as np
flows = np.array([i.mol for i in streams])
Ts = np.array([i.T for i in streams])
np.save('flows_steady_state', flows)
np.save('temperatures_steady_state', Ts)