# -*- coding: utf-8 -*-
"""
Created on Sat Mar 16 13:38:11 2024

@author: cortespea
"""
from matplotlib import pyplot as plt
from numpy.testing import assert_allclose
import biosteam as bst
import numpy as np
import pandas as pd
from biosteam.utils import colors as c
from colorpalette import Color
from scipy import interpolate
from scipy.ndimage import gaussian_filter
from math import sqrt
from warnings import catch_warnings, filterwarnings
from typing import NamedTuple
from . import systems
import os
import pickle
import thermosteam as tmo
from scipy.integrate import romb
from scipy.optimize import minimize
from warnings import warn
from biosteam import report
import openpyxl
import json
from thermosteam.phenomenode import colors as phenomena_colors
import matplotlib.animation as animation
import matplotlib.ticker as ticker

__all__ = (
    'get_variable_profiles',
    'plot_profile',
    'register',
    'all_systems',
    'test_convergence',
    # 'save_stream_tables_and_specifications',
    'create_graph_theoretic_convergence_gif',
    'create_convergence_plot',
    'get_phenomegraph',
    'create_graphs',
    'create_convergence_gif',
    'plot_oscillating_variables',
    'get_LLE_partition_coefficients',
)

# %% Make simulation rigorous to achieve lower tolerances

bst.System.strict_convergence = False
bst.System.default_maxiter = 200
bst.System.default_molar_tolerance = 1e-6
bst.System.default_relative_molar_tolerance = 1e-6
bst.MultiStageEquilibrium.default_maxiter = 20
bst.MultiStageEquilibrium.max_attempts = 0 # 0 for tracking
bst.MultiStageEquilibrium.default_S_tolerance = 1e-9
bst.MultiStageEquilibrium.default_relative_S_tolerance = 1e-12
bst.MultiStageMixerSettlers.default_maxiter = 30 # 50
bst.MultiStageMixerSettlers.max_attempts = 0
bst.PhasePartition.B_relaxation_factor = 0.5
# bst.MultiStageEquilibrium.default_fallback = ('sequential', None)
# bst.MultiStageMixerSettlers.default_fallback = None
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
    xtol=1e-16, rtol=1e-16, maxiter=200, checkiter=False, 
    checkconvergence=False, convergenceiter=20,
)
tmo.LLE.pseudo_equilibrium_inner_loop_options = dict(
    xtol=1e-12, rtol=1e-16, maxiter=200, checkiter=False,
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

# %% System creation

try:
    images_folder = os.path.join(os.path.dirname(__file__), 'images')
    simulations_folder = os.path.join(os.path.dirname(__file__), 'simulations')
    graphs_folder = os.path.join(os.path.dirname(__file__), 'graphs')
except:
    images_folder = os.path.join(os.getcwd(), 'images')
    simulations_folder = os.path.join(os.getcwd(), 'simulations')
    graphs_folder = os.path.join(os.getcwd(), 'graphs')
    
stages_file = os.path.join(simulations_folder, 'system_stages.json')
all_systems = {}
system_titles = {}
system_convergence_times = {}
system_tickmarks = {} # Key: closeup-bool, name-str
system_labels = {}
system_yticks = {}
units = set()
try:
    with open(stages_file) as f: system_stages = json.load(f)
except:
    system_stages = {}

def register(name, title, time, tickmarks, label, 
             yticks=None, 
             ideal_tickmarks=None, 
             closeup_ideal_tickmarks=None,
             closeup_nonideal_tickmarks=None,
             f=None,
             unit=False):
    if f is None: f = getattr(systems, 'create_system_' + name, None) or getattr(systems, 'create_' + name + '_system')
    if yticks is None: yticks = [(-10, -5, 0, 5), (-10, -5, 0, 5)]
    if name not in system_stages: 
        sys = f(None)
        system_stages[name] = len(sys.stages)
    create_system = f
    name_fs = [
        (name, f),
        (name + '_ideal', lambda *args, **kwargs: create_system(*args, ideal=True, **kwargs)),
    ]
    for i, f in name_fs:
        all_systems[i] = f
        system_titles[i] = title
        system_convergence_times[i] = time
        system_labels[i] = label
        system_yticks[i] = yticks
        if unit: units.add(i)
    if ideal_tickmarks is None: ideal_tickmarks = tickmarks
    if closeup_ideal_tickmarks is None: closeup_ideal_tickmarks = ideal_tickmarks
    if closeup_nonideal_tickmarks is None: closeup_nonideal_tickmarks = tickmarks
    system_tickmarks[False, name] = tickmarks
    system_tickmarks[False, name + '_ideal'] = ideal_tickmarks
    system_tickmarks[True, name] = closeup_nonideal_tickmarks
    system_tickmarks[True, name + '_ideal'] = closeup_ideal_tickmarks

def get_LLE_partition_coefficients(sys):
    partition_coefficients = {}
    with catch_warnings():
        filterwarnings('ignore')
        f_sys = all_systems[sys]
        bst.F.set_flowsheet(sys)
        po = f_sys('Sequential modular')
        po.set_tolerance(rmol=1e-6, mol=1e-6, subsystems=True, method='fixed-point', 
                         maxiter=200)
        po.simulate(design_and_cost=False)
        for i in po.units:
            if getattr(i, 'phases', None) == ('L', 'l'):
                if hasattr(i, 'stages'):
                    partition_coefficients[i.ID] = np.array([i.K for i in i.stages])
                else:
                    partition_coefficients[i.ID] = i.K
    return partition_coefficients

register(
    'acetic_acid_complex', 'Rigorous system',
    90, [0, 17, 34, 51, 68], 'AcOH\nindustrial\ndewatering', 
    [(-15, -10, -5, 0, 5), (-15, -10, -5, 0, 5)],
    ideal_tickmarks=[0, 6, 12, 18, 24],
    closeup_ideal_tickmarks=[0, 0.5, 1, 1.5, 2, 2.5],
    closeup_nonideal_tickmarks=[0, 5, 10, 15, 20],
    # [(-5, -2.5, 0, 2.5, 5), (-8, -5, -2, 1, 4)],
)
register(
    'acetic_acid_simple', 'Subsystem',
    10, 
    [0, 1.5, 3, 4.5, 6], 
    'AcOH\npartial\ndewatering',
    [(-15, -10, -5, 0, 5, 10), (-15, -10, -5, 0, 5, 10)],
    ideal_tickmarks=[0, 0.4, 0.8, 1.2, 1.6],
    closeup_ideal_tickmarks=[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6],
    closeup_nonideal_tickmarks=[0, 0.3, 0.6, 0.9, 1.2, 1.5],
    # [(-15, -10, -5, 0, 5, 10), (-15, -10, -5, 0, 5, 10)],
)
register(
    'acetic_acid_complex_decoupled', 'Shortcut system',
    6, [0, 0.2, 0.4, 0.6, 0.8, 1],'AcOH\nshortcut\ndewatering',
    [(-15, -10, -5, 0, 5), (-15, -10, -5, 0, 5)],
    ideal_tickmarks=[0.1, 0.2, 0.3, 0.4],
)
# register(
#     'alcohol_narrow_flash', 'Alcohol flash narrow',
#     0.05, [0.01, 0.02, 0.03, 0.04, 0.05], 'Alcohol\nflash\nnarrow'
# )
# register(
#     'alcohol_wide_flash', 'Alcohol flash wide',
#     0.05, [0.01, 0.02, 0.03, 0.04, 0.05], 'Alcohol\nflash\nwide'
# )

register(
    'butanol_purification', 'Butanol purification',
    6, [0, 0.1, 0.2, 0.3, 0.4], 'BtOH\nseparation',
    [(-15, -10, -5, 0, 5), (-15, -10, -5, 0, 5)],
)
register(
    'ethanol_purification', 'Ethanol purification',
    5, [0, 0.03, 0.06, 0.09, 0.12], 'EtOH\nseparation',
    [(-15, -10, -5, 0, 5), (-15, -10, -5, 0, 5)],
)
# register(
#     'hydrocarbon_narrow_flash', 'Hydrocarbon flash narrow',
#     0.05, [0.01, 0.02, 0.03, 0.04, 0.05], 'Alkane\nflash\nnarrow'
# )
# register(
#     'hydrocarbon_wide_flash', 'Hydrocarbon flash wide',
#     0.1, [0.02, 0.04, 0.06, 0.08, 0.10], 'Alkane\nflash\nwide'
# )
register(
    'haber_bosch_process', 'Haber-Bosch',
    1, [0, 0.006, 0.012, 0.018, 0.024], 'Haber-Bosch\nammonia\nproduction',
    [[-10, -7.5, -5, -2.5, 0], [-10, -7.5, -5, -2.5, 0]],
    closeup_nonideal_tickmarks=[0, 0.005, 0.01, 0.015, 0.02],
)

register(
    'stripper', 'Stripper',
    1, [0, 0.1, 0.2, 0.3, 0.4, 0.5], 'Stripper',
    [(-15, -10, -5, 0, 5), (-15, -10, -5, 0, 5)],
    unit=True
)

register(
    'flash', 'Flash',
    0.05, [0, 0.01, 0.02, 0.03, 0.04, 0.05], 'Flash',
    [(-15, -10, -5, 0, 5), (-15, -10, -5, 0, 5)],
    unit=True
)

with open(stages_file, 'w') as file: json.dump(system_stages, file)

# %% Testing

def test_convergence(systems=None, algs=None, maxiter=None):
    if systems is None: systems = list(all_systems)
    elif isinstance(systems, str): systems = [systems]
    outs = []
    if algs is None: algs = ('po', 'sm')
    names = {
        'po': 'Phenomena based',
        'pm': 'phenomena modular',
        'sm': 'sequential modular',
    }
    with catch_warnings():
        filterwarnings('ignore')
        time = bst.Timer()
        for sys in systems:
            f_sys = all_systems[sys]
            new = []
            print(sys)
            for alg in algs:
                name = names[alg]
                bst.F.set_flowsheet(alg)
                po = f_sys(name)
                po.track_convergence()
                if alg == 'po':
                    maxiter = 500
                else:
                    maxiter = 100
                po.set_tolerance(rmol=1e-6, mol=1e-6, subsystems=True, method='fixed-point', 
                                 maxiter=maxiter if maxiter is None else maxiter)
                time.start()
                po.simulate(design_and_cost=False)
                print(f'- {name.capitalize()}', time.measure())
                po.show()
                new.append(po)
            outs.append(new)
            left, right = new
            for s_sm, s_dp in zip(left.streams, right.streams):
                # assert s_sm.node_tag == s_dp.node_tag
                actual = s_sm.mol
                value = s_dp.mol
                try:
                    assert_allclose(actual, value, rtol=0.01, atol=0.01)
                except:
                    if s_sm.source:
                        print(f'- Multiple steady stages for {s_sm}, {s_sm.source}-{s_sm.source.outs.index(s_sm)}: L-{actual} R-{value}')
                    else:
                        print(f'- Multiple steady stages for {s_sm}, {s_sm.sink.ins.index(s_sm)}-{s_sm.sink}: L-{actual} R-{value}')
    return outs

# %% Stream tables and specifications

# def save_stream_tables_and_specifications(systems=None):
#     if systems is None: systems = list(all_systems)
#     for sys in systems:
#         po = Tracker(sys, 'Phenomena based')
#         for i in range(50): po.run()
#         with po.system.stage_configuration() as conf:
#             for i in conf.streams:
#                 if i.ID == '': i.ID = ''
#         name = system_labels[sys].replace('\n', ' ')
#         file_name = f'{name}.xlsx'
#         file = os.path.join(simulations_folder, file_name)
#         po.system.save_report(
#             file=file,
#             sheets={
#                 'Flowsheet',
#                 'Stream table',
#                 'Specifications',
#                 'Reactions',
#             },
#             stage=True,
#         )


# %% Convergence time

def convergence_time(sm, po):
    ysm = np.array(sm)
    ypo = np.array(po)
    cutoff = ysm.min() + 1
    sm_index = sum(ysm > cutoff)
    po_index = sum(ypo > cutoff)
    return sm['Time'][sm_index], po['Time'][po_index]


# %% Graph-theoretic representation of convergence

def get_phenomegraph(
        system=None, algorithm=None, load=True, save=True,
        simulate=True, bipartite=True, new_graph=False, mock=True,
        convergence_rate=25,
    ):
    bst.MultiStageEquilibrium.max_attempts = 0 # For consistent tracking
    if system is None: system = 'acetic_acid_simple'
    if algorithm is None: algorithm = 'Phenomena based'
    alg = algorithm.replace(' ', '_').replace('-', '_')
    file = os.path.join(graphs_folder, f"{system}_{alg}")
    if not bipartite: file += '_monopartite'
    if new_graph:
        dotfile = None
    else:
        dotfile = tmo.DotFile(os.path.join(graphs_folder, system))
    if load:
        try:
            with open(file, 'rb') as f: phenomena_graph = pickle.load(f)
        except Exception as e:
            print(e)
            return get_phenomegraph(
                system=system, algorithm=algorithm,
                load=False, save=save, bipartite=bipartite,
                new_graph=new_graph, simulate=simulate,
            )
    else:
        bst.F.clear()
        sys = all_systems[system](algorithm)
        if mock:
            sys.mock_convergence(convergence_rate=convergence_rate)
        else:
            sys.set_tolerance(
                mol=1e-2,
                rmol=1e-2,
                maxiter=500,
                method='fixed-point',
                subsystems=True,
            )
            if simulate:
                sys.track_convergence(
                    True, 
                )
                sys.simulate()
        phenomena_graph = sys.get_phenomegraph(
            dotfile=dotfile, bipartite=bipartite,
        )
        if save:
            if not simulate: file += '_nosim'
            sys.track_convergence(False)
            with open(file, 'wb') as f: pickle.dump(phenomena_graph, f)
    return phenomena_graph

def create_graphs(system=None, bipartite=True, algorithms=None):
    if system is None: system = 'acetic_acid_simple'
    if algorithms is None: 
        if system in units:
            algorithms = ('Phenomena based',)
        else:
            algorithms = ('sequential modular', 'Phenomena based')
    elif isinstance(algorithms, str): algorithms = [algorithms]
    header = system
    if not bipartite: header += '_monopartite'
    first = True
    for algorithm in algorithms:
        graph = get_phenomegraph(system, algorithm=algorithm, simulate=False, new_graph=False)
        if first:
            file = os.path.join(graphs_folder, f"{header}_graph.png")
            graph.write(file)
        for subgraph in graph.subgraphs:
            _, name = subgraph.name.split('.')
            file = os.path.join(graphs_folder, f"{header}_subgraph_{name}.png")
            subgraph.write(file)
        
def create_convergence_plot(
        system=None, algorithm=None, load=True, save=True,
    ):
    if system is None: system = 'acetic_acid_simple'
    if algorithm is None: algorithm = 'Phenomena based'
    alg = algorithm.replace(' ', '_').replace('-', '_')
    phenomena_graph = get_phenomegraph(
        system=system, algorithm=algorithm, 
        load=load, save=save, 
    )
    file = os.path.join(graphs_folder, f"{system}_{alg}")
    return phenomena_graph.plot_convergence_profile(file=file)

def create_graph_theoretic_convergence_gif(
        system=None, algorithm=None, load=True, save=True, 
        fps=3, new_graph=False, mock=True, convergence_rate=25,
    ):
    if system is None: system = 'acetic_acid_simple'
    if algorithm is None: algorithm = 'Phenomena based'
    alg = algorithm.replace(' ', '_').replace('-', '_')
    bipartite_graph = get_phenomegraph(
        system=system, algorithm=algorithm, new_graph=new_graph,
        load=load, save=True, bipartite=True, mock=mock, convergence_rate=convergence_rate,
    )
    aggregated_graph = get_phenomegraph(
        system=system, algorithm=algorithm, new_graph=new_graph,
        load=load, save=save, bipartite=False, mock=mock, convergence_rate=convergence_rate,
    )
    reference = {eq.name: n for n, nodes in enumerate(aggregated_graph.phenomenodes) for eq in nodes.equations}
    file = os.path.join(graphs_folder, f"{system}_{alg}_graph")
    
    bipartite_graph.convergence_gif(
        time=aggregated_graph.time,
        profiles=aggregated_graph.profiles,
        reference=reference,
        file=file, fps=fps, 
    )
    for i, pg in enumerate(bipartite_graph.subgraphs):
        pg.convergence_gif(file=f'{file}_{i}', fps=fps, 
                           profiles=aggregated_graph.profiles,
                           time=aggregated_graph.time,
                           reference=reference)
    
def create_convergence_gif(
        system=None, algorithm=None, load=False, 
        fps=3, new_graph=False,
        mock=True, convergence_rate=25,
    ):
    plt.style.use('dark_background')
    plt.rcParams.update({
        "figure.facecolor":  (0.0, 0.0, 0.0, 0),  # red   with alpha = 30%
        "axes.facecolor":    (0.0, 0.0, 0.0, 1),  # green with alpha = 50%
        "savefig.facecolor": (0.0, 0.0, 0.0, 0),  # blue  with alpha = 20%
    })
    fs = 9
    bst.set_font(fs)
    n_rows = 1
    n_cols = 1
    bst.set_figure_size(aspect_ratio=1.1 / n_cols, width='half')
    fig, ax = plt.subplots(n_rows, n_cols)
    profiles = get_variable_profiles(system, algorithm, number=100, load=load, mock=True, convergence_rate=convergence_rate)
    time = profiles['time']
    time = time - time[0] 
    categories = [i for i in profiles.keys() if i != 'time']
    yticks = [0, 20, 40, 60, 80, 100]
    vline, = plt.plot([0, 0], [0, 102], color='#E7E6E6')
    for category in categories:
        y = profiles[category]
        y -= y.min()
        y *= 100 / y.max()
        color = phenomena_colors[category.split('_')[0]]
        plt.sca(ax)
        plt.step(time, y, '-', color=color, lw=2, where='post', alpha=1)
        # plt.step(time, y, lw=0, marker='.', color=color, markersize=1, where='post')
    plt.xlim(time[0], time[-1])
    plt.ylim(yticks[0], yticks[-1] + 2)
    ax.xaxis.set_major_locator(ticker.MaxNLocator(4, integer=True))
    plt.ylabel('MSE [%]')
    plt.xlabel('Iteration')
    bst.utils.style_axis(
        ax, yticks=yticks, 
    )
    left = 0.1
    top = 0.85
    if n_rows == 2:
        left = 0.25
        bottom = 0.1
    elif n_rows == 1:
        bottom = 0.15
        if n_cols == 1:
            left = 0.25
            top = 0.90
            bottom = 0.2
    else:
        bottom = 0.08
    plt.subplots_adjust(right=0.96, left=left, bottom=bottom, top=top, hspace=0, wspace=0)
    N_frames = max(5 * len(time), 200)
    scale = len(time) / N_frames
    def update(frame):
        x = frame * scale
        vline.set_data([x, x], [0, 102])
        return (vline,)
    
    fps /= scale
    ani = animation.FuncAnimation(
        fig=fig, 
        func=update, 
        frames=N_frames, 
        interval=fps,
    )
    plt.show()
    file = os.path.join(images_folder, f"{system}_{algorithm}")
    duration = [1000 / fps] * N_frames
    duration[-1] += 5000 
    print(sum(duration) / 1000)
    writter = animation.PillowWriter(fps=fps)
    writter.dpi = 300
    def finish():
        writter._frames[0].save(
            writter.outfile, format='GIF', save_all=True, append_images=writter._frames[1:],
            duration=duration, loop=1, optimize=False,
        )
    writter.finish = finish
    ani.save(
        filename=file + '.gif', 
        writer=writter, 
    )
    return ani
    
    
    # for i in ('svg', 'png'):
    #     name = f'profile.{i}'
    #     file = os.path.join(images_folder, name)
    #     plt.savefig(file, dpi=900, transparent=True)
    # for i in ('svg', 'png'):
    #     system_names = '_'.join(systems)
    #     name = f'{system_names}_profile.{i}'
    #     file = os.path.join(images_folder, name)
    #     plt.savefig(file, dpi=900, transparent=True)


# %% Profiling and benchmarking utilities

categories = (
    'material', 'energy', 'vle', 'lle',
    'reaction', 'compression', 'shortcut'
)

class ErrorPoint(NamedTuple):
    time: float; error: float

default_steady_state_cutoff = 1

def steady_state_error(profile, steady_state_cutoff=None):
    if steady_state_cutoff is None: steady_state_cutoff = default_steady_state_cutoff
    minimum_error = np.log10(10 ** profile['Component flow rate error'][-1] + 10 ** profile['Temperature error'][-1]) 
    return minimum_error + steady_state_cutoff
    
def benchmark(profile, steady_state_error=None):
    if steady_state_error is None: steady_state_error = steady_state_error(profile)
    time = profile['Time']
    error = np.log10(10 ** profile['Component flow rate error'] + 10 ** profile['Temperature error'])
    if np.isnan(error).any(): breakpoint()
    time = interpolate.interp1d(error, time, bounds_error=False)(steady_state_error)
    if np.isnan(error).any(): breakpoint()
    return time 

def create_tracked_system(name, algorithm, time=None):
    bst.F.set_flowsheet('temp', new=True)
    sys = all_systems[name](algorithm)
    sys._setup_units()
    sys.track_convergence()
    sys.set_tolerance(
        maxtime=time or system_convergence_times[name], 
        method='fixed-point',
        mol=0, rmol=0,
        T=0, rT=0,
        maxiter=int(1e16),
    )
    return sys
    
def get_steady_state(name, load=True):
    file_name = f"{name}_sequential modular_0"
    file = os.path.join(simulations_folder, file_name)
    bst.F.clear()
    if load:
        try:
            with open(file, 'rb') as f: 
                time, variable_profiles = pickle.load(f)
        except: 
            return get_steady_state(name, False)
    else:
        sys = create_tracked_system(name, 'sequential modular') 
        sys.simulate(design_and_cost=False)
        sys.assert_tracking_ok()
        sys.flowsheet.clear()
        results = time, variable_profiles = sys.get_variable_profiles()
        with open(file, 'wb') as f: pickle.dump(results, f)
    steady_state = variable_profiles.iloc[-1]
    return steady_state

def plot_oscillating_variables(name, algorithm, time=None): 
    bst.F.clear()
    sys = create_tracked_system(name, algorithm, time) 
    sys._setup()
    sys.converge()
    sys.assert_tracking_ok()
    sys.flowsheet.clear()
    time, variable_profiles = sys.get_variable_profiles(time=True)
    variable_profiles.index = time
    oscillating_columns = []
    past_names = set()
    for eq in sys.equation_nodes: 
        names = [i.name for i in eq.outputs]
        eq_name = eq.name
        if 'phenomena' not in eq_name: continue
        for name in names:
            if name in past_names: continue
            e = variable_profiles[name]
            def add(name, e):
                mean = e.mean()
                e0 = e.iloc[0]
                N = len(e)
                crossed_mean = 0
                for i in range(1, N-1):
                    e1 = e.iloc[i]
                    crossed_mean += (
                        e1 >= mean and e0 < mean 
                        or
                        e1 <= mean and e0 > mean 
                    )
                    e1 = e0
                if crossed_mean > 5:
                    oscillating_columns.append([name, e])
            if e.ndim == 2: 
                for subname in e:
                    add([name, subname], e[subname])
            else:
                add(name, e)
        past_names.update(names)
    for name, e in oscillating_columns:
        plt.figure()
        e.plot.line()
        plt.ylabel(name)
    # oscillating_profiles = variable_profiles[oscillating_columns]
    # oscillating_profiles.index = time
    # oscillating_profiles.plot.line()
    return oscillating_columns

def get_variable_profiles(
        name, algorithm, 
        number=0, load=True, mock=False, convergence_rate=None,
        grouped_error=True,
    ): 
    file_name = f"{name}_{algorithm}_{number}"
    file = os.path.join(simulations_folder, file_name)
    bst.F.clear()
    if not mock and load:
        try:
            with open(file, 'rb') as f: time, variable_profiles = pickle.load(f)
            sys = all_systems[name](algorithm)
            sys._setup_units()
            sys.flowsheet.clear()
            sys.track_convergence()
        except:
            return get_variable_profiles(name, algorithm, number, False, mock)
    else:
        if mock:
            sys = all_systems[name](algorithm)
            sys.mock_convergence(convergence_rate=convergence_rate or 20)
            sys.flowsheet.clear()
        else:
            sys = create_tracked_system(name, algorithm) 
            sys.simulate(design_and_cost=False)
            sys.assert_tracking_ok()
            sys.flowsheet.clear()
        results = time, variable_profiles = sys.get_variable_profiles()
        with open(file, 'wb') as f: pickle.dump(results, f)
    if not grouped_error: return time, variable_profiles
    if mock:
        errors = variable_profiles.copy()
    else:
        steady_state = get_steady_state(name, True)
        steady_state.values[steady_state < 1e-6] = 0
        steady_state.values[steady_state > 1e6] = 1e6
        variable_profiles.values[variable_profiles < 1e-6] = 0
        variable_profiles.values[variable_profiles > 1e6] = 1e6
        errors = (variable_profiles - steady_state)
        # breakpoint()
        # ms = [i for i in errors if 'U1' in i[0] and 'Phi' in i[0]]
        # print(errors[ms].sum())
    profiles = {}
    past_names = set()
    for eq in sys.equation_nodes: 
        names = [i.name for i in eq.tracked_outputs]
        eq_name = eq.name
        for name in names:
            if name in past_names: continue
            e = errors[name].values
            if e.ndim == 2: e = e.sum(axis=1)
            squared_error = e * e # squared error
            for category in categories:
                if category in eq_name: break
            else:
                category = 'vle_phenomena' # Just assume distillation is vle related
            if category in profiles:
                profiles[category].append(squared_error)
            else:
                profiles[category] = [squared_error]
        past_names.update(names)
    for key, lst in profiles.items(): profiles[key] = sum(lst) / len(lst)
    profiles['time'] = time
    return profiles

# %%  Profile plot

def dct_mean_profile(dcts: list[dict], categories):
    size = np.min([len(i['time']) for i in dcts])
    t = np.array([i['time'][:size] for i in dcts]).mean(axis=0)
    names = [i for i in categories if i in dcts[0]]
    mean = {i: np.array([dct[i][:size] for dct in dcts]).mean(axis=0) for i in names}
    mean['time'] = t
    for i in names:
        values = mean[i]
        mask = np.isnan(values)
        values[mask] = values[~mask].max() # Maximum error
        values[values < 1e-15] = 1e-15
        values[:] = np.log10(values)
    return mean

def dct_mean_std(dcts: list[dict], keys: list[str]):
    n = len(dcts)
    values = {i: np.zeros(n) for i in keys}
    for i, dct in enumerate(dcts):
        for key in keys: values[key][i] = dct[key]
    return {i: (values[i].mean(), values[i].std()) for i in keys}

def plot_profile(
        systems=None, algorithms=None, N=1, load=True, save=True, dark=False,
        closeup=False, fs=None, width=None, aspect_ratio=None, label=True,
    ):
    if dark: 
        plt.style.use('dark_background')
        plt.rcParams.update({
            "figure.facecolor":  (0.0, 0.0, 0.0, 0),  # red   with alpha = 30%
            "axes.facecolor":    (0.0, 0.0, 0.0, 1),  # green with alpha = 50%
            "savefig.facecolor": (0.0, 0.0, 0.0, 0),  # blue  with alpha = 20%
        })
    if algorithms is None: algorithms = ('Phenomena based', 'phenomena modular', 'sequential modular', )
    if isinstance(systems, str): systems = [systems]
    if systems is None: systems = list(all_systems)
    if fs is None: fs = 10
    bst.set_font(fs)
    n_rows = len(algorithms)
    n_cols = len(systems)
    if width is None: width = 'full' if n_cols > 1 else 'half'
    if aspect_ratio is None: aspect_ratio = 0.325 * n_rows if n_rows > 1 else 0.65
    bst.set_figure_size(aspect_ratio=aspect_ratio, width=width)
    fig, all_axes = plt.subplots(n_rows, n_cols)
    all_axes = np.reshape(np.array(all_axes), [n_rows, n_cols])
    for m, sys in enumerate(systems):
        time = system_convergence_times[sys]
        axes = all_axes[:, m]
        profiles = {
            alg: [get_variable_profiles(sys, alg, i, load=load) for i in range(N)]
            for alg in algorithms
        }
        mean_profiles = {i: dct_mean_profile(j, categories) for i, j in profiles.items()}
        for dct in mean_profiles.values():
            time = dct['time']
            time[:] -= time[0]
        # yticks_list = system_yticks[sys]
        for n, (alg, mean_profile) in enumerate(mean_profiles.items()):
            ax = axes[n]
            plt.sca(ax)
            xticks = system_tickmarks[closeup, sys]
            for category in categories:
                # if alg == 'sequential modular' and category == 'lle': breakpoint()
                if category not in mean_profile: continue
                y = np.array(mean_profile[category])
                # y = gaussian_filter(y, 0.2)
                t = np.array(mean_profile['time'])
                c = phenomena_colors[category.split('_')[0]]
                # print(alg, category)
                # print(y)
                if alg == 'Phenomena based' and sys == 'acetic_acid_complex':
                    time_, variable_profiles = get_variable_profiles(
                        'acetic_acid_complex', 'Phenomena based',
                        grouped_error=False,
                    )
                    single_phase = variable_profiles['settler.K', 'Water'] == 1
                    time_[:] -= time_[0]
                    time_ = time_[single_phase]
                    for i in time_: plt.axvline(i, c='darkgray', lw=0.01)
                plt.step(t, y, '-', color=c, lw=1, alpha=1, where='post')
            
            yticklabels = m == 0
            # yticks = yticks_list[n]
            yticks = [-12, -6, 0, 6, 12]
            if yticklabels: yticklabels = [str(i) for i in yticks]
            # if yticklabels: yticklabels = [r'$\mathrm{10}^{' f'{i}' '}$' for i in yticks]
            xticklabels = xtick0 = n == n_rows-1
            xtickf = True
            ytick0 = n == n_rows-1
            ytickf = n == 0
            plt.xlim(0, xticks[-1])
            plt.ylim(yticks[0], yticks[-1])
            # if m == 0: plt.ylabel(f'{alg}\nlog$_{{10}}$ error')
            bst.utils.style_axis(
                ax, xticks=xticks, yticks=yticks, 
                xtick0=xtick0, xtickf=xtickf, ytick0=ytick0, ytickf=ytickf,
                xticklabels=xticklabels,
                yticklabels=yticklabels,
            )
    if label:
        for axis in all_axes[-1]:
            plt.sca(axis)
            plt.xlabel('Time [s]')
        if n_rows == 3:
            central_axis = all_axes[1, 0]
            plt.sca(central_axis)
            plt.ylabel('log$_{10}$Error')
        else:
            for axis in all_axes[:, 0]:
                plt.sca(axis)
                plt.ylabel('log$_{10}$Error')
    left = 0.2
    top = 0.85
    if n_rows == 2:
        left = 0.25
        bottom = 0.1
    elif n_rows == 1:
        bottom = 0.15
        if n_cols == 1:
            left = 0.25
            top = 0.90
            bottom = 0.2
    else:
        bottom = 0.15
    plt.subplots_adjust(right=0.96, left=left, bottom=bottom, top=top, hspace=0, wspace=0.1)
    for i in ('svg', 'png'):
        name = f'profile.{i}'
        file = os.path.join(images_folder, name)
        plt.savefig(file, dpi=900, transparent=True)
    for i in ('svg', 'png'):
        system_names = '_'.join(systems)
        name = f'{system_names}_profile.{i}'
        file = os.path.join(images_folder, name)
        plt.savefig(file, dpi=900, transparent=True)
    # return fig, all_axes, sms, pos
