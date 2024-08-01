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
import os
import pickle
import thermosteam as tmo

__all__ = (
    'BenchmarkModel',
    'run_monte_carlo',
    'plot_monte_carlo',
    'plot_benchmark',
    'plot_profile',
    'register',
    'test_convergence'
)

# %% Uncertainty analysis

class BenchmarkModel:
    
    def __init__(self, system):
        model = bst.Model(None, specification=lambda: None)
        parameter = model.parameter
        tracker = Tracker(system, 'phenomena oriented')
        benchmark = tracker.benchmark()
        recycles = tracker.system.get_all_recycles()
        self.material_flows = np.array([i.mol for i in recycles])
        self.Ts = np.array([i.T for i in recycles])
        chemicals = [i.ID for i in bst.settings.chemicals]
        self.model = model
        @parameter(units='-', name="log absolute tolerance", bounds=[-9, -1])
        def set_tolerance(atol):
            atol = 10 ** atol
            # atol = 1e-12
            self.sm_tracker = Tracker(system, 'sequential modular', atol=atol)
            self.sm_system = self.sm_tracker.system
            self.sm_recycles = self.sm_system.get_all_recycles()
            self.po_tracker = Tracker(system, 'phenomena oriented', atol=atol)
            self.po_system = self.po_tracker.system
            self.po_recycles = self.po_system.get_all_recycles()
        
        # for i, (flows, T) in enumerate(zip(self.material_flows, self.Ts)):
        #     @parameter(units='K', name=f"Recycle {i}", element='temperature', bounds=[T - 15, T + 15])
        #     def set_temperature(T, index=i):
        #         return
        #         self.sm_recycles[index].T = T
        #         self.po_recycles[index].T = T
                
        #     for chemical, flow in zip(chemicals, flows):
        #         @parameter(units='-', name=chemical, element=f"Recycle {i}", bounds=[-3, 3])
        #         def set_chemical_flow(flow, index=i, chemical=chemical, flow_max=flow):
        #             flow = flow_max * 10 ** flow
        #             self.sm_recycles[index].imol[chemical] = flow
        #             self.po_recycles[index].imol[chemical] = flow
        
        @model.metric(units='%')
        def relative_simulation_time():
            recycles = self.sm_recycles
            material_flows = np.array([i.mol for i in recycles])
            Ts = np.array([i.T for i in recycles])
            sm_b = np.mean([self.sm_tracker.benchmark()['Time'] for i in range(5)])
            po_b = np.mean([self.po_tracker.benchmark()['Time'] for i in range(5)])
            T_diff = (self.Ts - Ts)
            material_diff = (self.material_flows - material_flows)
            nmd = 0.5 * np.log((material_diff * material_diff).sum())
            ntd = 0.5 * np.log((T_diff * T_diff).sum())
            self._distance = (
                 nmd + ntd
            ) / 2
            self._time_po = po_b
            self._time_sm = sm_b
            return 100 * self._time_po / self._time_sm
        
        @model.metric(units='s')
        def sequential_modular_simulation_time():
            return self._time_sm
        
        @model.metric(units='s')
        def phenomena_oriented_simulation_time():
            return self._time_po
        
        @model.metric(units='-')
        def distance():
            return self._distance

def run_monte_carlo(N=50, system='alcohol_wide_flash', autosave=True, autoload=True):
    bm = BenchmarkModel(system)
    samples = bm.model.sample(N, rule='L', seed=0)
    bm.model.load_samples(samples)
    bm.model.evaluate(
        notify=10,
        autosave=autosave, autoload=autoload,
        file=os.path.join(simulations_folder, f'{system}_MC{N}')
    )
    bm.model.table.to_excel(os.path.join(simulations_folder, f'{system}_MC{N}.xlsx'))
    bm.model.table.to_excel(os.path.join(simulations_folder, f'{system}_MC.xlsx'))

def plot_monte_carlo(system='alcohol_wide_flash'):
    file = os.path.join(simulations_folder, f'{system}_MC.xlsx')
    df = pd.read_excel(file, header=[0, 1], index_col=[0])
    plt.scatter(df.iloc[:, 0], df.iloc[:, -4])
    plt.show()

# %% System creation and testing

try:
    images_folder = os.path.join(os.path.dirname(__file__), 'images')
    simulations_folder = os.path.join(os.path.dirname(__file__), 'simulations')
except:
    images_folder = os.path.join(os.getcwd(), 'images')
    simulations_folder = os.path.join(os.getcwd(), 'simulations')

all_systems = {}
system_titles = {}
system_convergence_times = {}
system_tickmarks = {}
system_labels = {}

def _register(f, name, title, time, tickmarks, label):
    all_systems[name] = f
    system_titles[name] = title
    system_convergence_times[name] = time
    system_tickmarks[name] = tickmarks
    system_labels[name] = label
    return f

def register(name, title, time, tickmarks, label):
    return lambda f: _register(f, name, title, time, tickmarks, label)

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

def test_convergence(systems=None):
    if systems is None: systems = list(all_systems)
    with catch_warnings():
        filterwarnings('ignore')
        time = bst.TicToc()
        for sys in systems:
            f_sys = all_systems[sys]
            sm = f_sys('sequential modular')
            sm.flatten()
            sm.set_tolerance(rmol=1e-5, mol=1e-5, subsystems=True, method='fixed-point', maxiter=500)
            time.tic()
            sm.simulate()
            print(sys)
            print('- Sequential modular', time.toc())
            po = f_sys('phenomena oriented')
            po.flatten()
            po.set_tolerance(rmol=1e-5, mol=1e-5, subsystems=True, method='fixed-point', maxiter=500)
            time.tic()
            po.simulate()
            print('- Phenomena oriented', time.toc())
            for s_sm, s_dp in zip(sm.streams, po.streams):
                actual = s_sm.mol
                value = s_dp.mol
                try:
                    assert_allclose(actual, value, rtol=0.01, atol=0.01)
                except:
                    print(f'- Multiple steady stages for {s_sm}: sm-{actual} po-{value}')
                    break

# %% Profiling and benchmarking utilities

class Tracker:
    __slots__ = ('system', 'run', 'streams', 
                 'adiabatic_stages', 'stages',
                 'profile_time', 'rtol', 'atol')
    cutoff_time = 200
    
    def __init__(self, name, algorithm, rtol=1e-16, atol=1e-6):
        sys = all_systems[name](algorithm)
        sys.flatten()
        sys._setup_units()
        bst.F.clear()
        self.system = sys
        if sys.algorithm == 'Sequential modular':
            self.run = sys.run_sequential_modular
        elif sys.algorithm == 'Phenomena oriented':
            self.run = sys.run_phenomena
        else:
            raise ValueError('invalid algorithm')
        self.rtol = rtol
        self.atol = atol
        self.profile_time = system_convergence_times[name]
        self.streams, self.adiabatic_stages, self.stages, = streams_and_stages(sys)
        
    def benchmark(self):
        f = self.run
        streams = self.streams
        adiabatic_stages = self.adiabatic_stages
        stages = self.stages
        total_time = self.cutoff_time
        time = bst.TicToc()
        net_time = 0
        temperatures = np.array([i.T for i in streams])
        flows = np.array([i.mol for i in streams])
        while net_time < total_time:
            time.tic()
            f()
            net_time += time.toc()
            new_temperatures = np.array([i.T for i in streams])
            new_flows = np.array([i.mol for i in streams])
            dF = np.abs(flows - new_flows)
            nonzero_index = dF > 0
            if nonzero_index.any(): 
                dF = dF[nonzero_index]
                F_max = np.maximum.reduce([abs(flows[nonzero_index]), abs(new_flows[nonzero_index])])
                rdF = dF / F_max
                rdF_max = rdF.max()
                dF_max = dF.max()
            else:
                dF_max = rdF_max = 0
            dT = np.abs(temperatures - new_temperatures)
            T_max = np.maximum.reduce([abs(temperatures), abs(new_temperatures)])
            dT_max = dT.max()
            rdT = dT / T_max
            rdT_max = rdT.max()
            if (rdF_max < self.rtol or dF_max < self.atol) and (rdT_max < self.rtol or dT_max < self.atol): break
            flows = new_flows
            temperatures = new_temperatures
        dM = max([abs(i.mass_balance_error()) for i in stages])
        dEs = [abs(dT_error(i)) for i in adiabatic_stages]
        dE = max(dEs) if dEs else 0
        return {
            'Time': net_time, 
            'Stream temperature': rdT_max,
            'Component flow rate': rdF_max, 
            'Absolute stream temperature': dT_max,
            'Absolute component flow rate': dF_max, 
            'Energy balance': dE, 
            'Material balance': dM,
            'Steady state flows': flows,
            'Steady state temperatures': temperatures,
        }
        
    def profile(self):
        f = self.run
        streams = self.streams
        adiabatic_stages = self.adiabatic_stages
        stages = self.stages
        total_time = self.profile_time
        time = bst.TicToc()
        flow_error = []
        energy_error = []
        material_error = []
        temperature_error = []
        record = []
        net_time = 0
        temperatures = np.array([i.T for i in streams])
        flows = np.array([i.mol for i in streams])
        while net_time < total_time:
            time.tic()
            f()
            net_time += time.toc()
            new_temperatures = np.array([i.T for i in streams])
            new_flows = np.array([i.mol for i in streams])
            record.append(net_time)
            dF = np.abs(flows - new_flows).sum()
            dT = np.abs(temperatures - new_temperatures).sum()
            flow_error.append(
                np.log10(dF + 1e-15)
            )
            temperature_error.append(
                np.log10(dT + 1e-15)
            )
            energy_error.append(
                np.log10(sum([abs(dT_error(i)) for i in adiabatic_stages]) + 1e-15)
            )
            material_error.append(
                np.log10(sum([abs(i.mass_balance_error()) for i in stages]) + 1e-15)
            )
            flows = new_flows
            temperatures = new_temperatures
        return {
            'Time': record, 
            'Stream temperature': temperature_error,
            'Component flow rate': flow_error, 
            'Energy balance': energy_error, 
            'Material balance': material_error,
        }

def streams_and_stages(sys):
    all_stages = []
    adiabatic_stages = []
    streams = []
    for i in sys.units:
        if hasattr(i, 'stages'):
            all_stages.extend(i.stages)
            for j in i.stages:
                streams.extend(j.outs + j.ins)
                if j.B_specification is None and j.T_specification is None:
                    adiabatic_stages.append(j)
        else:
            try:
                if i.B_specification is None and i.T_specification is None:
                    adiabatic_stages.append(j)
            except:
                pass
            all_stages.append(i)
            streams.extend(i.outs + i.ins)
    return (streams, adiabatic_stages, all_stages)

def dT_error(stage):
    if all([i.isempty() for i in stage.outs]): 
        return 0
    else:
        return (
            sum([i.H for i in stage.outs]) - sum([i.H for i in stage.ins])
        ) / sum([i.C for i in stage.outs])

def uncertainty_percent(xdx, ydy):
    x, dx = xdx
    y, dy = ydy
    if x == 0 and y == 0:
        return [1, 0]
    else:
        z = x / y
        if x == 0:
            dxx = 0
        else:
            dxx = dx / x
        if y == 0:
            dyy = 0
        else:
            dyy = dy / y
        return [100 * z, 100 * z * sqrt(dxx*dxx + dyy*dyy)]

# %% Benchmark plot

def plot_benchmark(systems=None, exclude=None, N=3, load=True, save=True):
    if systems is None: systems = list(all_systems)
    if exclude is not None: systems = [i for i in systems if i not in exclude]
    n_systems = len(systems)
    results = np.zeros([n_systems, 3])
    for m, sys in enumerate(systems):
        time = system_convergence_times[sys]
        sms = []
        pos = []
        if load:
            try:
                sm_name = f'sm_{time}_{sys}_benchmark_{N}.npy'
                file = os.path.join(simulations_folder, sm_name)
                with open(file, 'rb') as f: sms = pickle.load(f)
                po_name = f'po_{time}_{sys}_benchmark_{N}.npy'
                file = os.path.join(simulations_folder, po_name)
                with open(file, 'rb') as f: pos = pickle.load(f)
            except:
                for i in range(N):
                    sm = Tracker(sys, 'sequential modular').benchmark()
                    sms.append(sm)
                for i in range(N):
                    po = Tracker(sys, 'phenomena oriented').benchmark()
                    pos.append(po)
                if save:
                    sm_name = f'sm_{time}_{sys}_benchmark_{N}.npy'
                    file = os.path.join(simulations_folder, sm_name)
                    with open(file, 'wb') as f: pickle.dump(sms, f)
                    po_name = f'po_{time}_{sys}_benchmark_{N}.npy'
                    file = os.path.join(simulations_folder, po_name)
                    with open(file, 'wb') as f: pickle.dump(pos, f)
        else:
            sm = Tracker(sys, 'sequential modular').benchmark()
            for i in range(N):
                sm = Tracker(sys, 'sequential modular').benchmark()
                sms.append(sm)
            po = Tracker(sys, 'phenomena oriented').benchmark()
            for i in range(N):
                po = Tracker(sys, 'phenomena oriented').benchmark()
                pos.append(po)
            if save:
                sm_name = f'sm_{time}_{sys}_benchmark_{N}.npy'
                file = os.path.join(simulations_folder, sm_name)
                with open(file, 'wb') as f: pickle.dump(sms, f)
                po_name = f'po_{time}_{sys}_benchmark_{N}.npy'
                file = os.path.join(simulations_folder, po_name)
                with open(file, 'wb') as f: pickle.dump(pos, f)
        data = np.array([j['Time'] / i['Time'] for i, j in zip(sms, pos)])[1:]
        mean = np.mean(data)
        sm_better = mean > 1
        if sm_better:
            mean = 1 / mean
            data = 1 / data
        results[m] = [100 * mean, 100 * np.std(data), sm_better]
        # sm = dct_mean_std(sms, keys)
        # po = dct_mean_std(pos, keys)
        # system_results.append((sm, po))
    # Assume only time matters from here on
    # for i, (sm, po) in enumerate(system_results):
    #     results[i] = uncertainty_percent(po['Time'], sm['Time'])
    n_rows = 1
    n_cols = 1
    fs = 11
    bst.set_font(fs)
    bst.set_figure_size(aspect_ratio=0.9)
    fig, ax = plt.subplots(n_rows, n_cols)
    red = Color(fg='#f1777f').RGBn
    blue = Color(fg='#5fc1cf').RGBn
    black = Color(fg='#7b8580').RGBn
    # csm = Color(fg='#33BBEE').RGBn
    # cpo = Color(fg='#EE7733').RGBn
    yticks = (0, 25, 50, 75, 100)
    yticklabels = [f'{i}%' for i in yticks]
    xticks = list(range(n_systems))
    xticklabels = [system_labels[sys] for sys in systems]
    sm_index, = np.where(results[:, -1])
    po_index, = np.where(~results[:, -1].astype(bool))
    plt.errorbar([xticks[i] for i in sm_index], results[sm_index, 0], 
                 yerr=results[sm_index, 1], color=red, marker='x', linestyle='', ms=10,
                 capsize=5, capthick=1.5, ecolor=black)
    plt.errorbar([xticks[i] for i in po_index], results[po_index, 0], yerr=results[po_index, 1], color=blue,
                 marker='s', linestyle='', capsize=5, capthick=1.5, ecolor=black)
    plt.axhline(y=100, color='grey', ls='--', zorder=-1)
    plt.ylabel('Simulation time [%]')
    bst.utils.style_axis(
        ax, xticks=xticks, yticks=yticks,
        xticklabels=xticklabels,
        yticklabels=yticklabels,
    )
    plt.subplots_adjust(right=0.96, left=0.2, hspace=0, wspace=0)
    for i in ('svg', 'png'):
        name = f'PO_SM_{time}_benchmark_{N}.{i}'
        file = os.path.join(images_folder, name)
        plt.savefig(file, dpi=900, transparent=True)
    return fig, ax, sms, pos


# %%  Profile plot

def dct_mean_profile(dcts: list[dict], keys: list[str], ub: float, n=50):
    mean = {i: np.zeros(n) for i in keys}
    tmin = np.min([np.min(i['Time']) for i in dcts])
    t = np.linspace(0, ub, n)
    goods = [np.zeros(n) for i in range(len(keys))]
    for dct in dcts:
        for i, good in zip(keys, goods): 
            x = dct['Time']
            mean['Time'] = t
            y = dct[i]
            values = interpolate.interp1d(x[:len(y)], y, bounds_error=False)(t)
            mask = ~np.isnan(values)
            mean[i][mask] += values[mask]
            mean[i][t < tmin] = np.nan
            good[t > tmin] += 1
    for i, j in zip(keys, goods): mean[i][j > 0] /= j[j > 0]
    return mean

def dct_mean_std(dcts: list[dict], keys: list[str]):
    n = len(dcts)
    values = {i: np.zeros(n) for i in keys}
    for i, dct in enumerate(dcts):
        for key in keys: values[key][i] = dct[key]
    return {i: (values[i].mean(), values[i].std()) for i in keys}

def plot_profile(
        systems=None, N=3, load=True, save=True
    ):
    if systems is None: systems = list(all_systems)
    fs = 9
    bst.set_font(fs)
    keys = (
        'Component flow rate',
        'Stream temperature',
        # 'Stripping factor',
        'Material balance',
        'Energy balance',
    )
    units = (
        r'$[\mathrm{mol} \cdot \mathrm{hr}^{\mathrm{-1}}]$',
        r'$[\mathrm{K}]$',
        # r'$[-]$',
        r'$[\mathrm{mol} \cdot \mathrm{hr}^{\mathrm{-1}}]$',
        r'$[\mathrm{K}]$',
    )
    n_rows = len(units)
    n_cols = len(systems)
    if n_rows == 4:
        bst.set_figure_size(aspect_ratio=1.1)
    elif n_rows == 2:
        bst.set_figure_size(aspect_ratio=0.6)
    else:
        bst.set_figure_size(aspect_ratio=1)
    fig, all_axes = plt.subplots(n_rows, n_cols)
    if n_cols == 1:
        all_axes = np.reshape(all_axes, [n_rows, n_cols]) 
    yticks = [-15, -10, -5, 0, 5, 10]
    for m, sys in enumerate(systems):
        time = system_convergence_times[sys]
        axes = all_axes[:, m]
        sms = []
        pos = []
        if load:
            sm_name = f'sm_{time}_{sys}_profile.npy'
            file = os.path.join(simulations_folder, sm_name)
            with open(file, 'rb') as f: sms = pickle.load(f)
            po_name = f'po_{time}_{sys}_profile.npy'
            file = os.path.join(simulations_folder, po_name)
            with open(file, 'rb') as f: pos = pickle.load(f)
        else:
            for i in range(N):
                sm = Tracker(sys, 'sequential modular').profile()
                sms.append(sm)
            for i in range(N):
                po = Tracker(sys, 'phenomena oriented').profile()
                pos.append(po)
            if save:
                sm_name = f'sm_{time}_{sys}_profile.npy'
                file = os.path.join(simulations_folder, sm_name)
                with open(file, 'wb') as f: pickle.dump(sms, f)
                po_name = f'po_{time}_{sys}_profile.npy'
                file = os.path.join(simulations_folder, po_name)
                with open(file, 'wb') as f: pickle.dump(pos, f)
        tub = system_tickmarks[sys][-1]
        tub = min(tub, min([dct['Time'][-1] for dct in sms]), min([dct['Time'][-1] for dct in pos]))
        sm = dct_mean_profile(sms, keys, tub)
        po = dct_mean_profile(pos, keys, tub)
        csm = Color(fg='#33BBEE').RGBn
        cpo = Color(fg='#EE7733').RGBn
        yticklabels = m == 0
        if yticklabels:
            yticklabels = [r'$\mathrm{10}^{' f'{i}' '}$' for i in yticks]
        labels = {
            'Stripping factor': 'Stripping factor\nconvergence error',
            'Component flow rate': 'Flow rate\nconvergence error',
            'Stream temperature': 'Temperature\nconvergence error',
            'Material balance': 'Stage material\nbalance error',
            'Energy balance': 'Stage energy\nbalance error',
        }
        for n, (i, ax, u) in enumerate(zip(keys, axes, units)):
            plt.sca(ax)
            if n == n_rows-1: plt.xlabel('Time [s]')
            label = labels[i]
            if m == 0: plt.ylabel(f'{label}\n{u}')
            ysm = gaussian_filter(np.array(sm[i]), 0.2)
            ypo = gaussian_filter(np.array(po[i]), 0.2)
            tsm = np.array(sm['Time'])
            tpo = np.array(po['Time'])
            ysm[ysm < -14] = -14
            ypo[ypo < -14] = -14
            # plt.plot(xsm, ysm, lw=0, marker='s', color=csm, markersize=2.5)
            plt.plot(tsm, ysm, '--', color=csm, lw=1.5)
            # plt.plot(xpo, ypo, lw=0, marker='d', color=cpo, markersize=2.5)
            plt.plot(tpo, ypo, '-', color=cpo, lw=1.5)
            if m == 0 and n == 0:
                index = int(len(tsm) * 0.55)
                xy = x, y = (tsm[index], ysm[index])
                ax.annotate('Sequential\nmodular',
                    xy=xy, 
                    xytext=(x+0.1*tub, y+1),
                    arrowprops=dict(arrowstyle="->", color=csm),
                    color=csm,
                    fontsize=fs,
                    fontweight='bold',
                )
                index = int(len(tpo) * 0.45)
                xy = x, y = (tpo[index], ypo[index])
                ax.annotate('Phenomena\noriented',
                    xy=xy, 
                    xytext=(x-0.1*tub, y-5),
                    arrowprops=dict(arrowstyle="->", color=cpo),
                    ha='right',
                    color=cpo,
                    fontsize=fs,
                    fontweight='bold',
                )
            xticks = system_tickmarks[sys]
            xticklabels = xtick0 = n == n_rows-1
            xtickf = m == n_cols-1
            ytick0 = n == n_rows-1
            ytickf = n == 0
            plt.xlim(0, xticks[-1])
            plt.ylim(-15, 10)
            bst.utils.style_axis(
                ax, xticks=xticks, yticks=yticks, 
                xtick0=xtick0, xtickf=xtickf, ytick0=ytick0, ytickf=ytickf,
                xticklabels=xticklabels,
                yticklabels=yticklabels,
            )
    plt.subplots_adjust(right=0.96, left=0.2, hspace=0, wspace=0)
    letter_color = c.neutral.shade(25).RGBn
    titles = [system_titles[i] for i in systems]
    for ax, letter in zip(all_axes[0], titles):
        plt.sca(ax)
        ylb, yub = plt.ylim()
        xlb, xub = plt.xlim()
        plt.text((xlb + xub) * 0.5, ylb + (yub - ylb) * 1.2, letter, color=letter_color,
                  horizontalalignment='center',verticalalignment='center',
                  fontsize=10, fontweight='bold')
    for i in ('svg', 'png'):
        name = f'PO_SM_profile.{i}'
        file = os.path.join(images_folder, name)
        plt.savefig(file, dpi=900, transparent=True)
    return fig, all_axes, sms, pos
