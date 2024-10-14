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

__all__ = (
    'BenchmarkModel',
    'run_monte_carlo',
    'plot_monte_carlo',
    'plot_benchmark',
    'plot_profile',
    'register',
    'test_convergence',
    'Tracker'
)

# %% Make simulation rigorous to achieve lower tolerances

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
# bst.units.stage.PhasePartition.S_relaxation_factor = 0
# bst.units.stage.PhasePartition.B_relaxation_factor = 0
# bst.units.stage.PhasePartition.K_relaxation_factor = 0
# bst.units.stage.PhasePartition.T_relaxation_factor = 0

# %% Uncertainty analysis

class BenchmarkModel:
    
    def __init__(self):
        model = bst.Model(None, specification=lambda: None)
        parameter = model.parameter
        self.model = model
        
        def load_trackers():
            system = 'acetic_acid_complex'
            self.sm_tracker = Tracker(
                system, 'sequential modular',
                extraction_stages=self.extraction_stages, 
                raffinate_distillation_stages=self.raffinate_distillation_stages,
                extract_distillation_stages=self.raffinate_distillation_stages,
                atol=1e-10
            )
            self.sm_system = self.sm_tracker.system
            self.sm_recycles = self.sm_system.get_all_recycles()
            self.po_tracker = Tracker(
                system, 'phenomena oriented', 
                extraction_stages=self.extraction_stages, 
                raffinate_distillation_stages=self.raffinate_distillation_stages,
                extract_distillation_stages=self.raffinate_distillation_stages,
                atol=1e-10
            )
            self.po_system = self.po_tracker.system
            self.po_recycles = self.po_system.get_all_recycles()
        
        @parameter(units='-', bounds=[1, 12], hook=lambda x: int(round(x)))
        def set_extraction_stages(extraction_stages):
            self.extraction_stages = extraction_stages
        
        @parameter(units='-', bounds=[1, 12], hook=lambda x: int(round(x)))
        def set_raffinate_distillation_stages(raffinate_distillation_stages):
            self.raffinate_distillation_stages = raffinate_distillation_stages
            
        @parameter(units='-', bounds=[1, 12], hook=lambda x: int(round(x)))
        def set_extract_distillation_stages(extract_distillation_stages):
            self.extract_distillation_stages = extract_distillation_stages
        
        # benchmark = tracker.benchmark()
        # chemicals = [i.ID for i in bst.settings.chemicals]
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
        
        @model.metric(units='%', element='Sequential modular')
        def distance_from_steady_state():
            sm_b = self.sm_tracker.benchmark()
            po_b = self.po_tracker.benchmark()
            sm_e = sm_b['Error history']
            po_e = po_b['Error history']
            t_lb = max(po_e[0].time, sm_e[0].time)
            t_ub = min(po_e[0].time, sm_e[0].time)
            x, y = zip(*sm_e)
            t = np.linspace(t_lb, t_ub, 2**10 + 1)
            sm_y = interpolate.interp1d(x, y, bounds_error=False)(t)
            x, y = zip(*sm_e)
            po_y = interpolate.interp1d(x, y, bounds_error=False)(t)
            return romb(po_y / sm_y, t[1] - t[0]) / (t_ub - t_lb)
          
        @model.metric(units='%', element='Phenomena-oriented')
        def distance_from_steady_state():
            sm_b = self.sm_tracker.benchmark()
            po_b = self.po_tracker.benchmark()
            sm_e = sm_b['Error history']
            po_e = po_b['Error history']
            t_lb = max(po_e[0].time, sm_e[0].time)
            t_ub = min(po_e[0].time, sm_e[0].time)
            x, y = zip(*sm_e)
            t = np.linspace(t_lb, t_ub, 2**10 + 1)
            sm_y = interpolate.interp1d(x, y, bounds_error=False)(t)
            x, y = zip(*sm_e)
            po_y = interpolate.interp1d(x, y, bounds_error=False)(t)
            return romb(po_y / sm_y, t[1] - t[0]) / (t_ub - t_lb)
        
        @model.metric(units='%')
        def average_relative_distance_from_steady_state():
            sm_b = self.sm_tracker.benchmark()
            po_b = self.po_tracker.benchmark()
            sm_e = sm_b['Error history']
            po_e = po_b['Error history']
            t_lb = max(po_e[0].time, sm_e[0].time)
            t_ub = min(po_e[0].time, sm_e[0].time)
            x, y = zip(*sm_e)
            t = np.linspace(t_lb, t_ub, 2**10 + 1)
            sm_y = interpolate.interp1d(x, y, bounds_error=False)(t)
            x, y = zip(*sm_e)
            po_y = interpolate.interp1d(x, y, bounds_error=False)(t)
            return romb(po_y / sm_y, t[1] - t[0]) / (t_ub - t_lb)
        
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

# %% System creation

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

def register(name, title, time, tickmarks, label):
    f = getattr(systems, 'create_system_' + name, None) or getattr(systems, 'create_' + name + '_system')
    all_systems[name] = f
    system_titles[name] = title
    system_convergence_times[name] = time
    system_tickmarks[name] = tickmarks
    system_labels[name] = label

register(
    'acetic_acid_reactive_purification', 'Acetic acid reactive purification',
    10, [2, 4, 6, 8, 10], 'AA\nr. sep.'
)
register(
    'acetic_acid_simple', 'Subsystem',
    40, [0, 6, 12, 18, 24], 'AcOH sep.\nsubsystem'
)
register(
    'acetic_acid_complex_decoupled', 'Shortcut system',
    8, [0, 2, 4, 6, 8], 'AcOH sep.\nshortcut',
)
register(
    'acetic_acid_complex', 'Rigorous system',
    400, [0, 60, 120, 180, 240, 300], 'AcOH sep.\nrigorous'
)
register(
    'alcohol_narrow_flash', 'Alcohol flash narrow',
    0.05, [0.01, 0.02, 0.03, 0.04, 0.05], 'Alcohol\nflash\nnarrow'
)
register(
    'alcohol_wide_flash', 'Alcohol flash wide',
    0.05, [0.01, 0.02, 0.03, 0.04, 0.05], 'Alcohol\nflash\nwide'
)
register(
    'butanol_purification', 'Butanol purification',
    0.4, [0.1, 0.2, 0.3, 0.4], 'BtOH\nsep.'
)
register(
    'ethanol_purification', 'Ethanol purification',
    0.25, [0.05, 0.1, 0.15, 0.20, 0.25], 'EtOH\nsep.',
)
register(
    'haber_bosch_process', 'Haber-Bosch process',
    10, [2, 4, 6, 8, 10], 'HB\nprocess.'
)
register(
    'hydrocarbon_narrow_flash', 'Hydrocarbon flash narrow',
    0.05, [0.01, 0.02, 0.03, 0.04, 0.05], 'Alkane\nflash\nnarrow'
)
register(
    'hydrocarbon_wide_flash', 'Hydrocarbon flash wide',
    0.1, [0.02, 0.04, 0.06, 0.08, 0.10], 'Alkane\nflash\nwide'
)
register(
    'lactic_acid_purification', 'Lactic acid purification',
    10, [2, 4, 6, 8, 10], 'LA\nsep.'
)
# %% Testing

def test_convergence(systems=None):
    if systems is None: systems = list(all_systems)
    with catch_warnings():
        filterwarnings('ignore')
        time = bst.TicToc()
        for sys in systems:
            f_sys = all_systems[sys]
            bst.F.set_flowsheet('SM')
            sm = f_sys('sequential modular')
            sm.flatten()
            sm.set_tolerance(rmol=1e-5, mol=1e-5, subsystems=True, method='fixed-point', maxiter=500)
            time.tic()
            sm.simulate()
            print(sys)
            print('- Sequential modular', time.toc())
            bst.F.set_flowsheet('PO')
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
                    if s_sm.source:
                        print(f'- Multiple steady stages for {s_sm}, {s_sm.source}-{s_sm.source.outs.index(s_sm)}: sm-{actual} po-{value}')
                    else:
                        print(f'- Multiple steady stages for {s_sm}, {s_sm.sink.ins.index(s_sm)}-{s_sm.sink}: sm-{actual} po-{value}')

# %% Profiling and benchmarking utilities

class ErrorPoint(NamedTuple):
    time: float; error: float

default_steady_state_cutoff = 2
    
def steady_state_error(profile, steady_state_cutoff=None):
    if steady_state_cutoff is None: steady_state_cutoff = default_steady_state_cutoff
    minimum_error = np.log10(10 ** profile['Component flow rate error'][-1] + 10 ** profile['Temperature error'][-1])
    return minimum_error + steady_state_cutoff
    
def benchmark(profile, steady_state_error=None):
    if steady_state_error is None: steady_state_error = steady_state_error(profile)
    time = profile['Time']
    error = np.log10(10 ** profile['Component flow rate error'] + 10 ** profile['Temperature error'])
    return interpolate.interp1d(error, time, bounds_error=False)(steady_state_error)
        

class Tracker:
    __slots__ = ('system', 'run', 'streams', 
                 'adiabatic_stages', 'stages',
                 'profile_time', 'rtol', 'atol',
                 'kwargs', 'name', 'algorithm')
    
    def __init__(self, name, algorithm, rtol=1e-16, atol=1e-6, **kwargs):
        sys = all_systems[name](algorithm, **kwargs)
        sys.flatten()
        sys._setup_units()
        sys.flowsheet.clear()
        self.system = sys
        self.name = name
        self.kwargs = kwargs
        if sys.algorithm == 'Sequential modular':
            self.run = sys.run_sequential_modular
        elif sys.algorithm == 'Phenomena oriented':
            self.run = sys.run_phenomena
        else:
            raise ValueError('invalid algorithm')
        self.algorithm = algorithm
        self.rtol = rtol
        self.atol = atol
        self.profile_time = system_convergence_times[name]
        self.streams, self.adiabatic_stages, self.stages, = streams_and_stages(sys)
        
    def estimate_steady_state(self, load=True): # Uses optimization to estimate steady state
        options = '_'.join(['{i}_{j}' for i, j in self.kwargs.items()])
        if options:
            file_name = f"{self.name}_{options}_steady_state"
        else:
            file_name = f"{self.name}_steady_state"
        file = os.path.join(simulations_folder, file_name)
        if load:
            try:
                with open(file, 'rb') as f: return pickle.load(f)
            except: pass
        sys = all_systems[self.name]('phenomena-oriented', **self.kwargs)
        sys.flatten()
        try: sys.simulate()
        except: pass
        for i in range(100): sys.run()
        sys.flowsheet.clear()
        try: sys.simulate()
        except: pass
        p = bst.units.stage.PhasePartition
        settings = (
            p.S_relaxation_factor, p.B_relaxation_factor, 
            p.K_relaxation_factor, p.T_relaxation_factor
        )
        p.S_relaxation_factor = 0
        p.B_relaxation_factor = 0
        p.K_relaxation_factor = 0
        p.T_relaxation_factor = 0
        try:
            results = self._estimate_steady_state(sys)
            with open(file, 'wb') as f: pickle.dump(results, f)
        finally:
            (p.S_relaxation_factor, p.B_relaxation_factor, 
             p.K_relaxation_factor, p.T_relaxation_factor) = settings
        return results
        
    def _estimate_steady_state(self, system):
        po = system
        streams, adiabatic_stages, all_stages = streams_and_stages(po)
        flows = np.array([i.mol for i in streams])
        Ts = np.array([i.T for i in streams])
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
            elif not isinstance(i, (bst.Mixer, bst.Separator, bst.SinglePhaseStage)):
                raise ValueError(f'cannot optimize {i!r}')
        all_stages = stages + shortcuts
        N_chemicals = po.units[0].chemicals.size
        S = np.array([i.B * i.K for i in all_stages])
        S_full = S
        S_index = S_index = [
            i for i, j in enumerate(all_stages)
            if not (getattr(j, 'B_specification', None) == 0 or getattr(j, 'B_specification', None) == np.inf)
        ]
        lle_stages = []
        vle_stages = []
        for i in S_index:
            s = all_stages[i]
            if not isinstance(s, bst.StageEquilibrium): continue
            phases = getattr(s, 'phases', None)
            if phases == ('g', 'l'):
                vle_stages.append(s)
            elif phases == ('L', 'l'):
                lle_stages.append(s)
        S = S[S_index]
        N_S_index = len(S_index)
        lnS = np.log(S).flatten()
        count = [0]
        
        energy_error = lambda stage: abs(
            (sum([i.H for i in stage.outs]) - sum([i.H for i in stage.ins])) / sum([i.C for i in stage.outs])
        )
        
        def B_error(stage):
            B = getattr(stage, 'B_specification', None)
            if B is not None:
                top, bottom = stage.partition.outs
                return 100 * (top.F_mol / bottom.F_mol - B)
            else:
                return 0
        
        def T_equilibrium_error(stage):
            if len(stage.outs) == 2:
                top, bottom = stage.outs
            else:
                top, bottom = stage.partition.outs
            return (top.dew_point_at_P().T - bottom.T)
        
        def lnS_objective(lnS):
            S_original = np.exp(lnS)
            S = S_full
            S[S_index] = S_original.reshape([N_S_index, N_chemicals])
            for i in S_index:
                s = all_stages[i]
                if getattr(s, 'B_specification', None) is not None:
                    s.K = S[i] / s.B_specification
                else:
                    s.B = 1 # Work around S = K * B
                    s.K = S[i]
            with po.stage_configuration(aggregated=False) as conf:
                conf._solve_material_flows(composition_sensitive=False)
            # breakpoint()
            for i in vle_stages:
                partition = i.partition
                # print('----')
                # print(i.K * i.B)
                partition._run_decoupled_KTvle()
                if i.B_specification is None: partition._run_decoupled_B()
                # print(i.K * i.B)
                T = partition.T
                for i in (partition.outs + i.outs): i.T = T
            for i in shortcuts:
                for s in i.outs:
                    if s.phase == 'l': 
                        bp = s.bubble_point_at_P()
                        s.T = bp.T
                    elif s.phase == 'g': 
                        dp = s.dew_point_at_P()
                        s.T = dp.T
            # Ts = np.array([i.T for i in lle_stages])
            # for i in range(10):
            #     Ts_new = np.array([i.T for i in lle_stages])
            #     with po.stage_configuration(aggregated=False) as conf:
            #         conf.solve_energy_departures(temperature_only=True)
            #     for i in lle_stages:
            #         for j in i.outs: j.T = i.T
            #     print(np.abs(Ts_new - Ts).sum())   
            #     if np.abs(Ts_new - Ts).sum() < 1e-9: break
            #     Ts = Ts_new
            for i in lle_stages:
                i.partition._run_lle(update=False)
            total_energy_error = sum([energy_error(stage) for stage in stages if stage.B_specification is None and stage.T_specification is None])
            specification_errors = np.array([B_error(all_stages[i]) for i in S_index])
            # temperature_errors = np.array([T_equilibrium_error(i) for i in vle_stages])
            for i in shortcuts: i._run()
            S_new = np.array([(all_stages[i].K * all_stages[i].B) for i in S_index]).flatten()
            splits_new = S_new / (S_new + 1)
            splits = S_original / (S_original + 1)
            diff = (splits_new - splits)
            total_split_error = (diff * diff).sum()
            total_specification_error = (specification_errors * specification_errors).sum()
            # total_temperature_error = (temperature_errors * temperature_errors).sum()
            total_error = total_split_error + total_energy_error + total_specification_error #+ total_temperature_error
            err = np.sqrt(total_error)
            if not count[0] % 100:
                print(err)
                print(total_split_error, total_energy_error, total_specification_error)
                po.show()
            count[0] += 1
            return err
        benchmark = lnS_objective(lnS)
        result = minimize(
            lnS_objective, 
            lnS,
            tol=0.5, 
            method='CG',
            options=dict(maxiter=80),
        )
        optimization = lnS_objective(result.x)
        if benchmark > optimization:
            streams, adiabatic_stages, all_stages = streams_and_stages(po)
            flows = np.array([i.mol for i in streams])
            Ts = np.array([i.T for i in streams])
            return flows, Ts, optimization
        else:
            return flows, Ts, benchmark
        
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
        net_time = 0
        temperatures = np.array([i.T for i in streams])
        flows = np.array([i.mol for i in streams])
        temperature_history = []
        flow_history = []
        record = []
        steady_state_flows, steady_state_temperatures, _ = self.estimate_steady_state()
        while net_time < total_time:
            time.tic()
            f()
            net_time += time.toc()
            new_temperatures = np.array([i.T for i in streams])
            new_flows = np.array([i.mol for i in streams])
            record.append(net_time)
            flow_history.append(new_flows)
            temperature_history.append(new_temperatures)
            dF = np.abs(flows - new_flows).sum()
            dT = np.abs(temperatures - new_temperatures).sum()
            flow_error.append(
                np.log10(dF + 1e-25)
            )
            temperature_error.append(
                np.log10(dT + 1e-25)
            )
            energy_error.append(
                np.log10(sum([abs(dT_error(i)) for i in adiabatic_stages]) + 1e-25)
            )
            material_error.append(
                np.log10(sum([abs(i.mass_balance_error()) for i in stages]) + 1e-25)
            )
            flows = new_flows
            temperatures = new_temperatures
        cfe = np.log10([np.abs(steady_state_flows - i).sum() + 1e-25 for i in flow_history])
        te = np.log10([np.abs(steady_state_temperatures - i).sum() + 1e-25 for i in temperature_history])
        return {
            'Time': record, 
            'Component flow rate error': cfe,
            'Temperature error': te,
            'Stream temperature': temperature_error,
            'Component flow rate': flow_error, 
            'Energy balance': energy_error, 
            'Material balance': material_error,
        }

def streams_and_stages(sys):
    all_stages = []
    adiabatic_stages = []
    streams = []
    past_streams = set()
    # print(f'----{sys.algorithm}----')
    for i in sorted(sys.unit_path, key=lambda x: x.ID):
        # print(i)
        if hasattr(i, 'stages'):
            all_stages.extend(i.stages)
            for j in i.stages:
                new_streams = [i for i in (j.outs + j.ins) if i.imol not in past_streams]
                streams.extend(new_streams)
                past_streams.update([i.imol for i in new_streams])
                if j.B_specification is None and j.T_specification is None:
                    adiabatic_stages.append(j)
        else:
            try:
                if i.B_specification is None and i.T_specification is None:
                    adiabatic_stages.append(j)
            except:
                pass
            all_stages.append(i)
            new_streams = [j for j in (i.outs + i.ins) if j.imol not in past_streams]
            streams.extend(new_streams)
            past_streams.update([i.imol for i in new_streams])
    return (streams, adiabatic_stages, all_stages)

def dT_error(stage):
    if all([i.isempty() for i in stage.outs]): 
        return 0
    else:
        return (
            sum([i.H for i in stage.outs]) - sum([i.H for i in stage.ins])
        ) / sum([i.C for i in stage.outs])

def division_mean_std(xdx, ydy):
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
        return [z, z * sqrt(dxx*dxx + dyy*dyy)]

# %% Benchmark plot

def plot_benchmark(systems=None, exclude=None, N=3, load=True, save=True):
    if systems is None: systems = list(all_systems)
    if exclude is not None: systems = [i for i in systems if i not in exclude]
    n_systems = len(systems)
    results = np.zeros([n_systems, 3])
    Ns = N * np.ones(n_systems, int)
    for m, sys in enumerate(systems):
        N = Ns[m]
        time = system_convergence_times[sys]
        sms = []
        pos = []
        if load:
            try:
                sm_name = f'sm_{time}_{sys}_profile.npy'
                file = os.path.join(simulations_folder, sm_name)
                with open(file, 'rb') as f: sms = pickle.load(f)
                po_name = f'po_{time}_{sys}_profile.npy'
                file = os.path.join(simulations_folder, po_name)
                with open(file, 'rb') as f: pos = pickle.load(f)
            except:
                for i in range(N):
                    sm = Tracker(sys, 'sequential modular').profile()
                    sms.append(sm)
                for i in range(N):
                    po = Tracker(sys, 'phenomena oriented').profile()
                    pos.append(po)
                if save:
                    sm_name = f'sm_{time}_{sys}_benchmark_{N}.npy'
                    file = os.path.join(simulations_folder, sm_name)
                    with open(file, 'wb') as f: pickle.dump(sms, f)
                    po_name = f'po_{time}_{sys}_benchmark_{N}.npy'
                    file = os.path.join(simulations_folder, po_name)
                    with open(file, 'wb') as f: pickle.dump(pos, f)
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
        cutoff = max([steady_state_error(i) for i in sms + pos])
        sms = np.array([benchmark(i, cutoff) for i in sms])
        pos = np.array([benchmark(i, cutoff) for i in pos])
        pos_mean_std = [np.mean(pos), np.std(pos)]
        sms_mean_std = [np.mean(sms), np.std(sms)]
        mean, std = division_mean_std(pos_mean_std, sms_mean_std)
        sm_better = False
        # sm_better = mean > 1
        # if sm_better: mean = 1 / mean
        results[m] = [100 * mean, 100 * std, sm_better]
        # sm = dct_mean_std(sms, keys)
        # po = dct_mean_std(pos, keys)
        # system_results.append((sm, po))
    # Assume only time matters from here on
    # for i, (sm, po) in enumerate(system_results):
    #     results[i] = uncertainty_percent(po['Time'], sm['Time'])
    n_rows = 1
    n_cols = 1
    fs = 10
    bst.set_font(fs)
    bst.set_figure_size(aspect_ratio=0.9, width='half')
    fig, ax = plt.subplots(n_rows, n_cols)
    red = Color(fg='#f1777f').RGBn
    blue = Color(fg='#5fc1cf').RGBn
    black = Color(fg='#7b8580').RGBn
    # csm = Color(fg='#33BBEE').RGBn
    # cpo = Color(fg='#EE7733').RGBn
    yticks = (0, 25, 50, 75, 100, 125)
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
    return fig, ax, results


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
            values = interpolate.interp1d(x, y, bounds_error=False)(t)
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
    fs = 12
    bst.set_font(fs)
    keys = (
        'Component flow rate error',
        # 'Temperature error',
        # 'Component flow rate',
        # 'Stream temperature',
        # 'Stripping factor',
        # 'Material balance',
        # 'Energy balance',
    )
    units = (
        r'$[\mathrm{mol} \cdot \mathrm{hr}^{\mathrm{-1}}]$',
        # r'$[\mathrm{K}]$',
        # r'$[\mathrm{mol} \cdot \mathrm{hr}^{\mathrm{-1}}]$',
        # r'$[\mathrm{K}]$',
        # r'$[-]$',
        # r'$[\mathrm{mol} \cdot \mathrm{hr}^{\mathrm{-1}}]$',
        # r'$[\mathrm{K}]$',
    )
    n_rows = len(units)
    n_cols = len(systems)
    if n_cols >= 2: 
        width = 'full'
    else:
        width = 'half'
    if n_rows == 4:
        bst.set_figure_size(aspect_ratio=1.1 / n_cols, width=width)
    elif n_rows == 2:
        if n_cols >= 2:
            aspect_ratio = 0.75
        else:
            aspect_ratio = 1.5
        bst.set_figure_size(aspect_ratio=aspect_ratio, width=width)
    else:
        if n_cols >= 2:
            aspect_ratio = 0.75 / 2
        elif n_rows == 1:
            aspect_ratio = 1.4 / 2
        else:
            aspect_ratio = 1.5 / 2
        bst.set_figure_size(aspect_ratio=aspect_ratio, width=width)
    fig, all_axes = plt.subplots(n_rows, n_cols)
    if n_rows == 1:
        if n_cols == 1:
            all_axes = np.array([[all_axes]])
        else:
            all_axes = all_axes.reshape([n_rows, n_cols])
    if n_cols == 1:
        all_axes = np.reshape(all_axes, [n_rows, n_cols]) 
    for m, sys in enumerate(systems):
        time = system_convergence_times[sys]
        axes = all_axes[:, m]
        sms = []
        pos = []
        for i in range(N):
            sm_name = f'sm_{time}_{sys}_profile_{i}.npy'
            file = os.path.join(simulations_folder, sm_name)
            if load:
                try:
                    with open(file, 'rb') as f: sm = pickle.load(f)
                except:
                    sm = Tracker(sys, 'sequential modular').profile()
            if save:
                sm_name = f'sm_{time}_{sys}_profile_{i}.npy'
                file = os.path.join(simulations_folder, sm_name)
                with open(file, 'wb') as f: pickle.dump(sm, f)
            sms.append(sm)
        for i in range(N):
            po_name = f'po_{time}_{sys}_profile_{i}.npy'
            file = os.path.join(simulations_folder, po_name)
            if load:
                try:
                    with open(file, 'rb') as f: po = pickle.load(f)
                except:
                    po = Tracker(sys, 'phenomena oriented').profile()
            if save:
                po_name = f'po_{time}_{sys}_profile_{i}.npy'
                file = os.path.join(simulations_folder, po_name)
                with open(file, 'wb') as f: pickle.dump(po, f)
            pos.append(po)
        tub = system_tickmarks[sys][-1]
        tub = min(tub, min([dct['Time'][-1] for dct in sms]), min([dct['Time'][-1] for dct in pos]))
        sm = dct_mean_profile(sms, keys, tub)
        po = dct_mean_profile(pos, keys, tub)
        csm = Color(fg='#33BBEE').RGBn
        cpo = Color(fg='#EE7733').RGBn
        labels = {
            'Component flow rate error': 'Flow rate error',
            'Temperature error': 'Temperature error',
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
            if m == 0: plt.ylabel(f'{label} {u}')
            ysm = np.array(sm[i])
            ypo = np.array(po[i])
            ysm = gaussian_filter(ysm, 0.2)
            ypo = gaussian_filter(ypo, 0.2)
            tsm = np.array(sm['Time'])
            tpo = np.array(po['Time'])
            plt.plot(tsm, ysm, '--', color=csm, lw=1.5, alpha=0.5)
            plt.plot(tsm, ysm, lw=0, marker='.', color=csm, markersize=2.5)
            plt.plot(tpo, ypo, lw=0, marker='.', color=cpo, markersize=2.5)
            plt.plot(tpo, ypo, '-', color=cpo, lw=1.5, alpha=0.5)
            ypo_finite = ypo[~np.isnan(ypo)]
            yticklabels = m == 0
            lb = -15
            ub = 10
            if yticklabels:
                # lb = ypo_finite.min()
                # ub = ypo_finite.max()
                step = 5
                # if (ub - lb) < 18:
                #     step = 3    
                # elif (ub - lb) < 20:
                #     step = 5
                # else:
                #     step = 10
                # lb = int(np.floor(ypo_finite.min() / step) * step)
                # ub = int(np.ceil(ypo_finite.max() / step) * step)
                yticks = [i for i in range(lb, ub+1, step)]
                yticklabels = [r'$\mathrm{10}^{' f'{i}' '}$' for i in yticks]
            if m == 0 and n == 0:
                index = int(len(tsm) * 0.5)
                xy = x, y = (tsm[index], ysm[index])
                ax.annotate('Sequential\nmodular',
                    xy=xy, 
                    xytext=(x-0.01*tub, y+1),
                    # arrowprops=dict(arrowstyle="->", color=csm),
                    color=csm,
                    fontsize=fs,
                    fontweight='bold',
                )
                index = int(len(tpo) * 0.5)
                xy = x, y = (tpo[index], ypo[index])
                ax.annotate('Phenomena\noriented',
                    xy=xy, 
                    xytext=(x+0.025*tub, y-5),
                    # arrowprops=dict(arrowstyle="->", color=cpo),
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
            plt.ylim(lb, ub)
            bst.utils.style_axis(
                ax, xticks=xticks, yticks=yticks, 
                xtick0=xtick0, xtickf=xtickf, ytick0=ytick0, ytickf=ytickf,
                xticklabels=xticklabels,
                yticklabels=yticklabels,
            )
    letter_color = c.neutral.shade(25).RGBn
    titles = [system_titles[i] for i in systems]
    # for ax, letter in zip(all_axes[0], titles):
    #     plt.sca(ax)
    #     ylb, yub = plt.ylim()
    #     xlb, xub = plt.xlim()
    #     plt.text((xlb + xub) * 0.5, ylb + (yub - ylb) * 1.1, letter, color=letter_color,
    #               horizontalalignment='center',verticalalignment='center',
    #               fontsize=fs, fontweight='bold')
    left = 0.1
    top = 0.85
    if n_rows == 2:
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
    for i in ('svg', 'png'):
        name = f'PO_SM_profile.{i}'
        file = os.path.join(images_folder, name)
        plt.savefig(file, dpi=900, transparent=True)
    for i in ('svg', 'png'):
        system_names = '_'.join(systems)
        name = f'PO_SM_{system_names}_profile.{i}'
        file = os.path.join(images_folder, name)
        plt.savefig(file, dpi=900, transparent=True)
    return fig, all_axes, sms, pos
