# -*- coding: utf-8 -*-
"""
Created on Sat Mar 16 13:38:11 2024

@author: cortespea
"""
from matplotlib import pyplot as plt
from numpy.testing import assert_allclose
import biosteam as bst
import numpy as np
from biosteam.utils import colors as c
from colorpalette import Color
from scipy import interpolate
from scipy.ndimage import gaussian_filter
from math import sqrt
import os
import pickle

__all__ = (
    'plot_benchmark',
    'plot_profile',
    'register',
)

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

def _register(f, name, title, time, tickmarks):
    all_systems[name] = f
    system_titles[name] = title
    system_convergence_times[name] = time
    system_tickmarks[name] = tickmarks
    return f

def register(name, title, time, tickmarks):
    return lambda f: _register(f, name, title, time, tickmarks)

def profile_sequential_modular(system, total_time):
    sm = all_systems[system]('sequential modular')
    sm.flatten()
    sm._setup_units()
    data = profile(sm.run_sequential_modular, *streams_and_stages(sm), total_time, init=sm.run_sequential_modular)
    bst.F.clear()
    return data

def benchmark_sequential_modular(system, total_time):
    sm = all_systems[system]('sequential modular')
    sm.flatten()
    sm._setup_units()
    data = benchmark(sm.run_sequential_modular, *streams_and_stages(sm), total_time, init=sm.run_sequential_modular)
    bst.F.clear()
    return data

def profile_phenomena_oriented(system, total_time):
    po = all_systems[system]('phenomena oriented')
    po.flatten()
    po._setup_units()
    data = profile(po.run_phenomena, *streams_and_stages(po), total_time,
                   init=po.run_sequential_modular)
    bst.F.clear()
    return data

def benchmark_phenomena_oriented(system, total_time):
    po = all_systems[system]('phenomena oriented')
    po.flatten()
    po._setup_units()
    data = benchmark(po.run_phenomena, *streams_and_stages(po), total_time,
                   init=po.run_sequential_modular)
    bst.F.clear()
    return data

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
                all_stages.append(i)
            except:
                pass
            streams.extend(i.outs + i.ins)
    return (streams, adiabatic_stages, all_stages)

def dT_error(stage):
    if all([i.isempty() for i in stage.outs]): 
        return 0
    else:
        return (
            sum([i.H for i in stage.outs]) - sum([i.H for i in stage.ins])
        ) / sum([i.C for i in stage.outs])

def benchmark(f, streams, adiabatic_stages, stages, total_time, init):
    time = bst.TicToc()
    net_time = 0
    init()
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
        dF = dF[nonzero_index]
        dF_max = np.maximum.reduce([abs(flows[nonzero_index]), abs(new_flows[nonzero_index])])
        rdF = dF / dF_max
        dT = np.abs(temperatures - new_temperatures).sum()
        dT_max = np.maximum.reduce([abs(temperatures), abs(new_temperatures)])
        rdT = dT / dT_max
        if (rdF.max() < 1e-9 or dF_max.max() < 1e-6) and (rdT.max() < 1e-9 or dT_max.max() < 1e-6): break
        flows = new_flows
        temperatures = new_temperatures
    dM = sum([abs(i.mass_balance_error()) for i in stages])
    dE = sum([abs(dT_error(i)) for i in adiabatic_stages])
    return {
        'Time': net_time, 
        'Stream temperature': dT,
        'Component flow rate': dF, 
        'Energy balance': dE, 
        'Material balance': dM,
    }

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
        return [100 * z, 100 * z * sqrt(dxx + dyy)]

def plot_benchmark(systems=None, N=10, load=True, save=True):
    if systems is None: systems = list(all_systems)
    fs = 9
    bst.set_font(fs)
    keys = (
        'Time',
    )
    units = (
        r'$[s]$',
    )
    n_rows = 1
    n_cols = 1
    bst.set_figure_size(aspect_ratio=1)
    fig, ax = plt.subplots(n_rows, n_cols)
    yticks = [-15, -10, -5, 0, 5, 10]
    system_results = []
    n_systems = len(systems)
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
                    sm = benchmark_sequential_modular(sys, time)
                    po = benchmark_phenomena_oriented(sys, time)
                    sms.append(sm)
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
                sm = benchmark_sequential_modular(sys, time)
                po = benchmark_phenomena_oriented(sys, time)
                sms.append(sm)
                pos.append(po)
            if save:
                sm_name = f'sm_{time}_{sys}_benchmark_{N}.npy'
                file = os.path.join(simulations_folder, sm_name)
                with open(file, 'wb') as f: pickle.dump(sms, f)
                po_name = f'po_{time}_{sys}_benchmark_{N}.npy'
                file = os.path.join(simulations_folder, po_name)
                with open(file, 'wb') as f: pickle.dump(pos, f)
        sm = dct_mean_std(sms, keys)
        po = dct_mean_std(pos, keys)
        system_results.append((sm, po))
    # Assume only time matters from here on
    results = np.zeros([n_systems, 2])
    for i, (sm, po) in enumerate(system_results):
        results[i] = uncertainty_percent(po['Time'], sm['Time'])
    csm = Color(fg='#33BBEE').RGBn
    cpo = Color(fg='#EE7733').RGBn
    mean, std = results['Time'].T
    yticks = (0, 50, 100, 150, 200)
    yticklabels = [f'{i}%' for i in yticks]
    xticks = list(range(n_systems))
    xticklabels = []
    for sys in systems:
        names = sys.upper().split('_')
        xticklabels.append(
            ''.join([i[0] for i in names])
        )
    plt.ylabel('Time [% Sequential modular]')
    plt.scatter()
    bst.utils.style_axis(
        ax, xticks=xticks, yticks=yticks,
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
        name = f'PO_SM_benchmark_AA_purification.{i}'
        file = os.path.join(images_folder, name)
        plt.savefig(file, dpi=900, transparent=True)
    return fig, all_axes, sms, pos

def profile(f, streams, adiabatic_stages, stages, total_time, init):
    time = bst.TicToc()
    KB_error = []
    flow_error = []
    energy_error = []
    material_error = []
    temperature_error = []
    record = []
    net_time = 0
    init()
    KBs = np.array([i.K * i.B for i in stages if hasattr(i, 'K')])
    temperatures = np.array([i.T for i in streams])
    flows = np.array([i.mol for i in streams])
    while net_time < total_time:
        time.tic()
        f()
        net_time += time.toc()
        new_KBs = np.array([i.K * i.B for i in stages if hasattr(i, 'K')])
        new_temperatures = np.array([i.T for i in streams])
        new_flows = np.array([i.mol for i in streams])
        record.append(net_time)
        dF = np.abs(flows - new_flows).sum()
        dT = np.abs(temperatures - new_temperatures).sum()
        KB_error.append(
            np.log10(np.abs(new_KBs - KBs).sum() + 1e-15)
        )
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
        KBs = new_KBs
        flows = new_flows
        temperatures = new_temperatures
    return {
        'Time': record, 
        'Stripping factor': KB_error,
        'Stream temperature': temperature_error,
        'Component flow rate': flow_error, 
        'Energy balance': energy_error, 
        'Material balance': material_error,
    }

# %%  Run

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
        systems=None, N=10, load=True, save=True
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
                sm = profile_sequential_modular(sys, time)
                po = profile_phenomena_oriented(sys, time)
                sms.append(sm)
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
        name = f'PO_SM_benchmark_AA_purification.{i}'
        file = os.path.join(images_folder, name)
        plt.savefig(file, dpi=900, transparent=True)
    return fig, all_axes, sms, pos
