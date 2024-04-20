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
import os
import pickle

try:
    images_folder = os.path.join(os.path.dirname(__file__), 'images')
    simulations_folder = os.path.join(os.path.dirname(__file__), 'simulations')
except:
    images_folder = os.path.join(os.getcwd(), 'images')
    simulations_folder = os.path.join(os.getcwd(), 'simulations')

thermo = bst.Thermo(['Water', 'AceticAcid', 'EthylAcetate'], cache=True)
bst.settings.set_thermo(thermo)

def create_system_simple(alg):
    solvent_feed_ratio = 1
    with bst.System(algorithm=alg) as sys:
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
            solvent.imass['EthylAcetate'] = max(
                0, broth * solvent_feed_ratio - EtAc_recycle
            )

        @solvent.equation('material')
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
    return sys

def create_system_complex(alg):
    solvent_feed_ratio = 1.5
    chemicals = bst.settings.chemicals
    acetic_acid_broth = bst.Stream(
        ID='acetic_acid_broth', AceticAcid=6660, Water=43600, units='kg/hr',
    )
    ethyl_acetate = bst.Stream(
        ID='fresh_solvent', EthylAcetate=15000, units='kg/hr',
    )
    glacial_acetic_acid = bst.Stream(
        'glacial_acetic_acid', 
    )
    wastewater = bst.Stream(
        'wastewater',
    )
    solvent_recycle = bst.Stream(
        'solvent_rich', 
    )
    reflux = bst.Stream('reflux')
    water_rich = bst.Stream('water_rich')
    distillate = bst.Stream('distillate')
    distillate_2 = bst.Stream('distillate_2')
    with bst.System(algorithm=alg) as sys:
        # @ethyl_acetate.equation('material')
        # def fresh_solvent_flow_rate():
        #     f = np.ones(chemicals.size)
        #     r = np.zeros(chemicals.size)
        #     v = r.copy()
        #     index = chemicals.index('EthylAcetate')
        #     r[index] = 1
        #     v[index] = solvent_feed_ratio * acetic_acid_broth.F_mass / chemicals.EthylAcetate.MW
        #     return (
        #         {ethyl_acetate: f,
        #          ED-0: r,
        #          distillate: r,
        #          distillate_2: r},
        #          v
        #     )
        @ethyl_acetate.equation('material')
        def fresh_solvent_flow_rate():
            f = np.ones(chemicals.size)
            r = np.zeros(chemicals.size)
            v = r.copy()
            index = chemicals.index('EthylAcetate')
            r[index] = 1
            v[index] = solvent_feed_ratio * acetic_acid_broth.F_mass / chemicals.EthylAcetate.MW
            return (
                {ethyl_acetate: f,
                  solvent_recycle: r},
                  v
            )
        extractor = bst.MultiStageMixerSettlers(
            'extractor', 
            ins=(acetic_acid_broth, solvent_recycle, ethyl_acetate), 
            outs=('extract', 'raffinate'),
            top_chemical='EthylAcetate',
            feed_stages=(0, -1, -1),
            N_stages=15,
            collapsed_init=False,
            use_cache=True,
            thermo=thermo,
        )
        # @extractor.add_specification(run=True)
        # def run_settler_first():
        #     settler.run()
        @extractor.add_specification(run=True)
        def adjust_fresh_solvent_flow_rate():
            broth = acetic_acid_broth.F_mass
            EtAc_recycle = solvent_recycle.imass['EthylAcetate']
            ethyl_acetate.imass['EthylAcetate'] = max(
                0, broth * solvent_feed_ratio - EtAc_recycle
            )
        HX = bst.StageEquilibrium(
            'HX_extract',
            ins=[extractor.extract], 
            phases=('g', 'l'),
            B=1,
        )
        ED = bst.MESHDistillation(
            'extract_distiller',
            ins=[HX-0, HX-1, reflux],
            outs=['vapor', ''],
            LHK=('EthylAcetate', 'AceticAcid'),
            N_stages=15,
            feed_stages=(7, 7, 0),
            reflux=None,
            boilup=3,
            use_cache=True,
        )
        settler = bst.StageEquilibrium(
            'settler',
            ins=(ED-0, distillate, distillate_2), 
            outs=(solvent_recycle, water_rich, ''),
            phases=('L', 'l'),
            top_chemical='EthylAcetate',
            top_split=0.4,
            T=310,
            # partition_data={
            #     'K': np.array([ 0.253,  2.26 , 40.816]),
            #     'IDs': ('Water', 'AceticAcid', 'EthylAcetate'),
            # },
            thermo=thermo,
        )
        HX = bst.StageEquilibrium(
            'HX_reflux',
            ins=[settler-2], 
            outs=['', reflux],
            phases=('g', 'l'),
            B=0,
        )
        # @settler.add_specification(run=True)
        # def adjust_fresh_solvent_flow_rate():
        #     broth = acetic_acid_broth.F_mass
        #     EtAc_recycle = sum([i.imass['EthylAcetate'] for i in (ED-0, distillate, distillate_2)])
        #     ethyl_acetate.imass['EthylAcetate'] = max(
        #         0, broth * solvent_feed_ratio - EtAc_recycle
        #     )
        # settler.coupled_KL = True
        AD = bst.ShortcutColumn(
            'acetic_acid_distiller',
            LHK=('EthylAcetate', 'AceticAcid'),
            ins=ED-1,
            outs=[distillate_2, glacial_acetic_acid],
            partial_condenser=False,
            Lr=0.999,
            Hr=0.999,
            k=1.5,
        )
        HX = bst.StageEquilibrium(
            'HX',
            ins=[water_rich, extractor.raffinate], 
            phases=('g', 'l'),
            B=0,
        )
        AD.check_LHK = False
        RD = bst.MESHDistillation(
            'raffinate_distiller',
            LHK=('EthylAcetate', 'Water'),
            ins=[HX-0, HX-1],
            outs=['', wastewater, distillate],
            full_condenser=True,
            N_stages=10,
            feed_stages=(1, 2),
            reflux=1,
            boilup=2,
        )
    return sys

def profile_sequential_modular(total_time, kind='simple', N=None):
    if kind == 'simple':
        sm = create_system_simple('sequential modular')
    else:
        sm = create_system_complex('sequential modular')
    sm.flatten()
    sm._setup_units()
    data = profile(sm.run_sequential_modular, *streams_and_stages(sm), total_time, init=sm.run_sequential_modular)
    bst.F.clear()
    return data

def profile_phenomena_oriented(total_time, kind='simple', consolidated=None):
    if kind == 'simple':
        po = create_system_simple('phenomena oriented')
    else:
        po = create_system_complex('phenomena oriented')
        for i in po.units: i.consolidated = consolidated
    po.flatten()
    po._setup_units()
    data = profile(po.run_phenomena, *streams_and_stages(po), total_time,
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
            sum([i.H + i.Hf for i in stage.outs]) - sum([i.H + i.Hf for i in stage.ins])
        ) / sum([i.C for i in stage.outs])

def profile(f, streams, adiabatic_stages, stages, total_time, init=None):
    time = bst.TicToc()
    KB_error = []
    recycle_error = []
    energy_error = []
    material_error = []
    temperature_error = []
    record = []
    net_time = 0
    init()
    KBs = np.array([i.K * i.B for i in stages if hasattr(i, 'K')])
    recycle_temperatures = np.array([i.T for i in streams])
    recycle_values = np.array([i.mol for i in streams])
    while net_time < total_time:
        time.tic()
        f()
        net_time += time.toc()
        new_KBs = np.array([i.K * i.B for i in stages if hasattr(i, 'K')])
        new_recycle_temperatures = np.array([i.T for i in streams])
        new_recycle_values = np.array([i.mol for i in streams])
        record.append(net_time)
        KB_error.append(
            np.log(np.abs(new_KBs - KBs).sum())
        )
        recycle_error.append(
            np.log10(np.abs(recycle_values - new_recycle_values).sum() + 1e-15)
        )
        temperature_error.append(
            np.log10(np.abs(recycle_temperatures - new_recycle_temperatures).sum() + 1e-15)
        )
        energy_error.append(
            np.log10(sum([abs(dT_error(i)) for i in adiabatic_stages]) + 1e-15)
        )
        material_error.append(
            np.log10(sum([abs(i.mass_balance_error()) for i in stages]) + 1e-15)
        )
        KBs = new_KBs
        recycle_values = new_recycle_values
        recycle_temperatures = new_recycle_temperatures
    return {
        'Time': record, 
        'Stripping factor': KB_error,
        'Stream temperature': temperature_error,
        'Component flow rate': recycle_error, 
        'Energy balance': energy_error, 
        'Material balance': material_error,
    }

cache = {}

# %%  Run


def dct_mean(dcts: list[dict], keys: list[str], ub: float, n=50):
    mean = {i: np.zeros(n) for i in keys}
    tmin = np.min([np.min(i['Time']) for i in dcts])
    t = np.linspace(0, ub, n)
    goods = [np.zeros(n) for i in range(len(keys))]
    for dct in dcts:
        for i, good in zip(keys, goods): 
            x = dct['Time']
            mean['Time'] = t
            values = interpolate.interp1d(x, dct[i], bounds_error=False)(t)
            mask = ~np.isnan(values)
            mean[i][mask] += values[mask]
            mean[i][t < tmin] = np.nan
            good[t > tmin] += 1
    for i, j in zip(keys, goods): mean[i][j > 0] /= j[j > 0]
    return mean

def plot_profile(kinds=('simple', 'complex',), times=[10, 100], N=10, load=True, save=True):
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
    if n_rows == 4:
        bst.set_figure_size(aspect_ratio=1.1)
    elif n_rows == 2:
        bst.set_figure_size(aspect_ratio=0.6)
    else:
        bst.set_figure_size(aspect_ratio=1)
    fig, all_axes = plt.subplots(n_rows, 2)
    yticks = [-15, -10, -5, 0, 5, 10]
    for m, (kind, time) in enumerate(zip(kinds, times)):
        axes = all_axes[:, m]
        sms = []
        pos = []
        if load:
            sm_name = f'sm_{time}_{kind}.npy'
            file = os.path.join(simulations_folder, sm_name)
            with open(file, 'rb') as f: sms = pickle.load(f)
            po_name = f'po_{time}_{kind}.npy'
            file = os.path.join(simulations_folder, po_name)
            with open(file, 'rb') as f: pos = pickle.load(f)
        else:
            for i in range(N):
                key = (kind, time, i)
                if key in cache:
                    sm, po = cache[key]
                else:
                    sm = profile_sequential_modular(time, kind=kind)
                    po = profile_phenomena_oriented(time, kind=kind)
                    cache[key] = (sm, po)
                sms.append(sm)
                pos.append(po)
            if save:
                sm_name = f'sm_{time}_{kind}.npy'
                file = os.path.join(simulations_folder, sm_name)
                with open(file, 'wb') as f: pickle.dump(sms, f)
                po_name = f'po_{time}_{kind}.npy'
                file = os.path.join(simulations_folder, po_name)
                with open(file, 'wb') as f: pickle.dump(pos, f)
        if kind == 'simple':
            tub = 10
        else:
            tub = 80
        sm = dct_mean(sms, keys, tub)
        po = dct_mean(pos, keys, tub)
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
            ysm[ysm < -10] = -10
            ypo[ypo < -10] = -10
            # plt.plot(xsm, ysm, lw=0, marker='s', color=csm, markersize=2.5)
            plt.plot(tsm, ysm, '--', color=csm, lw=1.5)
            # plt.plot(xpo, ypo, lw=0, marker='d', color=cpo, markersize=2.5)
            plt.plot(tpo, ypo, '-', color=cpo, lw=1.5)
            if m == 0 and n == 0:
                index = int(len(tsm) * 0.55)
                xy = x, y = (tsm[index], ysm[index])
                ax.annotate('Sequential\nmodular',
                    xy=xy, 
                    xytext=(x+1, y+1),
                    arrowprops=dict(arrowstyle="->", color=csm),
                    color=csm,
                    fontsize=fs,
                    fontweight='bold',
                )
                index = int(len(tpo) * 0.45)
                xy = x, y = (tpo[index], ypo[index])
                ax.annotate('Phenomena\noriented',
                    xy=xy, 
                    xytext=(x-1, y-5),
                    arrowprops=dict(arrowstyle="->", color=cpo),
                    ha='right',
                    color=cpo,
                    fontsize=fs,
                    fontweight='bold',
                )
            if kind == 'simple':
                xticks = [i for i in range(0, tub + 1, 2)]
            else:
                xticks = [*range(0, tub + 1, 10)]
            xticklabels = xtick0 = n == n_rows-1
            xtickf = m == 1
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
    titles = ['Coupled acetic acid\ndistillation & liquid extraction', 'Glacial acetic acid\npurification']
    letter_color = c.neutral.shade(25).RGBn
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

fig, axes, sm, po = plot_profile()