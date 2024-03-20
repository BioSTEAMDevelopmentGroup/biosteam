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
            ins=[HX-1, reflux],
            outs=['vapor', ''],
            LHK=('EthylAcetate', 'AceticAcid'),
            N_stages=10,
            feed_stages=(5, 0),
            reflux=None,
            boilup=4.37,
            use_cache=True,
        )
        settler = bst.StageEquilibrium(
            'settler',
            ins=(ED-0, distillate, distillate_2), 
            outs=(solvent_recycle, water_rich, ''),
            phases=('L', 'l'),
            top_chemical='EthylAcetate',
            top_split=0.1,
            T=340,
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
        HX = bst.StageEquilibrium(
            'HX_raffinate',
            ins=[extractor.raffinate, water_rich], 
            phases=('g', 'l'),
            B=0,
        )
        AD = bst.ShortcutColumn(
            'acetic_acid_distiller',
            LHK=('Water', 'AceticAcid'),
            ins=ED-1,
            outs=[distillate_2, glacial_acetic_acid],
            partial_condenser=False,
            Lr=0.999,
            Hr=0.999,
            k=1.5,
        )
        RD = bst.MESHDistillation(
            'raffinate_distiller',
            LHK=('EthylAcetate', 'Water'),
            ins=[HX-1],
            outs=['', wastewater, distillate],
            full_condenser=True,
            N_stages=15,
            feed_stages=(2,),
            reflux=0.1,
            boilup=1,
        )
    return sys

def profile_sequential_modular(total_time, kind='simple', N=None):
    if kind == 'simple':
        sm = create_system_simple('sequential modular')
    else:
        sm = create_system_complex('sequential modular')
    sm.flatten()
    sm._setup_units()
    data = profile(sm.run_sequential_modular, *streams_and_stages(sm), total_time)
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
                   po.run_sequential_modular)
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

def profile(f, streams, adiabatic_stages, stages, total_time, g=None, kind=None):
    time = bst.TicToc()
    recycle_values = np.array([i.mol for i in streams])
    recycle_temperatures = np.array([i.T for i in streams])
    recycle_error = []
    energy_error = []
    material_error = []
    temperature_error = []
    record = []
    net_time = 0
    while net_time < total_time:
        time.tic()
        f()
        net_time += time.toc()
        record.append(net_time)
        new_recycle_temperatures = np.array([i.T for i in streams])
        new_recycle_values = np.array([i.mol for i in streams])
        recycle_error.append(
            np.log10(np.abs(recycle_values - new_recycle_values).sum() + 1e-15)
        )
        temperature_error.append(
            np.log10(np.abs(recycle_temperatures - new_recycle_temperatures).sum() + 1e-15)
        )
        recycle_values = new_recycle_values
        recycle_temperatures = new_recycle_temperatures
        energy_error.append(
            np.log10(sum([abs(dT_error(i)) for i in adiabatic_stages]) + 1e-15)
        )
        material_error.append(
            np.log10(sum([abs(i.mass_balance_error()) for i in stages]) + 1e-15)
        )
    return {
        'Time': record, 
        'Stream temperature': temperature_error,
        'Component flow rate': recycle_error, 
        'Energy balance': energy_error, 
        'Material balance': material_error,
    }

cache = {}

# %%  Run

def dct_mean(dcts: list[dict], keys: list[str], ub: float, n=50):
    N = len(dcts)
    mean = {i: np.zeros(n) for i in keys}
    tmin = np.min([np.min(i['Time']) for i in dcts])
    t = np.linspace(0, ub, n)
    for dct in dcts:
        for i in keys: 
            x = dct['Time']
            mean['Time'] = t
            values = interpolate.interp1d(x, dct[i], bounds_error=False)(t)
            mask = ~np.isnan(values)
            mean[i][mask] += values[mask]
            mean[i][t < tmin] = np.nan
    for i in keys: mean[i] /= N
    return mean

def plot_profile(kinds=('simple', 'complex'), times=[10, 40], N=4):
    bst.set_font()
    bst.set_figure_size(7)
    fig, all_axes = plt.subplots(4, 2)
    keys = (
        'Component flow rate',
        'Material balance',
        'Stream temperature',
        'Energy balance',
    )
    units = (
        r'$[\mathrm{log}_{\mathrm{10}}\ \mathrm{mol} \cdot \mathrm{hr}^{\mathrm{-1}}]$',
        r'$[\mathrm{log}_{\mathrm{10}}\ \mathrm{mol} \cdot \mathrm{hr}^{\mathrm{-1}}]$',
        r'$[\mathrm{log}_{\mathrm{10}}\ \mathrm{K}]$',
        r'$[\mathrm{log}_{\mathrm{10}}\ \mathrm{K}]$',
    )
    for m, (kind, time) in enumerate(zip(kinds, times)):
        axes = all_axes[:, m]
        key = (kind, time, N)
        if key in cache:
            sms, pos = cache[key]
        else:
            # profile_sequential_modular(0.1, kind=kind)
            # profile_phenomena_oriented(0.1, kind=kind)
            sms = [profile_sequential_modular(time, kind=kind) for i in range(N)]
            pos = [profile_phenomena_oriented(time, kind=kind) for i in range(N)]
            cache[key] = (sms, pos)
        if kind == 'simple':
            tub = 10
        else:
            tub = 40
        sm = dct_mean(sms, keys, tub)
        po = dct_mean(pos, keys, tub)
        csm = Color(fg='#33BBEE').RGBn
        cpo = Color(fg='#EE7733').RGBn
        yticklabels = m == 0
        for n, (i, ax, u) in enumerate(zip(keys, axes, units)):
            plt.sca(ax)
            if n == 3: plt.xlabel('Time [s]')
            if i == 'Stream temperature':
                label = 'Temperature'
            elif i == 'Component flow rate':
                label = 'Flow rate'
            else:
                label = i
            if m == 0: plt.ylabel(f'{label}\n{u}')
            # for sm, po in zip(sms, pos):
            ysm = np.array(sm[i])
            ypo = np.array(po[i])
            tsm = np.array(sm['Time'])
            tpo = np.array(po['Time'])
            # ysm[ysm < -10] = -10
            # ypo[ypo < -10] = -10
            # plt.plot(xsm, ysm, lw=0, marker='s', color=csm, markersize=2.5)
            plt.plot(tsm, ysm, '-', color=csm, lw=1.5)
            # plt.plot(xpo, ypo, lw=0, marker='d', color=cpo, markersize=2.5)
            plt.plot(tpo, ypo, '-', color=cpo, lw=1.5)
            if kind == 'simple':
                step = 1
            else:
                step = tub // 10
            xticks = [*range(0, tub + 1, step)]
            xticklabels = xtick0 = n == 3
            ytick0 = n == 3
            plt.xlim(0, xticks[-1])
            plt.ylim(-16, 8)
            bst.utils.style_axis(
                ax, xticks=xticks, yticks=[*range(-16, 8, 6)], 
                xtick0=xtick0, ytick0=ytick0, xticklabels=xticklabels,
                yticklabels=yticklabels,
            )
    plt.subplots_adjust(hspace=0, wspace=0.15)
    return fig, all_axes, sms, pos

fig, axes, sm, po = plot_profile()