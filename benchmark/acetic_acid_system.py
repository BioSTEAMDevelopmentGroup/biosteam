# -*- coding: utf-8 -*-
"""
Created on Sat Mar 16 13:38:11 2024

@author: cortespea
"""
import biosteam as bst
import numpy as np
from .profile import register

__all__ = (
    'create_acetic_acid_simple_system',
    'create_acetic_acid_complex_system',
)

@register(
    'acetic_acid_simple', 'Acetic acid\ndistillation & liquid extraction',
    10, [0, 2, 4, 6, 8, 10]
)
def create_acetic_acid_simple_system(alg):
    solvent_feed_ratio = 1
    thermo = bst.Thermo(['Water', 'AceticAcid', 'EthylAcetate'], cache=True)
    bst.settings.set_thermo(thermo)
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

@register(
    'acetic_acid_complex', 'Glacial acetic acid\npurification',
    80, [0, 10, 20, 30, 40, 50, 60, 70, 80],
)
def create_acetic_acid_complex_system(alg):
    thermo = bst.Thermo(['Water', 'AceticAcid', 'EthylAcetate'], cache=True)
    bst.settings.set_thermo(thermo)
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

