# -*- coding: utf-8 -*-
"""
Created on Sat Mar 16 13:38:11 2024

@author: cortespea
"""
import biosteam as bst
import numpy as np

__all__ = (
    'create_acetic_acid_simple_system',
    'create_acetic_acid_complex_system',
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
    return sys

def create_acetic_acid_complex_system(
        alg, 
        extractor_stages=15, 
        raffinate_distillation_stages=10,
        extract_distillation_stages=15,
    ):
    thermo = bst.Thermo(['Water', 'AceticAcid', 'EthylAcetate'], cache=True)
    bst.settings.set_thermo(thermo)
    solvent_feed_ratio = 1.5
    chemicals = bst.settings.chemicals
    acetic_acid_broth = bst.Stream(
        ID='acetic_acid_broth', AceticAcid=6660, Water=43600, units='kg/hr',
        # T=310,
    )
    ethyl_acetate = bst.Stream(
        ID='fresh_solvent', EthylAcetate=15000, units='kg/hr', 
        # T=310,
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
    # bst.units.MultiStageEquilibrium.default_maxiter = 30
    # bst.units.MultiStageEquilibrium.default_attempts = 5
    # bst.units.MultiStageEquilibrium.default_molar_tolerance = 1e-9
    # bst.units.MultiStageEquilibrium.default_relative_molar_tolerance = 1e-9
    with bst.System(algorithm=alg) as sys:
        # @ethyl_acetate.material_balance
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
        @ethyl_acetate.material_balance
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
            N_stages=extractor_stages,
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
            EtAc_required = broth * solvent_feed_ratio
            if EtAc_required < EtAc_recycle:
                solvent_recycle.F_mass *= EtAc_required / EtAc_recycle
                EtAc_recycle = solvent_recycle.imass['EthylAcetate']
            EtAc_fresh = EtAc_required - EtAc_recycle
            ethyl_acetate.imass['EthylAcetate'] = max(
                0, EtAc_fresh
            )
            
        HX = bst.StageEquilibrium(
            'HX_extract',
            ins=[extractor.extract], 
            phases=('g', 'l'),
            B=1,
        )
        extract_feed_stage = int(extract_distillation_stages / 2)
        ED = bst.MESHDistillation(
            'extract_distillation',
            ins=[HX-0, HX-1, reflux],
            outs=['vapor', ''],
            LHK=('EthylAcetate', 'AceticAcid'),
            N_stages=extract_distillation_stages,
            feed_stages=(extract_feed_stage, extract_feed_stage, 0),
            
            reflux=None,
            boilup=3,
            use_cache=True,
        )
        mixer = bst.Mixer('distillation_recycle', ins=[ED-0, distillate_2])
        settler = bst.StageEquilibrium(
            'settler',
            ins=(mixer-0, distillate), 
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
            'acetic_acid_distillation',
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
            'raffinate_distillation',
            LHK=('EthylAcetate', 'Water'),
            ins=[HX-0, HX-1],
            outs=['', wastewater, distillate],
            full_condenser=True,
            N_stages=raffinate_distillation_stages,
            feed_stages=(1, 2),
            reflux=1,
            boilup=2,
            use_cache=True,
        )
    return sys

def test_simple_acetic_acid_purification_system():
    pass

def test_complex_acetic_acid_purification_system():
    po = create_acetic_acid_complex_system('phenomena-oriented')
    po.flatten()
    po.set_tolerance(mol=1e-3, rmol=1e-3, maxiter=20)
    sm = create_acetic_acid_complex_system('sequential modular')
    sm.flatten()
    sm.set_tolerance(mol=1e-3, rmol=1e-3, maxiter=20)
    po.simulate()
    