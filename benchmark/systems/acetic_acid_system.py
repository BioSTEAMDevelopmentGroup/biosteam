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
    'create_acetic_acid_complex_decoupled_system',
)

def create_acetic_acid_simple_system(alg, ideal=False):
    solvent_feed_ratio = 1
    IDs = ['Water', 'AceticAcid', 'EthylAcetate']
    bst.settings.set_thermo(IDs, cache=True, ideal=ideal)
    with bst.System('acetic_acid_simple', algorithm=alg) as sys:
        feed = bst.Stream('feed', AceticAcid=6660, Water=43600)
        solvent = bst.Stream('solvent', EthylAcetate=65000)
        solvent.F_mass = feed.F_mass * solvent_feed_ratio
        recycle = bst.Stream('recycle')
        LE = bst.MultiStageEquilibrium(
            N_stages=5, ins=[feed, solvent, recycle],
            feed_stages=(0, -1, -1),
            phases=('L', 'l'),
            maxiter=200,
            use_cache=True,
            method='fixed-point',
            partition_data={
                'IDs': IDs,
                'K': np.array([
                    [0.314,  1.787, 12.835],
                    [ 0.3  ,  1.841, 14.721],
                    [ 0.285,  1.906, 17.297],
                    [ 0.271,  1.985, 20.962],
                    [ 0.258,  2.083, 26.435]
                ]),
            } if ideal else None,
        )
        # DAA = bst.MultiStageEquilibrium(N_stages=6, ins=[LE-0], feed_stages=[3],
        #     outs=['vapor', 'liquid'],
        #     stage_specifications={0: ('Reflux', 0.673), -1: ('Boilup', 2.57)},
        #     maxiter=200,
        #     phases=('g', 'l'),
        #     method='fixed-point',
        #     use_cache=True,
        # )
        DEA = bst.MESHDistillation(N_stages=6, ins=[LE-1], feed_stages=[3],
            outs=['', 'bottoms', recycle],
            # stage_specifications={0: ('Reflux', 0.673), -1: ('Boilup', 2.57)},
            reflux=0.673,
            boilup=2.57,
            full_condenser=True,
            maxiter=200,
            method='fixed-point',
            use_cache=True,
            LHK=('EthylAcetate', 'AceticAcid'),
        )
        # HX = bst.SinglePhaseStage(ins=DEA-0, outs=recycle, T=320, phase='l')
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
        ideal=False
    ):
    IDs = ['Water', 'AceticAcid', 'EthylAcetate']
    bst.settings.set_thermo(IDs, cache=True, ideal=ideal)
    solvent_feed_ratio = 1.0
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
    with bst.System('acetic_acid_complex', algorithm=alg) as sys:
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
            partition_data={
                'IDs': IDs,
                'K': np.array([
                    [ 0.198,  2.109, 31.506],
                    [ 0.198,  2.109, 31.511],
                    [ 0.198,  2.109, 31.52 ],
                    [ 0.198,  2.109, 31.533],
                    [ 0.198,  2.109, 31.555],
                    [ 0.198,  2.11 , 31.592],
                    [ 0.198,  2.111, 31.652],
                    [ 0.198,  2.113, 31.753],
                    [ 0.198,  2.115, 31.922],
                    [ 0.199,  2.12 , 32.205],
                    [ 0.199,  2.128, 32.683],
                    [ 0.201,  2.141, 33.49 ],
                    [ 0.204,  2.163, 34.845],
                    [ 0.21 ,  2.198, 37.094],
                    [ 0.222,  2.248, 40.681]
                ]),
            } if ideal else None,
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
        extract_feed_stage = int(extract_distillation_stages / 2)
        ED = bst.MESHDistillation(
            'extract_distillation',
            ins=[extractor.extract, reflux],
            outs=['', '', ''],
            LHK=('EthylAcetate', 'AceticAcid'),
            N_stages=extract_distillation_stages,
            feed_stages=(extract_feed_stage, 1),
            full_condenser=True,
            boilup=3,
            use_cache=True,
        )
        # ED.fallback = ('sequential modular', 'fixed-point', 5)
        mixer = bst.Mixer('distillation_recycle', ins=[ED-2, distillate, distillate_2])
        settler = bst.StageEquilibrium(
            'settler',
            ins=(mixer-0), 
            outs=(solvent_recycle, water_rich, reflux),
            phases=('L', 'l'),
            top_chemical='EthylAcetate',
            top_split=0.7,
            partition_data={
                'IDs': IDs,
                'K': np.array([ 0.248,  2.301, 45.828]),
            } if ideal else None,
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
            ins=[HX-1],
            outs=['', wastewater, distillate],
            full_condenser=True,
            N_stages=raffinate_distillation_stages,
            feed_stages=(2,),
            reflux=1,
            boilup=2,
            use_cache=True,
        )
        # RD.fallback = ('sequential modular', 'fixed-point', 5)
        # RD = bst.ShortcutColumn(
        #     'raffinate_distillation',
        #     LHK=('EthylAcetate', 'Water'),
        #     ins=HX-1,
        #     outs=[distillate, wastewater],
        #     partial_condenser=False,
        #     k=1.4,
        #     Lr=0.999,
        #     Hr=0.95,
        # )
    return sys

def create_acetic_acid_complex_decoupled_system(
        alg, 
        extractor_stages=15, 
        ideal=False
    ):
    IDs = ['Water', 'AceticAcid', 'EthylAcetate']
    bst.settings.set_thermo(IDs, cache=True, ideal=ideal)
    solvent_feed_ratio = 1.0
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
    raffinate_distillate = bst.Stream('raffinate_distillate')
    extract_distillate = bst.Stream('extract_distillate')
    # bst.units.MultiStageEquilibrium.default_maxiter = 30
    # bst.units.MultiStageEquilibrium.default_attempts = 5
    # bst.units.MultiStageEquilibrium.default_molar_tolerance = 1e-9
    # bst.units.MultiStageEquilibrium.default_relative_molar_tolerance = 1e-9
    with bst.System('acetic_acid_complex_decoupled', algorithm=alg) as sys:
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
            partition_data={
                'IDs': IDs,
                'K': np.array([
                    [ 0.198,  2.109, 31.51 ],
                    [ 0.198,  2.109, 31.52 ],
                    [ 0.198,  2.109, 31.535],
                    [ 0.198,  2.109, 31.559],
                    [ 0.198,  2.11 , 31.596],
                    [ 0.198,  2.111, 31.656],
                    [ 0.198,  2.112, 31.75 ],
                    [ 0.198,  2.114, 31.902],
                    [ 0.198,  2.118, 32.146],
                    [ 0.198,  2.124, 32.545],
                    [ 0.198,  2.133, 33.2  ],
                    [ 0.199,  2.149, 34.287],
                    [ 0.201,  2.175, 36.107],
                    [ 0.205,  2.217, 39.17 ],
                    [ 0.215,  2.279, 44.289]
                ]),
            } if ideal else None,
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
        ED = bst.ShortcutColumn(
            'extract_distillation',
            ins=extractor.extract,
            outs=[extract_distillate, glacial_acetic_acid],
            LHK=('EthylAcetate', 'AceticAcid'),
            partial_condenser=False,
            Lr=0.999,
            Hr=0.9999,
            k=1.4,
        )
        ED.check_LHK = False
        mixer = bst.Mixer('distillate_mixer', ins=[extract_distillate, raffinate_distillate], outs=solvent_recycle)
        
        # @settler.add_specification(run=True)
        # def adjust_fresh_solvent_flow_rate():
        #     broth = acetic_acid_broth.F_mass
        #     EtAc_recycle = sum([i.imass['EthylAcetate'] for i in (ED-0, distillate, distillate_2)])
        #     ethyl_acetate.imass['EthylAcetate'] = max(
        #         0, broth * solvent_feed_ratio - EtAc_recycle
        #     )
        # settler.coupled_KL = True
        RD = bst.ShortcutColumn(
            'raffinate_distillation',
            LHK=('EthylAcetate', 'Water'),
            ins=extractor.raffinate,
            outs=[raffinate_distillate, wastewater],
            partial_condenser=False,
            k=1.4,
            Lr=0.999,
            Hr=0.95,
        )
    return sys

def test_simple_acetic_acid_purification_system():
    po = create_acetic_acid_simple_system('phenomena-based')
    po.set_tolerance(mol=1e-3, rmol=1e-3, maxiter=20)
    sm = create_acetic_acid_simple_system('sequential modular')
    sm.set_tolerance(mol=1e-3, rmol=1e-3, maxiter=20)
    po.simulate()
    sm.simulate()

def test_complex_acetic_acid_purification_system():
    # po = create_acetic_acid_complex_system('phenomena-based')
    # po.set_tolerance(mol=1e-3, rmol=1e-3, maxiter=20)
    sm = create_acetic_acid_complex_system('sequential modular')
    sm.set_tolerance(mol=1e-3, rmol=1e-3, maxiter=20)
    # po.simulate()
    sm.simulate()
    
def test_complex_decoupled_acetic_acid_purification_system():
    po = create_acetic_acid_complex_decoupled_system('phenomena-oriented')
    po.set_tolerance(mol=1e-3, rmol=1e-3, maxiter=20)
    sm = create_acetic_acid_complex_decoupled_system('sequential modular')
    sm.set_tolerance(mol=1e-3, rmol=1e-3, maxiter=20)
    po.simulate()
    sm.simulate()
    
# if __name__ == "__main__":
#     test_simple_acetic_acid_purification_system()
#     test_complex_acetic_acid_purification_system()