# -*- coding: utf-8 -*-
"""
"""
import biosteam as bst
import numpy as np
from numpy.testing import assert_allclose

def test_trivial_lle_case():
    import biosteam as bst
    import numpy as np
    from numpy.testing import assert_allclose
    # lle
    with bst.System(algorithm='phenomena oriented') as sys:
        bst.settings.set_thermo(['Water', 'Octanol', 'Methanol'], cache=True)
        feed = bst.Stream(Water=50, Methanol= 5, units='kg/hr')
        solvent = bst.Stream(Octanol=50, units='kg/hr')
        stage = bst.StageEquilibrium(
            phases=('L', 'l'), ins=[feed, solvent], top_chemical='Octanol',
        )
    sys.simulate()
    streams = stage.ins + stage.outs
    actuals = [i.mol.copy() for i in streams]
    sys.run_phenomena()
    values = [i.mol for i in streams]
    for actual, value in zip(actuals, values):
        assert_allclose(actual, value)

def test_trivial_vle_case():
    import biosteam as bst
    import numpy as np
    from numpy.testing import assert_allclose
    # vle
    with bst.System(algorithm='phenomena oriented') as sys:
        bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
        liquid = bst.Stream(Water=50, Ethanol=10, units='kg/hr')
        vapor = bst.Stream(Ethanol=10, Water=50, phase='g', units='kg/hr')
        stage = bst.StageEquilibrium(phases=('g', 'l'), ins=[liquid, vapor])
    sys.simulate()
    streams = stage.ins + stage.outs
    actuals = [i.mol.copy() for i in streams]
    sys.run_phenomena()
    values = [i.mol for i in streams]
    for actual, value in zip(actuals, values):
        assert_allclose(actual, value)
 
def test_trivial_liquid_extraction_case():
    # liquid extraction
    with bst.System(algorithm='phenomena oriented') as sys:
        bst.settings.set_thermo(['Water', 'Methanol', 'Octanol'], cache=True)
        feed = bst.Stream(Water=500, Methanol=50)
        solvent = bst.Stream(Octanol=500, T=330)
        MSE = bst.MultiStageEquilibrium(N_stages=2, ins=[feed, solvent], phases=('L', 'l'))
    sys.simulate()
    extract, raffinate = MSE.outs
    actual = extract.imol['Methanol'] / feed.imol['Methanol']
    T_actual = extract.T
    sys.run_phenomena()
    value = extract.imol['Methanol'] / feed.imol['Methanol']
    T = extract.T
    assert_allclose(actual, value, rtol=0.01, atol=1e-3)
    assert_allclose(T_actual, T, rtol=0.01, atol=1e-3)

def test_trivial_distillation_case():    
    # distillation
    with bst.System(algorithm='phenomena oriented') as sys:
        bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
        feed = bst.Stream(Ethanol=80, Water=100, T=80.215 + 273.15)
        MSE = bst.MultiStageEquilibrium(N_stages=5, ins=[feed], feed_stages=[2],
            outs=['vapor', 'liquid'],
            stage_specifications={0: ('Reflux', 0.673), -1: ('Boilup', 2.57)},
            phases=('g', 'l'),
        )
    sys.simulate()
    vapor, liquid = MSE.outs
    actual = round(vapor.imol['Ethanol'] / feed.imol['Ethanol'], 2)
    sys.run_phenomena()
    assert round(vapor.imol['Ethanol'] / feed.imol['Ethanol'], 2) == actual
    
def test_simple_acetic_acid_separation_no_recycle():
    bst.settings.set_thermo(['Water', 'AceticAcid', 'EthylAcetate'], cache=True)
    @bst.SystemFactory
    def system(ins, outs):
        feed = bst.Stream(AceticAcid=6660, Water=43600)
        solvent = bst.Stream(EthylAcetate=65000)
        LE = bst.MultiStageEquilibrium(
            N_stages=6, ins=[feed, solvent], phases=('L', 'l'),
            maxiter=200,
            use_cache=True,
        )
        # DAA = bst.MultiStageEquilibrium(N_stages=6, ins=[LE-0], feed_stages=[3],
        #     outs=['vapor', 'liquid'],
        #     stage_specifications={0: ('Reflux', 0.673), -1: ('Boilup', 2.57)},
        #     phases=('g', 'l'),
        #     method='anderson',
        #     use_cache=True,
        # )
        DEA = bst.MultiStageEquilibrium(N_stages=6, ins=[LE-1], feed_stages=[3],
            outs=['vapor', 'liquid'],
            stage_specifications={0: ('Reflux', 0.673), -1: ('Boilup', 2.57)},
            phases=('g', 'l'),
            use_cache=True,
        )
    init_sys = system()
    init_sys.simulate()
    po = system(algorithm='phenomena oriented', 
                    molar_tolerance=1e-6,
                    relative_molar_tolerance=1e-6,
                    method='fixed-point')
    sm = system(algorithm='sequential modular',
                    molar_tolerance=1e-6,
                    relative_molar_tolerance=1e-6,
                    method='fixed-point')
    
    for i in range(2): 
        po.simulate()
        po.run_phenomena()
    for i in range(2): 
        sm.simulate()
    
    for s_sm, s_dp in zip(sm.streams, po.streams):
        actual = s_sm.mol
        value = s_dp.mol
        assert_allclose(actual, value, rtol=1e-3, atol=1e-3)

def test_simple_acetic_acid_separation_with_recycle():
    import biosteam as bst
    import numpy as np
    from numpy.testing import assert_allclose
    bst.settings.set_thermo(['Water', 'AceticAcid', 'EthylAcetate'], cache=True)
    solvent_feed_ratio = 1
    @bst.SystemFactory
    def system(ins, outs):
        feed = bst.Stream('feed', AceticAcid=6660, Water=43600)
        solvent = bst.Stream('solvent', EthylAcetate=65000)
        recycle = bst.Stream('recycle')
        LE = bst.MultiStageEquilibrium(
            N_stages=6, ins=[feed, solvent, recycle],
            feed_stages=(0, -1, -1),
            phases=('L', 'l'),
            maxiter=200,
            use_cache=True,
        )
        # DAA = bst.MultiStageEquilibrium(N_stages=6, ins=[LE-0], feed_stages=[3],
        #     outs=['vapor', 'liquid'],
        #     stage_specifications={0: ('Reflux', 0.673), -1: ('Boilup', 2.57)},
        #     maxiter=200,
        #     phases=('g', 'l'),
        #     use_cache=True,
        # )
        DEA = bst.MultiStageEquilibrium(N_stages=6, ins=[LE-1], feed_stages=[3],
            outs=['vapor', 'liquid'],
            stage_specifications={0: ('Reflux', 0.673), -1: ('Boilup', 2.57)},
            phases=('g', 'l'),
            maxiter=200,
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
        
    init_sys = system()
    init_sys.simulate()
    po = system(algorithm='phenomena oriented', 
                    molar_tolerance=1e-9,
                    relative_molar_tolerance=1e-9,
                    method='fixed-point')
    sm = system(algorithm='sequential modular',
                    molar_tolerance=1e-9,
                    relative_molar_tolerance=1e-9,
                    method='fixed-point')
    time = bst.TicToc()
    
    time.tic()
    for i in range(1): po.simulate()
    t_phenomena = time.toc()
    
    time.tic()
    for i in range(1): sm.simulate()
    t_sequential = time.toc()
    
    print('SM', t_sequential)
    print('PO', t_phenomena)
    
    for s_sm, s_dp in zip(sm.streams, po.streams):
        actual = s_sm.mol
        value = s_dp.mol
        assert_allclose(actual, value, rtol=1e-1, atol=1e-3)

# def test_ethanol_purification_system():
#     import biosteam as bst
#     import numpy as np
#     bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
#     feed = bst.Stream(
#         ID='feed', phase='l', T=373.83, P=607950, 
#         Water=1.11e+04, Ethanol=545.8, units='kmol/hr'
#     )

# def test_integrated_settler_acetic_acid_separation_system(): # integratted settler
#     from numpy.testing import assert_allclose
#     import biosteam as bst
#     import numpy as np
#     thermo = bst.Thermo(['Water', 'AceticAcid', 'EthylAcetate'], cache=True)
#     bst.settings.set_thermo(thermo)
#     def create_system(alg):
#         solvent_feed_ratio = 1.5
#         chemicals = bst.settings.chemicals
#         acetic_acid_broth = bst.Stream(
#             ID='acetic_acid_broth', AceticAcid=6660, Water=43600, units='kg/hr',
#         )
#         ethyl_acetate = bst.Stream(
#             ID='fresh_solvent', EthylAcetate=15000, units='kg/hr',
#         )
#         glacial_acetic_acid = bst.Stream(
#             'glacial_acetic_acid', 
#         )
#         wastewater = bst.Stream(
#             'wastewater',
#         )
#         solvent_recycle = bst.Stream(
#             'solvent_rich', 
#         )
#         reflux = bst.Stream('reflux')
#         water_rich = bst.Stream('water_rich')
#         distillate = bst.Stream('distillate')
#         distillate_2 = bst.Stream('distillate_2')
#         with bst.System(algorithm=alg) as sys:
#             # @ethyl_acetate.equation('material')
#             # def fresh_solvent_flow_rate():
#             #     f = np.ones(chemicals.size)
#             #     r = np.zeros(chemicals.size)
#             #     v = r.copy()
#             #     index = chemicals.index('EthylAcetate')
#             #     r[index] = 1
#             #     v[index] = solvent_feed_ratio * acetic_acid_broth.F_mass / chemicals.EthylAcetate.MW
#             #     return (
#             #         {ethyl_acetate: f,
#             #          ED-0: r,
#             #          distillate: r,
#             #          distillate_2: r},
#             #          v
#             #     )
#             @ethyl_acetate.equation('material')
#             def fresh_solvent_flow_rate():
#                 f = np.ones(chemicals.size)
#                 r = np.zeros(chemicals.size)
#                 v = r.copy()
#                 index = chemicals.index('EthylAcetate')
#                 r[index] = 1
#                 v[index] = solvent_feed_ratio * acetic_acid_broth.F_mass / chemicals.EthylAcetate.MW
#                 return (
#                     {ethyl_acetate: f,
#                       solvent_recycle: r},
#                       v
#                 )
#             extractor = bst.MultiStageMixerSettlers(
#                 'extractor', 
#                 ins=(acetic_acid_broth, solvent_recycle, ethyl_acetate), 
#                 outs=('extract', 'raffinate'),
#                 top_chemical='EthylAcetate',
#                 feed_stages=(0, -1, -1),
#                 N_stages=15,
#                 collapsed_init=False,
#                 use_cache=True,
#                 thermo=thermo,
#             )
#             # @extractor.add_specification(run=True)
#             # def run_settler_first():
#             #     settler.run()
#             @extractor.add_specification(run=True)
#             def adjust_fresh_solvent_flow_rate():
#                 broth = acetic_acid_broth.F_mass
#                 EtAc_recycle = solvent_recycle.imass['EthylAcetate']
#                 ethyl_acetate.imass['EthylAcetate'] = max(
#                     0, broth * solvent_feed_ratio - EtAc_recycle
#                 )
#             HX = bst.StageEquilibrium(
#                 'HX_extract',
#                 ins=[extractor.extract], 
#                 phases=('g', 'l'),
#                 B=1,
#             )
#             ED = bst.MESHDistillation(
#                 'extract_distiller',
#                 ins=[HX-0, HX-1, reflux],
#                 outs=['vapor', ''],
#                 LHK=('EthylAcetate', 'AceticAcid'),
#                 N_stages=15,
#                 feed_stages=(7, 7, 0),
#                 reflux=None,
#                 boilup=3,
#                 use_cache=True,
#                 # algorithm='optimize',
#                 # method='CG',
#             )
#             settler = bst.StageEquilibrium(
#                 'settler',
#                 ins=(ED-0, distillate, distillate_2), 
#                 outs=(solvent_recycle, water_rich, ''),
#                 phases=('L', 'l'),
#                 top_chemical='EthylAcetate',
#                 top_split=0.4,
#                 T=310,
#                 # partition_data={
#                 #     'K': np.array([ 0.253,  2.26 , 40.816]),
#                 #     'IDs': ('Water', 'AceticAcid', 'EthylAcetate'),
#                 # },
#                 thermo=thermo,
#             )
#             HX = bst.StageEquilibrium(
#                 'HX_reflux',
#                 ins=[settler-2], 
#                 outs=['', reflux],
#                 phases=('g', 'l'),
#                 B=0,
#             )
#             # @settler.add_specification(run=True)
#             # def adjust_fresh_solvent_flow_rate():
#             #     broth = acetic_acid_broth.F_mass
#             #     EtAc_recycle = sum([i.imass['EthylAcetate'] for i in (ED-0, distillate, distillate_2)])
#             #     ethyl_acetate.imass['EthylAcetate'] = max(
#             #         0, broth * solvent_feed_ratio - EtAc_recycle
#             #     )
#             # settler.coupled_KL = True
#             AD = bst.ShortcutColumn(
#                 'acetic_acid_distiller',
#                 LHK=('EthylAcetate', 'AceticAcid'),
#                 ins=ED-1,
#                 outs=[distillate_2, glacial_acetic_acid],
#                 partial_condenser=False,
#                 Lr=0.999,
#                 Hr=0.999,
#                 k=1.5,
#             )
#             HX = bst.StageEquilibrium(
#                 'HX',
#                 ins=[water_rich, extractor.raffinate], 
#                 phases=('g', 'l'),
#                 B=0,
#             )
#             AD.check_LHK = False
#             RD = bst.MESHDistillation(
#                 'raffinate_distiller',
#                 LHK=('EthylAcetate', 'Water'),
#                 ins=[HX-0, HX-1],
#                 outs=['', wastewater, distillate],
#                 full_condenser=True,
#                 N_stages=10,
#                 feed_stages=(1, 2),
#                 reflux=1,
#                 boilup=2,
#                 # algorithm='optimize',
#                 # method='CG',
#             )
#         return sys
#     time = bst.TicToc()
#     time.tic()
#     bst.F.set_flowsheet('SM')
#     sm = create_system('sequential modular')
#     sm.flatten()
#     sm.set_tolerance(
#         rmol=1e-3, mol=1e-3, subsystems=True,
#         method='fixed-point', maxiter=300,
#     )
#     sm.simulate()
#     t_sequential = time.toc()
#     # for i in sm.units: print(i, i.mass_balance_error())
#     # return
#     time.tic()
#     bst.F.set_flowsheet('PO')
#     po = create_system('phenomena oriented')
#     po.flatten()
#     po.set_tolerance(rmol=1e-3, mol=1e-3, 
#                      subsystems=True,
#                      method='fixed-point',
#                      maxiter=300)
#     po.simulate()
#     t_phenomena = time.toc()

#     print('SM', t_sequential)
#     print('PO', t_phenomena)

#     # TODO: Figure out why it doesn't match
#     # for s_sm, s_dp in zip(sm.streams, po.streams):
#     #     actual = s_sm.mol
#     #     value = s_dp.mol
#     #     assert_allclose(actual, value, rtol=0.2, atol=2)


if __name__ == '__main__':
    pass
    # Art Westerberg
    # Lorence Bigler
    test_trivial_lle_case()
    test_trivial_vle_case()
    test_trivial_liquid_extraction_case()
    test_trivial_distillation_case()
    test_simple_acetic_acid_separation_no_recycle()
    test_simple_acetic_acid_separation_with_recycle()
    test_integrated_settler_acetic_acid_separation_system()
    # test_complex_acetic_acid_separation_system()