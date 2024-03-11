# -*- coding: utf-8 -*-
"""
"""
import biosteam as bst
import numpy as np
from numpy.testing import assert_allclose

# @bst.SystemFactory(
#     ins=(),
#     outs=('vapor', 'liquid'),
# )
# def distillation(ins, outs, N_stages, R, B, feed_stages):
#     top, bottom = outs
#     midins = list(ins)
#     feed_stages = list(feed_stages)
#     if R is None:
#         N_stages -= 1
#         feed_stages = [i - 1 for i in feed_stages]
#         condenser_feed = top
#     else:
#         reflux = bst.Stream()
#         condenser_feed = bst.Stream()
#         bst.StageEquilibrium(
#             '', ins=condenser_feed, outs=[top, reflux], B=1/R, phases=('g', 'l')
#         )
#         midins.append(reflux)
#         feed_stages.append(0)
#     if B is None:
#         N_stages -= 1
#     else:
#         reboiler_feed = bst.Stream()
#         boilup = bst.Stream()
#         bst.StageEquilibrium(
#             '', ins=reboiler_feed, outs=[boilup, bottom], B=B, phases=('g', 'l')
#         )
#         midins.append(boilup)
#         feed_stages.append(-1)
#     bst.MultiStageEquilibrium(
#         '', N_stages=N_stages, ins=midins, feed_stages=feed_stages,
#         outs=[],
#         phases=('g', 'l'),
#         method='fixed-point',
#         use_cache=True,
#     )

def test_trivial_lle_case():
    # lle
    with bst.System(algorithm='phenomena oriented') as sys:
        bst.settings.set_thermo(['Water', 'Octanol', 'Methanol'], cache=True)
        feed = bst.Stream(Water=50, Methanol= 5, units='kg/hr')
        solvent = bst.Stream(Octanol=50, units='kg/hr')
        stage = bst.StageEquilibrium(phases=('L', 'l'), ins=[feed, solvent])
    sys.simulate()
    streams = stage.ins + stage.outs
    actuals = [i.mol.copy() for i in streams]
    sys.run_phenomena()
    values = [i.mol for i in streams]
    for actual, value in zip(actuals, values):
        assert_allclose(actual, value)

def test_trivial_vle_case():
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
    assert_allclose(actual, value, rtol=0.01)
    assert_allclose(T_actual, T, rtol=0.01)

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
            method='fixed-point',
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
            method='fixed-point',
            use_cache=True,
        )
    init_sys = system()
    init_sys.simulate()
    dp_sys = system(algorithm='phenomena oriented', 
                    molar_tolerance=1e-6,
                    relative_molar_tolerance=1e-6)
    sm_sys = system(algorithm='sequential modular',
                    molar_tolerance=1e-6,
                    relative_molar_tolerance=1e-6,
                    method='wegstein')
    
    for i in range(2): 
        dp_sys.simulate()
        dp_sys.run_phenomena()
    for i in range(2): 
        sm_sys.simulate()
    
    for s_sm, s_dp in zip(sm_sys.streams, dp_sys.streams):
        actual = s_sm.mol
        value = s_dp.mol
        assert_allclose(actual, value, rtol=0.01)

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
        
    init_sys = system()
    init_sys.simulate()
    dp_sys = system(algorithm='phenomena oriented', 
                    molar_tolerance=1e-6,
                    relative_molar_tolerance=1e-6,
                    method='fixed-point')
    sm_sys = system(algorithm='sequential modular',
                    molar_tolerance=1e-6,
                    relative_molar_tolerance=1e-6,
                    method='fixed-point')
    time = bst.TicToc()
    
    time.tic()
    for i in range(1): dp_sys.simulate()
    t_phenomena = time.toc()
    
    time.tic()
    for i in range(1): sm_sys.simulate()
    t_sequential = time.toc()
    
    print('SM', t_sequential)
    print('PO', t_phenomena)
    
    for s_sm, s_dp in zip(sm_sys.streams, dp_sys.streams):
        actual = s_sm.mol
        value = s_dp.mol
        assert_allclose(actual, value, rtol=0.05, atol=1e-3)

def test_complex_acetic_acid_separation_system():
    import biosteam as bst
    import numpy as np
    bst.settings.set_thermo(['Water', 'AceticAcid', 'EthylAcetate'], cache=True)
    ideal_thermo = bst.settings.thermo.ideal()
    def create_system(alg):
        solvent_feed_ratio = 1
        chemicals = bst.settings.chemicals
        acetic_acid_broth = bst.Stream(
            ID='acetic_acid_broth', AceticAcid=6660, Water=43600, units='kg/hr',
            thermo=ideal_thermo
        )
        ethyl_acetate = bst.Stream(
            ID='fresh_solvent', EthylAcetate=15000, units='kg/hr',
            thermo=ideal_thermo
        )
        glacial_acetic_acid = bst.Stream('glacial_acetic_acid', thermo=ideal_thermo)
        wastewater = bst.Stream('wastewater', thermo=ideal_thermo)
        solvent_recycle = bst.Stream('solvent_rich', thermo=ideal_thermo)
        with bst.System(algorithm=alg) as sys:
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
                ins=(acetic_acid_broth, ethyl_acetate, solvent_recycle), 
                outs=('extract', 'raffinate'),
                top_chemical='EthylAcetate',
                feed_stages=(0, -1, -1),
                N_stages=5,
                collapsed_init=False,
                use_cache=True,
            )
            
            @extractor.add_specification(run=True)
            def adjust_fresh_solvent_flow_rate():
                broth = acetic_acid_broth.F_mass
                EtAc_recycle = solvent_recycle.imass['EthylAcetate']
                ethyl_acetate.imass['EthylAcetate'] = max(
                    0, broth * solvent_feed_ratio - EtAc_recycle
                )
            reflux = bst.Stream()
            ED = bst.MESHDistillation(
                'extract_distiller',
                ins=[extractor.extract, reflux],
                outs=['vapor', glacial_acetic_acid],
                LHK=('EthylAcetate', 'AceticAcid'),
                N_stages=10,
                feed_stages=(7, 0),
                reflux=None,
                boilup=4.37,
                thermo=ideal_thermo,
                use_cache=True,
            )
            condenser = bst.SinglePhaseStage(
                'condenser',
                ins=ED-0, 
                outs='condensate',
                phase='l',
                T=310,
                thermo=ideal_thermo,
            )
            water_rich = bst.Stream('water_rich')
            distillate = bst.Stream('distillate')
            settler = bst.StageEquilibrium(
                'settler',
                ins=(condenser-0, distillate), 
                outs=(solvent_recycle, water_rich, reflux),
                phases=('L', 'l'),
                top_chemical='EthylAcetate',
                top_split=0.273,
            )
            # settler.coupled_KL = True
            HX = bst.StageEquilibrium(
                'HX',
                ins=[extractor.raffinate, water_rich], 
                phases=('g', 'l'),
                B=0,
            )
            # RD = bst.ShortcutColumn(
            #     'raffinate_distiller',
            #     LHK=('EthylAcetate', 'Water'),
            #     ins=HX.outs[1],
            #     outs=[distillate, wastewater],
            #     partial_condenser=False,
            #     Lr=0.999,
            #     Hr=0.999,
            #     k=1.5,
            # )
            RD = bst.MESHDistillation(
                'raffinate_distiller',
                LHK=('EthylAcetate', 'Water'),
                ins=[HX-1],
                outs=['', wastewater, distillate],
                full_condenser=True,
                N_stages=8,
                feed_stages=(3, 3),
                reflux=1.5,
                boilup=2.,
                thermo=ideal_thermo,
            )
        return sys
    time = bst.TicToc()

    time.tic()
    sm = create_system('sequential modular')
    sm.flatten()
    sm.set_tolerance(
        rmol=1e-3, mol=1e-6, subsystems=True,
        method='fixed-point', maxiter=200,
    )
    sm.simulate()
    t_sequential = time.toc()
    # for i in sm.units: print(i, i.mass_balance_error())
    # return
    time.tic()
    po = create_system('phenomena oriented')
    po.set_tolerance(rmol=1e-3, mol=1e-6, 
                     subsystems=True,
                     method='fixed-point',
                     maxiter=200)
    po.simulate()
    t_phenomena = time.toc()
    
    print('SM', t_sequential)
    print('PO', t_phenomena)
    
    for s_sm, s_dp in zip(sm.streams, po.streams):
        actual = s_sm.mol
        value = s_dp.mol
        assert_allclose(actual, value, rtol=0.1)

# def test_complex_acetic_acid_separation_system():
#     import biosteam as bst
#     import numpy as np
#     @bst.SystemFactory
#     def system(ins, outs):
#         solvent_feed_ratio = 1.5
#         bst.settings.set_thermo(['Water', 'AceticAcid', 'EthylAcetate'], cache=True)
#         chemicals = bst.settings.chemicals
#         acetic_acid_broth = bst.Stream(
#             ID='acetic_acid_broth', AceticAcid=6660, Water=43600,
#         )
#         ethyl_acetate = bst.Stream(
#             ID='fresh_solvent', EthylAcetate=15000, units='kg/hr'
#         )
#         glacial_acetic_acid = bst.Stream(ID='glacial_acetic_acid')
#         wastewater = bst.Stream('wastewater')
#         solvent_recycle = bst.Stream('solvent_rich')
           
#         @ethyl_acetate.equation('material')
#         def fresh_solvent_flow_rate():
#             f = np.ones(chemicals.size)
#             r = np.zeros(chemicals.size)
#             v = r.copy()
#             index = chemicals.index('EthylAcetate')
#             r[index] = 1
#             v[index] = solvent_feed_ratio * acetic_acid_broth.F_mass / chemicals.EthylAcetate.MW
#             return (
#                 {ethyl_acetate: f,
#                   solvent_recycle: r},
#                   v
#             )
        
#         water_rich = bst.Stream('water_rich')
#         steam = bst.Stream('steam', Water=100, phase='g', T=390)
#         vapor_extract = bst.Stream('vapor_extract', phase='g')
#         liquid_extract = bst.Stream('liquid_extract', phase='l')
#         extractor = bst.MultiStageMixerSettlers(
#             'extractor', 
#             ins=(acetic_acid_broth, ethyl_acetate, solvent_recycle), 
#             outs=('extract', 'raffinate'),
#             feed_stages=(0, -1, -1),
#             N_stages=10,
#             use_cache=True,
#         )
        
#         @extractor.add_specification(run=True)
#         def adjust_fresh_solvent_flow_rate():
#             broth = acetic_acid_broth.F_mass
#             EtAc_recycle = solvent_recycle.imass['EthylAcetate']
#             ethyl_acetate.imass['EthylAcetate'] = max(
#                 0, broth * solvent_feed_ratio - EtAc_recycle
#             )
        
#         water_heater = bst.SinglePhaseStage('heater',
#             ins=[water_rich, extractor.raffinate],
#             outs=['carrier'],
#             phase='l',
#             T=360,
#         )
        
#         stripper = bst.Stripper('adiabatic_column',
#             N_stages=3, ins=[water_heater-0, steam], 
#             outs=['to_distillation', wastewater],
#             solute="AceticAcid", 
#             use_cache=True,
#         )
        
#         @steam.equation('material')
#         def steam_flow_rate():
#             feed, steam = stripper.ins
#             f = np.zeros(chemicals.size)
#             s = np.ones(chemicals.size)
#             v = f.copy()
#             index = chemicals.index('Water')
#             f[index] = -1
#             v[index] = 0
#             return (
#                 {feed: f,
#                   steam: s},
#                   v
#             )
        
#         @stripper.add_specification(run=True)
#         def adjust_steam_flow_rate():
#             feed, steam = stripper.ins
#             steam.imass['Water'] = feed.imass['Water']
        
#         reflux = bst.Stream('reflux')
        
#         distillation = bst.MESHDistillation(
#             'distillation',
#             N_stages=20,
#             ins=[vapor_extract, liquid_extract, reflux],
#             feed_stages=[8, 12, 0],
#             outs=['', glacial_acetic_acid, 'distillate'],
#             full_condenser=True,
#             reflux=1.0,
#             boilup=3.5,
#             use_cache=True,
#             LHK=('Water', 'AceticAcid'),
#             collapsed_init=False,
#         )
        
#         @distillation.add_specification
#         def skip_sim():
#             if not distillation.ins[0].isempty():
#                 distillation._run()
        
#         hx0 = bst.SinglePhaseStage(
#             'cooler',
#             ins=[distillation.outs[2], stripper.vapor],
#             outs=['cooled_distillate'],
#             T=315,
#             phase='l',
#         )
#         flash = bst.StageEquilibrium(
#             'flash',
#             ins=extractor.extract,
#             outs=[vapor_extract, liquid_extract],
#             B=5,
#             phases=('g', 'l')
#         )
#         settler = bst.StageEquilibrium(
#             'settler',
#             ins=hx0-0, 
#             outs=(solvent_recycle, water_rich, reflux),
#             phases=('L', 'l'),
#             top_chemical='EthylAcetate',
#             top_split=0.5,
#         )
#     sys = system()
#     sys.diagram()
#     sys.simulate()

if __name__ == '__main__':
    test_trivial_lle_case()
    test_trivial_vle_case()
    test_trivial_liquid_extraction_case()
    test_trivial_distillation_case()
    test_simple_acetic_acid_separation_no_recycle()
    test_simple_acetic_acid_separation_with_recycle()
    # test_complex_acetic_acid_separation_system()