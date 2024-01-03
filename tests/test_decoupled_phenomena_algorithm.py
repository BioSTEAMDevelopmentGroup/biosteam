# -*- coding: utf-8 -*-
"""
"""
import biosteam as bst
import numpy as np
from numpy.testing import assert_allclose

def test_trivial_lle_case():
    # lle
    with bst.System(algorithm='decoupled phenomena') as sys:
        bst.settings.set_thermo(['Water', 'Octanol', 'Methanol'], cache=True)
        feed = bst.Stream(Water=50, Methanol= 5, units='kg/hr')
        solvent = bst.Stream(Octanol=50, units='kg/hr')
        stage = bst.StageEquilibrium(phases=('L', 'l'), ins=[feed, solvent])
    sys.simulate()
    streams = stage.ins + stage.outs
    actuals = [i.mol.copy() for i in streams]
    sys.run_decoupled_phenomena()
    values = [i.mol for i in streams]
    for actual, value in zip(actuals, values):
        assert_allclose(actual, value)

def test_trivial_vle_case():
    # vle
    with bst.System(algorithm='decoupled phenomena') as sys:
        bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
        liquid = bst.Stream(Water=50, Ethanol=10, units='kg/hr')
        vapor = bst.Stream(Ethanol=10, Water=50, phase='g', units='kg/hr')
        stage = bst.StageEquilibrium(phases=('g', 'l'), ins=[liquid, vapor])
    sys.simulate()
    streams = stage.ins + stage.outs
    actuals = [i.mol.copy() for i in streams]
    sys.run_decoupled_phenomena()
    values = [i.mol for i in streams]
    for actual, value in zip(actuals, values):
        assert_allclose(actual, value)
 
def test_trivial_liquid_extraction_case():
    # liquid extraction
    with bst.System(algorithm='decoupled phenomena') as sys:
        bst.settings.set_thermo(['Water', 'Methanol', 'Octanol'], cache=True)
        feed = bst.Stream(Water=500, Methanol=50)
        solvent = bst.Stream(Octanol=500, T=330)
        MSE = bst.MultiStageEquilibrium(N_stages=2, ins=[feed, solvent], phases=('L', 'l'))
    sys.simulate()
    extract, raffinate = MSE.outs
    actual = extract.imol['Methanol'] / feed.imol['Methanol']
    T_actual = extract.T
    sys.run_decoupled_phenomena()
    value = extract.imol['Methanol'] / feed.imol['Methanol']
    T = extract.T
    assert_allclose(actual, value, rtol=0.001)
    assert_allclose(T_actual, T, rtol=0.001)

def test_trivial_distillation_case():    
    # distillation
    with bst.System(algorithm='decoupled phenomena') as sys:
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
    sys.run_decoupled_phenomena()
    assert round(vapor.imol['Ethanol'] / feed.imol['Ethanol'], 2) == actual
    
def test_simple_acetic_acid_separation_no_recycle():
    with bst.System(algorithm='decoupled phenomena') as sys:
        bst.settings.set_thermo(['Water', 'AceticAcid', 'EthylAcetate'], cache=True)
        feed = bst.Stream(AceticAcid=6660, Water=43600)
        solvent = bst.Stream(EthylAcetate=65000)
        LE = bst.MultiStageEquilibrium(
            N_stages=6, ins=[feed, solvent], phases=('L', 'l'),
            maxiter=200,
            use_cache=True,
            method='fixed-point',
        )
        DAA = bst.MultiStageEquilibrium(N_stages=6, ins=[LE-0], feed_stages=[3],
            outs=['vapor', 'liquid'],
            stage_specifications={0: ('Reflux', 0.673), -1: ('Boilup', 2.57)},
            phases=('g', 'l'),
            method='anderson',
            use_cache=True,
        )
        DEA = bst.MultiStageEquilibrium(N_stages=6, ins=[LE-1], feed_stages=[3],
            outs=['vapor', 'liquid'],
            stage_specifications={0: ('Reflux', 0.673), -1: ('Boilup', 2.57)},
            phases=('g', 'l'),
            method='excitingmixing',
            use_cache=True,
        )
    sys.simulate()
    for i in sys.units: print(i.ID, i.iter)
    streams = [*sys.ins, *sys.outs]
    actuals = [i.mol.copy() for i in streams]
    sys.run_decoupled_phenomena()
    values = [i.mol for i in streams]
    for actual, value in zip(actuals, values):
        assert_allclose(actual, value, rtol=0.001)

def test_simple_acetic_acid_separation_with_recycle():
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
            outs=['', 'liquid', recycle],
            stage_specifications={0: ('Reflux', float('inf')), -1: ('Boilup', 2.57)},
            bottom_side_draws={0: 0.673 / (1 + 0.673)},
            phases=('g', 'l'),
            maxiter=200,
            method='fixed-point',
            use_cache=True,
        )
        chemicals = bst.settings.chemicals
        
        @LE.add_specification(run=True)
        def fresh_solvent_flow_rate():
            broth = feed.F_mass
            EtAc_recycle = recycle.imass['EthylAcetate']
            solvent.imass['EthylAcetate'] = max(
                0, broth * solvent_feed_ratio - EtAc_recycle
            )
        
        @solvent.equation('mol')
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
    dp_sys = system(algorithm='decoupled phenomena', 
                    molar_tolerance=1e-6,
                    relative_molar_tolerance=1e-6)
    sm_sys = system(algorithm='sequential modular',
                    molar_tolerance=1e-6,
                    relative_molar_tolerance=1e-6)
    time = bst.TicToc()
    
    time.tic()
    for i in range(1): dp_sys.simulate()
    time.toc()
    
    time.tic()
    for i in range(1): sm_sys.simulate()
    time.toc()
    print(time.record)
    
    dp_sys.show()
    sm_sys.show()
    
    values = [i.mol for i in dp_sys.streams]
    actuals = [i.mol for i in sm_sys.streams]
    
    for actual, value in zip(actuals, values):
        assert_allclose(actual, value, rtol=0.001)

def test_complex_acetic_acid_separation_system(solvent_feed_ratio=1):
    with bst.System(algorithm='decoupled phenomena') as sys:
        bst.settings.set_thermo(['Water', 'AceticAcid', 'EthylAcetate'], cache=True)
        chemicals = bst.settings.chemicals
        acetic_acid_broth = bst.Stream(
            ID='acetic_acid_broth', AceticAcid=500, Water=10000, units='kg/hr'
        )
        ethyl_acetate = bst.Stream(
            ID='fresh_solvent', EthylAcetate=15000, units='kg/hr'
        )
        glacial_acetic_acid = bst.Stream(ID='glacial_acetic_acid')
        wastewater = bst.Stream('wastewater')
        solvent_recycle = bst.Stream('solvent_rich')

        @ethyl_acetate.equation('mol')
        def fresh_solvent_flow_rate():
            f = np.ones(chemicals.size)
            r = np.zeros(chemicals.size)
            v = r.copy()
            index = chemicals.index('EthylAcetate')
            r[index] = 1
            v[index] = solvent_feed_ratio * acetic_acid_broth.F_mass
            return (
                {ethyl_acetate: f,
                 solvent_recycle: r},
                 v
            )
            
        ideal_thermo = bst.settings.thermo.ideal()
        water_rich = bst.Stream('water_rich')
        steam = bst.Stream('steam', Water=100, phase='g', T=390)
        vapor_extract = bst.Stream('vapor_extract', phase='g', thermo=ideal_thermo)
        liquid_extract = bst.Stream('liquid_extract', phase='l', thermo=ideal_thermo)
        extractor = bst.MultiStageMixerSettlers(
            'extractor', 
            ins=(acetic_acid_broth, ethyl_acetate, solvent_recycle), 
            outs=('extract', 'raffinate'),
            feed_stages=(0, -1, -1),
            N_stages=3,
            use_cache=True,
        )
        
        water_heater = bst.SinglePhaseStage(
            ins=[water_rich, extractor.raffinate],
            outs=['carrier'],
            phase='l',
            T=360,
        )
        
        absorber = bst.Absorber(
            N_stages=3, ins=[water_heater-0, steam], 
            outs=['to_distillation', wastewater],
            solute="AceticAcid", 
            use_cache=True,
        )
        
        @steam.equation
        def steam_flow_rate():
            feed, steam = absorber.ins
            f = np.zeros(chemicals.size)
            s = np.ones(chemicals.size)
            v = s.copy()
            index = chemicals.index('Water')
            f[index] = -1
            v[index] = 0
            return (
                {feed: f,
                 steam: s},
                 v
            )
        
        absorber.line = 'Absorber'
        
        distillation = bst.MESHDistillation(
            N_stages=10,
            ins=[vapor_extract, liquid_extract, absorber.vapor],
            feed_stages=[4, 6, 0],
            outs=['', glacial_acetic_acid, 'distillate'],
            full_condenser=True,
            reflux=1.0,
            boilup=3.5,
            use_cache=True,
            LHK=('Water', 'AceticAcid'),
            collapsed_init=False,
        )
        hx0 = bst.SinglePhaseStage(
            ins=[distillation.outs[2]],
            outs=['cooled_distillate'],
            T=320,
            phase='l',
            thermo=ideal_thermo,
        )
        flash = bst.StageEquilibrium(
            ins=extractor.extract,
            outs=[vapor_extract, liquid_extract],
            B=5,
            thermo=ideal_thermo,
            phases=('g', 'l')
        )
        settler = bst.StageEquilibrium(
            ins=hx0-0, 
            outs=(solvent_recycle, water_rich),
            phases=('L', 'l'),
            solvent='EthylAcetate',
        )
    sys.diagram()
    sys.simulate()

if __name__ == '__main__':
    # test_trivial_lle_case()
    # test_trivial_vle_case()
    # test_trivial_liquid_extraction_case()
    # test_trivial_distillation_case()
    # test_simple_acetic_acid_separation_no_recycle()
    test_simple_acetic_acid_separation_with_recycle()
    # test_complex_acetic_acid_separation_system()