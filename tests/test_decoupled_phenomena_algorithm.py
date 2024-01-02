# -*- coding: utf-8 -*-
"""
"""
import biosteam as bst
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
    
def test_acetic_acid_separation_no_recycle():
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

def test_acetic_acid_separation_system(solvent_feed_ratio=1.5):
    with bst.System(algorithm='decoupled phenomena') as sys:
        bst.settings.set_thermo(['Water', 'AceticAcid', 'EthylAcetate'], cache=True)
        chemicals = bst.settings.chemicals
        acetic_acid_broth = bst.Stream(
            ID='acetic_acid_broth', AceticAcid=500, Water=10000, units='kg/hr'
        )
        ethyl_acetate = bst.Stream(
            ID='acetic_acid_broth', EthylAcetate=15000, units='kg/hr'
        )
        glacial_acetic_acid = bst.Stream(ID='glacial_acetic_acid')
        wastewater = bst.Stream
        solvent_recycle = bst.Stream('solvent_rich')
        solvent = bst.Stream(ID='wastewater')
        # TODO:
        # @bst.mass_balance
        # def f():
        #     return [(
        #         {(fresh, 'EthylAcetate'): 1,
        #          (recycle, 'EthylAcetate'): 1},
        #          solvent_feed_ratio * acetic_acid_broth.F_mass
        #     )]
            
        solvent_mixer = bst.Mixer(ins=[ethyl_acetate, solvent_recycle], outs=solvent)
        solvent_mixer.solvent_feed_ratio = solvent_feed_ratio
        ideal_thermo = bst.settings.thermo.ideal()
        water_rich = bst.Stream('water_rich')
        steam = bst.Stream('steam', Water=100, phase='g', T=390)
        warm_extract = bst.Stream('warm_extract', thermo=ideal_thermo)
        hot_extract = bst.MultiStream('hot_extract', phases=('g', 'l'), thermo=ideal_thermo)
        extractor = bst.MultiStageMixerSettlers(
            'extractor', 
            ins=(acetic_acid_broth, solvent), 
            outs=('extract', 'raffinate'),
            N_stages=6,
            use_cache=True,
        )
        
        
        water_heat_integration = bst.HXprocess(
            ins=[extractor.raffinate, water_rich],
            outs=[wastewater, 'carrier']
        )
        
        absorber = bst.Absorber(
            N_stages=3, ins=[water_heat_integration-1, steam], 
            solute="AceticAcid", outs=['vapor', 'liquid'],
            use_cache=True,
        )
        @absorber.add_specification(run=False)
        def adjust_steam():
            feed, steam = absorber.ins
            if feed.isempty(): return
            steam.F_mass = feed.F_mass
            # speed_up = all([not i.isempty() for i in absorber.outs])
            absorber.run()
            # if speed_up:
            #     absorber.partition_data.update({
            #         'K': gmean([i.partition.K for i in absorber.stages]),
            #         'IDs': absorber.stages[0].partition.IDs,
            #         'phi': 0.5 # Initial phase fraction guess. This is optional.
            #     })
            #     absorber._setup()
        
        absorber.line = 'Absorber'
        # mixer = bst.Mixer(ins=[absorber.vapor, hot_extract], rigorous=True, thermo=ideal_thermo)
        # distillation = bst.ShortcutColumn(
        #     ins=[mixer-0],
        #     outs=['distillate', glacial_acetic_acid],
        #     thermo=ideal_thermo,
        #     k=1.2,
        #     LHK=('Water', 'AceticAcid'),
        #     y_top=0.99,
        #     x_bot=0.01,
        #     partial_condenser=False,
        # )
        distillation = bst.MESHDistillation(
            N_stages=10,
            ins=[hot_extract, absorber.vapor],
            feed_stages=[5, 0],
            outs=['', glacial_acetic_acid, 'distillate'],
            full_condenser=True,
            reflux=1.0,
            boilup=3.5,
            use_cache=True,
            LHK=('Water', 'AceticAcid'),
        )
        distillation.collapsed_init = False
        
        hx0 = bst.HXprocess(
            ins=[distillation.outs[2], extractor.extract],
            outs=['cooled_distillate', warm_extract],
            thermo=ideal_thermo,
        )
        hx1 = bst.HXutility(
            ins=hx0-1,
            outs=hot_extract,
            V=0.95,
            rigorous=True,
            heat_only=True,
            thermo=ideal_thermo,
        )
        settler = bst.LLESettler(
            ins=hx0-0, 
            outs=(solvent_recycle, water_rich)


if __name__ == '__main__':
    # test_trivial_lle_case()
    # test_trivial_vle_case()
    # test_trivial_liquid_extraction_case()
    # test_trivial_distillation_case()
    test_acetic_acid_separation_no_recycle()