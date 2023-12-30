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
        feed = bst.Stream(AceticAcid=6660, Water=23600)
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
            method='anderson',
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

    
if __name__ == '__main__':
    # test_trivial_lle_case()
    # test_trivial_vle_case()
    # test_trivial_liquid_extraction_case()
    # test_trivial_distillation_case()
    test_acetic_acid_separation_no_recycle()