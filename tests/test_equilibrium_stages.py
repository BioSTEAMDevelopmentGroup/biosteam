# -*- coding: utf-8 -*-
"""
"""
import pytest
import biosteam as bst
from numpy.testing import assert_allclose

def test_multi_stage_adiabatic_vle():
    bst.settings.set_thermo(['AceticAcid', 'EthylAcetate', 'Water', 'MTBE'], cache=True)
    feed = bst.Stream('feed', Water=75, AceticAcid=5, MTBE=20, T=320)
    steam = bst.Stream('steam', Water=100, phase='g', T=390)
    MSE = bst.MultiStageEquilibrium(N_stages=2, ins=[feed, steam], feed_stages=[0, -1],
        outs=['vapor', 'liquid'],
        phases=('g', 'l'),
    )
    MSE.simulate()
    with bst.System() as sys:
        vapor_recycle = bst.Stream()
        M0 = bst.Mixer(ins=[feed.copy(), vapor_recycle], conserve_phases=True)
        stage0 = bst.Flash(ins=M0-0, P=101325, Q=0)
        M1 = bst.Mixer(ins=[stage0.liquid, steam.copy()], conserve_phases=True)
        stage1 = bst.Flash(ins=M1-0,
                           outs=[vapor_recycle, 'liquid'], 
                           P=101325, Q=0)
    sys.set_tolerance(mol=0.01, rmol=0.001)
    sys.simulate()
    assert_allclose(
        MSE.stages[0].vapor.mol, 
        stage0.vapor.mol,
        atol=1,
        rtol=0.01,
    )
    assert_allclose(
        MSE.stages[0].liquid.mol, 
        stage0.liquid.mol,
        atol=1,
        rtol=0.01,
    )
    assert_allclose(
        MSE.stages[1].vapor.mol, 
        stage1.vapor.mol,
        atol=1,
        rtol=0.01,
    )
    assert_allclose(
        MSE.stages[1].liquid.mol, 
        stage1.liquid.mol,
        atol=1,
        rtol=0.01,
    )
    assert_allclose(
        MSE.stages[1].vapor.T, 
        stage1.vapor.T,
        atol=1,
        rtol=0.01,
    )
    assert_allclose(
        MSE.stages[1].liquid.T, 
        stage1.liquid.T,
        atol=1,
        rtol=0.01,
    )
   
def test_distillation():
    import biosteam as bst
    bst.settings.set_thermo(
        ['Water', 'AceticAcid', 'EthylAcetate'],
        cache=True
    )
    hot_extract = bst.MultiStream(
        phases=('g', 'l'), T=358.05, P=101325,
        g=[('Water', 20.29), 
           ('AceticAcid', 3.872), 
           ('EthylAcetate', 105.2)],
        l=[('Water', 1.878), 
           ('AceticAcid', 0.6224), 
           ('EthylAcetate', 4.311)]
    )
    distillation = bst.MESHDistillation(
        N_stages=10,
        ins=[hot_extract],
        feed_stages=[5],
        outs=['', 'glacial_acetic_acid', 'distillate'],
        full_condenser=True,
        reflux=1.0,
        boilup=3.5,
        use_cache=True,
        LHK=('Water', 'AceticAcid'),
        method='anderson',
    )
    distillation.simulate()
    flows = (
        [0.0, 0.0, 0.0],
        [1.3946134158640209, 3.72342260007752, 32.960792036663534],
        [20.773386584135977, 0.7709773999224796, 76.55020796333646]
    )
    for i, j in zip(distillation.outs, flows):    
        assert_allclose(i.mol, j)
    
if __name__ == '__main__':
    test_multi_stage_adiabatic_vle()
    test_distillation()



