# -*- coding: utf-8 -*-
"""
"""
import pytest
import biosteam as bst
from numpy.testing import assert_allclose

def test_multi_stage_vle():
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
    sys.set_tolerance(mol=0.1, rmol=0.001, maxiter=100)
    sys.simulate()
    assert_allclose(
        MSE.stages[0].vapor.mol, 
        stage0.vapor.mol,
        atol=0.1,
        rtol=0.01,
    )
    assert_allclose(
        MSE.stages[0].liquid.mol, 
        stage0.liquid.mol,
        atol=0.1,
        rtol=0.01,
    )
    assert_allclose(
        MSE.stages[1].vapor.mol, 
        stage1.vapor.mol,
        atol=0.1,
        rtol=0.01,
    )
    assert_allclose(
        MSE.stages[1].liquid.mol, 
        stage1.liquid.mol,
        atol=0.1,
        rtol=0.01,
    )
    
if __name__ == '__main__':
    test_multi_stage_vle()



