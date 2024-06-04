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
   
def test_reactive_distillation():
    import biosteam as bst
    from thermosteam.constants import R
    from math import exp
    from numpy.testing import assert_allclose
    bst.settings.set_thermo(['EthylLactate', 'LacticAcid', 'H2O', 'Ethanol'], cache=True)
    
    class Esterification(bst.KineticReaction):
        
        def volume(self, stream):
            return 0.01 # Kg of catalyst
        
        def rate(self, stream):
            T = stream.T
            if T > 365: return 0 # Prevents multiple steady states.
            kf = 6.52e3 * exp(-4.8e4 / (R * T))
            kr = 2.72e3 * exp(-4.8e4 / (R * T))
            LaEt, La, H2O, EtOH = stream.mol / stream.F_mol
            return 3600 * (kf * La * EtOH - kr * LaEt * H2O) # kmol / kg-catalyst / hr
    
    rxn = Esterification('LacticAcid + Ethanol -> H2O + EthylLactate')
    stream = bst.Stream(
        H2O=10, Ethanol=10, LacticAcid=2, T=355,
    )
    stream.vle(V=0, P=101325, liquid_conversion=rxn)
    distillation = bst.MESHDistillation(
        N_stages=5,
        ins=[stream],
        feed_stages=[2],
        outs=['distillate', 'bottoms_product'],
        full_condenser=False,
        reflux=1.0,
        boilup=2.0,
        use_cache=True,
        LHK=('H2O', 'EthylLactate'),
        stage_reactions={i: rxn for i in range(2)},
        method='fixed-point',
        maxiter=100,
    )
    distillation.simulate()
    flows = [
        [6.575895653800856e-06, 3.6651391022703107e-09, 
         3.1210737204600005, 8.043951319613521],
        [0.019035234263940305, 2.0190418064944553, 
         6.897968089699594, 1.9750904905460722],
    ]
    for i, j in zip(distillation.outs, flows):    
        assert_allclose(i.mol, j, rtol=1e-5, atol=1e-3)
    sys = bst.System.from_units(units=[distillation])
    sys._setup()
    sys.run_phenomena()
    for i, j in zip(distillation.outs, flows):    
        assert_allclose(i.mol, j, rtol=1e-5, atol=1e-3)
    
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
        outs=['', 'bottoms_product', 'distillate'],
        full_condenser=True,
        reflux=1.0,
        boilup=3.5,
        use_cache=True,
        LHK=('Water', 'AceticAcid'),
        method='fixed-point',
    )
    distillation.simulate()
    flows = [
        [0.0, 0.0, 0.0],
        [0.1120748575014966, 4.350505208368632, 22.358805264239038],
        [22.055925142498495, 0.14389479163136612, 87.1521947357609],
    ]
    for i, j in zip(distillation.outs, flows):    
        assert_allclose(i.mol, j, rtol=1e-5, atol=1e-3)
    
    
if __name__ == '__main__':
    test_multi_stage_adiabatic_vle()
    test_distillation()
    test_reactive_distillation()



