# -*- coding: utf-8 -*-
"""
"""
import pytest
import biosteam as bst
import thermosteam as tmo
from math import exp, log
import numpy as np
from thermosteam.constants import R
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
   
def test_lactic_acid_ethanol_reactive_distillation():
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
    
    rxn = Esterification('LacticAcid + Ethanol -> H2O + EthylLactate', reactant='LacticAcid')
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
        [6.585111268569887e-06, 3.6465771244802687e-09,
         3.110437134095077, 8.015523262914513],
        [0.019007540635984595, 1.9809858706061696,
         6.908576991652181, 1.9654626113382356],
    ]
    for i, j in zip(distillation.outs, flows):    
        assert_allclose(i.mol, j, rtol=1e-5, atol=1e-3)
    sys = bst.System.from_units(units=[distillation])
    sys._setup()
    sys.run_phenomena()
    for i, j in zip(distillation.outs, flows):    
        assert_allclose(i.mol, j, rtol=1e-5, atol=1e-3)
    
def test_acetic_acid_reactive_distillation():
    import biosteam as bst
    import thermosteam as tmo
    from thermosteam.constants import R
    import numpy as np
    from math import exp
    from numpy.testing import assert_allclose
    bst.settings.set_thermo(['Water', 'AceticAcid', 'MethylAcetate', 'Methanol'], cache=True)
    Gamma = tmo.equilibrium.DortmundActivityCoefficients(bst.settings.chemicals)
    class Esterification(bst.KineticReaction):
        def volume(self, stream): # kg of catalyst
            rho_cat = 770 # kg / m3
            liquid_volume = self.liquid_volume
            catalyst_volume = 0.5 * liquid_volume
            catalyst_mass = catalyst_volume * rho_cat
            return catalyst_mass
        
        def rate(self, stream): # kmol/kg-catalyst/hr
            F_mol = stream.F_mol
            if not F_mol: return 0
            K_AcOH = 3.15
            K_MeOH = 5.64
            K_MeOAc = 4.15
            K_H2O = 5.24
            K_adsorption = np.array([K_H2O, K_AcOH, K_MeOAc, K_MeOH])
            MWs = stream.chemicals.MW
            k0f = 8.497e6
            EAf = 60.47
            k0r = 6.127e5
            EAr = 63.73
            T = stream.T
            mol = stream.mol
            x = mol / F_mol
            gamma = Gamma(x, T)
            activity = gamma * x
            a = K_adsorption * activity / MWs
            H2O, LA, BuLA, BuOH = activity
            kf = k0f * exp(-EAf / (R * T))
            kr = k0r * exp(-EAr / (R * T))
            a_sum = sum(a)
            return 3600 * (kf * LA * BuOH - kr * BuLA * H2O) / (a_sum * a_sum) # kmol / kg-catalyst / hr
        
    methanol = bst.Stream(
        Methanol=30,
    )
    methanol.vle(V=1, P=101325)
    acetic_acid = bst.Stream(
        AceticAcid=10, Water=200
    )
    acetic_acid.vle(V=0, P=101325)
    stage_reactions = {
        i: Esterification('AceticAcid + Methanol -> MethylAcetate + Water', reactant='AceticAcid')
        for i in range(3, 9)
    }
    distillation = bst.MESHDistillation(
        N_stages=10,
        ins=[methanol, acetic_acid],
        feed_stages=[3, 8],
        outs=['distillate', 'bottoms_product'],
        full_condenser=False,
        reflux=1.0,
        boilup=2.0,
        use_cache=True,
        LHK=('MethylAcetate', 'Water'),
        stage_reactions=stage_reactions,
        method='fixed-point',
        maxiter=100,
    )
    distillation.simulate()
    flows = [
        [6.585111268569887e-06, 3.6465771244802687e-09,
         3.110437134095077, 8.015523262914513],
        [0.019007540635984595, 1.9809858706061696,
         6.908576991652181, 1.9654626113382356],
    ]
    for i, j in zip(distillation.outs, flows):    
        assert_allclose(i.mol, j, rtol=1e-5, atol=1e-3)
    sys = bst.System.from_units(units=[distillation])
    sys._setup()
    sys.run_phenomena()
    for i, j in zip(distillation.outs, flows):    
        assert_allclose(i.mol, j, rtol=1e-5, atol=1e-3)
    
    
def test_lactic_acid_butanol_reactive_distillation():
    import biosteam as bst
    from thermosteam.constants import R
    from math import exp
    from numpy.testing import assert_allclose
    bst.settings.set_thermo(['Water', 'LacticAcid', 'ButylLactate', 'Butanol'], cache=True)
    
    class Esterification(bst.KineticReaction):
        
        def volume(self, stream): # kg of catalyst
            rho_cat = 770 # kg / m3
            liquid_volume = self.liquid_volume
            catalyst_volume = 0.5 * liquid_volume
            catalyst_mass = catalyst_volume * rho_cat
            return catalyst_mass
        
        def rate(self, stream): # kmol/kg-catalyst/hr
            T = stream.T
            # if T > 370: return 0 # Prevents multiple steady states
            kf = 2.59e4 * exp(-5.340e4 / (R * T))
            kr = 3.80e3 * exp(-5.224e4 / (R * T))
            H2O, LA, BuLA, BuOH = stream.mol / stream.F_mol
            return 3600 * (kf * LA * BuOH - kr * BuLA * H2O) # kmol / kg-catalyst / hr
        
    stream = bst.Stream(
        H2O=10, Butanol=10, LacticAcid=2, T=355,
    )
    stream.vle(V=0, P=101325)
    stage_reactions = {
        i: Esterification('LacticAcid + Butanol -> Water + ButylLactate', reactant='LacticAcid')
        for i in range(2)
    }
    distillation = bst.MESHDistillation(
        N_stages=5,
        ins=[stream],
        feed_stages=[2],
        outs=['distillate', 'bottoms_product'],
        full_condenser=False,
        reflux=1.0,
        boilup=2.0,
        use_cache=True,
        LHK=('Butanol', 'ButylLactate'),
        stage_reactions=stage_reactions,
        method='fixed-point',
        maxiter=100,
    )
    distillation.simulate()
    flows = [
        [6.585111268569887e-06, 3.6465771244802687e-09,
         3.110437134095077, 8.015523262914513],
        [0.019007540635984595, 1.9809858706061696,
         6.908576991652181, 1.9654626113382356],
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
        [0.0, 0.0, 0.0]
        [0.11220194378011927, 4.351310579539389, 22.35906337238964]
        [22.06172189576452, 0.14410640814689596, 87.14686571087641]
    ]
    for i, j in zip(distillation.outs, flows):    
        assert_allclose(i.mol, j, rtol=1e-5, atol=1e-3)
    
    
if __name__ == '__main__':
    test_multi_stage_adiabatic_vle()
    test_distillation()
    test_lactic_acid_ethanol_reactive_distillation()



