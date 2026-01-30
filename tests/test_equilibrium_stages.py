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
    
def test_multi_stage_adiabatic_vle_non_condensables():
    import biosteam as bst
    bst.settings.set_thermo(['Water', 'Ethanol', 'CO2', 'N2'], cache=True)
    liq = bst.Stream('feed', Water=100, units='kg/hr', T=320)
    gas = bst.Stream('steam', CO2=100, N2=100, Ethanol=40, units='kg/hr', phase='g', T=305.15)
    MSE = bst.MultiStageEquilibrium(N_stages=10, ins=[liq, gas], feed_stages=[0, -1],
        outs=['vapor', 'liquid'],
        phases=('g', 'l'),
    )
    MSE.simulate()

def test_small_esterification_column():
    import biosteam as bst
    bst.settings.set_thermo(['Water', 'LacticAcid', 'MethylLactate', 'Methanol', 'SuccinicAcid'], cache=True)
    from math import exp, log
    from thermosteam.constants import R
    
    class Esterification(bst.KineticReaction):
        catalyst_fraction = 0.5 
        
        def volume(self, stream): # kg of catalyst
            rho_cat = 770 # kg / m3
            liquid_volume = self.liquid_volume
            catalyst_volume = self.catalyst_fraction * liquid_volume
            catalyst_mass = catalyst_volume * rho_cat
            return catalyst_mass
        
        def rate(self, stream):
            T = stream.T
            kf = 2.98e4 * exp(-51300 / (R * T))
            kr = 2.74e2 * exp(-47200 / (R * T))
            H2O, La, LaMe, MeOH, _ = stream.mol / stream.F_mol
            return 3600 * (kf * La * MeOH - kr * LaMe * H2O) # kmol / kg-catalyst / hr
    
    
    feed = bst.Stream(
        'feed',
        LacticAcid=4.174,
        Water=5.460,
        SuccinicAcid=0.531,
        MethylLactate=1e-6,
        total_flow=10.165,
        P=0.106 * 101325,
        T=72.5 + 273.15
    )
    feed.T = feed.bubble_point_at_P(P=0.2 * 101325).T
    methanol = bst.Stream('methanol', Methanol=16.695, P=0.2 * 101325, phase='g')
    methanol.T = methanol.dew_point_at_P().T
    esterification = bst.MESHDistillation(
        'esterification',
        ins=(feed, methanol), 
        outs=('esterification_distillate', 'bottoms'),
        N_stages=3,
        feed_stages=(1, 1),
        reflux=1, # Unknown (update later)
        boilup=1,
        # bottoms_product_to_feed=0.0334, # Fraction of feed.
        full_condenser=True,
        stage_reactions={
            i: Esterification('LacticAcid + Methanol -> Water + MethylLactate', reactant='LacticAcid')
            for i in range(1, 2)
        },
        maxiter=50,
        LHK=('Methanol', 'MethylLactate'),
        P=0.2 * 101325,
        use_cache=True,
        algorithms=('inside out',)
    )
    esterification.simulate()

def test_esterification_column():
    import biosteam as bst
    bst.settings.set_thermo(['Water', 'LacticAcid', 'MethylLactate', 'Methanol', 'SuccinicAcid'], cache=True)
    from math import exp, log
    from thermosteam.constants import R
    
    class Esterification(bst.KineticReaction):
        catalyst_fraction = 0.5 
        
        def volume(self, stream): # kg of catalyst
            rho_cat = 770 # kg / m3
            liquid_volume = self.liquid_volume
            catalyst_volume = self.catalyst_fraction * liquid_volume
            catalyst_mass = catalyst_volume * rho_cat
            return catalyst_mass
        
        def rate(self, stream):
            T = stream.T
            kf = 2.98e4 * exp(-51300 / (R * T))
            kr = 2.74e2 * exp(-47200 / (R * T))
            H2O, La, LaMe, MeOH, _ = stream.mol / stream.F_mol
            rate = 3600 * (kf * La * MeOH - kr * LaMe * H2O) # kmol / kg-catalyst / hr
            print(rate)
            return 0
    
    
    feed = bst.Stream(
        'feed',
        LacticAcid=4.174,
        Water=5.460,
        SuccinicAcid=0.531,
        MethylLactate=1e-6,
        total_flow=10.165,
        P=0.106 * 101325,
        T=72.5 + 273.15
    )
    feed.T = feed.bubble_point_at_P(P=0.2 * 101325).T
    methanol = bst.Stream('methanol', Methanol=16.695, P=0.2 * 101325, phase='g')
    methanol.T = methanol.dew_point_at_P().T
    esterification = bst.MESHDistillation(
        'esterification',
        ins=(feed, methanol), 
        outs=('esterification_distillate', 'bottoms'),
        N_stages=21,
        feed_stages=(1, 19),
        reflux=1, # Unknown (update later)
        boilup=1,
        # bottoms_product_to_feed=0.0334, # Fraction of feed.
        full_condenser=True,
        stage_reactions={
            i: Esterification('LacticAcid + Methanol -> Water + MethylLactate', reactant='LacticAcid')
            for i in range(1, 19)
        },
        maxiter=50,
        LHK=('Methanol', 'MethylLactate'),
        P=0.2 * 101325,
        use_cache=True,
        algorithms=('inside out',)
    )
    esterification.simulate()
    
def test_esterification_column_no_reaction():
    import biosteam as bst
    bst.settings.set_thermo(
        ['Water', 'LacticAcid', 'MethylLactate', 'Methanol', 'SuccinicAcid'],
        pkg='ideal gas',
        cache=True
    )
    
    feed = bst.Stream(
        'feed',
        LacticAcid=4.174,
        Water=5.460,
        SuccinicAcid=0.531,
        MethylLactate=1e-6,
        total_flow=10.165,
        P=0.106 * 101325,
        T=72.5 + 273.15
    )
    feed.T = feed.bubble_point_at_P(P=0.2 * 101325).T
    methanol = bst.Stream('methanol', Methanol=16.695, P=0.2 * 101325, phase='g')
    methanol.T = methanol.dew_point_at_P().T
    esterification = bst.MESHDistillation(
        'esterification',
        ins=(feed, methanol), 
        outs=('esterification_distillate', 'bottoms'),
        N_stages=21,
        feed_stages=(1, 19),
        reflux=1, # Unknown (update later)
        boilup=1,
        # bottoms_product_to_feed=0.0334, # Fraction of feed.
        full_condenser=True,
        maxiter=50,
        LHK=('Methanol', 'MethylLactate'),
        P=0.2 * 101325,
        use_cache=True,
        algorithms=('inside out',)
    )
    esterification.simulate()
    
   
# def test_lactic_acid_ethanol_reactive_distillation():
#     import biosteam as bst
#     from thermosteam.constants import R
#     from math import exp
#     from numpy.testing import assert_allclose
#     bst.settings.set_thermo(['EthylLactate', 'LacticAcid', 'H2O', 'Ethanol'], cache=True)
    
#     class Esterification(bst.KineticReaction):
        
#         def volume(self, stream):
#             return 0.01 # Kg of catalyst
        
#         def rate(self, stream):
#             T = stream.T
#             if T > 365: return 0 # Prevents multiple steady states.
#             kf = 6.52e3 * exp(-4.8e4 / (R * T))
#             kr = 2.72e3 * exp(-4.8e4 / (R * T))
#             LaEt, La, H2O, EtOH = stream.mol / stream.F_mol
#             return 3600 * (kf * La * EtOH - kr * LaEt * H2O) # kmol / kg-catalyst / hr
    
#     rxn = Esterification('LacticAcid + Ethanol -> H2O + EthylLactate', reactant='LacticAcid')
#     stream = bst.Stream(
#         H2O=10, Ethanol=10, LacticAcid=2, T=355,
#     )
#     stream.vle(V=0, P=101325, liquid_conversion=rxn)
#     distillation = bst.MESHDistillation(
#         N_stages=5,
#         ins=[stream],
#         feed_stages=[2],
#         outs=['distillate', 'bottoms_product'],
#         full_condenser=False,
#         reflux=1.0,
#         boilup=2.0,
#         use_cache=True,
#         LHK=('H2O', 'EthylLactate'),
#         stage_reactions={i: rxn for i in range(1,3)},
#         method='fixed-point',
#         maxiter=100,
#     )
#     distillation.simulate()
#     flows = [
#         [6.779963298101173e-06, 3.6172767647177665e-09, 3.1106893353339236, 8.015420206297735],
#         [0.01955718389933403, 1.9804360325200911, 6.908874628528711, 1.9650158298396327],
#     ]
#     for i, j in zip(distillation.outs, flows):    
#         assert_allclose(i.mol, j, rtol=1e-5, atol=1e-3)
#     sys = bst.System.from_units(units=[distillation])
#     sys._setup()
#     sys.run_phenomena()
#     for i, j in zip(distillation.outs, flows):    
#         assert_allclose(i.mol, j, rtol=1e-5, atol=1e-3)
    
# def test_acetic_acid_reactive_distillation():
#     import biosteam as bst
#     import thermosteam as tmo
#     from thermosteam.constants import R
#     import numpy as np
#     from math import exp
#     from numpy.testing import assert_allclose
#     bst.settings.set_thermo(['Water', 'AceticAcid', 'MethylAcetate', 'Methanol'], cache=True)
#     Gamma = tmo.equilibrium.DortmundActivityCoefficients(bst.settings.chemicals)
#     class Esterification(bst.KineticReaction):
#         def volume(self, stream): # kg of catalyst
#             rho_cat = 770 # kg / m3
#             liquid_volume = self.liquid_volume
#             catalyst_volume = 0.5 * liquid_volume
#             catalyst_mass = catalyst_volume * rho_cat
#             return catalyst_mass
        
#         def rate(self, stream): # kmol/kg-catalyst/hr
#             F_mol = stream.F_mol
#             if not F_mol: return 0
#             K_AcOH = 3.15
#             K_MeOH = 5.64
#             K_MeOAc = 4.15
#             K_H2O = 5.24
#             K_adsorption = np.array([K_H2O, K_AcOH, K_MeOAc, K_MeOH])
#             MWs = stream.chemicals.MW
#             k0f = 8.497e6
#             EAf = 60.47
#             k0r = 6.127e5
#             EAr = 63.73
#             T = stream.T
#             mol = stream.mol
#             x = mol / F_mol
#             gamma = Gamma(x, T)
#             activity = gamma * x
#             a = K_adsorption * activity / MWs
#             H2O, LA, BuLA, BuOH = activity
#             kf = k0f * exp(-EAf / (R * T))
#             kr = k0r * exp(-EAr / (R * T))
#             a_sum = sum(a)
#             return 3600 * (kf * LA * BuOH - kr * BuLA * H2O) / (a_sum * a_sum) # kmol / kg-catalyst / hr
        
#     methanol = bst.Stream(
#         Methanol=30,
#     )
#     methanol.vle(V=1, P=101325)
#     acetic_acid = bst.Stream(
#         AceticAcid=10, Water=200
#     )
#     acetic_acid.vle(V=0, P=101325)
#     stage_reactions = {
#         i: Esterification('AceticAcid + Methanol -> MethylAcetate + Water', reactant='AceticAcid')
#         for i in range(3, 9)
#     }
#     distillation = bst.MESHDistillation(
#         N_stages=10,
#         ins=[methanol, acetic_acid],
#         feed_stages=[3, 8],
#         outs=['distillate', 'bottoms_product'],
#         full_condenser=False,
#         reflux=1.0,
#         boilup=2.0,
#         use_cache=True,
#         LHK=('MethylAcetate', 'Water'),
#         stage_reactions=stage_reactions,
#         method='fixed-point',
#         maxiter=100,
#     )
#     distillation.simulate()
#     flows = [
#         [6.585111268569887e-06, 3.6465771244802687e-09,
#          3.110437134095077, 8.015523262914513],
#         [0.019007540635984595, 1.9809858706061696,
#          6.908576991652181, 1.9654626113382356],
#     ]
#     for i, j in zip(distillation.outs, flows):    
#         assert_allclose(i.mol, j, rtol=1e-5, atol=1e-3)
#     sys = bst.System.from_units(units=[distillation])
#     sys._setup()
#     sys.run_phenomena()
#     for i, j in zip(distillation.outs, flows):    
#         assert_allclose(i.mol, j, rtol=1e-5, atol=1e-3)
    
    
# def test_lactic_acid_butanol_reactive_distillation():
#     import biosteam as bst
#     from thermosteam.constants import R
#     from math import exp
#     from numpy.testing import assert_allclose
#     bst.settings.set_thermo(['Water', 'LacticAcid', 'ButylLactate', 'Butanol'], cache=True)
    
#     class Esterification(bst.KineticReaction):
        
#         def volume(self, stream): # kg of catalyst
#             rho_cat = 770 # kg / m3
#             liquid_volume = self.liquid_volume
#             catalyst_volume = 0.5 * liquid_volume
#             catalyst_mass = catalyst_volume * rho_cat
#             return catalyst_mass
        
#         def rate(self, stream): # kmol/kg-catalyst/hr
#             T = stream.T
#             # if T > 370: return 0 # Prevents multiple steady states
#             kf = 2.59e4 * exp(-5.340e4 / (R * T))
#             kr = 3.80e3 * exp(-5.224e4 / (R * T))
#             H2O, LA, BuLA, BuOH = stream.mol / stream.F_mol
#             return 3600 * (kf * LA * BuOH - kr * BuLA * H2O) # kmol / kg-catalyst / hr
        
#     stream = bst.Stream(
#         H2O=10, Butanol=10, LacticAcid=2, T=355,
#     )
#     stream.vle(V=0, P=101325)
#     stage_reactions = {
#         i: Esterification('LacticAcid + Butanol -> Water + ButylLactate', reactant='LacticAcid')
#         for i in range(2)
#     }
#     distillation = bst.MESHDistillation(
#         N_stages=5,
#         ins=[stream],
#         feed_stages=[2],
#         outs=['distillate', 'bottoms_product'],
#         full_condenser=False,
#         reflux=1.0,
#         boilup=2.0,
#         use_cache=True,
#         LHK=('Butanol', 'ButylLactate'),
#         stage_reactions=stage_reactions,
#         method='fixed-point',
#         maxiter=100,
#     )
#     distillation.simulate()
#     flows = [
#         [6.585111268569887e-06, 3.6465771244802687e-09,
#          3.110437134095077, 8.015523262914513],
#         [0.019007540635984595, 1.9809858706061696,
#          6.908576991652181, 1.9654626113382356],
#     ]
#     for i, j in zip(distillation.outs, flows):    
#         assert_allclose(i.mol, j, rtol=1e-5, atol=1e-3)
#     sys = bst.System.from_units(units=[distillation])
#     sys._setup()
#     sys.run_phenomena()
#     for i, j in zip(distillation.outs, flows):    
#         assert_allclose(i.mol, j, rtol=1e-5, atol=1e-3)
    
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
        outs=['distillate', 'bottoms_product'],
        full_condenser=True,
        reflux=1.0,
        boilup=3.5,
        use_cache=True,
        LHK=('Water', 'AceticAcid'),
    )
    distillation.simulate()
    flows = [
        [22.055922550575424, 0.1438997998268773, 87.15249997005776],
        [0.11207744942457795, 4.350500200173122, 22.35850002994225],
    ]
    for i, j in zip(distillation.outs, flows):    
        assert_allclose(i.mol, j, rtol=1e-3, atol=1e-3)
    
if __name__ == '__main__':
    test_multi_stage_adiabatic_vle()
    test_distillation()
    # test_esterification_column()


