# -*- coding: utf-8 -*-
"""
"""
import biosteam as bst
import thermosteam as tmo
import numpy as np
from numpy.testing import assert_allclose

# bst.MultiStageEquilibrium.optimize_result = False

def test_trivial_lle_case():
    import biosteam as bst
    import numpy as np
    from numpy.testing import assert_allclose
    # lle
    with bst.System(algorithm='phenomena based') as sys:
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
    with bst.System(algorithm='phenomena based') as sys:
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
    with bst.System(algorithm='phenomena based') as sys:
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
    assert_allclose(actual, value, rtol=1e-3, atol=1e-3)
    assert_allclose(T_actual, T, rtol=1e-3, atol=1e-3)

def test_trivial_distillation_case():   
    import biosteam as bst
    import thermosteam as tmo
    import numpy as np
    from numpy.testing import assert_allclose
    # distillation
    with bst.System(algorithm='phenomena based') as sys:
        bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
        feed = bst.Stream(Ethanol=80, Water=100, T=353.455)
        MSE = bst.MultiStageEquilibrium(N_stages=5, ins=[feed], feed_stages=[2],
            outs=['vapor', 'liquid'],
            stage_specifications={0: ('Reflux', 0.673), -1: ('Boilup', 2.57)},
            phases=('g', 'l'),
            use_cache=True,
            maxiter=200,
        )
        MSE.molar_tolerance = 1e-9
        MSE.relative_molar_tolerance = 1e-9
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
    po = system(algorithm='phenomena based', 
                    molar_tolerance=1e-9,
                    relative_molar_tolerance=1e-9,
                    method='fixed-point')
    sm = system(algorithm='sequential modular',
                    molar_tolerance=1e-9,
                    relative_molar_tolerance=1e-9,
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
            EtAc_required = broth * solvent_feed_ratio
            if EtAc_required < EtAc_recycle:
                recycle.F_mass *= EtAc_required / EtAc_recycle
                EtAc_recycle = recycle.imass['EthylAcetate']
            EtAc_fresh = EtAc_required - EtAc_recycle
            solvent.imass['EthylAcetate'] = max(
                0, EtAc_fresh
            )
        
        @solvent.material_balance
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
        
    init_sys = system(method='fixed-point')
    init_sys.simulate()
    po = system(algorithm='phenomena based', 
                    molar_tolerance=1e-9,
                    relative_molar_tolerance=1e-9,
                    maxiter=2000,
                    method='fixed-point')
    sm = system(algorithm='sequential modular',
                    molar_tolerance=1e-9,
                    relative_molar_tolerance=1e-9,
                    method='fixed-point')
    time = bst.Timer()
    
    time.start()
    for i in range(1): po.simulate()
    t_phenomena = time.elapsed_time
    
    time.start()
    for i in range(1): sm.simulate()
    t_sequential = time.elapsed_time
    
    for s_sm, s_dp in zip(sm.streams, po.streams):
        actual = s_sm.mol
        value = s_dp.mol
        assert_allclose(actual, value, rtol=1e-6, atol=1e-6)

def test_vlle_case():
    import numpy as np
    import biosteam as bst
    import thermosteam as tmo
    from numpy.testing import assert_allclose
    tmo.settings.set_thermo(['Water', 'Ethanol', 'Octane'], cache=True)
    chemicals = tmo.settings.chemicals
    @bst.SystemFactory
    def system(ins, outs):
        T = 351
        P = 101325
        extract = tmo.Stream('extract')
        raffinate = tmo.Stream('raffinate')
        vapor = tmo.Stream('vapor')
        feed = tmo.Stream('feed', Water=1, Ethanol=0.5, Octane=2, T=T, P=P)
        lle_stage = bst.StageEquilibrium(
            phases=('L', 'l'), 
            ins=[extract, raffinate, feed],
            T=T, P=P,
        )
        vle_L_stage = bst.StageEquilibrium(
            phases=('g', 'l'), 
            ins=[vapor, lle_stage.extract],
            outs=['', extract],
            T=T, P=P,
        )
        vle_l_stage = bst.StageEquilibrium(
            phases=('g', 'l'), 
            ins=[vle_L_stage.vapor, lle_stage.raffinate],
            outs=[vapor, raffinate],
            T=T, P=P,
        )
        total_components = feed.mol.to_array()
        @lle_stage.add_specification(run=True)
        def spec():
            recycles = (extract.mol, raffinate.mol, vapor.mol)
            total_recycles = sum(recycles)
            if total_recycles.all():
                feed.mol[:] = 0
                factor = total_components / total_recycles
                for i in recycles: i *= factor
            else:
                mol = (
                    total_components 
                    - extract.mol
                    - raffinate.mol
                    - vapor.mol
                )
                mol[mol < 0] = 0
                feed.mol = mol
        
        original = lle_stage._create_material_balance_equations
        def spec(composition_sensitive):
            feed.mol[:] = 0
            N = chemicals.size
            e = np.ones(N)
            r = np.ones(N)
            v = np.ones(N)
            new_material_balance = (
                {extract: e,
                 raffinate: r,
                 vapor: v},
                 total_components
            )
            if composition_sensitive:
                return original(composition_sensitive)
            else:
                balances = original(composition_sensitive)
                balances[0] = new_material_balance
                return balances
            
        lle_stage._create_material_balance_equations = spec
    
    init_sys = system(method='fixed-point')
    init_sys.simulate()
    po = system(algorithm='phenomena based', 
                    molar_tolerance=1e-9,
                    relative_molar_tolerance=1e-9,
                    maxiter=2000,
                    method='fixed-point')
    sm = system(algorithm='sequential modular',
                    molar_tolerance=1e-9,
                    relative_molar_tolerance=1e-9,
                    method='fixed-point')
    time = bst.Timer()
    
    time.start()
    for i in range(1): po.simulate()
    t_phenomena = time.elapsed_time
    
    time.start()
    for i in range(1): sm.simulate()
    t_sequential = time.elapsed_time
    
    for s_sm, s_dp in zip(sm.streams, po.streams):
        actual = s_sm.mol
        value = s_dp.mol
        assert_allclose(actual, value, rtol=1e-4, atol=1e-4)
    
    
    # assert_allclose(s.mol, [1, 0.5, 2]) # mass balance
    # total = s.F_mol
    # xl = s.imol['l'].sum() / total
    # xL = s.imol['L'].sum() / total if 'L' in s.phases else 0.
    # xg = s.imol['g'].sum() / total
    # assert_allclose(xl, 0.2748928033836762, atol=2e-3, rtol=2e-3)
    # assert_allclose(xL, 0.6370268977038833, atol=2e-3, rtol=2e-3) # mass balance
    # assert_allclose(xg, 0.08808029891244049, atol=2e-3, rtol=2e-3) # mass balance

if __name__ == '__main__':
    test_trivial_lle_case()
    test_trivial_vle_case()
    test_trivial_liquid_extraction_case()
    test_trivial_distillation_case()
    test_simple_acetic_acid_separation_no_recycle()
    test_simple_acetic_acid_separation_with_recycle()
    test_vlle_case()