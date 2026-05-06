# -*- coding: utf-8 -*-
"""
"""
import biosteam as bst
import thermosteam as tmo
import numpy as np
from numpy.testing import assert_allclose

def test_flash():
    import biosteam as bst
    import numpy as np
    from numpy.testing import assert_allclose
    bst.settings.set_thermo(['Water', 'Ethanol'], ideal=True)
    feed_liquid = bst.Stream('feed_liquid', Water=50, Ethanol=50, phase='l', T = 400)
    feed_gas = bst.Stream('feed_gas', Water=50, Ethanol=50, phase='g', T=400)
    flash_1 = bst.StageEquilibrium('flash_1', ins=[feed_liquid, feed_gas], Q=0, P=101325, phases=('g', 'l'))
    flashes = [flash_1]
    for flash in flashes:
        flash.T = 350
        flash.K = np.array([0.5, 1.5])
        flash.B = 1.
        flash.outs[0].mol[:] = [50, 50]
        flash.outs[1].mol[:] = [50, 50]
    sys = bst.System('sys', path=flashes, algorithm='phenomena based')
    sys.track_convergence()
    for i in range(10):
        sys.run_phenomena()
    phenomena_graph = sys.get_phenomegraph()
    shape = phenomena_graph.variable_profiles.shape
    assert shape == (30, 8)
        

def test_2_stage_flash():
    import biosteam as bst
    import numpy as np
    from numpy.testing import assert_allclose
    bst.settings.set_thermo(['Water', 'Ethanol'], ideal=True)
    feed_liquid = bst.Stream('feed_liquid', Water=50, Ethanol=50, phase='l', T = 400)
    feed_gas = bst.Stream('feed_gas', Water=50, Ethanol=50, phase='g', T=400)
    flash_1 = bst.StageEquilibrium('flash_1', ins=[feed_liquid, feed_gas], Q=0, P=101325, phases=('g', 'l'))
    flash_2 = bst.StageEquilibrium('flash_2', ins=[flash_1-0, feed_liquid], B=1, P=101325, phases=('g', 'l'))
    flashes = [flash_1, flash_2]
    for flash in flashes:
        flash.T = 350
        flash.K = np.array([0.5, 1.5])
        flash.B = 1.
        flash.outs[0].mol[:] = [50, 50]
        flash.outs[1].mol[:] = [50, 50]
    sys = bst.System('sys', path=flashes, algorithm='phenomena based')
    sys.track_convergence()
    for i in range(10):
        sys.run_phenomena()
    phenomena_graph = sys.get_phenomegraph()
    shape = phenomena_graph.variable_profiles.shape
    assert shape == (30, 15)

def test_column():
    import biosteam as bst
    import numpy as np
    from numpy.testing import assert_allclose
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    feed = bst.Stream('feed', Ethanol=80, Water=100, T=80.215 + 273.15)
    D1 = bst.MESHDistillation('D1', N_stages=5, ins=[feed], feed_stages=[2],
        outs=['vapor', 'liquid'],
        reflux=0.673, boilup=2.57,
        LHK=('Ethanol', 'Water') )
    D1.set_flow_rates(D1.hot_start())
    sys = bst.System('sys', path=[D1])
    sys.track_convergence(True)
    for i in range(10):
        sys.run_phenomena()
    phenomena_graph = sys.get_phenomegraph()
    shape = phenomena_graph.variable_profiles.shape
    assert shape == (30, 38)
 
def test_system():
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
    sys = system(algorithm='phenomena oriented', 
                    molar_tolerance=1e-9,
                    relative_molar_tolerance=1e-9,
                    maxiter=2000,
                    method='fixed-point')
    sys.track_convergence(True)
    for i in range(10):
        sys.run_phenomena()
    phenomena_graph = sys.get_phenomegraph()
    shape = phenomena_graph.variable_profiles.shape
    assert shape == (37, 133)
 
if __name__ == '__main__':
    test_flash()
    test_2_stage_flash()
    test_column()
    test_system()