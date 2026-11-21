# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import pytest
import biosteam as bst
import numpy as np
from numpy.testing import assert_allclose

def ethanol_water_column(N_stages=3, ID=None, **kwargs):
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    feed = bst.Stream(Ethanol=80, Water=100, T=80.215 + 273.15)
    kwargs.setdefault('outs', ['', ''])
    kwargs.setdefault('stage_specifications', {0: ('Reflux', 0.673), -1: ('Boilup', 2.57)})
    MSE = bst.MultiStageEquilibrium(
        ID, N_stages=N_stages, ins=[feed], feed_stages=[1],
        phases=('g', 'l'), **kwargs
    )
    return MSE

def check_column(MSE):
    vapor, liquid = MSE.outs[:2]
    feed = MSE.ins[0]
    assert_allclose(
        sum([i.mol for i in MSE.outs]), feed.mol, rtol=1e-2, atol=1e-3
    )
    y = vapor.imol['Ethanol'] / vapor.F_mol
    x = liquid.imol['Ethanol'] / liquid.F_mol
    z = feed.imol['Ethanol'] / feed.F_mol
    assert y > z > x
    Tb_ethanol = 351.3
    Tb_water = 373.2
    for stage in MSE.stages:
        assert Tb_ethanol - 2 < stage.T < Tb_water + 2

def test_mse_vle_algorithms():
    bst.F.set_flowsheet('mse_vle_algorithms')
    for algorithm, decomposition, method in [
            ('inside out', 'bubble point', 'fixed-point'),
            ('inside out', 'sum rates', 'wegstein'),
            ('inside out', 'modular', 'fixed-point'),
            ('phenomena', 'bubble point', 'wegstein'),
            ('phenomena', 'sum rates', 'fixed-point'),
            ('phenomena modular', None, 'fixed-point'),
            ('sequential modular', None, 'fixed-point'),
        ]:
        MSE = ethanol_water_column(
            algorithms=algorithm, methods=method,
            vle_decomposition=decomposition, maxiter=3,
        )
        MSE.simulate()
        vapor, liquid = MSE.outs
        assert np.isfinite(vapor.mol).all() and (vapor.mol >= 0).all()
        assert np.isfinite(liquid.mol).all() and (liquid.mol >= 0).all()
        assert np.isfinite([i.T for i in MSE.stages]).all()

def test_mse_simultaneous_correction():
    bst.F.set_flowsheet('mse_simultaneous_correction')
    for method in ('hybr', 'trf'):
        MSE = ethanol_water_column(
            algorithms=('simultaneous correction',), methods=method, maxiter=20,
        )
        MSE.simulate()
        vapor, liquid = MSE.outs
        assert np.isfinite(vapor.mol).all()
        assert np.isfinite(liquid.mol).all()
    MSE = ethanol_water_column(
        algorithms=('simultaneous correction',), methods='invalid',
    )
    with pytest.raises(ValueError):
        MSE.simulate()

def test_mse_invalid_inputs():
    bst.F.set_flowsheet('mse_invalid_inputs')
    MSE = ethanol_water_column(algorithms='phenomena', methods='invalid')
    with pytest.raises(ValueError):
        MSE.simulate()
    with pytest.raises(RuntimeError):
        bst.MultiStageEquilibrium(
            N_stages=2, ins=[bst.Stream(Water=1)], feed_stages=[0],
            outs=['', '', ''], phases=('g', 'l'),
        )

def test_mse_default_and_utilities():
    bst.F.set_flowsheet('mse_default_and_utilities')
    MSE = ethanol_water_column(maxiter=20)
    MSE.simulate()
    check_column(MSE)
    df = MSE.stream_table()
    assert df.shape[0] == MSE.N_stages
    MSE.correct_overall_mass_balance()
    assert_allclose(
        sum([i.mol for i in MSE.outs]), MSE.ins[0].mol, rtol=1e-6, atol=1e-6
    )
    assert np.isfinite(MSE.residual())
    assert MSE.composition_sensitive is False
    stage = MSE.stages[0]
    assert stage.vapor is stage.outs[0]
    assert stage.liquid is stage.outs[1]
    assert stage.top is stage.outs[0]
    assert stage.bottom is stage.outs[1]
    assert stage.extract is stage.outs[0]
    assert stage.raffinate is stage.outs[1]
    assert stage.vapor_side_draw is None
    assert stage.liquid_side_draw is None
    assert stage.top_side_draw is None
    assert stage.bottom_side_draw is None
    assert stage.extract_side_draw is None
    assert stage.raffinate_side_draw is None
    assert len(stage.K) == 2 and (np.asarray(stage.K) > 0).all()
    assert_allclose(np.sum(stage.x), 1, rtol=1e-2)
    assert_allclose(np.sum(stage.y), 1, rtol=1e-2)
    # Empty feed
    MSE.ins[0].empty()
    MSE.simulate()
    assert all([i.isempty() for i in MSE.outs])

def test_mse_temperature_and_duty_specifications():
    bst.F.set_flowsheet('mse_temperature_and_duty_specifications')
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    feed = bst.Stream(Ethanol=80, Water=100, T=355)
    MSE = bst.MultiStageEquilibrium(
        N_stages=3, ins=[feed], feed_stages=[1], outs=['', ''],
        phases=('g', 'l'), T=355, maxiter=10,
        algorithms=('phenomena',)
    )
    MSE.simulate()
    for stage in MSE.stages: assert_allclose(stage.T, 355, rtol=1e-3)
    feed = bst.Stream(Ethanol=80, Water=100, T=355)
    MSE = bst.MultiStageEquilibrium(
        N_stages=3, ins=[feed], feed_stages=[1], outs=['', ''],
        phases=('g', 'l'), maxiter=10,
        stage_specifications={0: ('Temperature', 352), 2: ('Duty', 1e6)},
        algorithms=('phenomena',)
    )
    MSE.simulate()
    assert np.isfinite([i.T for i in MSE.stages]).all()

def test_mse_flow_specification_and_pressure_profile():
    bst.F.set_flowsheet('mse_flow_specification_and_pressure_profile')
    MSE = ethanol_water_column(
        stage_specifications={0: ('Reflux', 1), -1: ('Flow', 0.5)},
        P=[101325, 105000], maxiter=20,
    )
    assert len(MSE.P) == 3
    MSE.simulate()
    vapor, liquid = MSE.outs
    assert_allclose(
        vapor.mol + liquid.mol, MSE.ins[0].mol, rtol=1e-2, atol=1e-3
    )

def test_mse_side_draws():
    bst.F.set_flowsheet('mse_side_draws')
    MSE = ethanol_water_column(
        N_stages=4,
        outs=['', '', '', ''],
        top_side_draws={1: 0.2}, bottom_side_draws=[(2, 0.2)], maxiter=20,
    )
    MSE.simulate()
    assert_allclose(
        sum([i.mol for i in MSE.outs]), MSE.ins[0].mol, rtol=1e-2, atol=1e-3
    )
    stage = MSE.stages[1]
    assert stage.vapor_side_draw is stage.outs[2]
    assert stage.top_side_draw is stage.outs[2]
    stage = MSE.stages[2]
    assert stage.liquid_side_draw is stage.outs[-1]
    assert stage.bottom_side_draw is stage.outs[-1]

def test_mse_lle():
    bst.F.set_flowsheet('mse_lle')
    bst.settings.set_thermo(['Water', 'Methanol', 'Octanol'], cache=True)
    feed = bst.Stream('feed', Water=500, Methanol=50)
    solvent = bst.Stream('solvent', Octanol=500)
    MSE = bst.MultiStageEquilibrium(
        N_stages=2, ins=[feed, solvent], phases=('L', 'l'), maxiter=20,
    )
    MSE.simulate()
    extract, raffinate = MSE.outs
    assert MSE.composition_sensitive
    assert extract.imol['Methanol'] / feed.imol['Methanol'] > 0.5
    assert extract.imol['Octanol'] / solvent.imol['Octanol'] > 0.9
    assert_allclose(
        extract.mol + raffinate.mol, feed.mol + solvent.mol, rtol=1e-2, atol=1e-3
    )
    # Stage-wise partition coefficients
    feed = bst.Stream('feed', Water=500, Methanol=50)
    solvent = bst.Stream('solvent', Octanol=500)
    K = np.array([[1.38], [1.38]])
    MSE = bst.MultiStageEquilibrium(
        N_stages=2, ins=[feed, solvent], phases=('L', 'l'),
        partition_data={
            'K': K,
            'IDs': ('Methanol',),
            'raffinate_chemicals': ('Water',),
            'extract_chemicals': ('Octanol',),
        },
    )
    MSE.simulate()
    extract, raffinate = MSE.outs
    assert extract.imol['Methanol'] / feed.imol['Methanol'] > 0.5
    assert_allclose(extract.imol['Octanol'], 500, rtol=1e-6)

def test_mse_convergence_analysis():
    bst.F.set_flowsheet('mse_convergence_analysis')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    for algorithm, method in [
            (None, None),
            ('simultaneous correction', 'hybr'),
            ('simultaneous correction', 'trf'),
            ('phenomena', 'fixed-point'),
            ('phenomena', 'wegstein'),
            ('sequential modular', None),
            ('phenomena modular', None),
            ('inside out', None),
        ]:
        MSE = ethanol_water_column(maxiter=4)
        profile = MSE.convergence_analysis(
            iterations=3, algorithm=algorithm, method=method, 
            plot=True, verbose=False, yticks=[-5, 0, 5], xticks=[0, 2, 4],
        )
        assert np.isfinite(profile.log_residual).all()
        assert len(profile.iteration) == len(profile.log_residual)
        plt.close('all')
    MSE = ethanol_water_column(maxiter=4)
    for algorithm, method in [
            ('simultaneous correction', 'invalid'),
            ('phenomena', 'invalid'),
        ]:
        with pytest.raises(ValueError):
            MSE.convergence_analysis(
                iterations=2, algorithm=algorithm, method=method, plot=False
            )

def test_single_phase_stage():
    bst.F.set_flowsheet('single_phase_stage')
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    water = bst.Stream(Water=10, T=300)
    ethanol = bst.Stream(Ethanol=5, T=310)
    stage = bst.SinglePhaseStage(ins=[water, ethanol], T=320, P=101325, phase='l')
    stage.simulate()
    outlet = stage.outs[0]
    assert_allclose(outlet.T, 320)
    assert_allclose(outlet.mol, water.mol + ethanol.mol)
    stage = bst.SinglePhaseStage(ins=[water, ethanol], Q=1e5)
    stage.simulate()
    outlet = stage.outs[0]
    assert_allclose(outlet.H, water.H + ethanol.H + 1e5, rtol=1e-3)
    assert outlet.T > 300

def test_phenomena_based_vle_stages():
    bst.F.set_flowsheet('phenomena_based_vle_stages')
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    for spec in [dict(B=1), dict(T=355), dict(Q=0)]:
        with bst.System(algorithm='phenomena based') as sys:
            feed = bst.Stream(Water=50, Ethanol=20, T=300)
            H1 = bst.SinglePhaseStage(ins=feed, T=350, phase='l')
            H2 = bst.SinglePhaseStage(ins=H1-0, Q=0, phase='l')
            S1 = bst.StageEquilibrium(ins=H2-0, phases=('g', 'l'), **spec)
            S2 = bst.StageEquilibrium(ins=S1-1, phases=('g', 'l'), B=0.5)
        sys.simulate()
        sys.run_phenomena()
        products = [S1.outs[0], *S2.outs]
        assert_allclose(
            sum([i.mol for i in products]), feed.mol, rtol=1e-2, atol=1e-3
        )
        for i in products: assert (i.mol >= 0).all()

def test_phenomena_based_lle_stages():
    bst.F.set_flowsheet('phenomena_based_lle_stages')
    bst.settings.set_thermo(['Water', 'Ethanol', 'Octanol'], cache=True)
    for spec in [dict(T=310), dict(Q=0)]:
        with bst.System(algorithm='phenomena based') as sys:
            feed = bst.Stream(Water=50, Ethanol=5, units='kg/hr')
            solvent = bst.Stream(Octanol=50, units='kg/hr')
            solvent2 = bst.Stream(Octanol=20, units='kg/hr')
            S1 = bst.StageEquilibrium(
                phases=('L', 'l'), ins=[feed, solvent], 
                top_chemical='Octanol', **spec
            )
            S2 = bst.StageEquilibrium(
                phases=('L', 'l'), ins=[S1-1, solvent2], top_chemical='Octanol'
            )
        sys.simulate()
        sys.run_phenomena()
        products = [S1.outs[0], *S2.outs]
        assert_allclose(
            sum([i.mol for i in products]), 
            feed.mol + solvent.mol + solvent2.mol, 
            rtol=1e-2, atol=1e-3
        )
        assert S1.outs[0].imol['Octanol'] > S1.outs[1].imol['Octanol']

def test_mse_lle_algorithms():
    bst.F.set_flowsheet('mse_lle_algorithms')
    bst.settings.set_thermo(['Water', 'Methanol', 'Octanol'], cache=True)
    for kwargs in [
            dict(algorithms='sequential modular'),
            dict(algorithms='inside out'), # Defaults to phenomena
            dict(algorithms='phenomena', methods='wegstein'),
            dict(algorithms='phenomena', T=300),
            dict(algorithms='phenomena', stage_specifications={1: ('Temperature', 305)}),
        ]:
        feed = bst.Stream(Water=500, Methanol=50)
        solvent = bst.Stream(Octanol=500)
        MSE = bst.MultiStageEquilibrium(
            N_stages=3, ins=[feed, solvent], phases=('L', 'l'), maxiter=5,
            **kwargs
        )
        MSE.simulate()
        extract, raffinate = MSE.outs
        assert extract.imol['Methanol'] / feed.imol['Methanol'] > 0.5
        assert extract.imol['Octanol'] / solvent.imol['Octanol'] > 0.9
        assert raffinate.imol['Water'] > extract.imol['Water']

def test_mse_vle_with_nonvolatile():
    bst.F.set_flowsheet('mse_vle_with_nonvolatile')
    bst.settings.set_thermo(['Water', 'Ethanol', 'Glycerol'], cache=True)
    for kwargs in [
            dict(algorithms='inside out'),
            dict(algorithms='phenomena'),
            dict(algorithms='simultaneous correction'),
            dict(algorithms='phenomena', vle_decomposition='sum rates',
                 stage_specifications={0: ('Reflux', 1), 1: ('Temperature', 360), -1: ('Boilup', 2)}),
        ]:
        feed = bst.Stream(Ethanol=80, Water=100, Glycerol=5, T=355)
        kwargs.setdefault('stage_specifications', {0: ('Reflux', 1), -1: ('Boilup', 2)})
        MSE = bst.MultiStageEquilibrium(
            N_stages=3, ins=[feed], feed_stages=[1], outs=['', ''], 
            phases=('g', 'l'), maxiter=5, **kwargs
        )
        MSE.simulate()
        vapor, liquid = MSE.outs
        for i in MSE.outs: assert np.isfinite(i.mol).all()
        assert np.isfinite([i.T for i in MSE.stages]).all()
    
def test_mse_damping_and_line_search():
    bst.F.set_flowsheet('mse_damping_and_line_search')
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    for name, value in [('damping', 0.3), ('_line_search', True)]:
        MSE = ethanol_water_column(algorithms='phenomena', maxiter=5)
        setattr(MSE, name, value)
        MSE.simulate()
        check_column(MSE)

def test_stage_equilibrium_specifications():
    bst.F.set_flowsheet('stage_equilibrium_specifications')
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    feed = bst.Stream(Water=50, Ethanol=20, T=355)
    stage = bst.StageEquilibrium(ins=feed, phases=('g', 'l'), B=1)
    stage.simulate()
    vapor, liquid = stage.outs
    assert_allclose(vapor.F_mol / liquid.F_mol, 1, rtol=1e-2)
    assert_allclose(vapor.mol + liquid.mol, feed.mol, rtol=1e-6)
    assert vapor.imol['Ethanol'] / vapor.F_mol > liquid.imol['Ethanol'] / liquid.F_mol
    for name, value in [('Q', 1.), ('B', 2.), ('W', 2.), ('F', 0.5), 
                        ('T', 300.), ('P', 1e5), ('K', np.array([1., 2.])),
                        ('x', np.array([0.5, 0.5])), ('y', np.array([0.5, 0.5]))]:
        setattr(stage, name, value)
        assert_allclose(getattr(stage, name), value)
    stage.specify_variables(Reflux=0)
    assert stage.specified_variable == 'B'
    assert stage.B == float('inf')
    stage.specify_variables(Reflux=2, wt=True)
    assert stage.specified_variable == 'W'
    stage.specify_variables(Boilup=2, wt=True)
    assert stage.specified_variable == 'W'
    stage.specify_variables(Pressure=1e5)
    assert stage.specified_variable == 'Q'
    stage.specify_variables(Q=0, P=1e5)
    assert stage.specified_variable == 'Q'
    with pytest.raises(ValueError):
        stage.specify_variables(Q=0, B=1)
    with pytest.raises(ValueError):
        stage.specify_variables(Q=0, B=1, P=1e5)
    with pytest.raises(RuntimeError):
        stage.specify_variables(Foo=1)

if __name__ == '__main__':
    test_mse_vle_algorithms()
    test_mse_simultaneous_correction()
    test_mse_invalid_inputs()
    test_mse_default_and_utilities()
    test_mse_temperature_and_duty_specifications()
    test_mse_flow_specification_and_pressure_profile()
    test_mse_side_draws()
    test_mse_lle()
    test_mse_convergence_analysis()
    test_single_phase_stage()
    test_phenomena_based_vle_stages()
    test_phenomena_based_lle_stages()
    test_mse_lle_algorithms()
    test_mse_vle_with_nonvolatile()
    test_mse_damping_and_line_search()
    test_stage_equilibrium_specifications()
