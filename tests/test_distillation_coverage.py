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

def methanol_water_feed(**kwargs):
    bst.settings.set_thermo(['Water', 'Methanol', 'Glycerol'], cache=True)
    feed = bst.Stream(Water=80, Methanol=100, Glycerol=5, **kwargs)
    feed.T = feed.bubble_point_at_P().T
    return feed

def test_binary_distillation_recovery_specification():
    bst.F.set_flowsheet('binary_distillation_recovery_specification')
    feed = methanol_water_feed()
    D1 = bst.BinaryDistillation(
        ins=feed, LHK=('Methanol', 'Water'), Lr=0.95, Hr=0.95, k=1.5,
        is_divided=True,
    )
    D1.simulate()
    distillate, bottoms = D1.outs
    assert_allclose(distillate.imol['Methanol'], 95, rtol=1e-3)
    assert_allclose(bottoms.imol['Water'], 0.95 * 80, rtol=1e-3)
    assert_allclose(distillate.mol + bottoms.mol, feed.mol, rtol=1e-6)
    assert D1.line == 'Divided Distillation Column'
    assert D1.design_results['Rectifier stages'] > 0
    assert D1.design_results['Stripper stages'] > 0
    assert D1.distillate is distillate
    assert D1.bottoms_product is bottoms
    assert D1.product_specification_format == 'Recovery'
    with pytest.raises(AssertionError):
        D1.y_top = 0.9
    with pytest.raises(AssertionError):
        D1.Lr = 1.5
    with pytest.raises(ValueError):
        D1.Hr = 1.5
    D1.product_specification_format = 'Composition'
    assert D1.Lr is None and D1.Hr is None
    with pytest.raises(AssertionError):
        D1.Lr = 0.9
    with pytest.raises(ValueError):
        D1.Hr = 0.9
    with pytest.raises(AttributeError):
        D1.product_specification_format = 'invalid'
    D1.product_specification_format = 'Recovery'
    with pytest.raises(ValueError):
        bst.BinaryDistillation(
            ins=feed.copy(), LHK=('Methanol', 'Water'), Lr=0.95, k=1.5,
        )
    with pytest.raises(ValueError):
        bst.BinaryDistillation(
            ins=feed.copy(), LHK=('Methanol', 'Water'), Lr=0.95, Hr=0.95, k=1.5,
            product_specification_format='invalid',
        )

def test_binary_distillation_plot_and_rigorous_column():
    bst.F.set_flowsheet('binary_distillation_plot_and_rigorous_column')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    feed = methanol_water_feed()
    D1 = bst.BinaryDistillation(
        ins=feed, LHK=('Methanol', 'Water'), y_top=0.95, x_bot=0.05, k=1.5,
    )
    with pytest.raises(RuntimeError):
        D1._plot_stages()
    D1.simulate()
    D1.plot_stages()
    plt.close('all')
    kwargs = D1.to_rigorous_column(kwargs_only=True, print_code=True)
    assert kwargs['N_stages'] == int(D1.design_results['Theoretical stages'])
    kwargs = D1.to_rigorous_column(kwargs_only=True, dP=0)
    assert kwargs['P'] == D1.P
    D2 = D1.to_rigorous_column(maxiter=3)
    assert isinstance(D2, bst.MESHDistillation)

def test_binary_distillation_infeasible_specifications():
    bst.F.set_flowsheet('binary_distillation_infeasible_specifications')
    # Ethanol-water azeotrope cannot be broken
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    feed = bst.Stream(Water=80, Ethanol=20)
    feed.T = feed.bubble_point_at_P().T
    D1 = bst.BinaryDistillation(
        ins=feed, LHK=('Ethanol', 'Water'), y_top=0.95, x_bot=0.01, k=1.5,
    )
    with pytest.raises(RuntimeError):
        D1.simulate()
    # Intermediate volatile chemical
    bst.settings.set_thermo(['Water', 'Ethanol', 'Methanol'], cache=True)
    feed = bst.Stream(Water=80, Ethanol=50, Methanol=20)
    feed.T = feed.bubble_point_at_P().T
    D1 = bst.BinaryDistillation(
        ins=feed, LHK=('Methanol', 'Water'), Lr=0.9, Hr=0.9, k=1.5,
    )
    with pytest.raises(RuntimeError):
        D1.simulate()
    D1.check_LHK = False
    D1.simulate()
    distillate, bottoms = D1.outs
    assert_allclose(distillate.mol + bottoms.mol, feed.mol, rtol=1e-6)

def test_shortcut_column_extras():
    bst.F.set_flowsheet('shortcut_column_extras')
    feed = methanol_water_feed()
    D1 = bst.ShortcutColumn(
        ins=feed, LHK=('Methanol', 'Water'), y_top=0.95, x_bot=0.05, k=1.5,
    )
    D1.simulate()
    distillate, bottoms = D1.outs
    assert_allclose(distillate.mol + bottoms.mol, feed.mol, rtol=1e-6)
    assert distillate.imol['Methanol'] > bottoms.imol['Methanol']
    with pytest.raises(TypeError):
        D1.plot_stages()
    assert D1.S_node is D1.S_node
    D1.reset_cache()
    D1.simulate()
    assert_allclose(distillate.mol + bottoms.mol, feed.mol, rtol=1e-6)

def test_mesh_distillation_configurations():
    bst.F.set_flowsheet('mesh_distillation_configurations')
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    for kwargs in [
            dict(full_condenser=True, boilup=2),
            dict(full_condenser=True, reflux=1, boilup=2),
            dict(reflux=None, boilup=2, LHK=None),
            dict(reflux=1, bottoms_product_to_feed=0.5),
        ]:
        feed = bst.Stream(Ethanol=20, Water=80, T=360)
        kwargs.setdefault('LHK', ('Ethanol', 'Water'))
        D1 = bst.MESHDistillation(
            ins=feed, N_stages=4, feed_stages=[1], maxiter=10, **kwargs,
        )
        D1.simulate()
        distillate = D1.outs[0]
        bottoms = D1.outs[1]
        assert_allclose(
            sum([i.mol for i in D1.outs]), feed.mol, rtol=1e-2, atol=1e-3
        )
        assert bottoms.imol['Water'] > 0
        assert D1.purchase_cost > 0
    D1.reflux = 2
    assert D1.reflux == 2
    D1.boilup = 3
    assert D1.boilup == 3

def test_mesh_distillation_vacuum():
    bst.F.set_flowsheet('mesh_distillation_vacuum')
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    feed = bst.Stream(Ethanol=20, Water=80, T=330, P=0.3 * 101325)
    D1 = bst.MESHDistillation(
        ins=feed, N_stages=3, feed_stages=[1], reflux=1, boilup=2,
        LHK=('Ethanol', 'Water'), P=0.3 * 101325, maxiter=10,
    )
    with pytest.warns(RuntimeWarning):
        D1.simulate()
    assert_allclose(
        sum([i.mol for i in D1.outs]), feed.mol, rtol=1e-2, atol=1e-3
    )

def test_absorber():
    bst.F.set_flowsheet('absorber')
    bst.settings.set_thermo(['Water', 'Ethanol', 'CO2'], cache=True)
    gas = bst.Stream(CO2=100, Ethanol=5, phase='g', T=310)
    water = bst.Stream(Water=200, T=300)
    A1 = bst.Absorber(
        N_stages=2, ins=[water, gas], solute='Ethanol', outs=['vapor', 'liquid'],
    )
    A1.simulate()
    vapor, liquid = A1.outs
    assert A1.line == 'Absorber'
    assert liquid.imol['Ethanol'] > vapor.imol['Ethanol']
    assert_allclose(
        vapor.mol + liquid.mol, gas.mol + water.mol, rtol=1e-2, atol=1e-3
    )
    A1.reset_cache()
    A1.simulate()
    assert A1.purchase_cost > 0

def test_shortcut_distillation_phenomena_based():
    bst.F.set_flowsheet('shortcut_distillation_phenomena_based')
    for cls in (bst.BinaryDistillation, bst.ShortcutColumn):
        for spec in (dict(y_top=0.95, x_bot=0.05), dict(Lr=0.95, Hr=0.95)):
            bst.settings.set_thermo(['Water', 'Methanol', 'Glycerol'], cache=True)
            with bst.System(algorithm='phenomena based') as sys:
                feed = bst.Stream(Water=80, Methanol=100, Glycerol=5)
                feed.T = feed.bubble_point_at_P().T
                D1 = cls(ins=feed, LHK=('Methanol', 'Water'), k=1.5, **spec)
            sys.simulate()
            distillate, bottoms = D1.outs
            recovery = distillate.imol['Methanol'] / feed.imol['Methanol']
            sys.run_phenomena()
            assert_allclose(
                distillate.imol['Methanol'] / feed.imol['Methanol'], 
                recovery, rtol=1e-2
            )
            assert_allclose(distillate.mol + bottoms.mol, feed.mol, rtol=1e-2, atol=1e-3)
            assert_allclose(bottoms.imol['Glycerol'], 5, rtol=5e-2)

if __name__ == '__main__':
    test_binary_distillation_recovery_specification()
    test_binary_distillation_plot_and_rigorous_column()
    test_binary_distillation_infeasible_specifications()
    test_shortcut_column_extras()
    test_mesh_distillation_configurations()
    test_mesh_distillation_vacuum()
    test_absorber()
    test_shortcut_distillation_phenomena_based()
