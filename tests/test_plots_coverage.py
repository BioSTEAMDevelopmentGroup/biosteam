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
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import biosteam as bst
from numpy.testing import assert_allclose
from biosteam.plots import plots
from biosteam.plots import sankey
from biosteam.plots.utils import MetricBar

def get_unit_groups(values=((1., 2.), (3., -1.))):
    groups = []
    for n, (a, b) in enumerate(values):
        group = bst.UnitGroup(f'Group {n}', units=[])
        group.metric(lambda a=a: a, 'Cost', 'USD')
        group.metric(lambda b=b: b, 'Duty', 'kW')
        group.metric(lambda a=a: 2 * a, 'Power', 'kW')
        group.metric(lambda b=b: 2 * b, 'Flow', 'kg/hr')
        groups.append(group)
    return groups

def test_rounded_tickmarks():
    ticks = plots.rounded_tickmarks_from_data(
        [np.array([0.2, 3.5]), [np.array([1., 7.4])]], 5,
        step_min=1, lb_max=0, ub_min=8, ub_max=20, lb_min=-1,
    )
    assert len(ticks) == 5
    assert ticks[0] <= 0 and ticks[-1] >= 7.4
    ticks = plots.rounted_tickmarks_from_range(0.1, 0.9, 5, f=lambda x: round(x, 2))
    assert len(ticks) == 5
    ticks = plots.rounded_linspace(0.13, 0.91, 4, p=0.05)
    assert len(ticks) == 4
    assert_allclose(ticks[0], 0.1, rtol=1e-2)
    ticks = plots.rounded_linspace(-3, 7, 3, center=0)
    assert 0 in ticks
    assert plots.rounded_linspace(1, 1, 3) == [0, 1]

def test_colors_hatches_and_subplots():
    colors, hatches = plots.default_colors_and_hatches(20, None, None)
    assert len(colors) >= 20 and len(hatches) >= 20
    colors, hatches = plots.default_colors_and_hatches(3, ['r'], ['x'])
    assert len(colors) == len(hatches) == 3
    for N, shape in [(3, (3,)), (4, (2, 2)), (5, (5,))]:
        fig, axes = plots.subplots(N)
        assert axes.shape == shape
        plt.close('all')
    fig, axes = plots.contour_subplots(2, 2)
    assert axes.shape == (2, 3)
    fig, axes, ax_colorbar = plots.contour_subplots(1, 2, single_colorbar=True, wbar=2)
    assert axes.shape == (1, 3)
    plt.close('all')

def test_modify_stacked_bars():
    fig, axes = plt.subplots(2, 2)
    plots.modify_stacked_bars(axes, 1, ['a'], ['r'], ['x'])
    plt.close('all')

def test_plot_heatmap():
    data = np.arange(6.).reshape([2, 3])
    labels = np.array([['a', 'b', 'c'], ['d', 'e', 'f']])
    im, cbar = plots.plot_heatmap(data, cell_labels=labels,
                                  xlabels=['x1', 'x2', 'x3'], ylabels=['y1', 'y2'])
    assert cbar is None
    fig, ax = plt.subplots()
    metric_bar = MetricBar('Metric', 'USD', plt.cm.viridis, [0, 5])
    im, cbar = plots.plot_heatmap(data, ax=ax, metric_bar=metric_bar)
    assert cbar is not None
    plt.close('all')

def test_plot_stacked_bar():
    data = np.array([[1., 2.], [3., 4.], [0.5, 0.2]])
    names = pd.MultiIndex.from_tuples([('a', 'A'), ('b', 'B'), ('c', 'C')])
    fig, axes = plots.plot_stacked_bar(
        data, names, ['Cost [USD]', 'Duty [kW]'], ylabel='Value [USD/yr]'
    )
    assert axes.shape == (1,)
    plt.close('all')

def test_plot_unit_groups():
    bst.F.set_flowsheet('plot_unit_groups')
    groups = get_unit_groups()
    fig, axes = plots.plot_unit_groups(groups)
    assert axes.size == 1
    plt.close('all')
    joint_group = bst.UnitGroup('Joint', units=[])
    joint_group.metric(lambda: 4., 'Cost', 'USD')
    joint_group.metric(lambda: 1., 'Duty', 'kW')
    joint_group.metric(lambda: 8., 'Power', 'kW')
    joint_group.metric(lambda: 2., 'Flow', 'kg/hr')
    fig, axes = plots.plot_unit_groups(groups, fraction=True, joint_group=joint_group)
    assert axes.size == 1
    plt.close('all')
    fig, axes = plots.plot_unit_groups(groups, fraction=True, joint_group=joint_group,
                                       bold_label=True)
    plt.close('all')

def test_plot_unit_groups_across_coordinate():
    bst.F.set_flowsheet('plot_unit_groups_across_coordinate')
    state = {'x': 1.}
    groups = []
    for n, sign in enumerate((1., -1.)):
        group = bst.UnitGroup(f'Group {n}', units=[])
        group.metric(lambda: 10 * state['x'], 'Cost', 'USD')
        group.metric(lambda sign=sign: 10 * sign * state['x'], 'Duty', 'kW')
        groups.append(group)
    plots.plot_unit_groups_across_coordinate(
        lambda x: state.__setitem__('x', x), [1., 2.], 'X', groups,
        colors=['r', 'b'], hatches=['', 'x'],
    )
    assert len(plt.gcf().axes) >= 2
    plt.close('all')

def test_plot_uncertainty_boxes_and_scatter():
    np.random.seed(0)
    x = np.random.rand(20)
    y = np.random.rand(20)
    boxes = plots.plot_uncertainty_boxes(np.random.rand(20, 2), hatch='x',
                                         xmarks=['a', 'b'])
    assert len(boxes) == 2
    boxes = plots.plot_uncertainty_boxes(x, vertical=False, outliers=True,
                                         light_color=np.array([0.5, 0.5, 0.5]),
                                         dark_color=np.array([0.1, 0.1, 0.1]))
    assert len(boxes) == 1
    plt.close('all')
    box_kwargs = dict(fill=np.array([0.2, 0.3, 0.4]), edge=np.array([0., 0., 0.]))
    for xbox, ybox in [(None, None), (False, None), (None, False), (False, False)]:
        fig, ax, axes = plots.plot_scatter(x, y, xbox=xbox, ybox=ybox,
                                           xbox_kwargs=box_kwargs, ybox_kwargs=box_kwargs,
                                           xlabel='x', ylabel='y')
        assert ax is not None
        plt.close('all')
    fig, ax, axes = plots.plot_scatter(
        [x, x + 1], [y, y + 1], xbox_kwargs=[box_kwargs, box_kwargs],
        ybox_kwargs=[box_kwargs, box_kwargs], xticks=[0, 1, 2], yticks=[0, 1, 2],
    )
    plt.close('all')
    fig, ax = plt.subplots()
    plots.plot_scatter(x, y, ax=ax, fig=fig)
    plt.close('all')
    fig, axes = plt.subplots(2, 2)
    box = plots.Box(axes[0, 0], **box_kwargs)
    assert box.light is box.fill and box.dark is box.edge
    plots.plot_scatter(x, y, axes=axes, fig=fig, xbox=box, ybox=box)
    box.reset()
    plt.close('all')
    fig, ax, axes = plots.plot_kde(x, y)
    plt.close('all')
    fig, ax, axes = plots.plot_kde([x, x], [y, np.ones(20)])
    plt.close('all')
    fig, ax, axes = plots.plot_kde(np.ones(20), y)
    plt.close('all')
    with pytest.raises(ValueError):
        plots.plot_scatter([x, x], [y])
    with pytest.raises(ValueError):
        plots.plot_scatter([x], [y], z=[x, x])
    assert plots.show_box([False, True])
    assert not plots.show_box(0)

def test_plot_scatter_with_z():
    np.random.seed(0)
    x = np.random.rand(10)
    y = np.random.rand(10)
    z = np.random.rand(10)
    fig, ax = plt.subplots()
    plots.plot_scatter(x, y, ax=ax, fig=fig, z=z, zbar=False)
    plt.close('all')

def test_plot_scatter_1d_and_kde_1d():
    np.random.seed(0)
    xs = [np.random.rand(10), np.random.rand(10)]
    ys = [np.random.rand(10), np.random.rand(10)]
    fig, axes = plots.plot_kde_1d(xs, ys, autobox=False)
    assert axes.shape == (1, 2)
    plt.close('all')
    fig, axes = plots.plot_scatter_1d(xs[:1], ys[:1], autobox=False, xlabel='x',
                                      xticklabels=True, xbox_kwargs=[{}])
    assert axes.shape == (1, 1)
    plt.close('all')
    zs = np.random.rand(2, 10)
    fig, axes = plots.plot_scatter_1d(xs, ys, autobox=False, zs=zs, zlabel='z')
    plt.close('all')
    fig, axes = plots.plot_scatter_1d(xs, ys, zs=zs, zlabel='z', xticklabels=False,
                                      xlabel='x', ylabel='y', xticks=[0, 0.5, 1], 
                                      yticks=[0, 0.5, 1])
    plt.close('all')

def test_plot_kde_2d():
    np.random.seed(0)
    xs = np.random.rand(2, 10)
    ys = np.random.rand(2, 10)
    zs = np.random.rand(2, 10)
    fig, axes = plots.plot_kde_2d(xs, ys, zs=zs, xlabel='x', ylabel='y',
                                  xticks=[0, 0.5, 1], yticks=[0, 0.5, 1],
                                  xticklabels=True, yticklabels=True,
                                  xbox_kwargs={}, ybox_kwargs={}, zbar=False)
    assert axes.shape == (1, 2)
    plt.close('all')
    fig, axes = plots.plot_kde_2d(xs, ys, zs=zs, autobox=False, zbar=False)
    assert axes.shape == (1, 2)
    plt.close('all')
    fig, axes = plots.plot_kde_2d(xs, ys, zs=zs, xbox_kwargs=[{}], ybox_kwargs=[{}], zbar=False)
    plt.close('all')

def test_generate_contour_data():
    import os, tempfile
    f = lambda x, y: np.array([x + y, x * y])
    X, Y, Z = plots.generate_contour_data(f, (0, 1), (0, 1), n=3, smooth=False)
    assert X.shape == Y.shape == (3, 3)
    assert Z.shape == (3, 3, 2)
    assert_allclose(Z[..., 0], X + Y, rtol=1e-2)
    with tempfile.TemporaryDirectory() as tmp:
        file = os.path.join(tmp, 'data.npy')
        X, Y, Z = plots.generate_contour_data(
            lambda x, y, a: np.array([a * x, y]), (0, 1), (0, 1), n=3, file=file,
            load=False, smooth=False, args=(2.,), filterwarnings=False,
        )
        X, Y, Z2 = plots.generate_contour_data(f, (0, 1), (0, 1), n=3, file=file)
        assert_allclose(Z, Z2)
    X, Y, Z = plots.generate_contour_data(
        lambda X, Y: np.stack([X, Y], axis=-1), (0, 1), (0, 1), n=3,
        vectorize=False, smooth=False,
    )
    assert Z.shape == (3, 3, 2)

def test_quadrants_and_titles():
    np.random.seed(0)
    x = np.random.rand(10) - 0.5
    y = np.random.rand(10) - 0.5
    fig, ax = plt.subplots()
    plt.xlim(-1, 1); plt.ylim(-1, 1)
    text = ['A()', 'B', 'C()', 'D()']
    labeled = plots.plot_quadrants(text=text, data=(x, y), fs=8)
    assert labeled == [True, True, True, True]
    plt.close('all')
    fig, ax = plt.subplots()
    plt.xlim(-1, 1); plt.ylim(-1, 1)
    text = ['A', 'B()', 'C', 'D']
    labeled = plots.plot_quadrants(
        text=text, data=(x, y), x=(-0.1, 0.1), y=(-0.1, 0.1), rotate=1,
        quadrant_color=['r', 'g', 'b', 'k'], text_color=[None, 'k', None, 'k'],
    )
    assert len(labeled) == 4
    plt.close('all')
    fig, ax = plt.subplots()
    plt.xlim(-1, 1); plt.ylim(-1, 1)
    labeled = plots.plot_quadrants(text=None, x=0.2, y=0.3)
    assert labeled == 4 * [False]
    plt.close('all')
    fig, ax = plt.subplots()
    plt.xlim(-1, 1); plt.ylim(-1, 1)
    labels = plots.label_quadrants(x=None, y=None, xr=0.1, yr=(-0.1, 0.1),
                                   text=['A', 'B', 'C', 'D'], color=4*['k'])
    assert labels == 4 * [True]
    plots.add_titles([ax], ['A'], 'k')
    plt.close('all')

def test_sankey_helpers():
    bst.F.set_flowsheet('sankey_helpers')
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    feed = bst.Stream('feed', Water=10, Ethanol=5, price=0.1)
    M1 = bst.Mixer('M1', [feed, 'feed_2'])
    S1 = bst.Splitter('S1', M1-0, ['product', 'other'], split=0.5)
    group = sankey.asugroup(M1)
    assert isinstance(group, bst.UnitGroup)
    assert sankey.asugroup(group) is group
    with pytest.raises(ValueError):
        sankey.asugroup(1)
    assert sankey.stream_cost(feed) >= 0
    assert sankey.stream_cost(None) == 0.
    assert sankey.abbreviate_word('Ethanol') == 'Eth.'
    assert sankey.abbreviate_word('Gas') == 'Gas'
    assert sankey.abbreviate_word('aeiou') == 'aeiou'
    assert sankey.abbreviate_name('Water treatment') == 'Wat. treatm.'
    group_1 = bst.UnitGroup('Mixing', [M1])
    group_2 = bst.UnitGroup('Splitting', [S1])
    assert sankey.get_unit_group(feed, [group_1, group_2]) is group_1
    s1 = bst.Stream('s_dup', Water=1)
    s2 = bst.Stream('s_dup', Water=2)
    M2 = bst.Mixer('M2', [s1])
    M3 = bst.Mixer('M3', [s2])
    group_3 = bst.UnitGroup('Mixing 2', [M2, M3])
    streams = sankey.reduced_streams([s1, s2, feed], [group_1, group_3])
    assert len(streams) == 2
    assert_allclose(streams[0].imol['Water'], 3, rtol=1e-2)
    M4 = bst.Mixer('M4', [s1.copy('s_dup')])
    streams = sankey.reduced_streams([s1, M4.ins[0]], [group_1, group_3, bst.UnitGroup('Other', [M4])])
    assert len(streams) == 2

if __name__ == '__main__':
    test_rounded_tickmarks()
    test_colors_hatches_and_subplots()
    test_modify_stacked_bars()
    test_plot_heatmap()
    test_plot_stacked_bar()
    test_plot_unit_groups()
    test_plot_unit_groups_across_coordinate()
    test_plot_uncertainty_boxes_and_scatter()
    test_plot_scatter_with_z()
    test_plot_scatter_1d_and_kde_1d()
    test_plot_kde_2d()
    test_generate_contour_data()
    test_quadrants_and_titles()
    test_sankey_helpers()
