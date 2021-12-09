# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import biosteam as bst
from biosteam.utils import colors as c, CABBI_colors
from .utils import style_axis, style_plot_limits, fill_plot, set_axes_labels
from thermosteam.units_of_measure import format_units, reformat_units
import matplotlib.patches as mpatches
from math import floor, ceil
from matplotlib.ticker import MultipleLocator

__all__ = (
    'rounted_tickmarks_from_range',
    'rounded_tickmarks_from_data',
    'plot_unit_groups',
    'plot_unit_groups_across_coordinate',
    'plot_montecarlo', 
    'plot_montecarlo_across_coordinate',
    'plot_scatter_points', 
    'plot_single_point_sensitivity',
    'plot_spearman', 
    'plot_spearman_1d',
    'plot_spearman_2d',
    'plot_horizontal_line',
    'plot_bars', 
    'plot_vertical_line', 
    'plot_scatter_points',
    'plot_contour_1d', 
    'plot_contour_2d', 
    'plot_contour_single_metric',
    'plot_contour_across_coordinate',
    'plot_contour_2d_curves',
    'plot_heatmap',
)

# %% Utilities

default_light_color = c.brown_tint.RGBn
default_dark_color = c.brown_shade.RGBn

def plot_horizontal_line(y, color='grey', **kwargs): # pragma: no coverage
    """Plot horizontal line."""
    plt.axhline(y=y, color=color, **kwargs) 

def plot_vertical_line(x, color='grey', **kwargs): # pragma: no coverage
    """Plot vertical line."""
    plt.axvline(x=x, color=color, **kwargs) 

def plot_scatter_points(xs, ys, color=None, s=50, zorder=1e6, edgecolor='black', marker='o', **kwargs): # pragma: no coverage
    """Plot scatter points and return patch artist."""
    if xs is None: xs = tuple(range(len(ys)))
    if color is None: color = default_dark_color
    return plt.scatter(xs, ys, marker=marker, s=s, color=color,
                       zorder=zorder, edgecolor=edgecolor, **kwargs) 

def rounded_tickmarks_from_data(data, N_ticks, step_min=None, 
                                lb_max=None, ub_min=None, expand=None, f=None,
                                center=None, lb_min=None, ub_max=None):
    get_max = lambda x: max([i.max() for i in x]) if isinstance(x, list) else x.max()
    get_min = lambda x: min([i.min() for i in x]) if isinstance(x, list) else x.min()
    lb = min([get_min(i) for i in data])
    ub = max([get_max(i) for i in data])
    return rounted_tickmarks_from_range(lb, ub, N_ticks, step_min, lb_max, ub_min, expand, f, center,
                                        lb_min, ub_max)

def rounted_tickmarks_from_range(lb, ub, N_ticks, step_min=None, lb_max=None, ub_min=None,
                                 expand=None, f=None, center=None, lb_min=None, ub_max=None):
    if lb_max is not None: lb = min(lb, lb_max)
    if expand is None: expand = 0.10
    diff = expand * (ub - lb)
    ub += diff
    if ub_min is not None: ub = max(ub, ub_min)
    if ub_max is not None: ub = min(ub, ub_max)
    if lb_min is not None: lb = max(lb, lb_min)
    return rounded_linspace(lb, ub, N_ticks, step_min, f, center)

def rounded_linspace(lb, ub, N, step_min, f=None, center=None):
    if step_min is not None:
        lb = floor(lb / step_min) * step_min
        ub = ceil(ub / step_min) * step_min
    step = (ub - lb) / (N - 1)
    if f is None:
        f = int
        if int(step) == step: step = int(step)
        if int(lb) == lb: lb = int(lb)
    else:
        step = f(step)
        lb = f(lb)
    values = [0, 1] if step == 0 else [lb + step * i for i in range(N)]
    if center is not None:
        offset = min(values, key=lambda x: abs(center - x))
        values = [i - offset for i in values[0:-1]]
        values = [values[0] - step, *values, values[-1] + step]
    return [f(i) for i in values]
        
def default_colors_and_hatches(length, colors, hatches):
    if colors is None:
        colors = [i.HEX for i in CABBI_colors]
    if hatches is None:
        hatches = ['', 'x', '-', '/', '\\', '+', '/|', r'\\',
                   '//', '\\|', '.', 'o', '*']
    N_hatches = len(hatches)
    N_colors = len(colors)
    colors *= int(np.ceil(length / N_colors))
    hatches *= int(np.ceil(length / N_hatches))
    return colors, hatches

def subplots(N_axes):
    if N_axes % 3 == 0:
        N_rows = int(N_axes / 3)
        N_cols = 3
    elif N_axes % 2 == 0:
        N_rows = int(N_axes / 2)
        N_cols = 2
    else:
        N_rows = N_axes
        N_cols = 1
    return plt.subplots(nrows=N_rows, ncols=N_cols)

def contour_subplots(N_rows, N_cols, single_colorbar=False):
    widths = np.ones(N_cols + 1)
    widths[-1] /= 4
    gs_kw = dict(width_ratios=widths)
    fig, axes = plt.subplots(nrows=N_rows, ncols=N_cols + 1, gridspec_kw=gs_kw)
    axes = axes.reshape([N_rows, N_cols + 1])
    if single_colorbar:
        gs = axes[0, 0].get_gridspec()
        for ax in axes[:, -1]: ax.remove()
        ax_colorbar = fig.add_subplot(gs[:, -1])
        return fig, axes, ax_colorbar
    else:
        return fig, axes
    
def modify_stacked_bars(axes, N_marks, names, colors, hatches, legend=True, 
                        loc='upper left', bbox_to_anchor=(1.05, 1),
                        labelspacing=1.5, handlelength=4, scale=1.0, **kwargs):
    for ax in axes.flatten(): 
        plt.sca(ax)
        lg = ax.get_legend()
        if lg is not None: lg.remove()
        bars = ax.patches
        color_sets = sum([[i]*N_marks for i in colors], [])
        hatch_sets = sum([[i]*N_marks for i in hatches], [])
        for bar, hatch, color in zip(bars, hatch_sets, color_sets):
            bar.set_facecolor(color)
            bar.set_hatch(hatch)
    if legend:
        if axes.ndim == 1: 
            ax = axes[-1]
        elif axes.ndim == 2:
            ax = axes[0, -1]
        else:
            raise ValueError('axes dimensions must 1 or 2')
        plt.sca(ax)
        patches = [mpatches.Patch(facecolor=i, hatch=j, label=k, edgecolor='k') for i,j,k in 
               zip(colors, hatches, names)]
        patches.reverse()
        leg = plt.legend(
            handles=patches, loc=loc, bbox_to_anchor=bbox_to_anchor,
            labelspacing=labelspacing, handlelength=handlelength,
            **kwargs,
        )
        leg.get_frame().set_linewidth(0.0)
        for patch in leg.get_patches():
            patch.set_height(22 * scale)
            patch.set_width(22 * scale)
            patch.set_y(-6 * scale)


# %% General

def plot_bars(scenarios, ys, colors, edgecolors, labels, positions=None): # pragma: no coverage
    barwidth = 0.50
    N_scenarios = len(scenarios)
    N_labels = len(labels)
    if positions is None: positions = N_labels * np.arange(N_scenarios, dtype=float)
    data = (ys, colors, edgecolors, labels)
    for y, color, edgecolor, label in zip(*data):
        plt.bar(positions, y, barwidth,
                align='center', label=label,
                color=color, edgecolor=edgecolor)
        positions += barwidth
    plt.xticks(positions-barwidth*(N_labels+1)/2, scenarios)
    plt.tight_layout()
    plt.legend()

def plot_heatmap(
        data, ax=None, cell_labels=None,
        metric_bar=None, ax_cbar=None, xlabels=None, ylabels=None, **kwargs
    ):
    if ax is None: fig, ax = plt.subplots()
    else: fig = ax.figure    
    if ax_cbar is None: ax_cbar = ax
    if metric_bar is None:
        cbar = None
        cmap = None
    else:
        cmap = metric_bar.cmap
    data = np.asarray(data)
    nrows, ncols = data.shape
    im = ax.imshow(data, cmap=cmap, **kwargs)
    if metric_bar is not None: cbar = metric_bar.colorbar(fig, ax_cbar, im)
    if cell_labels is not None:
        for i in range(nrows):
            for j in range(ncols):
                ax.text(j, i, cell_labels[i, j], ha="center", va="center", color="k")
    xticks = np.arange(ncols)
    yticks = np.arange(nrows)
    style_axis(ax, xticks, yticks, xlabels, ylabels, offset_xticks=True, 
               offset_yticks=True)
    ax.set_aspect('auto')
    return im, cbar

# %% Plot unit groups

def plot_unit_groups_across_coordinate(f, x, name, unit_groups,
        colors=None, hatches=None, fraction=False, legend_kwargs=None,
        **kwargs):
    df = bst.UnitGroup.df_from_groups_across_coordinate(unit_groups, f, x)
    metrics = unit_groups[0].metrics
    N_metrics = len(metrics)
    fig, axes = plt.subplots(nrows=N_metrics, ncols=1)
    axes_flat = axes.flatten()
    for i in range(N_metrics):
        metric = metrics[i]
        col = metric.name_with_units
        df_metric = df[col]
        ax = axes_flat[i]
        plt.sca(ax)
        df_metric.T.plot(kind='bar', stacked=True, edgecolor='k',
                         ax=ax, **kwargs)
        plt.ylabel(reformat_units(col).replace(' [', '\n['))
    fig.align_ylabels(axes)
    xticks = list(range(len(x)))
    for ax in axes:
        plt.sca(ax)
        plt.xticks(xticks, (), rotation=0)
    plt.sca(ax)
    plt.xticks(xticks, x, rotation=0)
    plt.xlabel(name)
    
    data = [df[i.name_with_units].values for i in metrics]
    data_ub = [np.where(i > 0, i, 0.) for i in data]
    data_lb = [np.where(i < 0, i, 0.) for i in data]
    ubs = np.array([i.sum(axis=0).max() for i in data_ub])
    lbs = np.array([i.sum(axis=0).min() for i in data_lb])
    tickmarks = [rounted_tickmarks_from_range(lb, ub, 5, 5, 0., f=int) for lb, ub in
                 zip(lbs, ubs)]
    for i in range(N_metrics):
        ax = axes_flat[i]
        plt.sca(ax)
        yticks = tickmarks[i]
        plt.ylim([yticks[0], yticks[-1]])
        if yticks[0] < 0:
            plot_horizontal_line(0, zorder=0,
                color=CABBI_colors.black.RGBn, linestyle='--'
            )
        style_axis(ax,  
            yticks=yticks,
            ytickf=False,
        )
    if legend_kwargs is None: legend_kwargs = {}
    modify_stacked_bars(axes, len(x), [i.name for i in unit_groups],
                        colors, hatches, legend_kwargs)
    plt.subplots_adjust(hspace=0.1, wspace=0.4)

def plot_unit_groups(unit_groups, colors=None,
                     hatches=None, fraction=False, joint_group=None, 
                     format_total=None, bold_label=False, legend=True, 
                     legend_kwargs=None, **kwargs):
    """Plot unit groups as a stacked bar chart."""
    colors, hatches = default_colors_and_hatches(len(unit_groups), colors, hatches)
    df = bst.UnitGroup.df_from_groups(
        unit_groups, fraction,
    )
    names = [i.name for i in unit_groups]
    if fraction:
        if joint_group is None:
            units = sum([i.units for i in unit_groups], [])
            joint_group = bst.UnitGroup(None, units)
            joint_group.autofill_metrics()
        N_metrics = len(joint_group.metrics)
        if format_total is None: format_total = lambda x: format(x, '.3g')
        if bold_label:
            bar_labels = [r"$\mathbf{" f"{format_total(i())}" "}$" "\n"
                           "$\mathbf{[" f"{format_units(i.units, '', False)}" "]}$"
                          for i in joint_group.metrics]
        else:
            bar_labels = [f"{format_total(i())}\n[{format_units(i.units)}]"
                          for i in joint_group.metrics]
        # bar_labels = [r"$\mathbf{" + i + "}$" for i in bar_labels]
        df.T.plot(kind='bar', stacked=True, edgecolor='k', **kwargs)
        locs, labels = plt.xticks()
        plt.xticks(locs, ['\n['.join(i.get_text().split(' [')) for i in labels])
        if legend: plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
        plt.xticks(rotation=0)
        fig = plt.gcf()
        ax = plt.gca()
        ax.set_ylabel('Cost and Utility Breakdown [%]')
        values = df.values
        negative_values = np.where(values < 0., values, 0.).sum(axis=0)
        lb = min(0., 20 * floor(negative_values.min() / 20))
        plt.ylim(lb, 100)
        plot_horizontal_line(0, color=CABBI_colors.black.RGBn, linestyle='--')
        style_axis(top=False, yticks=np.arange(lb, 101, 20))
        xticks, _ = plt.xticks()
        xlim = plt.xlim()
        y_twin = ax.twiny()
        plt.sca(y_twin)
        y_twin.tick_params(axis='x', top=True, direction="in", length=0)
        y_twin.zorder = 2
        plt.xlim(xlim)
        if len(xticks) != len(bar_labels): xticks = xticks[1:]
        plt.xticks(xticks, bar_labels, va='baseline')
        N_marks = N_metrics
    else:
        metrics = unit_groups[0].metrics
        N_metrics = len(metrics)
        fig, axes = subplots(N_metrics)
        N_rows, N_cols = axes.shape
        axes_flat = axes.flatten()
        for i, col in enumerate(df):
            ax = axes_flat[i]
            plt.sca(ax)
            data = df[[col]]
            data.T.plot(kind='bar', stacked=True, edgecolor='k', ax=ax, **kwargs)
            plt.ylabel(reformat_units(col))
            plt.xticks((), (), rotation=0)
        plt.subplots_adjust(hspace=0.1, wspace=0.4)
        for i in range(N_cols): fig.align_ylabels(axes[:, i])
        N_marks = 1
    axes = np.array([ax])
    if legend_kwargs is None: legend_kwargs = {}
    modify_stacked_bars(axes, N_marks, names, colors, hatches, legend, **legend_kwargs)
    return fig, axes
    
# %% Sensitivity analysis

def plot_spearman(rhos, top=None, name=None, color_wheel=None, index=None):
    if isinstance(rhos, list): 
        return plot_spearman_2d(rhos, top, name, color_wheel=color_wheel, index=index)
    else:
        return plot_spearman_1d(rhos, top, name, color=color_wheel, index=index)

def format_spearman_plot(ax, index, name, yranges):
    plot_vertical_line(0, color=c.neutral_shade.RGBn, lw=1)
    ax.xaxis.set_major_locator(MultipleLocator(0.5))
    ax.xaxis.set_major_formatter('{x:.2f}')
    ax.xaxis.set_minor_locator(MultipleLocator(0.25))
    yticks = [i[0]+i[1]/2 for i in yranges]
    ax.set_xlim(-1, 1)
    if name: ax.set_xlabel(f"Spearman's correlation with {name}")
    ax.set_yticks(yticks)
    ax.tick_params(axis='y', right=False, direction="inout", length=4)
    ax.tick_params(which='both', axis='x', direction="inout", length=4)
    ax.set_yticklabels(index)
    ax.grid(False)
    ylim = plt.ylim()
    
    ax2 = ax.twinx()
    plt.sca(ax2)
    plt.yticks(yticks, [])
    plt.ylim(*ylim)
    ax2.zorder = 1000
    ax2.tick_params(direction="in")
    
    ax3 = ax.twiny()
    plt.sca(ax3)
    ax3.tick_params(which='both', direction="in", labeltop=False, bottom=False, length=2)
    # ax3.xaxis.set_major_locator(MultipleLocator(0.5))
    # ax3.xaxis.set_major_formatter('{x:.2f}')
    # ax3.xaxis.set_minor_locator(MultipleLocator(0.25))
    ax3.zorder = 1000

def format_single_point_sensitivity_plot(center, diff, ax, index, name, yranges):
    plot_vertical_line(center, color=c.neutral_shade.RGBn, lw=1)
    yticks = [i[0]+i[1]/2 for i in yranges]
    yranges = np.array(yranges)
    ax.set_xlim(center - diff, center + diff)
    ax.set_xlabel(name)
    ax.set_yticks(yticks)
    ax.tick_params(axis='y', right=False, direction="inout", length=4)
    ax.tick_params(which='both', axis='x', direction="inout", length=4)
    ax.set_yticklabels(index)
    ax.grid(False)
    ylim = plt.ylim()
    
    ax2 = ax.twinx()
    plt.sca(ax2)
    plt.yticks(yticks, [])
    plt.ylim(*ylim)
    ax2.zorder = 1000
    ax2.tick_params(direction="in")
    
    ax3 = ax.twiny()
    plt.sca(ax3)
    ax3.tick_params(which='both', direction="in", labeltop=False, bottom=False, length=2)
    # ax3.xaxis.set_major_locator(MultipleLocator(0.5))
    # ax3.xaxis.set_major_formatter('{x:.2f}')
    # ax3.xaxis.set_minor_locator(MultipleLocator(0.25))
    ax3.zorder = 1000

def plot_single_point_sensitivity(baseline, lb, ub, 
        top=None, name=None, colors=None, w=1., s=1., offset=0., style=True, 
        fig=None, ax=None, sort=True, index=None
    ): # pragma: no coverage
    """
    Display Spearman's rank correlation plot.
    
    Parameters
    ----------
    baseline : float
        Baseline metric value.
    lb : pandas.Series
        Metric values at parameter lower bounds.
    bb : pandas.Series
        Metric values at parameter upper bounds.
    top=None : float, optional
        Number of parameters to plot (from highest values).
    
    Returns
    -------
    fig : matplotlib Figure
    ax : matplotlib AxesSubplot
    """
    if isinstance(lb, pd.Series):
        if index is None: index = lb.index
        if name is None: name = lb.name
        lb = lb.values
    if isinstance(ub, pd.Series):
        if index is None: index = ub.index
        if name is None: name = ub.name
        ub = ub.values
    if index is None:
        raise ValueError('must pass index if lb or ub is not a pandas Series object')
    # Sort parameters for plot
    diff = np.abs(ub - lb)
    if sort:
        number = list(sorted(range(len(index)), key=lambda x: diff[x]))
        index = [index[i] for i in number]
        lb = [lb[i] for i in number]
        ub = [ub[i] for i in number]
    if top:
        lb = lb[-top:]
        ub = ub[-top:]
        index = index[-top:]
    
    yranges = [(offset + s*i, w) for i in range(len(index))]
    
    # Plot bars one by one
    if ax is None:
        fig, ax = plt.subplots()
    if colors is None: 
        colors = c.blue_tint.RGBn, c.red_tint.RGBn
    color_left, color_right = colors
    for i, y in enumerate(yranges):
        xlb = [baseline, lb[i] - baseline]
        xub = [baseline, ub[i] - baseline]
        ax.broken_barh([xlb], y, facecolors=color_left,
                       edgecolors=c.blue_dark.RGBn)
        ax.broken_barh([xub], y, facecolors=color_right,
                       edgecolors=c.blue_dark.RGBn)
    plot_vertical_line(baseline)
    if style:
        diff = 1.05 * max(np.abs(ub - baseline).max(), np.abs(lb - baseline).max())
        format_single_point_sensitivity_plot(baseline, diff, ax, index, name, yranges)
    return fig, ax    

def plot_spearman_1d(rhos, top=None, name=None, color=None, 
                     w=1., s=1., offset=0., style=True, 
                     fig=None, ax=None, sort=True, index=None,
                     cutoff=None): # pragma: no coverage
    """
    Display Spearman's rank correlation plot.
    
    Parameters
    ----------
    rhos : pandas.Series
         Spearman's rank correlation coefficients to be plotted.
    top=None : float, optional
        Number of parameters to plot (from highest values).
    
    Returns
    -------
    fig : matplotlib Figure
    ax : matplotlib AxesSubplot
    """
    # Sort parameters for plot
    abs_ = abs
    if isinstance(rhos, pd.Series):
        if index is None: index = rhos.index
        if name is None: name = rhos.name
        rhos = rhos.values
    if cutoff:
        cutoff_index = np.where(np.any(np.abs(rhos) > cutoff))
        index = index[cutoff_index]
        rhos = rhos[cutoff_index]
    if sort:
        rhos, index = zip(*sorted(zip(rhos, index),
                                  key=lambda x: abs_(x[0])))

    if top:
        rhos = rhos[-top:]
        index = index[-top:]
    
    xranges = [(0, i) for i in rhos]
    yranges = [(offset + s*i, w) for i in range(len(rhos))]
    
    # Plot bars one by one
    if ax is None:
        fig, ax = plt.subplots()
    if color is None: color = c.blue_tint.RGBn
    for x, y in zip(xranges, yranges):
        ax.broken_barh([x], y, facecolors=color,
                       edgecolors=c.blue_dark.RGBn)
    
    if style:
        if index is None:
            raise ValueError('must pass index if rhos is not a pandas Series object')
        format_spearman_plot(ax, index, name, yranges)
    return fig, ax

def plot_spearman_2d(rhos, top=None, name=None, color_wheel=None, index=None,
                     cutoff=None, sort=True): # pragma: no coverage
    """
    Display Spearman's rank correlation plot.
    
    Parameters
    ----------
    rhos : list[pandas.Series]
         Spearman's rank correlation coefficients to be plotted.
    top=None : float, optional
        Number of parameters to plot (from highest values).
    
    Returns
    -------
    fig : matplotlib Figure
    ax : matplotlib AxesSubplot
    """
    if name is None: name = rhos[0].name
    if index is None: index = rhos[0].index
    rhos = list(reversed(rhos))
    values = np.array([i.values for i in rhos])
    indices = list(range(values.shape[1]))
    if cutoff:
        cutoff_index, = np.where(np.any(np.abs(rhos) > cutoff, axis=0))
        indices = [indices[i] for i in cutoff_index]
    if sort:
        rhos_mean = np.abs(values).max(axis=0)
        indices.sort(key=lambda x: rhos_mean[x])
    if top is not None: indices = indices[-top:]
    rhos = [[rho[i] for i in indices] for rho in values]
    index = [index[i] for i in indices]
    N = len(rhos)
    s = N + 1
    if not color_wheel: color_wheel = CABBI_colors.wheel()
    fig, ax = plt.subplots()
    for i, rho in enumerate(rhos):
        plot_spearman_1d(rho, color=color_wheel[N - i - 1].RGBn, s=s, offset=i,
                         fig=fig, ax=ax, style=False, sort=False, top=None)
    # Plot central line
    yranges = [(s/2 + s*i - 1., 1.) for i in range(len(rhos[0]))]
    format_spearman_plot(ax, index, name, yranges)
    return fig, ax

# %% Monte Carlo

def plot_montecarlo(data, 
                    light_color=None,
                    dark_color=None,
                    positions=None,
                    xmarks=None,
                    transpose=None): # pragma: no coverage
    """
    Return box plot of Monte Carlo evaluation.
    
    Parameters
    ----------
    data : numpy.ndarray or pandas.DataFrame
        Metric values with uncertainty. Each row represents a sample and 
        each column represents a metric. 
    light_colors : Iterable(numpy.ndarray)
        RGB normalized to 1. Defaults to brown.
    dark_colors : Iterable(numpy.ndarray)
        RGB normalized to 1. Defaults to brown.
    transpose : bool, optional 
        If True, data will be transposed. If False, data will not be transposed. 
        If no values is given, data will be transposed when the number of columns 
        is greater than the number of rows.
    
    Returns
    -------
    bx : Patch
    
    """
    if isinstance(data, pd.DataFrame): data = data.values
    if transpose is None and data.ndim == 2:
        N_rows, N_cols = data.shape
        if N_cols > N_rows: data = data.transpose()
    elif transpose:
        data = data.transpose()
    if not positions:
        if data.ndim == 1: 
            positions = (0,)
        else:
            positions = list(range(data.shape[1]))
    if light_color is None: light_color = default_light_color
    if dark_color is None: dark_color = default_dark_color
    bx = plt.boxplot(x=data, positions=positions, patch_artist=True,
                     widths=0.8, whis=[5, 95],
                     boxprops={'facecolor':light_color,
                               'edgecolor':dark_color},
                     medianprops={'color':dark_color,
                                  'linewidth':1.5},
                     flierprops = {'marker':'D',
                                   'markerfacecolor': light_color,
                                   'markeredgecolor': dark_color,
                                   'markersize':6})
    if xmarks: plt.xticks(positions, xmarks)
    return bx

def plot_montecarlo_across_coordinate(xs, ys, 
                                      light_color=None,
                                      dark_color=None): # pragma: no coverage
    """
    Plot Monte Carlo evaluation across a coordinate.
    
    Parameters
    ----------
    xs : numpy.ndarray(ndim=1)
        Coordinate values for each column in ``ys``.
    ys : numpy.ndarray(ndim=2)
        Metric values with uncertainty. Each row represents a sample and each 
        column represent a metric along the x-coordinate.
    light_color : numpy.ndarray
        RGB normalized to 1. Defaults to brown.
    dark_color : numpy.ndarray
        RGB normalized to 1. Defaults to brown.
    
    Returns
    -------
    percentiles : numpy.ndarray(ndim=2)
        5, 25, 50, 75 and 95th percentiles by row (5 rows total).
    
    """
    if light_color is None: light_color = default_light_color
    if dark_color is None: dark_color = default_dark_color
    q05, q25, q50, q75, q95 = percentiles = np.percentile(ys, [5,25,50,75,95], axis=0)

    plt.plot(xs, q50, '-',
             color=dark_color,
             linewidth=1.5) # Median
    plt.fill_between(xs, q25, q75,
                     color=light_color,
                     linewidth=1.0)
    plt.plot(xs, q05, '-.',
             color=dark_color,
             linewidth=1.0) # Lower whisker
    plt.plot(xs, q95, '-.',
             color=dark_color,
             linewidth=1.0) # Upper whisker
    
    return percentiles
  
# %% Contours
  
def plot_contour_1d(X_grid, Y_grid, data, 
                    xlabel, ylabel, xticks, yticks, 
                    metric_bars, fillcolor=None, label=False, **styleaxiskw): # pragma: no coverage
    """Create contour plots and return the figure and the axes."""
    n = len(metric_bars)
    assert data.shape == (*X_grid.shape, n), (
        "data shape must be (X, Y, M), where (X, Y) is the shape of both X_grid and Y_grid, "
        "and M is the number of metrics"
    )
    gs_kw = dict(height_ratios=[1, 0.25])
    fig, axes = plt.subplots(ncols=n, nrows=2, gridspec_kw=gs_kw)
    if styleaxiskw is None: styleaxiskw = {}
    cps = np.zeros([n], dtype=object)
    linecolor = c.neutral_shade.RGBn
    for i in range(n):
        metric_bar = metric_bars[i]
        ax = axes[0, i]
        plt.sca(ax)
        style_plot_limits(xticks, yticks)
        yticklabels = i == 0
        xticklabels = True
        if fillcolor is not None: fill_plot(fillcolor)
        cp = plt.contourf(X_grid, Y_grid, data[:, :, i],
                          levels=metric_bar.levels,
                          cmap=metric_bar.cmap)
        if label:
            cs = plt.contour(cp, zorder=1e6,
                             linestyles='dashed', linewidths=1.,
                             norm=metric_bar.norm,
                             levels=metric_bar.levels, colors=[linecolor])
            clabels = ax.clabel(cs, levels=[i for i in cs.levels if i!=metric_bar.levels[-1]], inline=True, fmt=metric_bar.fmt,
                      colors=['k'], zorder=1e6)
            for clabel in clabels: clabel.set_rotation(0)
        cps[i] = cp
        style_axis(ax, xticks, yticks, xticklabels, yticklabels)
        cbar_ax = axes[1, i]
        plt.sca(cbar_ax)
        cb = metric_bar.colorbar(fig, cbar_ax, cp, shrink=0.8, orientation='horizontal')
        plt.axis('off')
    set_axes_labels(axes[:-1], xlabel, ylabel)
    plt.subplots_adjust(hspace=0.1, wspace=0.1)
    return fig, axes, cps, cb

def plot_contour_2d(X_grid, Y_grid, Z_1d, data, 
                    xlabel, ylabel, xticks, yticks, 
                    metric_bars, Z_label=None,
                    Z_value_format=lambda Z: str(Z),
                    fillcolor=None, styleaxiskw=None,
                    label=False): # pragma: no coverage
    """Create contour plots and return the figure and the axes."""
    nrows = len(metric_bars)
    ncols = len(Z_1d)
    assert data.shape == (*X_grid.shape, nrows, ncols), (
        "data shape must be (X, Y, M, Z), where (X, Y) is the shape of both X_grid and Y_grid, "
        "M is the number of metrics, and Z is the number of elements in Z_1d"
    )
    fig, axes = contour_subplots(nrows, ncols)
    if styleaxiskw is None: styleaxiskw = {}
    cps = np.zeros([nrows, ncols], dtype=object)
    cbs = np.zeros([nrows], dtype=object)
    linecolor = c.neutral_shade.RGBn
    for row in range(nrows):
        metric_bar = metric_bars[row]
        for col in range(ncols):
            ax = axes[row, col]
            plt.sca(ax)
            style_plot_limits(xticks, yticks)
            yticklabels = col == 0
            xticklabels = row == nrows - 1
            if fillcolor is not None: fill_plot(fillcolor)
            metric_data = data[:, :, row, col]
            cp = plt.contourf(X_grid, Y_grid, metric_data,
                              levels=metric_bar.levels,
                              cmap=metric_bar.cmap)
            if label:
                cs = plt.contour(cp, zorder=1e16,
                                 linestyles='dashed', linewidths=1.,
                                 levels=cp.levels, colors=[linecolor])
                clabels = ax.clabel(cs, levels=[i for i in cs.levels[::2] if i!=metric_bar.levels[-1]], inline=True, fmt=metric_bar.fmt,
                          colors=['k'], zorder=1e16)
                for i in clabels: i.set_rotation(0)
            cps[row, col] = cp
            style_axis(ax, xticks, yticks, xticklabels, yticklabels, **styleaxiskw)
        cbar_ax = axes[row, -1]
        cbs[row] = metric_bar.colorbar(fig, cbar_ax, cp, shrink=0.8)
        # plt.clim()
    for col in range(ncols):
        if not col and Z_label:
            title = f"{Z_label}: {Z_value_format(Z_1d[col])}"
        else:
            title = Z_value_format(Z_1d[col])
        ax = axes[0, col]
        ax.set_title(title)
    for ax in axes[:, -1]:
        plt.sca(ax)
        plt.axis('off')
    set_axes_labels(axes[:, :-1], xlabel, ylabel)
    plt.subplots_adjust(hspace=0.1, wspace=0.1)
    return fig, axes, cps, cbs
       
def plot_contour_single_metric(X_grid, Y_grid, data, 
                    xlabel, ylabel, xticks, yticks, metric_bar,
                    titles=None, fillcolor=None, styleaxiskw=None,
                    label=False): # pragma: no coverage
    """Create contour plots and return the figure and the axes."""
    *_, nrows, ncols = data.shape
    assert data.shape == (*X_grid.shape, nrows, ncols), (
        "data shape must be (X, Y, M, N), where (X, Y) is the shape of both X_grid and Y_grid"
    )
    fig, axes, ax_colorbar = contour_subplots(nrows, ncols, single_colorbar=True)
    if styleaxiskw is None: styleaxiskw = {}
    cps = np.zeros([nrows, ncols], dtype=object)
    linecolor = c.neutral_shade.RGBn
    for row in range(nrows):
        for col in range(ncols):
            ax = axes[row, col]
            plt.sca(ax)
            style_plot_limits(xticks, yticks)
            yticklabels = col == 0
            xticklabels = row == nrows - 1
            if fillcolor is not None: fill_plot(fillcolor)
            metric_data = data[:, :, row, col]
            cp = plt.contourf(X_grid, Y_grid, metric_data,
                              levels=metric_bar.levels,
                              cmap=metric_bar.cmap,
                              norm=metric_bar.norm)
            for i in cp.collections: i.set_edgecolor('face') # For svg background
            if label:
                cs = plt.contour(cp, zorder=1,
                                 linestyles='dashed', linewidths=1.,
                                 norm=metric_bar.norm,
                                 levels=metric_bar.levels, colors=[linecolor])
                clabels = ax.clabel(cs, levels=[i for i in cs.levels if i!=metric_bar.levels[-1]], inline=True, fmt=metric_bar.fmt,
                          colors=['k'], zorder=1)
                for i in clabels: i.set_rotation(0)
            cps[row, col] = cp
            style_axis(ax, xticks, yticks, xticklabels, yticklabels, **styleaxiskw)
    cb = metric_bar.colorbar(fig, ax_colorbar, cp, fraction=0.5)
    plt.sca(ax_colorbar)
    plt.axis('off')
    if titles:
        for col, title in enumerate(titles):
            ax = axes[0, col]
            ax.set_title(title)
    set_axes_labels(axes[:, :-1], xlabel, ylabel)
    plt.subplots_adjust(hspace=0.1, wspace=0.1)
    return fig, axes, cps, cb

def plot_contour_2d_curves(X_grid, Y_grid, Z_1d, data, 
                    xlabel, ylabel, xticks, yticks, 
                    metric_bars, Z_label=None,
                    Z_value_format=lambda Z: str(Z),
                    fillcolor=None, styleaxiskw=None): # pragma: no coverage
    """Create contour curve plots and return the figure and the axes."""
    nrows = len(metric_bars)
    ncols = len(Z_1d)
    assert data.shape == (*X_grid.shape, nrows, ncols), (
        "data shape must be (X, Y, M, Z), where (X, Y) is the shape of both X_grid and Y_grid, "
        "M is the number of metrics, and Z is the number of elements in Z_1d"
    )
    widths = np.ones(ncols)
    gs_kw = dict(width_ratios=widths)
    fig, axes = plt.subplots(ncols=ncols, nrows=nrows, gridspec_kw=gs_kw)
    axes = axes.reshape([nrows, ncols])
    if styleaxiskw is None: styleaxiskw = {}
    cps = np.zeros([nrows, ncols], dtype=object)
    for row in range(nrows):
        metric_bar = metric_bars[row]
        for col in range(ncols):
            ax = axes[row, col]
            plt.sca(ax)
            style_plot_limits(xticks, yticks)
            yticklabels = col == 0
            xticklabels = row == nrows - 1
            if fillcolor is not None: fill_plot(fillcolor)
            metric_data = data[:, :, row, col]
            cp = plt.contour(X_grid, Y_grid, metric_data,
                              levels=metric_bar.levels,
                              cmap=metric_bar.cmap)
            clabels = ax.clabel(cp, levels=cp.levels, inline=True, fmt=lambda x: f'{round(x):,}',
                      colors=['k'], zorder=1e16)
            for i in clabels: i.set_rotation(0)
            cps[row, col] = cp
            style_axis(ax, xticks, yticks, xticklabels, yticklabels, **styleaxiskw)
    for col in range(ncols):
        if not col and Z_label:
            title = f"{Z_label}: {Z_value_format(Z_1d[col])}"
        else:
            title = Z_value_format(Z_1d[col])
        ax = axes[0, col]
        ax.set_title(title)
    for ax in axes[:, -1]:
        plt.sca(ax)
        plt.axis('off')
    set_axes_labels(axes[:, :-1], xlabel, ylabel)
    plt.subplots_adjust(hspace=0.1, wspace=0.1)
    return fig, axes, cps

def plot_contour_across_coordinate(X_grid, Y_grid, Z_1d, data, 
                                   xlabel, ylabel, xticks, yticks, 
                                   metric_bar, Z_label=None,
                                   Z_value_format=lambda Z: str(Z),
                                   fillcolor=None): # pragma: no coverage
    """Create contour plots and return the figure and the axes."""
    ncols = len(Z_1d)
    assert data.shape == (*X_grid.shape, ncols), (
        "data shape must be (X, Y, Z), where (X, Y) is the shape of both X_grid and Y_grid, "
        "and Z is the number of elements in Z_1d"
    )
    widths = np.ones(ncols + 1)
    widths[-1] *= 0.38196601125
    gs_kw = dict(width_ratios=widths)
    fig, axes = plt.subplots(ncols=ncols + 1, nrows=1, gridspec_kw=gs_kw)
    xticklabels = True
    for col in range(ncols):
        ax = axes[col]
        plt.sca(ax)
        style_plot_limits(xticks, yticks)
        yticklabels = col == 0
        if fillcolor is not None: fill_plot(fillcolor)
        cp = plt.contourf(X_grid, Y_grid, data[:, :, col],
                          levels=metric_bar.levels,
                          cmap=metric_bar.cmap)
        style_axis(ax, xticks, yticks, xticklabels, yticklabels)
    cbar_ax = axes[-1]
    metric_bar.colorbar(fig, cbar_ax, cp, fraction=0.35, pad=0.15)
    for col in range(ncols):
        if not col and Z_label:
            title = f"{Z_label}: {Z_value_format(Z_1d[col])}"
        else:
            title = Z_value_format(Z_1d[col])
        ax = axes[col]
        ax.set_title(title)
    plt.sca(axes[-1])
    style_plot_limits(xticks, yticks)
    plt.axis('off')
    set_axes_labels(axes[np.newaxis, :-1], xlabel, ylabel)
    plt.subplots_adjust(hspace=0.1, wspace=0.1)
    return fig, axes
            
# def plot_contour_across_metric(X_grid, Y_grid, data, 
#                                xlabel, ylabel, xticks, yticks, 
#                                metric_bars, Z_value_format=lambda Z: str(Z),
#                                fillcolor=None):
#     """Create contour plots and return the figure and the axes."""
#     ncols = len(metric_bars)
#     assert data.shape == (*X_grid.shape, ncols), (
#         "data shape must be (X, Y, M), where (X, Y) is the shape of both X_grid and Y_grid, "
#         "and M is the number of metric bars"
#     )
#     widths = np.ones(ncols + 1)
#     widths[-1] /= 4
#     gs_kw = dict(width_ratios=widths)
#     fig, axes = plt.subplots(ncols=ncols + 1, nrows=1, gridspec_kw=gs_kw)
#     xticklabels = True
#     for col in range(ncols):
#         ax = axes[col]
#         plt.sca(ax)
#         style_plot_limits(xticks, yticks)
#         yticklabels = col == 0
#         style_axis(ax, xticks, yticks, xticklabels, yticklabels)
#         if fillcolor is not None: fill_plot(fillcolor)
#         cp = plt.contourf(X_grid, Y_grid, data[:, :, col],
#                           cmap=metric_bar.cmap)
#     cbar_ax = axes[-1]
#     metric_bar.colorbar(fig, cbar_ax, cp)
#     for col in range(ncols):
#         title = Z_value_format(Z_1d[col])
#         ax = axes[col]
#         ax.set_title(title)
#     plt.sca(axes[-1])
#     plt.axis('off')
#     set_axes_labels(axes[:-1], xlabel, ylabel)
#     plt.subplots_adjust(hspace=0.1, wspace=0.1)
#     return fig, axes
