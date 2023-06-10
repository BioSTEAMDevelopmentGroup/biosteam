# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
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
from .utils import style_axis, style_plot_limits, fill_plot, set_axes_labels, MetricBar, closest_index
from thermosteam.units_of_measure import format_units, reformat_units
import matplotlib.patches as mpatches
from math import floor, ceil
from matplotlib.ticker import MultipleLocator
from scipy.stats import kde
from collections import deque

__all__ = (
    'rounded_linspace',
    'rounted_tickmarks_from_range',
    'rounded_tickmarks_from_data',
    'annotate_line',
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
    'plot_contour_2d', 
    'plot_contour_single_metric',
    'plot_heatmap',
    'plot_kde_2d',
    'plot_kde',
    'plot_quadrants',
    'plot_stacked_bar',
    'generate_contour_data',
)

# %% Utilities

default_light_color = c.orange_tint.RGBn
default_dark_color = c.orange_shade.RGBn

def annotate_line(text, x, xs, ys, dy=0.2, dy_text=0.22, position='under', 
                  color=None): # pragma: no coverage
    """
    Annotate line with text and arrow pointing to text.
    
    Parameters
    ----------
    text : str
    x : float
        Arrow position
    xs : numpy.ndarray(dim=1)
    ys : numpy.ndarray(dim=1)
    dy : float
        Length of arrow to y-position.
    dy_text : float
        Distance of text to arrow.
    position : {'under', 'over'}
        Relative position of text to line.
    color : numpy.ndarray
        RGB normalized to 1. Defaults to brown.
    
    """
    index = closest_index(x, xs)
    x = xs[index]
    y = ys[index]
    if position == 'under':
        y *= 0.998
        y_text = y - dy - dy_text
    elif position == 'over':
        y *= 1.002
        y_text = y + dy + dy_text
    else:
        raise ValueError(f"position must be either 'over' or 'under', not '{position}'")
    dx = 0
    if color is None: color = default_dark_color
    color = 0.60 * color
    plt.arrow(x, y, dx, dy, linestyle='-', alpha=0.8, color=color, linewidth=1)
    plt.text(x, y_text, text, color=0.75*color, horizontalalignment='center', fontsize=12)
    

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
                                center=None, lb_min=None, ub_max=None, p=None,
                                f_max=None, f_min=None):
    get_max = lambda x: max([(f_max or np.max)(i) for i in x]) if isinstance(x, list) else ((f_max or np.max)(x))
    get_min = lambda x: min([(f_min or np.min)(i) for i in x]) if isinstance(x, list) else (f_min or np.min)(x)
    lb = min([get_min(i) for i in data])
    ub = max([get_max(i) for i in data])
    return rounted_tickmarks_from_range(lb, ub, N_ticks, step_min, lb_max, ub_min, expand, f, center,
                                        lb_min, ub_max, p)

def rounted_tickmarks_from_range(lb, ub, N_ticks, step_min=None, lb_max=None, ub_min=None,
                                 expand=None, f=None, center=None, lb_min=None, ub_max=None, p=None):
    if lb_max is not None: lb = min(lb, lb_max)
    if expand is None: expand = 0.10
    diff = expand * (ub - lb)
    ub += diff
    if ub_min is not None: ub = max(ub, ub_min)
    if ub_max is not None: ub = min(ub, ub_max)
    if lb_min is not None: lb = max(lb, lb_min)
    return rounded_linspace(lb, ub, N_ticks, step_min, f, center, p)

def rounded_linspace(lb, ub, N, step_min=None, f=None, center=None, p=None):
    if step_min is not None:
        lb = floor(lb / step_min) * step_min
        ub = ceil(ub / step_min) * step_min
    step = (ub - lb) / (N - 1)
    if p:
        x = lb % p
        if x: lb -= x
        step = (ub - lb) / (N - 1)
        x = step % p
        if x: step += p - x
        f = lambda x: x
    elif f is None:
        f = int
        step = int(ceil(step))
        lb = int(floor(lb))
    else:
        step = f(step)
        lb = f(lb)
    values = [0, 1] if step == 0 else [lb + step * i for i in range(N)]
    if center is not None:
        offset = min(values, key=lambda x: abs(center - x))
        values = [i - offset for i in values[0:-1]]
        values = [*values, values[-1] + step]
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

def contour_subplots(N_rows, N_cols, single_colorbar=False, wbar=None):
    if wbar is None: wbar = 1
    widths = np.ones(N_cols + 1)
    widths[-1] *= wbar / 4
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

# %% Plot data tables

def plot_stacked_bar(data, names, xlabels, colors=None, hatches=None, legend=True, 
                    # format_total=None, bold_label=False, fraction=False, 
                    legend_kwargs=None, ylabel=None, **kwargs):
    """Plot data table as a stacked bar chart."""
    colors, hatches = default_colors_and_hatches(data.shape[0], colors, hatches)
    N_metrics = len(xlabels)
    if isinstance(names, pd.MultiIndex): names = [i[-1] for i in names]
    if ylabel:
        ylabel, *other = ylabel.split('[')
        units = other[-1].split(']')[0]
        units = format_units(units)
        ylabel += f"[{units}]"
    df = pd.DataFrame(data, index=names, columns=xlabels)
    df.T.plot(kind='bar', stacked=True, edgecolor='k', **kwargs)
    locs, labels = plt.xticks()
    plt.xticks(locs, ['\n['.join(i.get_text().split(' [')) for i in labels])
    if legend: plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
    plt.xticks(rotation=0)
    fig = plt.gcf()
    ax = plt.gca()
    if ylabel is not None: ax.set_ylabel(ylabel)
    xticks, _ = plt.xticks()
    xlim = plt.xlim()
    y_twin = ax.twiny()
    plt.sca(y_twin)
    y_twin.tick_params(axis='x', top=True, direction="in", length=0)
    y_twin.zorder = 2
    plt.xlim(xlim)
    plt.xticks(xticks, ['' for i in xticks], va='baseline')
    N_marks = N_metrics
    axes = np.array([ax])
    if legend_kwargs is None: legend_kwargs = {}
    modify_stacked_bars(axes, N_marks, names, colors, hatches, legend, **legend_kwargs)
    return fig, axes

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

def format_spearman_plot(ax, index, name, yranges, xlabel=None):
    plot_vertical_line(0, color=c.neutral_shade.RGBn, lw=1)
    ax.xaxis.set_major_locator(MultipleLocator(0.5))
    ax.xaxis.set_major_formatter('{x:.2f}')
    ax.xaxis.set_minor_locator(MultipleLocator(0.25))
    yticks = [i[0]+i[1]/2 for i in yranges]
    ax.set_xlim(-1, 1)
    if name: 
        ax.set_xlabel(f"Spearman's correlation with {name}" if xlabel is None else xlabel)
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
        colors = c.red_tint.RGBn, c.blue_tint.RGBn
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
                     cutoff=None, xlabel=None): # pragma: no coverage
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
        format_spearman_plot(ax, index, name, yranges, xlabel)
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
        rhos_max = np.abs(values).max(axis=0)
        indices.sort(key=lambda x: rhos_max[x])
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

class Box:
    __slots__ = ('axis', 'light', 'dark', '_position', '_baseline_position')
    
    def __init__(self, axis, position=0, light=None, dark=None):
        self.axis = axis
        self.light = light
        self.dark = dark
        self._baseline_position = self._position = position
        
    def get_position(self, shift=1):
        self._position += shift
        return self._position
    
    def reset(self):
        self._position = self._baseline_position

def plot_montecarlo(data, 
                    light_color=None,
                    dark_color=None,
                    positions=None,
                    xmarks=None,
                    transpose=None,
                    vertical=True,
                    outliers=True): # pragma: no coverage
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
    if outliers: 
        flierprops = {'marker':'D',
                      'markerfacecolor': light_color,
                      'markeredgecolor': dark_color,
                      'markersize':6}
    else:
        flierprops = {'marker': ''}
    bx = plt.boxplot(x=data, positions=positions, patch_artist=True,
                     widths=0.8, whis=[5, 95], vert=vertical,
                     boxprops={'facecolor':light_color,
                               'edgecolor':dark_color},
                     medianprops={'color':dark_color,
                                  'linewidth':1.5},
                     flierprops=flierprops)
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

# %% KDE

def plot_kde(x, y, nbins=100, ax=None,
             xticks=None, yticks=None, xticklabels=None, yticklabels=None,
             xtick0=True, ytick0=True, xtickf=True, ytickf=True,
             xbox=None, ybox=None, xbox_kwargs=None, ybox_kwargs=None, 
             aspect_ratio=1.25, cmaps=None, **kwargs):
    axis_not_given = ax is None
    if axis_not_given:
        grid_kw = dict(height_ratios=[1, 8], width_ratios=[8, aspect_ratio])
        fig, all_axes = plt.subplots(
            ncols=2, nrows=2, 
            gridspec_kw=grid_kw,
        )
        ax_empty = all_axes[0, 1]
        ax = all_axes[1, 0]
        xbox_ax = all_axes[0, 0]
        ybox_ax = all_axes[1, 1]
        xbox = Box(xbox_ax, **(xbox_kwargs or {}))
        ybox = Box(ybox_ax, **(ybox_kwargs or {}))
    xs = x if isinstance(x, tuple) else (x,)
    ys = y if isinstance(y, tuple) else (y,)
    if cmaps is None: cmaps = len(xs) * [None]
    for x, y, cmap in zip(xs, ys, cmaps):
        # Evaluate a gaussian kde on a regular grid of nbins x nbins over data extents
        k = kde.gaussian_kde([x, y])
        z = k(np.vstack([x, y]))
        # Sort the points by density, so that the densest points are plotted last
        idx = z.argsort()
        try:
            x, y, z = x[idx], y[idx], z[idx]
        except:
            breakpoint()
        
        # 2D Density with shading
        plt.sca(ax)
        
        plt.scatter(x, y, c=z, s=1., cmap=cmap, **kwargs)
        if xbox:
            plt.sca(xbox.axis)
            plot_montecarlo(x, xbox.light, xbox.dark, positions=(xbox.get_position(-1),), vertical=False, outliers=False)
        if ybox:
            plt.sca(ybox.axis)
            plot_montecarlo(y, ybox.light, ybox.dark, positions=(ybox.get_position(),), vertical=True, outliers=False)
    style_axis(ax, xticks, yticks, xticklabels, yticklabels, trim_to_limits=True,
               xtick0=xtick0, ytick0=ytick0, xtickf=xtickf, ytickf=ytickf)
    plt.sca(ax)
    if axis_not_given:
        if xticks is None:
            x0, xf = plt.xlim()
        else:
            x0 = xticks[0]
            xf = xticks[-1]
        if yticks is None:
            y0, yf = plt.ylim()
        else:
            y0 = yticks[0]
            yf = yticks[-1]
        plt.sca(ax_empty); plt.axis('off')
        plt.sca(xbox.axis); plt.axis('off')
        plt.xlim([x0, xf])
        plt.sca(ybox.axis); plt.axis('off')
        plt.ylim([y0, yf])
        plt.subplots_adjust(
            hspace=0.05, wspace=0.05,
            top=0.95, bottom=0.12,
            left=0.1, right=0.96,
        )
    return ax
    
def plot_kde_2d(xs, ys, nbins=100, axes=None, xboxes=None, yboxes=None,
                xticks=None, yticks=None, xticklabels=None, yticklabels=None,
                autobox=True, xbox_kwargs=None, ybox_kwargs=None, aspect_ratio=1.,
                **kwargs):
    N_rows, N_cols, *_ = xs.shape
    if axes is None:
        if autobox:
            grid_kw = dict(height_ratios=[1/N_cols, *N_rows*[4]], width_ratios=[*N_cols*[4], aspect_ratio/N_rows])
            fig, all_axes = plt.subplots(
                ncols=N_cols + 1, nrows=N_rows + 1, 
                gridspec_kw=grid_kw,
            )
            ax_empty = all_axes[0, -1]
            axes = all_axes[1:, :-1]
            xbox_axes = all_axes[0, :-1]
            ybox_axes = all_axes[1:, -1]
            if xbox_kwargs is None: xbox_kwargs = N_cols*[{}]
            if ybox_kwargs is None: ybox_kwargs = N_rows*[{}]
            xboxes = [Box(xbox_axes[i], **xbox_kwargs[i]) for i in range(N_cols)]
            yboxes = [Box(ybox_axes[i], **ybox_kwargs[i]) for i in range(N_rows)]
            plt.sca(ax_empty)
            plt.axis('off')
            for i, box in enumerate(xboxes):
                plt.sca(box.axis)
                plt.axis('off')
                if xticks is not None: plt.xlim([xticks[i][0], xticks[i][-1]])
            for i, box in enumerate(yboxes):
                plt.sca(box.axis)
                plt.axis('off')
                if yticks is not None: plt.ylim([yticks[i][0], yticks[i][-1]])
        else:
            fig, axes = plt.subplots(ncols=N_cols, nrows=N_rows)
            axes = axes.reshape([N_rows, N_cols])
    for i in range(N_rows):
        for j in range(N_cols):
            x = xs[i, j]
            y = ys[i, j]
            ax = axes[i, j]
            xbox = None if xboxes is None else xboxes[j]
            ybox = None if yboxes is None else yboxes[i]
            xticksj = None if xticks is None else xticks[j]
            yticksi = None if yticks is None else yticks[i]
            xticklabelsj = None if xticklabels is None else xticklabels[j]
            yticklabelsi = None if yticklabels is None else yticklabels[i]
            if xticksj is not None and i != N_rows - 1:
                xticklabelsj = len(xticksj) * ['']
            if yticksi is not None and j != 0:
                yticklabelsi = len(yticksi) * ['']
            plot_kde(x, y, nbins=nbins, ax=ax, 
                     xbox=xbox,
                     ybox=ybox,
                     xticks=xticksj,
                     yticks=yticksi,
                     xticklabels=xticklabelsj,
                     yticklabels=yticklabelsi,
                     xtick0=j==0,
                     ytick0=i==N_rows-1,
                     xtickf=j==N_cols-1,
                     ytickf=i==0,
                     **kwargs)
    if xboxes: 
        for i in xboxes: i.reset()
    if yboxes: 
        for i in yboxes: i.reset()
    plt.subplots_adjust(hspace=0, wspace=0)
    return axes

    
# %% Contours

def generate_contour_data(
        z_at_xy, xlim, ylim, n=5, file=None, load=True, save=True,
        strict_convergence=None, filterwarnings=True, smooth=True, 
        vectorize=True, args=(),
    ):
    if strict_convergence is not None: 
        bst.System.strict_convergence = strict_convergence
    x0, xf = xlim
    y0, yf = ylim
    x = np.linspace(x0, xf, n)
    y = np.linspace(y0, yf, n)
    X, Y = np.meshgrid(x, y)
    if file and load:
        Z = np.load(file, allow_pickle=True)
    else:
        if filterwarnings:
            from warnings import filterwarnings
            filterwarnings('ignore')
        data0 = z_at_xy(x0, y0, *args)
        if vectorize:
            N_args = len(args)
            Z_at_XY = np.vectorize(
                z_at_xy, signature=f'(),()->{data0.shape}',
                excluded=tuple(range(2, 2 + N_args)),
            )
        else:
            Z_at_XY = z_at_xy
        Z = Z_at_XY(X, Y, *args)
        if smooth: # Smooth curves due to avoid discontinuities
            from scipy.ndimage.filters import gaussian_filter
            *_, M, N = Z.shape
            for i in range(M):
                for j in range(N):
                    Z[:, :, i, j] = gaussian_filter(Z[:, :, i, j], smooth)
    if file and save and not load: np.save(file, Z)
    return X, Y, Z

def plot_contour_2d(X, Y, Z, 
                    xlabel, ylabel, xticks, yticks, 
                    metric_bars, titles=None, 
                    fillcolor=None, styleaxiskw=None,
                    label=False, wbar=1): # pragma: no coverage
    """Create contour plots and return the figure and the axes."""
    if isinstance(metric_bars[0], MetricBar):
        nrows = len(metric_bars)
        ncols = Z.shape[-1] if titles is None else len(titles)
        row_bars = True
    else:
        nrows = len(metric_bars)
        ncols = len(metric_bars[0])
        row_bars = False
    assert Z.shape == (*X.shape, nrows, ncols), (
       f"Z was shape {Z.shape}, but expeted shape {(*X.shape, nrows, ncols)}; "
        "Z.shape must be (X, Y, M, N), where (X, Y) is the shape of both X and Y, "
        "M is the number of metrics, and N is the number of elements in titles (if given)"  
    )
    if row_bars:
        fig, axes = contour_subplots(nrows, ncols, wbar=wbar)
        cbs = np.zeros([nrows], dtype=object)
    else:
        wr = np.ones(ncols)
        wr[-1] += wbar / (3 * ncols)
        gs = dict(
            width_ratios=wr
        )
        fig, axes = plt.subplots(nrows=nrows, ncols=ncols, gridspec_kw=gs)
        cbs = np.zeros([nrows, ncols], dtype=object)
    if styleaxiskw is None: styleaxiskw = {}
    cps = np.zeros([nrows, ncols], dtype=object)
    linecolor = c.neutral_shade.RGBn
    other_axes = [[] for i in range(nrows)]
    for row in range(nrows):
        metric_row = metric_bars[row]
        for col in range(ncols):
            if row_bars:
                metric_bar = metric_row
            else:
                metric_bar = metric_row[col]
            ax = axes[row, col]
            plt.sca(ax)
            style_plot_limits(xticks, yticks)
            yticklabels = col == 0
            xticklabels = row == nrows - 1
            if fillcolor is not None: fill_plot(fillcolor)
            metric_data = Z[:, :, row, col]
            cp = plt.contourf(X, Y, metric_data,
                              levels=metric_bar.levels,
                              cmap=metric_bar.cmap)
            if label:
                cs = plt.contour(cp, zorder=1,
                                 linestyles='dashed', linewidths=0.5,
                                 levels=cp.levels, colors=[linecolor])
                clabels = ax.clabel(cs, levels=[i for i in cs.levels[::2] if i!=metric_bar.levels[-1]], inline=True, fmt=metric_bar.fmt,
                          colors=['k'], zorder=1)
                for i in clabels: i.set_rotation(0)
            cps[row, col] = cp
            if not row_bars:
                clabel = col == ncols - 1
                if clabel: 
                    pad = 0.175
                else:
                    pad = 0.05
                cbs[row, col] = metric_bar.colorbar(fig, ax, cp, shrink=metric_bar.shrink, label=clabel, pad=pad)
            other_axes[row].append(
                style_axis(ax, xticks, yticks, xticklabels, yticklabels, **styleaxiskw)
            )
        if row_bars:
            cbar_ax = axes[row, -1]
            cbs[row] = metric_bar.colorbar(fig, cbar_ax, cp, fraction=0.5, shrink=metric_bar.shrink,)
        
        # plt.clim()
    for col in range(ncols):
        title = titles[col]
        ax = axes[0, col]
        ax.set_title(title)
    if row_bars:
        for ax in axes[:, -1]:
            plt.sca(ax)
            plt.axis('off')
        set_axes_labels(axes[:, :-1], xlabel, ylabel)
    else:
        set_axes_labels(axes, xlabel, ylabel)
    plt.subplots_adjust(hspace=0.1, wspace=0.1)
    return fig, axes, cps, cbs, other_axes
       
def plot_contour_single_metric(
        X, Y, Z, xlabel, ylabel, xticks, yticks, metric_bar,
        titles=None, fillcolor=None, styleaxiskw=None, label=False
    ): # pragma: no coverage
    """Create contour plots and return the figure and the axes."""
    *_, nrows, ncols = Z.shape
    assert Z.shape == (*X.shape, nrows, ncols), (
        "Z.shape must be (X, Y, M, N), where (X, Y) is the shape of both X and Y"
    )
    fig, axes, ax_colorbar = contour_subplots(nrows, ncols, single_colorbar=True)
    if styleaxiskw is None: styleaxiskw = {}
    cps = np.zeros([nrows, ncols], dtype=object)
    linecolor = c.neutral_shade.RGBn
    other_axes = []
    for row in range(nrows):
        for col in range(ncols):
            ax = axes[row, col]
            plt.sca(ax)
            style_plot_limits(xticks, yticks)
            yticklabels = col == 0
            xticklabels = row == nrows - 1
            if fillcolor is not None: fill_plot(fillcolor)
            metric_data = Z[:, :, row, col]
            cp = plt.contourf(X, Y, metric_data,
                              levels=metric_bar.levels,
                              cmap=metric_bar.cmap,
                              norm=metric_bar.norm)
            for i in cp.collections:
                i.set_edgecolor('face') # For svg background
            if label:
                cs = plt.contour(cp, zorder=1,
                                 linestyles='dashed', linewidths=1.,
                                 norm=metric_bar.norm,
                                 levels=metric_bar.levels, colors=[linecolor])
                clabels = ax.clabel(
                    cs, levels=[i for i in cs.levels if i!=metric_bar.levels[-1]], inline=True, fmt=metric_bar.fmt,
                    colors=['k'], zorder=1
                )
                for i in clabels: i.set_rotation(0)
            cps[row, col] = cp
            dct = style_axis(ax, xticks, yticks, xticklabels, yticklabels, **styleaxiskw)
            other_axes.append(dct)
    cb = metric_bar.colorbar(fig, ax_colorbar, cp, fraction=0.5)
    plt.sca(ax_colorbar)
    plt.axis('off')
    if titles:
        for col, title in enumerate(titles):
            ax = axes[0, col]
            ax.set_title(title)
    set_axes_labels(axes[:, :-1], xlabel, ylabel)
    plt.subplots_adjust(hspace=0.1, wspace=0.1)
    return fig, axes, cps, cb, other_axes
            
def color_quadrants(color=None, x=None, y=None, xlim=None, ylim=None, 
                    line_color=None, linewidth=1.0):
    if x is None: x = 0
    if y is None: y = 0
    if line_color is None: line_color = c.grey.RGBn
    if xlim is None: xlim = plt.xlim()
    if ylim is None: ylim = plt.ylim()
    x0, x1 = xlim
    y0, y1 = ylim
    top_left, top_right, bottom_left, bottom_right = color
    # Top left
    if top_left is not None:
        plt.fill_between([x0, x], y, y1,
                         color=top_left,
                         linewidth=linewidth,
                         zorder=0)
    # Top right
    if top_right is not None:
        plt.fill_between([x, x1], y, y1,
                         color=top_right,
                         linewidth=linewidth,
                         zorder=0)
    # Bottom left
    if bottom_left is not None:
        plt.fill_between([x0, x], y0, y,
                         color=bottom_left,
                         linewidth=linewidth,
                         zorder=0)
    # Bottom right
    if bottom_right is not None:
        plt.fill_between([x, x1], y0, y,
                         color=bottom_right,
                         linewidth=linewidth,
                         zorder=0)
    plot_vertical_line(x, line_color, zorder=0)
    plot_horizontal_line(y, line_color, zorder=0)

def label_quadrants(
        x=None, y=None, text=None, color=None,
    ):
    xlb, xub = plt.xlim()
    ylb, yub = plt.ylim()
    data_given = not (x is None or y is None)
    if data_given:
        y_mt_0 = y > 0
        y_lt_0 = y < 0
        x_mt_0 = x > 0
        x_lt_0 = x < 0
    xpos = lambda x: xlb + (xub - xlb) * x
    ypos = lambda y: ylb + (yub - ylb) * y
    xleft = 0.02
    xright = 0.98
    ytop = 0.94
    ybottom = 0.02
    labeled = 4 * [False]
    top_left, top_right, bottom_left, bottom_right = text
    top_left_color, top_right_color, bottom_left_color, bottom_right_color = color
    if yub > 0. and xlb < 0. and top_left:
        if data_given and top_left.endswith('()'):
            p = (y_mt_0 & x_lt_0).sum() / y.size
            top_left = f"{p:.0%} {top_left.strip('()')}"
        plt.text(xpos(xleft), ypos(ytop), top_left, color=top_left_color,
                 horizontalalignment='left', verticalalignment='top',
                 fontsize=10, fontweight='bold', zorder=10)
        labeled[0] = True
    if yub > 0. and xub > 0. and top_right:
        if data_given and top_right.endswith('()'):
            p = (y_mt_0 & x_mt_0).sum() / y.size
            top_right = f"{p:.0%} {top_right.strip('()')}"
        plt.text(xpos(xright), ypos(ytop), top_right, color=top_right_color,
                 horizontalalignment='right', verticalalignment='top',
                 fontsize=10, fontweight='bold', zorder=10)
        labeled[1] = True
    if ylb < 0. and xlb < 0. and bottom_left:
        if data_given and bottom_left.endswith('()'):
            p = (y_lt_0 & x_lt_0).sum() / y.size
            bottom_left = f"{p:.0%} {bottom_left.strip('()')}"
        plt.text(xpos(xleft), ypos(ybottom), bottom_left, color=bottom_left_color,
                 horizontalalignment='left', verticalalignment='bottom',
                 fontsize=10, fontweight='bold', zorder=10)
        labeled[2] = True
    if ylb < 0. and xub > 0. and bottom_right:
        if data_given and bottom_right.endswith('()'):
            p = (y_lt_0 & x_mt_0).sum() / y.size
            bottom_right = f"{p:.0%} {bottom_right.strip('()')}"
        plt.text(xpos(xright), ypos(ybottom), bottom_right, color=bottom_right_color,
                 horizontalalignment='right', verticalalignment='bottom',
                 fontsize=10, fontweight='bold', zorder=10)
        labeled[3] = True
    return labeled

def plot_quadrants(
        text, data=None, x=None, y=None, rotate=0,
    ):
    quadrant_color = deque([
        (*c.CABBI_teal.tint(90).RGBn, 0.9), None,
        None, (*c.red.tint(80).RGBn, 0.9)
    ])
    text_color = deque([
        CABBI_colors.teal.shade(50).RGBn, CABBI_colors.grey.shade(75).RGBn, 
        CABBI_colors.grey.shade(75).RGBn, c.red.shade(50).RGBn
    ])
    quadrant_color.rotate(rotate)
    text_color.rotate(rotate)
    labeled_quadrants = format_quadrants(
        data, x, y, text, text_color, quadrant_color,
    )
    for i, labeled in enumerate(labeled_quadrants):
        if labeled: text[i] = '()' if text[i].endswith('()') else None
    
    return labeled_quadrants

def format_quadrants(
        data=None,
        x=None, y=None, 
        text=None,
        text_color=None,
        quadrant_color=None,
        xlim=None, ylim=None,
    ):
    if data is None: data = (None, None) 
    color_quadrants(quadrant_color, x, y, xlim, ylim)
    return label_quadrants(
        *data, text, text_color,
    )
    
def add_titles(axes, titles, color):
    if titles:
        plt.subplots_adjust(
            top=0.90,
        )
        for ax, letter in zip(axes, titles):
            plt.sca(ax)
            ylb, yub = plt.ylim()
            xlb, xub = plt.xlim()
            plt.text((xlb + xub) * 0.5, ylb + (yub - ylb) * 1.17, letter, color=color,
                      horizontalalignment='center', verticalalignment='center',
                      fontsize=12, fontweight='bold')
