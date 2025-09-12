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
from biosteam.utils import colors as c, CABBI_colors, GG_colors
from colorpalette import Color, ColorWheel
from .utils import style_axis, style_plot_limits, fill_plot, set_axes_labels, MetricBar, closest_index
from thermosteam.units_of_measure import format_units, reformat_units
import matplotlib.patches as mpatches
from math import floor, ceil
from matplotlib.ticker import MultipleLocator
from scipy.stats.kde import gaussian_kde
from collections import deque
from itertools import product
from matplotlib import colormaps
import matplotlib.colors as clr

__all__ = (
    'rounded_linspace',
    'rounted_tickmarks_from_range',
    'rounded_tickmarks_from_data',
    'annotate_point',
    'annotate_line',
    'plot_unit_groups',
    'plot_unit_groups_across_coordinate',
    'plot_uncertainty_boxes',
    'plot_montecarlo', 
    'plot_uncertainty_across_coordinate',
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
    'plot_contour',
    'plot_contour_2d', 
    'plot_contour_single_metric',
    'plot_heatmap',
    'plot_kde_2d',
    'plot_kde_1d',
    'plot_kde',
    'plot_quadrants',
    'plot_stacked_bar',
    'generate_contour_data',
    'default_colors_and_hatches',
    'modify_stacked_bars',
    'title_color',
    'contour_subplots',
)

# %% Utilities

plt.rcParams['figure.dpi'] = 300 # High DPI (default is 100; so low!)
default_light_color = c.orange_tint.RGBn
default_dark_color = c.orange_shade.RGBn
title_color = c.neutral.shade(25).RGBn
default_color_wheel = GG_colors.wheel(keys=['red', 'blue', 'purple', 'orange', 'green'])

def annotate_point(
        text, x, y, dx=0, dy=0.2, dx_text=0, dy_text=0.22,
        textcolor=None, linecolor=None, fontsize=12, 
        horizontalalignment='center', 
        verticalalignment='bottom',
        arrowkwargs=None, 
        textkwargs=None,
        xlinestyle='-', 
        ylinestyle='-',
    ): # pragma: no coverage
    """
    Annotate point with text and arrow pointing to text.
    
    Parameters
    ----------
    text : str
    x : float
        Arrow position
    dx : float
        Length of arrow to x-position.
    dy : float
        Length of arrow to y-position.
    dy_text : float
        Distance of text to arrow.
    color : numpy.ndarray
        RGB normalized to 1. Defaults to brown.
    
    """
    xtext = x + dx + dx_text
    ytext = y + dy + dy_text
    if textcolor is None: textcolor = 0.45 * default_dark_color
    if linecolor is None: linecolor = 0.60 * default_dark_color
    if arrowkwargs is None: arrowkwargs = {}
    if textkwargs is None: textkwargs = {}
    plt.arrow(x, y, dx, dy, linestyle=xlinestyle, alpha=0.8, color=linecolor, 
              **arrowkwargs)
    plt.text(xtext, ytext, text, color=textcolor, 
             horizontalalignment=horizontalalignment,
             verticalalignment=verticalalignment,
             fontsize=fontsize, **textkwargs)

def annotate_line(text, x, xs, ys, *args, **kwargs): # pragma: no coverage
    """
    Annotate line with text and arrow pointing to text.
    
    Parameters
    ----------
    text : str
    x : float
        Arrow position
    xs : numpy.ndarray(dim=1)
    ys : numpy.ndarray(dim=1)
    dx : float
        Length of arrow to x-position.
    dy : float
        Length of arrow to y-position.
    dy_text : float
        Distance of text to arrow.
    color : numpy.ndarray
        RGB normalized to 1. Defaults to brown.
    
    """
    index = closest_index(x, xs)
    x = xs[index]
    y = ys[index]
    annotate_point(text, x, y, **kwargs)

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
    return values
        
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
                          r"$\mathbf{[" f"{format_units(i.units, '', False)}" "]}$"
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
        if lb < 0: plot_horizontal_line(0, color=CABBI_colors.black.RGBn, linestyle='--')
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

def plot_spearman(rhos, top=None, name=None, color_wheel=None, index=None, **kwargs):
    if isinstance(rhos, list): 
        return plot_spearman_2d(rhos, top, name, color_wheel=color_wheel, index=index, **kwargs)
    else:
        return plot_spearman_1d(rhos, top, name, color=color_wheel, index=index, **kwargs)

def format_spearman_plot(ax, index, name, yranges, xlabel=None):
    plot_vertical_line(0, color=c.neutral_shade.RGBn, lw=1)
    ax.xaxis.set_major_locator(MultipleLocator(0.5))
    ax.xaxis.set_major_formatter('{x:.2f}')
    ax.xaxis.set_minor_locator(MultipleLocator(0.25))
    yticks = [i[0]+i[1]/2 for i in yranges]
    ax.set_xlim(-1, 1)
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    elif name is not None:
        ax.set_xlabel(f"Spearman's correlation with {name}")
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
        top=None, name=None, colors=None, w=0.8, s=1., offset=0., style=True, 
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
    mask = np.isnan(lb) | np.isnan(ub)
    lb[mask] = baseline
    ub[mask] = baseline
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
                     w=None, s=1., offset=0., style=True, 
                     fig=None, ax=None, sort=True, index=None,
                     cutoff=None, xlabel=None, edgecolors=None): # pragma: no coverage
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
    if edgecolors is None: edgecolors = 'k'
    if w is None: w = 0.8
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
                       edgecolors=edgecolors)
    
    if style:
        if index is None:
            raise ValueError('must pass index if rhos is not a pandas Series object')
        format_spearman_plot(ax, index, name, yranges, xlabel)
    return fig, ax

def plot_spearman_2d(rhos, top=None, name=None, color_wheel=None, index=None,
                     cutoff=None, sort=True, xlabel=None, sort_index=None, w=None,
                     edgecolors=None): # pragma: no coverage
    """
    Display Spearman's rank correlation plot.
    
    Parameters
    ----------
    rhos : list[pandas.Series]|array
         Spearman's rank correlation coefficients to be plotted.
    top=None : float, optional
        Number of parameters to plot (from highest values).
    color_wheel=None: Iterable, optional
        Iterable of colorpalette.Color objects or RGBn values.
    
    Returns
    -------
    fig : matplotlib Figure
    ax : matplotlib AxesSubplot
    """
    rhos = list(reversed(rhos))
    if name is None: name = rhos[0].name
    if index is None: index = rhos[0].index
    values = np.array([i.values if hasattr(i, 'values') else i for i in rhos])
    indices = list(range(values.shape[-1]))
    if cutoff:
        cutoff_index, = np.where(np.any(np.abs(rhos) > cutoff, axis=0))
        indices = [indices[i] for i in cutoff_index]
    if sort:
        if sort_index is None:
            rhos_max = np.abs(values).max(axis=0)
        else:
            rhos_max = np.abs(values)[sort_index]
        indices.sort(key=lambda x: rhos_max[x])
    if top is not None: indices = indices[-top:]
    rhos = [[rho[i] for i in indices] for rho in values]
    index = [index[i] for i in indices]
    N = len(rhos)
    s = N + 1
    if not color_wheel: color_wheel = CABBI_colors.wheel()
    fig, ax = plt.subplots()
    for i, rho in enumerate(rhos):
        if edgecolors:
            edgecolor = edgecolors[N - i - 1]
            if hasattr(edgecolor, 'RGBn'): edgecolor = edgecolor.RGBn
        else:
            edgecolor = None
        color = color_wheel[N - i - 1]
        if hasattr(color, 'RGBn'): 
            if edgecolor is None:
                edgecolor = color.shade(50).RGBn
            color = color.RGBn
        plot_spearman_1d(rho, color=color, s=s, offset=i, w=w, edgecolors=edgecolor,
                         fig=fig, ax=ax, style=False, sort=False, top=None)
    # Plot central line
    yranges = [(s/2 + s*i - 1., 1.) for i in range(len(rhos[0]))]
    format_spearman_plot(ax, index, name, yranges, xlabel)
    return fig, ax

# %% Monte Carlo

class Box:
    __slots__ = ('axis', 'fill', 'edge', '_position', '_baseline_position')
    
    def __init__(self, axis, position=None, fill=None, edge=None, light=None, dark=None):
        if light is not None: fill = light # For backwards compatibility
        if dark is not None: edge = dark # For backwards compatibility
        if position is None: position = [0]
        self.axis = axis
        self.fill = fill
        self.edge = edge
        self._baseline_position = position[0]
        self._position = position
        
    @property # For backwards compatibility
    def light(self): return self.fill
    @property # For backwards compatibility
    def dark(self): return self.edge
    
    def get_position(self, shift=1):
        self._position[0] += shift
        return self._position[0]
    
    def reset(self):
        self._position = self._baseline_position

def plot_uncertainty_boxes(
        data, 
        fill_color=None,
        edge_color=None,
        positions=None,
        xmarks=None,
        transpose=None,
        vertical=True,
        outliers=False,
        bounds=True,
        width=None,
        hatch=None,
        light_color=None,
        dark_color=None,
    ): # pragma: no coverage
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
    if transpose is None and hasattr(data, 'ndim') and data.ndim == 2:
        N_rows, N_cols = data.shape
        if N_cols < N_rows: data = data.transpose()
    elif transpose:
        data = data.transpose()
    if hasattr(data, 'ndim') and data.ndim == 1:
        data = [data]
    if not positions:
        if hasattr(data, 'ndim'):
            if data.ndim != 2: 
                positions = (0,)
            else:
                positions = list(range(data.shape[0]))
        else:
            positions = list(range(len(data)))
    N_positions = len(positions)
    if width is None: width = 0.8
    if light_color is not None: fill_color = light_color
    if dark_color is not None: edge_color = dark_color
    if fill_color is None: fill_color = default_light_color
    if edge_color is None: edge_color = default_dark_color
    if isinstance(fill_color[0], float):
        fill_color = N_positions * [fill_color]
    if isinstance(edge_color[0], float):
        edge_color = N_positions * [edge_color]
    boxes = []
    for x, pos, fill, edge in zip(data, positions, fill_color, edge_color):
        if outliers: 
            flierprops = {'marker':'D',
                          'markerfacecolor': fill,
                          'markeredgecolor': edge,
                          'markersize':6}
        else:
            flierprops = {'marker': ''}
        bx = plt.boxplot(x=[x], positions=[pos], patch_artist=True,
                         widths=width, whis=[5, 95], vert=vertical,
                         boxprops={'facecolor':fill,
                                   'edgecolor':edge},
                         medianprops={'color':edge,
                                      'linewidth':1.5},
                         capprops=dict(color=edge),
                         whiskerprops=dict(color=edge),
                         manage_ticks=False,
                         flierprops=flierprops)
        boxes.extend(bx['boxes'])
    # bx = plt.boxplot(x=data, positions=positions, patch_artist=True,
    #                  widths=width, whis=[5, 95], orientation='vertical' if vertical else 'horizontal')
    # boxprops = {'facecolor':fill_color,
    #           'edgecolor':edge_color}
    # medianprops = {'color':edge_color,
    #              'linewidth':1.5}
    # capprops = dict(color=edge_color)
    # whiskerprops = dict(color=edge_color)
    # whiskers = bx['whiskers']
    # flierprops = {'marker':'D',
    #               'markerfacecolor': fill_color,
    #               'markeredgecolor': edge_color,
    #               'markersize':6} if outliers else {'marker': ''}
    # caps = bx['caps']
    # boxes = bx['boxes']
    # medians = bx['medians']
    # fliers = bx['fliers']
    # for i in range(N_positions):
    #     for name, value in boxprops.items():
    #         getattr(boxes[i], 'set_' + name)(value if isinstance(value, (str, float)) else value[i])
    #     for name, value in medianprops.items():
    #         getattr(medians[i], 'set_' + name)(value if isinstance(value, (str, float)) else value[i])
    #     for name, value in capprops.items():
    #         getattr(caps[i], 'set_' + name)(value if isinstance(value, (str, float)) else value[i])
    #     for name, value in flierprops.items():
    #         getattr(fliers[i], 'set_' + name)(value if isinstance(value, (str, float)) else value[i])
    #     for name, value in whiskerprops.items():
    #         getattr(whiskers[i], 'set_' + name)(value if isinstance(value, (str, float)) else value[i])
    # boxes = bx['boxes']
    if bounds:
        if vertical:
            plt.scatter(x=positions, y=[i.min() for i in data], marker='1', c=edge_color)
            plt.scatter(x=positions, y=[i.max() for i in data], marker='2', c=edge_color)
        else:
            plt.scatter(x=[i.min() for i in data], y=positions, marker='1', c=edge_color)
            plt.scatter(x=[i.max() for i in data], y=positions, marker='2', c=edge_color)
    if xmarks: plt.xticks(positions, xmarks)
    if hatch:
        for box in boxes: box.set(hatch = hatch)
    return boxes
    # return bx

plot_montecarlo = plot_uncertainty_boxes # Backwards compatibility

def plot_uncertainty_across_coordinate(
        xs, ys, 
        p5_color=None,
        fill_color=None,
        median_color=None,
        smooth=0
    ): # pragma: no coverage
    """
    Plot Monte Carlo evaluation across a coordinate.
    
    Parameters
    ----------
    xs : numpy.ndarray(ndim=1)
        Coordinate values for each column in ``ys``.
    ys : numpy.ndarray(ndim=2)
        Metric values with uncertainty. Each row represents a sample and each 
        column represent a metric along the x-coordinate.
    p5_color : numpy.ndarray
        RGB normalized to 1. Defaults to brown.
    fill_color : numpy.ndarray
        RGB normalized to 1. Defaults to brown.
    median_color : numpy.ndarray
        RGB normalized to 1. Defaults to brown.
    
    Returns
    -------
    percentiles : numpy.ndarray(ndim=2)
        5, 25, 50, 75 and 95th percentiles by row (5 rows total).
    
    """
    if fill_color is None: fill_color = default_light_color
    if median_color is None: median_color = default_dark_color
    if p5_color is None: p5_color = 0.5 * (default_light_color + default_dark_color)
    q05, q25, q50, q75, q95 = percentiles = np.percentile(ys, [5,25,50,75,95], axis=0)

    if smooth:
        from scipy.ndimage.filters import gaussian_filter
        for i, p in enumerate(percentiles):
            percentiles[i] = gaussian_filter(p, smooth)

    plt.plot(xs, q50, '-',
             color=median_color,
             linewidth=1.5) # Median
    plt.fill_between(xs, q25, q75,
                     color=fill_color,
                     linewidth=1.0)
    plt.plot(xs, q05, '-.',
             color=p5_color,
             linewidth=1.0) # Lower whisker
    plt.plot(xs, q95, '-.',
             color=p5_color,
             linewidth=1.0) # Upper whisker
    
    return percentiles

plot_montecarlo_across_coordinate = plot_uncertainty_across_coordinate

# %% KDE

def plot_kde(
        x, y, axes=None, fig=None,
        xticks=None, yticks=None, xticklabels=None, yticklabels=None,
        xtick0=True, ytick0=True, xtickf=True, ytickf=True,
        xbox=None, ybox=None, xbox_kwargs=None, ybox_kwargs=None, 
        aspect_ratio=1.25, colors=None, xbox_width=None,
        ybox_width=None, zorders=None, xlabel=None, ylabel=None,
        ax=None, kde= None, transparency=None, color_wheel=None, **kwargs
    ):
    if kde is None: kde = True
    if color_wheel is None: color_wheel = default_color_wheel
    axis_not_given = axes is None and ax is None
    xs = x if isinstance(x, (tuple, list)) or x.ndim == 2 else (x,)
    ys = y if isinstance(y, (tuple, list)) or y.ndim == 2 else (y,)
    N_xs = len(xs)
    N_ys = len(ys)
    if axis_not_given:
        # grid_kw = dict(height_ratios=[0.4 * N_internal_x, 4], width_ratios=[4, 0.4 * N_internal_y])
        grid_kw = dict(height_ratios=[1, 8], width_ratios=[8, aspect_ratio])
        fig, axes = plt.subplots(
            ncols=2, nrows=2, 
            gridspec_kw=grid_kw,
        )
        ax_empty = axes[0, 1]
        ax = axes[1, 0]
        xbox_ax = axes[0, 0]
        ybox_ax = axes[1, 1]
    elif ax:
        ax_empty = xbox_ax = ybox_ax = None
    else:
        ax_empty = axes[0, 1]
        ax = axes[1, 0]
        xbox_ax = axes[0, 0]
        ybox_ax = axes[1, 1]
    if xbox is None: 
        if xbox_kwargs is None: xbox_kwargs = {}
        position = [0]
        if isinstance(xbox_kwargs, dict):
            xboxes = [Box(xbox_ax, position, **xbox_kwargs)] * N_xs
        else:
            xboxes = [Box(xbox_ax, position, **i) for i in xbox_kwargs]
    elif isinstance(xbox, Box):
        xboxes = [xbox] * N_xs
    else:
        xboxes = xbox
    if ybox is None:
        if ybox_kwargs is None: ybox_kwargs = {}
        position = [0]
        if isinstance(ybox_kwargs, dict):
            yboxes = [Box(ybox_ax, position, **ybox_kwargs)] * N_xs
        else:
            yboxes = [Box(ybox_ax, position, **i) for i in ybox_kwargs]
    elif isinstance(ybox, Box):
        yboxes = [ybox] * N_ys
    else:
        yboxes = ybox
    if colors is None: 
        if kde:
            if N_xs == 1:
                colors = [colormaps['viridis']]
            else:
                colors = [
                    clr.LinearSegmentedColormap.from_list(
                        (c:=color_wheel[i]).ID,
                        [*[c.shade(60 - 20 * j).RGBn for j in range(3)],
                         *[c.tint(20 * j).RGBn for j in range(3)]],
                        N=256
                    )
                    for i in range(N_xs)
                ]
        else:
            if transparency is None: transparency = 1
            colors = [color_wheel[i].RGBn for i in range(len(xs))]
    if zorders is None: zorders = len(xs) * [5]
    for x, y, color, zorder, xbox, ybox in zip(xs, ys, colors, zorders, xboxes, yboxes):
        plt.sca(ax)
        scatter_kwargs = kwargs.copy()
        if kde:
            # Evaluate a gaussian kde on a regular grid of nbins x nbins over data extents
            k = gaussian_kde([x, y])
            z = k(np.vstack([x, y]))
            
            # Sort the points by density, so that the densest points are plotted last
            idx = z.argsort()
            x, y, z = x[idx], y[idx], z[idx]
            
            # 2D Density with shading
            scatter_kwargs['cmap'] = color
            scatter_kwargs['c'] = z
        else:
            scatter_kwargs['c'] = color
        if 's' not in scatter_kwargs and 'size' not in scatter_kwargs:
            scatter_kwargs['s'] = 1.
        plt.scatter(x, y, zorder=zorder, **scatter_kwargs)
        if xbox:
            plt.sca(xbox.axis)
            if xbox.fill is None and xbox.edge is None:
                if kde:
                    fill = xbox.fill
                    edge = xbox.edge
                else:
                    fill = color
                    edge = Color(fg=color).shade(60).RGBn
            else:
                fill = xbox.fill
                edge = xbox.edge
            plot_uncertainty_boxes(x, fill, edge, positions=(xbox.get_position(-1),), vertical=False, outliers=False, width=xbox_width, bounds=True)
        if ybox:
            plt.sca(ybox.axis)
            if ybox.fill is None and ybox.edge is None:
                if kde:
                    fill = ybox.fill
                    edge = ybox.edge
                else:
                    fill = color
                    edge = Color(fg=color).shade(60).RGBn
            else:
                fill = ybox.fill
                edge = ybox.edge
            plot_uncertainty_boxes(y, fill, edge, positions=(ybox.get_position(-1),), vertical=True, outliers=False, width=ybox_width, bounds=True)
    # if xboxes:
    #     plt.sca(xbox_ax)
    #     fill_colors = []
    #     edge_colors = []
    #     for xbox, color in zip(xboxes, colors):
    #         if xbox.fill is None and xbox.edge is None:
    #             if kde:
    #                 fill = xbox.fill
    #                 edge = xbox.edge
    #             else:
    #                 fill = color
    #                 edge = Color(fg=color).shade(60).RGBn
    #         fill_colors.append(fill)
    #         edge_colors.append(edge)
    #     plot_uncertainty_boxes(
    #         xs, fill_colors, edge_colors, 
    #         positions=[*range(0, -len(xs), -1)],
    #         vertical=False, outliers=False, 
    #         width=xbox_width, bounds=True
    #     )
    # if yboxes:
    #     fill_colors = []
    #     edge_colors = []
    #     for ybox, color in zip(yboxes, colors):
    #         if ybox.fill is None and ybox.edge is None:
    #             if kde:
    #                 fill = ybox.fill
    #                 edge = ybox.edge
    #             else:
    #                 fill = color
    #                 edge = Color(fg=color).shade(60).RGBn
    #         fill_colors.append(fill)
    #         edge_colors.append(edge)
    #     plt.sca(ybox_ax)
    #     plot_uncertainty_boxes(
    #         ys, fill_colors, edge_colors, 
    #         positions=[*range(len(ys))],
    #         vertical=True, outliers=False, 
    #         width=ybox_width, bounds=True
    #     )
    style_axis(ax, xticks, yticks, xticklabels, yticklabels, trim_to_limits=True,
               xtick0=xtick0, ytick0=ytick0, xtickf=xtickf, ytickf=ytickf)
    plt.sca(ax)
    if xlabel is not None: plt.xlabel(xlabel)
    if ylabel is not None: plt.ylabel(ylabel)
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
        if xbox:
            plt.sca(xbox.axis); plt.axis('off')
            plt.xlim([x0, xf])
        if ybox:
            plt.sca(ybox.axis); plt.axis('off')
            plt.ylim([y0, yf])
        plt.subplots_adjust(
            hspace=0.05, wspace=0.05,
            top=0.95, bottom=0.12,
            left=0.1, right=0.96,
        )
    return fig, ax, axes
  
def plot_kde_1d(
        xs, ys, axes=None, xboxes=None, yboxes=None,
        xticks=None, yticks=None, xticklabels=None, yticklabels=None,
        autobox=True, xbox_kwargs=None, ybox_kwargs=None, aspect_ratio=1.,
        xlabel=None, ylabel=None, fs=None, colors=None, kde=None, color_wheel=None, **kwargs
    ):
    if kde is None: kde = True
    N_cols = len(xs)
    N_rows = 1
    if xticklabels and not hasattr(xticklabels, '__len__'):
        xticklabels = N_cols * [xticklabels]
    if xticks is not None and len(xticks) != N_cols:
        xticks = N_cols * [xticks]
    if color_wheel is None: color_wheel = default_color_wheel
    if axes is None:
        if autobox:
            N_internal_x = sum([(1 if hasattr(x, 'ndim') and x.ndim == 1 else len(x))
                                for x in xs]) / N_cols
            N_internal_y = sum([(1 if hasattr(y, 'ndim') and y.ndim == 1 else len(y))
                                for y in ys])
            grid_kw = dict(height_ratios=[0.3 * N_internal_x, *N_rows*[4]], width_ratios=[*N_cols*[4], 0.3 * N_internal_y])
            fig, all_axes = plt.subplots(
                ncols=N_cols + 1, nrows=N_rows + 1, 
                gridspec_kw=grid_kw,
            )
            ax_empty = all_axes[0, -1]
            axes = all_axes[1:, :-1]
            xbox_axes = all_axes[0, :-1]
            ybox_axis = all_axes[1, -1]
            if xbox_kwargs is None: 
                xbox_kwargs = N_cols*[{}]
            elif isinstance(xbox_kwargs, dict):
                xbox_kwargs = N_cols*[xbox_kwargs]
            elif len(xbox_kwargs) == 1: 
                xbox_kwargs = N_cols * xbox_kwargs
            if ybox_kwargs is None:
                ybox_kwargs = {}
            xboxes = [Box(xbox_axes[i], **xbox_kwargs[i]) for i in range(N_cols)]
            ybox = Box(ybox_axis, **ybox_kwargs)
            plt.sca(ax_empty)
            plt.axis('off')
            for i, box in enumerate(xboxes):
                plt.sca(box.axis)
                plt.axis('off')
                if xticks is not None: 
                    try: plt.xlim([xticks[i][0], xticks[i][-1]])
                    except: pass
            plt.sca(ybox.axis)
            plt.axis('off')
            if yticks is not None: 
                try: plt.ylim([yticks[0], yticks[-1]])
                except: pass
        else:
            fig, axes = plt.subplots(ncols=N_cols, nrows=N_rows)
            axes = axes.reshape([N_rows, N_cols])
    for i in range(N_cols):
        x = xs[i]
        y = ys[i]
        ax = axes[0, i]
        xticksi = None if xticks is None else xticks[i]
        xticklabelsi = None if xticklabels is None else xticklabels[i]
        plot_kde(x, y, ax=ax, fig=fig,
                 xbox=[False] if hasattr(x, 'ndim') and x.ndim == 1 else len(x) * [False],
                 ybox=[False] if hasattr(y, 'ndim') and y.ndim == 1 else len(y) * [False],
                 xticks=xticksi,
                 yticks=yticks,
                 xticklabels=xticklabelsi,
                 yticklabels=yticklabels if i == 0 else False,
                 xtick0=True,
                 ytick0=True,
                 xtickf=i==N_cols-1,
                 ytickf=True,
                 colors=colors,
                 kde=kde,
                 color_wheel=color_wheel,
                 **kwargs)
    if xboxes is not None:
        for i, xbox in enumerate(xboxes):
            x = xs[i]
            plt.sca(xbox.axis)
            x = x if isinstance(x, (tuple, list)) or x.ndim == 2 else (x,)
            if xbox.fill is None and xbox.edge is None:
                fill = []
                edge = []
                for i in range(len(x)):
                    color = color_wheel[i]
                    fill.append(color.RGBn)
                    edge.append(color.shade(60).RGBn)
            else:
                fill = xbox.fill
                edge = xbox.edge
            plot_uncertainty_boxes(
                x, fill, edge,
                positions=[*range(0, -len(x), -1)],
                vertical=False, outliers=False, 
                bounds=True
            )
        for i in range(N_cols):
            ax = axes[0, i]
            plt.sca(ax)
            lim = plt.xlim()
            plt.sca(xbox.axis)
            plt.xlim(lim)
    if ybox is not None:
        y = []
        fill = []
        edge = []
        for values in ys:
            if hasattr(values, 'ndim') and values.ndim == 1:
                y.append(values)
                fill.append(xbox.fill)
                edge.append(xbox.edge)
            else:
                y.extend(values)
                if ybox.fill is None and ybox.edge is None:
                    for i in range(len(values)):
                        color = color_wheel[i]
                        fill.append(color.RGBn)
                        edge.append(color.shade(60).RGBn)
                else:
                    fill.append(xbox.fill)
                    edge.append(xbox.edge)
        plt.sca(ybox.axis)
        plot_uncertainty_boxes(
            y, fill, edge,
            positions=[*range(len(y))],
            vertical=True, outliers=False, 
            bounds=True
        )
        ax = axes[0, 0]
        plt.sca(ax)
        ylim = plt.ylim()
        plt.sca(ybox.axis)
        plt.ylim(ylim)
    if xlabel is not None: fig.supxlabel(xlabel, fontsize=fs)
    plt.subplots_adjust(hspace=0, wspace=0)
    plt.sca(axes[0, 0])
    if ylabel is not None: plt.ylabel(ylabel, fontsize=fs)
    return fig, axes
  
def plot_kde_2d(
        xs, ys, axes=None, xboxes=None, yboxes=None,
        xticks=None, yticks=None, xticklabels=None, yticklabels=None,
        autobox=True, xbox_kwargs=None, ybox_kwargs=None, aspect_ratio=1.,
        xlabel=None, ylabel=None, fs=None, **kwargs
    ):
    xs = np.asarray(xs)
    ys = np.asarray(ys)
    if xs.ndim == 2:
        N, M = xs.shape
        xs = xs.reshape([1, N, M])
    if ys.ndim == 2:
        N, M = ys.shape
        ys = ys.reshape([1, N, M])
    N_rows, N_cols, *_ = xs.shape
    if xticklabels and not hasattr(xticklabels, '__len__'):
        xticklabels = N_cols * [xticklabels]
    if yticklabels and not hasattr(yticklabels, '__len__'):
        yticklabels = N_rows * [yticklabels]
    if xticks is not None and len(xticks) != N_cols:
        xticks = N_cols * [xticks]
    if yticks is not None and len(yticks) != N_rows:
        yticks = N_rows * [yticks]
    if axes is None:
        if autobox:
            grid_kw = dict(height_ratios=[0.5 * N_rows, *N_rows*[4]], width_ratios=[*N_cols*[4], 0.5 * N_cols])
            fig, all_axes = plt.subplots(
                ncols=N_cols + 1, nrows=N_rows + 1, 
                gridspec_kw=grid_kw,
            )
            ax_empty = all_axes[0, -1]
            axes = all_axes[1:, :-1]
            xbox_axes = all_axes[0, :-1]
            ybox_axes = all_axes[1:, -1]
            if xbox_kwargs is None: 
                xbox_kwargs = N_cols*[{}]
            elif isinstance(xbox_kwargs, dict):
                xbox_kwargs = N_cols*[xbox_kwargs]
            elif len(xbox_kwargs) == 1: 
                xbox_kwargs = N_cols * xbox_kwargs
            if ybox_kwargs is None:
                ybox_kwargs = N_rows*[{}]
            elif isinstance(ybox_kwargs, dict):
                ybox_kwargs = N_rows*[ybox_kwargs]
            elif len(ybox_kwargs) == 1: 
                ybox_kwargs = N_rows * ybox_kwargs
            xboxes = [Box(xbox_axes[i], **xbox_kwargs[i]) for i in range(N_cols)]
            yboxes = [Box(ybox_axes[i], **ybox_kwargs[i]) for i in range(N_rows)]
            plt.sca(ax_empty)
            plt.axis('off')
            for i, box in enumerate(xboxes):
                plt.sca(box.axis)
                plt.axis('off')
                if xticks is not None: 
                    try: plt.xlim([xticks[i][0], xticks[i][-1]])
                    except: pass
            for i, box in enumerate(yboxes):
                plt.sca(box.axis)
                plt.axis('off')
                if yticks is not None: 
                    try: plt.ylim([yticks[i][0], yticks[i][-1]])
                    except: pass
        else:
            fig, axes = plt.subplots(ncols=N_cols, nrows=N_rows)
            axes = axes.reshape([N_rows, N_cols])
    for i in range(N_rows):
        for j in range(N_cols):
            x = xs[i, j]
            y = ys[i, j]
            ax = axes[i, j]
            xticksj = None if xticks is None else xticks[j]
            yticksi = None if yticks is None else yticks[i]
            xticklabelsj = None if xticklabels is None else xticklabels[j]
            yticklabelsi = None if yticklabels is None else yticklabels[i]
            xticklabelsj = i == N_rows - 1
            yticklabelsi = j == 0
            plot_kde(x, y, ax=ax, fig=fig,
                     xbox=[False],
                     ybox=[False],
                     xticks=xticksj,
                     yticks=yticksi,
                     xticklabels=xticklabelsj,
                     yticklabels=yticklabelsi,
                     xtick0=True,
                     ytick0=i==N_rows-1,
                     xtickf=j==N_cols-1,
                     ytickf=i==0,
                     **kwargs)
    if xboxes is not None:
        for i, xbox in enumerate(xboxes):
            x = xs[:, i]
            plt.sca(xbox.axis)
            plot_uncertainty_boxes(x, xbox.light, xbox.dark,
                            positions=[*range(0, -len(x), -1)],
                            vertical=False, outliers=False, 
                            bounds=True)
            ax = axes[0, i]
            plt.sca(ax)
            xlim = plt.xlim()
            plt.sca(xbox.axis)
            plt.xlim(xlim)
    if yboxes is not None:
        for i, (ybox, y) in enumerate(zip(yboxes, ys)):
            plt.sca(ybox.axis)
            plot_uncertainty_boxes(y, ybox.light, ybox.dark,
                            positions=[*range(len(y))],
                            vertical=True, outliers=False, 
                            bounds=True)
            ax = axes[i, 0]
            plt.sca(ax)
            ylim = plt.ylim()
            plt.sca(ybox.axis)
            plt.ylim(ylim)
    if xlabel is not None: fig.supxlabel(xlabel, fontsize=fs)
    if ylabel is not None: fig.supylabel(ylabel, fontsize=fs)
    plt.subplots_adjust(hspace=0, wspace=0)
    plt.sca(ax)
    return fig, axes

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
        data0 = np.asarray(z_at_xy(x0, y0, *args))
        shape = data0.shape
        if len(shape) == 1:
            shape = f"({shape[0]})"
        if vectorize:
            N_args = len(args)
            Z_at_XY = np.vectorize(
                z_at_xy, signature=f'(),()->{shape}',
                excluded=tuple(range(2, 2 + N_args)),
            )
        else:
            Z_at_XY = z_at_xy
        Z = Z_at_XY(X, Y, *args)
        if smooth: # Smooth curves due to avoid discontinuities
            from scipy.ndimage.filters import gaussian_filter
            A, B, *other = Z.shape
            for index in product(*[range(i) for i in other]):
                Z[(..., *index)] = gaussian_filter(Z[(..., *index)], smooth)
    if file and save and not load: np.save(file, Z)
    return X, Y, Z

def plot_contour(
        X, Y, Z, 
        xlabel, ylabel, xticks, yticks, 
        metric_bars, titles=None, 
        fillcolor=None, styleaxiskw=None, label_size=10,
        label=False, wbar=1, contour_label_interval=2,
        highlight_levels=None, 
        highlight_color=None
    ): # pragma: no coverage
    """Create contour plots and return the figure and the axes."""
    if isinstance(metric_bars, MetricBar):
        metric_bars = [metric_bars]
        nrows = 1
        row_bars = True
    elif isinstance(metric_bars[0], MetricBar):
        nrows = len(metric_bars)
        row_bars = True
    else:
        nrows = len(metric_bars)
        ncols = len(metric_bars[0])
        row_bars = False
    A, B = Z.shape[:2]
    if Z.ndim == 2:
        ncols = 1
        Z = Z[:, :, None, None]
    elif Z.ndim == 3:
        if Z.shape == (*X.shape, nrows):
            Z = Z[:, :, :, None]
        else:
            Z = Z[:, :, None, :]
    assert Z.shape == (*X.shape, nrows, ncols), (
       f"Z was shape {Z.shape}, but expeted shape {(*X.shape, nrows, ncols)}; "
        "Z.shape must be (A, B, M, N), where (A, B) is the shape of both X and Y, "
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
        axes = axes.reshape([nrows, ncols])
        cbs = np.zeros([nrows, ncols], dtype=object)
    if styleaxiskw is None: styleaxiskw = ncols * [{}]
    if isinstance(styleaxiskw, dict): styleaxiskw = ncols * [styleaxiskw]
    cps = np.zeros([nrows, ncols], dtype=object)
    linecolor = np.array([*c.neutral_shade.RGBn, 0.1])
    other_axes = [[] for i in range(nrows)]
    if highlight_color is None: highlight_color = 'r'
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
                cs = plt.contour(cp, zorder=1, linewidths=0.8,
                                 levels=cp.levels, colors=[linecolor])
                clevels = [i for i in cp.levels[:-1][::contour_label_interval]]
                clabels = ax.clabel(
                    cs, levels=clevels,
                    inline=True, fmt=metric_bar.fmt,
                    colors=['k'], zorder=1
                )
                for i in clabels: i.set_rotation(0)
            if highlight_levels:
                cs = plt.contour(cp, zorder=1, linewidths=0.8,
                                 levels=highlight_levels, colors=[highlight_color])
                clabels = ax.clabel(
                    cs, levels=highlight_levels,
                    inline=True, fmt=metric_bar.fmt,
                    colors=[highlight_color], zorder=1
                )
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
                style_axis(ax, xticks, yticks, xticklabels, yticklabels, **styleaxiskw[col])
            )
        if row_bars:
            cbar_ax = axes[row, -1]
            cbs[row] = metric_bar.colorbar(fig, cbar_ax, cp, shrink=metric_bar.shrink,)
        
        # plt.clim()
    if titles is not None:
        for col in range(ncols):
            title = titles[col]
            ax = axes[0, col]
            ax.set_title(title, color=title_color, fontsize=10, fontweight='bold')
    if row_bars:
        for ax in axes[:, -1]:
            plt.sca(ax)
            plt.axis('off')
        # set_axes_labels(axes[:, :-1], xlabel, ylabel)
    # else:
        # set_axes_labels(axes, xlabel, ylabel)
    fig.supxlabel(xlabel, size=label_size)
    fig.supylabel(ylabel, size=label_size)
    plt.subplots_adjust(hspace=0.1, wspace=0.1)
    return fig, axes, cps, cbs, other_axes
plot_contour_2d = plot_contour       

def plot_contour_single_metric(
        X, Y, Z, xlabel, ylabel, xticks, yticks, metric_bar,
        titles=None, fillcolor=None, styleaxiskw=None, label=False,
        contour_label_interval=2, highlight_levels=None, highlight_color=None,
        label_fs=None,
    ): # pragma: no coverage
    """Create contour plots and return the figure and the axes."""
    if Z.ndim < 4:
        if Z.ndim == 3:
            Z = Z[:, :, :, None]
        elif Z.ndim == 2:
            Z = Z[:, :, None, None]
    *_, nrows, ncols = Z.shape
    assert Z.shape == (*X.shape, nrows, ncols), (
        "the first 2 dimensions of Z must have the same shape as X and Y"
    )
    fig, axes, ax_colorbar = contour_subplots(nrows, ncols, single_colorbar=True)
    if styleaxiskw is None: styleaxiskw = dict(xtick0=False, ytick0=False)
    cps = np.zeros([nrows, ncols], dtype=object)
    linecolor = np.array([*c.neutral_shade.RGBn, 0.1])
    other_axes = []
    if highlight_color is None: highlight_color = 'r'
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
            cp.set_edgecolors('face') # For svg background
            if label:
                cs = plt.contour(cp, zorder=1, linewidths=0.8,
                                 levels=cp.levels, colors=[linecolor])
                levels = levels=[i for i in cp.levels[:-1][::contour_label_interval]]
                clabels = ax.clabel(
                    cs, levels=levels,
                    inline=True, fmt=metric_bar.fmt,
                    colors=['k'], zorder=1, fontsize=label_fs,
                )
                for i in clabels: i.set_rotation(0)
            cps[row, col] = cp
            
            if highlight_levels:
                cs = plt.contour(cp, zorder=1, linewidths=0.8,
                                 levels=highlight_levels, colors=[highlight_color])
                clabels = ax.clabel(
                    cs, levels=highlight_levels,
                    inline=True, fmt=metric_bar.fmt,
                    colors=[highlight_color], zorder=1
                )
                for i in clabels: i.set_rotation(0)
            
            if row == nrows - 1 and 'ytick0' not in styleaxiskw:
                sak = styleaxiskw.copy()
                sak['ytick0'] = True
            else:
                sak = styleaxiskw
            if col == 0 and 'xtick0' not in styleaxiskw:
                sak = sak.copy()
                sak['xtick0'] = True
            dct = style_axis(ax, xticks, yticks, xticklabels, yticklabels, **sak)
            other_axes.append(dct)
    cb = metric_bar.colorbar(fig, ax_colorbar, cp, fraction=0.5)
    plt.sca(ax_colorbar)
    plt.axis('off')
    if titles:
        for col, title in enumerate(titles):
            ax = axes[0, col]
            ax.set_title(title, color=title_color, fontsize=10, fontweight='bold')
    fig = plt.gcf()
    if xlabel: fig.supxlabel(xlabel)
    if ylabel: fig.supylabel(ylabel)
    plt.subplots_adjust(hspace=0.1, wspace=0.1, bottom=0.2)
    return fig, axes, cps, cb, other_axes
            
def color_quadrants(color=None, x=None, y=None, xlim=None, ylim=None, 
                    line_color=None, fill_color=None, linewidth=1.0):
    if x is None: 
        xlb = 0
        xub = 0
    elif hasattr(x, '__len__'):
        xlb, xub = x
    else:
        xlb = xub = x
    if y is None: 
        ylb = 0
        yub = 0
    elif hasattr(y, '__len__'):
        ylb, yub = y
    else:
        ylb = yub = y
    if line_color is None: line_color = c.grey.RGBn
    if fill_color is None: fill_color = CABBI_colors.green_dirty.tint(85).RGBn
    if xlim is None: xlim = plt.xlim()
    if ylim is None: ylim = plt.ylim()
    x0, x1 = xlim
    y0, y1 = ylim
    top_left, top_right, bottom_left, bottom_right = color
    # Top left
    if top_left is not None:
        plt.fill_between([x0, xlb], yub, y1,
                         color=top_left,
                         linewidth=linewidth,
                         zorder=0)
    # Top right
    if top_right is not None:
        plt.fill_between([xub, x1], yub, y1,
                         color=top_right,
                         linewidth=linewidth,
                         zorder=0)
    # Bottom left
    if bottom_left is not None:
        plt.fill_between([x0, xlb], y0, ylb,
                         color=bottom_left,
                         linewidth=linewidth,
                         zorder=0)
    # Bottom right
    if bottom_right is not None:
        plt.fill_between([xub, x1], y0, ylb,
                         color=bottom_right,
                         linewidth=linewidth,
                         zorder=0)
    if yub == ylb: 
        plot_horizontal_line(ylb, line_color, zorder=0)
    else:
        plt.fill_between([x0, x1], ylb, yub,
                         color=fill_color,
                         linewidth=linewidth,
                         zorder=0)
        plot_horizontal_line(ylb, line_color, zorder=0)
        plot_horizontal_line(yub, line_color, zorder=0)
    if xub == xlb: 
        plot_vertical_line(xlb, line_color, zorder=0)
    else:
        plt.fill_between([xlb, xub], y0, y1,
                         color=fill_color,
                         linewidth=linewidth,
                         zorder=0)
        plot_vertical_line(xlb, line_color, zorder=0)
        plot_vertical_line(xub, line_color, zorder=0)
    
def label_quadrants(
        x=None, y=None, xr=None, yr=None, text=None, color=None,
        fs=None,
    ):
    if xr is None: 
        xrlb = 0
        xrub = 0
    elif hasattr(xr, '__len__'):
        xrlb, xrub = xr
    else:
        xrlb = xrub = xr
    if yr is None: 
        yrlb = 0
        yrub = 0
    elif hasattr(yr, '__len__'):
        yrlb, yrub = yr
    else:
        yrlb = yrub = yr
    xlb, xub = plt.xlim()
    ylb, yub = plt.ylim()
    data_given = not (x is None or y is None)
    if data_given:
        y_mt_0 = y > yrub
        y_lt_0 = y < yrlb
        x_mt_0 = x > xrub
        x_lt_0 = x < xrlb
    xpos = lambda x: xlb + (xub - xlb) * x
    ypos = lambda y: ylb + (yub - ylb) * y
    xleft = 0.02
    xright = 0.98
    ytop = 0.94
    ybottom = 0.02
    labeled = 4 * [False]
    top_left, top_right, bottom_left, bottom_right = text
    top_left_color, top_right_color, bottom_left_color, bottom_right_color = color
    if yub > yrub and xlb < xrlb and top_left:
        if data_given and top_left.endswith('()'):
            p = (y_mt_0 & x_lt_0).sum() / y.size
            top_left = f"{p:.0%} {top_left.strip('()')}"
        plt.text(xpos(xleft), ypos(ytop), top_left, color=top_left_color,
                 horizontalalignment='left', verticalalignment='top',
                 fontsize=fs, fontweight='bold', zorder=10)
        labeled[0] = True
    if yub > yrub and xub > xrub and top_right:
        if data_given and top_right.endswith('()'):
            p = (y_mt_0 & x_mt_0).sum() / y.size
            top_right = f"{p:.0%} {top_right.strip('()')}"
        plt.text(xpos(xright), ypos(ytop), top_right, color=top_right_color,
                 horizontalalignment='right', verticalalignment='top',
                 fontsize=fs, fontweight='bold', zorder=10)
        labeled[1] = True
    if ylb < yrlb and xlb < xrlb and bottom_left:
        if data_given and bottom_left.endswith('()'):
            p = (y_lt_0 & x_lt_0).sum() / y.size
            bottom_left = f"{p:.0%} {bottom_left.strip('()')}"
        plt.text(xpos(xleft), ypos(ybottom), bottom_left, color=bottom_left_color,
                 horizontalalignment='left', verticalalignment='bottom',
                 fontsize=fs, fontweight='bold', zorder=10)
        labeled[2] = True
    if ylb < yrlb and xub > xrub and bottom_right:
        if data_given and bottom_right.endswith('()'):
            p = (y_lt_0 & x_mt_0).sum() / y.size
            bottom_right = f"{p:.0%} {bottom_right.strip('()')}"
        plt.text(xpos(xright), ypos(ybottom), bottom_right, color=bottom_right_color,
                 horizontalalignment='right', verticalalignment='bottom',
                 fontsize=fs, fontweight='bold', zorder=10)
        labeled[3] = True
    return labeled

def plot_quadrants(
        text, data=None, x=None, y=None, rotate=0, fs=None,
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
        data, x, y, text, text_color, quadrant_color, fs
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
        fs=None,
    ):
    if data is None: data = (None, None) 
    color_quadrants(quadrant_color, x, y, xlim, ylim)
    return label_quadrants(
        *data, x, y, text, text_color, fs
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
