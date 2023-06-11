# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from thermosteam.utils import plots
from thermosteam.utils.plots import *
from biosteam.utils import colors
from typing import NamedTuple, Iterable, Callable
from matplotlib.colors import Colormap, TwoSlopeNorm
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LinearSegmentedColormap
from matplotlib.cm import ScalarMappable

__all__ = ('CABBI_green_colormap', 'MetricBar',
           'expand', 'closest_index', 'color_bar', *plots.__all__)

# %% Data classes

class MetricBar(NamedTuple): # pragma: no coverage
    name: str = None
    units: str = None
    cmap: Colormap = None
    ticks: Iterable[float] = None
    N_levels: int = 20
    N_decimals: int = 0
    units_dlim: str = '\n'
    ub: bool = False
    lb: bool = False
    forced_size: float = None
    ylabelkwargs: dict = {}
    center: bool = None
    title_position: str = 'top'
    shrink: bool = 0.8
    
    def fmt(self, x):
        value = f'{round(x, self.N_decimals):,}'
        value = value.rstrip('0')
        if value[-1] == '.':
            value = value[:-1]
        return value
    
    @property
    def levels(self):
        ticks = self.ticks 
        if ticks is None: 
            return None
        else:
            return np.linspace(ticks[0], ticks[-1], self.N_levels)
    
    @property
    def title(self):
        if self.units:
            return f'{self.name}{self.units_dlim}[{self.units}]'
        else:
            return self.name
    
    @property
    def norm(self):
        if self.center is not None:
            return TwoSlopeNorm(vmin=self.ticks[0], 
                                vcenter=self.center, 
                                vmax=self.ticks[-1])
    
    def colorbar(self, fig, ax, colorplot, label=True, **cbarkwargs):
        if self.forced_size is not None:
            cbarkwargs['fraction'] = self.forced_size
        cbar = fig.colorbar(colorplot, ax=ax, ticks=self.ticks, **cbarkwargs)
        cbar_ax = cbar.ax
        # cbar_ax.locator_params(nbins=self.N_ticks)
        if label:
            if self.title_position == 'top':
                cbar_ax.set_title(self.title, self.ylabelkwargs)
            else:
                cbar_ax.set_ylabel(self.title, self.ylabelkwargs)
        try:
            ylabels = [y.get_text() for y in cbar_ax.get_yticklabels()]
            ylabels = [(i if i[0].isdigit() else '-'+i[1:]) for i in ylabels]
            if self.ub:
                ylabels[-1] = '>' + ylabels[-1]
            if self.lb:
                ylabels[0] = '<' + ylabels[0]
            cbar_ax.set_yticklabels(ylabels)
        except:
            pass
        return cbar
        
# %% Helpful functions

def CABBI_green_colormap(N_levels=25): # pragma: no coverage
    """
    Return a matplotlib.colors.LinearSegmentedColormap object
    that serves as CABBI's green colormap theme for contour plots.
    
    """
    CABBI_colors = (colors.CABBI_yellow.RGBn,
                    colors.CABBI_green.RGBn,
                    colors.CABBI_teal_green.shade(75).RGBn)
    return LinearSegmentedColormap.from_list('CABBI', CABBI_colors, N_levels)

def color_bar(RGBn: list, vmin=0, vmax=100, ax=None,
              label=None, orientation='vertical', N_levels=25): # pragma: no coverage
    cmap = LinearSegmentedColormap.from_list(label, RGBn, N_levels)
    norm = Normalize(vmin=vmin, vmax=vmax)
    if not ax: 
        if orientation == 'vertical':
            shape = (0.5, 5)
        elif orientation == 'horizontal':
            shape = (5, 0.5)
        else:
            raise ValueError("orientation must be either 'vertical' or 'horizonta'; not %s" %orientation)
        fig, ax = plt.subplots(figsize=shape)
    return plt.colorbar(ScalarMappable(norm=norm, cmap=cmap), 
                        orientation=orientation, 
                        cax=ax, label=label)
    

def expand(lower, upper, lb, ub): # pragma: no coverage
    dx = (upper - lower)/12
    return max(lower-dx, lb), min(upper+dx, ub)

def closest_index(x, xs): # pragma: no coverage
    if xs[0] < xs[-1]:
        for i, xi in enumerate(xs):
            if x < xi: break
    else:
        for i, xi in enumerate(xs):
            if x > xi: break
    return i

del plots