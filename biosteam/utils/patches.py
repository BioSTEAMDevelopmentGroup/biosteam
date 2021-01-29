# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from thermosteam.utils import colors
import numpy as np

__all__ = ('DoubleColorCircle', 'DoubleColorLegend')

# %% Legend handlers

def make_patch(pos, width, height, Patch, handlepatch, kwargs): # pragma: no cover
    if issubclass(Patch , patches.Circle): 
        patch = Patch(pos+[width/2,height/2], height/2,
                      transform=handlepatch.get_transform(), **kwargs)    
    elif issubclass(Patch, patches.Rectangle): 
        patch = Patch(pos, width, height,
                      transform=handlepatch.get_transform(), **kwargs)    
    else:
        raise NotImplementedError(f"patch '{Patch}' not implemented")
    handlepatch.add_artist(patch)
    return patch

class DoubleColor: # pragma: no cover
    __slots__ = ('leftpatches', 'rightpatches', 'Patch', 'patches', 'both')
    
    def __init__(self, leftpatches, rightpatches, morepatches=None,
                 shape='Rectangle', both={}):
        self.leftpatches = leftpatches
        self.rightpatches = rightpatches
        self.Patch = getattr(patches, shape)
        self.patches = morepatches or []
        self.both = both
        
    def legend_artist(self, legend, handle, fontsize, patch):
        x0, y0 = patch.xdescent, patch.ydescent
        width, height = patch.width/2, patch.height
        leftpos = np.array([x0, y0])
        rightpos = np.array([x0+width, y0])
        Patch = self.Patch
        for l in self.leftpatches:
            lpatch = make_patch(leftpos, width, height, Patch, patch, l)
        for r in self.rightpatches:
            rpatch = make_patch(rightpos, width, height, Patch, patch, r)
            patch.add_artist(rpatch)
        for patch in self.patches:
            patch.legend_artist(legend, handle, fontsize, patch)
        both = self.both
        if both:
            if 'fill' not in both:
                both['fill'] = False
            if 'edgecolor' not in both:
                both['edgecolor'] = 'k'
            patch.add_artist(patches.Rectangle(leftpos, width*2, height,
                transform=patch.get_transform(), **both))
        return lpatch


class DoubleColorCircle: # pragma: no cover
    __slots__ = ('left', 'right', 'both')
    
    def __init__(self, left, right, both):
        self.left = left
        self.right = right
        self.both = both
        
    def legend_artist(self, legend, handle, fontsize, patch):
        x0, y0 = patch.xdescent, patch.ydescent
        width, height = patch.width, patch.height
        pos = np.array([x0+width/2, y0+height/2])
        leftpatch = patches.Arc(pos, width/2, width/2,
                                angle=.0, theta1= 90.0, theta2=-90.0,
                                **self.left, hatch = 20*'o',
                                transform=patch.get_transform())
        rightpatch = patches.Arc(pos, width/2, width/2,
                                angle=.0, theta1=-90.0, theta2= 90.0,
                                **self.right, hatch = 20*'o',
                                transform=patch.get_transform())
        patch.add_artist(leftpatch)
        patch.add_artist(rightpatch)
        both = self.both
        if both:
            if 'fill' not in both:
                both['fill'] = False
            if 'edgecolor' not in both:
                both['edgecolor'] = 'k'
            patch.add_artist(patches.Arc(pos, width/2, width/2,
                                         transform=patch.get_transform(),
                                         **both))
        return leftpatch


class TwoColorArrow: # pragma: no cover
    __slots__ = ('left', 'right', 'pos')
    
    def __init__(self, left, right, pos='mid'):
        self.left = left
        self.right = right
        self.pos = pos
        
    def legend_artist(self, legend, handle, fontsize, patch):
        x0, y0 = patch.xdescent, patch.ydescent
        width = patch.width
        height = patch.height
        pos = self.pos
        if pos=='mid':
            y0 += height/2
        elif pos=='bot':
            pass
        elif pos=='top':
            y0 += height
        else:
            ValueError('invalid pos in ArrowPatch object')
        leftpatch = patches.Arrow(x0, y0, width/2, 0,
                                  **self.left)
        rightpatch = patches.Arrow(x0+width/2, y0, width/2, 0,
                                   **self.right)
        patch.add_artist(leftpatch)
        patch.add_artist(rightpatch)
        return leftpatch


class Legend: # pragma: no cover
    __slots__ = ('handler_map')
    
    def __init__(self, handler_map=None):
        self.handler_map = handler_map or {}
    
    def legend(self, loc='upper right', **kwargs):
        hmap = self.handler_map
        names = hmap.keys()
        plt.legend(names, names, handler_map=self.handler_map, prop=kwargs,
                   loc=loc)


class DoubleColorLegend(Legend): # pragma: no cover
    __slots__ = ()
    
    def __init__(self, key_leftrightpatches=None, handler_map=None):
        self.handler_map = handler_map or {}
        if key_leftrightpatches:
            for key, left_right in key_leftrightpatches.items():
                handler_map[key] = DoubleColor(*left_right)
    
    def add_box(self, name,
                leftcolor=colors.orange_tint.RGBn, 
                rightcolor=colors.blue_tint.RGBn,
                both={}):  
        self.handler_map[name] = DoubleColor([{'facecolor': leftcolor}],
                                             [{'facecolor': rightcolor}],
                                             None, 'Rectangle', both=both)

    def add_circle(self, name,
                   leftcolor=colors.orange_shade.RGBn,
                   rightcolor=colors.blue_shade.RGBn,
                   both={}):
        self.handler_map[name] = DoubleColorCircle(
                   {'color': leftcolor,
                    'edgecolor': 'black'},
                   {'color': rightcolor,
                    'edgecolor': 'black'},
                    both=both)
    
    # def add2key(self, key, patch, side='left'):
    #     db = self.handler_map[key]
    #     sidepatches = getattr(db, side+'patches')
    #     sidepatches.append(patch)
    
    # def newkey(self, key, leftpatches, rightpatches, patches=None, shape='Rectangle', **kwargs):
    #     if isinstance(leftpatches, dict):
    #         leftpatches = [leftpatches]
    #     if isinstance(rightpatches, dict):
    #         rightpatches = [rightpatches]
    #     self.handler_map[key] = DoubleColor(leftpatches, rightpatches, patches, shape, **kwargs)