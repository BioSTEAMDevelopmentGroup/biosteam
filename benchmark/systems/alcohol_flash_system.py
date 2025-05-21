# -*- coding: utf-8 -*-
"""
Created on Sat Mar 16 13:38:11 2024

@author: cortespea
"""
import biosteam as bst

__all__ = (
    'create_system_alcohol_narrow_flash',
    'create_system_alcohol_wide_flash',
)

def create_system_alcohol_narrow_flash(alg):
    bst.settings.set_thermo(['heptanol', 'octanol'], cache=True, Gamma=bst.IdealActivityCoefficients)
    feed = bst.Stream('feed', heptanol=100, octanol=100)
    recycle = bst.Stream('liquid_recycle')
    liquid_product = bst.Stream('liquid_product')
    vapor_product = bst.Stream('vapor_product')
    stage = bst.StageEquilibrium(
        'stage', ins=[feed, recycle], outs=[vapor_product, recycle, liquid_product], 
        B=1, bottom_split=0.4, phases=('g', 'l')
    )
    sys = bst.System.from_units('sys', [stage])
    return sys

def create_system_alcohol_wide_flash(alg):
    bst.settings.set_thermo(['propanol', 'octanol'], cache=True, Gamma=bst.IdealActivityCoefficients)
    feed = bst.Stream('feed', propanol=100, octanol=100)
    recycle = bst.Stream('liquid_recycle')
    liquid_product = bst.Stream('liquid_product')
    vapor_product = bst.Stream('vapor_product')
    stage = bst.StageEquilibrium(
        'stage', ins=[feed, recycle], outs=[vapor_product, recycle, liquid_product], 
        B=1, bottom_split=0.4, phases=('g', 'l')
    )
    sys = bst.System.from_units('sys', [stage])
    return sys
