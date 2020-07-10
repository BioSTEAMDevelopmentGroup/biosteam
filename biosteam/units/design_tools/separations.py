# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
This module contains unit operation separation methods.

"""
from thermosteam.exceptions import InfeasibleRegion

__all__ = ('lle', 'vle', 'split', 'mix_and_split',
           'adjust_moisture_content', 
           'mix_and_split_with_moisture_content')

def mix_and_split_with_moisture_content(ins, retentate, permeate,
                                        split, moisture_content):
    """
    Run splitter mass and energy balance with mixing all input streams and 
    and ensuring retentate moisture content.
    
    Parameters
    ----------
    ins : Iterable[Stream]
        Inlet fluids with solids.
    retentate : Stream
    permeate : Stream
    split : array_like
        Component splits to the retentate.
    moisture_content : float
        Fraction of water in retentate.

    """
    mix_and_split(ins, retentate, permeate, split)
    adjust_moisture_content(retentate, permeate, moisture_content)

def adjust_moisture_content(retentate, permeate, moisture_content):
    """
    Adjust retentate moisture content with water from permeate.
    
    Parameters
    ----------
    retentate : Stream
    permeate : Stream
    moisture_content : float
        Fraction of water in retentate.

    """
    F_mass = retentate.F_mass
    mc = moisture_content
    retentate.imol['7732-18-5'] = water = (F_mass * mc/(1-mc))/18.01528
    permeate.imol['7732-18-5'] -= water
    if permeate.imol['7732-18-5'] < 0:
        raise InfeasibleRegion('not enough water; permeate moisture content')

def mix_and_split(ins, top, bottom, split):
    """
    Run splitter mass and energy balance with mixing all input streams.
    
    Parameters
    ----------
    ins : Iterable[Stream]
        All inlet fluids.
    top : Stream
        Top inlet fluid.
    bottom : Stream
        Bottom inlet fluid
    split : array_like
        Component-wise split of feed to the top stream.
    
    """
    top.mix_from(ins)
    bottom.copy_like(top)
    top_mol = top.mol
    top_mol[:] *= split
    bottom.mol[:] -= top_mol

def split(feed, top, bottom, split):
    """
    Run splitter mass and energy balance with mixing all input streams.
    
    Parameters
    ----------
    feed : Stream
        Inlet fluid.
    top : Stream
        Top inlet fluid.
    bottom : Stream
        Bottom inlet fluid
    split : array_like
        Component-wise split of feed to the top stream.
    
    """
    if feed is not top: top.copy_like(feed)
    bottom.copy_like(top)
    top_mol = top.mol
    top_mol[:] *= split
    bottom.mol[:] -= top_mol

def lle(feed, top, bottom, top_chemical, efficiency, multi_stream=None):
    """Run LLE mass and energy balance."""
    if multi_stream:
        ms = multi_stream
        ms.copy_like(feed)
    else:
        ms = feed.copy()
    ms.lle(feed.T, top_chemical=top_chemical)
    top_phase = 'l'
    bottom_phase = 'L'
    if not top_chemical:
        rho_l = ms['l'].rho
        rho_L = ms['L'].rho
        top_L = rho_L < rho_l
        if top_L:
            top_phase = 'L'
            bottom_phase = 'l'
    top.mol[:] = ms.imol[top_phase]
    bottom.mol[:] = ms.imol[bottom_phase]
    top.T = bottom.T = feed.T
    top.P = bottom.P = feed.P
    if efficiency < 1.:
        top.mol *= efficiency
        bottom.mol *= efficiency
        mixing = (1 - efficiency) / 2 * feed.mol
        top.mol += mixing
        bottom.mol += mixing
        
def vle(feed, vap, liq, T=None, P=None, V=None, Q=None, x=None, y=None,
        multi_stream=None):
    """Run VLE mass and energy balance."""
    if multi_stream:
        ms = multi_stream
        ms.copy_like(feed)
    else:
        ms = feed.copy()
    H = feed.H + Q if Q is not None else None
    ms.vle(P=P, H=H, T=T, V=V, x=x, y=y)

    # Set Values
    vap.phase = 'g'
    liq.phase = 'l'
    vap.mol[:] = ms.imol['g']
    liq.mol[:] = ms.imol['l']
    vap.T = liq.T = ms.T
    vap.P = liq.P = ms.P