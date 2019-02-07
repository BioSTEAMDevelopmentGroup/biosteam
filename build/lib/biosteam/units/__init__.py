#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 15 20:39:46 2018

@author: Yoel Rene Cortes-Pena
"""
import numpy as np
from biosteam.graphics import Graphics
from .mixer import Mixer
from .splitter import Splitter, InvSplitter
from .pump import Pump
from .hx import HX, HXutility, HXprocess
from .flash import Flash, Flash_PQin, Flash_PV, P69_flash
from .multi_effect_evaporator import MultiEffectEvaporator
from .centrifuge_LLE import Centrifuge_LLE, Centrifuge_LLE_Lazy
from .distillation import Dist, Distillation, Stripper
from .tank import Tank, StorageTank, MixTank, PCT
from .transesterification import Transesterification
from .fermentation import Fermentation
from .enzyme_treatment import EnzymeTreatment
from .settler import Settler
from .clarifier import Clarifier
from .yeast_centrifuge import YeastCentrifuge
from .crushing_mill import CrushingMill
from .RVF import RVF
from .juice_screener import JuiceScreener
from .screener import Screener
from .molecular_sieve import MolecularSieve, MolSieve
from .separator import Separator
from .reactor import Reactor, BatchReactor
from .balance import MassBalance, EnergyBalance

# %% All units

__all__ = ['Mixer', 'Splitter', 'InvSplitter', 'Tank', 'MixTank', 'StorageTank', 'Separator', 'HX', 'HXutility', 'HXprocess', 'Pump', 'Flash', 'Flash_PQin', 'Flash_PV', 'P69_flash', 'Distillation', 'Stripper', 'Transesterification', 'Fermentation', 'Centrifuge_LLE', 'Centrifuge_LLE_Lazy', 'MultiEffectEvaporator', 'EnzymeTreatment', 'CrushingMill', 'RVF', 'JuiceScreener', 'Screener', 'MolecularSieve', 'MolSieve', 'PCT', 'Settler', 'YeastCentrifuge', 'Clarifier', 'Reactor', 'BatchReactor', 'MassBalance', 'EnergyBalance']


# %% Enhance Graphics

# Flash
edge_out = Flash._Graphics.edge_out
edge_out[0]['tailport'] = 'n'
edge_out[1]['tailport'] = 's'
node = Flash._Graphics.node
node['width'] = '1'
node['height'] = '1.1'

# Distillation
edge_out = Dist._Graphics.edge_out
edge_out[0]['tailport'] = 'n'
edge_out[1]['tailport'] = 's'
node = Dist._Graphics.node
node['width'] = '1'
node['height'] = '1.2'

# Single stream heat exchanger
HXutility._Graphics = Graphics()
HXutility._Graphics.node['shape'] = 'circle'
HXutility._Graphics.node['color'] = 'none'
HXutility._Graphics.node['margin'] = '0'

def HXutility_node(hx):
    si = hx.ins[0]
    so = hx.outs[0]
    gi = si.phase == 'g'
    li = so.phase == 'l'
    go = so.phase == 'g'
    lo = so.phase == 'l'
    Ti = si.T
    To = so.T
    if Ti > To or (gi and lo):
        hx._Graphics.node['fillcolor'] = '#cfecf0'
    elif Ti < To or (li and go):
        hx._Graphics.node['fillcolor'] = '#fad6d8'
    else:
        hx._Graphics.node['fillcolor'] = '#cfecf0:#fad6d8'

HXutility._Graphics.node_function = HXutility_node

# Double stream heat exchanger
HXprocess._Graphics = Graphics()
HXprocess._Graphics.node['shape'] = 'circle'
HXprocess._Graphics.node['color'] = 'none'
HXprocess._Graphics.node['margin'] = '0'
HXprocess._Graphics.node['gradientangle'] = '90'
HXprocess._Graphics.node['fillcolor'] = '#cfecf0:#fad6d8'

# Mixer
Mixer._Graphics.node['shape'] = 'triangle'
Mixer._Graphics.node['orientation'] = '270'
Mixer._Graphics.edge_out[0]['tailport'] = 'e'

# Splitter
Splitter._Graphics.node['shape'] = 'triangle'
Splitter._Graphics.node['orientation'] = '90'
Splitter._Graphics.node['fillcolor'] = "#bfbfbf:white"
Splitter._Graphics.edge_in[0]['headport'] = 'w'

# Balance
MBG = MassBalance._Graphics
MBG.node['shape'] = 'note'
MBG.node['fillcolor'] = '#F0F0F0'
MBG.in_system = False