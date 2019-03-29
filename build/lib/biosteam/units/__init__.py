#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 15 20:39:46 2018

@author: Yoel Rene Cortes-Pena
"""
import numpy as np
from . import designtools
from . import flash
from . import centrifuge_LLE

from biosteam.graphics import Graphics
from .designtools import *
from .mixer import Mixer
from .splitter import Splitter, InvSplitter
from .pump import Pump
from .hx import HX, HXutility, HXprocess
from .flash import *
from .flash import Flash
from .multi_effect_evaporator import MultiEffectEvaporator
from .centrifuge_LLE import *
from .distillation import Dist, Distillation, Stripper
from .tank import Tank, StorageTank, MixTank, PCT
from .transesterification import Transesterification
from .fermentation import Fermentation
from .enzyme_treatment import EnzymeTreatment
from .settler import Settler
from .clarifier import Clarifier
from .yeast_centrifuge import YeastCentrifuge
from .crushing_mill import CrushingMill
from .rvf import RVF
from .juice_screener import JuiceScreener
from .screener import Screener
from .molecular_sieve import MolecularSieve, MolSieve
from .separator import Separator
from .reactor import Reactor, BatchReactor
from .balance import MassBalance, EnergyBalance
from .conveying_belt import ConveyingBelt
from .shredder import Shredder
from .magnetic_separator import MagneticSeparator
from .screw_feeder import ScrewFeeder
from .vibrating_screen import VibratingScreen

# %% All units

__all__ = ['Mixer', 'Splitter', 'InvSplitter', 'Tank', 'MixTank', 'StorageTank', 'Separator', 'HX', 'HXutility', 'HXprocess', 'Pump', 'Distillation', 'Stripper', 'Transesterification', 'Fermentation', 'Centrifuge_LLE', 'MultiEffectEvaporator', 'EnzymeTreatment', 'CrushingMill', 'RVF', 'JuiceScreener', 'Screener', 'MolecularSieve', 'MolSieve', 'PCT', 'Settler', 'YeastCentrifuge', 'Clarifier', 'Reactor', 'BatchReactor', 'MassBalance', 'EnergyBalance', 'ConveyingBelt', 'Shredder', 'MagneticSeparator', 'ScrewFeeder', 'VibratingScreen']

__all__.extend(designtools.__all__)
__all__.extend(flash.__all__)

# %% Enhance Graphics

# Flash
edge_out = Flash._graphics.edge_out
edge_out[0]['tailport'] = 'n'
edge_out[1]['tailport'] = 's'
node = Flash._graphics.node
node['width'] = '1'
node['height'] = '1.1'

# Distillation
edge_out = Dist._graphics.edge_out
edge_out[0]['tailport'] = 'n'
edge_out[1]['tailport'] = 's'
node = Dist._graphics.node
node['width'] = '1'
node['height'] = '1.2'

# Single stream heat exchanger
HXutility._graphics = graphics = Graphics()
graphics.node['shape'] = 'circle'
graphics.node['color'] = 'none'
graphics.node['margin'] = '0'

def HXutility_node(hx):
    try:
        si = hx.ins[0]
        so = hx.outs[0]
        gi = si.phase == 'g'
        li = so.phase == 'l'
        go = so.phase == 'g'
        lo = so.phase == 'l'
        Ti = si.T
        To = so.T
        graphics = hx._graphics
        if Ti > To or (gi and lo):
            graphics.node['fillcolor'] = '#cfecf0'
            graphics.name = 'Cooling'
        elif Ti < To or (li and go):
            graphics.node['fillcolor'] = '#fad6d8'
            graphics.name = 'Heating'
        else:
            graphics.node['fillcolor'] = '#cfecf0:#fad6d8'
    except:
        graphics = hx._graphics
        graphics.name = 'Heat exchange'

graphics.node_function = HXutility_node

# Double stream heat exchanger
HXprocess._graphics = graphics = Graphics()
graphics.name = 'HXprocess'
graphics.node['shape'] = 'circle'
graphics.node['color'] = 'none'
graphics.node['margin'] = '0'
graphics.node['gradientangle'] = '90'
graphics.node['fillcolor'] = '#cfecf0:#fad6d8'

# Mixer
Mixer._graphics.node['shape'] = 'triangle'
Mixer._graphics.node['orientation'] = '270'
Mixer._graphics.edge_out[0]['tailport'] = 'e'

# Splitter
Splitter._graphics.node['shape'] = 'triangle'
Splitter._graphics.node['orientation'] = '90'
Splitter._graphics.node['fillcolor'] = "#bfbfbf:white"
Splitter._graphics.edge_in[0]['headport'] = 'w'

# Balance
MBG = MassBalance._graphics
MBG.node['shape'] = 'note'
MBG.node['fillcolor'] = '#F0F0F0'
MBG.in_system = False