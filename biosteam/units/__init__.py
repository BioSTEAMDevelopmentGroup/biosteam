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
from .flash import Flash, Evaporator_PQin, Evaporator_PV, LazyFlash
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
from .boiler_turbogenerator import BoilerTurbogenerator
from .vacuum import vacuum_system

# %% All units

__all__ = ['Mixer', 'Splitter', 'InvSplitter', 'Tank', 'MixTank', 'StorageTank', 'Separator', 'HX', 'HXutility', 'HXprocess', 'Pump', 'Flash', 'Evaporator_PQin', 'Evaporator_PV', 'LazyFlash', 'Distillation', 'Stripper', 'Transesterification', 'Fermentation', 'Centrifuge_LLE', 'Centrifuge_LLE_Lazy', 'MultiEffectEvaporator', 'EnzymeTreatment', 'CrushingMill', 'RVF', 'JuiceScreener', 'Screener', 'MolecularSieve', 'MolSieve', 'PCT', 'Settler', 'YeastCentrifuge', 'Clarifier', 'Reactor', 'BatchReactor', 'MassBalance', 'EnergyBalance', 'BoilerTurbogenerator', 'vacuum_system']


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