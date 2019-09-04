#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 21 11:05:24 2017

The oil and sugar separation (pretreatment) section for the baseline lipid cane biorefinery is defined here as System objects. The systems include all streams and units starting from enzyme treatment to purification of the sugar solution and the oil stream.

@author: Yoel
"""
import numpy as np
from biosteam import System, Stream
from biosteam.units import Mixer, EnzymeTreatment, CrushingMill, \
    HXutility, RVF, SplitFlash, VibratingScreen, \
    MagneticSeparator, Clarifier, MixTank, \
    Shredder, ConveyingBelt, Splitter, \
    SplitCentrifuge_LLE, Pump, StorageTank
from biosteam.biorefineries.lipidcane.species import pretreatment_species
from biosteam.biorefineries.lipidcane.process_settings import price

__all__ = ('pretreatment_sys', 'Lipid_cane')

# %% Species

Stream.species = pretreatment_species
psp = ('Ash', 'CaO', 'Cellulose', 'Ethanol', 'Flocculant',
       'Glucose', 'Hemicellulose', 'Lignin', 'Lipid',
       'Solids', 'H3PO4', 'Sucrose', 'Water')
psp1 = ('Ash', 'Cellulose', 'Glucose', 'Hemicellulose',
        'Lignin', 'Lipid', 'Solids', 'Sucrose', 'Water')
psp2 = ('Ash', 'CaO', 'Cellulose', 'Flocculant', 'Glucose',
        'Hemicellulose', 'Lignin', 'Lipid',
        'H3PO4', 'Sucrose', 'Water')

# %% Streams

f1 = (2000.042, 26986.69 , 2007.067, 15922.734, 14459.241,
      10035.334, 5017.667, 22746.761, 234157.798)
Lipid_cane = Stream('Lipid cane', f1, psp1, units='kg/hr',
                    price=price['Lipid cane'])

enzyme = Stream('Enzyme', Cellulose=100, Water=900, units='kg/hr',
                price=price['Protease'])

imbibition_water = Stream('Imbibition water',
                          Water=87023.35,
                          T = 338.15, units='kg/hr')

H3PO4 = Stream('H3PO4', H3PO4=74.23, Water=13.10, units='kg/hr',
               price=price['H3PO4'])  # to P9

lime = Stream('Lime', CaO=333.00, Water=2200.00, units='kg/hr',
              price=price['Lime'])  # to P5

polymer = Stream('Polymer', Flocculant=0.83, units='kg/hr',
                 price=price['Polymer'])  # to P68

wash_water = Stream('Water 1', Water=16770, units='kg/hr')  # to P14
wash_water.T = 363.15

lipid_wash = Stream('Water 2', Water=1350, units='kg/hr')  # to P75
lipid_wash.T = 358.15

S254 = Stream('S254', Ash=1, units='kg/hr')  # to P46

# %% Units

# Feed the shredder
F1 = ConveyingBelt('F1', ins=Lipid_cane)
F1.cost_items['Conveying belt'].ub = 2500

# Separate metals
F2 = MagneticSeparator('F2', ins=F1.outs)

# Shred fresh cane
F3 = Shredder('F3', ins=F2.outs)

# Hydrolyze starch
P137 = EnzymeTreatment('P137', T=323.15)  # T=50

# Finely crush lipid cane
Mill = CrushingMill('Mill',
                    split=(0.92, 0.92, 0.04, 0.92, 0.92, 0.04, 0.1, 1),
                    order=('Ash', 'Cellulose', 'Glucose', 'Hemicellulose',
                           'Lignin', 'Sucrose', 'Lipid', 'Solids'),
                    moisture_content=0.5)

# Convey out bagasse
F4 = ConveyingBelt('F4', ins=Mill.outs[0], outs='Bagasse')

# Mix in water
P21 = Mixer('P21')

# Screen out fibers
P56 = VibratingScreen('P56',
                      split=(0.35, 0.35, 0.88, 0.35, 0.35, 0.88, 0, 0.88, 0.88),
                      order=psp1)

# Store juice before treatment
P1 = StorageTank('P1')
P1.tau = 12

# Heat up before adding acid
P3 = HXutility('P3', T=343.15)

# Mix in acid
P9 = MixTank('P9')

# Pump acid solution
F5 = Pump('F5')

# Mix lime solution
F6 = MixTank('F6')
F6.tau = 1
F7 = Pump('F7')

# Blend acid lipid solution with lime
P5 = MixTank('P5')

# Mix recycle
P4 = Mixer('P4')

# Heat before adding flocculant
P7 = HXutility('P7', T=372.15)

# Mix in flocculant
P68 = MixTank('P68')
P68.tau = 1/4

# Separate residual solids
P12 = Clarifier('P12',
                split=(0, 0, 0, 0.522, 0.522, 0, 0,
                       0.98, 0.522, 0.522, 0.522),
                order=psp2)

# Remove solids as filter cake
P14 = RVF('P14', 
          moisture_content=0.80,
          split=(0.85, 0.85, 0.85, 0.01, 0.85, 0.85, 0.01),
          order=('Ash', 'CaO', 'Cellulose', 'Glucose',
                 'Hemicellulose', 'Lignin', 'Sucrose'))
F8 = Pump('F8')

# Get filter cake
P46 = Mixer('P46', outs='Filter cake')

# Separate oil and sugar
P10 = MixTank('P10', outs=('', ''))
split = np.zeros(len(pretreatment_species), float)
index = pretreatment_species.indices(('Lipid', 'Water'))
split[index] = (1, 0.0001)
P10._split = split
P10._run = lambda : Splitter._run(P10)
del split, index

# Cool the oil
P49 = HXutility('P49', T=343.15)

# Screen out small fibers from sugar stream
P50 = VibratingScreen('P50', outs=('Sugar', 'Fiber fines'),
                      split=1-np.array((0, 0, 0, 1, 0.002, 0, 0,0, 0, 0.002, 0.002)),
                      order=psp2)
Sugar = P50-0
P50.mesh_opening = 2

# Add distilled water to wash lipid
P75 = MixTank('P75')
P75.tau = 2

# Centrifuge out water
P76 = SplitCentrifuge_LLE('P76',
                          split=(0.99, 0.01),
                          order=('Lipid', 'Water'))

# Vacume out water
P69 = SplitFlash('P69', T=347.15, P=2026.5,
                 split=(0.0001, 0.999), order=('Lipid', 'Water'))
Lipid = P69.outs[1]

# %% Process specifications

# Specifications dependent on lipid cane flow rate
_enzyme_mass = enzyme.mass[[9, 12]]
_CaO_Water_mass = lime.mass[[7, 12]]
_H3PO4_Water_mass = H3PO4.mass[[1, 12]]
last_lipidcane_massnet = int(Lipid_cane.massnet)
def correct_flows():
    global last_lipidcane_massnet
    massnet = Lipid_cane.massnet
    if int(massnet) != last_lipidcane_massnet:
        # correct enzyme, lime, phosphoric acid, and imbibition water
        _enzyme_mass[:] = 0.003 * massnet * np.array([0.1, 0.9])
        _CaO_Water_mass[:] = 0.001 * massnet * np.array([0.046, 0.954])
        _H3PO4_Water_mass[:] = 0.00025 * massnet
        imbibition_water_mass.value = 0.25* massnet
        last_lipidcane_massnet = int(massnet)

# Specifications within a system
def correct_lipid_wash_water():
    lipid_wash.mol[12] = P49.outs[0].mol[-2]*100/11

solids_index = Stream.indices(['Ash', 'CaO', 'Cellulose', 'Hemicellulose', 'Lignin'])
def correct_wash_water():
    solids = solidsmol[solids_index].sum()
    wash_water.mol[12] = 0.0574*solids

imbibition_water_mass = imbibition_water.mass.item(12)


# %% Pretreatment system set-up

(P137-0, P21-0)-Mill-1-P56-0-P1
(P56-1, imbibition_water)-P21

crushing_mill_recycle = System('crushing_mill_recycle',
                               network=(Mill, P56, P21),
                               recycle=P21-0)

(F7-0, F8-0)-P4-P7
(P7-0, polymer)-P68-P12
(P12-1, wash_water)-P14-1-F8
clarification_recycle = System('clarification_recycle',
                               network=(P4, P7, P68, P12, P14, F8),
                               recycle=P14-1)

(F3-0, enzyme)-P137
P1-0-P3
(P3-0, H3PO4)-P9-F5
(F5-0, lime-F6-0)-P5-F7
P12-0-P10-0-P49
(P49-0, lipid_wash)-P75-P76-0-P69
(P14-0, S254)-P46
P10-1-P50

pretreatment_sys = System('pretreatment_sys',
                          network=(F1, F2, F3,
                                   correct_flows, P137,
                                   crushing_mill_recycle,
                                   F4, P1, P3, P9,
                                   F5, F6, P5, F7,
                                   correct_wash_water,
                                   clarification_recycle,
                                   P10, P49,
                                   correct_lipid_wash_water,
                                   P75, P76,
                                   P69, P46, P50))

solidsmol = F7.outs[0].mol
