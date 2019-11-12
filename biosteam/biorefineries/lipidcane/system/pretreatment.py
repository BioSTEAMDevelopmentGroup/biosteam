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

__all__ = ('pretreatment_sys', 'lipid_cane', 'lipidcane', 'area_100', 'area_200')

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
lipidcane = lipid_cane = Stream('lipid_cane', f1, psp1, units='kg/hr',
                                price=price['Lipid cane'])

enzyme = Stream('enzyme', Cellulose=100, Water=900, units='kg/hr',
                price=price['Protease'])

imbibition_water = Stream('imbibition_water',
                          Water=87023.35,
                          T = 338.15, units='kg/hr')

H3PO4 = Stream('H3PO4', H3PO4=74.23, Water=13.10, units='kg/hr',
               price=price['H3PO4'])  # to T203

lime = Stream('lime', CaO=333.00, Water=2200.00, units='kg/hr',
              price=price['Lime'])  # to P5

polymer = Stream('polymer', Flocculant=0.83, units='kg/hr',
                 price=price['Polymer'])  # to T205

rvf_wash_water = Stream('rvf_wash_water',
                        Water=16770, units='kg/hr',
                        T=363.15)  # to C202

oil_wash_water = Stream('oil_wash_water',
                        Water=1350, units='kg/hr',
                        T=358.15)  # to T207

# %% Units

Stream.default_ID = 'd'
Stream.default_ID_number = 0
# Stream.default_ID_number = 100

# Feed the shredder
U101 = ConveyingBelt('U101', ins=lipid_cane)
U101.cost_items['Conveying belt'].ub = 2500

# Separate metals
U102 = MagneticSeparator('U102', ins=U101.outs)

# Shredded cane
U103 = Shredder('U103', ins=U102.outs)

# Stream.default_ID_number = 200

# Hydrolyze starch
T201 = EnzymeTreatment('T201', T=323.15)  # T=50

# Finely crush lipid cane
U201 = CrushingMill('U201',
                    split=(0.92, 0.92, 0.04, 0.92, 0.92, 0.04, 0.1, 1),
                    order=('Ash', 'Cellulose', 'Glucose', 'Hemicellulose',
                           'Lignin', 'Sucrose', 'Lipid', 'Solids'),
                    moisture_content=0.5)

# Convey out bagasse
U202 = ConveyingBelt('U202', ins=U201.outs[0], outs='Bagasse')

# Mix in water
M201 = Mixer('M201')

# Screen out fibers
S201 = VibratingScreen('S201',
                       split=(0.35, 0.35, 0.88, 0.35,
                              0.35, 0.88, 0, 0.88, 0.88),
                       order=psp1)

# Store juice before treatment
T202 = StorageTank('T202')
T202.tau = 12

# Heat up before adding acid
H201 = HXutility('H201', T=343.15)

# Mix in acid
T203 = MixTank('T203')

# Pump acid solution
P201 = Pump('P201')

# Mix lime solution
T204 = MixTank('T204')
T204.tau = 1
P202 = Pump('P202')

# Blend acid lipid solution with lime
T205 = MixTank('T205')

# Mix recycle
M202 = Mixer('M202')

# Heat before adding flocculant
H202 = HXutility('H202', T=372.15)

# Mix in flocculant
T206 = MixTank('T206')
T206.tau = 1/4

# Separate residual solids
C201 = Clarifier('C201',
                split=(0, 0, 0, 0.522, 0.522, 0, 0,
                       0.98, 0.522, 0.522, 0.522),
                order=psp2)

# Remove solids as filter cake
C202 = RVF('C202', 
           outs=('filte_cake', ''),
           moisture_content=0.80,
           split=(0.85, 0.85, 0.85, 0.01, 0.85, 0.85, 0.01),
           order=('Ash', 'CaO', 'Cellulose', 'Glucose',
                  'Hemicellulose', 'Lignin', 'Sucrose'))
P203 = Pump('P203')

# Separate oil and sugar
T207 = MixTank('T207', outs=('', ''))
split = np.zeros(len(pretreatment_species), float)
index = pretreatment_species.indices(('Lipid', 'Water'))
split[index] = (1, 0.0001)
T207._split = split
T207._run = lambda : Splitter._run(T207)
del split, index

# Cool the oil
H203 = HXutility('H203', T=343.15)

# Screen out small fibers from sugar stream
S202 = VibratingScreen('S202', outs=('', 'fiber_fines'),
                      split=1-np.array((0, 0, 0, 1, 0.002, 0, 0,0, 0, 0.002, 0.002)),
                      order=psp2)
sugar = S202-0
S202.mesh_opening = 2

# Add distilled water to wash lipid
T208 = MixTank('T208')
T208.tau = 2

# Centrifuge out water
C203 = SplitCentrifuge_LLE('C203',
                           split=(0.99, 0.01),
                           order=('Lipid', 'Water'))

# Vacume out water
F201 = SplitFlash('F201', T=347.15, P=2026.5,
                 split=(0.0001, 0.999), order=('Lipid', 'Water'))
lipid = F201.outs[1]

# %% Process specifications

# Specifications dependent on lipid cane flow rate
_enzyme_mass = enzyme.mass[[9, 12]]
_CaO_Water_mass = lime.mass[[7, 12]]
_H3PO4_Water_mass = H3PO4.mass[[1, 12]]
last_lipidcane_massnet = int(lipid_cane.massnet)
def correct_flows():
    global last_lipidcane_massnet
    massnet = lipid_cane.massnet
    if int(massnet) != last_lipidcane_massnet:
        # correct enzyme, lime, phosphoric acid, and imbibition water
        _enzyme_mass[:] = 0.003 * massnet * np.array([0.1, 0.9])
        _CaO_Water_mass[:] = 0.001 * massnet * np.array([0.046, 0.954])
        _H3PO4_Water_mass[:] = 0.00025 * massnet
        imbibition_water_mass.value = 0.25* massnet
        last_lipidcane_massnet = int(massnet)

# Specifications within a system
def correct_lipid_wash_water():
    oil_wash_water.mol[12] = H202.outs[0].mol[-2]*100/11

solids_index = Stream.indices(['Ash', 'CaO', 'Cellulose', 'Hemicellulose', 'Lignin'])
def correct_wash_water():
    solids = solidsmol[solids_index].sum()
    rvf_wash_water.mol[12] = 0.0574*solids

imbibition_water_mass = imbibition_water.mass.item(12)


# %% Pretreatment system set-up

(U103-0, enzyme)-T201
(T201-0, M201-0)-U201-1-S201-0-T202
(S201-1, imbibition_water)-M201
crushing_mill_recycle_sys = System('crushing_mill_recycle_sys',
                               network=(U201, S201, M201),
                               recycle=M201-0)

T202-0-H201
(H201-0, H3PO4)-T203-P201
(P201-0, lime-T204-0)-T205-P202
(P202-0, P203-0)-M202-H202
(H202-0, polymer)-T206-C201
(C201-1, rvf_wash_water)-C202-1-P203
clarification_recycle_sys = System('clarification_recycle_sys',
                                   network=(M202, H202, T206, C201, C202, P203),
                                   recycle=C202-1)

C201-0-T207-0-H203
(H203-0, oil_wash_water)-T208-C203-0-F201
T207-1-S202

pretreatment_sys = System('pretreatment_sys',
                          network=(U101, U102, U103,
                                   correct_flows, T201,
                                   crushing_mill_recycle_sys,
                                   U202, T202, H201, T203,
                                   P201, T204, T205, P202,
                                   correct_wash_water,
                                   clarification_recycle_sys,
                                   T207, H203, S202,
                                   correct_lipid_wash_water,
                                   T208, C203, F201,))

solidsmol = P202.outs[0].mol

area_100 = System('area_100', network=(U101, U102, U103))
units = pretreatment_sys.units.copy()
for i in area_100.network: units.discard(i)
area_200_network = sorted(units, key=lambda x: x.ID)
area_200 = System('area_200', network=area_200_network)
