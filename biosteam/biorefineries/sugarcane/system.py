# -*- coding: utf-8 -*-
"""
Created on Tue May 28 11:11:38 2019

@author: yoelr
"""
import numpy as np
from biosteam import System, Stream, Species, TEA, find
from biosteam.units import Mixer, EnzymeTreatment, CrushingMill, \
    HXutility, RVF, VibratingScreen, \
    MagneticSeparator, Clarifier, MixTank, \
    Shredder, ConveyingBelt, \
    Pump, StorageTank, Splitter, HXprocess, \
    MultiEffectEvaporator, Fermentation, \
    Distillation, SolidsCentrifuge, MolecularSieve, \
    VentScrubber
from biosteam.units.facilities import CoolingTower, \
                                      BoilerTurbogenerator, \
                                      ChilledWaterPackage, \
                                      ProcessWaterCenter
from biosteam.biorefineries.lipidcane.tea import LipidcaneTEA
from biosteam.biorefineries.sugarcane.species import sugarcane_species
from biosteam.biorefineries.sugarcane.process_settings import price
# import warnings
# warnings.filterwarnings('ignore')

# %% Streams
Stream.species = sugarcane_species
sugar_cane_composition = {'Glucose': 0.0120811,
                          'Lignin': 0.0327653,
                          'Solids': 0.015,
                          'Sucrose': 0.136919,
                          'Ash': 0.006,
                          'Cellulose': 0.0611531,
                          'Hemicellulose': 0.036082,
                          'Water': 0.7}
Sugar_cane = Stream('sugar_cane', units='kg/hr',
                    **{i:j*333334 for i,j in sugar_cane_composition.items()},
                    price=price['Sugar cane'])

enzyme = Stream('enzyme', Cellulose=100, Water=900, units='kg/hr',
                price=price['Protease'])

imbibition_water = Stream('imbibition_water',
                          Water=87023.35,
                          T = 338.15, units='kg/hr')

H3PO4 = Stream('H3PO4', H3PO4=74.23, Water=13.10, units='kg/hr',
               price=price['H3PO4'])  # to P9

lime = Stream('lime', CaO=333.00, Water=2200.00, units='kg/hr',
              price=price['Lime'])  # to P5

polymer = Stream('polymer', Flocculant=0.83, units='kg/hr',
                 price=price['Polymer'])  # to P68

wash_water = Stream('water_1', Water=16770, units='kg/hr')  # to P14
wash_water.T = 363.15

S254 = Stream('S254', Ash=1, units='kg/hr')  # to P46

# %% Units

# Feed the shredder
F1 = ConveyingBelt('F1', ins=Sugar_cane)
F1.cost_items['Conveying belt'].ub = 2500

# Separate metals
F2 = MagneticSeparator('F2', ins=F1.outs)

# Shred fresh cane
F3 = Shredder('F3', ins=F2.outs)

# Hydrolyze starch
P137 = EnzymeTreatment('P137', 'S102', T=323.15)  # T=50

# Finely crush cane
Mill = CrushingMill('Mill', outs=('S107', ''),
                    split=(0.92, 0.92, 0.04, 0.92, 0.92, 0.04, 1),
                    order=('Ash', 'Cellulose', 'Glucose', 'Hemicellulose',
                           'Lignin', 'Sucrose', 'Solids'),
                    moisture_content=0.5)

# Convey out bagasse
F4 = ConveyingBelt('F4', ins=Mill.outs[0], outs='Bagasse')

# Mix in water
P21 = Mixer('P21', 'S171')

# Screen out fibers
P56 = VibratingScreen('P56', ('S113', 'S105'),
                      split=(0.35, 0.35, 0.88, 0.35, 0.35, 0.88, 0.88),
                      order=('Ash', 'Cellulose', 'Glucose', 'Hemicellulose',
                             'Lignin', 'Sucrose', 'Water'))
# Store juice before treatment
P1 = StorageTank('P1')
P1.tau = 12

# Heat up before adding acid
P3 = HXutility('P3', 'S108', T=343.15)

# Mix in acid
P9 = MixTank('P9', 'S109')

# Pump acid solution
F5 = Pump('F5')

# Mix lime solution
F6 = MixTank('F6')
F6.tau = 1
F7 = Pump('F7')

# Blend acid solution with lime
P5 = MixTank('P5', 'S118')

# Mix recycle
P4 = Mixer('P4', 'S122')

# Heat before adding flocculant
P7 = HXutility('P7', 'S110', T=372.15)

# Mix in flocculant
P68 = MixTank('P68', 'S111')
P68.tau = 1/4
P68.electricity_rate = 1

# Separate residual solids
P12 = Clarifier('P12', ('S123', 'S104'), 
                split=(0.522, 0.522, 0.522, 0.522, 0.522),
                order=('Flocculant', 'Glucose', 'H3PO4', 'Sucrose', 'Water'))

# Remove solids as filter cake
P14 = RVF('P14', ('S124', 'S116'), 
          moisture_content=0.80,
          split=(0.85, 0.85, 0.85, 0.01, 0.85, 0.85, 0.01),
          order=('Ash', 'CaO', 'Cellulose', 'Glucose',
                 'Hemicellulose', 'Lignin', 'Sucrose'))
F8 = Pump('F8')

# Get filter cake
P46 = Mixer('P46', 'filter_cake')

# Screen out small fibers from sugar stream
P50 = VibratingScreen('P50', ('sugar', 'fiber_fines'),
                      split=(0.998, 0.998, 0.998, 1),
                      order=('Glucose', 'Sucrose', 'Water', 'Flocculant'))
Sugar = P50-0
P50.mesh_opening = 2


# %% Process specifications

# Specifications dependent on cane flow rate
_enzyme_mass = enzyme.mass[enzyme.indices(['Cellulose', 'Water'])]
_CaO_Water_mass = lime.mass[lime.indices(['CaO', 'Water'])]
_H3PO4_Water_mass = H3PO4.mass[H3PO4.indices(['H3PO4', 'Water'])]
last_sugarcane_massnet = int(Sugar_cane.massnet)
def correct_flows():
    global last_sugarcane_massnet
    massnet = Sugar_cane.massnet
    if int(massnet) != last_sugarcane_massnet:
        # correct enzyme, lime, phosphoric acid, and imbibition water
        _enzyme_mass[:] = 0.003 * massnet * np.array([0.1, 0.9])
        _CaO_Water_mass[:] = 0.001 * massnet * np.array([0.046, 0.954])
        _H3PO4_Water_mass[:] = 0.00025 * massnet
        imbibition_water_mass.value = 0.25* massnet
        last_sugarcane_massnet = int(massnet)

def correct_wash_water():
    solids = solidsmass.sum()
    wash_water.mol[12] = 0.0574*solids

imbibition_water_mass = imbibition_water.mass.item(7)


# %% Pretreatment system set-up

(P137-0, P21-0)-Mill-1-P56-0-P1
(P56-1, imbibition_water)-P21

crushing_mill_recycle = System('crushing_mill_recycle',
                               network=(Mill, P56, P21),
                               recycle=P21-0)

(F7-0, F8-0)-P4-P7
(P7-0, polymer)-P68-P12-0-P50
(P12-1, wash_water)-P14-1-F8
clarification_recycle = System('clarification_recycle',
                               network=(P4, P7, P68, P12, P14, F8),
                               recycle=P14-1)

(F3-0, enzyme)-P137
P1-0-P3
(P3-0, H3PO4)-P9-F5
(F5-0, lime-F6-0)-P5-F7
(P14-0, S254)-P46

pretreatment_sys = System('pretreatment_sys',
                          network=(F1, F2, F3,
                                   correct_flows, P137,
                                   crushing_mill_recycle,
                                   F4, P1, P3, P9,
                                   F5, F6, P5, F7,
                                   correct_wash_water,
                                   clarification_recycle,
                                   P46, P50))
solidsmass = F7.outs[0].mass[Stream.indices(
        ['Ash', 'CaO', 'Cellulose', 'Hemicellulose', 'Lignin'])]

# %% Ethanol function

Ethanol_MW = sugarcane_species.Ethanol.MW
Water_MW = sugarcane_species.Water.MW

def Ethanol_molfrac(e: 'Ethanol mass fraction'):
    """Return ethanol mol fraction in a ethanol water mixture"""
    return e/Ethanol_MW / (e/Ethanol_MW + (1-e)/Water_MW)


# %% Input Streams

# Fresh water
S134 = Stream('S134', Water=5000, units='kg/hr')

# Gasoline
denature = Stream('denaturant', Octane=230.69,
                  units='kg/hr', price=price['Gasoline'])

# Yeast
S144 = Stream('yeast', Water=24700, DryYeast=10300,
              units='kg/hr')

# Ethanol product
ethanol = Stream('ethanol', price=price['Ethanol'])


# %% Units

# Split sugar solution
P16 = Splitter('P16',
               split=0.265)

# Concentrate sugars
P15 = MultiEffectEvaporator('P15',
                            component='Water', # component being evaporated
                            P=(101325, 73581, 50892, 32777, 20000),
                            V=0.95248) # fraction evaporated
P15.components['condenser'].U = 1.85
# Note: value of steam ~ 6.86 for the following 
# (101325, 73580.467, 50891.17, 32777.406, 19999.925, 11331.5),

# Mix sugar solutions
P17 = Mixer('P17')

# Cool for fermentation
P51 = HXutility('P51', T=295.15)

# Ethanol Production
P24 = Fermentation('P24', outs=('', 'CO2'),
                   tau=10, efficiency=0.90, N=8) 
T1 = StorageTank('T1')
T1.tau = 4

scrubber = VentScrubber('Scrubber', ins=(S134, P24-0), gas=('CO2',))

# Separate 99% of yeast
P19 = SolidsCentrifuge('P19', outs=('', 'recycle_yeast'),
                       split=(1, 0.99999, 1, 0.96, 0.01),
                       order=('Ethanol', 'Glucose', 'H3PO4', 
                              'Water', 'DryYeast'),
                       solids=('DryYeast',))

# Mix in Water
P23 = Mixer('P23', 'S147')
Q1 = Pump('Q1')

# Heat up before beer column
# Exchange heat with stillage
P32 = HXprocess('P32', outs=('', 'stillage'),
                fluid_type='ss', U=1.28)

# Beer column
xbot = Ethanol_molfrac(0.00001)
ytop = Ethanol_molfrac(0.574)
P25 = Distillation('P25', P=101325,
                   y_top=ytop, x_bot=xbot, k=1.20,
                   LHK=('Ethanol', 'Water'))
P25.tray_material = 'Stainless steel 316'
P25.vessel_material = 'Stainless steel 316'
P25._boiler.U = 1.85
Q2 = Pump('Q2')

# Mix ethanol Recycle (Set-up)
P28 = Mixer('P28', 'S150')

ytop = Ethanol_molfrac(0.9061726)
P30 = Distillation('P30', P=101325,
                    y_top=ytop, x_bot=xbot, k=1.20,
                    LHK=('Ethanol', 'Water'))
P30.tray_material = 'Stainless steel 316'
P30.vessel_material = 'Stainless steel 316'
P30.is_divided = True
P30._boiler.U = 1.85
Q3 = Pump('Q3')

# Superheat vapor for mol sieve
P34 = HXutility('P34', T=115+273.15, V=1)

# Molecular sieve
P33 = MolecularSieve('P33',
                     split=(2165.14/13356.04, 1280.06/1383.85),
                     order=('Ethanol', 'Water'))

# Condense ethanol product
P41 = HXutility('P41', 'S149', V=0, T=350.)
T2 = StorageTank('T2')
T2.line = 'Ethanol storage'
T2.tau = 6*24
Q4 = Pump('Q4')

# Storage for gasoline
T3 = StorageTank('T3')
T3.tau = 6*24
Q5 = Pump('Q5')

# Denatured ethanol product
T4 = MixTank('T4', outs=ethanol)
T4.tau = 0.10

# Add denaturant
P39 = Mixer('P39')

# Recycle water to Process Condensate Tank
P35_0 = Mixer('P34', outs='Water_4')

# Yeast mixing
T5 = MixTank('T5')
T5.tau = 0.1
S144-T5

# Multi-effect evaporator pumps
Q6 = Pump('Q6')


# %% Set up EtOH system

P50-0-P16-1-P15-0-Q6
(P16-0, Q6-0)-P17-P51
(P51-0, S144-T5-0)-P24-1-T1-0-P19
(P19-0, scrubber-1)-P23-Q1
(Q1-0, Q2-0)-P32-0-P25-1-Q2
EtOH_start_network = (P16, P15, Q6, P17, P51, T5, P24, T1,
                      P19, scrubber, P23, Q1, P32, P25, Q2, P32)

(P25-0, P33-0)-P28-0-P30-0-P34-P33
P30-1-Q3
purification_recycle = System('purification_recycle',
                              network=(P28, P30, P34, P33),
                              recycle=P28-0)

gas_index = sugarcane_species.IDs.index('Octane')
def adjust_denaturant():
    denature.mol[gas_index] = 0.021*Q4.outs[0].massnet/114.232
    
P33-1-P41-0-T2-0-Q4
denature-T3-Q5
(Q5-0, Q4-0)-P39-T4
EtOH_end_network=(Q3, P41, T2, Q4, adjust_denaturant, T3, Q5, P39, T4)

(Q3-0, P15-1)-P35_0
EtOH_process_water_network = (P35_0,)

ethanol_sys = System('ethanol_sys',
                     network=(EtOH_start_network
                            + (purification_recycle,)
                            + EtOH_end_network
                            + EtOH_process_water_network))


# %% Facilities
emission = Stream('Emission')
water = Species('Water',)
stream = find.stream

BT = BoilerTurbogenerator('BT',
                          ins=F4-0,
                          outs=emission,
                          boiler_efficiency=0.80,
                          turbogenerator_efficiency=0.85)

Stream.species = water
CT = CoolingTower('CT')
process_water_streams = (stream.water_1, stream.imbibition_water, stream.S134)
recycle_water_streams = ()
makeup_water = Stream('makeup_water', price=0.000254)
process_water = Stream('process_water')
pws_indices = [i.index('Water') for i in process_water_streams]
pws_flow_index_pairs = tuple(zip([i._mol for i in process_water_streams], pws_indices))
def update_water():
    process_water._mol[0] = sum([i[j] for i,j in pws_flow_index_pairs])

CWP = ChilledWaterPackage('CWP')
PWC = ProcessWaterCenter('PWC',
                         ins=('', makeup_water),
                         outs=(process_water,))

# %% Set up system
sugarcane_sys = System('sugar_cane_system',
                       network=(pretreatment_sys.network
                              + ethanol_sys.network),
                       facilities=(CWP, BT, CT, update_water, PWC))
units = sugarcane_sys.units


# %% Perform TEA
sugarcane_sys.reset_names()
sugarcane_tea = LipidcaneTEA(system=sugarcane_sys, IRR=0.15, duration=(2018, 2038),
                             depreciation='MACRS7', income_tax=0.35,
                             operating_days=200, lang_factor=3,
                             construction_schedule=(0.4, 0.6), WC_over_FCI=0.05,
                             labor_cost=2.5e6, fringe_benefits=0.4,
                             property_tax=0.001, property_insurance=0.005,
                             supplies=0.20, maintenance=0.01, administration=0.005)
sugarcane_sys.simulate()
sugarcane_tea.IRR = sugarcane_tea.solve_IRR()


