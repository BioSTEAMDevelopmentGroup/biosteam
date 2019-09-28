#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 18 16:17:00 2017

The biodiesel production section for the baseline lipid cane biorefinery is defined here as System objects. The systems include all streams and units starting from transesterification to purification of biodiesel and glycerol products.

@author: yoel
"""
from biosteam import System, Stream, MixedStream
from biosteam.units import Transesterification, \
    SplitFlash, MixTank, Distillation, RatioCentrifuge_LLE, \
    InvSplitter, Pump, StorageTank, \
    HXutility, HXprocess, Flash, \
    SplitCentrifuge_LLE, MassBalance
from biosteam.biorefineries.lipidcane.species import biodiesel_species
from biosteam.biorefineries.lipidcane.process_settings import price
from numpy import array


__all__ = ('area_400',)

# %% Stream Settings

# Set species library
Stream.species = biodiesel_species
# Stream.default_ID_number = 400


# %% Fresh Streams

# Fresh degummed oil
oil = Stream('oil', Lipid=8853.49, Water=1.00e-02,
             units='kg/hr', T=347.15)

# Fresh methanol
methanol = Stream('methanol', Methanol=1,
                  price=price['Methanol'])

# Catalyst
catalyst = Stream('catalyst', NaOCH3=0.25,
                  Methanol=0.75, units='kg/hr',
                  price=price['NaOCH3'])
                  # price=0.25*price['NaOCH3'] + 0.75*Methanol.price)

# Water to remove glycerol
biodiesel_wash_water = Stream('biodiesel_wash_water', Water=13.6, T=273.15+60, 
               price=price['Water'])

# Acid to neutralize catalyst after second centrifuge
HCl1 = Stream('HCl1', HCl=0.21, Water=0.79,
              price=price['HCl'])
              # price=0.21*price['HCl'] + 0.79*Water.price) # 35% HCl by mass

# Acid to remove soaps after first centrifuge
HCl2 = Stream('HCl2', HCl=0.21, Water=0.79,
              price=HCl1.price)

# Base to neutralize acid before distillation
NaOH = Stream('NaOH', NaOH=1, price=price['NaOH'])

# Products
crude_glycerol = Stream('crude_glycerol',
                        price=price['Crude glycerol'])
biodiesel = Stream('biodiesel',
                   price=price['Biodiesel'])

# Waste
waste = Stream('waste', price=price['Waste'])

# %% Units

"""Biodiesel Transesterification Section"""

# Aparently reactors are adiabatic, meoh coming in at 40C, lipid at 60C

# From USDA Biodiesel model
x_cat = 1.05e-05 # Catalyst molar fraction in methanol feed

# Mix Recycle Methanol and Fresh Methanol
T401 = StorageTank('T401')
P401 = Pump('P401')

# Storage Tank for Catalyst
T402 = StorageTank('T402')
P402 = Pump('P402')

# Tank for oil
T403 = StorageTank('T403')
T403.tau = 4
P403 = Pump('P403')

# Mix Methanol and Catalyst stream
T404 = MixTank('T404')
P404 = Pump('P404')

# Mass Balance for Methanol, Recycle Methanol, and Catalyst stream
B401 = MassBalance('B401', species=('Methanol', 'NaOCH3'), streams=(0, 1))

# Split Methanol/Catalyst to reactors
S401 = InvSplitter('S401')

# First Reactor
R401 = Transesterification('R401', efficiency=0.90, methanol2lipid=6, T=333.15,
                         catalyst_molfrac=x_cat)

# Centrifuge to remove glycerol
C401 = SplitCentrifuge_LLE('C401',
                         split=(0.99,  # Lipid
                                0.40,  # Methanol
                                0.06,  # Glycerol
                                0.999, # Biodiesel
                                0.40,  # Water
                                0,     # NaOH
                                0,     # HCl
                                0.40)) # NaOCH3

P405 = Pump('P405')

# Second Reactor
R402 = Transesterification('R402', efficiency=0.90, methanol2lipid=6, T=333.15,
                         catalyst_molfrac=x_cat) 

# Centrifuge to remove glycerol
C402 = SplitCentrifuge_LLE('C402',
                         split=(0.90,  # Lipid
                                0.10,  # Methanol
                                0.05,  # Glycerol
                                0.999, # Biodiesel
                                0.10,  # Water
                                0,     # NaOH
                                0,     # HCl
                                0.10)) # NaOCH3

# Acids and bases per catalyst by mol
catalyst_index = oil.index('NaOCH3')
k1 = 0.323/1.5; k2 = 1.060/1.5; k3 = 0.04505/1.5;
catalyst_mol = T404.outs[0].mol[7]
def adjust_acid_and_base():
    global catalyst_mol
    # Adjust according to USDA biodiesel model
    new = T404.outs[0].mol[7]
    if new != catalyst_mol:
        catalyst_mol = new
        NaOH._mol[5] = k1 * new
        HCl1._mol[6] = k2 * new
        HCl2._mol[6] = k3 * new

"""Biodiesel Purification Section"""

# Water wash
T405 = MixTank('T405')
T405.tau = 0.25
P406 = Pump('P406')

# Centrifuge out water
C403 = RatioCentrifuge_LLE('C403',
                         Kspecies=('Methanol', 'Glycerol'),
                         Ks=array((0.382, 0.183)),
                         top_solvents=('Biodiesel',),
                         top_split=(0.999,),
                         bot_solvents=('Water', 'Lipid', 'NaOH', 'HCl'),
                         bot_split=(0.999, 1, 1, 1))

# Vacuum dry biodiesel
# Consider doing this by split, and keeping this unit constant
# 290 Pa, 324 K according to USDA Biodiesel Model
F401 = SplitFlash('F401',
                order=('Water', 'Methanol', 'Biodiesel'),
                split=(0.9999, 0.9999, 0.00001),
                P=2026.5, T=331.5, Q=0)
F401.line = 'Vacuum dryer'
F401.material = 'Stainless steel 316'
P407 = Pump('P407')

"""Glycerol Purification Section"""

# Condense vacuumed methanol to recycle
H401 = HXutility('H401', V=0, T=295)
P408 = Pump('P408')

# Mix recycled streams and HCl
T406 = MixTank('T406')
P409 = Pump('P409')

# Centrifuge out waste fat
# assume all the lipid, free_lipid and biodiesel is washed out
C404 = SplitCentrifuge_LLE('C404', outs=('', waste),
                         order=('Methanol', 'Glycerol', 'Water'),
                         split=(0.999, 0.999, 0.999))

# Add and mix NaOH
T407 = MixTank('T407')
P410 = Pump('P410')

# Methanol/Water distillation column
D401 = Distillation('D401',
                  LHK=('Methanol', 'Water'), P=101325,
                  y_top=0.99999, x_bot=0.0001, k=2.5)
D401.is_divided = True
D401.tray_material = 'Stainless steel 316'
D401.vessel_material = 'Stainless steel 316'

# Condense recycle methanol
H402 = HXutility('H402', V=0, T=315)
P411 = Pump('P411')

# Glycerol/Water flash (not a distillation column)
w = 0.20/Stream.species.Water.MW
g = 0.80/Stream.species.Glycerol.MW
x = w/(w+g)
F402 = Flash('F402',
           species_IDs=('Water', 'Glycerol'),
           P=101325,
           x=(x, 1-x),
           HNK=('Biodiesel',),
           LNK=('Methanol',))
F402.HasGlycolGroups = True
F402.vessel_material = 'Stainless steel 316'

# Condense water to recycle
H403 = HXprocess('H403', fluid_type='ll', outs=(MixedStream(), ''),
                     species_IDs=('Methanol', 'Water'), 
                     HNK=('Glycerol', 'Biodiesel', 'NaOH', 'HCl'))
    
# Use heat from glycerol
H404 = HXprocess('H404', fluid_type='ls',
                     species_IDs=('Methanol', 'Water'), 
                     HNK=('Glycerol', 'Biodiesel', 'NaOH', 'HCl'))
run_bot = H404._run

# Storage tank for glycerol
T408 = StorageTank('T408', outs=crude_glycerol)
H404-1-T408

# Storage tank for biodiesel
T409 = StorageTank('T409', outs=biodiesel)
F401-1-P407-0-T409

# %% Set up systems

# Biodiesel Transesterification Section
oil-T403-P403
(P403-0, S401-0)-R401-0-C401
(C401-0, S401-1)-R402-0-C402
transesterification_network = (T403, P403, R401, C401, P405, R402, C402,
                               adjust_acid_and_base)

"""
Specs for product https://www.afdc.energy.gov/fuels/biodiesel_specifications.html
minimum spec requirements:
 0.050 wt % water (lower to 0.045)
 0.2 wt % meoh (lower to 0.15)
 0.02 wt % glycerol (lower to 0.015)
 0.4 wt % free lipids (lower to 0.35)
"""

# Find Water Flow
(C402-0, H403-0, biodiesel_wash_water, HCl1)-T405-P406-C403
def adjust_water_flow():
    H403.outs[0].liquid_mol[4] = 800*C402.outs[0]._mol[2]

# Glycerol recycle and purification section
C403-0-F401
F401-0-H401-P408
C401-1-P405
(P405-0, C402-1, C403-1, P408-0, HCl2)-T406-P409-C404
(C404-0, NaOH)-T407-P410
(H403-1, F402-1)-H404
H404-0-D401-1-F402
(F402-0, P410-0)-H403
glycerol_recycle_sys = System('glycerol_recycle_sys',
                              network=(adjust_water_flow, 
                                       T405, P406, C403, F401, H401,
                                       P408, P407, T406, P409, C404,
                                       T407, P410, H403, H404, D401,
                                       F402, H403, H404),
                              recycle=F402-0)                           
# Find Fresh Methanol Flow
D401-0-H402-P411    # Recycle methanol
methanol-T401-P401  # Mix fresh and recycled methanol
catalyst-T402-P402  # Fresh catalyst
(P411-0, P401-0, P402-0)-T404-P404-S401  # Mix Catalyst with Methanol
meoh_network = (H402, P411, T401, P401, T402, P402, T404, P404, S401)

# Set connections for mass balance proxy
(catalyst, methanol, D401-0)-B401
B401**(1**R401, 1**R402)

# Complete System
area_400 = System('area_400',
                  network=transesterification_network
                          + (glycerol_recycle_sys, B401)
                          + meoh_network
                          + (T408, T409))

# Initial conditions
index = oil.indices(['Methanol', 'Glycerol', 'Water'])
F402.outs[0].mol[index] = (3.483e-02, 8.284e-03, 2.094e+01)
F402.outs[0].T = 374.986
H403.outs[0].T = 374.986
H403.outs[0].liquid_mol[index] = (3.483e-02, 8.284e-03, 2.094e+01)

