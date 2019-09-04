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


__all__ = ('biodiesel_sys',)

# %% Stream Settings

# Set species library
Stream.species = biodiesel_species


# %% Fresh Streams

# Fresh degummed oil
Oil = Stream('Oil', Lipid=8853.49, Water=1.00e-02,
             units='kg/hr', T=347.15)

# Fresh Methanol
Methanol = Stream('Methanol', Methanol=1,
                  price=price['Methanol'])

# Catalyst
Catalyst = Stream('Catalyst', NaOCH3=0.25,
                  Methanol=0.75, units='kg/hr',
                  price=price['NaOCH3'])
                  # price=0.25*price['NaOCH3'] + 0.75*Methanol.price)

# Water to remove glycerol
Water = Stream('Water 3', Water=13.6, T=273.15+60, 
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
Crude_glycerol = Stream('Crude glycerol',
                        price=price['Crude glycerol'])
Biodiesel = Stream('Biodiesel',
                   price=price['Biodiesel'])

# Waste
Waste = Stream('Waste', price=price['Waste'])

# %% Units

"""Biodiesel Transesterification Section"""

# Aparently reactors are adiabatic, meoh coming in at 40C, lipid at 60C

# From USDA Biodiesel model
x_cat = 1.05e-05 # Catalyst molar fraction in methanol feed

# Mix Recycle Methanol and Fresh Methanol
T1 = StorageTank('T1')
P1 = Pump('P1')

# Storage Tank for Catalyst
T2 = StorageTank('T2')
P2 = Pump('P2')

# Tank for Oil
T3 = StorageTank('T3')
T3.tau = 4
P3 = Pump('P3')

# Mix Methanol and Catalyst stream
M1 = MixTank('M1')
P0 = Pump('P0')

# Mass Balance for Methanol, Recycle Methanol, and Catalyst stream
MB1 = MassBalance('MB1', species=('Methanol', 'NaOCH3'), streams=(0, 1))

# Split Methanol/Catalyst to reactors
S1 = InvSplitter('S1')

# First Reactor
R1 = Transesterification('R1', efficiency=0.90, methanol2lipid=6, T=333.15,
                         catalyst_molfrac=x_cat)

# Centrifuge to remove glycerol
C1 = SplitCentrifuge_LLE('C1',
                         split=(0.99,  # Lipid
                                0.40,  # Methanol
                                0.06,  # Glycerol
                                0.999, # Biodiesel
                                0.40,  # Water
                                0,     # NaOH
                                0,     # HCl
                                0.40)) # NaOCH3

P4 = Pump('P4')

# Second Reactor
R2 = Transesterification('R2', efficiency=0.90, methanol2lipid=6, T=333.15,
                         catalyst_molfrac=x_cat) 

# Centrifuge to remove glycerol
C2 = SplitCentrifuge_LLE('C2',
                         split=(0.90,  # Lipid
                                0.10,  # Methanol
                                0.05,  # Glycerol
                                0.999, # Biodiesel
                                0.10,  # Water
                                0,     # NaOH
                                0,     # HCl
                                0.10)) # NaOCH3

# Acids and bases per catalyst by mol
catalyst_index = Oil.index('NaOCH3')
k1 = 0.323/1.5; k2 = 1.060/1.5; k3 = 0.04505/1.5;
catalyst_mol = M1.outs[0].mol[7]
def adjust_acid_and_base():
    global catalyst_mol
    # Adjust according to USDA biodiesel model
    new = M1.outs[0].mol[7]
    if new != catalyst_mol:
        catalyst_mol = new
        NaOH._mol[5] = k1 * new
        HCl1._mol[6] = k2 * new
        HCl2._mol[6] = k3 * new

"""Biodiesel Purification Section"""

# Water wash
W1 = MixTank('W1')
W1.tau = 0.25
X0 = Pump('X0')

# Centrifuge out water
C3 = RatioCentrifuge_LLE('C3',
                         Kspecies=('Methanol', 'Glycerol'),
                         Ks=array((0.382, 0.183)),
                         top_solvents=('Biodiesel',),
                         top_split=(0.999,),
                         bot_solvents=('Water', 'Lipid', 'NaOH', 'HCl'),
                         bot_split=(0.999, 1, 1, 1))

# Vacuum dry biodiesel
# Consider doing this by split, and keeping this unit constant
# 290 Pa, 324 K according to USDA Biodiesel Model
V1 = SplitFlash('V1', outs=('', 'biodiesel'),
                order=('Water', 'Methanol', 'Biodiesel'),
                split=(0.9999, 0.9999, 0.00001),
                P=2026.5, T=331.5, Q=0)
V1.line = 'Vacuum dryer'
V1.material = 'Stainless steel 316'
P6 = Pump('P6')

"""Glycerol Purification Section"""

# Condense vacuumed methanol to recycle
HX1 = HXutility('HX1', V=0, T=295)
P11 = Pump('P11')

# Mix recycled streams and HCl
T5 = MixTank('T5')
P7 = Pump('P7')

# Centrifuge out waste fat
# assume all the lipid, free_lipid and biodiesel is washed out
C4 = SplitCentrifuge_LLE('C4', outs=('', Waste),
                         order=('Methanol', 'Glycerol', 'Water'),
                         split=(0.999, 0.999, 0.999))

# Add and mix NaOH
T6 = MixTank('T6')
P8 = Pump('P8')

# Methanol/Water distillation column
D1 = Distillation('D1',
                  LHK=('Methanol', 'Water'), P=101325,
                  y_top=0.99999, x_bot=0.0001, k=2.5)
D1.is_divided = True
D1.tray_material = 'Stainless steel 316'
D1.vessel_material = 'Stainless steel 316'

# Condense recycle methanol
HXD1 = HXutility('HXD1', V=0, T=315)
P9 = Pump('P9')

# Glycerol/Water flash (not a distillation column)
w = 0.20/Stream.species.Water.MW
g = 0.80/Stream.species.Glycerol.MW
x = w/(w+g)
D2 = Flash('D2',
           species_IDs=('Water', 'Glycerol'),
           P=101325,
           x=(x, 1-x),
           HNK=('Biodiesel',),
           LNK=('Methanol',))
D2.HasGlycolGroups = True
D2.vessel_material = 'Stainless steel 316'

# Condense water to recycle
HXD2_top = HXprocess('HXD2_top', fluid_type='ll', outs=(MixedStream(), ''),
                     species_IDs=('Methanol', 'Water'), 
                     HNK=('Glycerol', 'Biodiesel', 'NaOH', 'HCl'))
    
# Use heat from glycerol
HXD2_bot = HXprocess('HXD2_bot', fluid_type='ls',
                     species_IDs=('Methanol', 'Water'), 
                     HNK=('Glycerol', 'Biodiesel', 'NaOH', 'HCl'))
run_bot = HXD2_bot._run

# Storage tank for glycerol
T7 = StorageTank('T7', outs=Crude_glycerol)
HXD2_bot-1-T7

# Storage tank for biodiesel
T8 = StorageTank('T8', outs=Biodiesel)
V1-1-P6-0-T8

# %% Set up systems

# Biodiesel Transesterification Section
Oil-T3-P3
(P3-0, S1-0)-R1-0-C1
(C1-0, S1-1)-R2-0-C2
transesterification_network = (T3, P3, R1, C1, P4, R2, C2,
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
(C2-0, HXD2_top-0, Water, HCl1)-W1-X0-C3
def adjust_water_flow():
    HXD2_top.outs[0].liquid_mol[4] = 800*C2.outs[0]._mol[2]

# Glycerol recycle and purification section
C3-0-V1
V1-0-HX1-P11
C1-1-P4
(P4-0, C2-1, C3-1, P11-0, HCl2)-T5-P7-C4
(C4-0, NaOH)-T6-P8
(HXD2_top-1, D2-1)-HXD2_bot
HXD2_bot-0-D1-1-D2
(D2-0, P8-0)-HXD2_top
glycerol_recycle_sys = System('glycerol_recycle_sys',
                              network=(adjust_water_flow, 
                                       W1, X0, C3, V1, HX1, P11, P6, T5,
                                       P7, C4, T6, P8, HXD2_top, HXD2_bot, 
                                       D1, D2, HXD2_top, HXD2_bot),
                              recycle=D2-0)                           
# Find Fresh Methanol Flow
D1-0-HXD1-P9    # Recycle methanol
Methanol-T1-P1  # Mix fresh and recycled methanol
Catalyst-T2-P2  # Fresh catalyst
(P9-0, P1-0, P2-0)-M1-P0-S1  # Mix Catalyst with Methanol
meoh_network = (HXD1, P9, T1, P1, T2, P2, M1, P0, S1)

# Set connections for mass balance proxy
(Catalyst, Methanol, D1-0)-MB1
MB1**(1**R1, 1**R2)

# Complete System
biodiesel_sys = System('biodiesel_sys', network=transesterification_network
                                              + (glycerol_recycle_sys, MB1)
                                              + meoh_network
                                              + (T7, T8))

# Initial conditions
index = Oil.indices(['Methanol', 'Glycerol', 'Water'])
D2.outs[0].mol[index] = (3.483e-02, 8.284e-03, 2.094e+01)
D2.outs[0].T = 374.986
HXD2_top.outs[0].T = 374.986
HXD2_top.outs[0].liquid_mol[index] = (3.483e-02, 8.284e-03, 2.094e+01)

