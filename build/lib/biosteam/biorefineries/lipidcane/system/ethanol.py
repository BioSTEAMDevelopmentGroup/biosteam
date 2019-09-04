#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 20 17:38:57 2017

The ethanol production section for the baseline lipid cane biorefinery is defined here as System objects. The systems include all streams and units starting from sugar concentration to denaturing the ethanol product.

@author: yoel
"""
# import numpy as np
from biosteam import Stream, System
from biosteam.biorefineries.lipidcane.species import ethanol_species
from biosteam.units import Mixer, Splitter, HXutility, HXprocess, \
     MultiEffectEvaporator, Fermentation, StorageTank, Pump, \
     Distillation, SolidsCentrifuge, MolecularSieve, MixTank, VentScrubber
from biosteam.biorefineries.lipidcane.process_settings import price

__all__ = ('ethanol_sys',)


# %% Species

Stream.species = ethanol_species
sp_c = ('Glucose', 'H3PO4', 'Sucrose', 'Water')
sp_r = ('Ethanol', 'Glucose', 'H3PO4', 'Water', 'DryYeast')

# %% Other

Ethanol_MW = ethanol_species.Ethanol.MW
Water_MW = ethanol_species.Water.MW

def Ethanol_molfrac(e: 'Ethanol mass fraction'):
    """Return ethanol mol fraction in a ethanol water mixture"""
    return e/Ethanol_MW / (e/Ethanol_MW + (1-e)/Water_MW)


# %% Input Streams

# Fresh water
S134 = Stream('S134', Water=5000, units='kg/hr')

# Gasoline
denature = Stream('Denaturant', Octane=230.69,
                  units='kg/hr', price=price['Gasoline'])

# Yeast
S144 = Stream('Yeast', Water=24700, DryYeast=10300,
              units='kg/hr')

# Sugar Stream (from Pretreatment section)
Sugar_solution = Stream('Sugar_solution', Glucose=1888.23, H3PO4=0,
                        Sucrose=21399.94, Water=264523.53,
                        units='kg/hr', T=99+273.15)

# Ethanol product
Ethanol = Stream('Ethanol', price=price['Ethanol'])


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
P24 = Fermentation('P24', outs=('', 'CO2'), tau=10, efficiency=0.90, N=8) 
T1 = StorageTank('T1')
T1.tau = 4
T1.line = 'Beer tank'

scrubber = VentScrubber('Scrubber', ins=(S134, P24-0), gas=('CO2',))

# Separate 99% of yeast
P19 = SolidsCentrifuge('P19', outs=('', 'Recycle yeast'),
                       split=(1, 0.99999, 1, 0.96, 0.01),
                       order=('Ethanol', 'Glucose', 'H3PO4', 
                              'Water', 'DryYeast'),
                       solids=('DryYeast',))

# Mix in Water
P23 = Mixer('P23', 'S147')
Q1 = Pump('Q1')

# Heat up before beer column
# Exchange heat with stillage
P32 = HXprocess('P32', outs=('', 'Stillage'),
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
P34 = HXutility('P34', T=115+273.15, V=0)

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
T4 = MixTank('T4', outs=Ethanol)
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

Sugar_solution-P16-1-P15-0-Q6
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

gas_index = ethanol_species.IDs.index('Octane')
def adjust_denaturant():
    denature.mol[gas_index] = 0.011*Q4.outs[0].massnet/114.232
    
P33-1-P41-0-T2-0-Q4
denature-T3-Q5
(Q5-0, Q4-0)-P39-T4
EtOH_end_network=(Q3, P41, T2, Q4, adjust_denaturant, T3, Q5, P39, T4)

(Q3-0, P15-1)-P35_0
EtOH_process_water_network=(P35_0,)    

ethanol_sys = System('ethanol_sys',
                     network=(EtOH_start_network
                            + (purification_recycle,)
                            + EtOH_end_network
                            + EtOH_process_water_network))
