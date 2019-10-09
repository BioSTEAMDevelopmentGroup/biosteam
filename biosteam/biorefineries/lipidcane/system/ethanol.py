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

__all__ = ('area_300',)


# %% Species

# Stream.default_ID_number = 300
Stream.species = ethanol_species
sp_c = ('Glucose', 'H3PO4', 'Sucrose', 'Water')
sp_r = ('Ethanol', 'Glucose', 'H3PO4', 'Water', 'DryYeast')

# %% Other

MW_etoh = ethanol_species.Ethanol.MW
MW_water = ethanol_species.Water.MW

def Ethanol_molfrac(e: 'Ethanol mass fraction'):
    """Return ethanol mol fraction in a ethanol water mixture"""
    return e/MW_etoh / (e/MW_etoh + (1-e)/MW_water)


# %% Input Streams

# Fresh water
stripping_water = Stream('stripping_water', Water=5000, units='kg/hr')

# Gasoline
denaturant = Stream('denaturant', Octane=230.69,
                  units='kg/hr', price=price['Gasoline'])

# Yeast
yeast = Stream('yeast', Water=24700, DryYeast=10300,
               units='kg/hr')

# Sugar Stream (from Pretreatment section)
sugar_solution = Stream('sugar_solution', Glucose=1888.23, H3PO4=0,
                        Sucrose=21399.94, Water=264523.53,
                        units='kg/hr', T=99+273.15)

# Ethanol product
ethanol = Stream('ethanol', price=price['Ethanol'])


# %% Units

# Split sugar solution
S301 = Splitter('S301',
               split=0.265)

# Concentrate sugars
F301 = MultiEffectEvaporator('F301',
                            component='Water', # component being evaporated
                            P=(101325, 73581, 50892, 32777, 20000),
                            V=0.95248) # fraction evaporated
F301.components['condenser'].U = 1.85
# Note: value of steam ~ 6.86 for the following 
# (101325, 73580.467, 50891.17, 32777.406, 19999.925, 11331.5),

# Mix sugar solutions
M301 = Mixer('M301')

# Cool for fermentation
H301 = HXutility('H301', T=295.15)

# Ethanol Production
R301 = Fermentation('R301', outs=('CO2', ''), tau=10, efficiency=0.90, N=8) 
T301 = StorageTank('T301')
T301.tau = 4
T301.line = 'Beer tank'

D301 = VentScrubber('D301', ins=(stripping_water, R301-0), gas=('CO2',))

# Separate 99% of yeast
C301 = SolidsCentrifuge('C301', outs=('', 'recycle_yeast'),
                       split=(1, 0.99999, 1, 0.96, 0.01),
                       order=('Ethanol', 'Glucose', 'H3PO4', 
                              'Water', 'DryYeast'),
                       solids=('DryYeast',))

# Mix in Water
M302 = Mixer('M302')
P301 = Pump('P301')

# Heat up before beer column
# Exchange heat with stillage
H302 = HXprocess('H302', outs=('', 'stillage'),
                fluid_type='ss', U=1.28)

# Beer column
xbot = Ethanol_molfrac(0.00001)
ytop = Ethanol_molfrac(0.574)
D302 = Distillation('D302', P=101325,
                   y_top=ytop, x_bot=xbot, k=1.20,
                   LHK=('Ethanol', 'Water'))
D302.tray_material = 'Stainless steel 316'
D302.vessel_material = 'Stainless steel 316'
D302._boiler.U = 1.85
P302 = Pump('P302')

# Mix ethanol Recycle (Set-up)
M303 = Mixer('M303')

ytop = Ethanol_molfrac(0.9061726)
D303 = Distillation('D303', P=101325,
                    y_top=ytop, x_bot=xbot, k=1.20,
                    LHK=('Ethanol', 'Water'))
D303.tray_material = 'Stainless steel 316'
D303.vessel_material = 'Stainless steel 316'
D303.is_divided = True
D303._boiler.U = 1.85
P303 = Pump('P303')

# Superheat vapor for mol sieve
H303 = HXutility('H303', T=115+273.15, V=1)

# Molecular sieve
U301 = MolecularSieve('U301',
                     split=(2165.14/13356.04, 1280.06/1383.85),
                     order=('Ethanol', 'Water'))

# Condense ethanol product
H304 = HXutility('H304', 'S149', V=0, T=350.)
T302 = StorageTank('T302')
T302.line = 'Ethanol storage'
T302.tau = 6*24
P304 = Pump('P304')

# Storage for gasoline
T303 = StorageTank('T303')
T303.tau = 6*24
P305 = Pump('P305')

# denaturantd ethanol product
T304 = MixTank('T304', outs=ethanol)
T304.tau = 0.10

# Add denaturant
M304 = Mixer('M304')

# Recycle water to Process Condensate Tank
M305 = Mixer('M305', outs='recycle_water')

# Yeast mixing
T305 = MixTank('T305')
T305.tau = 0.1
yeast-T305

# Multi-effect evaporator pumps
P306 = Pump('P306')


# %% Set up EtOH system

sugar_solution-S301-1-F301-0-P306
(S301-0, P306-0)-M301-H301
(H301-0, yeast-T305-0)-R301-1-T301-0-C301
(C301-0, D301-1)-M302-P301
(P301-0, P302-0)-H302-0-D302-1-P302
EtOH_start_network = (S301, F301, P306, M301, H301, T305, R301, T301,
                      C301, D301, M302, P301, H302, D302, P302, H302)

(D302-0, U301-0)-M303-0-D303-0-H303-U301
D303-1-P303
ethanol_recycle_sys = System('ethanol_recycle_sys',
                           network=(M303, D303, H303, U301),
                           recycle=M303-0)

gas_index = ethanol_species.IDs.index('Octane')
def adjust_denaturant():
    denaturant.mol[gas_index] = 0.021*P304.outs[0].massnet/114.232
    
U301-1-H304-0-T302-0-P304
denaturant-T303-P305
(P305-0, P304-0)-M304-T304
EtOH_end_network=(P303, H304, T302, P304,
                  adjust_denaturant, T303,
                  P305, M304, T304)

(P303-0, F301-1)-M305
EtOH_process_water_network=(M305,)    

area_300 = System('area_300',
                  network=(EtOH_start_network
                         + (ethanol_recycle_sys,)
                         + EtOH_end_network
                         + EtOH_process_water_network))
