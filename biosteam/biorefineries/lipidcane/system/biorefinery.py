# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 22:18:21 2019

@author: Guest Group
"""
from biosteam.biorefineries.lipidcane.system.ethanol import area_300, sugar_solution
from biosteam.biorefineries.lipidcane.system.biodiesel import area_400, oil
from biosteam.biorefineries.lipidcane.system.pretreatment import \
    pretreatment_sys, pretreatment_species, U202, lipid, sugar
from biosteam.biorefineries.lipidcane.tea import LipidcaneTEA
from biosteam.units.facilities import CoolingTower, \
                                      BoilerTurbogenerator, \
                                      ChilledWaterPackage, \
                                      ProcessWaterCenter
from biosteam import System, Stream, find, Species
from biosteam.units import Junction, Splitter
from biosteam.biorefineries.lipidcane.utils import set_lipid_fraction
import warnings

warnings.filterwarnings('ignore')

__all__ = ('lipidcane_sys', 'lipidcane_tea', 'area_500', 'area_600', 'BT')

# %% Facilities
Stream.species = pretreatment_species
emission = Stream('emission')
water = Species('Water',)
stream = find.stream

# Stream.default_ID_number = 500

BT = BoilerTurbogenerator('BT',
                          ins=U202-0, # Bagasse from conveyor belt
                          outs=emission,
                          boiler_efficiency=0.80,
                          turbogenerator_efficiency=0.85)

Stream.default_ID_number = 600

Stream.species = water
CT = CoolingTower('CT')
process_water_streams = (stream.biodiesel_wash_water,
                         stream.cooling_tower_makeup_water,
                         stream.boiler_makeup_water)
makeup_water = Stream('makeup_water', price=0.000254)
process_water = Stream('process_water')
pws_indices = [i.index('Water') for i in process_water_streams]
pws_flow_index_pairs = tuple(zip([i._mol for i in process_water_streams], pws_indices))
def update_water():
    process_water._mol[0] = sum([i[j] for i,j in pws_flow_index_pairs])

CWP = ChilledWaterPackage('CWP')
PWC = ProcessWaterCenter('PWC',
                         ins=('recycle_water', makeup_water),
                         outs=process_water)
        
S601 = Splitter('S601', ins=process_water, outs=process_water_streams,
                split=(1,), order=('Water',))
UO = find.unit
area_500 = System('area_500', (BT,))
area_600 = System('area_600', (CT, CWP, PWC, S601))

# %% Set up system
connect_sugar = Junction(sugar, sugar_solution, ('Water', 'Glucose', 'Sucrose'))
connect_lipid = Junction(lipid, oil, ('Lipid',))

lipidcane_sys = System('lipidcane_sys',
                       network=pretreatment_sys.network
                             + (connect_sugar, connect_lipid)
                             + area_300.network
                             + area_400.network,
                       facilities=(CWP, BT, CT, update_water, PWC))
units = lipidcane_sys._costunits

# %% Perform TEA

lipidcane_tea = LipidcaneTEA(system=lipidcane_sys, IRR=0.15, duration=(2018, 2038),
                             depreciation='MACRS7', income_tax=0.35,
                             operating_days=200, lang_factor=3,
                             construction_schedule=(0.4, 0.6), WC_over_FCI=0.05,
                             labor_cost=2.5e6, fringe_benefits=0.4,
                             property_tax=0.001, property_insurance=0.005,
                             supplies=0.20, maintenance=0.01, administration=0.005)
set_lipid_fraction(0.10)
lipidcane_sys.simulate()
lipidcane_tea.IRR = lipidcane_tea.solve_IRR()
