# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 23:12:28 2019

@author: yoelr
"""
from biosteam import Stream, System
import biosteam as bst
from biosteam.biorefineries.cornstover.process_settings import price, ethanol_density_kggal
from biosteam.biorefineries.cornstover.species import species, get_grouped_species, species_groups
from biosteam.biorefineries.cornstover.tea import CornstoverTEA
from biosteam.biorefineries.cornstover import units
import biosteam.reaction as rxn
import numpy as np

bst.CE = 522
System.maxiter = 200
System.molar_tolerance = 1

Ethanol_MW = species.Ethanol.MW
Water_MW = species.H2O.MW

def Ethanol_molfrac(e: 'Ethanol mass fraction'):
    """Return ethanol mol fraction in a ethanol water mixture"""
    return e/Ethanol_MW / (e/Ethanol_MW + (1-e)/Water_MW)

def find_split(IDs, flow0, flow1):
    flow0 = np.asarray(flow0)
    splits = flow0/(flow0 + np.asarray(flow1))
    species = Stream.species
    array = np.zeros(len(species))  
    for ID, split in zip(IDs, splits):
        if ID in species_groups:
            array[species.indices(species_groups[ID])] = split
        else:
            array[species.index(ID)] = split
    return array


# %% Streams

bst.find.set_flowsheet(bst.Flowsheet('Cornstover'))
synonym = species.set_synonym 
synonym('CaO', 'Lime')
synonym('H2O', 'Water')
synonym('H2SO4', 'SulfuricAcid')
synonym('NH3', 'Ammonia')
synonym('Denaturant', 'Octane')
synonym('CO2', 'CarbonDioxide')
bst.Stream.species = pretreatment_species = bst.WorkingSpecies.subgroup(species, 
        ['Acetate', 'AceticAcid', 'Arabinan', 'Ash', 'Cellulase',
         'Ethanol', 'Extract', 'Furfural', 'Glucan', 'Glucose',
         'GlucoseOligomer', 'H2O', 'H2SO4', 'HMF', 'Lignin',
         'Mannan', 'NH3', 'Protein', 'SolubleLignin', 'Sucrose',
         'Xylan', 'Xylose', 'Arabinose', 'XyloseOligomer',
         'ArabinoseOligomer', 'Mannose', 'MannoseOligomer',
         'Galactan', 'Galactose', 'GalactoseOligomer'])

# feed flow
drycomposition = pretreatment_species.kwarray(
                 Glucan=0.3505, Xylan=0.1953, Lignin=0.1576,
                 Ash=0.0493, Acetate=0.0181, Protein=0.0310,
                 Extract=0.1465, Arabinan=0.0238, Galactan=0.0143,
                 Mannan=0.0060, Sucrose=0.0077)
moisture_content = pretreatment_species.kwarray(Water=0.20)
netflow = 104167.0
feedflow = netflow*(drycomposition*0.8 + moisture_content)

feed = Stream('Cornstover',
              feedflow,
              units='kg/hr',
              price=price['Feedstock'])
process_water212 = Stream('S212',
                          T=368.15,
                          P=4.7*101325,
                          Water=140000,
                          units='kg/hr')
rectifier_bottoms_product = Stream('S516',
                                   T=114+273.15,
                                   P=6.1*101325,
                                   Ethanol=18,
                                   Water=36629,
                                   Furfural=72,
                                   HMF=100,
                                   units='kg/hr')
sulfuric_acid = Stream('Sulfuric acid',
                       P=5.4*101325,
                       T=294.15,
                       Water=139,
                       SulfuricAcid=1842,
                       units='kg/hr',
                       price=price['Sulfuric acid'])
steam = Stream('Steam',
               phase='g',
               T=268+273.15,
               P=13*101325,
               Water=24534+3490,
               units='kg/hr')
ammonia = Stream('Ammonia',
                 Ammonia=1051,
                 units='kg/hr',
                 phase='l',
                 price=price['Ammonia'])
cellulase = Stream('Cellulase',
                   Cellulase=13836,
                   units='kg/hr',
                   price=price['Cellulase'])

# %% Pretreatment system

C100 = units.FeedStockHandling(ins=feed)
T201 = units.SulfuricAcidTank201(ins=sulfuric_acid)
A201 = units.SulfuricAcidMixer(ins=(rectifier_bottoms_product, T201-0))
M1 = bst.Mixer(ins=(A201-0, process_water212, C100-0))
M1_2 = units.SteamMixer(ins=(M1-0, steam), P=5.5*101325)
PRS200 = units.PretreatmentReactorSystem(ins=M1_2-0)
P1 = units.BlowdownDischargePump203(ins=PRS200-1)
T208 = units.OligomerConversionTank208(ins=P1-0)
T204 = units.Flash204(ins=T208-0, P=101325, Q=0)
M2 = bst.Mixer(ins=(PRS200-0, T204-0))
HX1 = units.WasteVaporCondenser244(ins=M2-0, T=99+273.15, V=0)
M210 = units.AmmoniaMixer210(ins=(ammonia, process_water212))
M3 = bst.Mixer(ins=(T204-1, M210-0))
T209 = units.AmmoniaAdditionTank209(ins=M3-0)
H301 = units.HydrolysateCooler301(ins=T209-0, T=48+273.15)
A308 = units.EnzymeHydrolysateMixer308(ins=(H301-0, cellulase))

solids_indices = Stream.indices(['Glucan', 'Xylan', 'Lignin',
                                 'Ash', 'Acetate', 'Protein',
                                 'Extract', 'Arabinan', 'Mannan',
                                 'Sucrose'])

sulfuric_acid_over_feed = sulfuric_acid.mol/feed.massnet
def update_sulfuric_acid_loading():
    # Also plant air
    feed_massnet = feed.massnet
    plant_air.mol[0] = 0.8 *feed_massnet
    sulfuric_acid.mol[:] = sulfuric_acid_over_feed * feed_massnet

hydrolyzate = T204.outs[1]
ammonia_over_hydrolyzate = ammonia.mol/310026.22446428984
def update_ammonia_loading():
    ammonia.mol[:] = ammonia_over_hydrolyzate * hydrolyzate.massnet

cooled_hydrolyzate = H301.outs[0]
cellulase_over_hydrolyzate = cellulase.mol/435396.4815414322
def update_cellulase_and_nutrient_loading():
    cooled_hydrolyzate_massnet = cooled_hydrolyzate.massnet
    cellulase.mol[:] = cellulase_over_hydrolyzate * cooled_hydrolyzate_massnet
    DAP310.mol[:] = DAP310_over_hydrolyzate * cooled_hydrolyzate_massnet
    DAP312.mol[:] = DAP312_over_hydrolyzate * cooled_hydrolyzate_massnet
    CSL309.mol[:] = CSL309_over_hydrolyzate * cooled_hydrolyzate_massnet
    CSL311.mol[:] = CSL311_over_hydrolyzate * cooled_hydrolyzate_massnet
    
ptsys = System('pretreatment',
               network=(C100, T201, A201, M1, M1_2,
                        PRS200, P1, T208, T204, M2,
                        HX1, M210, M3, T209, T209,
                        H301,
                        A308))

# %% Fermentation system

bst.Stream.species = fermentation_species = bst.WorkingSpecies.subgroup(species, 
        ['Acetate', 'AceticAcid', 'Arabinan', 'Ash', 'CO2', 'CSL',
         'Cellobiose', 'DAP', 'Denaturant', 'Enzyme', 'Ethanol',
         'Extract', 'Furfural', 'Glucan', 'Glucose', 'GlucoseOligomer', 
         'Glycerol', 'H2O', 'H2SO4', 'HMF', 'LacticAcid', 'Lignin',
         'Mannan', 'NH3', 'O2', 'Protein', 'SolubleLignin',
         'SuccinicAcid', 'Sucrose', 'Xylan', 'Xylitol', 'Xylose',
         'XyloseOligomer', 'Z_mobilis', 'Arabinose', 'Mannose',
         'Galactan', 'Galactose', 'GalactoseOligomer',
         'ArabinoseOligomer', 'MannoseOligomer'])

process_water274 = Stream('S274',
                          Water=150310,
                          P=5*101325,
                          T=33+273.15,
                          units='kg/hr')
DAP310 = Stream('DAP310',
                DAP=26,
                units='kg/hr',
                price=price['DAP'])
DAP312 = Stream('DAP312',
                DAP=116,
                units='kg/hr',
                price=price['DAP'])
DAP_storage = units.DapTank(ins=Stream('DAP_fresh'))
DAP_storage.link_streams()
IS0 = bst.InvSplitter(ins=DAP_storage-0, outs=(DAP310, DAP312))
CSL309 = Stream('CSL309',
                CSL=211,
                units='kg/hr',
                price=price['CSL'])
CSL311 = Stream('CSL311',
                CSL=948,
                units='kg/hr',
                price=price['CSL'])
CSL_storage = units.CslTank(ins=Stream('CLS_fresh'))
CSL_storage.link_streams()
IS1 = bst.InvSplitter(ins=CSL_storage-0, outs=(CSL309, CSL311))
denaturant = Stream('Denaturant',
                    Octane=230.69,
                    units='kg/hr',
                    price=price['Denaturant'])
process_water524 = Stream('Stripping_water',
                          Water=26836,
                          units='kg/hr')
DAP310_over_hydrolyzate = DAP310.mol/451077.22446428984
DAP312_over_hydrolyzate = DAP312.mol/451077.22446428984
CSL309_over_hydrolyzate = CSL309.mol/451077.22446428984
CSL311_over_hydrolyzate = CSL311.mol/451077.22446428984

J1 = bst.Junction(upstream=A308-0, downstream=Stream())
M4 = bst.Mixer(ins=(J1-0, None))
F300 = units.SaccharificationAndCoFermentation(ins=(M4-0, CSL311, DAP312))
M5 = bst.Mixer(ins=(F300-2, CSL309, DAP310))
F301 = units.SeedTrain(ins=M5-0)
T301 = units.SeedHoldTank301(ins=F301-1)
T301-0-1-M4

fmsys = System('fermentation',
               network=(J1, M4, F300, M5, F301, T301),
               recycle=M4-0)

# %% Ethanol purification

M6 = bst.Mixer(ins=(F300-1, None))
T306 = units.BeerTank306(ins=M6-0)

M7 = bst.Mixer(ins=(F301-0, F300-0))
T512 = bst.VentScrubber(ins=(process_water524, M7-0),
                        gas=('CO2', 'NH3', 'O2'))
T512-1-1-M6

# Heat up before beer column
# Exchange heat with stillage
P32 = bst.HXprocess('P32', ins=(T306-0, None),
                    fluid_type='ss', U=1.28)

# Beer column
xbot = Ethanol_molfrac(0.00001)
ytop = Ethanol_molfrac(0.50)
P25 = bst.Distillation('P25', ins=P32-0,
                       P=101325, y_top=ytop, x_bot=xbot,
                       k=1.25, LHK=('Ethanol', 'Water'))
P25.tray_material = 'Stainless steel 304'
P25.vessel_material = 'Stainless steel 304'
P25.BM = 2.4
P25._boiler.U = 1.85
Q2 = bst.Pump('Q2', ins=P25-1)
Q2-0-1-P32

# Mix ethanol Recycle (Set-up)
P28 = bst.Mixer('P28', ins=(P25-0, None))

ytop = Ethanol_molfrac(0.915)
P30 = bst.Distillation('P30', ins=P28-0,
                       P=101325, y_top=ytop, x_bot=xbot,
                       k=1.25, LHK=('Ethanol', 'Water'))
P30.tray_material = 'Stainless steel 304'
P30.vessel_material = 'Stainless steel 304'
P30.is_divided = True
P30._boiler.U = 1.85
P30.BM = 2.8
Q3 = bst.Pump('Q3', ins=P30-1)
Q3.link_streams()
JX = bst.Junction(Q3-0, 0-A201)

# Superheat vapor for mol sieve
P34 = bst.HXutility('P34', ins=P30-0, T=115+273.15, V=1)

# Molecular sieve
P33 = bst.MolecularSieve('P33', ins=P34-0,
                         split=(2165.14/13356.04, 1280.06/1383.85),
                         order=('Ethanol', 'Water'))

P33-0-1-P28

pure_recycle_sys = System('purification_recycle',
                          network=(P28, P30, P34, P33),
                          recycle=P28-0)

# Condense ethanol product
P41 = bst.HXutility('P41', ins=P33-1, V=0, T=350.)
T2 = bst.StorageTank('T2', ins=P41-0)
T2.line = 'Ethanol storage'
T2.tau = 7*24
Q4 = bst.Pump('Q4', ins=T2-0)

# Storage for gasoline
T3 = bst.StorageTank('T3', ins=denaturant)
T3.tau = 7*24
Q5 = bst.Pump('Q5', ins=T3-0)

# Mix in denaturant
P39 = bst.Mixer('P39', ins=(Q5-0, Q4-0))

# Denatured ethanol product
ethanol = Stream('Ethanol', price=price['Ethanol'])
T4 = bst.MixTank('T4', ins=P39-0, outs=ethanol)
T4.tau = 0.05

gas_index = denaturant.index('Octane')
def adjust_denaturant():
    denaturant.mol[gas_index] = 0.011*Q4.outs[0].massnet/114.232

Q2.BM = Q3.BM = Q4.BM = Q5.BM = 3.1
T2.BM = T3.BM = 1.7

vent_stream = M7-0
stripping_water_over_vent = process_water524.mol / 21202.490455845436
def update_stripping_water():
    process_water524.mol[:] = stripping_water_over_vent * vent_stream.massnet

puresys = System('purification',
                 network=(M7,
                          update_stripping_water,
                          T512, 
                          M6, T306,
                          P32, P25,
                          P32, Q2, P32,
                          pure_recycle_sys,
                          Q3, P41, T2, Q4,
                          adjust_denaturant,
                          T3, Q5, P39, T4, JX))

# %% Lignin Separation
bst.Stream.species = species

recycled_water = bst.Stream(Water=1,
                            T=47+273.15,
                            P=3.9*101325,
                            units='kg/hr')

splits = [('Glucose', 19, 502),
          ('Xylose', 40, 1022),
          ('OtherSugars', 81, 2175),
          ('SugarOligomers', 60, 1552),
          ('OrganicSolubleSolids', 612, 15808),
          ('InorganicSolubleSolids', 97, 2513),
          ('Furfurals', 19, 513),
          ('OtherOrganics', 52, 1348),
          ('Glucan', 1230, 25),
          ('Xylan', 415, 8),
          ('OtherStructuralCarbohydrates', 94, 2),
          ('Lignin', 12226, 250),
          ('Protein', 3376, 69),
          ('CellMass', 925, 19),
          ('OtherInsolubleSolids', 4489, 92)]

S505 = units.PressureFilter505(ins=('', recycled_water),
                               moisture_content=0.35,
                               split=find_split(*zip(*splits)))
J2 = bst.Junction(P32-1, 0-S505)

# %% Waste water treatment

def burn(reactant, O2=0, H2O=0, CO2=0, SO2=0, NO2=0, N2=0, Ash=0, NaOH=0):
    r = rxn.Reaction(f"{reactant} + {O2}O2 -> {H2O}H2O + {CO2}CO2 + {Ash}Ash + "
                     f"{SO2}SO2 + {NO2}NO2 + {N2}N2 + {NaOH}NaOH", reactant, 1.)
    cmp = getattr(species, reactant)
    species.H2O.P = 101325
    species.H2O.T = 298.15
    cmp.Hc = (cmp.Hf - (H2O*species.H2O.Hf
                        + CO2*species.CO2.Hf
                        + SO2*species.SO2.Hf
                        + NO2*species.NO2.Hf) - H2O*species.H2O.Hvapm)
    return r

combustion = rxn.ParallelReaction([    
        burn('Glucose', 6, 6, 6),
        burn('Xylose', 5, 5, 5),
        burn('Sucrose', 12, 11, 12),
        burn('Extract', 6, 6, 6),
        burn('Arabinose', 5, 5, 5),
        burn('Galactose', 6, 6, 6),
        burn('Mannose', 6, 6, 6),
        burn('GlucoseOligomer', 6, 5, 6),
        burn('Cellobiose', 12, 11, 12),
        burn('XyloseOligomer', 5, 4, 5),
        burn('MannoseOligomer', 6, 5, 6),
        burn('GalactoseOligomer', 6, 5, 6),
        burn('ArabinoseOligomer', 5, 4, 5),
        burn('Xylitol', 5.5, 6, 5),
        burn('SolubleLignin', 8.5, 4, 8),
        burn('Ethanol', 3, 3, 2),
        burn('Furfural', 5, 2, 5),
        burn('HMF', 6, 3, 6),
        burn('H2SO4', -0.5, 1, SO2=1),
        burn('CH4', 2, 2, 1),
        burn('NO', 0.5, NO2=1),
        burn('NH3', 0.75, 1.5, N2=0.5),
        burn('LacticAcid', 3, 3, 3),
        burn('AceticAcid', 2, 2, 2),
        burn('NH4SO4', 1, 4, N2=1, SO2=1),
        burn('AmmoniumAcetate', 2.75, 3.5, 2, N2=0.5),
        burn('Glycerol', 3.5, 4, 3),
        burn('SuccinicAcid', 3.5, 3, 4),
        burn('Denaturant', 12, 9, 8), # Octane
        burn('Oil', 25.5, 17, 18),
        burn('WWTsludge', 6, 17, 18),
        burn('CellulaseNutrients', 6, 17, 18),
        burn('H2S', 1.5, 1, SO2=1),
        burn('CO', 0.5, CO2=1),
        burn('HNO3', -1.75, 0.5, N2=0.5),
        burn('NaNO3', -1.25, N2=0.5, H2O=-0.5, NaOH=1),
        burn('Cellulose', 6, 5, 6),
        burn('Xylan', 5, 4, 5),
        burn('Lignin', 8.5, 4, 8),
        burn('Enzyme', 1.1975, 0.795, 1, N2=0.12, SO2=0.01),
        burn('DenaturedEnzyme', 1.1975, 0.795, 1, N2=0.12, SO2=0.01),
        burn('Biomass', 1.2185, 0.82, 1, N2=0.115, SO2=0.0035),
        burn('Z_mobilis', 1.2, 0.9, 1, N2=0.1),
        burn('Acetate', 2, 2, 2),
        burn('Arabinan', 5, 4, 5),
        burn('Mannan', 6, 5, 6),
        burn('Galactan', 6, 5, 6),
        burn('Tar', 5, 5, 5),
        burn('T_reesei', 1.19375, 0.8225, 1, N2=0.1025, SO2=0.005),
        burn('Protein', 1.2445, 0.785, 1, N2=0.145, SO2=0.007),
        burn('Graphite', 1, CO2=1),
        burn('Lime', H2O=1, Ash=1),
        burn('CaSO4', -0.5, SO2=1, Ash=1)])

def growth(reactant):
    f = getattr(species, reactant).MW / species.WWTsludge.MW
    return rxn.Reaction(f"{f}{reactant} -> WWTsludge", reactant, 1.)
    
organic_groups = ['OtherSugars', 'SugarOligomers', 'OrganicSolubleSolids',
                  'Furfurals', 'OtherOrganics', 'Protein', 'CellMass']
organics = sum([species_groups[i] for i in organic_groups],
               ['Ethanol', 'AceticAcid', 'Xylose', 'Glucose'])
organics.remove('WWTsludge')

P_sludge = 0.05/0.91/species.WWTsludge.MW
MW = np.array([species.CH4.MW, species.CO2.MW])
mass = np.array([0.51, 0.49])*MW
mass /= mass.sum()
mass *= 0.86/(0.91)
P_ch4, P_co2 = mass/MW
def anaerobic_rxn(reactant):
    MW = getattr(species, reactant).MW
    return rxn.Reaction(f"{1/MW}{reactant} -> {P_ch4}CH4 + {P_co2}CO2 + {P_sludge}WWTsludge",
                        reactant, 0.91)

# TODO: Revise this with Jeremy
anaerobic_digestion = rxn.ParallelReaction([anaerobic_rxn(i) for i in organics] + 
                                           [rxn.Reaction(f"H2SO4 -> H2S + 2O2", 'H2SO4', 1.)])


# Note, nitogenous species included here, but most of it removed in anaerobic digester
# TODO: Add ammonium to reaction, make sure it can be a liquid, possibly add Henry's constant
aerobic_digestion = rxn.ParallelReaction([i*0.74 + 0.22*growth(i.reactant)
                                          for i in combustion
                                          if (i.reactant in organics)])
aerobic_digestion.X[:] = 0.96

splits = [('Ethanol', 1, 15),
          ('Water', 27158, 356069),
          ('Glucose', 3, 42),
          ('Xylose', 7, 85),
          ('OtherSugars', 13, 175),
          ('SugarOligomers', 10, 130),
          ('OrganicSolubleSolids', 182, 2387),
          ('InorganicSolubleSolids', 8, 110),
          ('Ammonia', 48, 633),
          ('AceticAcid', 0, 5),
          ('Furfurals', 5, 70),
          ('OtherOrganics', 9, 113),
          ('Cellulose', 19, 6),
          ('Xylan', 6, 2),
          ('OtherStructuralCarbohydrates', 1, 0),
          ('Lignin', 186, 64),
          ('Protein', 51, 18),
          ('CellMass', 813, 280),
          ('OtherInsolubleSolids', 68, 23)]

well_water = Stream(Water=1, T=15+273.15)
M8 = bst.Mixer(ins=(S505-1, ''))
J3 = bst.Junction(HX1-0, 1-M8)

WWTC = units.WasteWaterSystemCost(ins=M8-0)
anaerobic = units.AnaerobicDigestion(ins=(WWTC-0, well_water),
                                     reactions=anaerobic_digestion,
                                     sludge_split=find_split(*zip(*splits)))

air = Stream('Air lagoon', O2=51061, N2=168162, phase='g', units='kg/hr')
caustic = Stream('WWT caustic', Water=2252, NaOH=2252,
                 units='kg/hr', price=price['Caustic']*0.5)
# polymer = Stream('WWT polymer') # Empty in humbird report :-/

M9 = bst.Mixer(ins=(anaerobic-1, None))

caustic_over_waste = caustic.mol / 2544300.6261793654
air_over_waste = air.mol / 2544300.6261793654
waste = M9-0
def update_aerobic_input_streams():
    waste_massnet = waste.massnet
    caustic.mol[:] = waste_massnet * caustic_over_waste
    air.mol[:] = waste_massnet * air_over_waste

aerobic = units.AerobicDigestion(ins=(waste, air, caustic),
                                 reactions=aerobic_digestion)

splits = [('Ethanol', 0, 1),
          ('Water', 381300, 2241169),
          ('Glucose', 0, 2),
          ('Xylose', 1, 3),
          ('OtherSugars', 1, 7),
          ('SugarOligomers', 1, 6),
          ('OrganicSolubleSolids', 79, 466),
          ('InorganicSolubleSolids', 4828, 28378),
          ('Ammonia', 3, 16),
          ('Furfurals', 0, 3),
          ('OtherOrganics', 1, 7),
          ('CarbonDioxide', 6, 38),
          ('O2', 3, 17),
          ('N2', 5, 32),
          ('Cellulose', 0, 194),
          ('Xylan', 0, 65),
          ('OtherStructuralCarbohydrates', 0, 15),
          ('Lignin', 0, 1925),
          ('Protein', 0, 90),
          ('CellMass', 0, 19778),
          ('OtherInsolubleSolids', 0, 707)]

S0 = bst.Splitter(ins=aerobic-1, split=find_split(*zip(*splits)))

S1 = bst.Splitter(ins=S0-1, split=0.96)

M10 = bst.Mixer(ins=(S1-0, None))
M10-0-1-M9

M11 = bst.Mixer(ins=(anaerobic-2, S1-1))

centrifuge_species = ('Water', 'Glucose', 'Xylose', 'OtherSugars', 'SugarOligomers', 'OrganicSolubleSolids', 'InorganicSolubleSolids', 'Furfurals', 'OtherOrganics', 'CO2', 'COxSOxNOxH2S', 'Cellulose', 'Xylan', 'OtherStructuralCarbohydrates', 'Lignin', 'Protein', 'CellMass', 'OtherInsolubleSolids')
S623_flow = np.array([7708, 0, 0, 1, 1, 13, 75, 3, 0, 1, 1, 2, 25, 8, 2, 250, 52, 1523, 92])
S616_flow = np.array([109098, 3, 6, 13, 9, 187, 1068, 46, 5, 8, 14, 31, 1, 0, 0, 13, 3, 80, 5])

S2 = bst.Splitter(ins=M11-0, outs=('', 'Sludge'),
                  split=find_split(centrifuge_species, S616_flow, S623_flow))
S2-0-1-M10

brine_species = ('Water',  'Xylose', 'OtherSugars', 'SugarOligomers', 'OrganicSolubleSolids', 'InorganicSolubleSolids', 'Furfurals', 'OtherOrganics', 'CO2', 'COxSOxNOxH2S')
S627_flow = np.array([4967, 1, 1, 1, 79, 4828, 1, 3, 44])
S626_flow = np.array([376324, 0, 0,0, 0, 0,    0,  0, 0])

S3 = bst.Splitter(ins=S0-0, outs=('treated_water', 'waste_brine'),
                  split=find_split(brine_species, S626_flow, S627_flow))

WWTsystem = System('aerobic digestion system',
                   network=(M9, aerobic, S0, S1, M11, S2, M10),
                   recycle=M9-0)

# %% Facilities
M12 = bst.Mixer(ins=(S2-1, S505-0))
BT = bst.facilities.BoilerTurbogenerator(ins=M12-0, 
                                         turbogenerator_efficiency=0.85)
BT.cost_items['Turbogenerator'].n = 0.6
CWP = bst.facilities.ChilledWaterPackage()
CT = bst.facilities.CoolingTower()
process_water = bst.Stream(ID='Process_water', species=bst.Species('Water'))

           # flow                # index
process_water_data = ((caustic.mol,          caustic.index('Water')),
                      (process_water524.mol, process_water524.index('Water')),
                      (process_water212.mol, process_water212.index('Water')),
                      (process_water274.mol, process_water274.index('Water')),
                      (steam.mol,         steam.index('Water')))

def update_water_loss():
    process_water.mol[0] = sum([i[j] for i, j in process_water_data])
        
water = bst.Species('Water')
Makeup_water = Stream('Makeup_water', species=water, price=price['Makeup water'])
BT_water_makeup = Stream('BT_water_makeup', species=water, price=price['Makeup water'])
CT_water_makeup = Stream('CT_water_makeup', species=water, price=price['Makeup water'])

PWC = bst.facilities.ProcessWaterCenter(ins=(S3-0, Makeup_water,
                                             BT_water_makeup, CT_water_makeup),
                                        outs=(process_water, 'BT_water', 'CT_water'),
                                        BT=BT, CT=CT)
PWC.CT_blowdown = 0.0015
PWC.boiler_blowdown = 0.03
PWC.RO_rejection = 0
M8.ins.append(PWC.outs[1])

substance = bst.Species('Substance', cls=bst.Substance)
ash = Stream('Ash', species=substance,
             price=price['Ash disposal'])
lime = Stream('Lime', species=substance,
             price=price['FGD lime'])
boilerchems = Stream('Boiler_chemicals', species=substance,
                     price=price['Boiler chems'])
boilerchems_mol = boilerchems._mol
ash_mol = ash._mol
lime_mol = lime._mol
emission_mol = BT.outs[0]._mol
ash_index = BT.outs[0].index('Ash')
def update_lime_boilerchems_and_ash():
    emission_ash = emission_mol[ash_index]
    lime_mol[0] = lime = emission_ash * 0.21
    ash_mol[0] = (emission_ash + lime) * 1.18 # Include lime and other stuff
    boilerchems_mol[0] = 0.24620/865 * lime

CIP = Stream('CIP', species=substance, flow=(126,))
CIP_package = units.CIPpackage(ins=CIP)

plant_air = Stream('Plant air', flow=(83333,), species=substance)

ADP = bst.facilities.AirDistributionPackage(ins=plant_air)
ADP.link_streams()


FW = units.FireWaterTank(ins=Stream('Fire water', flow=(8343,),species=substance))

# %% Complete system

cornstover_sys = System('corn stover system',
                        network=(ptsys, fmsys, puresys, J2, S505,
                                 J3, M8, WWTC, anaerobic,
                                 update_aerobic_input_streams,
                                 WWTsystem, S3),
                        facilities=(M12, CWP, BT, CT, update_water_loss,
                                    PWC, ADP, update_lime_boilerchems_and_ash,
                                    CIP_package, IS0, IS1, DAP_storage,
                                    CSL_storage, FW))
cornstover_sys.products.update((ash, boilerchems))
baghouse_bags = Stream(ID='Baghouse_bags', species=substance, flow=(1,), price=11.1)
cornstover_sys.feeds.add(lime)
cornstover_sys.feeds.add(baghouse_bags)
WWTsystem.converge_method = 'Wegstein'
cornstover_sys.simulate()
cornstover_sys.simulate()
cornstover_sys.simulate()
ethanol_tea = CornstoverTEA(
        system=cornstover_sys, IRR=0.10, duration=(2007, 2037),
        depreciation='MACRS7', income_tax=0.35, operating_days=350.4,
        lang_factor=None, construction_schedule=(0.08, 0.60, 0.32),
        startup_months=3, startup_FOCfrac=1, startup_salesfrac=0.5,
        startup_VOCfrac=0.75, WC_over_FCI=0.05,
        finance_interest=0.08, finance_years=10, finance_fraction=0.4,
        OSBL_units=(WWTC, CWP, CT, PWC, ADP), # BT not included
        warehouse=0.04, site_development=0.09, additional_piping=0.045,
        proratable_costs=0.10, field_expenses=0.10, construction=0.20,
        contingency=0.10, other_indirect_costs=0.10, labor_cost=2.5e6,
        labor_burden=0.90, property_insurance=0.007, maintenance=0.03)
ethanol_tea.units.remove(BT)

Area800 = bst.TEA.like(System('Area800', (BT,)),
                       ethanol_tea)
Area800.labor_cost = 0
Area800.depreciation = 'MACRS20'
Area800.OSBL_units = (BT,)
cornstover_tea = bst.CombinedTEA([ethanol_tea, Area800], IRR=0.10)
cornstover_sys._TEA = cornstover_tea
ethanol.price = cornstover_tea.solve_price(ethanol, ethanol_tea)
ethanol.price = cornstover_tea.solve_price(ethanol, ethanol_tea)
ethanol_price_USDgal = ethanol.price * ethanol_density_kggal

# %% Areas

Area100 = bst.TEA.like(System('Area100', (C100,)), ethanol_tea)
Area200 = bst.TEA.like(System('Area200',
                              (T201, A201, PRS200, P1,
                               T208, T204, HX1, M210, T209)),
                       ethanol_tea)
Area300 = bst.TEA.like(System('Area300', 
                              (H301, A308, F300,
                               F301, T301, T306)),
                       ethanol_tea)
Area400 = () # Cellulase production (empty for now)
Area500 = bst.TEA.like(System('Area500', 
                              (T512, P32, P25, Q2,
                               P28, P30, Q3, P34,
                               P33, P41, P39, S505)),
                       ethanol_tea)
Area600 = bst.TEA.like(System('Area600', (WWTC,)),
                       ethanol_tea)
Area700 = bst.TEA.like(System('Area700',
                              (T2, T3, Q4, Q5, FW,
                               CSL_storage, DAP_storage)),
                       ethanol_tea) 
Area900 = bst.TEA.like(System('Area900', (CWP, CT, PWC, ADP, CIP_package)),
                       ethanol_tea)
areas = (Area100, Area200, Area300, Area500,
         Area600, Area700, Area800, Area900)
installation_costs = {i.system.ID: i.installation_cost/1e6
                      for i in areas}
utility_costs = {i.system.ID: i.utility_cost/1e6
                 for i in areas}

def get_utility(units, ID, attr):
    out = 0
    for i in units:
        if i._N_heat_utilities:
            for j in i._heat_utilities:
                if j.ID == ID:
                    out += getattr(j, attr)
    return out

get_rate = lambda units: sum([i._power_utility.rate
                              for i in units
                              if i._has_power_utility])

get_ecost = lambda units: sum([i._power_utility.cost
                               for i in units
                               if i._has_power_utility])*24*350.4/1e6

cooling_water_uses = {i.system.ID: 18.01528*get_utility(i.units, 'Cooling water', 'flow')/1e6
                      for i in areas}
electricity_uses = {i.system.ID: get_rate(i.units)/1e3/41 for i in areas}

electricity_costs = {i.system.ID: get_ecost(i.units) for i in areas}

