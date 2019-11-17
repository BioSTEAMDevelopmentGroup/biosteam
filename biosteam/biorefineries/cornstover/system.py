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

def Ethanol_molfrac(e):
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

bst.find.set_flowsheet(bst.Flowsheet('cornstover'))
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

# bst.Stream.default_ID_number = 100
    
# feed flow
drycomposition = pretreatment_species.kwarray(
                 Glucan=0.3505, Xylan=0.1953, Lignin=0.1576,
                 Ash=0.0493, Acetate=0.0181, Protein=0.0310,
                 Extract=0.1465, Arabinan=0.0238, Galactan=0.0143,
                 Mannan=0.0060, Sucrose=0.0077)
moisture_content = pretreatment_species.kwarray(Water=0.20)
netflow = 104167.0
feedflow = netflow*(drycomposition*0.8 + moisture_content)

cornstover = Stream('cornstover',
                    feedflow,
                    units='kg/hr',
                    price=price['Feedstock'])
warm_process_water = Stream('warm_process_water',
                         T=368.15,
                         P=4.7*101325,
                         Water=140000,
                         units='kg/hr')
rectifier_bottoms_product = Stream('',
                                   T=114+273.15,
                                   P=6.1*101325,
                                   Ethanol=18,
                                   Water=36629,
                                   Furfural=72,
                                   HMF=100,
                                   units='kg/hr')
sulfuric_acid = Stream('sulfuric_acid',
                       P=5.4*101325,
                       T=294.15,
                       Water=139,
                       SulfuricAcid=1842,
                       units='kg/hr',
                       price=price['Sulfuric acid'])
steam = Stream('steam',
               phase='g',
               T=268+273.15,
               P=13*101325,
               Water=24534+3490,
               units='kg/hr')
ammonia = Stream('ammonia',
                 Ammonia=1051,
                 units='kg/hr',
                 phase='l',
                 price=price['Ammonia'])
cellulase = Stream('cellulase',
                   units='kg/hr',
                   price=price['Enzyme'])

# %% Pretreatment system

U101 = units.FeedStockHandling('U101', ins=cornstover)
U101.cost_items['System'].cost = 0

# bst.Stream.default_ID_number = 200

T201 = units.SulfuricAcidTank('T201', ins=sulfuric_acid)
M201 = units.SulfuricAcidMixer('M201', ins=(rectifier_bottoms_product, T201-0))
M202 = bst.Mixer('M202', ins=(M201-0, warm_process_water, U101-0))
M203 = units.SteamMixer('M203', ins=(M202-0, steam), P=5.5*101325)
R201 = units.PretreatmentReactorSystem('R201', ins=M203-0)
P201 = units.BlowdownDischargePump('P201', ins=R201-1)
T202 = units.OligomerConversionTank('T202', ins=P201-0)
F201 = units.PretreatmentFlash('F201', ins=T202-0, P=101325, Q=0)
M204 = bst.Mixer('M204', ins=(R201-0, F201-0))
H201 = units.WasteVaporCondenser('H201', ins=M204-0, T=99+273.15, V=0)
M210 = units.AmmoniaMixer('M210', ins=(ammonia, warm_process_water))
M205 = bst.Mixer('M205', ins=(F201-1, M210-0))
T203 = units.AmmoniaAdditionTank('T203', ins=M205-0)

# bst.Stream.default_ID_number = 300

H301 = units.HydrolysateCooler('H301', ins=T203-0, T=48+273.15)
M301 = units.EnzymeHydrolysateMixer('M301', ins=(H301-0, cellulase))

solids_indices = Stream.indices(['Glucan', 'Xylan', 'Lignin',
                                 'Ash', 'Acetate', 'Protein',
                                 'Extract', 'Arabinan', 'Mannan',
                                 'Sucrose'])

sulfuric_acid_over_feed = sulfuric_acid.mol/cornstover.massnet
def update_sulfuric_acid_loading():
    # Also plant air
    feed_massnet = cornstover.massnet
    plant_air.mol[0] = 0.8 *feed_massnet
    sulfuric_acid.mol[:] = sulfuric_acid_over_feed * feed_massnet

hydrolyzate = F201.outs[1]
ammonia_over_hydrolyzate = ammonia.mol/310026.22446428984
def update_ammonia_loading():
    ammonia.mol[:] = ammonia_over_hydrolyzate * hydrolyzate.massnet

cooled_hydrolyzate = H301.outs[0]
cellulose_index = cornstover.index('Glucan')
enzyme_over_cellulose = 20/1000 * 10 # (20 g enzyme / cellulose) / (100 g cellulase / 1 L enzyme)
water_cellulase_mass = cellulase.mass[cellulase.indices(('Water', 'Cellulase'))]
water_cellulase_to_cellulase = np.array([0.9, 0.1])
def update_cellulase_and_nutrient_loading():
    cooled_hydrolyzate_massnet = cooled_hydrolyzate.massnet
    cellulose = cornstover._mass[cellulose_index]
    # Note: An additional 10% is produced for the media glucose/sophorose mixture
    water_cellulase_mass[:] = (enzyme_over_cellulose
                               * water_cellulase_to_cellulase
                               * cellulose * 1.1) 
    DAP1.mol[:] = DAP1_over_hydrolyzate * cooled_hydrolyzate_massnet
    DAP2.mol[:] = DAP2_over_hydrolyzate * cooled_hydrolyzate_massnet
    CSL1.mol[:] = CSL1_over_hydrolyzate * cooled_hydrolyzate_massnet
    CSL2.mol[:] = CSL2_over_hydrolyzate * cooled_hydrolyzate_massnet
    
pretreatment_sys = System('pretreatment_sys',
               network=(U101, T201, M201, M202, M203,
                        R201, P201, T202, F201, M204,
                        H201, M210, M205, T203, T203,
                        H301,
                        M301))

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

DAP1 = Stream('DAP1',
                DAP=26,
                units='kg/hr',
                price=price['DAP'])
DAP2 = Stream('DAP2',
                DAP=116,
                units='kg/hr',
                price=price['DAP'])
DAP_storage = units.DAPTank('DAP_storage', ins=Stream('DAP_fresh'), outs='DAP')
DAP_storage.link_streams()
S301 = bst.InvSplitter('S301', ins=DAP_storage-0, outs=(DAP1, DAP2))
CSL1 = Stream('CSL1',
                CSL=211,
                units='kg/hr',
                price=price['CSL'])
CSL2 = Stream('CSL2',
                CSL=948,
                units='kg/hr',
                price=price['CSL'])
CSL_storage = units.CSLTank('CSL_storage', ins=Stream('CSL_fresh'), outs='CSL')
CSL_storage.link_streams()
S302 = bst.InvSplitter('S302', ins=CSL_storage-0, outs=(CSL1, CSL2))
denaturant = Stream('denaturant',
                    Octane=230.69,
                    units='kg/hr',
                    price=price['Denaturant'])
stripping_water = Stream('stripping_water',
                          Water=26836,
                          units='kg/hr')
DAP1_over_hydrolyzate = DAP1.mol/451077.22446428984
DAP2_over_hydrolyzate = DAP2.mol/451077.22446428984
CSL1_over_hydrolyzate = CSL1.mol/451077.22446428984
CSL2_over_hydrolyzate = CSL2.mol/451077.22446428984

J1 = bst.Junction(upstream=M301-0, downstream=Stream())
M302 = bst.Mixer('M302', ins=(J1-0, None))
R301 = units.SaccharificationAndCoFermentation('R301', ins=(M302-0, CSL2, DAP2))
M303 = bst.Mixer('M303', ins=(R301-2, CSL1, DAP1))
R302 = units.SeedTrain('R302', ins=M303-0)
T301 = units.SeedHoldTank('T301', ins=R302-1)
T301-0-1-M302

fermentation_sys = System('fermentation_sys',
               network=(update_cellulase_and_nutrient_loading,
                        J1, M302, R301, M303, R302, T301),
               recycle=M302-0)

# %% Ethanol purification

M304 = bst.Mixer('M304', ins=(R302-0, R301-0))
T302 = units.BeerTank('T302')

# bst.Stream.default_ID_number = 400

M401 = bst.Mixer('M401', ins=(R301-1, None))
M401-0-T302

D401 = bst.VentScrubber('D401', ins=(stripping_water, M304-0),
                        gas=('CO2', 'NH3', 'O2'))
D401-1-1-M401

# Heat up before beer column
# Exchange heat with stillage
H401 = bst.HXprocess('H401', ins=(T302-0, None),
                    fluid_type='ss', U=1.28)

# Beer column
xbot = Ethanol_molfrac(0.00001)
ytop = Ethanol_molfrac(0.50)
D402 = bst.Distillation('D402', ins=H401-0,
                       P=101325, y_top=ytop, x_bot=xbot,
                       k=1.25, LHK=('Ethanol', 'Water'))
D402.tray_material = 'Stainless steel 304'
D402.vessel_material = 'Stainless steel 304'
D402.BM = 2.4
D402._boiler.U = 1.85
P401 = bst.Pump('P401', ins=D402-1)
P401-0-1-H401

# Mix ethanol Recycle (Set-up)
M402 = bst.Mixer('M402', ins=(D402-0, None))

ytop = Ethanol_molfrac(0.915)
D403 = bst.Distillation('D403', ins=M402-0,
                       P=101325, y_top=ytop, x_bot=xbot,
                       k=1.25, LHK=('Ethanol', 'Water'))
D403.tray_material = 'Stainless steel 304'
D403.vessel_material = 'Stainless steel 304'
D403.is_divided = True
D403._boiler.U = 1.85
D403.BM = 2.8
P402 = bst.Pump('P402', ins=D403-1)
P402.link_streams()
JX = bst.Junction(P402-0, 0-M201)

# Superheat vapor for mol sieve
H402 = bst.HXutility('H402', ins=D403-0, T=115+273.15, V=1)

# Molecular sieve
U401 = bst.MolecularSieve('U401', ins=H402-0,
                         split=(2165.14/13356.04, 1280.06/1383.85),
                         order=('Ethanol', 'Water'))

U401-0-1-M402

ethanol_recycle_sys = System('ethanol_recycle_sys',
                             network=(M402, D403, H402, U401),
                             recycle=M402-0)

# Condense ethanol product
H403 = bst.HXutility('H403', ins=U401-1, V=0, T=350.)

# IDnum_400 = bst.Stream.default_ID_number

# bst.Stream.default_ID_number = 700
T701 = bst.StorageTank('T701', ins=H403-0)
T701.line = 'Ethanol storage'
T701.tau = 7*24
P701 = bst.Pump('P701', ins=T701-0)

# Storage for gasoline
T702 = bst.StorageTank('T702', ins=denaturant)
T702.tau = 7*24
P702 = bst.Pump('P702', ins=T702-0)

# Mix in denaturant
ethanol = Stream('ethanol', price=price['Ethanol'])
M701 = bst.MixTank('M701', ins=(P702-0, P701-0), outs=ethanol)
M701.line = 'Mixer'
M701.tau = 0.05

gas_index = denaturant.index('Octane')
def adjust_denaturant():
    denaturant.mol[gas_index] = 0.022*P701.outs[0].massnet/114.232

P401.BM = P402.BM = P701.BM = P702.BM = 3.1
T701.BM = T702.BM = 1.7

vent_stream = M304-0
stripping_water_over_vent = stripping_water.mol / 21202.490455845436
def update_stripping_water():
    stripping_water.mol[:] = stripping_water_over_vent * vent_stream.massnet

puresys = System('purification',
                 network=(M304,
                          update_stripping_water,
                          D401, 
                          M401, T302,
                          H401, D402,
                          H401, P401, H401,
                          ethanol_recycle_sys,
                          P402, H403, T701, P701,
                          adjust_denaturant,
                          T702, P702, M701, JX))

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

# bst.Stream.default_ID_number = IDnum_400

S401 = units.PressureFilter('S401', ins=('', recycled_water),
                            moisture_content=0.35,
                            split=find_split(*zip(*splits)))
J2 = bst.Junction(H401-1, 0-S401)

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


# Note, nitogenous species included here, but most of it removed in R601 digester
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

# bst.Stream.default_ID_number = 600

well_water = Stream('well_water', Water=1, T=15+273.15)
M601 = bst.Mixer(ins=(S401-1, '', '', ''))
J3 = bst.Junction(H201-0, 1-M601)

WWTC = units.WasteWaterSystemCost('WWTC', ins=M601-0)
R601 = units.AnaerobicDigestion('R601', ins=(WWTC-0, well_water),
                                 reactions=anaerobic_digestion,
                                 sludge_split=find_split(*zip(*splits)))

air = Stream('air_lagoon', O2=51061, N2=168162, phase='g', units='kg/hr')
caustic = Stream('WWT_caustic', Water=2252, NaOH=2252,
                 units='kg/hr', price=price['Caustic']*0.5)
# polymer = Stream('WWT polymer') # Empty in humbird report :-/

M602 = bst.Mixer(ins=(R601-1, None))

caustic_over_waste = caustic.mol / 2544300.6261793654
air_over_waste = air.mol / 2544300.6261793654
waste = M602-0
def update_aerobic_input_streams():
    waste_massnet = waste.massnet
    caustic.mol[:] = waste_massnet * caustic_over_waste
    air.mol[:] = waste_massnet * air_over_waste

R602 = units.AerobicDigestion('R602', ins=(waste, air, caustic),
                              outs=('evaporated_water', ''),
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

S601 = bst.Splitter('S601', ins=R602-1, split=find_split(*zip(*splits)))

S602 = bst.Splitter('S602', ins=S601-1, split=0.96)

M603 = bst.Mixer('M603', ins=(S602-0, None))
M603-0-1-M602

M604 = bst.Mixer('M604', ins=(R601-2, S602-1))

centrifuge_species = ('Water', 'Glucose', 'Xylose', 'OtherSugars', 'SugarOligomers', 'OrganicSolubleSolids', 'InorganicSolubleSolids', 'Furfurals', 'OtherOrganics', 'CO2', 'COxSOxNOxH2S', 'Cellulose', 'Xylan', 'OtherStructuralCarbohydrates', 'Lignin', 'Protein', 'CellMass', 'OtherInsolubleSolids')
S623_flow = np.array([7708, 0, 0, 1, 1, 13, 75, 3, 0, 1, 1, 2, 25, 8, 2, 250, 52, 1523, 92])
S616_flow = np.array([109098, 3, 6, 13, 9, 187, 1068, 46, 5, 8, 14, 31, 1, 0, 0, 13, 3, 80, 5])

S603 = bst.Splitter('S603', ins=M604-0, outs=('', 'sludge'),
                  split=find_split(centrifuge_species, S616_flow, S623_flow))
S603-0-1-M603

brine_species = ('Water',  'Xylose', 'OtherSugars', 'SugarOligomers', 'OrganicSolubleSolids', 'InorganicSolubleSolids', 'Furfurals', 'OtherOrganics', 'CO2', 'COxSOxNOxH2S')
S627_flow = np.array([4967, 1, 1, 1, 79, 4828, 1, 3, 44])
S626_flow = np.array([376324, 0, 0,0, 0, 0,    0,  0, 0])

S604 = bst.Splitter('S604', ins=S601-0, outs=('treated_water', 'waste_brine'),
                  split=find_split(brine_species, S626_flow, S627_flow))

aerobic_digestion_sys = System('aerobic_digestion_sys',
                               network=(M602, R602, S601, S602, M604, S603, M603),
                               recycle=M602-0)

# %% Facilities

# bst.Stream.default_ID_number = 500

M501 = bst.Mixer(ins=(S603-1, S401-0))
BT = bst.facilities.BoilerTurbogenerator('BT', ins=M501-0, 
                                         turbogenerator_efficiency=0.85)
BT.outs[1].T = 373.15
BT.cost_items['Turbogenerator'].n = 0.6

# bst.Stream.default_ID_number = 700

CWP = bst.facilities.ChilledWaterPackage('CWP')
CT = bst.facilities.CoolingTower('CT')
CT.outs[1].T = 273.15 + 28
process_water = bst.Stream(ID='process_water',
                           species=bst.Species('Water'))

                       # flow                   # index
process_water_data = ((caustic.mol,             caustic.index('Water')),
                      (stripping_water.mol,     stripping_water.index('Water')),
                      (warm_process_water.mol,  warm_process_water.index('Water')),
                      (steam.mol,               steam.index('Water')),
                      (BT.outs[1].mol,          BT.outs[1].index('Water')),
                      (CT.outs[1].mol,          CT.outs[1].index('Water')))

def update_water_loss():
    process_water.mol[0] = sum([i[j] for i, j in process_water_data])
        
water = bst.Species('Water')
makeup_water = Stream('makeup_water', species=water, price=price['Makeup water'])

PWC = bst.facilities.ProcessWaterCenter('PWC',
                                        ins=(S604-0, makeup_water),
                                        outs=(process_water,))
J4 = bst.Junction(BT.outs[1], 2**M601)
J5 = bst.Junction(CT.outs[1], 3**M601)

substance = bst.Species('substance', cls=bst.Substance)
ash = Stream('ash', species=substance,
             price=price['Ash disposal'])
lime = Stream('lime', species=substance,
             price=price['FGD lime'])
boilerchems = Stream('boiler_chemicals', species=substance,
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
CIP_package = units.CIPpackage('CIP_package', ins=CIP)

plant_air = Stream('plant_air', flow=(83333,), species=substance)

ADP = bst.facilities.AirDistributionPackage('ADP', ins=plant_air)
ADP.link_streams()


FW = units.FireWaterTank('FT',
                         ins=Stream('fire_water', flow=(8343,), species=substance))

# %% Complete system

cornstover_sys = System('cornstover_sys',
                        network=(pretreatment_sys, fermentation_sys, puresys, J2, S401,
                                 J3, J4, J5, M601, WWTC, R601,
                                 update_aerobic_input_streams,
                                 aerobic_digestion_sys, S604),
                        facilities=(M501, CWP, BT, CT, update_water_loss,
                                    PWC, ADP, update_lime_boilerchems_and_ash,
                                    CIP_package, S301, S302, DAP_storage,
                                    CSL_storage, FW))
cornstover_sys.products.update((ash, boilerchems))
baghouse_bags = Stream(ID='Baghouse_bags', species=substance, flow=(1,), price=11.1)
cornstover_sys.feeds.add(lime)
cornstover_sys.feeds.add(baghouse_bags)
aerobic_digestion_sys.converge_method = 'Fixed point'
for i in range(3): cornstover_sys.simulate()
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

Area700 = bst.TEA.like(System(None, (M501, BT,)),
                       ethanol_tea)
Area700.labor_cost = 0
Area700.depreciation = 'MACRS7'
Area700.OSBL_units = (BT,)
cornstover_tea = bst.CombinedTEA([ethanol_tea, Area700], IRR=0.10)
cornstover_sys._TEA = cornstover_tea
ethanol.price = cornstover_tea.solve_price(ethanol, ethanol_tea)
ethanol.price = cornstover_tea.solve_price(ethanol, ethanol_tea)
ethanol_price_USDgal = ethanol.price * ethanol_density_kggal

# %% Areas

Area100 = bst.TEA.like(System(None, (U101,)), ethanol_tea)
Area200 = bst.TEA.like(System(None,
                              (T201, M201, R201, P201,
                               T202, F201, H201, M210, T203)),
                       ethanol_tea)
Area300 = bst.TEA.like(System(None, 
                              (H301, M301, R301,
                               R302, T301, T302)),
                       ethanol_tea)
Area400 = bst.TEA.like(System(None, 
                              (D401, H401, D402, P401,
                               M402, D403, P402, H402,
                               U401, H403, M701, S401)),
                       ethanol_tea)
Area500 = bst.TEA.like(System(None, (WWTC,)),
                       ethanol_tea)
Area600 = bst.TEA.like(System(None,
                              (T701, T702, P701, P702, M701, FW,
                               CSL_storage, DAP_storage)),
                       ethanol_tea) 
Area800 = bst.TEA.like(System(None, (CWP, CT, PWC, ADP, CIP_package)),
                       ethanol_tea)
areas = (Area100, Area200, Area300, Area400,
         Area500, Area600, Area700, Area800)
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
                              if i._has_power_utility])/1e3

get_ecost = lambda units: sum([i._power_utility.cost
                               for i in units
                               if i._has_power_utility])*24*350.4/1e6

cooling_water_uses = {i.system.ID: get_utility(i.units, 'Cooling water', 'duty')/1e6/4.182
                      for i in areas}
electricity_uses = {i: get_rate(j.units)/41 for i,j in enumerate(areas)}
electricity_costs = {i.system.ID: get_ecost(i.units) for i in areas}

