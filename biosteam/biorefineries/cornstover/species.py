# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 23:12:04 2019

@author: yoelr
"""
import biosteam as bst
import biosteam.compounds as cp
from biosteam.biorefineries.lipidcane import species as lcspecies
import pandas as pd
from copy import copy as _copy

__all__ = ('species', 'get_grouped_species', 'species_groups')

# %%  Constants

# Common structural carbohydrates properties

# Assume heat capacity of lignin, cellulose, and hemicellulose
# and all components at 350 K are about the same [2,2].
Cp_cellulosic = 1.364

# Heat of combustions of lignocellulosic material are
# based on wood as an approximiation [4]. Also, assume
# structural carbohydrates have the same heat of combustion
# as cellulose.
# These properties match NREL's
Hc_lignin = 21e3 # (J/g)
Hc_cellulosic = 17000 # (J/g)
cal2joule = 4.184

# %% Initialize species object and define functions

lcspecies = bst.Species.tospecies([*lcspecies.biodiesel_species,
                                   *lcspecies.ethanol_species,
                                   *lcspecies.pretreatment_species])

sp = bst.Species()

def remove(ID):
    try: missing.remove(ID)
    except ValueError: raise ValueError(f"'{ID}' not a required species")

def addspecies(*IDs, cls=None, **kwargs):
    if cls is None: cls = bst.Chemical
    sp.extend([cls(i, **kwargs) for i in IDs])
    for ID in IDs: remove(ID)

def addchemical(ID, ChemID=None, phase=None, **properties): 
    setattr(sp, ID, cp.StaticChemical(ChemID or ID, phase=phase, **properties))
    remove(ID)
    
def addsubstance(ID, *args, **kwargs):
    setattr(sp, ID, cp.Substance(ID, *args, **kwargs))
    remove(ID)

def addfrom(species, ID, species_ID=None):
    if not species_ID: species_ID = ID
    cmp = _copy(getattr(species, species_ID))
    setattr(sp, ID, cmp)
    cmp.ID = species_ID
    remove(ID)

def copy(ID, ChemID):
    cmp = _copy(getattr(sp, ChemID))
    setattr(sp, ID, cmp)
    cmp.ID = ChemID
    remove(ID)

# Species to define
species_IDs = [
        'H2O', 'Ethanol', 'Glucose', 'Galactose',
        'Mannose', 'Xylose', 'Arabinose', 'Cellobiose',
        'Sucrose', 'GlucoseOligomer', 'GalactoseOligomer',
        'MannoseOligomer', 'XyloseOligomer', 'ArabinoseOligomer',
        'Extract','SolubleLignin','HMF', 'Furfural', 'AceticAcid',
        'LacticAcid', 'Xylitol', 'Glycerol', 'SuccinicAcid',
        'NH3', 'H2SO4', 'NH4SO4', 'AmmoniumAcetate', 'DAP',
        'HNO3', 'NaNO3', 'NaOH', 'CellulaseNutrients',
        'Denaturant', 'Oil', 'Cellulose', 'Galactan', 'Mannan', 'Glucan',
        'Xylan', 'Arabinan', 'Lignin', 'Acetate', 'Protein',
        'Ash', 'Enzyme', 'DenaturedEnzyme', 'Z_mobilis', 'T_reesei',
        'Biomass', 'Tar', 'CaO', 'CaSO4', 'Graphite', 'N2', 'O2', 'CO2',
        'CH4', 'H2S', 'SO2', 'NO', 'CO', 'AmmoniumSulfate', 'NO2', 'CSL',
        'WWTsludge', 'Cellulase'
]
missing = species_IDs.copy()

# %% Define species

# TODO: Add heats of combustion

# As is in data bank
addspecies('H2O', 'Ethanol', 'AceticAcid', 'Furfural', 'Glycerol',
           'H2SO4', 'LacticAcid', 'SuccinicAcid') 
sp.SuccinicAcid.Hf = -940.26e3 # kJ/mol
sp.LacticAcid.T = sp.LacticAcid.Tb = 122+273.15
sp.LacticAcid.H_int_l_T_ref_l_to_Tb = sp.LacticAcid.H
sp.LacticAcid.Hvap_Tbm = sp.LacticAcid.Hvap
sp.LacticAcid.Hf = -163122*cal2joule
addchemical('HNO3', 'NitricAcid')
addchemical('Denaturant', 'Octane')
addchemical('DAP', 'Diammonium Phosphate')
addchemical('AmmoniumAcetate')
addchemical('NH4SO4', 'AmmoniumSulfate')
addchemical('NaNO3', 'SodiumNitrate')
addchemical('Oil', 'Oleic Acid', phase='l')
addchemical('HMF')

# Will remain in the vapor phase
addspecies('N2', 'NH3', 'O2', 'CH4', 'H2S', 'SO2', cls=cp.StaticChemical)
addchemical('CO2', phase='g')

# Analagous vapors
addsubstance('NO2', MW=46.01, obj=sp.N2, Hf=7925*cal2joule)
addsubstance('NO', MW=30.01, obj=sp.N2, Hf=82.05)
addsubstance('CO', MW=28.01, obj=sp.N2, Hf=-110.522)

# Will remain as  solid
addspecies('Glucose', 'Xylose', 'Sucrose',
           cls=cp.StaticChemical, Cp=Cp_cellulosic)
sp.Glucose.Hf = -300428*cal2joule
sp.Xylose.Hf = -249440*cal2joule
sp.Sucrose.Hf = -480900*cal2joule
addspecies('CaSO4', 'Graphite', 'AmmoniumSulfate',
           cls=cp.StaticChemical, Cp=Cp_cellulosic)

# Analagous sugars
copy('Mannose', 'Glucose')
copy('Galactose', 'Glucose')
copy('Arabinose', 'Xylose')

# Other analogues
copy('CellulaseNutrients', 'Glucose')
copy('Extract', 'Glucose')
copy('Acetate', 'AceticAcid')
sp.Acetate.Hf = -103373
copy('Tar', 'Xylose')

# Species taken from previous study
addfrom(lcspecies, 'CaO')
addfrom(lcspecies, 'Ash')
addfrom(lcspecies, 'NaOH')
addsubstance('Lignin', Cp=Cp_cellulosic, Hf=-108248*cal2joule, MW=152.15)
copy('SolubleLignin', 'Lignin')

# Create structural carbohydrates
addsubstance('GlucoseOligomer', obj=sp.Glucose, Cp=Cp_cellulosic,
             MW=162.1424, Hf=-233200*cal2joule)
copy('GalactoseOligomer', 'GlucoseOligomer')
copy('MannoseOligomer', 'GlucoseOligomer')
addsubstance('XyloseOligomer', obj=sp.Xylose, Cp=Cp_cellulosic,
             MW=132.11612, Hf=-182100*cal2joule)
copy('ArabinoseOligomer', 'XyloseOligomer')

# Other
addsubstance('Z_mobilis', MW=24.6265, Hf=-31169.39*cal2joule)
addsubstance('T_reesei', MW=23.8204, Hf=-23200.01*cal2joule)
addsubstance('Biomass', MW=23.238, Hf=-23200.01*cal2joule)
addsubstance('Cellulose', MW=162.1406, Hf=-233200.06*cal2joule, Hc=Hc_cellulosic)
addsubstance('Protein', MW=22.8396, Hf=-17618*cal2joule)
addsubstance('Enzyme', MW=24.0156, Hf=-17618*cal2joule)
addsubstance('Glucan', MW=162.14, Hf=-233200*cal2joule)
addsubstance('Xylan', MW=132.12, Hf=-182100*cal2joule)
addsubstance('Xylitol', MW=152.15, Hf=-243145*cal2joule)
addsubstance('Cellobiose', MW=342.30, Hf=-480900*cal2joule)
addsubstance('CSL', MW=1, Hf=sp.Protein.Hf/4+sp.H2O.Hf/2+sp.LacticAcid.Hf/4)
copy('DenaturedEnzyme', 'Enzyme')
copy('Arabinan', 'Xylan')
copy('Mannan',   'Glucan')
copy('Galactan', 'Glucan')

# %% TODO: Maybe remove this

# missing.extend([
#     'OtherSugars', 'SugarOligomers', 
#     'InorganicSolubleSolids',
#     'Furfurals', 'OtherOrganics',
#     'OtherOrganics', 'OtherOrganics',
#     'COxSOxNOxH2S', 'OtherStructuralCarbs',
#     'CellMass',
#     'OtherInsolubleSolids',
#     'OrganicSolubleSolids',
# ])
# copy('OtherSugars', 'Arabinose')
# copy('SugarOligomers', 'GlucoseOligomer')
# copy('OrganicSolubleSolids', 'SolubleLignin')
# copy('InorganicSolubleSolids', 'DAP')
# copy('Furfurals', 'HMF')
# copy('OtherOrganics', 'Denaturant')
# copy('COxSOxNOxH2S', 'NO')
# copy('OtherStructuralCarbs', 'Arabinan')

# For waste water
copy('WWTsludge', 'Z_mobilis')
copy('Cellulase', 'Enzyme')
# copy('OtherInsolubleSolids', 'Tar')
bst.Stream.species = species = sp

# %% Grouped species

species_groups = dict(
    OtherSugars = ['Arabinose',
                   'Mannose',
                   'Galactose',
                   'Cellobiose',
                   'Sucrose'],
    SugarOligomers = ['GlucoseOligomer',
                      'XyloseOligomer',
                      'GalactoseOligomer',
                      'ArabinoseOligomer',
                      'MannoseOligomer'],
    OrganicSolubleSolids = ['AmmoniumAcetate',
                            'SolubleLignin',
                            'Extract', 
                            'LacticAcid', 
                            'Cellulase'],
    InorganicSolubleSolids = ['AmmoniumSulfate',
                            'DAP',
                            'NaOH',
                            'HNO3',
                            'NaNO3'],
    Furfurals = ['Furfural',
                 'HMF'],
    OtherOrganics = ['Glycerol',
                     'Denaturant',
                     'Oil',
                     'SuccinicAcid',
                     'Xylitol'],
    COxSOxNOxH2S = ['NO',
                    'NO2',
                    'SO2',
                    'CO',
                    'H2S'],
    Protein = ['Protein',
               'Enzyme',
               'DenaturedEnzyme'],
    CellMass = ['WWTsludge',
                'Z_mobilis',
                'T_reesei'],
    OtherInsolubleSolids = ['Tar',
                            'Ash',
                            'Graphite',
                            'Lime'],
    OtherStructuralCarbohydrates = ['Arabinan', 
                                    'Mannan', 
                                    'Galactan']
)

def get_grouped_species(stream, units='kmol/hr'):
    s = bst.Stream(species=species)
    s.setflow(flow=stream.mol, species=stream.species.IDs)
    return pd.Series({i:s.getflow(*j, units=units).sum() for i, j in species_groups.items()})

assert not missing, str(missing)

# Fix sugar properties


# cellulosic = ('Glucan', 'Xylan', 'GlucoseOligomer','XyloseOligomer',
#               'Acetate', 'SolubleLignin')

