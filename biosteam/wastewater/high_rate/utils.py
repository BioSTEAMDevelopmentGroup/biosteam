# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-2024, Yalin Li <mailto.yalin.li@gmail.com>
# Copyright (C) 2023, Yoel Rene Cortes-Pena <yoelcortes@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

"""
This module contains data and utility functions for the high-rate wastewater treatment system.
"""

import biosteam as bst
import numpy as np
from chemicals.elements import molecular_weight
from warnings import warn
from thermosteam import Chemical, Chemicals, Stream, settings, Thermo
from thermosteam.reaction import (
    Reaction as Rxn,
    ParallelReaction as PRxn
)
mix_tank_purchase_cost_algorithms = bst.design_tools.mix_tank_purchase_cost_algorithms
TankPurchaseCostAlgorithm = bst.design_tools.TankPurchaseCostAlgorithm
ExponentialFunctor = bst.utils.ExponentialFunctor
remove_undefined_chemicals = bst.utils.remove_undefined_chemicals
default_chemical_dict = bst.utils.default_chemical_dict

__all__ = (
    # Combustion
    'get_combustion_energy',
    # Construction
    'IC_purchase_cost_algorithms',
    'select_pipe',
    'cost_pump',
    # Chemicals
    'default_insolubles',
    'get_insoluble_IDs',
    'get_soluble_IDs',
    'append_wwt_chemicals',
    'create_missing_wwt_chemicals',
    # Digestion
    'get_BD_dct',
    'get_digestable_chemicals',
    'compute_stream_COD',
    'get_digestion_rxns',
    # Miscellaneous
    'format_str',
    'remove_undefined_chemicals',
    'get_split_dct',
    'kph_to_tpd',
    # TEA/LCA
    'prices',
    'GWP_CFs',
    )


# %%

# =============================================================================
# Combustion
# =============================================================================

combustion_IDs = ['O2', 'H2O', 'CO2', 'N2', 'SO2', 'P4O10', 'Ash']
def create_combustion_chemicals():
    chems = Chemicals(combustion_IDs[:-1])
    for chem in chems:
        if chem.ID != 'P4O10': chem.at_state('g')
        else: chem.at_state('s')
    chems.append([Chemical('Ash', phase='s', phase_ref='s', search_db=False, MW=1.)])
    return chems


def get_combustion_energy(stream, combustion_eff=0.8):
    """
    Estimate the amount of energy generated from combustion of incoming streams.
    """
    chems = stream.chemicals
    to_add = []
    
    for ID in combustion_IDs:
        if not hasattr(chems, ID): to_add.append(ID)
    if to_add:
        combustion_chems = create_combustion_chemicals()
        for ID in to_add:
            to_add.append(combustion_chems.pop(ID))
        new_chems = Chemicals(to_add)
        new_chems.compile()
        settings.set_thermo(new_chems)
        new_stream = Stream()
        new_stream.imass[chems.IDs] = stream.imass[chems.IDs]
    else:
        new_chems = chems
        new_stream = stream
    rxns = new_stream.chemicals.get_combustion_reactions()
    reacted = new_stream.copy()
    rxns.force_reaction(reacted.mol)
    reacted.imol['O2'] = 0
    H_net = new_stream.H + new_stream.HHV - reacted.H
    settings.set_thermo(chems)
    return H_net


# %%

# =============================================================================
# Construction
# =============================================================================

##### Internal Circulation Reactor #####
IC_purchase_cost_algorithms = mix_tank_purchase_cost_algorithms.copy()
conventional = IC_purchase_cost_algorithms['Conventional']
# The cost correlation might not hold for the ranges beyond
ic = TankPurchaseCostAlgorithm(
    ExponentialFunctor(A=conventional.f_Cp.A,
                       n=conventional.f_Cp.n),
    V_min=np.pi/4*1.5**2*16, # 1.5 and 16 are the lower bounds of the width and height ranges in Kontos
    V_max=np.pi/4*12**2*25, # 12 and 25 are the lower bounds of the width and height ranges in Kontos
    V_units='m^3',
    CE=conventional.CE,
    material='Stainless steel')

IC_purchase_cost_algorithms['IC'] = ic

##### Pipe Selection #####
# Based on ANSI (American National Standards Institute) pipe chart
# the original code has a bug (no data for 22) and has been fixed here
boundaries = np.concatenate([
    np.arange(1/8, 0.5, 1/8),
    np.arange(0.5, 1.5, 1/4),
    np.arange(1.5, 5,   0.5),
    np.arange(5,   12,  1),
    np.arange(12,  36,  2),
    np.arange(36,  54,  6)
    ])

size = boundaries.shape[0]

pipe_dct = {
    1/8 : (0.405,  0.049), # OD (outer diameter), t (wall thickness)
    1/4 : (0.540,  0.065),
    3/8 : (0.675,  0.065),
    1/2 : (0.840,  0.083),
    3/4 : (1.050,  0.083),
    1   : (1.315,  0.109),
    1.25: (1.660,  0.109),
    1.5 : (1.900,  0.109),
    2   : (2.375,  0.109),
    2.5 : (2.875,  0.120),
    3   : (3.500,  0.120),
    3.5 : (4.000,  0.120),
    4   : (4.500,  0.120),
    4.5 : (5.000,  0.120),
    5   : (5.563,  0.134),
    6   : (6.625,  0.134),
    7   : (7.625,  0.134),
    8   : (8.625,  0.148),
    9   : (9.625,  0.148),
    10  : (10.750, 0.165),
    11  : (11.750, 0.165),
    12  : (12.750, 0.180),
    14  : (14.000, 0.188),
    16  : (16.000, 0.199),
    18  : (18.000, 0.188),
    20  : (20.000, 0.218),
    22  : (22.000, 0.250),
    24  : (24.000, 0.250),
    26  : (26.000, 0.250),
    28  : (28.000, 0.250),
    30  : (30.000, 0.312),
    32  : (32.000, 0.312),
    34  : (34.000, 0.312),
    36  : (36.000, 0.312),
    42  : (42.000, 0.312),
    48  : (48.000, 0.312)
    }


def select_pipe(Q, v):
    """Select pipe based on Q (flow in ft3/s) and velocity (ft/s)"""
    A = Q / v # cross-section area
    d = (4*A/np.pi) ** 0.5 # minimum inner diameter, [ft]
    d *= 12 # minimum inner diameter, [in]
    d_index = np.searchsorted(boundaries, d, side='left') # a[i-1] < v <= a[i]
    d_index = d_index-1 if d_index==size else d_index # if beyond the largest size
    OD, t = pipe_dct[boundaries[d_index]]
    ID = OD - 2*t # inner diameter, [in]
    return OD, t, ID


##### Pumping #####
def cost_pump(unit):
    """
    Calculate the cost of the pump and pump building for a unit
    based on its `Q_mgd` (hydraulic flow in million gallons per day),
    `recir_ratio` (recirculation ratio) attributes.
    """

    Q_mgd, recir_ratio = unit.Q_mgd, unit.recir_ratio

    # Installed pump cost, this is a fitted curve
    pumps = 2.065e5 + 7.721*1e4*Q_mgd

    # Design capacity of intermediate pumps, gpm,
    # 2 is the excess capacity factor to handle peak flows
    GPMI = 2 * Q_mgd * 1e6 / 24 / 60

    # Design capacity of recirculation pumps, gpm
    GPMR = recir_ratio * Q_mgd * 1e6 / 24 / 60

    building = 0.
    for GPM in (GPMI, GPMR):
        if GPM == 0:
            N = 0
        else:
            N = 1 # number of buildings
            GPMi = GPM
            while GPMi > 80000:
                N += 1
                GPMi = GPM / N

        PBA = N * (0.0284*GPM+640) # pump building area, [ft]
        building += 90 * PBA

    return pumps, building


# %%

# =============================================================================
# Chemicals
# =============================================================================

#!!! This should be individually added when adding the WWT process
default_insolubles = {
    # Sugarcane
    'CaO', 'Flocculant', 'Solids', 'Yeast',
    # Corn
    'Fiber', 'InsolubleProtein',
    # Corn stover
    'Ash', 'CaSO4', 'Cellulose', 'DenaturedEnzyme', 'Enzyme', 'Hemicellulose',
    'Lignin', 'Lime', 'P4O10', 'Protein', 'Tar', 'T_reesei', 'WWTsludge', 'Z_mobilis',
    # Lactic acid
     'Arabinan', 'Acetate', 'BaghouseBag', 'CalciumDihydroxide', 'CoolingTowerChems',
     'FermMicrobe', 'Galactan', 'Glucan','Mannan', 'Polymer', 'Xylan',
    # Oilcane
    'Biomass', 'Cellmass', 'DryYeast',
    # Others
    'CellMass',
    }

def get_insoluble_IDs(chemicals, insolubles):
    chem_IDs = set([i.ID for i in chemicals])
    new_insolubles = set(insolubles).intersection(chem_IDs)
    return tuple(new_insolubles)

def get_soluble_IDs(chemicals, insolubles):
    return tuple([i.ID for i in chemicals if not i.ID in insolubles])


_cal2joule = 4.184 # auom('cal').conversion_factor('J')

def create_missing_wwt_chemicals(chemicals):
    new_chemicals = Chemicals([])
    def add_chemical(ID, **data):
        missing = not hasattr(chemicals, ID)
        if missing:
            chemical = Chemical(ID, **data)
            new_chemicals.append(chemical)
        return missing

    add_chemical('NH3', phase='g', Hf=-10963*_cal2joule)
    add_chemical('H2S', phase='g', Hf=-4927*_cal2joule)
    add_chemical('SO2', phase='g')
    add_chemical('NH4OH', search_ID='AmmoniumHydroxide', phase='l', Hf=-336719)
    add_chemical('H2SO4', phase='l')
    add_chemical('HCl', phase='l')
    add_chemical('HNO3', phase='l', Hf=-41406*_cal2joule)
    add_chemical('NaOH', phase='l'),
    add_chemical('NaNO3', phase='l', Hf=-118756*_cal2joule)
    add_chemical('Na2SO4', phase='l', Hf=-1356380)
    if add_chemical('CaSO4', phase='s', Hf=-342531*_cal2joule):
        new_chemicals.CaSO4.Cn.move_up_model_priority('LASTOVKA_S', 0)
    add_chemical('NaOCl', phase='l', Hf=-347*1e3) # https://en.wikipedia.org/wiki/Sodium_hypochlorite
    add_chemical('CitricAcid', phase='l', Hf=-1543.8*1e3) # https://en.wikipedia.org/wiki/Citric_acid
    add_chemical('Bisulfite', phase='l')
    add_chemical('WWTsludge', search_db=False, phase='s',
                formula='CH1.64O0.39N0.23S0.0035', Hf=-23200.01*_cal2joule)
    if add_chemical('Polymer', search_db=False, phase='s', MW=1, Hf=0, HHV=0, LHV=0):
        new_chemicals.Polymer.Cn.add_model(evaluate=0, name='Constant')
    for i in new_chemicals: i.default()
    return new_chemicals


def append_wwt_chemicals(chemicals, set_thermo=True):
    chemicals = chemicals.chemicals if isinstance(chemicals, Thermo) else chemicals
    required_chemicals = (
        'NH3', 'H2S', 'SO2', 'NH4OH', 'H2SO4', 'HCl', 'HNO3', 'NaOH', 'NaNO3',
        'Na2SO4', 'CaSO4', 'NaOCl', 'CitricAcid', 'Bisulfite', 'WWTsludge', 'Polymer'
    )
    if all([(i in chemicals) for i in required_chemicals]): return chemicals
    chems = Chemicals([*chemicals, *create_missing_wwt_chemicals(chemicals)])

    # Add aliases and groups
    if isinstance(chemicals, bst.CompiledChemicals):
        chems.compile()
        get = getattr
        for grp, comp in chemicals._group_mol_compositions.items():
            group_IDs = [chem.ID for chem in get(chemicals, grp)]
            chems.define_group(grp, group_IDs, comp)

    if set_thermo: settings.set_thermo(chems)
    return chems



# %%

# =============================================================================
# Digestion
# =============================================================================

def get_CHONSP(chemical):
    organic = True
    atoms = chemical.atoms

    CHONSP = []
    for atom in ('C', 'H', 'O', 'N', 'S', 'P'):
        CHONSP.append(atoms.get(atom) or 0.,)

    if CHONSP[0] <= 0 or CHONSP[1] <= 0: # does not have C or H
        if not (len(atoms) == 1 and CHONSP[1] == 2): # count H2 as organic
            organic = False

    if sum(v for v in atoms.values()) != sum(CHONSP): # contains other elements
        organic = False

    return CHONSP if organic else [0.]*6


def get_COD_stoichiometry(chemical):
    r"""
    Get the molar stoichiometry for the theoretical
    chemical oxygen demand (COD) of a given chemical as in:

    .. math::
        C_nH_aO_bN_cS_dP_e + \frac{2n+0.5a-b-1.5c+3d+2.5e}{2}O_2
        -> nCO_2 + \frac{a-3c-2d}{2}H_2O + cNH_3 + dH_2SO_4 + \frac{e}{4}P_4O_{10}
    """
    Xs = nC, nH, nO, nN, nS, nP = get_CHONSP(chemical)

    excluded = ('O2', 'CO2', 'H2O', 'NH3', 'H2SO4', 'P4O10',
                *default_insolubles)
    if chemical.ID in excluded or chemical.locked_state=='g':
        return dict.fromkeys(excluded, 0)

    dct = {
        chemical.ID: -1. if sum([abs(i) for i in Xs])!=0 else 0.,
        'O2': -(nC+0.25*nH-0.5*nO-0.75*nN+1.5*nS+1.25*nP),
        'CO2': nC,
        'H2O': 0.5*nH-1.5*nN-nS, # assume one water reacts with SO3 to H2SO4
        'NH3': nN,
        'H2SO4': nS,
        'P4O10': 0.25*nP
        }

    return dct


def get_digestable_chemicals(chemicals):
    chems = [chemicals[i.ID] for i in chemicals
             if get_COD_stoichiometry(i)['O2']!=0]
    return chems


def get_BMP_stoichiometry(chemical):
    r"""
    Compute the theoretical biochemical methane potential (BMP) in
    mol :math:`CH_4`/mol chemical of a given chemical using:

    .. math::
        C_vH_wO_xN_yS_z + \frac{4v-w-2x+3y+2z}{2}H2O ->
        \frac{4v+w-2x-3y-2z}{8}CH4 + \frac{(4v-w+2x+3y+2z)}{8}CO2 + yNH_3 + zH_2S
    """
    Xs = nC, nH, nO, nN, nS, nP = get_CHONSP(chemical)

    excluded = ('H2O', 'CH4', 'CO2', 'NH3', 'H2S',
                *default_insolubles)
    if chemical.ID in excluded or chemical.locked_state=='g':
        return dict.fromkeys(excluded, 0)

    dct = {
        chemical.ID: -1. if sum([abs(i) for i in Xs])!=0 else 0.,
        'H2O': -(nC-0.25*nH-0.5*nO+0.75*nN+0.5*nS),
        'CH4': 0.5*nC+0.125*nH-0.25*nO-0.375*nN-0.25*nS,
        'CO2': 0.5*nC-0.125*nH+0.25*nO+0.375*nN+0.25*nS,
        'NH3': nN,
        'H2S': nS,
        }

    return dct


# Note that these biodegradabilities will then be multiplied by the yield
# of biogas/cell mass (0.86 is the default)
def get_BD_dct(chemicals, default_BD=1, **kwargs):
    BD_dct = dict.fromkeys([i.ID for i in get_digestable_chemicals(chemicals)],
                           default_BD)

    # Based on Kontos
    BD_dct['AceticAcid'] = 1 # 0.87/0.86>1
    BD_dct['Arabinose'] = 0.2/0.86
    BD_dct['Glucose'] = 1 # 0.87/0.86>1
    BD_dct['GlucoseOligomer'] = BD_dct['Glucan'] = 0.81/0.86
    BD_dct['HMF'] = 0.85/0.86
    BD_dct['LacticAcid'] = 0.85/0.86
    BD_dct['Lignin'] = BD_dct['SolubleLignin'] = 0.001/0.86
    BD_dct['Tar'] = 0.
    BD_dct['Xylose'] = 0.8/0.86
    BD_dct['XyloseOligomer'] = BD_dct['Xylan'] = 0.75/0.86

    # Assume biodegradabilities based on glucose and xylose
    BD_dct['Galactose'] = BD_dct['Mannose'] = \
        BD_dct['Arabinose']/BD_dct['Xylose'] * BD_dct['Glucose']
    BD_dct['GalactoseOligomer'] = BD_dct['Galactan'] = \
        BD_dct['MannoseOligomer'] = BD_dct['Mannan'] = \
             BD_dct['Glucan']

    BD_dct['ArabinoseOligomer'] = BD_dct['Arabinan'] = BD_dct['Xylan']

    if kwargs: # other input biodegradabilities
        BD_dct.update(kwargs)

    return BD_dct


def compute_stream_COD(stream):
    r"""
    Compute the chemical oxygen demand (COD) of a given stream in kg-O2/m3
    by summing the COD of each chemical in the stream using:

    .. math::
        COD [\frac{kg}{m^3}] = mol_{chemical} [\frac{kmol}{m^3}] * \frac{g O_2}{mol chemical}
    """
    chems = stream.chemicals
    mol = stream.mol
    if stream.F_vol == 0: return 0
    iCOD = np.array([-get_COD_stoichiometry(i)['O2'] for i in chems])
    COD = (mol*iCOD).sum()*molecular_weight({'O': 2}) / stream.F_vol
    return COD


def get_digestion_rxns(stream, BD, X_biogas, X_growth, biomass_ID):
    biomass_MW = getattr(stream.chemicals, biomass_ID).MW
    chems = [i for i in stream.chemicals if i.ID!=biomass_ID]
    BD = 1. if not BD else BD
    if isinstance(BD, float):
        BD = dict.fromkeys([i.ID for i in chems], BD)

    if X_biogas+X_growth > 1:
        raise ValueError('Sum of `X_biogas`/`X_decomp` and `X_biogas` is '
                         f'{X_biogas+X_growth}, larger than 100%.')

    biogas_rxns = []
    growth_rxns = []
    for i in chems:
        X = BD.get(i.ID)
        if not X:
            continue # assume no entry means not biodegradable

        biogas_stoyk = get_BMP_stoichiometry(i)
        if not biogas_stoyk.get(i.ID): # no conversion of this chemical
            continue

        iX_biogas = X * X_biogas # the amount of chemical used for biogas production
        iX_growth = X * X_growth # the amount of chemical used for cell growth

        if iX_biogas: # do not check atomic balance as P will not be accounted for
            biogas_rxn = Rxn(reaction=biogas_stoyk, reactant=i.ID, X=iX_biogas,
                             check_atomic_balance=False)
            biogas_rxns.append(biogas_rxn)

        if iX_growth:
        # Cannot check atom balance since the substrate may not have the atom
            growth_rxn = Rxn(f'{i.ID} -> {i.MW/biomass_MW}{biomass_ID}',
                             reactant=i.ID, X=iX_growth,
                             check_atomic_balance=False)


            growth_rxns.append(growth_rxn)

    if len(biogas_rxns)+len(growth_rxns)>1:
        return PRxn(biogas_rxns+growth_rxns)

    return []


# %%

# =============================================================================
# Miscellaneous
# =============================================================================

def format_str(string):
    string = string.replace(' ', '_')
    string = string.replace('-', '_')
    return string


def get_split_dct(chemicals, **split):
    # Copied from the corn stover biorefinery,
    # which is based on the 2011 NREL report (Humbird et al.),
    # assume no insolubles go to permeate
    insolubles_dct = dict.fromkeys(default_insolubles, 0.)
    split_dct = dict(
        Water=0.1454,
        Glycerol=0.125,
        LacticAcid=0.145,
        SuccinicAcid=0.125,
        HNO3=0.1454,
        Denaturant=0.125,
        DAP=0.1454,
        AmmoniumAcetate=0.145,
        AmmoniumSulfate=0.1454,
        H2SO4=0.1454,
        NaNO3=0.1454,
        Oil=0.125,
        N2=0.1351,
        NH3=0.1579,
        O2=0.15,
        CO2=0.1364,
        Xylose=0.25,
        Sucrose=0.125,
        Mannose=0.125,
        Galactose=0.125,
        Arabinose=0.125,
        Extract=0.145,
        NaOH=0.1454,
        SolubleLignin=0.145,
        GlucoseOligomer=0.1429,
        GalactoseOligomer=0.1429,
        MannoseOligomer=0.1429,
        XyloseOligomer=0.1429,
        ArabinoseOligomer=0.1429,
        Xylitol=0.125,
        Cellobiose=0.125,
        Cellulase=0.145
        )
    split_dct.update(insolubles_dct)
    default_chemical_dict(split_dct, chemicals, 0.15, 0.125, 0) # 'g', 'l', 's'

    if split is not None:
        split_dct.update(split)

    remove_undefined_chemicals(split_dct, chemicals)

    return split_dct


def kph_to_tpd(stream):
    dry_mass = stream.F_mass - stream.imass['Water']
    factor = 0.026455471462185312 # auom('kg').conversion_factor('ton')/auom('hr').conversion_factor('day')
    return dry_mass*factor


# %%

# =============================================================================
# Related to techno-economic analysis (TEA)
# =============================================================================

# Renewable natural gas, RIN D3, average of 2015 (first year with year-round data)-2021
# https://www.epa.gov/fuels-registration-reporting-and-compliance-help/rin-trades-and-price-information
# D3 designation based on entry Q of the approved pathways
# https://www.epa.gov/renewable-fuel-standard-program/approved-pathways-renewable-fuel
RIN_per_gal = 1.843 # $/gal ethanol

# According to the statue:
# https://www.ecfr.gov/current/title-40/chapter-I/subchapter-C/part-80/subpart-M/section-80.1415
# https://www.ecfr.gov/current/title-40/chapter-I/subchapter-C/part-80/subpart-M/section-80.1426
# "(5) 77,000 Btu (lower heating value) of compressed natural gas (CNG) or liquefied natural gas (LNG) shall represent one gallon of renewable fuel with an equivalence value of 1.0."
# All heating values from the H2 tools, accessed 9/27/2022
# https://h2tools.org/hyarc/calculator-tools/lower-and-higher-heating-values-fuels
LHV_LNG = 48.62 # MJ/kg, liquefied natural gas
MJ_per_BTU = 0.001055
LHV_LNG /= MJ_per_BTU # Btu/kg
LHV_ethanol = 77000 # Btu/gal per the statue instead of the 76330 from the H2 tools website
LNG_to_ethanol = LHV_LNG/LHV_ethanol # gal ethanol/kg LNG
RIN_price = RIN_per_gal * LNG_to_ethanol # $/kg LNG

prices = {
    'bisulfite': 0.08, # $/L
    'citric_acid': 0.22, # $/L
    'naocl': 0.14, # $/L
    'RIN': RIN_price, # in addition to the natural gas price
    }

# =============================================================================
# # Related to life cycle assessment (LCA)
# =============================================================================
# 100-year global warming potential (GWP) in kg CO2-eq/kg dry material
# All ecoinvent entries are from v3.8, allocation at the point of substitution

# Ecoinvent, market for sodium hydrogen sulfite, GLO,
# converted to 38% solution
bisulfite_CF = 1.2871 * 0.38

# Ecoinvent, market for citric acid, GLO
citric_acid_CF = 5.9048

# Ecoinvent, market for sodium hypochlorite, without water, in 15% solution state, RoW,
# converted to 12.5 wt% solution (15 vol%)
naocl_CF = 2.4871 * 0.125

# GREET, North America, from shale and conventional recovery, average US
ng_CF = 0.3899 + 1/16.04246*44.0095

GWP_CFs = {
    'Bisulfite': bisulfite_CF,
    'CitricAcid': citric_acid_CF,
    'CH4': ng_CF,
    'H2SO4': 43.3831/1e3, # assumed to be for the solution
    'HCl': 1.9683, # in the US, assumed to be concentrated solution
    'NaOCl': naocl_CF,
    'NaOH': 2.0092,
    'NH3': 2.6355,
}


def add_CFs(stream_registry, unit_registry, stream_CF_dct):
    has_steam = False
    for ID, key_factor in stream_CF_dct.items():
        if ID == 'steam':
            has_steam = True
            continue
        stream = stream_registry.search(ID) if isinstance(ID, str) \
            else getattr(unit_registry.search(ID[0]), ID[1])[ID[2]]
        if stream: # some streams might only exist in exist/new systems
            try:
                iter(key_factor)
                if isinstance(key_factor, str): key_factor = (key_factor,)
            except:
                key_factor = (key_factor,)
            key, factor = (key_factor[0], 1.) if len(key_factor) == 1 else key_factor

            stream.characterization_factors['GWP'] = GWP_CFs[key]*factor
    bst.PowerUtility.set_CF('GWP', *GWP_CFs['Electricity'])
    if has_steam:
        steam = stream_registry.search('steam')
        MJ = steam.H / 1e3 # enthalpy in kJ/hr
        steam.characterization_factors['GWP'] = MJ/steam.F_mass * GWP_CFs['Steam'] 