# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2024, Yoel Cortes-Pena <yoelcortes@gmail.com>
#               2023-2024, Yalin Li <mailto.yalin.li@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
This module contains unit operations for wastewater treatment of a 
cellulosic ethanol biorefinery as in [1]_.

.. contents:: :local:

Data
----
.. autodata:: biosteam.wastewater.conventional.non_digestables
    
Unit operations
---------------
.. autoclass:: biosteam.wastewater.conventional.AnaerobicDigestion 
.. autoclass:: biosteam.wastewater.conventional.AerobicDigestion
.. autoclass:: biosteam.wastewater.conventional.ReverseOsmosis
.. autoclass:: biosteam.wastewater.conventional.WastewaterSystemCost

Utilities
---------
.. autofunction:: biosteam.wastewater.conventional.get_digestable_organic_chemicals

System factories
----------------
.. autofunction:: biosteam.wastewater.conventional.create_conventional_wastewater_treatment_system

References
----------
.. [1] Humbird, D., Davis, R., Tao, L., Kinchin, C., Hsu, D., Aden, A.,
    Dudgeon, D. (2011). Process Design and Economics for Biochemical 
    Conversion of Lignocellulosic Biomass to Ethanol: Dilute-Acid 
    Pretreatment and Enzymatic Hydrolysis of Corn Stover
    (No. NREL/TP-5100-47764, 1013269). https://doi.org/10.2172/1013269

"""
from biosteam.utils import remove_undefined_chemicals, default_chemical_dict
import biosteam as bst
import thermosteam as tmo
from thermosteam import (
    utils,
    settings,
    separations,
    Reaction as Rxn,
    ParallelReaction as PRxn,
)
cost = bst.decorators.cost

__all__ = (
    'AnaerobicDigestion', 
    'AerobicDigestion', 
    'ReverseOsmosis',
    'WastewaterSystemCost',
    'SludgeCentrifuge',
    'get_digestable_organic_chemicals',
    'create_conventional_wastewater_treatment_system',
)


# %% Functional utilities

#: list[str] IDs for non-digestible components in wastewater
non_digestables = ['WWTsludge', 'Cellulose', 'Xylan', 'CellulaseNutrients',
                   'Mannan', 'Lignin', 'Galactan', 'Glucan', 'Acetate',
                   'Biomass', 'Arabinan', 'Tar', 'CO', 'CO2', 'CH4']

def get_digestable_organic_chemicals(chemicals):
    """
    Return a list of digestible organic chemical IDs.

    Parameters
    ----------
    chemicals : :class:`~thermosteam.Chemicals`
        Digestible organic chemicals will be retrieve from this parameter.

    """
    non_digestable_chemicals = set([chemicals[i] for i in non_digestables if i in chemicals])
    digestables = [i for i in chemicals if i not in non_digestable_chemicals]
    return [i for i in digestables if i.locked_state != 'g' and 'C' in i.atoms]

def anaerobic_digestion_reactions(
        chemicals, MW_sludge,
        biogas_CH4_fraction=0.51, # mol-CH4 / mol-biogas
        organics_to_biogas=0.86, # g-biogas / g-reactant
        organics_to_biomass=0.05, # g-biomass / g-reactant
        thermo=None,
    ):
    # Defaults are based on P49 in Humbird et al., 91% of organic components is destroyed,	
    # of which 86% is converted to biogas and 5% is converted to sludge,	
    # and the biogas is assumed to be 51% CH4 and 49% CO2 on a dry molar basis	
    MW_CH4 = 16.04246
    MW_CO2 = 44.0095
    x_CH4 = biogas_CH4_fraction
    x_CO2 = (1. - x_CH4)
    MW_biogas = x_CH4 * MW_CH4 +  x_CO2 * MW_CO2 # g-biogas / mol-biogas
    conversion = organics_to_biogas + organics_to_biomass # g-reacted / g-reactant
    f_biogas = organics_to_biogas / (conversion * MW_biogas) # mol-biogas / g-reacted
    f_sludge = organics_to_biomass / (conversion * MW_sludge) # mol-biomass / g-reacted
    f_CH4 = x_CH4 * f_biogas # mol-CH4 / g-reacted
    f_CO2 = x_CO2 * f_biogas # mol-CO2 / g-reacted
    thermo = settings.get_default_thermo(thermo)
    parsable_name = thermo.chemicals.get_parsable_synonym
    isvalid = utils.is_valid_ID
    def anaerobic_rxn(chemical):
        reactant = chemical.ID
        if not isvalid(reactant): reactant = parsable_name(reactant)
        if reactant == 'H2SO4':
            return Rxn("H2SO4 -> H2S + 2O2", 'H2SO4', 1.)
        else:
            MW_inv = 1. / chemical.MW	
            return Rxn(f'{MW_inv}{reactant} -> {f_CH4}CH4 + {f_CO2}CO2 + {f_sludge}WWTsludge',	
                       reactant, 0.91, chemicals=thermo.chemicals)	    
    
    return PRxn([anaerobic_rxn(i) for i in chemicals])

def aerobic_digestion_reactions(chemicals, MW_sludge, X_combustion=0.74, X_growth=0.22, thermo=None):
    # Based on P49 in Humbird et al. Defaults assume 96% of remaining soluble 
    # organic matter is removed after aerobic digestion, of which 74% is 
    # converted to water and CO2 and 22% to cell mass
    isvalid = utils.is_valid_ID
    thermo = settings.get_default_thermo(thermo)
    parsable_name = thermo.chemicals.get_parsable_synonym
    def growth(chemical):
        f = chemical.MW / MW_sludge
        reactant = chemical.ID
        if not isvalid(reactant): reactant = parsable_name(reactant)
        return Rxn(f"{reactant} -> {f}WWTsludge", reactant, X_growth, chemicals=thermo.chemicals)
    return PRxn([i.get_combustion_reaction(conversion=X_combustion) + growth(i)
                 for i in chemicals])

# %% Unit operations

@cost('Flow rate', 'Wastewater system', units='kg/hr', CE=551,
      cost=50280080., n=0.6, BM=1, kW=7139/1.05, S=393100)
class WastewaterSystemCost(bst.Unit): 
    """
    Create a unit that estimates the capital cost and electricity demand
    of a wastewater treatment system.
    
    Parameters
    ----------
    ins : 
        Wastewater.
        
    """


class AnaerobicDigestion(bst.Unit):
    """
    Create an anaerobic digestion unit operation. The model is based on 
    stoichiometric reactions and a specified fraction of water evaporated.
    
    Parameters
    ----------
    reactions : ReactionSet, optional
        Anaerobic digestion reactions. Default assumes 91% of organic components 
        is destroyed, of which 86% is converted to biogas and 5% is converted to 
        sludge.	The biogas is assumed to be 51% CH4 and 49% CO2 on a dry molar 
        basis.
    sludge_split : Array, optional
        Split between wastewater and sludge.
    ins : 
        * [0] Wastewater
    outs : 
        * [0] Biogas
        * [1] Wastewater
        * [2] Sludge
    
    """
    purchase_cost = installation_cost = 0
    _N_ins = 1
    _N_outs = 3
    def _init(self, reactions=None, sludge_split=None):
        chemicals = self.chemicals
        if not reactions:
            digestables = get_digestable_organic_chemicals(chemicals)
            reactions = anaerobic_digestion_reactions(digestables, chemicals.WWTsludge.MW, thermo=self.thermo)
        self.reactions = reactions
        if sludge_split is None:
            sludge_split = dict(
                Water=0.07087,
                Ethanol=0.0625,
                Furfural=0.06667,
                Glycerol=0.07377,
                LacticAcid=0.07084,
                SuccinicAcid=0.07377,
                HNO3=0.0678,
                Denaturant=0.07377,
                DAP=0.0678,
                AmmoniumAcetate=0.07084,
                AmmoniumSulfate=0.0678,
                H2SO4=0.0678,
                NaNO3=0.0678,
                Oil=0.07377,
                HMF=0.06667,
                NH3=0.07048,
                Glucose=0.06667,
                Xylose=0.07609,
                Sucrose=0.06915,
                Mannose=0.06915,
                Galactose=0.06915,
                Arabinose=0.06915,
                Extract=0.07084,
                Tar=0.7473,
                CaO=0.7473,
                Ash=0.7473,
                NaOH=0.0678,
                Lignin=0.744,
                SolubleLignin=0.07084,
                GlucoseOligomer=0.07143,
                GalactoseOligomer=0.07143,
                MannoseOligomer=0.07143,
                XyloseOligomer=0.07143,
                ArabinoseOligomer=0.07143,
                Z_mobilis=0.7438,
                T_reesei=0.7438,
                Cellulose=0.76,
                Protein=0.7391,
                Enzyme=0.7391,
                Xylan=0.75,
                Xylitol=0.07377,
                Cellobiose=0.06915,
                DenaturedEnzyme=0.7391,
                Arabinan=1,
                Mannan=1,
                Galactan=1,
                WWTsludge=0.7438,
                Cellulase=0.07084
            )
            remove_undefined_chemicals(sludge_split, chemicals)
            default_chemical_dict(sludge_split, chemicals, 0.07087, 0.07087, 0.744)
        self.sludge_split = chemicals.isplit(sludge_split)
        self._load_components()
    
    def _load_components(self):
        self.multi_stream = tmo.MultiStream(thermo=self.thermo)
    
    def _run(self):
        feed, = self.ins
        biogas, waste, sludge = self.outs
        biogas.phase = 'g'
        biogas.T = waste.T = sludge.T = 35+273.15
        sludge.copy_flow(feed)
        self.reactions(sludge)
        self.multi_stream.empty()
        self.multi_stream.copy_flow(sludge)
        self.multi_stream['g'].receive_vent(self.multi_stream['l'], energy_balance=False)
        biogas.mol[:] = self.multi_stream.imol['g']
        liquid_mol = self.multi_stream.imol['l']
        sludge.mol[:] = liquid_mol * self.sludge_split.data
        waste.mol[:] = liquid_mol - sludge.mol
        biogas.receive_vent(waste)
        

class AerobicDigestion(bst.Unit):
    """
    Create an aerobic digestion unit operation. Model is based on 
    stoichiometric reactions and a specified fraction of water evaporated.
    
    Parameters
    ----------
    ins : 
        * [0] Wastewater    
        * [1] Air    
        * [2] Caustic    
    outs : 
        * [0] Vent    
        * [1] Treated wastewater
    reactions : ReactionSet, optional
        Aerobic digestion reactions. Defaults assume 96% of remaining soluble 
        organic matter is removed after aerobic digestion, of which 74% is 
        converted to water and CO2 and 22% to cell mass.
    evaporation : float, optional
        Fraction of water evaporated. Defaults to 0.0113.
        
    """    
    _N_ins = 3
    _N_outs = 2
    purchase_cost = installation_cost = 0
    
    def _init(self, reactions=None, evaporation=0.0113):
        if not reactions:
            chemicals = self.chemicals
            digestables = get_digestable_organic_chemicals(self.chemicals)
            reactions = aerobic_digestion_reactions(digestables, chemicals.WWTsludge.MW, thermo=self.thermo)        
        self.reactions = reactions
        self.evaporation = evaporation
    
    def _run(self):
        waste, air, caustic = self._ins
        vent, water = self.outs
        vent.phase = 'g'
        water.copy_like(waste)
        water.mol[:] += caustic.mol
        self.reactions.force_reaction(water)
        O2 = -water.imass['O2']
        N2 = 0.78 / 0.22 * O2
        air.empty()
        air.imass['O2', 'N2'] += [1.2 * O2, 1.2 * N2]
        water.imass['O2'] += 1.2 * O2
        water.imass['N2'] += 1.2 * N2
        vent.copy_flow(water, ('CO2', 'O2', 'N2'))
        water_index = self.chemicals.index('7732-18-5')
        vent.mol[water_index] = water.mol[water_index] * self.evaporation
        water.mol[:] -= vent.mol
        

class SludgeCentrifuge(bst.SolidsSeparator):
    """
    Create a centrifuge to separate sludge. The model is based on 
    component splits.
    
    Parameters
    ----------
    ins : 
        Inlet fluid to be split.
    outs : 
        * [0] Liquid
        * [1] Sludge
    split : Defaults to Should be one of the following
        * [float] The fraction of net feed in the 0th outlet stream
        * [array_like] Componentwise split of feed to 0th outlet stream
        * [dict] ID-split pairs of feed to 0th outlet stream
    order=None : Iterable[str], defaults to biosteam.settings.chemicals.IDs
        Chemical order of split.
    moisture_content : float, optional
        Moisture content of sludge. Defaults to 0.79 based on stream 623 in [1]_
        (or 20% for insolubles).
    
    """
    purchase_cost = installation_cost = 0
    def _init(self, split=None, order=None, moisture_content=None, strict_moisture_content=None):
        chemicals = self.chemicals
        ID_water = chemicals['7732-18-5'].ID # Water must be defined
        if split is None:
            split = dict(
                Furfural=1,
                Glycerol=0.8889,
                LacticAcid=0.935,
                SuccinicAcid=0.8889,
                HNO3=0.9344,
                Denaturant=0.8889,
                DAP=0.9344,
                AmmoniumAcetate=0.935,
                AmmoniumSulfate=0.9344,
                H2SO4=0.935,
                NaNO3=0.9344,
                Oil=0.8889,
                HMF=1,
                NH3=0.9388,
                H2S=0.9394,
                SO2=0.9394,
                CO2=0.9333,
                NO2=0.9394,
                NO=0.9394,
                CO=0.9394,
                Glucose=1,
                Xylose=1,
                Sucrose=0.9286,
                Mannose=0.9286,
                Galactose=0.9286,
                Arabinose=0.9286,
                Extract=0.935,
                Tar=0.05155,
                CaO=0.05155,
                Ash=0.05155,
                NaOH=0.9344,
                Lignin=0.04943,
                SolubleLignin=0.935,
                GlucoseOligomer=0.9,
                GalactoseOligomer=0.9,
                MannoseOligomer=0.9,
                XyloseOligomer=0.9,
                ArabinoseOligomer=0.9,
                Z_mobilis=0.04991,
                T_reesei=0.04991,
                Cellulose=0.03846,
                Protein=0.05455,
                Enzyme=0.05455,
                Xylitol=0.8889,
                Cellobiose=0.9286,
                DenaturedEnzyme=0.05455,
                WWTsludge=0.04991,
                Cellulase=0.935
            )
            remove_undefined_chemicals(split, chemicals)
            default_chemical_dict(split, chemicals, 0.9394, 0.9286, 0.04991)
            split.pop(ID_water) # Remove water from split
        if ID_water not in split and moisture_content is None:
            # Only set moisture content if water split is not given
            moisture_content = 0.79
        bst.SolidsSeparator._init(
            self, split=split, order=order,
            moisture_content=moisture_content, strict_moisture_content=strict_moisture_content
        )

    def _run(self):
        ins = self.ins
        retentate, permeate = self.outs # Filtrate, solids
        separations.mix_and_split(ins, retentate, permeate, self.split)
        separations.adjust_moisture_content(permeate, retentate, self.moisture_content, None, self.strict_moisture_content)
            
        
# TODO: Split values seem arbitrary in NREL 2011 model, perhaps work on a better model
class MembraneBioreactor(bst.Splitter):
    """
    Create a membrane bioreactor to clarify sludge. The model is based on 
    component splits.
    
    Parameters 
    ----------
    ins : 
        Inlet fluid to be split.
    outs : 
        * [0] Liquid
        * [1] Sludge
    split : Defaults to Should be one of the following
        * [float] The fraction of net feed in the 0th outlet stream
        * [array_like] Componentwise split of feed to 0th outlet stream
        * [dict] ID-split pairs of feed to 0th outlet stream
    order=None : Iterable[str], defaults to biosteam.settings.chemicals.IDs
        Chemical order of split.
    
    """
    purchase_cost = installation_cost = 0
    def _init(self, split=None, order=None):
        if split is None:
            chemicals = self.chemicals
            split = dict(
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
            remove_undefined_chemicals(split, chemicals)
            default_chemical_dict(split, chemicals, 0.15, 0.125, 0.145)
        bst.Splitter._init(self, split=split, order=order)
        
class ReverseOsmosis(bst.Unit):
    """
    Create a reverse osmosis unit operation for recovering water from brine.
    The model is based on a fraction of water recovered.
    
    Parameters
    ----------
    ins : 
        Inlet fluid to be split.
    outs : 
        * [0] Filtered water
        * [1] Brine
    water_recovery : float, optional
        Water recovered to 0th stream. Defaults to 0.987
    
    """
    _N_ins = 1
    _N_outs = 2
    def _init(self, water_recovery=0.987):
        self.water_recovery = water_recovery
    
    @property
    def RO_treated_water(self):
        return self.outs[0]
    
    def _run(self):
        feed, = self.ins
        water, brine = self.outs
        water.copy_thermal_condition(feed)
        brine.copy_like(feed)
        water_index = self.chemicals.index('7732-18-5')
        water_flow = brine.mol[water_index]
        water_recovered = self.water_recovery * water_flow
        water.mol[water_index] = water_recovered
        brine.mol[water_index] = water_flow - water_recovered

def _create_cellulosic_ethanol_chemicals():
    from biorefineries.cellulosic import create_cellulosic_ethanol_chemicals
    return create_cellulosic_ethanol_chemicals()

@bst.SystemFactory(
    ID='wastewater_treatment_sys',
    outs=[dict(ID='biogas'),
          dict(ID='sludge'),
          dict(ID='RO_treated_water'),
          dict(ID='brine')],
    fixed_ins_size=False,
    fthermo=_create_cellulosic_ethanol_chemicals,
)
def create_conventional_wastewater_treatment_system(ins, outs, 
        NaOH_price=None, autopopulate=None
    ):
    """
    Return a system for wastewater treatment as described in Humbird et al. [1]_
    The system includes anaerobic and aerobic  digestion reactors,
    a membrane bioreactor, a sludge centrifuge, and a reverse osmosis unit.
         
    Parameters
    ----------
    ins : 
        Wastewater streams (without solids). Defaults to all product streams
        at run time that are not sold and cannot generate energy through combustion
        (i.e. streams that have no sink, no price, and a LHV less that 1 kJ / g).
    outs : 
        * [0] biogas
        * [1] sludge
        * [2] RO_treated_water
        * [3] brine 
    NaOH_price : float, optional
        Price of NaOH in USD/kg. The default is 0.07476.
    autopopulate : bool, optional
        Whether to automatically add wastewater streams.

    Examples
    --------
    >>> from biosteam import Stream, create_conventional_wastewater_treatment_system, settings
    >>> settings.set_thermo(create_conventional_wastewater_treatment_system.fthermo())
    >>> feed = Stream(
    ...     ID='wastewater', 
    ...     Water=2.634e+04, 
    ...     Ethanol=0.07225, 
    ...     AceticAcid=24.67, 
    ...     Furfural=6.206, 
    ...     Glycerol=1.784, 
    ...     LacticAcid=17.7, 
    ...     SuccinicAcid=3.472, 
    ...     DAP=1.001, 
    ...     AmmoniumSulfate=17.63, 
    ...     HMF=2.366, 
    ...     Glucose=2.816, 
    ...     Xylose=6.953, 
    ...     Arabinose=12.78, 
    ...     Extract=65.98, 
    ...     Ash=83.52, 
    ...     Lignin=1.659, 
    ...     SolubleLignin=4.202, 
    ...     GlucoseOligomer=6.796, 
    ...     GalactoseOligomer=0.01718, 
    ...     MannoseOligomer=0.009008, 
    ...     XyloseOligomer=2.878, 
    ...     ArabinoseOligomer=0.3508, 
    ...     Z_mobilis=0.6668, 
    ...     Protein=2.569, 
    ...     Glucan=0.1555, 
    ...     Xylan=0.06121, 
    ...     Xylitol=4.88, 
    ...     Cellobiose=0.9419, 
    ...     Arabinan=0.02242, 
    ...     Mannan=0.06448, 
    ...     Galactan=0.01504, 
    ...     Cellulase=25.4, 
    ...     units='kmol/hr'
    ... )
    >>> wwt_sys = create_conventional_wastewater_treatment_system(ins=feed)
    >>> wwt_sys.simulate()
    >>> wwt_sys.show('cwt100')
    System: wastewater_treatment_sys
    Highest convergence error among components in recycle
    stream M604-0 after 16 loops:
    - flow rate   1.13e+03 kmol/hr (0.84%)
    - temperature 6.49e-05 K (2.1e-05%)
    ins...
    [0] wastewater
        phase: 'l', T: 298.15 K, P: 101325 Pa
        composition (%): Water              94.6
                         Ethanol            0.000664
                         AceticAcid         0.295
                         Furfural           0.119
                         Glycerol           0.0328
                         LacticAcid         0.318
                         SuccinicAcid       0.0818
                         DAP                0.0264
                         AmmoniumSulfate    0.465
                         HMF                0.0595
                         Glucose            0.101
                         Xylose             0.208
                         Arabinose          0.383
                         Extract            2.37
                         Ash                0.0167
                         Lignin             0.0503
                         SolubleLignin      0.128
                         GlucoseOligomer    0.244
                         GalactoseOligomer  0.000617
                         MannoseOligomer    0.000324
                         XyloseOligomer     0.0862
                         ArabinoseOligomer  0.0105
                         Z_mobilis          0.00328
                         Protein            0.0117
                         Glucan             0.00503
                         Xylan              0.00161
                         Xylitol            0.148
                         Cellobiose         0.0643
                         Arabinan           0.000591
                         Mannan             0.00209
                         Galactan           0.000486
                         Cellulase          0.122
                         -----------------  5.01e+05 kg/hr
    outs...
    [0] biogas
        phase: 'g', T: 307.76 K, P: 101325 Pa
        composition (%): Water         3.18
                         Ethanol       2.46e-05
                         AceticAcid    0.0019
                         Furfural      0.00296
                         Glycerol      1.43e-08
                         LacticAcid    8.25e-07
                         SuccinicAcid  2.88e-09
                         CH4           26.6
                         CO2           70.2
                         ------------  2.13e+04 kg/hr
    [1] sludge
        phase: 'l', T: 307.88 K, P: 101325 Pa
        composition (%): Water              46
                         Ethanol            5.26e-05
                         AceticAcid         0.0267
                         Glycerol           0.0048
                         LacticAcid         0.0262
                         SuccinicAcid       0.012
                         DAP                0.0976
                         AmmoniumSulfate    1.72
                         SO2                0.000601
                         Arabinose          0.0338
                         Extract            0.195
                         Ash                2.56
                         NaOH               1.58
                         Lignin             7.54
                         SolubleLignin      0.0105
                         GlucoseOligomer    0.0312
                         GalactoseOligomer  7.88e-05
                         MannoseOligomer    4.13e-05
                         XyloseOligomer     0.011
                         ArabinoseOligomer  0.00134
                         Z_mobilis          0.0406
                         Protein            0.143
                         Glucan             0.0229
                         Xylan              0.0226
                         Xylitol            0.0217
                         Cellobiose         0.00568
                         Arabinan           0.0103
                         Mannan             0.0363
                         Galactan           0.00847
                         WWTsludge          39.8
                         Cellulase          0.01
                         -----------------  2.57e+03 kg/hr
    [2] RO_treated_water
        phase: 'l', T: 307.79 K, P: 101325 Pa
        composition (%): Water  100
                         -----  4.16e+05 kg/hr
    [3] brine
        phase: 'l', T: 307.79 K, P: 101325 Pa
        composition (%): Water              48.7
                         Ethanol            1.35e-05
                         AceticAcid         0.00609
                         Furfural           0.00244
                         Glycerol           0.000675
                         LacticAcid         0.00763
                         SuccinicAcid       0.00169
                         DAP                1.11
                         AmmoniumSulfate    19.6
                         HMF                0.00124
                         SO2                0.01
                         Glucose            0.00244
                         Xylose             0.00861
                         Arabinose          0.00791
                         Extract            0.0568
                         Ash                0.267
                         NaOH               23.5
                         Lignin             0.519
                         SolubleLignin      0.00306
                         GlucoseOligomer    0.00576
                         GalactoseOligomer  1.46e-05
                         MannoseOligomer    7.63e-06
                         XyloseOligomer     0.00203
                         ArabinoseOligomer  0.000248
                         Z_mobilis          1.99e-05
                         Protein            7.32e-05
                         Glucan             0.206
                         Xylan              0.0657
                         Xylitol            0.00305
                         Cellobiose         0.00133
                         Arabinan           0.0237
                         Mannan             0.0837
                         Galactan           0.0195
                         WWTsludge          5.9
                         Cellulase          0.00292
                         -----------------  1.12e+04 kg/hr
    
    """
    biogas, sludge, RO_treated_water, brine = outs
    RO_treated_water.register_alias('recycled_water')
    if NaOH_price is None: NaOH_price = 0.07476
    air = bst.Stream('air_lagoon', O2=51061, N2=168162, phase='g', units='kg/hr')
    caustic = bst.Stream('caustic', Water=2252, NaOH=2252,
                     units='kg/hr', price=NaOH_price)
    wastewater_mixer = bst.Mixer('M601', ins or [])
    wastewater_mixer.autopopulate = False if autopopulate is None else autopopulate
    
    @wastewater_mixer.add_specification(run=True)
    def autopopulate_waste_streams():
        if wastewater_mixer.autopopulate and not wastewater_mixer.ins: 
            sys = wastewater_mixer.system
            streams = bst.FreeProductStreams(sys.streams)
            wastewater_mixer.ins.extend(streams.noncombustible_slurries)
            for i in sys.facilities:
                if isinstance(i, bst.BlowdownMixer): wastewater_mixer.ins.append(i.outs[0])
            
    WWTC = WastewaterSystemCost('WWTC', wastewater_mixer-0)
    anaerobic_digestion = AnaerobicDigestion('R601', WWTC-0, (biogas, '', ''))
    recycled_sludge_mixer = bst.Mixer('M602', ('', anaerobic_digestion-1))
    
    caustic_over_waste = caustic.imol['Water', 'NaOH'] / 2544301
    air_over_waste = air.imol['O2', 'N2'] / 2544301
    air.mol[:] = 0.
    waste = recycled_sludge_mixer-0
    aerobic_digestion = AerobicDigestion('R602', (waste, air, caustic),
                                         outs=('evaporated_water', ''))
    
    @aerobic_digestion.add_specification
    def update_aerobic_input_streams():
        waste, air, caustic = aerobic_digestion.ins
        F_mass_waste = waste.F_mass
        caustic.imol['Water', 'NaOH'] = F_mass_waste * caustic_over_waste
        air.imol['O2', 'N2'] = F_mass_waste * air_over_waste
        aerobic_digestion._run()
    
    membrane_bioreactor = MembraneBioreactor('S601', aerobic_digestion-1)
    sludge_splitter = bst.Splitter('S602', membrane_bioreactor-1, split=0.96)
    fresh_sludge_mixer = bst.Mixer('M603', (anaerobic_digestion-2, sludge_splitter-1))
    sludge_centrifuge = SludgeCentrifuge('S603', fresh_sludge_mixer-0, outs=('', sludge))
    bst.Mixer('M604', [sludge_splitter-0, sludge_centrifuge-0], 0-recycled_sludge_mixer)
    ReverseOsmosis('S604', membrane_bioreactor-0, outs=(RO_treated_water, brine))


create_wastewater_treatment_system = create_conventional_wastewater_treatment_system