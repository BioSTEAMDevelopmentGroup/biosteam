# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
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
.. autodata:: biosteam.units.wastewater.non_digestables
    
Unit operations
---------------
.. autoclass:: biosteam.units.wastewater.AnaerobicDigestion 
.. autoclass:: biosteam.units.wastewater.AerobicDigestion
.. autoclass:: biosteam.units.wastewater.ReverseOsmosis
.. autoclass:: biosteam.units.wastewater.WastewaterSystemCost

Utilities
---------
.. autofunction:: biosteam.units.wastewater.get_digestable_organic_chemicals

System factories
----------------

.. autofunction:: biosteam.units.wastewater.create_wastewater_treatment_units
.. autodata:: biosteam.units.wastewater.create_wastewater_treatment_system

References
----------
.. [1] Humbird, D., Davis, R., Tao, L., Kinchin, C., Hsu, D., Aden, A.,
    Dudgeon, D. (2011). Process Design and Economics for Biochemical 
    Conversion of Lignocellulosic Biomass to Ethanol: Dilute-Acid 
    Pretreatment and Enzymatic Hydrolysis of Corn Stover
    (No. NREL/TP-5100-47764, 1013269). https://doi.org/10.2172/1013269

"""
from .._unit import Unit
from .decorators import cost
from .splitting import Splitter
from biosteam.utils import remove_undefined_chemicals, default_chemical_dict
import biosteam as bst
import thermosteam as tmo
from thermosteam import (
    utils,
    settings,
    Reaction as Rxn,
    ParallelReaction as PRxn,
)

__all__ = (
    'AnaerobicDigestion', 
    'AerobicDigestion', 
    'ReverseOsmosis',
    'WastewaterSystemCost',
    'get_digestable_organic_chemicals',
    'create_wastewater_treatment_system',
)


# %% Functional utilities

#: list[str] IDs for non-digestable components in wastewater
non_digestables = ['WWTsludge', 'Cellulose', 'Xylan', 'CellulaseNutrients',
                   'Mannan', 'Lignin', 'Galactan', 'Glucan', 'Acetate',
                   'Biomass', 'Arabinan', 'Tar', 'CO', 'CO2', 'CH4']

def get_digestable_organic_chemicals(chemicals):
    """
    Return a list of digestable organic chemical IDs.

    Parameters
    ----------
    chemicals : :class:`~thermosteam.Chemicals`
        Digestable organic chemicals will be retrieve from this parameter.

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
        f = MW_sludge / chemical.MW
        reactant = chemical.ID
        if not isvalid(reactant): reactant = parsable_name(reactant)
        return Rxn(f"{f}{reactant} -> WWTsludge", reactant, X_growth, chemicals=thermo.chemicals)
    return PRxn([i.get_combustion_reaction(conversion=X_combustion) + growth(i)
                 for i in chemicals])

# %% Unit operations

@cost('Flow rate', 'Wastewater system', units='kg/hr', CE=551,
      cost=50280080., n=0.6, BM=1, kW=7139/1.05, S=393100)
class WastewaterSystemCost(Unit): 
    """
    Create a unit that estimates the capital cost and electricity demand
    of a wastewater treatment system.
    
    Parameters
    ----------
    ins : stream
        Wastewater.
        
    """


class AnaerobicDigestion(Unit):
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
    ins : stream sequence
        * [0] Wastewater
        * [1] Cool well water
    outs : stream sequence
        * [0] Biogas
        * [1] Wastewater
        * [2] Sludge
        * [3] Hot well water
    
    """
    purchase_cost = installation_cost = 0
    _N_ins = 2
    _N_outs = 4
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 reactions=None, sludge_split=None):
        Unit.__init__(self, ID, ins, outs, thermo)
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
        feed, cool_water = self.ins
        biogas, waste, sludge, hot_water = self.outs
        biogas.phase = 'g'
        hot_water.link_with(cool_water, TP=False)
        biogas.T = waste.T = sludge.T = T = 35+273.15
        hot_water.T = feed.T - 5
        H_at_35C = feed.thermo.mixture.H('l', feed.mol, T, 101325)
        cool_water.mol[:] *= (feed.H - H_at_35C)/(hot_water.H - cool_water.H)
        sludge.copy_flow(feed)
        self.reactions(sludge)
        self.multi_stream.copy_flow(sludge)
        self.multi_stream.vle(P=101325, H=self.multi_stream.H)
        biogas.mol[:] = self.multi_stream.imol['g']
        liquid_mol = self.multi_stream.imol['l']
        sludge.mol[:] = liquid_mol * self.sludge_split.data
        waste.mol[:] = liquid_mol - sludge.mol
        biogas.receive_vent(waste)
        

class AerobicDigestion(Unit):
    """
    Create an aerobic digestion unit operation. Model is based on 
    stoichiometric reactions and a specified fraction of water evaporated.
    
    Parameters
    ----------
    ins : stream sequence
        * [0] Wastewater    
        * [1] Air    
        * [2] Caustic    
    outs : stream sequence
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
    
    def __init__(self, ID='', ins=None, outs=(), *, reactions=None, evaporation=0.0113):
        Unit.__init__(self, ID, ins, outs)
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
        water.mol[:] += air.mol + caustic.mol
        self.reactions.force_reaction(water)
        O2 = float(water.imass['O2'])
        if O2 < 0:
            N2 = 0.78 / 0.22 * O2
            air.imass['O2', 'N2'] += [O2, N2]
            water.imol['O2'] = 0.
            water.imass['N2'] += N2
        vent.copy_flow(water, ('CO2', 'O2', 'N2'))
        water_index = self.chemicals.index('7732-18-5')
        vent.mol[water_index] = water.mol[water_index] * self.evaporation
        water.mol[:] -= vent.mol
        

# TODO: Use moisture content
class SludgeCentrifuge(Splitter):
    """
    Create a centrifuge to separate sludge. The model is based on 
    component splits.
    
    Parameters
    ----------
    ins : stream
        Inlet fluid to be split.
    outs : stream sequence
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
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *, 
                 split=None, order=None):
        self._load_thermo(thermo)
        if split is None:
            chemicals = self.chemicals
            split = dict(
                Water=0.934,
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
        Splitter.__init__(self, ID, ins, outs, thermo, split=split, order=order)
        
# TODO: Split values seem arbitrary in NREL 2011 model, perhaps work on a better model
class MembraneBioreactor(Splitter):
    """
    Create a membrane bioreactor to clarify sludge. The model is based on 
    component splits.
    
    Parameters 
    ----------
    ins : stream
        Inlet fluid to be split.
    outs : stream sequence
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
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *, 
                 split=None, order=None):
        self._load_thermo(thermo)
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
        Splitter.__init__(self, ID, ins, outs, thermo, split=split, order=order)
        
class ReverseOsmosis(Unit):
    """
    Create a reverse osmosis unit operation for recovering water from brine.
    The model is based on a fraction of water recovered.
    
    Parameters
    ----------
    ins : stream
        Inlet fluid to be split.
    outs : stream sequence
        * [0] Filtered water
        * [1] Brine
    water_recovery : float, optional
        Water recovered to 0th stream. Defaults to 0.987
    
    """
    _N_ins = 1
    _N_outs = 2
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 water_recovery=0.987):
        Unit.__init__(self, ID, ins, outs, thermo)
        self.water_recovery = water_recovery
        
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


def create_wastewater_treatment_units(ins, outs, NaOH_price=0.15):
    """
    Create units for wastewater treatment, including anaerobic and aerobic 
    digestion reactors, a membrane bioreactor, a sludge centrifuge, 
    and a reverse osmosis unit.
    
    Parameters
    ----------
    ins : streams
        Wastewater streams (without solids).
    outs : stream sequence
        * [0] methane
        * [1] sludge
        * [2] treated_water
        * [3] waste_brine 
    NaOH_price : float, optional
        Price of NaOH in USD/kg. The default is 0.15.
    
    """
    methane, sludge, treated_water, waste_brine = outs
    well_water = bst.Stream('well_water', Water=1, T=15+273.15)
    air = bst.Stream('air_lagoon', O2=51061, N2=168162, phase='g', units='kg/hr')
    caustic = bst.Stream('caustic', Water=2252, NaOH=2252,
                     units='kg/hr', price=NaOH_price)
    
    wastewater_mixer = bst.Mixer('M601', ins)
    WWTC = WastewaterSystemCost('WWTC', wastewater_mixer-0)
    anaerobic_digestion = AnaerobicDigestion('R601', (WWTC-0, well_water), (methane, '', '', ''))
    recycled_sludge_mixer = bst.Mixer('M602', (anaerobic_digestion-1, None))
    
    caustic_over_waste = 2 * caustic.imol['Water', 'NaOH'] / 2544301
    air_over_waste = air.imol['O2', 'N2'] / 2544301
    air.mol[:] = 0.
    waste = recycled_sludge_mixer-0
    def update_aerobic_input_streams():
        waste, air, caustic = aerobic_digestion.ins
        F_mass_waste = waste.F_mass
        caustic.imol['Water', 'NaOH'] = F_mass_waste * caustic_over_waste
        air.imol['O2', 'N2'] = F_mass_waste * air_over_waste
        aerobic_digestion._run()
    
    aerobic_digestion = AerobicDigestion('R602', (waste, air, caustic),
                                         outs=('evaporated_water', ''))
    aerobic_digestion.specification = update_aerobic_input_streams
    membrane_bioreactor = MembraneBioreactor('S601', aerobic_digestion-1)
    sludge_splitter = bst.Splitter('S602', membrane_bioreactor-1, split=0.96)
    fresh_sludge_mixer = bst.Mixer('M603', (anaerobic_digestion-2, sludge_splitter-1))
    sludge_centrifuge = SludgeCentrifuge('S603', fresh_sludge_mixer-0, outs=('', sludge))
    sludge_centrifuge-0-1-recycled_sludge_mixer
    reverse_osmosis = ReverseOsmosis('S604', membrane_bioreactor-0,
                                     outs=(treated_water, waste_brine))

create_wastewater_treatment_system = bst.SystemFactory(
    f=create_wastewater_treatment_units,
    ID='wastewater_treatment_sys',
    outs=[dict(ID='methane'),
          dict(ID='sludge'),
          dict(ID='treated_water'),
          dict(ID='waste_brine')],
    fixed_ins_size=False,
)
"""
Return a system for wastewater treatment, which includes anaerobic and aerobic 
digestion reactors, a membrane bioreactor, a sludge centrifuge, and a reverse 
osmosis unit.
 
Parameters
----------
ID : str, optional
    Defaults to 'wastewater_treatment_sys'.
ins : streams, optional
    Wastewater streams (without solids).
outs : stream sequence
    * [0] methane
    * [1] sludge
    * [2] treated_water
    * [3] waste_brine 
mockup : bool, optional
    Whether to create a mock system.
area : int, optional
    If given, IDs of all units will follow the area naming convention as
    explained in :func:`~biosteam.process_tools.utils.rename_unit`.
udct : bool, optional
    Whether to also return fictionary of units with their original IDs as keys.
NaOH_price : float, optional
    Price of NaOH in USD/kg. The default is 0.15.

"""