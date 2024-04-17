# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-2024, Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

__all__ = (
    'create_high_rate_wastewater_treatment_system',    
)

import thermosteam as tmo, biosteam as bst
cost = bst.decorators.cost
CEPCI = bst.design_tools.CEPCI_by_year
from . import (
    append_wwt_chemicals, default_insolubles,
    InternalCirculationRx, AnMBR, PolishingFilter, BeltThickener, SludgeCentrifuge,
    get_combustion_energy, prices, GWP_CFs,
)

_mgd_to_cmh = 157.7255 # auom('gallon').conversion_factor('m3')*1e6/24
_gpm_to_cmh = 0.2271 # auom('gallon').conversion_factor('m3')*60
_Gcal_to_kJ = 4184000 # auom('kcal').conversion_factor('kJ')*1e6 # (also MMkcal/hr)
_kW_to_kJhr = 3600 # auom('kW').conversion_factor('kJ/hr')

Rxn = tmo.reaction.Reaction
ParallelRxn = tmo.reaction.ParallelReaction

__all__ = (
    'BiogasUpgrading',
    'create_high_rate_wastewater_treatment_system',
    'CHP',
    'ReverseOsmosis',
    'Skipped',
    )

@cost(basis='Volumetric flow', ID='Reactor', units='m3/hr',
      # 2.7 in million gallons per day (MGD)
      cost=2450000, S=2.7*_mgd_to_cmh, CE=CEPCI[2012], n=1, BM=1.8)
@cost(basis='Volumetric flow', ID='Evaporator', units='m3/hr',
      # 2.7 in million gallons per day (MGD)
       kW=1103.636, cost=5000000, S=2.7*_mgd_to_cmh, CE=CEPCI[2012], n=0.6, BM=1.6)
class ReverseOsmosis(bst.Unit):
    _N_ins = 1
    _N_outs = 2
    _units = {'Volumetric flow': 'm3/hr'}

    @property
    def RO_treated_water(self):
        return self.outs[0]
    
    @property
    def brine(self):
        return self.outs[1]

    def _run(self):
        influent = self.ins[0]
        water, brine = self.outs

        self.design_results['Volumetric flow'] = self.F_vol_in

        # Based on stream 626 and 627 in ref [1]
        water.imass['Water'] = 376324/(376324+4967) * influent.imass['Water']
        brine.mol = influent.mol - water.mol
        water.T = brine.T = influent.T


class Skipped(bst.Unit):
    """
    Copy ins[`main_in`] as ins[`main_out`].
    Can be also used to calculate the cost of wastewater treatment
    by clearing all WWT-related units and streams.

    Parameters
    ----------
    main_in : int
        Which influent will be copied to `main_out`.
    main_out : int
        Which effluent will be copied from `main_in`.
    wwt_units : Iterable
        Collection of units whose costs will be cleared when `clear_wwt` is True.
        ins and outs of the units will be emptied if its price isn't 0.
    wwt_streams : Iterable
        Collection of streams which will be emptied when `clear_wwt` is True.
        Usually should at least includes the biogas and sludge stream.
    clear_wwt : bool
        Whether to clear the costs of and select streams associated with
        `wwt_units`.
    """
    _ins_size_is_fixed = False
    _outs_size_is_fixed = False
    cache_dct = {} # to save some computation effort

    def _init(self, main_in=0, main_out=0, wwt_units=[], wwt_streams=[],
                 clear_wwt=False):
        self.main_in = main_in
        self.main_out = main_out
        self.wwt_units = wwt_units
        self.wwt_streams = wwt_streams
        self.clear_wwt = clear_wwt

    def _run(self):
        self.outs[self.main_out].copy_like(self.ins[self.main_in])

    def _cost(self):
        if not self.clear_wwt: return
        for u in self.wwt_units:
            u.baseline_purchase_costs.clear()
            u.installed_costs.clear()
            for hu in u.heat_utilities: hu.empty()
            u.power_utility(0)
            for s in u.ins+u.outs:
                if s.price: s.empty()
            u._utility_cost = 0.
        for s in self.wwt_streams: s.empty()

    @property
    def utility_cost(self): return 0.


class CHP(bst.Unit):
    """
    Used to estimate the cost of producing electricity as in [6]_.

    Parameters
    ----------
    eff : float
        Combined efficiency for combustion and power generation.
    unit_CAPEX : float
        Capital cost of the CHP per kW of power generated, $/kW.

    """
    _ins_size_is_fixed = False
    _F_BM_default = {'CHP': 1.}

    def _init(self, 
              eff=0.3375, # average of 40.5%, 27%, 36%, and 31.5% for the four types in ref [1]
              unit_CAPEX=1225):
        self.eff = eff
        self.unit_CAPEX = unit_CAPEX

    def _run(self):
        mixed, = self.outs
        mixed.mix_from(self.ins)
        H_net = get_combustion_energy(mixed, 1)
        self.H_for_power = H_net * self.eff


    def _design(self):
        kW = self.H_for_power / 3600 # kJ/hr to kW
        self.baseline_purchase_costs['CHP'] = kW*self.unit_CAPEX
        self.power_utility(-kW)


class BiogasUpgrading(bst.Unit):
    """
    Upgrade the biogas to renewable natural gas (RNG).
    Note that the second influent is a dummy stream to calculate the upgrading cost.

    Parameters
    ----------
    ratio : float
        How much of the incoming biogas will be upgraded to RNG.
    loss : float
        CH4 loss during upgrading.
    unit_upgrading_cost : float
        Unit cost for upgrading, in $/MMBtu.
    unit_upgrading_GWP : float
        Unit 100-yr global warming potential (GWP) for biogas upgrading, in kg CO2/MMBtu.

    References
    ----------
    Yang et al., Cost and Life-Cycle Greenhouse Gas Implications of
    Integrating Biogas Upgrading and Carbon Capture Technologies in Cellulosic Biorefineries.
    Environ. Sci. Technol. 2020.
    https://doi.org/10.1021/acs.est.0c02816.
    
    IEA. Outlook for Biogas and Biomethane: Prospects for Organic Growth; IEA: Paris, 2020.
    https://www.iea.org/reports/outlook-for-biogas-and-biomethane-prospects-for-organic-growth
    
    Rai et al., Comparative Life Cycle Evaluation of the Global Warming Potential (GWP)
    Impacts of Renewable Natural Gas Production Pathways. Environ. Sci. Technol. 2022.
    https://doi.org/10.1021/acs.est.2c00093.

    """

    # Default upgrading cost is calculated by
    # 0.18*0.914/35.3147*1055.056
    # where $0.18/NM is from ref 1,
    # 0.914 ft3/MJ is from U.S. EIA https://www.eia.gov/energyexplained/units-and-calculators/energy-conversion-calculators.php#natgascalc
    # 35.3147 ft3/m3 is auom('m3').conversion_factor('ft3')
    # 1055.056 MJ/MMBtu is auom('MMBtu').conversion_factor('MJ')
    # This cost is close to the general range of $2-4 for a
    # 3.5 million m3/yr biogas plant in ref 2.

    # Default upgrading GWP is from calculated by
    # 8.67/1e3*1055.056
    # where 8.67 g CO2/MJ is the average of upgrading GWP (5, 8, 13)
    # from Section S5 of ref 3

    _N_ins = 2
    _N_outs = 2

    def _init(self,
                 ratio=0, loss=0.05,
                 unit_upgrading_cost=4.92,
                 unit_upgrading_GWP=9.15):
        self.ratio = ratio
        self.loss = loss
        self.unit_upgrading_cost = unit_upgrading_cost
        self.unit_upgrading_GWP = unit_upgrading_GWP
        self.RIN_incentive = prices['RIN']
        # Credits from the displaced fossil natural gas
        self.FNG_price = bst.stream_prices.get('Natural gas', 0) # $/kg
        self.FNG_CF = GWP_CFs['CH4']


    def _run(self):
        biogas, foo = self.ins
        RNG, remained = self.outs
        RNG.empty()
        RNG.phase = remained.phase = 'g'
        RNG.imass['CH4'] = biogas.imass['CH4'] * self.ratio
        remained.mass = biogas.mass - RNG.mass # assume impurities left in the unused biogas
        RNG.mass *= (1 - self.loss) # lost methane not included in the unused biogas
        RNG.price = self.FNG_price + self.RIN_incentive
        RNG.characterization_factors['GWP'] = self.FNG_CF
        foo.copy_like(RNG)
        
        # Upgrading cost/GWP
        if not foo.isempty(): 
            foo.price = self.unit_upgrading_cost/1055.056*(foo.HHV/1e3)/foo.F_mass # HHV in kJ/hr
            foo.characterization_factors['GWP'] = self.unit_upgrading_GWP/1055.056*(foo.HHV/1e3)/foo.F_mass


# %%

# =============================================================================
# System function
# =============================================================================

@bst.SystemFactory(
    ID='wastewater_sys',
    outs=[dict(ID='RNG'), # renewable natural gas
          dict(ID='biogas'),
          dict(ID='sludge'),
          dict(ID='RO_treated_water'),
          dict(ID='brine')],
    fixed_ins_size=False,
    fthermo=append_wwt_chemicals,
)
def create_high_rate_wastewater_treatment_system(
        ins=None, outs=None, process_ID='6',
        flowsheet=None, autopopulate=False,
        skip_IC=False, IC_kwargs={},
        skip_AnMBR=False, AnMBR_kwargs={},
        skip_AeF=False, AeF_kwargs={}
        ):
    """
    Return a system for wastewater treatment (WWT) as described in Li et al. [1]_
    The system includes internal circulation (IC), anaerobic membrane bioreactors (AnMBR), 
    aerobic polishing filter (AeF),
    belt thickener and sludge centrifuge, and a reverse osmosis unit.
    
    Users can choose whether to skip IC
    (makes sense when influent COD is only on the order of grams per liter),
    AnMBR (not recommended unless the COD is very low),
    and/or AeF (not recommended if want to achieve low COD in the effluent).
    
    An optional biogas upgrading unit can be included to upgrade the biogas
    (from IC and AnMBR) as renewable natural gas (RNG) for sale with incentives,
    this is achieved through adjusting `BiogasUpgrading.ratio`.
         
    Parameters
    ----------
    ins : 
        Wastewater streams (without solids). Defaults to all product streams
        at run time that are not sold and cannot generate energy through combustion
        (i.e. streams that have no sink, no price, and a LHV less that 1 kJ / g).
    outs : 
        * [0] RNG
        * [1] biogas
        * [2] sludge
        * [3] RO_treated_water
        * [4] brine
    process_ID : float
        Number of the process.
        E.g., the default `process_ID` is 6,
        then the first mixer of this WWT system will be M601.
    flowsheet : Flowsheet, optional
        If provided, the WWT system will be added to the given flowsheet.
    autopopulate : bool, optional
        Whether to automatically add wastewater streams.
    skip_IC : bool
        Whether to skip the IC unit.
    IC_kwargs : dict
        kwargs to be passed to the IC unit (refer to the doc of the IC unit for details).
    skip_AnMBR : bool
        Whether to skip the AnMBR unit.
    AnMBR_kwargs : dict
        kwargs to be passed to the AnMBR unit (refer to the doc of the AnMBR unit for details).
    skip_AeF : bool
        Whether to skip the AeF unit.
    AeF_kwargs : dict
        kwargs to be passed to the AeF unit (refer to the doc of the AeF unit for details).

    Examples
    --------
    Check for PolishingFilter vent accumulation
    
    >>> from biosteam import Stream, create_high_rate_wastewater_treatment_system, settings
    >>> from biorefineries import cornstover as cs
    >>> settings.set_thermo(cs.create_chemicals())
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
    >>> wwt_sys = create_high_rate_wastewater_treatment_system(ins=feed)
    >>> wwt_sys.simulate()
    >>> wwt_sys.show('cwt100')
    System: wastewater_sys
    Highest convergence error among components in recycle
    stream M603-0 after 6 loops:
    - flow rate   1.57e+01 kmol/hr (0.066%)
    - temperature 1.19e-03 K (0.00039%)
    ins...
    [0] wastewater  
        phase: 'l', T: 298.15 K, P: 101325 Pa
        composition (%): Water              94.7
                         Ethanol            0.000664
                         AceticAcid         0.296
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
                         Lignin             0.0504
                         SolubleLignin      0.128
                         GlucoseOligomer    0.22
                         GalactoseOligomer  0.000556
                         MannoseOligomer    0.000291
                         XyloseOligomer     0.0759
                         ArabinoseOligomer  0.00925
                         Z_mobilis          0.00328
                         Protein            0.0117
                         Glucan             0.00503
                         Xylan              0.00161
                         Xylitol            0.148
                         Cellobiose         0.0643
                         Arabinan           0.000591
                         Mannan             0.00209
                         Galactan           0.000487
                         Cellulase          0.122
                         -----------------  5.01e+05 kg/hr
    outs...
    [0] RNG  
        phase: 'g', T: 298.15 K, P: 101325 Pa
        flow: 0
    [1] biogas  
        phase: 'g', T: 298.15 K, P: 101325 Pa
        composition (%): CH4  26.9
                         H2S  0.0348
                         CO2  73.1
                         ---  4.33e+04 kg/hr
    [2] sludge  
        phase: 'l', T: 307.85 K, P: 101325 Pa
        composition (%): Water              80
                         Ethanol            2.49e-05
                         AceticAcid         0.0111
                         Furfural           0.00447
                         Glycerol           0.00123
                         NH3                0.0448
                         LacticAcid         0.0117
                         SuccinicAcid       0.00307
                         DAP                0.0636
                         AmmoniumSulfate    1.12
                         HMF                0.00224
                         Glucose            0.0038
                         Xylose             0.00679
                         Arabinose          0.0144
                         Extract            0.0872
                         Ash                0.832
                         Lignin             2.51
                         SolubleLignin      0.00469
                         GlucoseOligomer    0.0081
                         GalactoseOligomer  2.05e-05
                         MannoseOligomer    1.07e-05
                         XyloseOligomer     0.0028
                         ArabinoseOligomer  0.000341
                         Z_mobilis          0.163
                         Protein            0.584
                         Glucan             0.251
                         Xylan              0.0805
                         Xylitol            0.00556
                         Cellobiose         0.00242
                         Arabinan           0.0295
                         Mannan             0.104
                         Galactan           0.0243
                         WWTsludge          14
                         Cellulase          0.00447
                         -----------------  1e+04 kg/hr
    [3] RO_treated_water  
        phase: 'l', T: 303.15 K, P: 101325 Pa
        composition (%): Water  100
                         -----  4.59e+05 kg/hr
    [4] brine  
        phase: 'l', T: 303.15 K, P: 101325 Pa
        composition (%): Water              71.2
                         Ethanol            2.18e-05
                         AceticAcid         0.00971
                         Furfural           0.00391
                         Glycerol           0.00108
                         NH3                1.07
                         LacticAcid         0.0121
                         SuccinicAcid       0.00269
                         DAP                1.48
                         AmmoniumSulfate    26
                         HMF                0.00195
                         Glucose            0.00332
                         Xylose             0.0133
                         Arabinose          0.0126
                         Extract            0.09
                         SolubleLignin      0.00484
                         GlucoseOligomer    0.00822
                         GalactoseOligomer  2.08e-05
                         MannoseOligomer    1.09e-05
                         XyloseOligomer     0.00284
                         ArabinoseOligomer  0.000346
                         Xylitol            0.00486
                         Cellobiose         0.00211
                         Cellulase          0.00462
                         -----------------  8.51e+03 kg/hr
    >>> u = wwt_sys.flowsheet.unit
    >>> print(round(u.R603.outs[3].F_mol, 2))
    50.46
    >>> wwt_sys.simulate()
    >>> print(round(u.R603.outs[3].F_mol, 2))
    50.46
    >>> wwt_sys.simulate()
    >>> print(round(u.R603.outs[3].F_mol, 2))
    50.46
    
    Check if system can finish simulating with a dilute influent stream.
    
    >>> from biosteam import Stream, create_high_rate_wastewater_treatment_system, settings
    >>> from biorefineries import cornstover as cs
    >>> settings.set_thermo(cs.create_chemicals())
    >>> feed = Stream(
    ...     ID='wastewater', 
    ...     Water=2.634e+05, 
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
    >>> wwt_sys = create_high_rate_wastewater_treatment_system(ins=feed)
    >>> wwt_sys.simulate()
    >>> wwt_sys.show('cwt100')
    System: wastewater_sys
    Highest convergence error among components in recycle
    stream M603-0 after 6 loops:
    - flow rate   5.52e+00 kmol/hr (0.027%)
    - temperature 4.36e-04 K (0.00014%)
    ins...
    [0] wastewater  
        phase: 'l', T: 298.15 K, P: 101325 Pa
        composition (%): Water              99.4
                         Ethanol            6.98e-05
                         AceticAcid         0.031
                         Furfural           0.0125
                         Glycerol           0.00344
                         LacticAcid         0.0334
                         SuccinicAcid       0.00859
                         DAP                0.00277
                         AmmoniumSulfate    0.0488
                         HMF                0.00625
                         Glucose            0.0106
                         Xylose             0.0219
                         Arabinose          0.0402
                         Extract            0.249
                         Ash                0.00175
                         Lignin             0.00529
                         SolubleLignin      0.0134
                         GlucoseOligomer    0.0231
                         GalactoseOligomer  5.84e-05
                         MannoseOligomer    3.06e-05
                         XyloseOligomer     0.00797
                         ArabinoseOligomer  0.000971
                         Z_mobilis          0.000344
                         Protein            0.00123
                         Glucan             0.000528
                         Xylan              0.000169
                         Xylitol            0.0156
                         Cellobiose         0.00676
                         Arabinan           6.21e-05
                         Mannan             0.000219
                         Galactan           5.11e-05
                         Cellulase          0.0128
                         -----------------  4.77e+06 kg/hr
    outs...
    [0] RNG  
        phase: 'g', T: 298.15 K, P: 101325 Pa
        flow: 0
    [1] biogas  
        phase: 'g', T: 298.15 K, P: 101325 Pa
        composition (%): CH4  26.9
                         H2S  0.0348
                         CO2  73.1
                         ---  4.33e+04 kg/hr
    [2] sludge  
        phase: 'l', T: 307.9 K, P: 101325 Pa
        composition (%): Water              79.9
                         Ethanol            2e-05
                         AceticAcid         0.00889
                         Furfural           0.00358
                         Glycerol           0.000986
                         NH3                0.0359
                         LacticAcid         0.00936
                         SuccinicAcid       0.00246
                         DAP                0.0511
                         AmmoniumSulfate    0.9
                         HMF                0.00179
                         Glucose            0.00304
                         Xylose             0.00544
                         Arabinose          0.0115
                         Extract            0.0698
                         Ash                0.838
                         Lignin             2.53
                         SolubleLignin      0.00375
                         GlucoseOligomer    0.00649
                         GalactoseOligomer  1.64e-05
                         MannoseOligomer    8.6e-06
                         XyloseOligomer     0.00224
                         ArabinoseOligomer  0.000273
                         Z_mobilis          0.165
                         Protein            0.589
                         Glucan             0.253
                         Xylan              0.0811
                         Xylitol            0.00445
                         Cellobiose         0.00193
                         Arabinan           0.0297
                         Mannan             0.105
                         Galactan           0.0245
                         WWTsludge          14.3
                         Cellulase          0.00358
                         -----------------  9.97e+03 kg/hr
    [3] RO_treated_water  
        phase: 'l', T: 303.15 K, P: 101325 Pa
        composition (%): Water  100
                         -----  4.67e+06 kg/hr
    [4] brine  
        phase: 'l', T: 303.15 K, P: 101325 Pa
        composition (%): Water              96.1
                         Ethanol            2.92e-06
                         AceticAcid         0.0013
                         Furfural           0.000522
                         Glycerol           0.000144
                         NH3                0.144
                         LacticAcid         0.00161
                         SuccinicAcid       0.000359
                         DAP                0.198
                         AmmoniumSulfate    3.49
                         HMF                0.000261
                         Glucose            0.000445
                         Xylose             0.00178
                         Arabinose          0.00168
                         Extract            0.012
                         SolubleLignin      0.000647
                         GlucoseOligomer    0.0011
                         GalactoseOligomer  2.78e-06
                         MannoseOligomer    1.46e-06
                         XyloseOligomer     0.000379
                         ArabinoseOligomer  4.62e-05
                         Xylitol            0.000651
                         Cellobiose         0.000282
                         Cellulase          0.000617
                         -----------------  6.42e+04 kg/hr
    
    """
    # Setup
    if flowsheet: bst.main_flowsheet.set_flowsheet(flowsheet)
    RNG, biogas, sludge, RO_treated_water, brine = outs
    RO_treated_water.register_alias('recycled_water')

    ##### Units #####
    # Mix waste liquids for treatment
    X = str(process_ID)
    MX01 = bst.Mixer(f'M{X}01', ins=ins or [])
    MX01.autopopulate = False if autopopulate is None else autopopulate
    
    @MX01.add_specification(run=True)
    def autopopulate_waste_streams():
        if MX01.autopopulate and not MX01.ins: 
            sys = MX01.system
            streams = bst.FreeProductStreams(sys.streams)
            MX01.ins.extend(streams.noncombustible_slurries)
            for i in sys.facilities:
                if isinstance(i, bst.BlowdownMixer): MX01.ins.append(i.outs[0])

    RX01_outs = (f'biogas_R{X}01', 'IC_eff', 'IC_sludge')
    if skip_IC:
        RX01 = Skipped(f'R{X}01', ins=MX01-0, outs=RX01_outs, main_in=0, main_out=1) # outs[0] is the gas phase
    else:
        RX01 = InternalCirculationRx(f'R{X}01', ins=MX01-0, outs=RX01_outs,
                                     T=35+273.15, **IC_kwargs)

    RX02_outs = (f'biogas_R{X}02', f'permeate_R{X}02', f'sludge_R{X}02', f'vent_R{X}02')
    if skip_AnMBR:
        RX02 = Skipped(f'R{X}02', ins=RX01-1, outs=RX02_outs, main_in=0, main_out=1)
    else:
        # Just setting the prices, flows will be updated upon simulation
        naocl_RX02 = tmo.Stream(f'naocl_R{X}02', NaOCl=0.125, Water=1-0.125, units='kg/hr')
        naocl_RX02.price = (naocl_RX02.F_mass/naocl_RX02.F_vol/1000)*prices['naocl'] # $/L to $/kg
        citric_RX02 = tmo.Stream(f'citric_R{X}02', CitricAcid=1, units='kg/hr')
        citric_RX02.price = (citric_RX02.F_mass/citric_RX02.F_vol/1000)*prices['citric_acid'] # $/L to $/kg
        bisulfite_RX02 = tmo.Stream(f'bisulfite_R{X}02', Bisulfite=0.38, Water=1-0.38, units='kg/hr')
        bisulfite_RX02.price = (bisulfite_RX02.F_mass/bisulfite_RX02.F_vol/1000)*prices['bisulfite'] # $/L to $/kg

        RX02 = AnMBR(f'R{X}02', ins=(RX01-1, '', naocl_RX02, citric_RX02,
                                     bisulfite_RX02, f'air_R{X}02'),
                     outs=RX02_outs,
                     reactor_type='CSTR',
                     membrane_configuration='cross-flow',
                     membrane_type='multi-tube',
                     membrane_material='ceramic',
                     include_aerobic_filter=False,
                     add_GAC=False,
                     include_degassing_membrane=True,
                     # Below include in the TEA
                     include_pump_building_cost=False,
                     include_excavation_cost=False, **AnMBR_kwargs)

    RX03_outs = (f'biogas_R{X}03', f'treated_R{X}03', f'sludge_R{X}03', f'vent_R{X}03')
    if skip_AeF:
        RX03 = Skipped(f'R{X}03', ins=(RX02-1, ''), outs=RX03_outs, main_in=0, main_out=1)
    else:
        RX03 = PolishingFilter(f'R{X}03', ins=(RX02-1, '', f'air_R{X}03'), outs=RX03_outs,
                              filter_type='aerobic',
                              include_degassing_membrane=False,
                              # Below include in the TEA
                              include_pump_building_cost=False,
                              include_excavation_cost=False,
                              **AeF_kwargs)

    MX02 = bst.Mixer(f'M{X}02', ins=(RX01-0, RX02-0, RX03-0))
    
    BiogasUpgrading('Upgrading', ins=(MX02-0, 'upgrading_placeholder'), outs=(RNG, biogas))

    # Recycled the majority of sludge (96%) to the aerobic filter,
    # 96% from the membrane bioreactor in ref [2]
    SX01 = bst.Splitter(f'S{X}01', ins=RX03-2, outs=(f'recycled_S{X}01', f'wasted_S{X}01'),
                              split=0.96)

    solubles = [i.ID for i in SX01.chemicals if not i.ID in default_insolubles]
    SX02 = BeltThickener(f'S{X}02', ins=(RX01-2, RX02-2, SX01-1),
                         outs=(f'eff_S{X}02', f'sludge_S{X}02'),
                         sludge_moisture=0.96, solubles=solubles)

    SX03 = SludgeCentrifuge(f'S{X}03', ins=SX02-1,
                            outs=(f'centrate_S{X}03', sludge),
                            sludge_moisture=0.8, solubles=solubles,
                            centrifuge_type='reciprocating_pusher')

    # Mix recycles to aerobic digestion
    bst.Mixer(f'M{X}03', ins=(SX01-0, SX02-0, SX03-0), outs=1-RX03)

    # Reverse osmosis to treat aerobically polished water
    ReverseOsmosis(f'S{X}04', ins=RX03-1, outs=(RO_treated_water, brine))