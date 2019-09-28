# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 00:55:03 2019

Equipment from the Humbird 2011 Report.
    
Humbird, D., Davis, R., Tao, L., Kinchin, C., Hsu, D., Aden, A., Dudgeon, D. (2011). Process Design and Economics for Biochemical Conversion of Lignocellulosic Biomass to Ethanol: Dilute-Acid Pretreatment and Enzymatic Hydrolysis of Corn Stover (No. NREL/TP-5100-47764, 1013269). https://doi.org/10.2172/1013269

@author: yoelr
"""
import os
import sys
from biosteam.units.factories import xl2mod
path = os.path.dirname(os.path.realpath(__file__)) + '\\'
xl2mod(path + '_humbird2011.xlsx', sys.modules[__name__])
del sys, xl2mod, os, path

from biosteam.utils.solvers import aitken_secant
from biosteam import Unit, MixedStream, Stream
from biosteam.units.decorators import cost
from biosteam.units.designtools import size_batch
import biosteam as bst
Rxn = bst.reaction.Reaction
ParallelRxn = bst.reaction.ParallelReaction
_gal2m3 = 0.003785
_gpm2m3hr = 0.227124
# _m3hr2gpm = 4.40287
_hp2kW = 0.7457
_Gcal2kJ = 4184e3

# %% Pretreatment

class SteamMixer(Unit):
    """
    **ins**
    
        [0] Feed
        
        [1] Steam
    
    **outs**
    
        [0] Mixed
    
    """
    _N_outs = 1
    _N_ins = 2
    _N_heat_utilities = 1
    def __init__(self, ID='', ins=None, outs=(), *, P):
        super().__init__(ID, ins, outs)
        self.P = P
    
    @staticmethod
    def _P_at_flow(mol_water, P, mol_array, index_water, mixed, ins):
        mol_array[index_water] = mol_water
        Stream.sum(mixed, ins)
        P_new = mixed.bubble_P()[0]
        return P - P_new
    
    def _run(self):
        ins = self._ins
        steam = ins[1]
        steam_mol = steam.molnet
        mixed = self.outs[0]
        index_water = mixed.index('7732-18-5')
        steam_mol = aitken_secant(self._P_at_flow,
                                  steam_mol, steam_mol+0.1, 
                                  1e-4, 1e-4,
                                  args=(self.P, steam.mol, index_water, mixed, ins))
        mixed.P = self.P
        hu = self._heat_utilities[0]
        hu.ID = 'Low pressure steam'
        hu.flow = steam_mol
        hu.cost = steam_mol*bst.HeatUtility.heating_agents['Low pressure steam'].price_kmol
    
    @property
    def installation_cost(self): return 0
    @property
    def purchase_cost(self): return 0
    def _design(self): pass
    def _cost(self): pass
        

@cost('Dry flow rate', units='kg/hr', S=83333, CE=522,
      cost=19812400, n=0.6, kW=4578, BM=1.5)
class PretreatmentReactorSystem(Unit):
    _N_ins = 1
    _N_outs = 2
    _graphics = bst.Flash._graphics
    def __init__(self, ID='', ins=None, outs=()):
        Unit.__init__(self, ID, ins, outs)
        self._mixedstream = MixedStream(None)
        self.reactions = ParallelRxn([
    #            Reaction definition                 Reactant    Conversion
    Rxn('Glucan + H2O -> Glucose',                   'Glucan',   0.0990),
    Rxn('Glucan + H2O -> GlucoseOligomer',           'Glucan',   0.0030),
    Rxn('Glucan -> HMF + 2 H2O',                     'Glucan',   0.0030),
    Rxn('Galactan + H2O -> GalactoseOligomer',       'Galactan', 0.0030),
    Rxn('Galactan -> HMF + 2 H2O',                   'Galactan', 0.0030),
    Rxn('Mannan + H2O -> MannoseOligomer',           'Mannan',   0.0030),
    Rxn('Mannan -> HMF + 2 H2O',                     'Mannan',   0.0030),
    Rxn('Sucrose -> HMF + Glucose + 2H2O',           'Sucrose',  0.0030),
    Rxn('Xylan + H2O -> Xylose',                     'Xylan',    0.9000),
    Rxn('Xylan + H2O -> XyloseOligomer',             'Xylan',    0.0024),
    Rxn('Xylan -> Furfural + 2 H2O',                 'Xylan',    0.0050),
    Rxn('Arabinan + H2O -> Arabinose',               'Arabinan', 0.9000),
    Rxn('Arabinan + H2O -> ArabinoseOligomer',       'Arabinan', 0.0024),
    Rxn('Arabinan -> Furfural + 2 H2O',              'Arabinan', 0.0050),
    Rxn('Acetate -> AceticAcid',                     'Acetate',  1.0000),
    Rxn('Lignin -> SolubleLignin',                   'Lignin',   0.0050)])
        vapor, liquid = self.outs
        vapor.phase = 'g'
    
    def _run(self):
        ms = self._mixedstream
        feed = self.ins[0]
        vapor, liquid = self.outs
        liquid.copylike(feed)
        self.reactions(liquid.mol) 
        ms.copylike(liquid)
        ms.VLE(T=130+273.15, Q=(liquid.Hf-feed.Hf))
        vapor.mol[:] = ms.vapor_mol
        liquid.mol[:] = ms.liquid_mol
        vapor.T = liquid.T = ms.T
        vapor.P = liquid.P = ms.P


@cost('Flow rate', 'Pumps',
      S=43149, CE=522, cost=24800, n=0.8, kW=40, BM=2.3)
@cost('Stage #1 reactor volume', 'Stage #1 reactors',
      cost=37700, S=20*_gal2m3, CE=522, n=0.7, BM=1.8)
@cost('Stage #2 reactor volume', 'Stage #2 reactors',
      cost=58300, S=200*_gal2m3, CE=522, n=0.7, BM=1.8)
@cost('Stage #3 reactor volume', 'Stage #3 reactors',
      cost=78800, S=2e3*_gal2m3, CE=522, n=0.7, BM=1.8)
@cost('Stage #4 reactor volume', 'Stage #4 reactors',
      cost=176e3, S=20e3*_gal2m3, CE=522, n=0.7, BM=1.8)
@cost('Stage #4 reactor volume', 'Stage #4 agitators',
      cost=26e3/2, S=20e3*_gal2m3, kW=7.5, CE=522, n=0.5, BM=1.5)
@cost('Stage #5 reactor volume', 'Stage #5 reactors',
      cost=590e3, S=200e3*_gal2m3, CE=522, n=0.7, BM=1.8)
@cost('Stage #5 reactor volume', 'Stage #5 agitators',
      cost=43e3/2, S=200e3*_gal2m3, kW=10, CE=522, n=0.5, BM=1.5)
class SeedTrain(Unit):
    _N_heat_utilities = 1
    
    _units= {'Flow rate': 'kg/hr',
             'Stage #1 reactor volume': 'm3',
             'Stage #2 reactor volume': 'm3',
             'Stage #3 reactor volume': 'm3',
             'Stage #4 reactor volume': 'm3',
             'Stage #5 reactor volume': 'm3'}
    
    @property
    def N_stages(self): 
        """Number of stages in series."""
        return 5
    
    #: Number of parallel seed trains
    N_trains = 2
    
    #: Cycle time for each batch (hr)
    tau_batch = 24
    
    @property
    def tau_turnover(self):
        """Turnover time (hr) calculated by batch time divided by number of trains."""
        return self.tau_batch/self.N_trains
    
    #: Operating temperature (K)
    T = 32+273.15
    
    # #: wt % media (e.g. corn steep liquor) in each stage 
    # media_loading = 0.50
    
    # #: Diammonium phosphate loading in g/L of fermentation broth
    # DAP = 0.67 
    
    def __init__(self, ID='', ins=None, outs=()):
        Unit.__init__(self, ID, ins, outs)
        self.reactions = ParallelRxn([
    #   Reaction definition                             Reactant    Conversion
    Rxn('Glucose -> 2 Ethanol + 2 CO2',                 'Glucose',   0.9500),
    Rxn('Glucose + 8.463 CSL + 0.018 DAP -> 6 Z_mobilis + 2.4 H2O',
                                                        'Glucose',   0.0200),
    Rxn('Glucose + 2 H2O -> 2 Glycerol + O2',           'Glucose',   0.0040),
    Rxn('Glucose + 2 CO2 -> 2 SuccinicAcid + O2',       'Glucose',   0.0060),
    Rxn('3 Xylose -> 5 Ethanol + 5 CO2',                'Xylose',    0.8000),
    Rxn('Xylose + 7.052 CSL + 0.015 DAP -> 5 Z_mobilis + 2 H2O',
                                                        'Xylose',    0.0400),
    Rxn('3 Xylose + 5 H2O -> 5 Glycerol + 2.5 O2',      'Xylose',    0.0030),
    Rxn('Xylose + H2O -> Xylitol + 0.5 O2',             'Xylose',    0.0460),
    Rxn('3 Xylose + 5 CO2 -> 5 SuccinicAcid + 2.5 O2',  'Xylose',    0.0090)])
    
    def _run(self):
        feed = self.ins[0]
        vent, effluent= self.outs
        effluent.copyflow(feed)
        self.reactions(effluent.mol)
        effluent.T = self.T
        vent.phase = 'g'
        vent.copyflow(effluent, ('CO2', 'NH3', 'O2'), remove=True)

    def _design(self): 
        maxvol = self.outs[1].volnet*self.tau_turnover
        vol = maxvol*10**-self.N_stages
        Design = self._Design
        for i in range(1, self.N_stages+1):
            Design[f'Stage #{i} reactor volume'] = vol
            vol *= 10 
        Design['Flow rate'] = sum([i.massnet for i in self.outs])
        self._heat_utilities[0](self._Hnet, self.T)

    def _cost(self):
        N = self.N_trains
        D = self._Design
        C = self._Cost
        kW = 0
        for i, x in self.cost_items.items():
            S = D[x._basis]
            q = S/x.S
            C[i] = N*bst.CE/x.CE*x.cost*q**x.n
            kW += N*x.kW*q
        self._power_utility(kW)
        

# %% Saccharification and fermentation

@cost('Flow rate', 'Recirculation pumps', kW=30, S=340*_gpm2m3hr,
      cost=47200, n=0.8, BM=2.3, CE=522, N='N_recirculation_pumps')
@cost('Reactor duty', 'Heat exchangers', CE=522, cost=23900,
      S=5*_Gcal2kJ, n=0.7, BM=2.2, N='N_reactors') # Based on a similar heat exchanger
@cost('Reactor volume', 'Agitators', CE=522, cost=52500,
      S=1e6*_gal2m3, n=0.5, kW=90, BM=1.5, N='N_reactors')
@cost('Reactor volume', 'Reactors', CE=522, cost=844000,
      S=1e6*_gal2m3, n=0.5, BM=1.5, N='N_reactors')
@cost('Flow rate', 'Transfer pumps', kW=58, S=352*_gpm2m3hr,
      cost=47200/5, CE=522, n=0.8, BM=2.3, N='N_transfer_pumps')
@cost('Tank volume', 'Tanks', cost=3840e3/8, S=250e3*_gal2m3, 
      CE=522, n=0.7, BM=2.0, N='N_tanks')
class SaccharificationAndCoFermentation(Unit):
    _N_ins = 3
    _N_outs = 3
    _N_heat_utilities = 2
    
    #: Saccharification temperature (K)
    T_saccharification = 48+273.15
    
    #: Fermentation temperature (K)
    T_fermentation = 32+273.15
    
    #: Residence time of countinuous saccharification tanks (hr)
    tau_tank = 24
    
    #: Saccharification time (hr)
    tau_saccharification = 60
    
    #: Co-Fermentation time (hr)
    tau_cofermentation = 36
    
    #: Unload and clean up time (hr)
    tau_0 = 4
    
    #: Working volume fraction (filled tank to total tank volume)
    V_wf = 0.9
    
    #: Number of reactors
    N_reactors = 12
    
    #: Number of continuous saccharification tanks
    N_tanks = 8
    
    #: Number of transfer pumps
    N_transfer_pumps = 5
    
    #: Number of recirculation pumps
    N_recirculation_pumps = 5
    
    _units = {'Flow rate': 'm3/hr',
              'Tank volume': 'm3',
              'Reactor volume': 'm3',
              'Reactor duty': 'kJ/hr'}
    
    # Split to outs[2]
    saccharified_slurry_split = 0.1
    
    def __init__(self, ID='', ins=None, outs=(), P=101325):
        Unit.__init__(self, ID, ins, outs)
        self.P = P
        #: [ParallelReaction] Enzymatic hydrolysis reactions including from downstream batch tank in co-fermentation.
        self.saccharification = ParallelRxn([
    #   Reaction definition                   Reactant    Conversion
    Rxn('Glucan -> GlucoseOligomer',          'Glucan',   0.0400),
    Rxn('Glucan + 0.5 H2O -> 0.5 Cellobiose', 'Glucan',   0.0120),
    Rxn('Glucan + H2O -> Glucose',            'Glucan',   0.9000),
    Rxn('Cellobiose + H2O -> Glucose',        'Cellobiose',  1.0000)])
    
        self.substrate_loss = ParallelRxn([
    # Losses
    Rxn('Glucose -> 2 LacticAcid',                      'Glucose',   0.0300),
    Rxn('3 Xylose -> 5 LacticAcid',                     'Xylose',    0.0300),
    Rxn('3 Arabinose -> 5 LacticAcid',                  'Arabinose', 0.0300),
    Rxn('Galactose -> 2 LacticAcid',                    'Galactose', 0.0300),
    Rxn('Mannose -> 2 LacticAcid',                      'Mannose',   0.0300),])
    
        self.cofermentation = ParallelRxn([
    #   Reaction definition                             Reactant    Conversion
    Rxn('Glucose -> 2 Ethanol + 2 CO2',                 'Glucose',   0.9500),
    Rxn('Glucose + 0.047 CSL + 0.018 DAP -> 6 Z_mobilis + 2.4 H2O',
                                                        'Glucose',   0.0200),
    Rxn('Glucose + 2 H2O -> 2 Glycerol + O2',           'Glucose',   0.0040),
    Rxn('Glucose + 2 CO2 -> 2 SuccinicAcid + O2',       'Glucose',   0.0060),
    Rxn('3 Xylose -> 5 Ethanol + 5 CO2',                'Xylose',    0.8500),
    Rxn('Xylose + 0.039 CSL + 0.015 DAP -> 5 Z_mobilis + 2 H2O',
                                                        'Xylose',    0.0190),
    Rxn('3 Xylose + 5 H2O -> 5 Glycerol + 2.5 O2',      'Xylose',    0.0030),
    Rxn('Xylose + H2O -> Xylitol + 0.5 O2',             'Xylose',    0.0460),
    Rxn('3 Xylose + 5 CO2 -> 5 SuccinicAcid + 2.5 O2',  'Xylose',    0.0090)])
    
        self.CSL2constituents = Rxn(
        'CSL -> 0.5 H2O + 0.25 LacticAcid + 0.25 Protein', 'CSL',    1.0000)
    
        self.saccharified_stream = bst.Stream(None)
    
    def _run(self):
        feed, CSL, DAP = self.ins
        vent, effluent, sidedraw = self.outs
        vent.P = effluent.P = sidedraw.P = self.P
        ss = self.saccharified_stream
        ss.T = sidedraw.T = self.T_saccharification
        vent.T = effluent.T = self.T_fermentation
        vent.phase = 'g'
        ss.copyflow(feed)
        self.saccharification(ss.mol)
        sidedraw.mol[:] = ss.mol*self.saccharified_slurry_split
        effluent.mol[:] = ss.mol - sidedraw.mol + CSL.mol + DAP.mol
        self.substrate_loss(effluent.mol)
        self.cofermentation(effluent.mol)
        self.CSL2constituents(effluent.mass)
        vent.copyflow(effluent, ('CO2', 'NH3', 'O2'), remove=True)
        vent.recieve_vent(effluent)
    
    def _design(self):
        effluent = self.outs[1]
        v_0 = effluent.volnet
        Design = self._Design
        Design['Tank volume'] = v_0*self.tau_tank/self.V_wf/self.N_tanks
        Design['Flow rate'] = v_0/self.N_transfer_pumps
        tau = self.tau_saccharification + self.tau_cofermentation
        Design.update(size_batch(v_0, tau, self.tau_0, self.N_reactors, self.V_wf))
        hu_cooling, hu_fermentation = self._heat_utilities
        hu_cooling(self.saccharified_stream.H_at(T=self.T_fermentation)
                   - self.saccharified_stream.H_at(T=self.T_saccharification),
                   self.T_fermentation)
        ei = effluent.index('Ethanol')
        ethanol = (sum([i._mol[ei] for i in self.outs])
                   - sum([i._mol[ei] for i in self.ins]))
        duty = ethanol*-5568
        hu_fermentation(duty, effluent.T)
        Design['Reactor duty'] = -duty
        
   
# %% Lignin separation

@cost('Flow rate', 'Flitrate tank agitator',
      cost=26e3, CE=551, kW=7.5*_hp2kW, S=31815, n=0.5, BM=1.5)
@cost('Flow rate', 'Discharge pump',
      cost=13040, CE=551, S=31815, n=0.8, BM=2.3)
@cost('Flow rate', 'Filtrate tank',
      cost=103e3, S=31815, CE=551, BM=2.0, n=0.7)
@cost('Flow rate', 'Feed pump', kW=74.57,
      cost= 18173, S=31815, CE=551, n=0.8, BM=2.3)
@cost('Flow rate', 'Stillage tank 531',
      cost=174800, CE=551, S=31815, n=0.7, BM=2.0)
@cost('Flow rate', 'Mafifold flush pump', kW=74.57,
      cost=17057, CE=551, S=31815, n=0.8, BM=2.3)
@cost('Flow rate', 'Recycled water tank',
      cost=1520, CE=551, S=31815, n=0.7, BM=3.0)
@cost('Flow rate', 'Lignin wet cake screw',  kW=15*_hp2kW,
      cost=2e4, CE=521.9, S=28630, n=0.8, BM=1.7)
@cost('Flow rate', 'Lignin wet cake conveyor', kW=10*_hp2kW,
      cost=7e4, CE=521.9, S=28630, n=0.8, BM=1.7)
@cost('Flow rate', 'Pressure filter',
      cost=3294700, CE=551, S=31815, n=0.8, BM=1.7)
@cost('Flow rate', 'Pressing air compressor reciever tank',
      cost=8e3, CE=551, S=31815, n=0.7, BM=3.1)
@cost('Flow rate', 'Cloth wash pump', kW=150*_hp2kW,
      cost=29154, CE=551, S=31815, n=0.8, BM=2.3)
@cost('Flow rate', 'Dry air compressor reciever tank',
      cost=17e3, CE=551, S=31815, n=0.7, BM=3.1)
@cost('Flow rate', 'Pressing air pressure filter',
      cost=75200, CE=521.9, S=31815, n=0.6, kW=112, BM=1.6)
@cost('Flow rate', 'Dry air pressure filter (2)',
      cost=405000, CE=521.9, S=31815, n=0.6, kW=1044, BM=1.6)
class PressureFilter(bst.SolidsSeparator):
    _units = {'Flow rate': 'kg/hr'}
    
    def _design(self):
        self._Design['Flow rate'] = self.outs[0].massnet


@cost('Flow rate', 'Waste water system', units='kg/hr', CE=551,
      cost=50280080., n=0.6, BM=1, kW=7139/1.05, S=393100)
class WasteWaterSystemCost(bst.Static): pass


class AnaerobicDigestion(bst.Unit):
    """Anaerobic digestion system as modeled by Humbird 2011
    
    **Parameters**
    
        **reactions:** [ReactionSet] Anaerobic digestion reactions.
        
        **sludge_split:** [Array] Split between waste water and sludge
        
    **ins**
    
        [0] Waste water
        
        [1] Cool well water
        
    **outs**
    
        [0] Biogas
        
        [1] Waste water
        
        [2] Sludge
        
        [3] Hot well water
    
    """
    purchase_cost = installation_cost = 0
    _N_ins = 2
    _N_outs = 4
    def __init__(self, ID='', ins=None, outs=(), *, reactions, sludge_split):
        Unit.__init__(self, ID, ins, outs)
        self.reactions = reactions
        self.sludge_split = sludge_split
        self.mixed_stream = bst.MixedStream()
    
    def _run(self):
        feed, cool_water = self.ins
        biogas, waste, sludge, hot_water = self.outs
        biogas.phase = 'g'
        hot_water.link = cool_water
        biogas.T = waste.T = sludge.T = T = 35+273.15
        hot_water.T = feed.T - 5
        cool_water.mol[:] *= (feed.H - feed.H_at(T=T))/(hot_water.H - cool_water.H)
        sludge.copyflow(feed)
        self.reactions(sludge.mol)
        self.mixed_stream.copyflow(sludge)
        self.mixed_stream.VLE(P=101325, Q=0)
        biogas.mol[:] = self.mixed_stream.vapor_mol
        liquid_mol = self.mixed_stream.liquid_mol
        sludge.mol[:] = liquid_mol * self.sludge_split
        waste.mol[:] = liquid_mol - sludge.mol
        biogas.recieve_vent(waste)
        
    
class AerobicDigestion(bst.Unit):
    """Anaerobic digestion system as modeled by Humbird 2011
    
    **Parameters**
    
        **reactions:** [ReactionSet] Anaerobic digestion reactions.
        
        **sludge_split:** [Array] Split between waste water and sludge
        
    **ins**
    
        [0] Waste water
        
        [1] Air
        
        [2] Caustic
        
    **outs**
    
        [0] Vent
        
        [1] Treated waste water

    """    
    _N_ins = 3
    _N_outs = 2
    purchase_cost = installation_cost = 0
    evaporation = 4/355
    
    def __init__(self, ID='', ins=None, outs=(), *, reactions):
        Unit.__init__(self, ID, ins, outs)
        self.reactions = reactions
    
    def _run(self):
        waste, air, caustic = self._ins
        vent, water = self.outs
        vent.phase = 'g'
        water.copylike(waste)
        water.mol[:] += air.mol
        water.mol[:] += caustic.mol
        self.reactions(water.mol)
        vent.copyflow(water, ('CO2', 'O2', 'N2'))
        wi = vent.index('Water')
        water_mol = water.mol[wi]
        vent.mol[wi] = water_mol * self.evaporation
        water.mol[:] -= vent.mol
        
@cost('Flow rate', units='kg/hr',
      S=63, cost=421e3, CE=522, BM=1.8, n=0.6)
class CIPpackage(bst.Facility):
    line = 'CIP Package'
    _N_ins = 1
    _N_outs = 1
    
        
    
# # %% Decorators 

# _massflow_units = {'Flow rate': 'kg/hr'}
# def _design_kg_hr(self):
#     self._results['Design']['Flow rate'] = self._ins[0].massnet

# def cost_kg_hr(name=None, *, cost, exp, S, kW=0, CE=CE[2009], N=1):
#     def decorator(cls):
#         cls._design = _design_kg_hr
#         cls._units = _massflow_units
#         return decorators.cost('Flow rate', cost=cost, exp=exp,
#                                CE=CE, S=S, kW=kW, N=N)(cls)
#     return decorator

# def _design_Gcal_hr(self):
#     duty = self._outs[0].H - self._ins[0].H
#     self._heat_utilities[0](duty, self.ins[0].T, self._kwargs['T'])
#     self._results['Design']['Duty'] = duty*2.39e-7  # Gcal/hr

# _heatflow_units = {'Duty': 'Gcal/hr'}
# def heat_utility(name=None, *, cost, S, exp=0.7, kW=0, N=1, CE=CE[2009], BM=2.2):
#     """Decorate class as a heat exchanger."""
#     def decorator(cls):
#         cls.BM = BM
#         cls._graphics = units.HXutility._graphics
#         cls._linkedstreams = True
#         cls._N_heat_utilities = 1
#         cls._units = _heatflow_units
#         cls._design = _design_Gcal_hr
#         return decorators.cost('Duty', name, cost=cost, exp=exp,
#                                CE=CE, S=S, kW=kW, N=N)(cls)
#     return decorator
    

# # %% Units

# @cost_kg_hr(cost=13329690, exp=0.6, S=94697, kW=511.321)
# class FeedStockHandling(Unit):
#     """All area 100 equipment:
#         * C101 Transfer Conveyor (2)
#         * C102 High Angle Transfer Converyor (2)  
#         * C103 Reversing Load-in Conveyor
#         * C104 Dome Reclaim System (2)
#         * C106 High Angle Transfer Conveyor
#         * C107 Elevated Transfer Conveyor
#         * M101 Truck Scale (2)
#         * M102 Truck Dumper (2)
#         * M103 Truck Dumper Hopper (2)
#         * M104 Concrete Feedstock Storage Dome (2)
#         * M105 Belt Scale (2)
#         * M106 Dust Collection System (6)
#     """
#     _linkedstreams = True
#     BM = 1.7
    
# @cost_kg_hr(cost=6000, exp=0.5, S=136260)
# class SulfuricAcidMixer(Unit):
#     """A-201"""
#     _linkedstreams = True
#     BM = 1.0
    
# @cost_kg_hr(cost=19812400, exp=0.6, S=83333, kW=5290)
# class PretreatmentReactorSystem(Unit):
#     """Includes the following:
#         * C201 Transfer Conveyor (2)
#         * C202 Distribution Conveyor (2)
#         * C203 Overfeed Conveyor (4)
#         * C204 Pressurized Heating Screw
#         * C205 Pressurized Pre-heater Discharge (2)
#         * C206 Pressurized Transport #1
#         * C207 Pressurized Transport #2
#         * M201 Doffing Roll Storage Bings (2)
#         * M202 Pin Drum Feeder (2)
#         * M203 Plug Screw Feeder (2)
#         * M204 Prehydrolysis / Vertical Preheater
#         * M205 Pin Drum Feeder (2)
#         * M206 Plug Screw Feeder  (2)
#         * M207 Pretreatment Reactor (3)
#     """
#     _linkedstreams = True
#     BM = 1.5
#     def _init(self):
#         self._water_mass = Stream.indices('Water')
    
#     def _design(self):
#         feed = self._ins[0]
#         self._results['Design']['Flow rate'] = feed.massnet - feed.mass[self._water_index]
        

# @heat_utility(cost=92e3, CE=CE[2010], S=-8, exp=0.70)
# class PretreatmentWaterHeater(Unit):
#     """H201"""
#     _kwargs = {'T': None}
    
#     def _run(self):
#         out = self._outs[0]
#         out.P = self._ins[0].P
#         out.T = self._kwargs['T']
    
    
# @heat_utility(cost=34e3, S=2, exp=0.70)
# class WasteVaporCondenser(Unit):
#     """H244"""
#     _kwargs = {'T': None}
#     def _setup(self):
#         self._outs[0].phase = 'l'
    
#     def _run(self):
#         feed = self._ins[0]
#         out = self._outs[0]
#         out.P = feed.P
#         out.T = self._kwargs['T'] or feed.T

# @decorators.cost('Flow rate', 'Discharge pump', kW=75*hp2kW, CE=CE[2009], cost=30e3, exp=0.80)
# @cost_kg_hr('Agitators', N=3, cost=90e3/3, S=252891, exp=0.5)
# class Flash204(units.Flash):
#     """Includes:
#         * Discharge pump
#         * Agitator
#     """
#     BM = {'Discharge pump': 1.5,
#           'Agitators': 2.3}