# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 00:55:03 2019

Equipment from the Humbird 2011 Report.
    
Humbird, D., Davis, R., Tao, L., Kinchin, C., Hsu, D., Aden, A., Dudgeon, D. (2011). Process Design and Economics for Biochemical Conversion of Lignocellulosic Biomass to Ethanol: Dilute-Acid Pretreatment and Enzymatic Hydrolysis of Corn Stover (No. NREL/TP-5100-47764, 1013269). https://doi.org/10.2172/1013269

@author: yoelr
"""

from .. import Unit, Stream, units
from .metaclasses import splitter
from . import decorators
from .designtools import CEPCI_by_year as CE

hp2kW = 0.7457

# %% Decorators 

_massflow_units = {'Flow rate': 'kg/hr'}
def _design_kg_hr(self):
    self._results['Design']['Flow rate'] = self._ins[0].massnet

def cost_kg_hr(name=None, *, cost, exp, S, kW=0, CE=CE[2009], N=1):
    def decorator(cls):
        cls._design = _design_kg_hr
        cls._units = _massflow_units
        return decorators.cost('Flow rate', cost=cost, exp=exp,
                               CE=CE, S=S, kW=kW, N=N)(cls)
    return decorator

def _design_Gcal_hr(self):
    duty = self._outs[0].H - self._ins[0].H
    self._heat_utilities[0](duty, self.ins[0].T, self._kwargs['T'])
    self._results['Design']['Duty'] = duty*2.39e-7  # Gcal/hr

_heatflow_units = {'Duty': 'Gcal/hr'}
def heat_utility(name=None, *, cost, S, exp=0.7, kW=0, N=1, CE=CE[2009], BM=2.2):
    """Decorate class as a heat exchanger."""
    def decorator(cls):
        cls.BM = BM
        cls._graphics = units.HXutility._graphics
        cls._linkedstreams = True
        cls._N_heat_utilities = 1
        cls._units = _heatflow_units
        cls._design = _design_Gcal_hr
        return decorators.cost('Duty', name, cost=cost, exp=exp,
                               CE=CE, S=S, kW=kW, N=N)(cls)
    return decorator
    

# %% Units

@cost_kg_hr(cost=13329690, exp=0.6, S=94697, kW=511.321)
class FeedStockHandling(Unit):
    """All area 100 equipment:
        * C101 Transfer Conveyor (2)
        * C102 High Angle Transfer Converyor (2)  
        * C103 Reversing Load-in Conveyor
        * C104 Dome Reclaim System (2)
        * C106 High Angle Transfer Conveyor
        * C107 Elevated Transfer Conveyor
        * M101 Truck Scale (2)
        * M102 Truck Dumper (2)
        * M103 Truck Dumper Hopper (2)
        * M104 Concrete Feedstock Storage Dome (2)
        * M105 Belt Scale (2)
        * M106 Dust Collection System (6)
    """
    _linkedstreams = True
    BM = 1.7
    
@cost_kg_hr(cost=6000, exp=0.5, S=136260)
class SulfuricAcidMixer(Unit):
    """A-201"""
    _linkedstreams = True
    BM = 1.0
    
@cost_kg_hr(cost=19812400, exp=0.6, S=83333, kW=5290)
class PretreatmentReactorSystem(Unit):
    """Includes the following:
        * C201 Transfer Conveyor (2)
        * C202 Distribution Conveyor (2)
        * C203 Overfeed Conveyor (4)
        * C204 Pressurized Heating Screw
        * C205 Pressurized Pre-heater Discharge (2)
        * C206 Pressurized Transport #1
        * C207 Pressurized Transport #2
        * M201 Doffing Roll Storage Bings (2)
        * M202 Pin Drum Feeder (2)
        * M203 Plug Screw Feeder (2)
        * M204 Prehydrolysis / Vertical Preheater
        * M205 Pin Drum Feeder (2)
        * M206 Plug Screw Feeder  (2)
        * M207 Pretreatment Reactor (3)
    """
    _linkedstreams = True
    BM = 1.5
    def _init(self):
        self._water_mass = Stream.indices('Water')
    
    def _design(self):
        feed = self._ins[0]
        self._results['Design']['Flow rate'] = feed.massnet - feed.mass[self._water_index]
        

@heat_utility(cost=92e3, CE=CE[2010], S=-8, exp=0.70)
class PretreatmentWaterHeater(Unit):
    """H201"""
    _kwargs = {'T': None}
    
    def _run(self):
        out = self._outs[0]
        out.P = self._ins[0].P
        out.T = self._kwargs['T']
    
    
@heat_utility(cost=34e3, S=2, exp=0.70)
class WasteVaporCondenser(Unit):
    """H244"""
    _kwargs = {'T': None}
    def _setup(self):
        self._outs[0].phase = 'l'
    
    def _run(self):
        feed = self._ins[0]
        out = self._outs[0]
        out.P = feed.P
        out.T = self._kwargs['T'] or feed.T

@decorators.cost('Flow rate', 'Discharge pump', kW=75*hp2kW, CE=CE[2009], cost=30e3, exp=0.80)
@cost_kg_hr('Agitators', N=3, cost=90e3/3, S=252891, exp=0.5)
class Flash204(units.Flash):
    """Includes:
        * Discharge pump
        * Agitator
    """
    BM = {'Discharge pump': 1.5,
          'Agitators': 2.3}
    
    
    
        