# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
.. contents:: :local:
    
Abstract Unit Operations
------------------------
.. autoclass:: biosteam.units.liquid_liquid_extraction.LLEUnit
.. autoclass:: biosteam.units.liquid_liquid_extraction.LiquidsCentrifuge
.. autoclass:: biosteam.units.liquid_liquid_extraction.LiquidsSettler

Centrifuges
-----------
.. autoclass:: biosteam.units.liquid_liquid_extraction.LiquidsSplitCentrifuge
.. autoclass:: biosteam.units.liquid_liquid_extraction.LLECentrifuge
.. autoclass:: biosteam.units.liquid_liquid_extraction.SLLECentrifuge
.. autoclass:: biosteam.units.liquid_liquid_extraction.SolidLiquidsSplitCentrifuge

Mixer-Settlers
--------------
.. autoclass:: biosteam.units.liquid_liquid_extraction.LiquidsMixingTank
.. autoclass:: biosteam.units.liquid_liquid_extraction.LLESettler
.. autoclass:: biosteam.units.liquid_liquid_extraction.LiquidsSplitSettler
.. autoclass:: biosteam.units.liquid_liquid_extraction.LiquidsPartitionSettler
.. autoclass:: biosteam.units.liquid_liquid_extraction.MixerSettler
.. autoclass:: biosteam.units.liquid_liquid_extraction.MultiStageMixerSettlers

References
----------
.. [1] Apostolakou, A. A.; Kookos, I. K.; Marazioti, C.; Angelopoulos, 
    K. C. Techno-Economic Analysis of a Biodiesel Production Process 
    from Vegetable Oils. Fuel Process. Technol. 2009, 90, 1023−1031
.. [2] Kwiatkowski, J. R.; McAloon, A. J.; Taylor, F.; Johnston, D. B. 
    Modeling the Process and Costs of Fuel Ethanol Production by the Corn 
    Dry-Grind Process. Industrial Crops and Products 2006, 23 (3), 288–296.
    https://doi.org/10.1016/j.indcrop.2005.08.004.

"""
import biosteam as bst
from .splitting import Splitter
from .design_tools import CEPCI_by_year, geometry, PressureVessel
from .decorators import cost, copy_algorithm
from .phase_equilibrium import MultiStageEquilibrium
from .._graphics import mixer_settler_graphics
from .. import Unit
from ._flash import RatioFlash
from thermosteam import separations as sep

__all__ = (
    'LLEUnit',
    'LiquidsCentrifuge',
    'LiquidsSplitCentrifuge', 
    'LiquidsRatioCentrifuge',
    'SLLECentrifuge', 
    'SolidLiquidsSplitCentrifuge',
    'LLECentrifuge',
    'LiquidsMixingTank',
    'LiquidsSettler', 
    'LLESettler', 
    'LiquidsSplitSettler',
    'LiquidsPartitionSettler',
    'MixerSettler',
    'MultiStageMixerSettlers',
)

# %% Abstract

class LLEUnit(bst.Unit, isabstract=True):
    r"""
    Abstract class for simulating liquid-liquid equilibrium.

    Parameters
    ----------
    ins : 
        Inlet fluid.
    outs : 
        * [0] Low density fluid.
        * [1] Heavy fluid.
    top_chemical : str, optional
        Identifier of chemical that will be favored in the low density phase.
    efficiency=1. : float, optional
        Fraction of feed in liquid-liquid equilibrium.
        The rest of the feed is divided equally between phases.
    cache_tolerance=1e-6 : float, optional
        Reuse previous partition coefficients to calculate LLE when 
        the change in molar fraction of all chemicals is below this 
        tolerance.
    forced_split_IDs : tuple[str], optional
        IDs of component with a user defined split.
    forced_split : 1d array, optional
        Component-wise split to 0th stream.
    
    Examples
    --------
    >>> from biorefineries.lipidcane import chemicals
    >>> from biosteam import units, settings, Stream
    >>> settings.set_thermo(chemicals['Methanol', 'Glycerol', 'Biodiesel', 'TAG'])
    >>> feed = Stream('feed', T=333.15,
    ...               TAG=0.996, Biodiesel=26.9,
    ...               Methanol=32.9, Glycerol=8.97)
    >>> C1 = units.LLEUnit('C1', ins=feed, outs=('light', 'heavy'))
    >>> C1.simulate()
    >>> C1.outs[0].show()
    Stream: light from <LLEUnit: C1>
    phase: 'l', T: 333.15 K, P: 101325 Pa
    flow (kmol/hr): Methanol   10.2
                    Glycerol   0.0239
                    Biodiesel  26.9
                    TriOlein   0.996
    
    """
    _N_outs = 2
    
    def _init(self, top_chemical=None, efficiency=1.0, cache_tolerance=1e-6,
              forced_split_IDs=None, forced_split=None):
        #: [str] Identifier of chemical that will be favored in the low density phase.
        self.top_chemical = top_chemical
        #: [float] Fraction of feed in liquid-liquid equilibrium.
        #: The rest of the feed is divided equally between phases.
        self.efficiency = efficiency 
        #: [float] The change in molar fraction of individual chemicals must be 
        #: below this tolerance to reuse partition coefficients.
        self.cache_tolerance = cache_tolerance
        #: array[float] Forced splits to 0th stream for given IDs. 
        self.forced_split = forced_split
        #: tuple[str] IDs corresponding to forced splits. 
        self.forced_split_IDs = forced_split_IDs
        self.multi_stream = bst.MultiStream(phases='lL', thermo=self.thermo)
        
    @property
    def solvent(self):
        return self.top_chemical
    @solvent.setter
    def solvent(self, solvent):
        self.top_chemical = solvent
        
    def _run(self):
        sep.lle(*self.ins, *self.outs, self.top_chemical, self.efficiency, self.multi_stream)
        IDs = self.forced_split_IDs
        if IDs:
            feed, = self.ins
            liquid, LIQUID = self.outs
            mol = feed.imol[IDs]
            liquid.imol[IDs] = mol_liquid = mol * self.forced_split
            LIQUID.imol[IDs] = mol - mol_liquid
   
    
# %% Centrifuge

# Electricity kW/(m3/hr) from USDA biosdiesel Super Pro model
# Possibly 1.4  kW/(m3/hr)
# https://www.sciencedirect.com/topics/engineering/disc-stack-centrifuge
# Microalgal fatty acids—From harvesting until extraction H.M. Amaro, et. al.,
# in Microalgae-Based Biofuels and Bioproducts, 2017

@cost('Flow rate', units='m^3/hr', CE=525.4, cost=28100,
      n=0.574, kW=1.4, ub=100, BM=2.03, N='Number of centrifuges')
class LiquidsCentrifuge(Unit, isabstract=True):
    r"""
    Abstract class for liquid centrifuges.

    Parameters
    ----------
    ins : 
        Inlet fluid.
    outs : 
        * [0] Low density fluid.
        * [1] Heavy fluid.
    
    Notes
    -----
    The f.o.b purchase cost is given by [1]_:

    .. math::
        
       C_{f.o.b}^{2007} = 28100 Q^{0.574} \ (Q < 100 \frac{m^3}{h})
    
    """
    _N_outs = 2
    line = 'Liquids centrifuge'


# TODO: Remove this in favor of partition coefficients
class LiquidsRatioCentrifuge(LiquidsCentrifuge):
    line = 'Liquids centrifuge'
    _init = RatioFlash._init
    _run = RatioFlash._run


class LiquidsSplitCentrifuge(LiquidsCentrifuge):
    r"""
    Create a liquids centrifuge simulated by component splits.

    Parameters
    ----------
    ins : 
        Inlet fluid.
    outs : 
        * [0] Low density fluid.
        * [1] Heavy fluid.
    split : Should be one of the following
            * [float] The fraction of net feed in the 0th outlet stream.
            * [array_like] Componentwise split of feed to 0th outlet stream.
            * [dict] ID-split pairs of feed to 0th outlet stream.
    order=None : Iterable[str], defaults to biosteam.settings.chemicals.IDs
        Chemical order of split.
    
    Notes
    -----
    The f.o.b purchase cost is given by [1]_:

    .. math::
        
       C_{f.o.b}^{2007} = 28100 Q^{0.574} (Q < 100 \frac{m^3}{h})
    
    """
    line = 'Liquids centrifuge'
    _init = Splitter._init
    _run = Splitter._run
    split = Splitter.split
    isplit = Splitter.isplit
    
class LLECentrifuge(LLEUnit, LiquidsCentrifuge):
    r"""
    Create a liquids centrifuge simulated by liquid-liquid equilibrium.

    Parameters
    ----------
    ins : 
        Inlet fluid.
    outs : 
        * [0] Low density fluid.
        * [1] Heavy fluid.
    top_chemical : str, optional
        Identifier of chemical that will be favored in the low density phase.
    efficiency : float,
        Fraction of feed in liquid-liquid equilibrium.
        The rest of the feed is divided equally between phases.
    
    Notes
    -----
    The f.o.b purchase cost is given by [1]_:
    
    .. math::
        
       C_{f.o.b}^{2007} = 28100 Q^{0.574} (Q < 100 \frac{m^3}{h})
    
    Examples
    --------
    >>> from biorefineries.lipidcane import chemicals
    >>> from biosteam import units, settings, Stream
    >>> settings.set_thermo(chemicals['Methanol', 'Glycerol', 'Biodiesel', 'TAG'])
    >>> feed = Stream('feed', T=333.15,
    ...               TAG=0.996, Biodiesel=26.9,
    ...               Methanol=32.9, Glycerol=8.97)
    >>> C1 = units.LLECentrifuge('C1', ins=feed, outs=('light', 'heavy'))
    >>> C1.simulate()
    >>> C1.outs[0].show()
    Stream: light from <LLECentrifuge: C1>
    phase: 'l', T: 333.15 K, P: 101325 Pa
    flow (kmol/hr): Methanol   10.2
                    Glycerol   0.0239
                    Biodiesel  26.9
                    TriOlein   0.996
    >>> C1.results()
    Liquids centrifuge                       Units       C1
    Electricity         Power                   kW     17.5
                        Cost                USD/hr     1.37
    Design              Flow rate           m^3/hr     12.5
    Purchase cost       Liquids centrifuge     USD 1.29e+05
    Total purchase cost                        USD 1.29e+05
    Utility cost                            USD/hr     1.37
    
    """
    line = 'Liquids centrifuge'
    

@cost('Flow rate', S=1725.61, units='L/min',
    CE=CEPCI_by_year[2007], cost=849000., n=0.6, kW=0.07, ub=2000., BM=2.03,
    N='Number of centrifuges')
class SLLECentrifuge(Unit):
    """
    Create a SLLECentrifuge object that separates the feed into solid, oil, and
    aqueous phases.
    
    Parameters
    ----------
    ins : 
        feed
    outs : 
        * [0] Low density fluid.
        * [1] Heavy fluid.
        * [2] Solids.
    solids_split : dict[str, float]
        Splits to 2nd outlet stream.
    top_chemical : str, optional
        Identifier of chemical that will be favored in the low density phase.
    efficiency : float, optional
        Fraction of feed in liquid-liquid equilibrium.
        The rest of the feed is divided equally between phases.
        Defaults to 1.0.
    moisture_content : float, optional
        Moisture content of solids. Defaults to 0.5.
    
    Notes
    -----
    Cost algorithm is based on a 3-phase decanter centrifuge from 
    a conventional dry-grind corn ethanol plant that separates
    aqueous, oil, and solid fractions (i.e. DDGS) from the bottoms product 
    of the beer column [2]_.
    
    Examples
    --------
    >>> import biosteam as bst
    >>> bst.settings.set_thermo(['Water', 'Hexane', bst.Chemical('Solids', search_db=False, default=True, phase='s')], cache=True)
    >>> feed = bst.Stream('feed', Water=100, Hexane=100, Solids=10)
    >>> C1 = bst.SLLECentrifuge('C1', feed, ['oil', 'aqueous', 'solids'], top_chemical='Hexane', solids_split={'Solids':1.0})
    >>> C1.simulate()
    >>> C1.show()
    SLLECentrifuge: C1
    ins...
    [0] feed
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Water   100
                        Hexane  100
                        Solids  10
    outs...
    [0] oil
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Water   0.791
                        Hexane  100
    [1] aqueous
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Water   98.7
                        Hexane  0.015
    [2] solids
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Water   0.555
                        Solids  10
    
    >>> C1.results()
    3-Phase decanter centrifuge                       Units       C1
    Electricity         Power                            kW   0.0101
                        Cost                         USD/hr 0.000792
    Design              Flow rate                     L/min      250
    Purchase cost       3-Phase decanter centrifuge     USD 2.88e+05
    Total purchase cost                                 USD 2.88e+05
    Utility cost                                     USD/hr 0.000792

    """
    line = '3-Phase decanter centrifuge'
    _N_ins = 1
    _N_outs = 3
    
    solvent = LLEUnit.solvent
    
    @property
    def solids_split(self):
        return self._solids_isplit.data
    @property
    def solids_isplit(self):
        return self._solids_isplit
    
    def _init(self, solids_split, top_chemical=None, efficiency=1.0,
              moisture_content=0.5):
        # [ChemicalIndexer] Splits to 0th outlet stream.
        self._solids_isplit = self.thermo.chemicals.isplit(solids_split)
        
        #: [str] Identifier of chemical that will be favored in the low density phase.
        self.top_chemical = top_chemical
        
        #: Fraction of feed in liquid-liquid equilibrium.
        #: The rest of the feed is divided equally between phases.
        self.efficiency = efficiency 
        
        #: Moisture content of retentate
        self.moisture_content = moisture_content
        assert self._solids_isplit['7732-18-5'] == 0, 'cannot define water split, only moisture content'

    def _run(self):
        feed = self.ins[0]
        top, bottom, solids = self.outs
        feed.split_to(solids, top, self.solids_split, energy_balance=False)
        sep.lle(top, top, bottom, self.top_chemical, self.efficiency)
        sep.adjust_moisture_content(solids, bottom, self.moisture_content)


@copy_algorithm(SLLECentrifuge, run=False)
class SolidLiquidsSplitCentrifuge(Unit):
    """
    Create a SolidLiquidsSplitCentrifuge object that separates the feed into solid, oil, and
    aqueous phases.
    
    Parameters
    ----------
    ins : 
        feed
    outs : 
        * [0] Low density fluid.
        * [1] Heavy fluid.
        * [2] Solids.
    aqueous_split : dict[str, float]
        Splits to [0] outlet stream.
    solids_split : dict[str, float]
        Splits to [2] outlet stream.
    moisture_content : float, optional
        Moisture content of solids. Defaults to 0.5.
    
    Notes
    -----
    Cost algorithm is based on a 3-phase decanter centrifuge from 
    a conventional dry-grind corn ethanol plant that separates
    aqueous, oil, and solid fractions (i.e. DDGS) from the bottoms product 
    of the beer column [2]_.
    
    The unit operation first splits the feed to the light and heavy fractions, 
    then fractionates the solids from the heavy phase.
    
    Examples
    --------
    >>> import biosteam as bst
    >>> bst.settings.set_thermo(['Water', 'Hexane', bst.Chemical('Solids', search_db=False, default=True, phase='s')], cache=True)
    >>> feed = bst.Stream('feed', Water=100, Hexane=100, Solids=10)
    >>> C1 = bst.SolidLiquidsSplitCentrifuge(
    ...     'C1', feed, ['oil', 'aqueous', 'solids'],
    ...     solids_split={'Solids':1.0},
    ...     aqueous_split={'Water':0.99, 'Hexane':0.01, 'Solids': 0.9}
    ... )
    >>> C1.simulate()
    >>> C1.show(flow='kg/hr')
    SolidLiquidsSplitCentrifuge: C1
    ins...
    [0] feed
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kg/hr): Water   1.8e+03
                      Hexane  8.62e+03
                      Solids  10
    outs...
    [0] oil
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kg/hr): Water   18
                      Hexane  8.53e+03
                      Solids  1
    [1] aqueous
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kg/hr): Water   1.77e+03
                      Hexane  86.2
    [2] solids
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kg/hr): Water   9
                      Solids  9
    
    >>> C1.results()
    3-Phase decanter centrifuge                       Units       C1
    Electricity         Power                            kW   0.0101
                        Cost                         USD/hr 0.000792
    Design              Flow rate                     L/min      250
    Purchase cost       3-Phase decanter centrifuge     USD 2.88e+05
    Total purchase cost                                 USD 2.88e+05
    Utility cost                                     USD/hr 0.000792
    
    """
    line = SLLECentrifuge.line
    _N_ins = 1
    _N_outs = 3
    
    @property
    def solids_split(self):
        return self._solids_isplit.data
    @property
    def solids_isplit(self):
        return self._solids_isplit
    
    @property
    def aqueous_split(self):
        return self._aqueous_isplit.data
    @property
    def aqueous_isplit(self):
        return self._aqueous_isplit
    
    def _init(self, aqueous_split, solids_split, moisture_content=0.5):
        # [ChemicalIndexer] Splits to 1st outlet stream (aqueous/heavy phase).
        self._aqueous_isplit = self.thermo.chemicals.isplit(aqueous_split)
        
        # [ChemicalIndexer] Splits to 2th outlet stream (solids) from aqueous/heavy phase.
        self._solids_isplit = self.thermo.chemicals.isplit(solids_split)
        
        #: Moisture content of retentate
        self.moisture_content = moisture_content
        assert self._solids_isplit['7732-18-5'] == 0, 'cannot define water split to solids, only moisture content'

    def _run(self):
        oil, aqueous, solids = self.outs
        self.ins[0].split_to(aqueous, oil, self.aqueous_split)
        aqueous.split_to(solids, aqueous, self.solids_split)
        sep.adjust_moisture_content(solids, aqueous, self.moisture_content)
        
        

# %% Mixing

# Cost base on table 16.32 of Seider's Product and Process Design Principles, 3rd edition
@cost('Power', 'Turbine agitator', N='Number of agitators',
      ub=60, CE=567, cost=3730, n=0.54, BM=2.25)
class LiquidsMixingTank(bst.Unit, PressureVessel):
    """
    Create a LiquidsMixingTank for mixing two liquid phases.
    
    Parameters
    ----------
    ins : 
        Inlet fluids to be mixed.
    outs : 
        Mixed outlet fluid.
    tau=0.022 : float
        Residence time [hr].
    agitator_kW_per_m3=1.0 : float
        Electricity consumption in kW / m3 of volume.
    vessel_material='Carbon steel' : str, optional
        Vessel construction material.
    vessel_type='Horizontal': 'Horizontal' or 'Vertical', optional
        Vessel type.
    length_to_diameter=1 : float
        Length to diameter ratio.
        
    """
    _units = {**PressureVessel._units,
              'Volume': 'm^3',
              'Power': 'hp'}
    _ins_size_is_fixed = False
    _N_ins = 3
    _N_outs = 1
    
    def _init(self, tau=0.022, agitator_kW_per_m3=1.0, length_to_diameter=1,
              vessel_material='Carbon steel', vessel_type='Vertical'):
        self.length_to_diameter = length_to_diameter
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type
        self.agitator_kW_per_m3 = agitator_kW_per_m3
        self.tau = tau
        
    def _run(self):
        self.outs[0].mix_from(self.ins)
        
    def _design(self):
        results = self.design_results
        results['Volume'] = volume = self.tau * self.outs[0].F_vol
        P = self.ins[0].get_property('P', 'psi')
        length_to_diameter = self.length_to_diameter
        results = self.design_results
        rate = self.agitator_kW_per_m3 * volume
        self.power_utility(rate)
        results['Power'] = 1.341 * rate # in hp
        D = geometry.cylinder_diameter_from_volume(volume, length_to_diameter)
        L = length_to_diameter * D
        results.update(self._vessel_design(P, D, L))
        
    def _cost(self):
        self._decorated_cost()
        D = self.design_results
        self.purchase_costs.update(
            self._vessel_purchase_cost(D['Weight'], D['Diameter'], D['Length'])
        )         
   

# %% Settling

class LiquidsSettler(bst.Unit, PressureVessel, isabstract=True):
    """
    Abstract Settler class for liquid-liquid extraction.
    
    Parameters
    ----------
    ins : 
        Inlet fluid with two liquid phases.
    outs : 
        * [0] Low density fluid.
        * [1] Heavy fluid.
    vessel_material='Carbon steel' : str, optional
        Vessel construction material.
    vessel_type='Horizontal': 'Horizontal' or 'Vertical', optional
        Vessel type.
    length_to_diameter=4 : float
        Length to diameter ratio.
    area_to_feed=0.1 : float
        Diameter * length per gpm of feed [ft2/gpm].
        
    """
    _N_ins = 1
    _N_outs = 2
    
    def _init(self, area_to_feed=0.1, 
              length_to_diameter=4,
              vessel_material='Carbon steel',
              vessel_type='Horizontal'):
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type
        self.length_to_diameter = length_to_diameter #: Length to diameter ratio
        self.area_to_feed = area_to_feed #: [ft2/gpm] Diameter * length per gpm of feed
    
    @staticmethod
    def _default_vessel_type():
        return 'Horizontal'
    
    def _design(self):
        feed = self.ins[0]
        F_vol_gpm = feed.get_total_flow('gpm')
        area = self.area_to_feed * F_vol_gpm
        length_to_diameter = self.length_to_diameter
        P = feed.get_property('P', 'psi')
        D = (area / length_to_diameter) ** 0.5
        L = length_to_diameter * D
        self.design_results.update(self._vessel_design(P, D, L))
        
    def _cost(self):
        D = self.design_results
        self.purchase_costs.update(
            self._vessel_purchase_cost(D['Weight'], D['Diameter'], D['Length'])
        )
    

class LLESettler(LLEUnit, LiquidsSettler):
    """
    Create a LLESettler object that rigorously simulates liquid-liquid extraction.
    
    Parameters
    ----------
    ins : 
        Inlet fluid with two liquid phases.
    outs : 
        * [0] Low density fluid.
        * [1] Heavy fluid.
    vessel_material='Carbon steel' : str, optional
        Vessel construction material.
    vessel_type='Horizontal': 'Horizontal' or 'Vertical', optional
        Vessel type.
    length_to_diameter=4 : float, optional
        Length to diameter ratio.
    area_to_feed=0.1 : float, optional
        Diameter * length per gpm of feed [ft2/gpm].
    top_chemical=None : str, optional
        Identifier of chemical that will be favored in the low density phase.
    efficiency=1.0 : float
        Fraction of feed in liquid-liquid equilibrium
    cache_tolerance=1e-6 : float, optional
        Molar tolerance of cached partition coefficients.
    
    """
    line = 'Settler'
    def _init(self, 
            area_to_feed=0.1, 
            length_to_diameter=4,
            vessel_material='Carbon steel',
            vessel_type='Horizontal',
            top_chemical=None,
            efficiency=1.0,
            cache_tolerance=1e-6,
        ):
        LLEUnit._init(self, top_chemical, efficiency)
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type
        self.length_to_diameter = length_to_diameter
        self.area_to_feed = area_to_feed
        self.cache_tolerance = cache_tolerance
        
        
class LiquidsSplitSettler(LiquidsSettler):
    """
    Create a LLESettler object that rigorously simulates liquid-liquid extraction.
    
    Parameters
    ----------
    ins : 
        Inlet fluid with two liquid phases.
    outs : 
        * [0] Low density fluid.
        * [1] Heavy fluid.
    split : Should be one of the following
        * [float] The fraction of net feed in the 0th outlet stream
        * [array_like] Componentwise split of feed to 0th outlet stream
        * [dict] ID-split pairs of feed to 0th outlet stream
    order=None : Iterable[str], defaults to biosteam.settings.chemicals.IDs
        Chemical order of split.
    vessel_material='Carbon steel' : str, optional
        Vessel construction material.
    vessel_type='Horizontal': 'Horizontal' or 'Vertical', optional
        Vessel type.
    length_to_diameter=4 : float, optional
        Length to diameter ratio.
    area_to_feed=0.1 : float, optional
        Diameter * length per gpm of feed [ft2/gpm].
    
    """
    line = 'Settler'
    def _init(self, 
            split, order=None,
            area_to_feed=0.1, 
            length_to_diameter=4,
            vessel_material='Carbon steel',
            vessel_type='Horizontal'
        ):
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type
        self.length_to_diameter = length_to_diameter
        self.area_to_feed = area_to_feed
        self._isplit = self.chemicals.isplit(split, order)        
    split = Splitter.split
    isplit = Splitter.isplit
    _run = Splitter._run
    
    
class LiquidsPartitionSettler(LiquidsSettler):
    """
    Create a LiquidsPartitionSettler object that simulates liquid-liquid 
    extraction by partition coefficients.
    
    Parameters
    ----------
    ins : 
        Inlet fluid with two liquid phases.
    outs : 
        * [0] Low density fluid.
        * [1] Heavy fluid.
    vessel_material='Carbon steel' : str, optional
        Vessel construction material.
    vessel_type='Horizontal': 'Horizontal' or 'Vertical', optional
        Vessel type.
    length_to_diameter=4 : float, optional
        Length to diameter ratio.
    area_to_feed=0.1 : float, optional
        Diameter * length per gpm of feed [ft2/gpm].
    partition_coefficients : 1d array, optional
        Partition coefficients of chemicals in equilibrium (molar 
        composition ratio of the top fluid over the bottom fluid). 
    partition_IDs: tuple[str], optional
        IDs of chemicals in equilibrium.
    
    """
    line = 'Settler'
    def _init(self, 
            partition_coefficients, partion_IDs, 
            area_to_feed=0.1, 
            length_to_diameter=4,
            vessel_material='Carbon steel',
            vessel_type='Horizontal'
        ):
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type
        self.length_to_diameter = length_to_diameter
        self.area_to_feed = area_to_feed
        self.partition_coefficients = partition_coefficients
        self.partion_IDs = partion_IDs
        self.reset_cache()
    
    def reset_cache(self, isdynamic=None):
        self._phi = None
    
    def _run(self):
        self._phi = sep.partition(*self.ins, *self.outs, 
                                  self.partion_IDs, self.partition_coefficients,
                                  self._phi)

# %% Mixer-settlers

class MixerSettler(bst.Unit):
    """
    Create a MixerSettler object that models liquid-liquid extraction using
    a mixing tank and a settler tank.
    
    Parameters
    ----------
    ins : 
        * [0] feed.
        * [1] solvent.
    outs : 
        * [0] extract.
        * [1] raffinate.
    solvent_ID : str, optional
        Name of main chemical in the solvent.
        Defaults to chemical with highest molar fraction in the solvent.
    mixer_data : dict, optional
        Arguments to initialize the "mixer" attribute, a :class:`~biosteam.units.LiquidsMixingTank` object.
    settler_data : dict, optional
        Arguments to initialize the "settler" attribute, a :class:`~biosteam.units.LiquidsSettler` object.
    
    Examples
    --------
    Simulate by rigorous LLE:
    
    >>> import biosteam as bst
    >>> bst.settings.set_thermo(['Water', 'Methanol', 'Octanol'], cache=True)
    >>> feed = bst.Stream('feed', Water=500, Methanol=50)
    >>> solvent = bst.Stream('solvent', Octanol=500)
    >>> MS1 = bst.MixerSettler('MS1', ins=(feed, solvent), outs=('extract', 'raffinate'))
    >>> MS1.simulate()
    >>> MS1.extract.imol['Methanol'] / MS1.feed.imol['Methanol']
    0.66
    >>> MS1.raffinate.imol['Water'] / MS1.feed.imol['Water']
    0.82
    >>> MS1.extract.imol['Octanol'] / MS1.solvent.imol['Octanol']
    0.99
    >>> MS1.results() # doctest: +SKIP
    Mixer settler                                             Units         MS1
    Power               Rate                                     kW        1.98
                        Cost                                 USD/hr       0.155
    Design              Mixer - Volume                          m^3        1.98
                        Mixer - Power                            hp        2.65
                        Mixer - Vessel type                            Vertical
                        Mixer - Length                           ft        1.36
                        Mixer - Diameter                         ft        1.36
                        Mixer - Weight                           lb        91.2
                        Mixer - Wall thickness                   in        0.25
                        Settler - Vessel type                        Horizontal
                        Settler - Length                         ft        12.6
                        Settler - Diameter                       ft        3.15
                        Settler - Weight                         lb    1.44e+03
                        Settler - Wall thickness                 in        0.25
    Purchase cost       Mixer - Turbine agitator                USD    6.32e+03
                        Mixer - Vertical pressure vessel        USD    4.59e+03
                        Mixer - Platform and ladders            USD         641
                        Settler - Horizontal pressure ve...     USD    9.99e+03
                        Settler - Platform and ladders          USD    2.66e+03
    Total purchase cost                                         USD    2.42e+04
    Utility cost                                             USD/hr       0.155
    
    Simulate with user defined partition coefficients:
    
    >>> import biosteam as bst
    >>> import numpy as np
    >>> bst.settings.set_thermo(['Water', 'Methanol', 'Octanol'])
    >>> feed = bst.Stream('feed', Water=500, Methanol=50)
    >>> solvent = bst.Stream('solvent', Octanol=500)
    >>> MS1 = bst.MixerSettler('MS1', 
    ...    ins=(feed, solvent), outs=('extract', 'raffinate'),
    ...    model='partition coefficients',
    ...    settler_data={
    ...        'partition_coefficients': np.array([1.451e-01, 1.380e+00, 2.958e+03]),
    ...        'partion_IDs': ('Water', 'Methanol', 'Octanol'),
    ...    },
    ... )
    >>> MS1.simulate()
    >>> MS1.extract.imol['Methanol'] / MS1.feed.imol['Methanol']
    0.66
    >>> MS1.raffinate.imol['Water'] / MS1.feed.imol['Water']
    0.82
    >>> MS1.extract.imol['Octanol'] / MS1.solvent.imol['Octanol']
    0.99
    >>> MS1.results() # doctest: +SKIP
    Mixer settler                                             Units         MS1
    Electricity         Power                                    kW        1.98
                        Cost                                 USD/hr       0.155
    Design              Mixer - Volume                          m^3        1.98
                        Mixer - Power                            hp        2.65
                        Mixer - Vessel type                            Vertical
                        Mixer - Length                           ft        1.36
                        Mixer - Diameter                         ft        1.36
                        Mixer - Weight                           lb        91.2
                        Mixer - Wall thickness                   in        0.25
                        Settler - Vessel type                        Horizontal
                        Settler - Length                                   12.6
                        Settler - Diameter                                 3.15
                        Settler - Weight                               1.44e+03
                        Settler - Wall thickness                           0.25
    Purchase cost       Mixer - Turbine agitator                USD    6.32e+03
                        Mixer - Vertical pressure vessel        USD    4.91e+03
                        Mixer - Platform and ladders            USD         686
                        Settler - Horizontal pressure ve...     USD    1.16e+04
                        Settler - Platform and ladders          USD    3.08e+03
    Total purchase cost                                         USD    2.65e+04
    Utility cost                                             USD/hr       0.155
    
    """
    _N_ins = 2
    _ins_size_is_fixed = False
    _N_outs = 2
    auxiliary_unit_names = ('mixer', 'settler')
    _graphics = mixer_settler_graphics
    _units = {}
    for i,j in LiquidsMixingTank._units.items(): 
        _units['Mixer - ' + i] = j
    for i,j in LiquidsSettler._units.items(): 
        _units['Settler - ' + i] = j
    
    def _init(self, solvent_ID=None, mixer_data={}, settler_data={}, model="LLE"):
        #: [LiquidsMixingTank] Mixer portion of the mixer-settler.
        #: All data and settings for the design of the mixing tank are stored here.
        self.mixer = mixer = LiquidsMixingTank(None, None, (None,),
                                               self.thermo, **mixer_data)
        self.multi_stream = multi_stream = mixer-0
        mixer._ins = self._ins
        model = model.lower()
        if model == 'lle':
            Settler = LLESettler
        elif model == 'split':
            Settler = LiquidsSplitSettler
        elif model == 'partition coefficients':
            Settler = LiquidsPartitionSettler
            
        #: [LiquidsSettler] Settler portion of the mixer-settler.
        #: All data and settings for the design of the settler are stored here.
        self.settler = Settler(None, multi_stream, thermo=self.thermo, **settler_data)
        
        #: [str] ID of carrier component in the feed.
        self.solvent_ID = solvent_ID 
    
    @property
    def feed(self):
        """[Stream] Feed with solute being extracted and carrier."""
        return self._ins[0]
    @feed.setter
    def feed(self, stream):
        self._ins[0] = stream

    @property
    def solvent(self):
        """[Stream] Solvent to extract solute."""
        return self._ins[1]
    @solvent.setter
    def solvent(self, stream):
        self._ins[1] = stream
        
    @property
    def raffinate(self):
        """[Stream] Raffinate after extraction."""
        return self._outs[1]
    @raffinate.setter
    def raffinate(self, stream):
        self._outs[1] = stream
        
    @property
    def extract(self):
        """[Stream] Extract with solvent."""
        return self._outs[0]
    @extract.setter
    def extract(self, stream):
        self._outs[0] = stream

    def _run(self):
        self.mixer._run()
        self.settler.solvent = self.solvent_ID or self.solvent.main_chemical
        self.settler._run()
        for i, j in zip([self.extract, self.raffinate], self.settler.outs): i.copy_like(j)
        
    def _design(self):
        mixer = self.mixer
        mixer._design()
        settler = self.settler
        settler._design()
        design_results = self.design_results
        for i,j in mixer.design_results.items():
            design_results['Mixer - ' + i] = j
        for i,j in settler.design_results.items():
            design_results['Settler - ' + i] = j
        
    def _cost(self):
        self.mixer._cost()
        self.settler._cost()
        
class MultiStageMixerSettlers(MultiStageEquilibrium):
    """
    Create a MultiStageMixerSettlers object that models a counter-current system
    of mixer-settlers for liquid-liquid extraction.
    
    Parameters
    ----------
    ins : 
        * [0] feed.
        * [1] solvent.
    outs : 
        * [0] extract
        * [1] raffinate
    N_stages : int
        Number of stages.
    partition_data : {'IDs': tuple[str], 'K': 1d array}, optional
        IDs of chemicals in equilibrium and partition coefficients (molar 
        composition ratio of the extract over the raffinate). If given,
        The mixer-settlers will be modeled with these constants. Otherwise,
        partition coefficients are computed based on temperature and composition.
    solvent_ID : str
        Name of main chemical in the solvent.
    mixer_data : dict
        Arguments to initialize the "mixer" attribute, a :class:`~biosteam.units.LiquidsMixingTank` object.
    settler_data : dict
        Arguments to initialize the "settler" attribute, a :class:`~biosteam.units.LiquidsSettler` object.
    
    Notes
    -----
    All mixer settlers are sized equally based on the assumption that the 
    volumetric flow rate of each phase does not change significantly across
    stages. 
    
    Examples
    --------
    Simulate by rigorous LLE:
    
    >>> import biosteam as bst
    >>> bst.settings.set_thermo(['Water', 'Methanol', 'Octanol'])
    >>> feed = bst.Stream('feed', Water=500, Methanol=50)
    >>> solvent = bst.Stream('solvent', Octanol=500)
    >>> MSMS1 = bst.MultiStageMixerSettlers('MSMS1', ins=(feed, solvent), outs=('extract', 'raffinate'), N_stages=2)
    >>> MSMS1.simulate()
    >>> MSMS1.extract.imol['Methanol'] / MSMS1.feed.imol['Methanol']
    0.83
    >>> MSMS1.raffinate.imol['Water'] / MSMS1.feed.imol['Water']
    0.82
    >>> MSMS1.extract.imol['Octanol'] / MSMS1.solvent.imol['Octanol']
    0.99
    >>> MSMS1.results() # doctest: +SKIP
    Multi stage mixer settlers                     Units       MSMS1
    Electricity         Power                         kW        3.96
                        Cost                      USD/hr       0.309
    Design              Mixer - Volume               m^3        1.98
                        Mixer - Power                 hp        2.65
                        Mixer - Vessel type                 Vertical
                        Mixer - Length                ft        1.36
                        Mixer - Diameter              ft        1.36
                        Mixer - Weight                lb        91.2
                        Mixer - Wall thickness        in        0.25
                        Settler - Vessel type             Horizontal
                        Settler - Length                        12.6
                        Settler - Diameter                      3.15
                        Settler - Weight                    1.44e+03
                        Settler - Wall thickness                0.25
    Purchase cost       Mixers and agitators         USD    1.12e+04
                        Settlers                     USD    2.93e+04
    Total purchase cost                              USD    4.05e+04
    Utility cost                                  USD/hr       0.309
    
    Simulate with user defined partition coefficients:
    
    >>> import biosteam as bst
    >>> import numpy as np
    >>> bst.settings.set_thermo(['Water', 'Methanol', 'Octanol'])
    >>> feed = bst.Stream('feed', Water=5000, Methanol=500)
    >>> solvent = bst.Stream('solvent', Octanol=5000)
    >>> MSMS1 = bst.MultiStageMixerSettlers('MSMS1', ins=(feed, solvent), outs=('extract', 'raffinate'), N_stages=10,
    ...     partition_data={
    ...         'K': np.array([0.1450, 1.380, 2957.]),
    ...         'IDs': ('Water', 'Methanol', 'Octanol'),
    ...         'phi': 0.590 # Initial phase fraction guess. This is optional.
    ...     }
    ... )
    >>> MSMS1.simulate()
    >>> MSMS1.extract.imol['Methanol'] / MSMS1.feed.imol['Methanol']
    0.99
    >>> MSMS1.raffinate.imol['Water'] / MSMS1.feed.imol['Water']
    0.82
    >>> MSMS1.extract.imol['Octanol'] / MSMS1.solvent.imol['Octanol']
    0.99
    >>> MSMS1.results() # doctest: +SKIP
    Multi stage mixer settlers                     Units       MSMS1
    Electricity         Power                         kW         198
                        Cost                      USD/hr        15.5
    Design              Mixer - Volume               m^3        19.8
                        Mixer - Power                 hp        26.5
                        Mixer - Vessel type                 Vertical
                        Mixer - Length                ft        2.93
                        Mixer - Diameter              ft        2.93
                        Mixer - Weight                lb         423
                        Mixer - Wall thickness        in        0.25
                        Settler - Vessel type             Horizontal
                        Settler - Length                        39.8
                        Settler - Diameter                      9.95
                        Settler - Weight                    2.52e+04
                        Settler - Wall thickness               0.438
    Purchase cost       Mixers and agitators         USD    1.15e+05
                        Settlers                     USD    6.15e+05
    Total purchase cost                              USD     7.3e+05
    Utility cost                                  USD/hr        15.5

    Because octanol and water do not mix well, it may be a good idea to assume
    that these solvents do not mix at all:
        
    >>> MSMS1 = bst.MultiStageMixerSettlers('MSMS1', ins=(feed, solvent), outs=('extract', 'raffinate'), N_stages=20,
    ...     partition_data={
    ...         'K': np.array([1.38]),
    ...         'IDs': ('Methanol',),
    ...         'raffinate_chemicals': ('Water',),
    ...         'extract_chemicals': ('Octanol',),
    ...     }
    ... )
    >>> MSMS1.simulate()
    >>> MSMS1.extract.imol['Methanol'] / feed.imol['Methanol'] # Recovery
    0.99
    >>> MSMS1.extract.imol['Octanol'] / solvent.imol['Octanol'] # Solvent stays in extract
    1.0
    >>> MSMS1.raffinate.imol['Water'] / feed.imol['Water'] # Carrier remains in raffinate
    1.0
       
    Simulate with a feed at the 4th stage:
    
    >>> dilute_feed = bst.Stream('dilute_feed', Water=100, Methanol=2)
    >>> MSMS1 = bst.MultiStageMixerSettlers('MSMS1', ins=(feed, dilute_feed, solvent), outs=('extract', 'raffinate'), N_stages=5,
    ...     feed_stages=[0, 3, -1], # Stage at which each inlet enters, respectively
    ...     partition_data={
    ...         'K': np.array([1.38]),
    ...         'IDs': ('Methanol',),
    ...         'raffinate_chemicals': ('Water',),
    ...         'extract_chemicals': ('Octanol',),
    ...     }
    ... )
    >>> MSMS1.simulate()
    >>> MSMS1.extract.imol['Methanol'] / (feed.imol['Methanol'] + dilute_feed.imol['Methanol']) # Recovery
    0.93
    
    Simulate with a 60% extract side draw at the 2nd stage and 10% raffinate side draw at the 3rd stage:
    
    >>> MSMS1 = bst.MultiStageMixerSettlers('MSMS1', ins=(feed, solvent), N_stages=4,  
    ...     # Extract side draws always come first, then raffinate side draws
    ...     outs=('extract', 'raffinate', 'extract_side_draw', 'raffinate_side_draw'),  
    ...     extract_side_draws=[(1, 0.6)], # Stage number and split fraction pairs
    ...     raffinate_side_draws=[(2, 0.10)], 
    ...     partition_data={
    ...         'K': np.array([1.38]),
    ...         'IDs': ('Methanol',),
    ...         'raffinate_chemicals': ('Water',),
    ...         'extract_chemicals': ('Octanol',),
    ...     }
    ... )
    >>> MSMS1.simulate()
    >>> (MSMS1.extract.imol['Methanol'] + MSMS1.outs[2].imol['Methanol']) / feed.imol['Methanol'] # Recovery
    1.0
    

    """
    _units = MixerSettler._units
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, 
                 feed_stages=None, extract_side_draws=None, raffinate_side_draws=None, 
                 **kwargs):
        if extract_side_draws is None: extract_side_draws = {}
        elif not isinstance(extract_side_draws, dict): extract_side_draws = dict(extract_side_draws)
        if raffinate_side_draws is None: raffinate_side_draws = {}
        elif not isinstance(raffinate_side_draws, dict): raffinate_side_draws = dict(raffinate_side_draws)
        if feed_stages is None: feed_stages = (0, -1)
        self._N_ins = len(feed_stages)
        self._N_outs = 2 + len(extract_side_draws) + len(raffinate_side_draws)
        Unit.__init__(self, ID, ins, outs, thermo, 
                      feed_stages=feed_stages,
                      extract_side_draws=extract_side_draws, 
                      raffinate_side_draws=raffinate_side_draws,
                      **kwargs)
    
    def _init(self, N_stages, feed_stages, extract_side_draws, 
              raffinate_side_draws, partition_data=None, solvent_ID=None,  
              mixer_data={}, settler_data={}, use_cache=None):
        bst.MultiStageEquilibrium._init(
            self, N_stages=N_stages, feed_stages=feed_stages, phases=('l', 'L'), P=101325,
            top_side_draws=extract_side_draws, bottom_side_draws=raffinate_side_draws,
            stage_specifications=None, partition_data=partition_data, 
            solvent=solvent_ID, use_cache=use_cache,
        )
        #: [LiquidsMixingTank] Used to design all mixing tanks. 
        #: All data and settings for the design of mixing tanks are stored here.
        self.mixer = mixer = LiquidsMixingTank(None, None, (None,),
                                               self.thermo, **mixer_data)
        mixer._ins = self._ins
        
        #: [LiquidsSettler] Used to design all settlers.
        #: All data and settings for the design of settlers are stored here.
        self.settler = LiquidsSettler(None, mixer-0, None, self.thermo, **settler_data)
        self.settler._outs = self._outs
        self.use_cache = use_cache
        self._last_args = (
            self.N_stages, self.feed_stages, self.extract_side_draws, self.use_cache,
            *self._ins, self.raffinate_side_draws, self.solvent_ID, self.partition_data
        )
        self.reset_cache()
    
    feed = MixerSettler.feed
    solvent = MixerSettler.solvent
    
    @property
    def partition_data(self):
        return self._partition_data
    @partition_data.setter
    def partition_data(self, partition_data):
        self._partition_data = partition_data
        self._last_args = None
    
    def _setup(self):
        super()._setup()
        args = (self.N_stages, self.feed_stages, self.extract_side_draws, self.use_cache,
                *self._ins, self.raffinate_side_draws, self.solvent_ID, self.partition_data)
        if args != self._last_args:
            MultiStageEquilibrium._init(
                self, N_stages=self.N_stages, feed_stages=self.feed_stages,
                phases=('l', 'L'), P=self.P,
                top_side_draws=self.extract_side_draws, 
                bottom_side_draws=self.raffinate_side_draws,
                stage_specifications=None, partition_data=self.partition_data, 
                solvent=self.solvent_ID, use_cache=self.use_cache, 
            )
            self.mixer._ins = self._ins
            self.settler._outs = self._outs
            self._last_args = args
    
    def reset_cache(self, isdynamic=None):
        self._last_args = None
        
    def _design(self):
        mixer = self.mixer
        mixer._run()
        mixer._design()
        settler = self.settler
        settler._design()
        design_results = self.design_results
        for i,j in mixer.design_results.items():
            design_results['Mixer - ' + i] = j
        for i,j in settler.design_results.items():
            design_results['Settler - ' + i] = j
        
    def _cost(self):
        N_stages = self.N_stages
        mixer = self.mixer
        settler = self.settler
        mixer._cost()
        settler._cost()
        self.power_utility.copy_like(mixer.power_utility)
        self.power_utility.scale(N_stages)
        purchase_costs = self.purchase_costs
        purchase_costs['Mixers and agitators'] = N_stages * mixer.purchase_cost
        purchase_costs['Settlers'] = N_stages * settler.purchase_cost
        baseline_purchase_costs = self.baseline_purchase_costs
        baseline_purchase_costs['Mixers and agitators'] = N_stages * mixer.purchase_cost
        baseline_purchase_costs['Settlers'] = N_stages * settler.purchase_cost
        installed_costs = self.installed_costs
        installed_costs['Mixers and agitators'] = N_stages * mixer.purchase_cost
        installed_costs['Settlers'] = N_stages * settler.purchase_cost
        
    @property
    def installed_cost(self):
        N_stages = self.N_stages
        return N_stages * (self.mixer.installed_cost + self.settler.installed_cost)
        
