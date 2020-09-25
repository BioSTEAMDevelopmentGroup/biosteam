# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import biosteam as bst
from ._lle_unit import LLEUnit
from ._splitter import Splitter
from .design_tools import PressureVessel
from thermosteam import separations

__all__ = ('LiquidsSettler', 'LLESettler', 'LiquidsSplitSettler')

class LiquidsSettler(bst.Unit, PressureVessel, isabstract=True):
    """
    Abstract Settler class for liquid-liquid extraction.
    
    Parameters
    ----------
    ins : stream
        Inlet fluid with two liquid phases.
    outs : stream sequence
        [0] Top fluid.
        [1] Bottom fluid.
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
    _N_heat_utilities = 0
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 area_to_feed=0.1, 
                 length_to_diameter=4,
                 vessel_material='Carbon steel',
                 vessel_type='Horizontal'):
        bst.Unit.__init__(self, ID, ins, outs, thermo)
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
    ins : stream
        Inlet fluid with two liquid phases.
    outs : stream sequence
        [0] Top fluid.
        [1] Bottom fluid.
    vessel_material='Carbon steel' : str, optional
        Vessel construction material.
    vessel_type='Horizontal': 'Horizontal' or 'Vertical', optional
        Vessel type.
    length_to_diameter=4 : float, optional
        Length to diameter ratio.
    area_to_feed=0.1 : float, optional
        Diameter * length per gpm of feed [ft2/gpm].
    top_chemical=None : str, optional
        Chemical selectively partitioned to the top phase
    efficiency=1.0 : float
        Fraction of feed in liquid-liquid equilibrium
    cache_tolerance=1e-6 : float, optional
        Molar tolerance of cached partition coefficients.
    
    """
    line = 'Settler'
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 area_to_feed=0.1, 
                 length_to_diameter=4,
                 vessel_material='Carbon steel',
                 vessel_type='Horizontal',
                 top_chemical=None,
                 efficiency=1.0,
                 cache_tolerance=1e-6,
        ):
        LLEUnit.__init__(self, ID, ins, outs, thermo, top_chemical, efficiency)
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
    ins : stream
        Inlet fluid with two liquid phases.
    outs : stream sequence
        [0] Top fluid.
        [1] Bottom fluid.
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
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 split, order=None,
                 area_to_feed=0.1, 
                 length_to_diameter=4,
                 vessel_material='Carbon steel',
                 vessel_type='Horizontal'):
        bst.Unit.__init__(self, ID, ins, outs, thermo)
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
    ins : stream
        Inlet fluid with two liquid phases.
    outs : stream sequence
        [0] Top fluid.
        [1] Bottom fluid.
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
    forced_split : 1d array, optional
        Componentwise split of feed to 0th outlet stream.
    forced_split_IDs=None : Iterable[str], optional
        Chemical order of `forced_split`.
    
    """
    line = 'Settler'
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 partition_coefficients, partion_IDs, 
                 forced_split_IDs=None, forced_split=None,
                 area_to_feed=0.1, 
                 length_to_diameter=4,
                 vessel_material='Carbon steel',
                 vessel_type='Horizontal'):
        bst.Unit.__init__(self, ID, ins, outs, thermo)
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type
        self.length_to_diameter = length_to_diameter
        self.area_to_feed = area_to_feed
        self.partition_coefficients = partition_coefficients
        self.partion_IDs = partion_IDs
        #: array[float] Forced splits to 0th stream for given IDs. 
        self.forced_split = forced_split
        #: tuple[str] IDs corresponding to forced splits. 
        self.forced_split_IDs = forced_split_IDs
        self.reset_cache()
    
    def reset_cache(self):
        self._phi = None
    
    def _run(self):
        self._phi = separations.partition(*self.ins, *self.outs, 
                                          self.partion_IDs, self.partition_coefficients,
                                          self._phi)