# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
This module contains unit operations for size enlargement.

.. contents:: :local:
    
Unit operations
---------------
.. autoclass:: biosteam.units.PelletMill

.. autoclass:: biosteam.units.BagassePelletMill

"""
from .decorators import cost
from .._unit import Unit
from warnings import warn

__all__ = ('PelletMill', 'BagassePelletMill')

# TODO: The pellet mill is carbon steel. Add material factors later
@cost('Flow rate', CE=567., lb=362.8739, ub=36287.39,
      S=1., units='kg/hr', n=0.11, cost=8659.2, BM=2.3)
class PelletMill(Unit):
    """
    Create a pellet mill for pelletizing solids.
    
    Parameters
    ----------
    ins : stream
    outs : stream 
    kW : float, optional
        Power requirement [kWh/kg]. Defaults to 0.01 based on a heuristic from 
        Perry's handbook.
    
    Notes
    -----
    Feed is pushed through holes in dies of various shapes. The friction of 
    the material in the die holes supplies the resistance necessary for 
    compaction. Adjustable knives shear the rodlike extrudates into pellets of 
    the desired length. Although several designs are in use, the most commonly 
    used pellet mills operate by applying power to the die and rotating it 
    around a freely turning roller with fixed horizontal or vertical axis. 
    
    """
    _N_ins = 1
    _N_outs = 1
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *, kW=0.01):
        Unit.__init__(self, ID, ins, outs, thermo)
        self.kW = kW
        
    def _cost(self):
        self._decorated_cost()
        self.power_utility.consumption = self.kW * self.design_results['Flow rate']
        

class BagassePelletMill(PelletMill):
    """
    Create a pellet mill for pelletizing dryed bagasse.
    
    Parameters
    ----------
    ins : stream
    outs : stream 
    kW : float, optional
        User enforced power requirement [kWh/kg]. Defaults according to
        the `capacity_to_power` class attribute.
    
    Notes
    -----
    Feed is pushed through holes in dies of various shapes. The friction of 
    the material in the die holes supplies the resistance necessary for 
    compaction. Adjustable knives shear the rodlike extrudates into pellets of 
    the desired length. Although several designs are in use, the most commonly 
    used pellet mills operate by applying power to the die and rotating it 
    around a freely turning roller with fixed horizontal or vertical axis. 
    
    """
    
    #: Tuple[Tuple[float, float], float] Capacity ranges [kg/hr] and default 
    #: power requirements [kW]. Values from industrial supplier:
    #: http://www.biopelletmill.com/large-biomass-pellet-mill.html
    #: http://www.biopelletmill.com/small-biomass-pellet-mill.html
    capacity_to_power = [
        ((50, 100), 5.5),
        ((80, 120), 7.5),
        ((120, 200), 11),
        ((160, 250), 15),
        ((250, 400), 22),
        ((400, 800), 61.5),
        ((800, 1500), 94.5),
        ((1500, 2000), 114.5),
    ]
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *, kW=None):
        Unit.__init__(self, ID, ins, outs, thermo)
        self.kW = kW
        
    def default_kW(self, F_mass):
        for (lb, ub), kW in self.capacity_to_power:
            if F_mass < ub: return kW
        return F_mass * kW / ub # Extrapolate
        
    def _cost(self):
        self._decorated_cost()
        F_mass = self.design_results['Flow rate']
        self.power_utility.consumption = self.default_kW(F_mass)
        
        

    