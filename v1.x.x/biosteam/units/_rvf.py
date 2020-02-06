# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 22:15:20 2018

@author: yoelr
"""
from ._solids_separator import SolidsSeparator
from .designtools import vacuum_system
import numpy as np
import biosteam as bst

__all__ = ('RotaryVacuumFilter', 'RVF')

class RotaryVacuumFilter(SolidsSeparator):
    """Create a RotaryVacuumFilter object.
    
    Parameters
    ----------
    ins
        [0] Feed
        
        [1] Wash water
    outs
        [0] Retentate
        
        [1] Permeate
    moisture_content : float
                       Fraction of water in retentate.
    
    """
    _has_power_utility = True
    BM = 2.32
    
    #: Revolutions per second
    rps = 20/3600
    
    #: Radius of the vessel (m)
    radius = 1
    
    #: Suction pressure (Pa)
    P_suction = 105
    
    #: For crystals (lb/day-ft^2)
    filter_rate = 6000
    _kwargs = {'moisture_content': 0.80} # fraction
    _bounds = {'Individual area': (10, 800)}
    _units = {'Area': 'ft^2',
              'Individual area': 'ft^2'}
    
    #: Efficiency of the vacuum pump
    power_efficiency = 0.9
    
    def _design(self):
        flow = sum(stream.massnet for stream in self.outs)
        self._Design['Area'] = self._calc_Area(flow, self.filter_rate)
        
    def _cost(self):
        Design = self._Design
        Area = Design['Area']
        ub = self._bounds['Individual area'][1]
        N_vessels = np.ceil(Area/ub)
        self._power(Area, N_vessels)
        iArea = Area/N_vessels # individual vessel
        Design['# RVF'] = N_vessels
        Design['Individual area'] = iArea
        logArea = np.log(iArea)
        Cost = np.exp(11.796-0.1905*logArea+0.0554*logArea**2)
        self._Cost['Cost of vessels'] = N_vessels*Cost*bst.CE/567
    
    def _power(self, area, N_vessels):
        s_cake, s_vacuumed = self.outs
        
        # # Weight of empty plate
        # mass_plates = 10*N_vessels
        
        # # Revolutions per s
        # rps = self.rps
        
        # # Cake volumetric flow meter per sec
        # Volcake = s_cake.volnet / 3600
        
        # # Thickness of cake layer, assumed uniform and constant
        # thCake = Volcake / rps;
        
        # # Mass of Cake
        # mass_cake = thCake * s_cake.rho
        radius = self.radius
        # cent_a = (2*np.pi*rps)**2*radius
        # cent_F = (mass_cake + mass_plates)*cent_a
        # work_rot = rps*2*np.pi*radius*cent_F
        Area = self._Design['Area']
        vol = radius*Area*0.0929/2 # m3
        
        # Assume same volume of air comes in as volume of liquid
        volflow = s_vacuumed.volnet
        massflow = volflow*1.2041 # multiply by density of air kg/m3 
        work_vacuum, self._Cost['Liquid-ring pump'] = vacuum_system(
                massflow, volflow, self.P_suction, vol)
        #power = work_rot/self.power_efficiency/1000 + work_vacuum # kW
        self._power_utility(work_vacuum)
    
    @staticmethod
    def _calc_Area(flow, filter_rate):
        """Return area in ft^2 given flow in kg/hr and filter rate in lb/day-ft^2."""
        return flow*52.91/filter_rate
    
        
RVF = RotaryVacuumFilter


