# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 22:15:20 2018

@author: yoelr
"""
from biosteam import Unit
from .splitter import Splitter
from .designtools import vacuum_system
import numpy as np

class RotaryVacuumFilter(Unit):
    _has_power_utility = True
    
    #: Revolutions per second
    rps = 20/3600
    
    #: Radius of the vessel (m)
    radius = 1
    
    #: For crystals (lb/day-ft^2)
    filter_rate = 6000
    _kwargs = {'split': None,
               'P_suction': 100} # Pa
    _bounds = {'Individual area': (10, 800)}
    _units = {'Area': 'ft^2',
              'Individual area': 'ft^2'}
    
    #: Efficiency of the vacuum pump
    power_efficiency = 0.9
    
    @staticmethod
    def _calc_Area(flow:'kg/hr', filter_rate:'lb/day-ft^2') -> 'ft^2':
        flow *= 52.91
        return flow/filter_rate
        
    _run = Splitter._run
    
    def _design(self):
        flow = sum(stream.massnet for stream in self.outs)
        Area = self._calc_Area(flow, self.filter_rate)
        self._results['Design']['Area'] = Area
        
    def _cost(self):
        results = self._results
        Design = results['Design']
        Area = Design['Area']
        ub = self._bounds['Individual area'][1]
        N_vessels = np.ceil(Area/ub)
        self._power(Area, N_vessels)
        iArea = Area/N_vessels # individual vessel
        Design['Individual area'] = iArea
        
        Cost = np.exp(11.796-0.1905*np.log(iArea)+0.0554*(np.log(iArea))**2)
        results['Cost']['Cost of vessels'] = N_vessels*Cost*self.CEPCI/567
    
    def _power(self, area, N_vessels) :
        r = self._results
        s_cake, s_vacuumed = self.outs
        P_suction = self._kwargs['P_suction'] # Max allowable pressure drop
        
        # Weight of empty plate
        mass_plates = 10*N_vessels
        
        # Revolutions per s
        rps = self.rps
        
        # Cake volumetric flow meter per sec
        Volcake = s_cake.volnet / 3600
        
        # Thickness of cake layer, assumed uniform and constant
        thCake = Volcake / rps;
        
        # Mass of Cake
        mass_cake = thCake * s_cake.rho
        radius = self.radius
        cent_a = (2*np.pi*rps)**2*radius
        cent_F = (mass_cake + mass_plates)*cent_a
        work_rot = rps*2*np.pi*radius*cent_F
        Area = r['Design']['Area']
        vol = radius*Area*0.0929/2 # m3
        
        # Assume same volume of air comes in as volume of liquid
        volflow = s_vacuumed.volnet
        massflow = volflow*1.2041 # multiply by density of air kg/m3 
        work_vacuum, r['Cost']['Liquid-ring pump'] = vacuum_system(massflow,
                                                                   volflow,
                                                                   P_suction, vol)
        power = work_rot/self.power_efficiency/1000 + work_vacuum # kW
        self._power_utility(power)
        
        
RVF = RotaryVacuumFilter


