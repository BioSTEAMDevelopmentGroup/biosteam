# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 22:17:05 2018

@author: yoelr
"""

from biosteam.units.splitter import Splitter
from biosteam import PowerUtility, Unit

__all__ = ['MolecularSieve', 'MolSieve']

class MolecularSieve(Unit):
    """Create an ethanol/water molecular sieve for bioethanol plants. The molecular sieve is modeled as a component wise separator. Costing is based on scaling by the 6/10ths rule from an NREL TEA report [1].
    
    **Parameters**
    
        **split:** *[array_like]* Component wise split in the 0th output stream

    **ins**
    
        [:] Input streams
        
    **outs**
    
        [0] Split stream
        
        [1] Remainder stream
    
    **References**
    
        [1] Process Design and Economics for Biochemical Conversion of Lignocellulosic Biomass to Ethanol Dilute-Acid Pretreatment and Enzymatic Hydrolysis of Corn Stover. D. Humbird, R. Davis, L. Tao, C. Kinchin, D. Hsu, and A. Aden (National Renewable Energy Laboratory Golden, Colorado). P. Schoen, J. Lukas, B. Olthof, M. Worley, D. Sexton, and D. Dudgeon (Harris Group Inc. Seattle, Washington and Atlanta, Georgia)
    
    **Examples**
    
        :doc:`MolecularSieve Example`
    
    """
    _kwargs = Splitter._kwargs
    _run = Splitter._run
    _has_power_utility = True
    
    #: Original Price (USD)
    C_0 = 2601000
    
    #: Original flow rate (kg/hr)
    V_0 = 22687 
    
    #: Scaling exponent
    exp = 0.60 
    
    def _cost(self):
        results = self._results
        Cost = results['Cost']
        CEPCI_old = 521.9 # Original CEPCI (2009)
        CEPCI_new = self.CEPCI
        F_CEPCI = CEPCI_new/CEPCI_old
        
        # Cost molecular sieve vessels with packing
        Cp_old = self.C_0 # Original purchase price (USD)
        flow_old = self.V_0 # Original flow rate (kg/hr)
        flow_new = self._massnet_in
        x = self.exp # Scaling exponent
        S = flow_new/flow_old # Size ratio
        Cp_new = F_CEPCI * Cp_old * S**x
        Cost['Molecular Sieve Cost'] = Cp_new
        
        # Cost pressure filter (pressing compressor)
        Cp_old = 75200
        Cost['Pressure Filter Pressing'] = F_CEPCI * Cp_old * S**x
        power0 = 112*S # kW
        
        # Cost pressure filter (drying compressor)
        N = 2
        Cp_old = 405000
        Cost['Pressure Filter Drying (2)'] = N * F_CEPCI * Cp_old * S**x
        power1 = N*522*S # kW
        
        # Find power utility
        total_power = power0 + power1
        self._power_utility(total_power)
        return Cost
        
MolSieve = MolecularSieve        
        


