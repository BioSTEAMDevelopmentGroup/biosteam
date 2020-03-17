# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 21:35:10 2019

@author: yoelr
"""
from ..._facility import Facility
from ..decorators import cost

__all__ = ('AirDistributionPackage',)

@cost('Flow rate', 'Plant air reciever',
      cost=16e3, CE=522, S=83333, n=0.6, BM=3.1)
@cost('Flow rate', 'Instrument air dryer',
      cost=15e3, CE=522, S=83333, n=0.6, BM=1.8)
@cost('Flow rate', 'Plant air compressor', units='kg/hr',
      cost=28e3, CE=551, S=83333, n=0.6, BM=1.6, kW=150*0.7457)
class AirDistributionPackage(Facility):
    """
    Create a AirDistributionPackage object that accounts for the capital cost 
    and power of air distribution based on flow rate correlations from [1]_.
    
    Parameters
    ----------
    ID : str, optional
        Unit ID.
    ins : stream
        Air flow to be distributed.
    outs : stream
        Distributed air flow.
    
    References
    ----------
    .. [1] Humbird, D., Davis, R., Tao, L., Kinchin, C., Hsu, D., Aden, A.,
        Dudgeon, D. (2011). Process Design and Economics for Biochemical 
        Conversion of Lignocellulosic Biomass to Ethanol: Dilute-Acid 
        Pretreatment and Enzymatic Hydrolysis of Corn Stover
        (No. NREL/TP-5100-47764, 1013269). https://doi.org/10.2172/1013269
    """
    network_priority = 0