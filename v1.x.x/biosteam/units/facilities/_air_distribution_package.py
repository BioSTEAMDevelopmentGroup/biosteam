# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 21:35:10 2019

@author: yoelr
"""

from .. import Static
from ..decorators import cost

__all__ = ('AirDistributionPackage',)

@cost('Flow rate', 'Plant air reciever',
      cost=16e3, CE=522, S=83333, n=0.6, BM=3.1)
@cost('Flow rate', 'Instrument air dryer',
      cost=15e3, CE=522, S=83333, n=0.6, BM=1.8)
@cost('Flow rate', 'Plant air compressor', units='kg/hr',
      cost=28e3, CE=551, S=83333, n=0.6, BM=1.6, kW=150*0.7457)
class AirDistributionPackage(Static): pass
    