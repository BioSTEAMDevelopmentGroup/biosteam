# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 21:23:56 2018

@author: yoelr
"""
from .. import Unit
from ._flash import RatioFlash
from ._splitter import Splitter
from .decorators import cost

__all__ = ('LiquidsCentrifuge', 'LiquidsSplitCentrifuge', 'LiquidsRatioCentrifuge')

# electricity kW/(m3/hr) from USDA biosdiesel Super Pro model
# Possibly 1.4  kW/(m3/hr)
# https://www.sciencedirect.com/topics/engineering/disc-stack-centrifuge
# Microalgal fatty acids—From harvesting until extraction H.M. Amaro, et. al.,
# in Microalgae-Based Biofuels and Bioproducts, 2017

@cost('Flow rate', units='m^3/hr', CE=525.4, cost=28100,
      n=0.574, kW=3.66, ub=100, BM=2.03, N='Number of centrifuges')
class LiquidsCentrifuge(Unit, isabstract=True):
    r"""Create a liquids centrifuge.

    Equation for f.o.b cost:

    :math:`C_{fob}^{2007} = 28100 Q^{0.574} (0.1 < Q < 100 \frac{m^3}{h})` 

    Parameters
    ----------
    ins
        [0] Input stream
        
    outs
        [0] 'liquid' phase stream
        
        [1] 'LIQUID' phase stream
    
    """
    _N_outs = 2
    _bounds = {'Flow rate': (0.1, 100)}


class LiquidsRatioCentrifuge(LiquidsCentrifuge):
    _N_heat_utilities = 0
    __init__ = RatioFlash.__init__
    _run = RatioFlash._run


class LiquidsSplitCentrifuge(LiquidsCentrifuge):
    __init__ = Splitter.__init__
    _run = Splitter._run
    split = Splitter.split
    
    
    