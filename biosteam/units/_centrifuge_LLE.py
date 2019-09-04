# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 21:23:56 2018

@author: yoelr
"""
from .. import Unit, MixedStream
from ._flash import RatioFlash
from ._splitter import Splitter
from .decorators import cost

__all__ = ('Centrifuge_LLE', 'SplitCentrifuge_LLE', 'RatioCentrifuge_LLE')

# electricity kW/(m3/hr) from USDA biosdiesel Super Pro model
# Possibly 1.4  kW/(m3/hr)
# https://www.sciencedirect.com/topics/engineering/disc-stack-centrifuge
# Microalgal fatty acids—From harvesting until extraction H.M. Amaro, et. al.,
# in Microalgae-Based Biofuels and Bioproducts, 2017

@cost('Flow rate', units='m^3/hr', CE=525.4, cost=28100,
      n=0.574, kW=3.66, ub=100, BM=2.03)
class Centrifuge_LLE(Unit):
    r"""Create an equlibrium based centrifuge with the option of having liquid non-keys and LIQUID non-keys completly separate into their respective phases.

    Equation for f.o.b cost:

    :math:`C_{fob}^{2007} = 28100 Q^{0.574} (0.1 < Q < 100 \frac{m^3}{h})` 

    Parameters
    ----------
    ins
        [0] Input stream
        
    outs
        [0] 'liquid' phase stream
        
        [1] 'LIQUID' phase stream
    species_IDs=None : tuple[str], optional
        IDs of equilibrium species
    split=None : tuple[float], optional
        Initial guess split fractions of each specie to the 'liquid'
    lNK=() : tuple[str], optional
        Species assumed to completely remain in the 'liquid' phase.
    LNK=() : tuple[str], optional
        Species assumed to completely remain in the 'LIQUID' phase.
    solvents=() : tuple[str], optional
        Species corresponding to specified solvent_split
    solvent_split=() : tuple[float], optional
        Split fractions of each specie to the 'liquid' phase.                
    
    **Examples**
    
        :doc:`notebooks/Centrifuge_LLE Example`
    
    """
    line = 'Liquids Centrifuge'
    _bounds = {'Flow rate': (0.1, 100)}

    def __init__(self, ID='', ins=None, outs=(),
                 species_IDs=None, split=None,
                 lNK=(), LNK=(), solvents=(),
                 solvent_split=()):
        Unit.__init__(ID, ins, outs)
        self.split = split
        self.species_IDs = species_IDs
        self.lNK = lNK
        self.LNK = LNK
        self.solvents = solvents
        self.solvent_split = solvent_split
        self._mixedstream =  MixedStream(None)

    def _run(self):
        liq, LIQ = self.outs
        feed = self.ins[0]
        ms = self._mixedstream
        ms.empty()
        ms.liquid_mol[:] = feed.mol
        ms.LLE(T=feed.T,  split=self.split, species_IDs=self.species_IDs,
               lNK=self.lNK, LNK=self.LNK, solvents=self.solvents,
               solvent_split=self.solvent_split)
        liq._mol[:] = ms.liquid_mol
        LIQ._mol[:] = ms.LIQUID_mol
        liq.T = LIQ.T = ms.T
        liq.P = LIQ.P = ms.P
        
# class PartitionCentrifuge_LLE(Centrifuge_LLE):
#     _N_heat_utilities = 0
#     kwargs = PartitionFlash._kwargs
#     _run = PartitionFlash._run 
#     _setup = PartitionFlash._setup
#     def _set_phases(self):
#         top, bot = self.outs
#         top.phase = 'l'
#         bot.phase = 'L'

class RatioCentrifuge_LLE(Centrifuge_LLE):
    _N_heat_utilities = 0
    __init__ = RatioFlash.__init__
    _run = RatioFlash._run


class SplitCentrifuge_LLE(Centrifuge_LLE):
    __init__ = Splitter.__init__
    _run = Splitter._run
    split = Splitter.split
    # def __init__(self, ID='', ins=None, outs=(), *, order=None, split):
    #     Splitter.__init__(self, ID, ins, outs, split=split, order=order)
    
    
    