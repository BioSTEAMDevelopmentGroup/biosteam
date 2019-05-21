# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 21:23:56 2018

@author: yoelr
"""
from .. import Unit, MixedStream
from ._flash import RatioFlash
from .metaclasses import splitter
from .decorators import cost, design

__all__ = ('Centrifuge_LLE', 'RatioCentrifuge_LLE', 'SplitCentrifuge_LLE')

# electricity kW/(m3/hr) from USDA biosdiesel Super Pro model
# Possibly 1.4  kW/(m3/hr)
# https://www.sciencedirect.com/topics/engineering/disc-stack-centrifuge
# Microalgal fatty acids—From harvesting until extraction H.M. Amaro, et. al.,
# in Microalgae-Based Biofuels and Bioproducts, 2017

@cost('Flow rate', CE=525.4, cost=28100, exp=0.574, kW=3.66)
@design('Flow rate', 'm^3/hr', lambda self: self._volnet_out)
class Centrifuge_LLE(Unit):
    r"""Create an equlibrium based centrifuge with the option of having liquid non-keys and LIQUID non-keys completly separate into their respective phases.

    Equation for f.o.b cost:

    :math:`C_{fob}^{2007} = 28100 Q^{0.574} (0.1 < Q < 100 \frac{m^3}{h})` 

    **Parameters**

        **species_IDs:** *tuple[str]* IDs of equilibrium species
    
        **split:** *tuple[float]* Initial guess split fractions of equilibrium specise to the 'liquid' phase.
    
        **lNK:** *tuple[str]* Species assumed to completely remain in the 'liquid' phase.
    
        **LNK:** *tuple[str]* Species assumed to completely remain in the 'LIQUID' phase.
    
        **solvents:** *tuple[str]* Species corresponding to specified solvent_split.
    
        **solvent_split:** *tuple[float]* Split fractions of each specie to the 'liquid' phase.
         
    **ins**
    
        [0] Input stream
        
    **outs**
    
        [0] 'liquid' phase stream
        
        [1] 'LIQUID' phase stream
    
    **Examples**
    
        :doc:`Centrifuge_LLE Example`
    
    """
    line = 'Liquids Centrifuge'
    
    _kwargs = {'species_IDs': None,
               'split': None,
               'lNK': (),
               'LNK': (),
               'solvents': (),
               'solvent_split': ()}
    _bounds = {'Flow rate': (0.1, 100)}

    def _setup(self):
        self._mixedstream =  MixedStream(None)

    def _run(self):
        liq, LIQ = self.outs
        feed = self.ins[0]
        ms = self._mixedstream
        ms.empty()
        ms.liquid_mol[:] = feed.mol
        ms.LLE(T=feed.T,  **self._kwargs)
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
    _kwargs = {'Kspecies': [],  # list of species that correspond to Ks
               'Ks': [],  # list of molar ratio partition coefficinets,
               # Ks = y/x, where y and x are molar ratios of two different phases
               'top_solvents': [],  # list of species that correspond to top_split
               'top_split': [],  # list of splits for top_solvents
               'bot_solvents': [],  # list of species that correspond to bot_split
               'bot_split': []}  # list of splits for bot_solvents
    _setup = Unit._setup
    _run = RatioFlash._run
    
class SplitCentrifuge_LLE(Centrifuge_LLE, metaclass=splitter): pass
    
    
    