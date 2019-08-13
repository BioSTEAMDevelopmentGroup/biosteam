# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 22:17:05 2018

@author: yoelr
"""
from ._splitter import Splitter
from .decorators import cost

__all__ = ('MolecularSieve',)


# @cost('Flow rate', 'Pressure filter drying (2)',
#       cost=405000, CE=521.9, S=22687, n=0.6, kW=1044)
# @cost('Flow rate', 'Pressure filter pressing',
#       cost=75200, CE=521.9, S=22687, n=0.6, kW=112)
@cost('Flow rate', 'Column', kW=151, BM=1.8,
      cost=2601000, CE=521.9, S=22687, n=0.6)
class MolecularSieve(Splitter):
    """Create an ethanol/water molecular sieve for bioethanol plants.
    The molecular sieve is modeled as a component wise separator. Costing
    is based on scaling by the 6/10ths rule from an NREL TEA report [1].
    
    **Parameters**
    
        **split:** [array_like] Componentwise split to the 0th output stream

    **ins**
    
        [0] Feed (gas)
        
    **outs**
    
        [0] Split stream (gas)
        
        [1] Remainder stream (gas)
    
    **References**
    
        [1] Process Design and Economics for Biochemical Conversion of
        Lignocellulosic Biomass to Ethanol Dilute-Acid Pretreatment and
        Enzymatic Hydrolysis of Corn Stover. D. Humbird, R. Davis, L.
        Tao, C. Kinchin, D. Hsu, and A. Aden (National Renewable Energy
        Laboratory Golden, Colorado). P. Schoen, J. Lukas, B. Olthof,
        M. Worley, D. Sexton, and D. Dudgeon (Harris Group Inc. Seattle,
        Washington and Atlanta, Georgia)
    
    **Examples**
    
        :doc:`MolecularSieve Example`
    
    """
    _N_heat_utilities = 2
    _units = {'Flow rate': 'kg/hr'}
    def __init__(self, ID='', ins=None, outs=(), *, order=None, split):
        Splitter.__init__(self, ID, ins, outs, order=order, split=split)
        self._outs[0]._phase = self._outs[1]._phase = 'g'

    def _design(self):
        self._Design['Flow rate'] = flow = self._outs[1].massnet
        T = self.ins[0].T
        self._heat_utilities[0](1429.65*flow, T)
        self._heat_utilities[1](-55.51*flow, T)

