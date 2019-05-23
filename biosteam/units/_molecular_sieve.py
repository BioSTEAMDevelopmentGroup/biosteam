# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 22:17:05 2018

@author: yoelr
"""
from .metaclasses import splitter
from .. import Unit
from .decorators import cost, design

__all__ = ('MolecularSieve',)


@cost('Flow rate', 'Pressure filter drying (2)', N=2, cost=405000, CE=521.9, S=22687, exp=0.6, kW=522)
@cost('Flow rate', 'Pressure filter pressing', cost=75200, CE=521.9, S=22687, exp=0.6, kW=112)
@cost('Flow rate', 'Column', cost=2601000, CE=521.9, S=22687, exp=0.6)
@design('Flow rate', 'kg/hr', lambda self: self._massnet_in)
class MolecularSieve(Unit, metaclass=splitter):
    """Create an ethanol/water molecular sieve for bioethanol plants. The molecular sieve is modeled as a component wise separator. Costing is based on scaling by the 6/10ths rule from an NREL TEA report [1].
    
    **Parameters**
    
        **split:** [array_like] Componentwise split to the 0th output stream

    **ins**
    
        [0] Feed (gas)
        
    **outs**
    
        [0] Split stream (gas)
        
        [1] Remainder stream (gas)
    
    **References**
    
        [1] Process Design and Economics for Biochemical Conversion of Lignocellulosic Biomass to Ethanol Dilute-Acid Pretreatment and Enzymatic Hydrolysis of Corn Stover. D. Humbird, R. Davis, L. Tao, C. Kinchin, D. Hsu, and A. Aden (National Renewable Energy Laboratory Golden, Colorado). P. Schoen, J. Lukas, B. Olthof, M. Worley, D. Sexton, and D. Dudgeon (Harris Group Inc. Seattle, Washington and Atlanta, Georgia)
    
    **Examples**
    
        :doc:`MolecularSieve Example`
    
    """
    def _init(self):
        self._outs[0]._phase = self._outs[1]._phase = 'g'


