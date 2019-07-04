# -*- coding: utf-8 -*-
"""
Created on Sun Jun 30 21:21:45 2019

@author: yoelr
"""
from .. import Unit
from .. reaction import Rxn, PRxn

class LignocellulosicPretreatment(Unit):
    
    
    def _init(self):
        self.conversion = PRxn([
            #   Reaction definition                 Reactant    Conversion
            Rxn('Glucan + H2O -> Glucose',          'Glucan',   0.0990),
            Rxn('Glucan + H2O -> Glucose Oligomer', 'Glucan',   0.0030),
            Rxn('Glucan -> HMF + 2 H2O',            'Glucan',   0.0030),
            Rxn('Sucrose -> HMF + Glucose + 2H2O',  'Sucrose',  0.0030),
            Rxn('Xylan + H2O -> Xylose',            'Xylan',    0.9000),
            Rxn('Xylan + H2O -> XyloseOligomer',    'Xylan',    0.0024),
            Rxn('Xylan -> Furfural + 2 H2O',        'Xylan',    0.0050),
            Rxn('Acetate -> AceticAcid',            'Acetate',  1.0000),
            Rxn('Lignin -> SolubleLignin',          'Lignin',   0.0050)])
        
    def _run(self):
        feed = self.ins[0]
        product = self.outs[0]
        prodmol = product.mol
        prodmol[:] = feed.mol
        prodmol += self.conversion(product.mol)
        
        