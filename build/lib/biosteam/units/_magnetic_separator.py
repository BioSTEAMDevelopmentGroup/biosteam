# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 11:32:01 2019

@author: yoelr
"""

from .. import Unit

class MagneticSeparator(Unit):
    _N_outs = 1
    _has_proxystream = True
    
    #: Original CEPCI
    CEPCI_0 = 576
    
    #: Original Price (USD)
    C_0 = 533571.0
    
    #: Original flow rate (kg/hr)
    V_0 = 333333
    
    #: Scaling exponent
    exp = 0.6
    
    def _run(self):
        self.outs[0].copylike(self.ins[0])
    
    def _cost(self):
        massflow = self.ins[0].massnet
        self._results['Cost']['Magnetic separator'] = self.CEPCI/self.CEPCI_0*self.C_0 * (massflow/self.V_0)**self.exp