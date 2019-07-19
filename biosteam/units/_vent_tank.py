# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 09:20:57 2019

@author: yoelr
"""
from .. import Unit

class VentTank:
        
    def _run(self):
        feed = self.ins[0]
        vent, effluent = self.outs
        effluent.copylike(feed)
        P_vapor = effluent.P_vapor
        
        
        