# -*- coding: utf-8 -*-
"""
Created on Mon Feb  4 19:38:37 2019

@author: Guest Group
"""

from . import pd

DataFrame = pd.DataFrame

class TEA:
    
    options = DataFrame()
    
    def __init__(self, system):
        self.system = system
        
    def operating_cost(self):
        pass
    
    def capital_cost(self):
        pass
    
    def cash_flow(self):
        pass
    
    def NPV(self):
        pass
    
    def MFSP(self, stream):
        pass
    
    