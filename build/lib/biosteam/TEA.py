# -*- coding: utf-8 -*-
"""
Created on Sat Aug 18 15:21:16 2018

@author: yoelr
"""


class TEA:
    def __init__(self, sys, project_length=15, IRR=0.10):
        self.sys = sys
        self.project_length = project_length
        self.IRR = IRR

    def operating_cost(self):
        pass

    def capital_cost(self):
        pass

    def NPV(self):
        op_cost = self.operating_cost()
        cap_cost = self.capital_cost()
        pass

    def TEA(self):
        self.run_results()
        return self.NPV()
