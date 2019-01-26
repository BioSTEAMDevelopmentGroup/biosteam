# -*- coding: utf-8 -*-
"""
Created on Wed Jan 16 17:49:26 2019

@author: yoelr
"""
from biosteam import Quantity

# def get_kg(dog):
#     return dog.weight.to('kg')    

# def set_kg(dog, weight_kg):
#     dog.weight = Quantity(weight_kg, 'kg')



class Dog:
    
    country = 'USA'
    nickname = 'doggy'
    
    def __init__(self, name, weight:'lbs'):
        self.name = name
        self.weight = Quantity(weight, 'lbs')

    @property
    def weight_kg(self):
        return self.weight.to('kg')  
    
    @weight_kg.setter
    def weight_kg(self, weight_kg):
        self.weight = Quantity(weight_kg, 'kg')

    def get_weight(self, units):
        return self.weight.to(units)
    
    def __repr__(self):
        return f'<{self.nickname}-{self.name}: weight={self.weight}>'
    
