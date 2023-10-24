# -*- coding: utf-8 -*-
"""
"""

__all__ = ('ProcessModel',)

class ProcessModel:
    """Abstract class for setting up an all-in-one handle with easy 
    access to all objects related to a process model, including chemicals, 
    streams, units, systems, and model parameters and metrics."""
    
    def load_system(self, system):
        self.system = system
        self.flowsheet = flowsheet = system.flowsheet
        self.__dict__.update(flowsheet.to_dict())
        for i in system.units:
            self.chemicals = i.chemicals
            break
            
    def load_model(self, model):
        self.model = model
        for i in model.parameters:
            setattr(self, i.setter.__name__, i)
        for i in model.metrics:
            setattr(self, i.getter.__name__, i)