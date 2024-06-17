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
    
    @property
    def parameters(self):
        return self.model._parameters
        
    @property
    def metrics(self):
        return self.model._metrics
            
    def _repr(self, m):
        clsname = type(self).__name__
        newline = "\n" + " "*(len(clsname)+2)
        return f'{clsname}: {newline.join([i.describe() for i in self.metrics])}'
    
    def __repr__(self):
        return f'<{type(self).__name__}: {len(self.parameters)}-parameters, {len(self.metrics)}-metrics>'
    
    def _info(self, p, m):
        return 'Process' + self.model._info(p, m)
    
    def show(self, p=None, m=None):
        """Return information on p-parameters and m-metrics."""
        print(self._info(p, m))
    _ipython_display_ = show
    