# -*- coding: utf-8 -*-
"""
"""
from dataclasses import dataclass
__all__ = ('ProcessModel', 'scenario')

def display_scenario(scenario):
    slots = scenario.__slots__
    units_of_measure = scenario.units_of_measure
    arguments = []
    for i in slots:
        j = getattr(scenario, i)
        if isinstance(j, (bool, str)):
            arg = f"{i}={j!r},"
        else:
            try:
                arg = f"{i}={j:.3g},"
            except:
                arg = f"{i}={j},"
        if i in units_of_measure:
            arg += ' # ' + units_of_measure[i]
        arguments.append(
            arg
        )
    if arguments:
        clsname = type(scenario).__name__
        N_spaces = len(clsname) +1
        if N_spaces > 4:
            arguments = '\n    '.join(arguments)
            print(
                f"{clsname}("
                f"\n    {arguments}\n"
                ")"
            
            )
        else:
            spaces = N_spaces * ' '
            arguments = f'\n{spaces}'.join(arguments)
            print(
                f"{clsname}({arguments})"
            
            )
    else:
        print(
            f"{type(scenario).__name__}()"
        )
    
def scenario(cls=None, **units_of_measure):
    if cls is None: return lambda cls: scenario(cls, **units_of_measure)
    cls = dataclass(cls,
        init=True, repr=True, eq=True, unsafe_hash=False, frozen=True,
        match_args=True, slots=True,
    )
    cls.units_of_measure = units_of_measure
    cls._ipython_display_ = cls.show = display_scenario
    return cls

class ProcessModel:
    """Abstract class for setting up an all-in-one handle with easy 
    access to all objects related to a process model, including chemicals, 
    streams, units, systems, and model parameters and metrics."""
    
    def baseline(self):
        sample = self.model.get_baseline_sample()
        return sample, self.model(sample)
    
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
    