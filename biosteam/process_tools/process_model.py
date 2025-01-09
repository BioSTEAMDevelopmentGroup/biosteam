# -*- coding: utf-8 -*-
"""
"""
from dataclasses import dataclass
from thermosteam.utils import AbstractMethod
import biosteam as bst

__all__ = ('ProcessModel', 'ScenarioComparison')

def copy(scenario, **kwargs):
    for i in scenario.__slots__:
        if i not in kwargs: kwargs[i] = getattr(scenario, i)
    return scenario.__class__(**kwargs)
    
def scenario_info(scenario):
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
            return (
                f"{clsname}("
                f"\n    {arguments}\n"
                ")"
            
            )
        else:
            spaces = N_spaces * ' '
            arguments = f'\n{spaces}'.join(arguments)
            return (
                f"{clsname}({arguments})"
            
            )
    else:
        return (
            f"{type(scenario).__name__}()"
        )
    
def display_scenario(scenario):
    print(scenario_info(scenario))
    
def iterate_scenario_data(scenario):
    for i in scenario.__slots__:
        yield (i, getattr(scenario, i))

def scenario_comparison(left, right):
    return ScenarioComparison(left, right)

@dataclass(
    init=True, repr=True, eq=True, 
    unsafe_hash=False, frozen=True,
    match_args=True, slots=True,
)
class ScenarioComparison:
    left: object
    right: object

    def show(self):
        parameters = ',\n'.join(['left=' + scenario_info(self.left), 'right=' + scenario_info(self.right)])
        parameters = parameters.replace('\n', '\n    ')
        return print(
            f'{type(self).__name__}(\n'
            f'    {parameters}\n'
            f')'
        )
    _ipython_display_ = show
        

class ProcessModel:
    """Abstract class for setting up an all-in-one handle with easy 
    access to all objects related to a process model, including chemicals, 
    streams, units, systems, and model parameters and metrics."""
    default_scenario = AbstractMethod
    create_system = AbstractMethod
    create_model = AbstractMethod
    as_scenario = AbstractMethod
    
    @classmethod
    def scenario_hook(cls, scenario, kwargs):
        if scenario is None:
            if cls.default_scenario:
                scenario = cls.default_scenario()
            else:
                try:
                    return cls.Scenario(**kwargs)
                except:
                    raise NotImplementedError('missing class method `default_scenario`')
        elif not isinstance(scenario, cls.Scenario):
            if cls.as_scenario:
                scenario = cls.as_scenario(scenario)
            else:
                raise NotImplementedError('missing class method `as_scenario`')
        if kwargs: scenario = scenario.copy(**kwargs)
        return scenario
    
    def __init_subclass__(cls):
        cls.cache = {}
        if '__new__' in cls.__dict__: return
        if not hasattr(cls, 'Scenario'):
            raise NotImplementedError(
                "ProcessModel sublcass missing a 'Scenario' attribute"
            )
        if 'Scenario' in cls.__dict__:
            Scenario = cls.Scenario
            Scenario.units_of_measure = units_of_measure = {}
            for i, j in tuple(Scenario.__dict__.items()):
                if i.startswith('__'): continue
                if isinstance(j, str):
                    if j[0] == '[' and j[-1] == ']':
                        delattr(Scenario, i)
                        units_of_measure[i] = j
                if isinstance(j, tuple) and len(j) == 2 and isinstance(j[-1], str):
                    value, units = j
                    if units[0] == '[' and units[-1] == ']':
                        setattr(Scenario, i, value)
                        units_of_measure[i] = units
                    
            # TODO: get defaults and separate units of measure
            cls.Scenario = Scenario = dataclass(cls.Scenario,
                init=True, repr=True, eq=True, unsafe_hash=False, frozen=True,
                match_args=True, slots=True,
            )
            Scenario.items = iterate_scenario_data
            Scenario._ipython_display_ = Scenario.show = display_scenario
            Scenario.copy = copy
            Scenario.__sub__ = scenario_comparison
            
    
    def __new__(cls, *, simulate=True, scenario=None, **kwargs):
        scenario = cls.scenario_hook(scenario, kwargs)
        if scenario in cls.cache: return cls.cache[scenario]
        self = super().__new__(cls)
        self.scenario = scenario
        self.flowsheet = bst.Flowsheet(repr(self))
        bst.main_flowsheet.set_flowsheet(self.flowsheet)
        system = self.create_system()
        self.load_system(system)
        model = self.create_model()
        self.load_model(model)
        if simulate: system.simulate()
        cls.cache[scenario] = self
        return self
    
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
            if i.baseline is not None: i.setter(i.baseline)
        for i in model.metrics:
            setattr(self, i.getter.__name__, i)
            
    @property
    def parameters(self):
        return self.model._parameters
        
    @property
    def metrics(self):
        return self.model._metrics
            
    def __repr__(self):
        scenario = self.scenario
        scenario_name = type(scenario).__name__
        process_name = type(self).__name__
        N = len(scenario_name)
        return process_name + repr(self.scenario)[N:]
    
    def show(self):
        """Print representation of process model."""
        scenario = self.scenario
        scenario_name = type(scenario).__name__
        process_name = type(self).__name__
        N = len(scenario_name)
        info = scenario_info(scenario)
        print(process_name + info[N:])
    
    _ipython_display_ = show
    