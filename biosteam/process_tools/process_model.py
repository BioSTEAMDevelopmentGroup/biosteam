# -*- coding: utf-8 -*-
"""
"""
from dataclasses import dataclass
from thermosteam.utils import AbstractMethod
from colorpalette import Color
import biosteam as bst

__all__ = ('ProcessModel', 'ScenarioComparison')

def copy(scenario, **kwargs):
    for i in scenario.__slots__:
        if i not in kwargs: kwargs[i] = getattr(scenario, i)
    return scenario.__class__(**kwargs)
    
def scenario_info(scenario, add_metadata):
    slots = scenario.__slots__
    if add_metadata: 
        grey = Color(fg='#878787')
        metadata = scenario.metadata
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
        if add_metadata and i in metadata:
            comment = '# ' + metadata[i]
            separate = len(arg + comment) > 50
            comment = grey(comment)
            if separate:
                arguments.append(comment)
                arguments.append(arg)
                continue
            else:
                arg += ' ' + comment
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
    
def display_scenario(scenario, metadata=True):
    print(scenario_info(scenario, metadata))
    
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
    """
    ProcessModel is an abstract class which has missing (i.e., "abstract") 
    attributes and methods. The user must define a `Scenario` dataclass which
    determines all inputs to the process model. Additionally, the user must 
    define the `as_scenario`, `create_thermo`, `create_system`, and 
    `create_model` methods for the process model to initialize its key components. 
    
    It may help to look at how ProcessModel objects are created (approximately):
    
    .. code-block:: python
    
        def __new__(cls, simulate=None, scenario=None, **kwargs):
            if scenario is None:
                self.scenario = cls.Scenario(**kwargs)
            else:
                # The Scenario object can be initialized through the `as_scenario` class method.
                self.scenario = cls.as_scenario(scenario)
            
            # No need to recreate a process model for repeated scenarios.
            if scenario in cls.cache: return cls.cache[scenario]
            self = super().__new__()
            
            # The thermodynamic property package is given by the `create_thermo` method.
            self.load_thermo(self.create_thermo())
            
            # If no system is returned by the `create_system` method, a new system is created from flowsheet units.
            self.flowsheet = bst.Flowsheet()
            system = self.create_system()
            if system is None: system = self.flowsheet.create_system()
            
            # This saves the system as self.system and all units/streams as attributes by ID.
            # For example, Stream('feedstock') will be stored as self.feestock.
            self.load_system(system) 
            
            # A Model object is loaded from the `create_model` method.
            # The model will be stored as self.model and all parameters and metrics as attributes by function name.
            # For example: 
            #
            # @model.metric
            # def MSP(): return self.tea.solve_price(self.product)
            #
            # ^ This becomes self.MSP.
            self.load_model(self.create_model())
            
            if simulate: self.system.simulate()
            self.cache[scenario] = self
            return self
    
    """
    #: **class-attribute** Class which defines arguments to the process model using
    #: the layout of a python dataclass: https://docs.python.org/3/library/dataclasses.html
    Scenario: type
    
    #: This method allows the process model to default the scenario.
    #: It should return a Scenario object.
    default_scenario = AbstractMethod
    
    #: **class method** This method allows the process model to interpret objects 
    #: (e.g., strings, numbers) as a Scenario.
    as_scenario = AbstractMethod
    
    #: This method should return a chemicals or thermo object.
    #: BioSTEAM will automatically set it as the thermodynmic property package.
    create_thermo = AbstractMethod
    
    #: This method should create unit operations and connect them.
    #: It can return a system object, optionally. Otherwise, BioSTEAM will 
    #: take care of creating the system from the units and saves 
    #: it as the self.system attribute.
    #: All streams and unit operations are also saved as attributes by their ID.
    create_system = AbstractMethod
    
    #: This method should return a model object. 
    #: The model will be saved as a self.model attribute. 
    #: All pareameters and metrics of the model object will also be saved as 
    #: attributes by their function names.
    create_model = AbstractMethod
    
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
        if not isinstance(scenario, cls.Scenario):
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
            Scenario.metadata = metadata = {}
            for i, j in tuple(Scenario.__dict__.items()):
                if i.startswith('__'): continue
                if isinstance(j, str):
                    if j[0] == '#':
                        delattr(Scenario, i)
                        metadata[i] = j.lstrip('# ')
                if isinstance(j, tuple) and len(j) == 2 and isinstance(j[-1], str):
                    value, j = j
                    if j[0] == '#':
                        setattr(Scenario, i, value)
                        metadata[i] = j.lstrip('# ')
                    
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
        thermo = self.create_thermo()
        if thermo is NotImplemented:
            raise NotImplementedError(f"{cls.__name__!r} is missing the `create_thermo` method")
        self.load_thermo(thermo)
        bst.main_flowsheet.set_flowsheet(self.flowsheet)
        unit_registry = self.flowsheet.unit
        unit_registry.open_context_level()
        system = self.create_system()
        if system is NotImplemented:
            raise NotImplementedError(f"{cls.__name__!r} is missing the `create_system` method")
        else:
            units = unit_registry.close_context_level()
            if system is None: system = bst.System.from_units(units=units)
            self.load_system(system)
        model = self.create_model()
        if model is NotImplemented:
            raise NotImplementedError(f"{cls.__name__!r} is missing the `create_model` method")
        elif model is None:
            raise RuntimeError('`create_model` must return a biosteam.Model object')
        self.load_model(model)
        if simulate: system.simulate()
        cls.cache[scenario] = self
        return self
    
    def load_thermo(self, thermo):
        bst.settings.set_thermo(thermo)
        thermo = bst.settings.get_thermo()
        self.chemicals = thermo.chemicals
        self.thermo = thermo
    
    def baseline(self):
        sample = self.model.get_baseline_sample()
        return sample, self.model(sample)
    
    def load_system(self, system):
        self.system = system
        self.flowsheet = flowsheet = system.flowsheet
        self.__dict__.update(flowsheet.to_dict())
            
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
    
    def show(self, metadata=True):
        """Print representation of process model."""
        scenario = self.scenario
        scenario_name = type(scenario).__name__
        process_name = type(self).__name__
        N = len(scenario_name)
        info = scenario_info(scenario, metadata)
        print(process_name + info[N:])
    
    _ipython_display_ = show
    