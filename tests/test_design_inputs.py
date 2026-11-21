# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Tests for design inputs (:func:`biosteam.design_input` and
:class:`biosteam.update_design_input`).
"""
import inspect
import pytest
import biosteam as bst
from inspect import Parameter
from biosteam._unit import (
    DesignInput,
    NumericalDesignInput,
    CategoricalDesignInput,
    LogicalDesignInput,
)

# %% Helpers

def create_vessel_class():
    class Vessel(bst.Unit):
        _F_BM_default = {'Vessel': 2.0}

        #: Residence time [hr].
        tau = bst.design_input('Residence time', 'hr', lb=0, default=1)

        #: Fractional efficiency.
        efficiency = bst.design_input(lb=0, ub=1, default=0.9)

        #: Construction material.
        material = bst.design_input(
            purchase_item='Vessel',
            catalog_type='Material',
            catalog={'Carbon steel': 1.0, 'Stainless steel': 2.5},
            default='Carbon steel',
        )

        #: Orientation (with design factor).
        orientation = bst.design_input(
            purchase_item='Vessel',
            catalog_type='Design',
            catalog={'Vertical': 1.0, 'Horizontal': 1.2},
            default='Vertical',
        )

        #: Type with no cost factor.
        vessel_type = bst.design_input(catalog=('Open', 'Closed'), default='Closed')

        #: Whether the vessel is insulated.
        insulated = bst.design_input(default=False)

        #: Number of vessels.
        N_vessels: int = bst.design_input(lb=1, default=1)

        @bst.design_input('Volume', 'm3', lb=0)
        def volume(self):
            return self.tau * self.feed.F_vol

        @volume.setter
        def volume(self, volume):
            self.tau = volume / self.feed.F_vol

        def _run(self):
            self.outs[0].copy_like(self.ins[0])

        def _cost(self):
            self.baseline_purchase_costs['Vessel'] = 1000.

    return Vessel

VESSEL_DESIGN_INPUTS = [
    'tau', 'efficiency', 'material', 'orientation', 'vessel_type',
    'insulated', 'N_vessels', 'volume'
]

@pytest.fixture
def Vessel():
    bst.settings.set_thermo(['Water'], cache=True)
    return create_vessel_class()

@pytest.fixture
def vessel(Vessel):
    feed = bst.Stream(None, Water=100)
    return Vessel(None, feed)

# %% Definition of design inputs

def test_design_input_factory_types():
    assert type(bst.design_input('Residence time', 'hr')) is NumericalDesignInput
    assert type(bst.design_input(lb=0)) is NumericalDesignInput
    assert type(bst.design_input(ub=1)) is NumericalDesignInput
    assert type(bst.design_input(catalog=('A', 'B'))) is CategoricalDesignInput
    assert type(bst.design_input(default=True)) is LogicalDesignInput
    assert type(bst.design_input(default=False)) is LogicalDesignInput
    assert type(bst.design_input()) is DesignInput
    assert type(bst.design_input(default='abc')) is DesignInput
    assert type(bst.design_input(default=1)) is DesignInput

    # Units default to empty string for numerical design inputs
    assert bst.design_input(lb=0).units == ''

    # Catalog type aliases
    for alias, key in [('Material', 'F_M'), ('material', 'F_M'), ('F_M', 'F_M'),
                       ('Design', 'F_D'), ('design', 'F_D'), ('F_D', 'F_D'),
                       ('Module', 'F_BM'), ('module', 'F_BM'), ('F_BM', 'F_BM')]:
        di = bst.design_input(catalog={'A': 1}, purchase_item='Vessel', catalog_type=alias)
        assert di.catalog_type == key
        assert di.purchase_item == 'Vessel'

def test_invalid_design_input_definitions():
    # Cost factor configuration requires a catalog
    with pytest.raises(ValueError):
        bst.design_input(catalog_type='Material')
    with pytest.raises(ValueError):
        bst.design_input(purchase_item='Vessel')

    # Both purchase_item and catalog_type (or neither) must be given
    with pytest.raises(ValueError):
        bst.design_input(catalog={'A': 1}, purchase_item='Vessel')
    with pytest.raises(ValueError):
        bst.design_input(catalog={'A': 1}, catalog_type='Material')
    with pytest.raises(ValueError):
        bst.design_input(catalog={'A': 1}, purchase_item='Vessel', catalog_type='Other')

    # Numerical options are not compatible with catalogs
    with pytest.raises(ValueError):
        bst.design_input(catalog=('A', 'B'), units='kg')
    with pytest.raises(ValueError):
        bst.design_input(catalog=('A', 'B'), lb=0)
    with pytest.raises(ValueError):
        bst.design_input(catalog=('A', 'B'), ub=1)

    # Numerical options are not compatible with logical design inputs
    with pytest.raises(ValueError):
        bst.design_input(default=True, units='kg')
    with pytest.raises(ValueError):
        bst.design_input(default=True, lb=0)
    with pytest.raises(ValueError):
        bst.design_input(default=True, ub=1)

    # Bare decorator does not accept other options
    with pytest.raises(ValueError):
        bst.design_input(lambda self: None, units='kg')

    # Both getter and setter must be defined (or none of them)
    with pytest.raises(RuntimeError):
        class GetterOnly(bst.Unit):
            @bst.design_input('Volume', 'm3')
            def volume(self): return 1.
    with pytest.raises(RuntimeError):
        class SetterOnly(bst.Unit):
            volume = bst.design_input('Volume', 'm3')
            @volume.setter
            def volume(self, volume): pass

    # Defaults are validated on class creation
    with pytest.raises(AttributeError):
        class OutOfBounds(bst.Unit):
            efficiency = bst.design_input(lb=0, ub=1, default=2)
    with pytest.raises(AttributeError):
        class NotInCatalog(bst.Unit):
            material = bst.design_input(catalog=('A', 'B'), default='C')

def test_design_input_registration(Vessel):
    # Design inputs are registered in order of definition
    assert Vessel._design_inputs == VESSEL_DESIGN_INPUTS
    assert all([isinstance(getattr(Vessel, i), DesignInput) for i in Vessel._design_inputs])

    # Design inputs are inherited and can be extended in subclasses
    class AgitatedVessel(Vessel):
        agitator_power = bst.design_input('Agitator power', 'kW', lb=0)

    assert AgitatedVessel._design_inputs == [*VESSEL_DESIGN_INPUTS, 'agitator_power']
    assert Vessel._design_inputs == VESSEL_DESIGN_INPUTS # Parent unchanged

    # Overriding a design input with a regular attribute removes it
    class FixedTimeVessel(Vessel):
        tau = 2

    assert 'tau' not in FixedTimeVessel._design_inputs

    # Design inputs from mixins are collected through the MRO
    class Mixin:
        pressure = bst.design_input('Pressure', 'Pa', lb=0)

    class PressureVessel(Mixin, Vessel): pass

    assert set(PressureVessel._design_inputs) == {*VESSEL_DESIGN_INPUTS, 'pressure'}

def test_design_input_names_and_metadata(Vessel):
    assert Vessel.tau.__name__ == 'tau'
    assert Vessel.tau.name == 'Residence time'
    assert Vessel.efficiency.name == 'Efficiency' # Default name from attribute
    assert Vessel.vessel_type.name == 'Vessel type'
    assert Vessel.N_vessels.name == 'N vessels'
    assert Vessel.volume.__name__ == 'volume'
    assert Vessel.volume.units == 'm3'
    assert Vessel.efficiency.units == ''
    assert (Vessel.efficiency.lb, Vessel.efficiency.ub) == (0, 1)
    assert Vessel.material.catalog_type == 'F_M'
    assert Vessel.orientation.catalog_type == 'F_D'
    assert Vessel.vessel_type.purchase_item is None
    assert Vessel.vessel_type.catalog_type is None

def test_design_input_docstrings(Vessel):
    assert Vessel.tau.__doc__ == 'Residence time [hr]. Must be positive.'
    assert Vessel.efficiency.__doc__ == 'Efficiency. Must be between 0 and 1.'
    assert Vessel.N_vessels.__doc__ == 'N vessels. Must be greater than 1.'
    assert Vessel.material.__doc__ == (
        "Material. Valid choices include 'Carbon steel' and 'Stainless steel'. "
        "The material will impact the material cost factor of 'Vessel'. "
        "Defaults to 'Carbon steel'."
    )
    assert Vessel.vessel_type.__doc__ == (
        "Vessel type. Valid choices include 'Open' and 'Closed'. "
        "Defaults to 'Closed'."
    )
    assert Vessel.insulated.__doc__ == 'Insulated. Defaults to False.'
    assert Vessel.volume.__doc__ == 'Volume [m3]. Must be positive.' # No docstring in getter

    class DocumentedVessel(bst.Unit):
        T_max = bst.design_input('Maximum temperature', 'K', ub=500)
        T_min = bst.design_input('Minimum temperature', 'K', ub=0)
        a = bst.design_input(catalog=('A',))
        abc = bst.design_input(
            catalog={'A': 1, 'B': 2, 'C': 3},
            purchase_item='Vessel', catalog_type='Module'
        )

        @bst.design_input('Volume', 'm3', lb=0)
        def volume(self) -> float:
            """Total volume."""
            return self._volume

        @volume.setter
        def volume(self, volume):
            self._volume = volume

        tau = bst.design_input('Residence time', 'hr', doc='Time in vessel.')

    assert DocumentedVessel.T_max.__doc__ == 'Maximum temperature [K]. Must be lower than 500.'
    assert DocumentedVessel.T_min.__doc__ == 'Minimum temperature [K]. Must be negative.'
    assert DocumentedVessel.a.__doc__ == "A. 'A' is the only valid choice."
    assert DocumentedVessel.abc.__doc__ == (
        "Abc. Valid choices include 'A', 'B', and 'C'. "
        "The abc will impact the bare module cost factor of 'Vessel'."
    )
    assert DocumentedVessel.volume.__doc__ == 'Total volume.' # Docstring from getter
    assert DocumentedVessel.tau.__doc__ == 'Time in vessel.' # Explicit docstring
    assert DocumentedVessel.volume.lb == 0 # Setter keeps configuration

    # Docs for parameters (used in Sphinx) include units and defaults
    assert Vessel.tau._make_doc(param=True) == 'Residence time [hr]. Defaults to 1.0 [hr].'
    assert Vessel.efficiency._make_doc(param=True) == 'Efficiency. Defaults to 0.9.'
    assert Vessel.material._make_doc(param=True) == "Material. Defaults to 'Carbon steel'."

def test_design_input_representation(Vessel):
    assert repr(Vessel.tau) == "NumericalDesignInput('Residence time', 'hr', default=1.0, lb=0)"
    assert repr(Vessel.efficiency) == "NumericalDesignInput('Efficiency', default=0.9, lb=0, ub=1)"
    assert repr(Vessel.insulated) == "LogicalDesignInput('Insulated', default=False)"
    assert repr(Vessel.material) == (
        "CategoricalDesignInput('Material', default='Carbon steel', "
        "catalog={'Carbon steel': 1.0, 'Stainless steel': 2.5}, "
        "purchase_item='Vessel', catalog_type='F_M')"
    )

    class Unit(bst.Unit):
        orientation = bst.design_input(catalog={'Vertical', 'Horizontal'})
        label = bst.design_input(default='abc')
        anything = bst.design_input()

    # Sets are sorted for a deterministic representation
    assert repr(Unit.orientation) == (
        "CategoricalDesignInput('Orientation', catalog={'Horizontal', 'Vertical'})"
    )
    assert repr(Unit.label) == "DesignInput('Label', default='abc')"
    assert repr(Unit.anything) == "DesignInput('Anything')"

def test_design_input_return_types():
    class Unit(bst.Unit):
        T: float = bst.design_input('Temperature', 'K', lb=0)
        N: int = bst.design_input(lb=1)
        flag = bst.design_input(default=True)
        option = bst.design_input(catalog=('A', 'B'))
        label = bst.design_input(default='abc')
        anything = bst.design_input()

        @bst.design_input('Volume', 'm3', lb=0)
        def volume(self) -> int:
            return self._volume

        @volume.setter
        def volume(self, volume):
            self._volume = volume

        @bst.design_input('Area', 'm2', lb=0)
        def area(self):
            return self._area

        @area.setter
        def area(self, area):
            self._area = area

    assert Unit.T._return_type is float
    assert Unit.N._return_type is int # From class annotation
    assert Unit.flag._return_type is bool
    assert Unit.option._return_type is str
    assert Unit.label._return_type is str # From default
    assert Unit.anything._return_type is None # No validation
    assert Unit.volume._return_type is int # From getter annotation
    assert Unit.area._return_type is float # Default for numerical inputs
    assert inspect.signature(Unit.area.fget).return_annotation is float

def test_string_annotations_are_resolved():
    namespace = {}
    exec(
        "from __future__ import annotations\n"
        "import biosteam as bst\n"
        "class Unit(bst.Unit):\n"
        "    N: int = bst.design_input(lb=1, default=2)\n"
        "    @bst.design_input(units='m3', lb=0)\n"
        "    def volume(self) -> int: return self._volume\n"
        "    @volume.setter\n"
        "    def volume(self, volume): self._volume = volume\n",
        namespace
    )
    Unit = namespace['Unit']
    assert Unit.N._return_type is int
    assert Unit.volume._return_type is int
    assert Unit.N.default == 2

def test_bare_design_input_decorator():
    bst.settings.set_thermo(['Water'], cache=True)

    class Unit(bst.Unit):
        @bst.design_input
        def label(self) -> str:
            """Label of unit."""
            return self._label

        @label.setter
        def label(self, label):
            self._label = label

    assert type(Unit.label) is DesignInput
    assert Unit.label.name == 'Label'
    assert Unit.label.__doc__ == 'Label of unit.'
    assert Unit._design_inputs == ['label']
    unit = Unit(None)
    unit.label = 'abc'
    assert unit.label == 'abc'
    with pytest.raises(AttributeError):
        unit.label = 1

# %% Setting and getting design inputs

def test_numerical_design_inputs(vessel):
    assert vessel.tau == 1
    vessel.tau = 2
    assert vessel.tau == 2
    assert type(vessel.tau) is float # Value is converted

    # Bounds are inclusive and enforced
    vessel.efficiency = 1
    vessel.efficiency = 0
    for value in (1.1, -0.1):
        with pytest.raises(AttributeError, match='efficiency'):
            vessel.efficiency = value
    with pytest.raises(AttributeError, match='cannot be below 0 hr'):
        vessel.tau = -1
    assert vessel.tau == 2 # Value unchanged after failed set

    # Type conversion must be lossless
    vessel.N_vessels = 2.0
    assert vessel.N_vessels == 2 and type(vessel.N_vessels) is int
    with pytest.raises(AttributeError, match='N_vessels must be a int'):
        vessel.N_vessels = 2.5
    with pytest.raises(AttributeError, match='N_vessels'):
        vessel.N_vessels = 0 # Below lower bound
    with pytest.raises(AttributeError, match='tau must be a float'):
        vessel.tau = 'long'
    with pytest.raises(AttributeError, match='tau must be a float'):
        vessel.tau = [1]
    assert vessel.tau == 2

    # None means the decision has not been made (no validation)
    vessel.efficiency = None
    assert vessel.efficiency is None

def test_design_inputs_with_getter_and_setter(vessel):
    feed = vessel.feed
    assert vessel.volume == vessel.tau * feed.F_vol
    vessel.volume = 4 * feed.F_vol
    assert vessel.tau == 4
    with pytest.raises(AttributeError):
        vessel.volume = -1 # Bounds are enforced before calling the setter
    assert vessel.tau == 4

def test_categorical_design_inputs(vessel):
    # Cost factors are set from defaults
    assert vessel.material == 'Carbon steel'
    assert vessel.orientation == 'Vertical'
    assert vessel.F_M['Vessel'] == 1.0
    assert vessel.F_D['Vessel'] == 1.0

    # Cost factors are updated with the design input
    vessel.material = 'Stainless steel'
    vessel.orientation = 'horizontal' # Capitalization is handled
    assert vessel.material == 'Stainless steel'
    assert vessel.orientation == 'Horizontal'
    assert vessel.F_M['Vessel'] == 2.5
    assert vessel.F_D['Vessel'] == 1.2

    # Invalid choices are rejected
    with pytest.raises(AttributeError, match="'Unobtanium' is not a recognised material"):
        vessel.material = 'Unobtanium'
    with pytest.raises(AttributeError):
        vessel.material = 1
    assert vessel.material == 'Stainless steel'
    assert vessel.F_M['Vessel'] == 2.5

    # Catalogs without cost factors do not affect costs
    vessel.vessel_type = 'Open'
    assert vessel.vessel_type == 'Open'
    with pytest.raises(AttributeError):
        vessel.vessel_type = 'Ajar'
    assert 'Vessel' not in vessel.F_BM or vessel.F_BM['Vessel'] == 2.0

    # Setting to None removes the cost factor
    vessel.material = None
    assert vessel.material is None
    assert 'Vessel' not in vessel.F_M

    # Cost factors from design inputs are used in costing
    vessel.material = 'Stainless steel'
    vessel.simulate()
    assert vessel.purchase_costs['Vessel'] == 1000. * 2.5 * 1.2
    # Installed cost = purchase cost * (F_BM + F_D * F_P * F_M - 1)
    assert vessel.installed_costs['Vessel'] == 1000. * (2.0 + 1.2 * 2.5 - 1)

def test_categorical_design_input_bare_module_factor():
    bst.settings.set_thermo(['Water'], cache=True)

    class Pump(bst.Unit):
        _F_BM_default = {'Pump': 3.0}
        seal = bst.design_input(
            catalog={'Packed': 1.0, 'Mechanical': 1.5},
            purchase_item='Pump', catalog_type='Module',
        )
        def _run(self): self.outs[0].copy_like(self.ins[0])
        def _cost(self): self.baseline_purchase_costs['Pump'] = 100.

    unit = Pump(None, bst.Stream(None, Water=1))
    assert unit.F_BM['Pump'] == 3.0
    unit.seal = 'Mechanical'
    assert unit.F_BM['Pump'] == 1.5
    unit.simulate()
    assert unit.installed_costs['Pump'] == 150.

def test_logical_and_generic_design_inputs(vessel):
    assert vessel.insulated is False
    vessel.insulated = True
    assert vessel.insulated is True
    for value in (1, 0, 'yes'):
        with pytest.raises(AttributeError, match='insulated must be either True or False'):
            vessel.insulated = value
    assert vessel.insulated is True

    class Unit(bst.Unit):
        label = bst.design_input(default='abc')
        anything = bst.design_input()

    unit = Unit(None)
    assert unit.label == 'abc'
    unit.label = 'xyz'
    with pytest.raises(AttributeError, match='label must be a str; int object is invalid'):
        unit.label = 1
    assert unit.label == 'xyz'

    # No return type means no validation
    for value in (1, 'a', [1, 2]):
        unit.anything = value
        assert unit.anything == value

def test_design_input_deletion(vessel):
    # Design inputs cannot be deleted unless a deleter is defined
    with pytest.raises(AttributeError, match='cannot delete'):
        del vessel.tau
    with pytest.raises(AttributeError, match='cannot delete'):
        del vessel.material

    class Unit(bst.Unit):
        @bst.design_input(units='m3', lb=0)
        def volume(self): return self._volume

        @volume.setter
        def volume(self, volume): self._volume = volume

        @volume.deleter
        def volume(self): del self._volume

    unit = Unit(None)
    unit.volume = 1
    del unit.volume
    assert not hasattr(unit, '_volume')

def test_design_input_unit_conversion(vessel):
    vessel.tau = 2
    assert vessel.get_design_input_units('tau') == 'hr'
    assert vessel.get_design_input('tau', 'min') == 120
    assert vessel.set_design_input('tau', 'day', 1) == 24
    assert vessel.tau == 24
    assert vessel.convert_design_input('tau', 'min', 60) == 1
    assert vessel.tau == 24 # Conversion does not set value
    with pytest.raises(AttributeError):
        vessel.set_design_input('tau', 'day', -1) # Bounds still enforced

# %% Initialization

def test_design_input_init_signature(Vessel):
    signature = inspect.signature(Vessel.__init__)
    parameters = signature.parameters

    # Design inputs are keyword-only arguments
    for name in VESSEL_DESIGN_INPUTS:
        assert parameters[name].kind is Parameter.KEYWORD_ONLY

    # Defaults are reflected in the signature
    assert parameters['tau'].default == 1
    assert parameters['efficiency'].default == 0.9
    assert parameters['material'].default == 'Carbon steel'
    assert parameters['insulated'].default is False
    assert parameters['volume'].default is Parameter.empty # No default

    # Unit parameters precede design inputs
    assert list(parameters)[:5] == ['self', 'ID', 'ins', 'outs', 'thermo']

def test_design_input_init_arguments(Vessel):
    feed = bst.Stream(None, Water=100)

    # Defaults are applied on initialization
    vessel = Vessel(None, feed)
    assert vessel.tau == 1
    assert vessel.efficiency == 0.9
    assert vessel.N_vessels == 1
    assert vessel.insulated is False

    # Design inputs can be passed as key-word arguments
    vessel = Vessel(
        None, feed, tau=3, efficiency=0.5, material='stainless steel',
        insulated=True, N_vessels=2,
    )
    assert vessel.tau == 3
    assert vessel.efficiency == 0.5
    assert vessel.material == 'Stainless steel'
    assert vessel.F_M['Vessel'] == 2.5
    assert vessel.insulated is True
    assert vessel.N_vessels == 2

    # Arguments are validated
    with pytest.raises(AttributeError):
        Vessel(None, feed, efficiency=2)
    with pytest.raises(AttributeError):
        Vessel(None, feed, material='Unobtanium')

def test_design_inputs_handled_by_init():
    bst.settings.set_thermo(['Water'], cache=True)

    class Unit(bst.Unit):
        a = bst.design_input(units='m', default=1.)
        b = bst.design_input(units='m')

        def _init(self, b, c=3):
            self.b = 2 * b
            self.c = c

    parameters = inspect.signature(Unit.__init__).parameters
    assert parameters['a'].default == 1.
    assert parameters['b'].default is Parameter.empty
    assert parameters['c'].default == 3

    # Defaults given by _init for design inputs show up in the signature
    class InitDefaultUnit(bst.Unit):
        a = bst.design_input(units='m')
        def _init(self, a=2.): self.a = a

    parameter = inspect.signature(InitDefaultUnit.__init__).parameters['a']
    assert parameter.default == 2.
    assert parameter.annotation is float
    assert InitDefaultUnit(None).a == 2.

    # Design inputs without return type are not annotated
    class UnannotatedUnit(bst.Unit):
        anything = bst.design_input()

    parameter = inspect.signature(UnannotatedUnit.__init__).parameters['anything']
    assert parameter.annotation is Parameter.empty
    assert Unit._init.parameters == frozenset(['b', 'c'])

    unit = Unit(None, b=1.)
    assert (unit.a, unit.b, unit.c) == (1., 2., 3)
    unit = Unit(None, a=4., b=2., c=5)
    assert (unit.a, unit.b, unit.c) == (4., 4., 5)

    # Defaults of design inputs handled by _init are passed to _init
    class DefaultedUnit(bst.Unit):
        a = bst.design_input(units='m', default=1.)

        def _init(self, a):
            self.a = 2 * a

    assert DefaultedUnit(None).a == 2.
    assert DefaultedUnit(None, a=5.).a == 10.

    # Defaults cannot be defined in both design input and _init
    with pytest.raises(RuntimeError):
        class ConflictingDefault(bst.Unit):
            a = bst.design_input(units='m', default=1.)
            def _init(self, a=2.): pass

# %% Updating inherited design inputs

def test_update_design_input(Vessel):
    feed = bst.Stream(None, Water=100)

    class LongVessel(Vessel):
        tau = bst.update_design_input(default=5., lb=1)
        material = bst.update_design_input(default='Stainless steel')
        efficiency = bst.update_design_input(name='Separation efficiency')

    # Updated copies are created; parent class is unchanged
    assert LongVessel.tau is not Vessel.tau
    assert repr(LongVessel.tau) == "NumericalDesignInput('Residence time', 'hr', default=5.0, lb=1)"
    assert repr(Vessel.tau) == "NumericalDesignInput('Residence time', 'hr', default=1.0, lb=0)"
    assert LongVessel.tau.__name__ == 'tau'
    assert LongVessel.tau.__doc__ == 'Residence time [hr]. Must be greater than 1.'
    assert LongVessel.efficiency.name == 'Separation efficiency'
    assert LongVessel.efficiency.default == 0.9
    assert LongVessel.material.catalog is Vessel.material.catalog
    assert LongVessel._design_inputs == Vessel._design_inputs

    vessel = LongVessel(None, feed)
    assert vessel.material == 'Stainless steel'
    assert vessel.F_M['Vessel'] == 2.5
    with pytest.raises(AttributeError):
        vessel.tau = 0.5 # Updated lower bound is enforced

    # Updates keep custom getters and setters
    class BigVessel(Vessel):
        volume = bst.update_design_input(lb=10)

    vessel = BigVessel(None, feed)
    vessel.volume = 20 * feed.F_vol
    assert vessel.tau == pytest.approx(20)
    with pytest.raises(AttributeError):
        vessel.volume = 1

    # Updates keep return types
    class ManyVessels(Vessel):
        N_vessels = bst.update_design_input(default=3)

    assert ManyVessels.N_vessels._return_type is int
    with pytest.raises(AttributeError):
        ManyVessels(None, feed).N_vessels = 2.5

    # Only inherited design inputs can be updated
    with pytest.raises(RuntimeError):
        class InvalidUpdate(Vessel):
            pressure = bst.update_design_input(default=101325)

    assert repr(bst.update_design_input(default=1)) == 'update_design_input(default=1)'

# %% Reporting design inputs

def test_design_inputs_in_results(vessel):
    vessel.tau = 2
    vessel.material = 'Stainless steel'
    vessel.simulate()
    df = vessel.results()
    for key in vessel._design_inputs:
        name = getattr(type(vessel), key).name
        assert ('Design inputs', name) in df.index
    series = vessel.results(with_units=False)
    assert series[('Design inputs', 'Residence time')] == 2
    assert series[('Design inputs', 'Material')] == 'Stainless steel'
    assert series[('Design inputs', 'Insulated')] == False
    assert df.loc[('Design inputs', 'Residence time'), 'Units'] == 'hr'

    # Undecided design inputs are not reported
    vessel.efficiency = None
    assert ('Design inputs', 'Efficiency') not in vessel.results().index

def test_unset_design_inputs_not_in_results():
    bst.settings.set_thermo(['Water'], cache=True)

    class Unit(bst.Unit):
        T = bst.design_input('Temperature', 'K', lb=0)
        P = bst.design_input('Pressure', 'Pa', lb=0, default=101325.)
        def _run(self): self.outs[0].copy_like(self.ins[0])

    unit = Unit(None, bst.Stream(None, Water=1))
    with pytest.raises(AttributeError):
        unit.T # Never set
    unit.simulate()
    index = unit.results().index
    assert ('Design inputs', 'Temperature') not in index
    assert ('Design inputs', 'Pressure') in index


if __name__ == '__main__':
    pytest.main([__file__])
