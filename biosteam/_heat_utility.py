# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2024, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from __future__ import annotations
from thermosteam.units_of_measure import (
    convert, DisplayUnits, UnitsOfMeasure, get_dimensionality,
    heat_utility_units_of_measure
)
from thermosteam.utils import unregistered, define_units_of_measure
from thermosteam import Thermo, Stream, ThermalCondition, settings
from .exceptions import DimensionError
from math import copysign
from collections import deque
from typing import Optional, TYPE_CHECKING, Iterable, Literal, Sequence
if TYPE_CHECKING: from biosteam import Unit

__all__ = ('HeatUtility', 'UtilityAgent')

# Costs from Table 17.1 in Warren D.  Seider et al.-Product and Process Design Principles Synthesis, Analysis and Evaluation-Wiley (2016)
# ^This table was made using data from Busche, 1995
# Entry temperature conditions of coolants taken from Table 12.1 in Warren, 2016

mol_basis_units = UnitsOfMeasure('kmol')
mass_basis_units = UnitsOfMeasure('kg')
energy_basis_units = UnitsOfMeasure('kJ')

# %% Utility agents

@unregistered
class UtilityAgent(Stream):
    """
    Create a UtilityAgent object that defines a utility option.
    
    Parameters
    ----------
    ID : 
        A unique identification. If ID is None, stream will not be registered.
        If no ID is given, stream will be registered with a unique ID.
    flow : 
        All flow rates corresponding to defined chemicals.
    phase : 
        'g' for gas, 'l' for liquid, and 's' for solid. Defaults to 'l'.
    T : 
        Temperature [K]. Defaults to 298.15.
    P : 
        Pressure [Pa]. Defaults to 101325.
    units : 
        Flow rate units of measure (only mass, molar, and
        volumetric flow rates are valid). Defaults to 'kmol/hr'.
    thermo : 
        Thermo object to initialize input and output streams. Defaults to
        :meth:`settings.thermo <thermosteam._settings.ProcessSettings.thermo>`.
    T_limit :
        Temperature limit of outlet utility streams [K]. If no limit is given,
        phase change is assumed. If utility agent heats up, `T_limit` is
        the maximum temperature. If utility agent cools down, `T_limit` is
        the minimum temperature. 
    heat_transfer_price :
        Price of transferred heat [USD/kJ]. Defaults to 1.
    regeneration_price :
        Price of regenerating the fluid for reuse [USD/kmol]. Defaults to 0.
    heat_transfer_efficiency :
        Fraction of heat transferred accounting for losses to the environment (must be between 0 to 1). Defaults to 1.
    isfuel :
        Whether to burn the agent as a isfuel for heat.
    dT :
        Minimum temperature change between inlet and outlet utility. A positive value
        prevents near infinite flows when utility agents use sensible heats.
    **chemical_flows : float
        ID - flow pairs.
        
    """
    __slots__ = ('T_limit', '_heat_transfer_price', 'utility_stream_dump',
                 '_regeneration_price', 'heat_transfer_efficiency', 'isfuel',
                 'dT')
    def __init__(self, 
                 ID: Optional[str]='',
                 phase: Optional[str]='l',
                 T: Optional[float]=298.15,
                 P: Optional[float]=101325.,
                 units: Optional[str]=None,
                 thermo: Optional[Thermo]=None, 
                 T_limit: Optional[float]=None,
                 heat_transfer_price: float=0.,
                 regeneration_price: float=0.,
                 heat_transfer_efficiency: float=1.,
                 isfuel: bool=False,
                 dT: Optional[float]=0,
                 **chemical_flows: float):
        self._thermal_condition = ThermalCondition(T, P)
        thermo = self._load_thermo(thermo)
        self._init_indexer(None, phase, thermo.chemicals, chemical_flows)
        if units is not None:
            name, factor = self._get_flow_name_and_factor(units)
            flow = getattr(self, name)
            flow[:] = self.mol / factor
        self.mol[:] /= self.mol.sum() # Total flow must be 1 kmol / hr
        self.mol.read_only = True # Flow rate cannot change anymore
        self._sink = self._source = None
        self.reset_cache()
        self._register(ID)
        self.T_limit = T_limit
        self.heat_transfer_price = heat_transfer_price
        self.regeneration_price = regeneration_price
        self.heat_transfer_efficiency = heat_transfer_efficiency
        self.isfuel = isfuel
        self.dT = dT
        self.utility_stream_dump = []
    
    def _get_property(self, name, flow=False, nophase=False, T=None, P=None):
        # Take advantage of how composition is always constant and temperatures
        # and pressures are not expected to vary much.
        property_cache = self._property_cache
        thermal_condition = self._thermal_condition
        imol = self._imol
        composition = imol.data
        if not T: T = thermal_condition._T
        if not P: P = thermal_condition._P
        if nophase:
            literal = (T, P)
        else:
            phase = imol._phase
            literal = (phase, T, P)
        key = (name, literal)
        if key in property_cache: return property_cache[key]
        calculate = getattr(self.mixture, name)
        if nophase:
            property_cache[name] = value = calculate(
                composition, T, P
            )
        else:
            property_cache[name] = value = calculate(
                phase, composition, T, P
            )
        if len(property_cache) > 100: property_cache.pop(property_cache.__iter__().__next__())
        return value
    
    @property
    def iscooling_agent(self) -> bool:
        """Whether the agent is a cooling agent."""
        T_limit = self.T_limit 
        return T_limit > self.T if T_limit else self.phase == 'l'
    
    @property
    def isheating_agent(self) -> bool:
        """Whether the agent is a heating agent."""
        T_limit = self.T_limit 
        return T_limit < self.T if T_limit else self.phase == 'g'
    
    @property
    def price(self):
        raise AttributeError(f"'{type(self).__name__}' object has no attribute 'price'")
    @price.setter
    def price(self, price):
        raise AttributeError(f"'{type(self).__name__}' object has no attribute 'price'")
    
    @property
    def cost(self):
        raise AttributeError(f"'{type(self).__name__}' object has no attribute 'cost'")
    @cost.setter
    def cost(self, cost):
        raise AttributeError(f"'{type(self).__name__}' object has no attribute 'cost'")
    
    def to_stream(self, ID: Optional[str]=None):
        """
        Return a copy as a :class:`~thermosteam.Stream` object.

        Examples
        --------
        >>> import biosteam as bst
        >>> bst.settings.set_thermo(['Water', 'Ethanol']) 
        >>> cooling_water = bst.HeatUtility.get_agent('cooling_water')
        >>> cooling_water_copy = cooling_water.to_stream('cooling_water_copy')
        >>> cooling_water_copy.show(flow='kg/hr')
        Stream: cooling_water_copy
        phase: 'l', T: 305.37 K, P: 101325 Pa
        flow (kg/hr): Water  18
        
        """
        new = Stream.__new__(Stream)
        new._sink = new._source = None
        new._thermo = self._thermo
        new._imol = self._imol.copy()
        new._thermal_condition = self._thermal_condition.copy()
        new.reset_cache()
        new._price = 0.
        new.characterization_factors = {}
        new._register(ID)
        return new
    
    @property
    def heat_transfer_price(self) -> float:
        """Price of transfered heat [USD/kJ]."""
        return self._heat_transfer_price
    @heat_transfer_price.setter
    def heat_transfer_price(self, price: float):
        assert price >= 0, "heat transfer price cannot be negative"
        self._heat_transfer_price = price
    
    @property
    def regeneration_price(self) -> float:
        """Price of regenerating the fluid for reuse [USD/kmol]."""
        return self._regeneration_price
    @regeneration_price.setter
    def regeneration_price(self, price: float):
        assert price >= 0, "regeneration price cannot be negative"
        self._regeneration_price = price
    
    def _info_phaseTP(self, phase, units, notation):
        T_notation = notation['T']
        P_notation = notation['P']
        T_units = units['T']
        P_units = units['P']
        if self.T_limit is None:
            T_limit = "None"
        else:
            T_limit = convert(self.T_limit, 'K', T_units)
            T_limit = f"{T_limit:.3g} K" if T_limit else "None"
        T = convert(self.T, 'K', T_units)
        P = convert(self.P, 'Pa', P_units)
        s = '' if isinstance(phase, str) else 's'
        ht_price = self.heat_transfer_price
        rg_price = self.regeneration_price
        ht_eff = self.heat_transfer_efficiency
        return (f"heat_transfer_efficiency: {ht_eff:.3f}\n"
                f"heat_transfer_price: {ht_price:.3g} USD/kJ\n"
                f"regeneration_price: {rg_price:.3g} USD/kmol\n"
                f"T_limit: {T_limit}\n"
                f"phase{s}: {repr(phase)}\n"
                f"T: {T:{T_notation}} {T_units}\n"
                f"P: {P:{P_notation}} {P_units}\n"
        )
    def __repr__(self):
        return f"<{type(self).__name__}: {self.ID}>"


# %%

@define_units_of_measure(heat_utility_units_of_measure)
class HeatUtility:
    """
    Create an HeatUtility object that can choose a utility stream and 
    calculate utility requirements. It can calculate required flow rate,
    temperature change, or phase change of utility. Calculations assume 
    counter current flow rate.
    
    Parameters
    ----------
    heat_transfer_efficiency : 
        Enforced fraction of heat transferred from utility (due
        to losses to environment).
    unit : 
        Parent unit using this heat utility.
    hxn_ok :
        Whether heat utility can be satisfied within a heat exchanger network.
    
    Examples
    --------
    Create a heat utility:
    
    >>> from biosteam import HeatUtility
    >>> hu = HeatUtility()
    >>> hu.show()
    HeatUtility: None
     duty: 0
     flow: 0
     cost: 0
    
    Calculate utility requirement by calling it with a duty (kJ/hr), and entrance and exit temperature (K):
     
    >>> hu(1000, 300, 350)
    >>> hu.show()
    HeatUtility: low_pressure_steam
     duty: 1.05e+03 kJ/hr
     flow: 0.0272 kmol/hr
     cost: 0.00647 USD/hr
   
    All results are accessible:
    
    >>> hu.ID, hu.duty, hu.flow, hu.cost
    ('low_pressure_steam', 1052.6315789473686, 0.02721274387089031, 0.006471190492497716)
   
    All results are accessible:
    
    >>> hu.ID, hu.duty, hu.flow, hu.cost
    ('low_pressure_steam', 1052.6315789473686, 0.02721274387089031, 0.006471190492497716)
           
    """
    __slots__ = (
        'inlet_utility_stream', 'outlet_utility_stream', 'agent', 
        'duty', 'flow', 'cost', 'heat_transfer_efficiency', 
        'unit', 'hxn_ok', 'unit_duty',
        'oxygen_rich_inlet', # Only for fuels
    )
    #: Minimum approach temperature difference. Used to assign
    #: the pinch temperature of the utility stream.
    dT: float = 5  
    
    #: Units of measure for IPython display.
    display_units: DisplayUnits = DisplayUnits(duty='kJ/hr', flow='kmol/hr', cost='USD/hr')

    #: Used broadly throughout utilities.
    thermo_water: Thermo = Thermo(['Water'])
    
    #: Used for refrigeration utility.
    thermo_propane: Thermo = Thermo(['Propane'])

    #: Used for refrigeration utility.
    thermo_propylene: Thermo = Thermo(['Propylene'])
    
    #: Used for refrigeration utility.
    thermo_ethylene: Thermo = Thermo(['Ethylene'])

    #: Used exclusively for furnaces.
    thermo_natural_gas: Thermo = Thermo(['Methane', 'N2', 'CO2', 'O2', 'H2O'])

    #: Characterization factor data (value and units) by agent ID and impact key.
    characterization_factors: dict[tuple[str, str], tuple[float, UnitsOfMeasure]] = {}
    
    #: All cooling utilities available.
    cooling_agents: list[UtilityAgent]
        
    #: All heating utilities available.
    heating_agents: list[UtilityAgent]

    @classmethod
    def set_CF(self, ID: str, key: str, value: float, basis: Optional[str]=None, units: Optional[str]=None):
        """
        Set the characterization factor of a utility agent for a given impact 
        key.

        Parameters
        ----------
        ID :
            ID of utility agent.
        key :
            Key of impact indicator.
        value :
            Characterization factor value.
        basis :
            Basis of characterization factor. Valid dimensions include weight, 
            molar, and energy (e.g. 'kg', 'kmol', 'kJ'). Defaults to 'kg'.
        units :
            Units of impact indicator. Before using this argument, the default units 
            of the impact indicator should be defined with 
            :meth:`settings.define_impact_indicator <thermosteam._settings.ProcessSettings.define_impact_indicator>`.
            Units must also be dimensionally consistent with the default units.
            
        Raises
        ------
        ValueError
            When the duty characterization factor for a cooling agent is positive (should be negative).
        DimensionError
            When characterization factor is not given in dimensions of
            weight, molar, or energy.

        See Also
        --------
        HeatUtility.get_CF
        
        Examples
        --------
        Set the GWP characterization factor for low pressure steam at 
        88.44 kg CO2e / mmBtu (GREET 2020; Steam Production via Small Boiler from North American Natural Gas):
        
        >>> import biosteam as bst
        >>> bst.HeatUtility.set_CF('low_pressure_steam', 'GWP [kg CO2e]', 88.44, basis='MMBtu')
        
        Retrieve the GWP characterization factor for low pressure steam on a
        Btu basis:
        
        >>> bst.HeatUtility.get_CF('low_pressure_steam', 'GWP [kg CO2e]', basis='Btu')
        8.844e-05

        """
        agent = self.get_agent(ID)
        if units is not None:
            original_units = settings.get_impact_indicator_units(key)
            value = original_units.unconvert(value, units)
        if basis is None:
            basis_units = mass_basis_units
        else:
            dim = get_dimensionality(basis)
            if dim == mol_basis_units.dimensionality: 
                basis_units = mol_basis_units
            elif dim == mass_basis_units.dimensionality:
                basis_units = mass_basis_units
            elif dim == energy_basis_units.dimensionality:
                basis_units = energy_basis_units
                if agent.iscooling_agent and value > 0.: 
                    raise ValueError(
                        'duty characterization factor must be negative for cooling agents'
                    )
            else:
                raise DimensionError(
                    "dimensions for characterization factors must be in a molar, "
                   f"mass or energy basis, not '{dim}'"
                )
            value *= basis_units.conversion_factor(basis)
        self.characterization_factors[agent.ID, key] = value, basis_units

    @classmethod
    def get_CF(self, ID: str, key: str, basis: Optional[str]=None, units: Optional[str]=None):
        """
        Return the characterization factor of a utility agent for a given impact 
        key.

        Parameters
        ----------
        ID :
            ID of utility agent.
        key :
            Key of impact indicator.
        basis :
            Basis of characterization factor. Valid dimensions include weight, 
            molar, and energy (e.g. 'kg', 'kmol', 'kJ'). Defaults to 'kg'.
        units :
            Units of impact indicator. Before using this argument, the default units 
            of the impact indicator should be defined with 
            :meth:`settings.define_impact_indicator <thermosteam._settings.ProcessSettings.define_impact_indicator>`.
            Units must also be dimensionally consistent with the default units.
            
        Raises
        ------
        DimensionalityError
            When the characterization factor cannot be converted to the given basis 
            due to inconsistent dimensions with the original basis.

        See Also
        --------
        HeatUtility.set_CF

        """
        try:
            value, basis_units = self.characterization_factors[ID, key]
        except KeyError:
            if basis is None:
                return 0., None
            else:
                return 0.
        if units is not None:
            original_units = settings.get_impact_indicator_units(key)
            value = original_units.convert(value, units)
        if basis is None:
            return value, basis_units.units
        else:
            return value / basis_units.conversion_factor(basis)

    def get_impact(self, key: str):
        agent = self.agent
        CF, basis_units = self.get_CF(agent.ID, key)
        if CF == 0.: return 0.
        if basis_units == 'kg':
            return self.flow * agent.MW * CF
        elif basis_units == 'mol':
            return self.flow * CF
        elif basis_units == 'kJ':
            return self.duty * CF
        else:
            raise RuntimeError("unknown error")

    def get_inventory(self, key: str):
        agent = self.agent
        CF, basis_units = self.get_CF(agent.ID, key)
        if basis_units == 'kg':
            return self.flow * agent.MW, basis_units + '/hr'
        elif basis_units == 'mol':
            return self.flow, basis_units + '/hr'
        elif basis_units == 'kJ':
            return self.duty, basis_units + '/hr'
        else:
            raise RuntimeError("unknown error")

    def __init__(self,
            heat_transfer_efficiency: Optional[float]=None, 
            unit: Optional[Unit]=None,
            hxn_ok: Optional[bool]=False,
        ):
        #: Enforced fraction of heat transferred from utility (due
        #: to losses to environment).
        self.heat_transfer_efficiency: float = heat_transfer_efficiency
        
        #: Parent unit using this heat utility.
        self.unit: Unit|None = unit
        
        #: Whether heat utility can be satisfied within a heat exchanger network.
        self.hxn_ok: bool = hxn_ok
        
        #: Total heat transferred from utility to both the process and the environment [kJ/hr].
        self.duty: float = 0.
        
        #: Effective heat transferred from utility to the unit operation [kJ/hr].
        self.unit_duty: float = 0.
        
        #: Flow rate of utility [kmol/hr].
        self.flow: float = 0.
        
        #: Cost of utility [USD/hr].
        self.cost: float = 0.
        
        #: Utility agent being used.
        self.agent: UtilityAgent|None = None
        
        #: Fresh utility stream
        self.inlet_utility_stream: Stream|None = None
        
        #: Used utility stream
        self.outlet_utility_stream: Stream|None = None

    def __bool__(self) -> bool:
        return bool(self.agent)

    @property
    def ID(self) -> str:
        """ID of utility agent being used."""
        agent = self.agent
        return agent.ID if agent else ""
    @ID.setter
    def ID(self, ID: str):
        """ID of utility agent being used."""
        self.agent = self.get_agent(ID)

    @classmethod
    def default_agents(cls):
        """Reset all agents back to BioSTEAM's defaults."""
        cls.default_heating_agents()
        cls.default_cooling_agents()

    @classmethod
    def default_heating_agents(cls):
        """Reset all heating agents back to BioSTEAM's defaults."""
        thermo_water = cls.thermo_water
        low_pressure_steam = UtilityAgent(
            'low_pressure_steam',
            Water=1, T=412.189, P=344738.0, phase='g',
            thermo=thermo_water,
            regeneration_price = 0.2378,
            heat_transfer_efficiency = 0.95,
        )
        medium_pressure_steam = UtilityAgent(
            'medium_pressure_steam',
            Water=1, T=454.770, P=1.041e+6, phase='g',
            thermo=thermo_water,
            regeneration_price = 0.2756,
            heat_transfer_efficiency = 0.90,
        ) 
        high_pressure_steam = UtilityAgent(
            'high_pressure_steam',
            Water=1, T=508.991, P=3.11e+6, phase='g', 
            thermo=thermo_water,
            regeneration_price = 0.3171,
            heat_transfer_efficiency = 0.85,
        ) 
        natural_gas = UtilityAgent(
            'natural_gas',
            Methane=1, T=298.15, P=200 * 101325, phase='g', 
            thermo=cls.thermo_natural_gas,
            heat_transfer_efficiency = 0.90, # Heat loss to environment
            regeneration_price=3.49672,
            T_limit=405, # Must be reasonably higher than the emission's dew point at 500 psig (401 K)
            isfuel=True,
        ) 
        cls.heating_agents = [low_pressure_steam, 
                              medium_pressure_steam, 
                              high_pressure_steam,
                              natural_gas]
        
    @classmethod
    def default_cooling_agents(cls):
        """Reset all cooling agents back to BioSTEAM's defaults."""
        thermo_water = cls.thermo_water
        cooling_water = UtilityAgent(
            'cooling_water',
            Water=1, T=305.372, P=101325,
            thermo=thermo_water,
            T_limit = 324.817,
            regeneration_price = 4.8785e-4,
            dT=2.,
        )
        chilled_water = UtilityAgent(
            'chilled_water',
            Water=1, T=280.372, P=101325,
            thermo=thermo_water,
            T_limit = 300.372,
            heat_transfer_price = 5e-6,
            dT=2.,
        )
        chilled_brine = UtilityAgent(
            'chilled_brine',
            Water=1, T=255.372, P=101325,
            thermo=thermo_water,
            T_limit = 275.372,
            heat_transfer_price = 8.145e-6,
            dT=2.,
        )
        propane = UtilityAgent(
            'propane',
            Propane=1,
            thermo=cls.thermo_propane,
            T=238.70,
            P=cls.thermo_propane.chemicals.Propane.Psat(238.70),
            heat_transfer_price = 13.17e-6,
            phase='l',
            dT=1.,
        )
        propylene = UtilityAgent(
            'propylene',
            Propylene=1,
            thermo=cls.thermo_propylene,
            T=227.59,
            P=cls.thermo_propylene.chemicals.Propylene.Psat(227.59),
            heat_transfer_price = 16.54e-6, # Lever rule with -30 and -90 F prices
            phase='l',
            dT=1.,
        )
        ethylene = UtilityAgent(
            'ethylene',
            Ethylene=1,
            thermo=cls.thermo_ethylene,
            T=172.04,
            P=cls.thermo_ethylene.chemicals.Ethylene.Psat(172.04),
            heat_transfer_price = 33.2e-06,
            phase='l',
            dT=1.,
        )
        cls.cooling_agents = [cooling_water, 
                              chilled_water, 
                              chilled_brine,
                              propane,
                              propylene,
                              ethylene]

    def copy(self) -> HeatUtility:
        hu = HeatUtility()
        hu.copy_like(self)
        return hu

    def copy_like(self, other: HeatUtility):
        """Copy all data from another heat utility."""
        if other.agent:
            self.load_agent(other.agent)
            self.inlet_utility_stream = other.inlet_utility_stream.copy()
            self.outlet_utility_stream = other.outlet_utility_stream.copy()
        self.flow = other.flow
        self.duty = other.duty
        self.unit_duty = other.unit_duty
        self.cost = other.cost
        self.heat_transfer_efficiency = other.heat_transfer_efficiency

    def scale(self, factor: float):
        """Scale utility data."""
        self.flow *= factor
        self.duty *= factor
        self.cost *= factor
        self.inlet_utility_stream.mol *= factor
        if self.agent.isfuel: 
            self.outlet_utility_stream.mol *= factor
            self.oxygen_rich_inlet.mol *= factor
        # No need to factor the outlet utility stream
        # because it shares the same flow rate data as the inlet

    def empty(self):
        """Remove utility requirements."""
        if self.agent:
            agent = self.agent
            dump = agent.utility_stream_dump
            if len(dump) < 1000:
                if agent.isfuel:
                    dump.append(
                        (self.inlet_utility_stream, self.outlet_utility_stream, self.oxygen_rich_inlet)
                    )
                else:
                    dump.append(
                        (self.inlet_utility_stream, self.outlet_utility_stream)
                    )
        self.cost = self.flow = self.duty = self.unit_duty = 0
        self.outlet_utility_stream = self.inlet_utility_stream = self.agent = None
        
    def set_utility_by_flow_rate(self, agent: Optional[UtilityAgent], F_mol: float):
        if F_mol == 0.: 
            self.empty()
            return
        self.load_agent(agent)
        heat_transfer_efficiency = self.heat_transfer_efficiency or agent.heat_transfer_efficiency
        if agent.T_limit or agent.isfuel: raise ValueError('agent must work by latent heat to set by flow rate')
        if self.inlet_utility_stream.phase == 'l' :
            self.outlet_utility_stream.phase = 'g'
            dh = -agent._get_property('Hvap', nophase=True)
        else:
            dh = agent._get_property('Hvap', nophase=True)
            self.outlet_utility_stream.phase = 'l'
        self.duty = duty = dh * F_mol
        self.unit_duty = duty * heat_transfer_efficiency
        self.outlet_utility_stream.mol[:] = F_mol
        self.flow = F_mol
        self.cost = agent._heat_transfer_price * abs(duty) + agent._regeneration_price * F_mol
        
    def __call__(self, 
            unit_duty: float, 
            T_in: float, 
            T_out: Optional[float]=None, 
            agent: Optional[UtilityAgent]=None):
        """
        Calculate utility requirements given the essential parameters.
        
        Parameters
        ----------
        unit_duty :
               Unit duty requirement [kJ/hr]
        T_in : 
               Inlet process stream temperature [K]
        T_out : 
               Outlet process stream temperature [K]
        agent : 
                Utility agent to use. Defaults to a suitable agent from 
                predefined heating/cooling utility agents.
        
        """
        if unit_duty == 0:
            self.empty()
            return
        T_out = T_out or T_in
        iscooling = unit_duty < 0. #: Whether the utility is cooling the process.
        
        # Note: These are pinch temperatures at the utility inlet and outlet. 
        # Not to be confused with the inlet and outlet of the process stream.
        # If cooling, T_pinch_in is the minimum utility temperature required.
        # If heating, T_pinch_in is the maximum utility temperature required.
        T_pinch_in, T_pinch_out = self.get_inlet_and_outlet_pinch_temperature(
            iscooling, T_in, T_out
        )
        ## Select heat transfer agent ##
        if agent:
            if agent is not self.agent: self.load_agent(agent)
        else:
            agent = (self.get_suitable_cooling_agent if iscooling else self.get_suitable_heating_agent)(T_pinch_in)
            self.load_agent(agent)
            
        ## Calculate utility requirement ##
        heat_transfer_efficiency = self.heat_transfer_efficiency or agent.heat_transfer_efficiency
        duty = unit_duty / heat_transfer_efficiency
        if agent.isfuel:
            T_emissions = self.get_outlet_temperature(
                T_pinch_out, agent.T_limit, iscooling
            )
            feed = self.inlet_utility_stream
            emissions = self.outlet_utility_stream
            emissions.copy_like(feed)
            emissions.P = 101325
            reactions = agent.chemicals.get_combustion_reactions()
            reactions.force_reaction(emissions)
            O2_consumption = -emissions.imol['O2']
            oxygen_rich_gas = self.oxygen_rich_inlet
            z_O2 = oxygen_rich_gas.imol['O2'] / oxygen_rich_gas.F_mol
            oxygen_rich_gas.F_mol = O2_consumption / z_O2
            emissions.mol += oxygen_rich_gas.mol
            F_emissions = emissions.F_mass
            z_CO2 = emissions.imass['CO2'] / F_emissions
            z_CO2_target = 0.055 # Usually between 4 - 7 for biomass and natural gas (https://www.sciencedirect.com/science/article/pii/S0957582021005127)
            F_emissions_new = z_CO2 * F_emissions / z_CO2_target
            dF_emissions = F_emissions_new - F_emissions
            oxygen_rich_gas.F_mass = F_mass_O2_new = oxygen_rich_gas.F_mass + dF_emissions
            emissions.mol += oxygen_rich_gas.mol * (dF_emissions / F_mass_O2_new)
            emissions.T = T_emissions
            emissions.P = 3548325.0 # 500 psig
            dh = feed.Hnet - emissions.Hnet
        elif agent.T_limit:
            # Temperature change
            self.outlet_utility_stream.T = T_outlet = self.get_outlet_temperature(
                T_pinch_out, agent.T_limit, iscooling
            )
            dh = agent._get_property('H') - agent._get_property('H', T=T_outlet)
        else:
            # Phase change
            if agent.phase == 'l':
                self.outlet_utility_stream.phase = 'g'
                dh = -agent._get_property('Hvap', nophase=True)
            else:
                dh = agent._get_property('Hvap', nophase=True)
                self.outlet_utility_stream.phase = 'l'
        
        # Update utility flow
        F_mol = duty / dh
        self.inlet_utility_stream.mol[:] *= F_mol
        if agent.isfuel: 
            emissions.mol[:] *= F_mol
            oxygen_rich_gas.mol[:] *= F_mol
        
        # Update results
        self.unit_duty = unit_duty
        self.flow = F_mol
        self.duty = duty
        self.cost = agent._heat_transfer_price * abs(duty) + agent._regeneration_price * F_mol

    @property
    def inlet_process_stream(self) -> Stream:
        """If a heat exchanger is available, this stream is the inlet 
        process stream to the heat exchanger."""
        heat_exchanger = self.unit
        if heat_exchanger:
            return heat_exchanger.inlet
        else:
            raise AttributeError('no heat exchanger available '
                                 'to retrieve process stream')

    @property
    def outlet_process_stream(self) -> Stream:
        """If a heat exchanger is available,
        this stream is the outlet process stream to the heat exchanger."""
        heat_exchanger = self.unit
        if heat_exchanger:
            return heat_exchanger.outlet
        else:
            raise AttributeError('no heat exchanger available '
                                 'to retrieve process stream')
    
    @staticmethod
    def heat_utilities_by_agent(heat_utilities: Iterable[HeatUtility]):
        """Return a dictionary of heat utilities sorted by agent ID."""
        heat_utilities = [i for i in heat_utilities if i.agent]
        heat_utilities_by_agent = {i.ID: [] for i in heat_utilities}
        for i in heat_utilities:
            heat_utilities_by_agent[i.agent.ID].append(i)
        return heat_utilities_by_agent

    @classmethod
    def sum(cls, heat_utilities: Iterable[HeatUtility]):
        """Return a HeatUtility object that reflects the sum of heat
        utilities."""
        heat_utility = cls()
        heat_utility.mix_from(heat_utilities)
        return heat_utility    

    @classmethod
    def sum_by_agent(cls, heat_utilities: Iterable[HeatUtility]):
        """Return a list of heat utilities that reflect the sum of heat utilities
        by agent."""
        heat_utilities_by_agent = cls.heat_utilities_by_agent(heat_utilities)
        return [cls.sum(i) for i in heat_utilities_by_agent.values()]

    @classmethod
    def get_agent(cls, ID: str):
        """Return utility agent with given ID."""
        for agent in cls.heating_agents + cls.cooling_agents:
            if agent.ID == ID: return agent
        raise LookupError(ID)

    @classmethod
    def get_heating_agent(cls, ID: str):
        """Return heating agent with given ID."""
        for agent in cls.heating_agents:
            if agent.ID == ID: return agent
        raise LookupError(ID)
  
    @classmethod
    def get_cooling_agent(cls, ID: str):
        """Return cooling agent with given ID."""
        for agent in cls.cooling_agents:
            if agent.ID == ID: return agent
        raise LookupError(ID)
    
    @classmethod
    def get_suitable_heating_agent(cls, T_pinch: float):
        """
        Return a heating agent that works at the pinch temperature.
        
        Parameters
        ----------
        T_pinch :
            Pinch temperature [K].
        
        """
        for agent in cls.heating_agents:
            if T_pinch < agent.T - agent.dT or agent.isfuel: return agent
        raise RuntimeError(f'no heating agent that can heat over {T_pinch} K')    

    @classmethod
    def get_suitable_cooling_agent(cls, T_pinch: float):
        """Return a cooling agent that works at the pinch temperature.
        
        Parameters
        ----------
        T_pinch : 
            Pinch temperature [K].
        
        """
        for agent in cls.cooling_agents:
            if T_pinch > agent.T + agent.dT: return agent
        raise RuntimeError(f'no cooling agent that can cool under {T_pinch} K')    

    def load_agent(self, agent: UtilityAgent):
        """Initialize utility streams with given agent."""
        if self.agent is agent: 
            self.inlet_utility_stream.copy_like(agent)
            if not agent.T_limit: self.outlet_utility_stream.copy_thermal_condition(agent)
            return
        elif self.agent: 
            self.empty()
        if agent.utility_stream_dump:
            if agent.isfuel:
                self.inlet_utility_stream, self.outlet_utility_stream, self.oxygen_rich_inlet = agent.utility_stream_dump.pop()
                self.oxygen_rich_inlet.reset_flow(O2=21, N2=79, phase='g', units='kg/hr')
            else:
                self.inlet_utility_stream, self.outlet_utility_stream = agent.utility_stream_dump.pop()
            # Prevent errors where utility streams are altered
            self.inlet_utility_stream.copy_like(agent) 
            if not agent.T_limit: self.outlet_utility_stream.copy_thermal_condition(agent)
        else:
            self.inlet_utility_stream = agent.to_stream()
            if agent.isfuel:
                self.outlet_utility_stream = agent.to_stream()
                self.oxygen_rich_inlet = Stream(O2=21, N2=79, phase='g', units='kg/hr', thermo=agent.thermo)
            else:
                self.outlet_utility_stream = self.inlet_utility_stream.flow_proxy()
        self.agent = agent

    def mix_from(self, heat_utilities: Iterable[HeatUtility]):
        """Mix all heat utilities to this heat utility."""
        heat_utilities = [i for i in heat_utilities if i.agent]
        N_heat_utilities = len(heat_utilities)
        if N_heat_utilities == 0:
            self.empty()
        elif N_heat_utilities == 1:
            self.copy_like(heat_utilities[0])
        else:
            heat_utility, *other_heat_utilities = heat_utilities
            agent = heat_utility.agent
            ID = agent.ID
            for i in other_heat_utilities:
                if i.agent.ID != ID:
                    raise ValueError(
                        "utility agent must be the same to mix heat utilities"
                    )
            self.load_agent(agent)
            self.flow = self.inlet_utility_stream.F_mol = sum([i.flow for i in heat_utilities])
            self.outlet_utility_stream.mix_from([i.outlet_utility_stream for i in heat_utilities])
            self.duty = sum([i.duty for i in heat_utilities])
            self.unit_duty = sum([i.unit_duty for i in heat_utilities])
            self.cost = sum([i.cost for i in heat_utilities])
            self.heat_transfer_efficiency = None

    def reverse(self):
        """Reverse direction of utility. If utility is being consumed,
        the utility is produced instead, and vice-versa."""
        self.flow *= -1
        self.duty *= -1
        self.unit_duty *= -1
        self.cost *= -1
        if self.duty:
            self.inlet_utility_stream, self.outlet_utility_stream = self.outlet_utility_stream, self.inlet_utility_stream

    # Subcalculations

    @staticmethod
    def get_outlet_temperature(T_pinch: float, T_limit: float, iscooling: bool):
        """
        Return outlet temperature of the utility in a counter current heat exchanger

        Parameters
        ----------
        T_pinch : 
                  Pinch temperature of utility stream [K].
        iscooling : 
                  True if utility is loosing energy.

        """
        if iscooling:
            return T_limit if T_limit and T_limit < T_pinch else T_pinch
        else:
            return T_limit if T_limit and T_limit > T_pinch else T_pinch
        
    @classmethod
    def get_inlet_and_outlet_pinch_temperature(cls, iscooling: bool, T_in: float, T_out: float):
        """Return pinch inlet and outlet temperature of utility."""
        dT = cls.dT
        if iscooling:
            if T_in + 1e-1 < T_out:
                raise ValueError("inlet must be hotter than outlet if cooling")
            T_pinch_in = T_out - dT
            T_pinch_out = T_in - dT
        else:
            if T_in > T_out + 1e-1:
                raise ValueError("inlet must be cooler than outlet if heating")
            T_pinch_in = T_out + dT
            T_pinch_out = T_in + dT
        return T_pinch_in, T_pinch_out

    # Representation
    def _info_data(self, duty: None, flow: None, cost: None):
        # Get units of measure
        su = self.display_units
        duty_units = duty or su.duty
        flow_units = flow or su.flow
        cost_units = cost or su.cost
        
        # Change units and return info string
        flow = self.inlet_utility_stream.get_total_flow(flow_units)
        duty = convert(self.duty, 'kJ/hr', duty_units)
        cost = convert(self.cost, 'USD/hr', cost_units)
        return duty, copysign(flow, self.flow), cost, duty_units, flow_units, cost_units
        
    def __repr__(self):
        if self.agent:
            duty, flow, cost, duty_units, flow_units, cost_units = self._info_data(None, None, None)
            return f'<{self.ID}: {self.duty:.3g} {duty_units}, {self.flow:.3g} {flow_units}, {self.cost:.3g} {cost_units}>'
        else:
            return f'<{type(self).__name__}: None>'
        
    def _info(self, duty: str|None, flow: str|None, cost: str|None):
        """Return string related to specifications"""
        if not self.agent:
            return (f'{type(self).__name__}: None\n'
                    +' duty: 0\n'
                    +' flow: 0\n'
                    +' cost: 0')
        else:
            (duty, flow, cost, duty_units,
             flow_units, cost_units) = self._info_data(duty, flow, cost)
            return (f'{type(self).__name__}: {self.ID}\n'
                   +f'duty:{duty: .3g} {duty_units}\n'
                   +f'flow:{flow: .3g} {flow_units}\n'
                   +f'cost:{cost: .3g} {cost_units}')
            

    def show(self, duty: Optional[str]=None, flow: Optional[str]=None, cost: Optional[str]=None):
        """Print all specifications"""
        print(self._info(duty, flow, cost))
    _ipython_display_ = show

    def __add__(self, other: HeatUtility) -> HeatUtility:
        if other == 0: return self # Special case to get Python built-in sum to work
        return self.__class__.sum([self, other])
        
    def __radd__(self, other: HeatUtility) -> HeatUtility:
        return self.__add__(other)

HeatUtility.default_agents()
settings_cls = settings.__class__
settings_cls.get_agent = HeatUtility.get_agent
settings_cls.get_cooling_agent = HeatUtility.get_cooling_agent
settings_cls.get_heating_agent = HeatUtility.get_heating_agent
settings_cls.set_utility_agent_CF = HeatUtility.set_CF
del settings_cls
