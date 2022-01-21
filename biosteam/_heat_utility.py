# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from thermosteam.units_of_measure import (
    convert, DisplayUnits, AbsoluteUnitsOfMeasure, get_dimensionality,
    heat_utility_units_of_measure
)
from thermosteam.utils import unregistered, units_of_measure
from thermosteam import Thermo, Stream, ThermalCondition, settings
from .exceptions import DimensionError
from math import copysign

__all__ = ('HeatUtility', 'UtilityAgent')

# Costs from Table 17.1 in Warren D.  Seider et al.-Product and Process Design Principles_ Synthesis, Analysis and Evaluation-Wiley (2016)
# ^This table was made using data from Busche, 1995
# Entry temperature conditions of coolants taken from Table 12.1 in Warren, 2016

mol_basis_units = AbsoluteUnitsOfMeasure('kmol')
mass_basis_units = AbsoluteUnitsOfMeasure('kg')
energy_basis_units = AbsoluteUnitsOfMeasure('kJ')

# %% Utility agents

@unregistered
class UtilityAgent(Stream):
    """
    Create a UtilityAgent object that defines a utility option.
    
    Parameters
    ----------
    ID='' : str
        A unique identification. If ID is None, utility agent will not be registered.
        If no ID is given, utility agent will be registered with a unique ID.
    flow=() : tuple
        All flow rates corresponding to chemical `IDs`.
    phase='l' : 'l', 'g', or 's'
        Either gas (g), liquid (l), or solid (s).
    T=298.15 : float
        Temperature [K].
    P=101325 : float
        Pressure [Pa].
    units='kmol/hr' : str
        Flow rate units of measure (only mass, molar, and
        volumetric flow rates are valid).
    thermo=None : Thermo
        Thermodynamic equilibrium package. Defaults to `biosteam.settings.get_thermo()`.
    T_limit : float, optional
        Temperature limit of outlet utility streams [K]. If no limit is given,
        phase change is assumed. If utility agent heats up, `T_limit` is
        the maximum temperature. If utility agent cools down, `T_limit` is
        the minumum temperature. 
    heat_transfer_price : float
        Price of transfered heat [USD/kJ].
    regeneration_price : float
        Price of regenerating the fluid for reuse [USD/kmol].
    heat_transfer_efficiency : float
        Fraction of heat transfered accounting for losses to the environment (must be between 0 to 1).
    **chemical_flows : float
        ID - flow pairs.
        
    """
    __slots__ = ('T_limit', '_heat_transfer_price',
                 '_regeneration_price', 'heat_transfer_efficiency')
    def __init__(self, ID='', flow=(), phase='l', T=298.15, P=101325., units='kmol/hr',
                 thermo=None, T_limit=None, heat_transfer_price=0.0,
                 regeneration_price=0.0, heat_transfer_efficiency=1.0,                  
                 **chemical_flows):
        self._thermal_condition = ThermalCondition(T, P)
        thermo = self._load_thermo(thermo)
        self._init_indexer(flow, phase, thermo.chemicals, chemical_flows)
        if units != 'kmol/hr':
            name, factor = self._get_flow_name_and_factor(units)
            flow = getattr(self, name)
            flow[:] = self.mol / factor
        self._islinked = False
        self._sink = self._source = None
        self.reset_cache()
        self._register(ID)
        self.T_limit = T_limit
        self.heat_transfer_price = heat_transfer_price
        self.regeneration_price = regeneration_price
        self.heat_transfer_efficiency = heat_transfer_efficiency
    
    def iscooling_agent(self):
        T_limit = self.T_limit 
        return T_limit > self.T if T_limit else self.phase == 'l'
    
    def isheating_agent(self):
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
    
    def to_stream(self, ID=None):
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
        new._islinked = False
        new._sink = new._source = None
        new._thermo = self._thermo
        new._imol = self._imol.copy()
        new._thermal_condition = self._thermal_condition.copy()
        new.reset_cache()
        new._price = 0.
        new.characterization_factors = {}
        new.ID = ID
        return new
    
    @property
    def heat_transfer_price(self):
        """Price of transfered heat [USD/kJ]."""
        return self._heat_transfer_price
    @heat_transfer_price.setter
    def heat_transfer_price(self, price):
        assert price >= 0, "heat transfer price cannot be negative"
        self._heat_transfer_price = price
    
    @property
    def regeneration_price(self):
        """Price of regenerating the fluid for reuse [USD/kmol]."""
        return self._regeneration_price
    @regeneration_price.setter
    def regeneration_price(self, price):
        assert price >= 0, "regeneration price cannot be negative"
        self._regeneration_price = price
    
    def _info_phaseTP(self, phase, T_units, P_units):
        T_limit = convert(self.T_limit, 'K', T_units)
        T_limit = f"{T_limit:.3g} K" if T_limit else "None"
        T = convert(self.T, 'K', T_units)
        P = convert(self.P, 'Pa', P_units)
        s = '' if isinstance(phase, str) else 's'
        ht_price = self.heat_transfer_price
        rg_price = self.regeneration_price
        ht_eff = self.heat_transfer_efficiency
        return (f" heat_transfer_efficiency: {ht_eff:.3f}\n"
                f" heat_transfer_price: {ht_price:.3g} USD/kJ\n"
                f" regeneration_price: {rg_price:.3g} USD/kmol\n"
                f" T_limit: {T_limit}\n"
                f" phase{s}: {repr(phase)}\n"
                f" T: {T:.5g} {T_units}\n"
                f" P: {P:.6g} {P_units}\n"
        )
    def __repr__(self):
        return f"<{type(self).__name__}: {self.ID}>"

       
# %%

@units_of_measure(heat_utility_units_of_measure)
class HeatUtility:
    """
    Create an HeatUtility object that can choose a utility stream and 
    calculate utility requirements. It can calculate required flow rate,
    temperature change, or phase change of utility. Calculations assume 
    counter current flow rate.
    
    Parameters
    ----------
    heat_transfer_efficiency=None : float, optional
        Enforced fraction of heat transfered from utility (due
        to losses to environment).
    
    Attributes
    ----------
    cooling_agents : list[UtilityAgent]
        All cooling utilities available.
    heating_agents : list[UtilityAgent]
        All heating utilities available.
    duty : float
        Total heat transfered from utility to both the process and the environment [kJ/hr].
    unit_duty : float
        Effective heat transfered from utility to the unit operation [kJ/hr].
    flow : float
        Flow rate of utility [kmol/hr].
    cost : float
        Cost of utility [USD/hr].
    heat_transfer_efficiency : float
        Enforced fraction of heat transfered from utility (due
        to losses to environment).
    T_pinch : float
        If cooling, the minimum utility temperature required.
        If heating, the maximum utility temperature required.
    iscooling : float
        Whether the utility is cooling the process.
    agent : UtilityAgent
        Utility agent being used.
    inlet_utility_stream : :class:`~thermosteam.Stream`
    outlet_utility_stream : :class:`~thermosteam.Stream`
    
    Examples
    --------
    Create a heat utility:
        
    .. code-block:: python
    
       >>> from biosteam import HeatUtility
       >>> hu = HeatUtility()
       >>> hu.show()
       HeatUtility: None
        duty: 0
        flow: 0
        cost: 0
    
    Calculate utility requirement by calling it with a duty (kJ/hr), and entrance and exit temperature (K):
        
    .. code-block:: python
    
       >>> hu(1000, 300, 350)
       >>> hu.show()
       HeatUtility: low_pressure_steam
        duty: 1.05e+03 kJ/hr
        flow: 0.0271 kmol/hr
        cost: 0.00645 USD/hr
   
    All results are accessible:
        
    .. code-block:: python
    
       >>> hu.ID, hu.duty, hu.flow, hu.cost
       ('low_pressure_steam', 1052., 0.02711, 0.006448)
           
    """
    __slots__ = ('inlet_utility_stream', 'outlet_utility_stream', 'duty',
                 'flow', 'cost', 'heat_transfer_efficiency', 'T_pinch',
                 'heat_exchanger', 'iscooling', 'agent', 'unit_duty')
    #: [float] Minimum approach temperature difference. Used to assign
    #: the pinch temperature of the utility stream.
    dT = 5  
    
    #: [DisplayUnits] Units of measure for IPython display.
    display_units = DisplayUnits(duty='kJ/hr', flow='kmol/hr', cost='USD/hr')

    #: [Thermo] Used broadly throughout BioSTEAM utilities.
    thermo_water = Thermo(['Water'])
    
    #: [Thermo] Used for BioSTEAM refrigeration utility.
    thermo_propane = Thermo(['Propane'])

    #: dict[tuple[str, str], tuple[float, AbsoluteUnitsOfMeasure]] 
    #: Characterization factor data (value and units) by agent ID and impact key.
    characterization_factors = {}
    
    @classmethod
    def set_CF(self, ID, key, value, basis=None, units=None):
        """
        Set the cacharacterization factor of a utility agent for a given impact 
        key.

        Parameters
        ----------
        ID : str
            ID of utility agent.
        key : str
            Key of impact indicator.
        value : float
            Characterization factor value.
        basis : str, optional
            Basis of characterization factor. Valid dimensions include weight, 
            molar, and energy (e.g. 'kg', 'kmol', 'kJ'). Defaults to 'kg'.
        units : str, optional
            Units of impact indicator. Before using this argument, the default units 
            of the impact indicator should be defined with 
            thermosteam.settings.define_impact_indicator.
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
            original_units = settings.get_impact_indicator_units(units)
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
                if agent.iscooling_agent() and value > 0.: 
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
    def get_CF(self, ID, key, basis=None, units=None):
        """
        Return the characterization factor of a utility agent for a given impact 
        key.

        Parameters
        ----------
        ID : str
            ID of utility agent.
        key : str
            Key of impact indicator.
        basis : str, optional
            Basis of characterization factor. Valid dimensions include weight, 
            molar, and energy (e.g. 'kg', 'kmol', 'kJ'). Defaults to 'kg'.
        units : str, optional
            Units of impact indicator. Before using this argument, the default units 
            of the impact indicator should be defined with 
            thermosteam.settings.define_impact_indicator.
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
            original_units = settings.get_impact_indicator_units(units)
            value = original_units.convert(value, units)
        if basis is None:
            return value, basis_units._units
        else:
            return value / basis_units.conversion_factor(basis)

    def get_impact(self, key):
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

    def get_inventory(self, key):
        agent = self.agent
        CF, basis_units = self.get_CF(agent.ID, key)
        if basis_units == 'kg':
            return self.flow * agent.MW, basis_units
        elif basis_units == 'mol':
            return self.flow, basis_units
        elif basis_units == 'kJ':
            return self.duty, basis_units
        else:
            raise RuntimeError("unknown error")

    def __init__(self, heat_transfer_efficiency=None, heat_exchanger=None):
        self.heat_transfer_efficiency = heat_transfer_efficiency
        self.heat_exchanger = heat_exchanger
        self.empty()

    def __bool__(self):
        return bool(self.agent)

    @property
    def ID(self):
        """[str] ID of utility agent being used."""
        agent = self.agent
        return agent.ID if agent else ""
    @ID.setter
    def ID(self, ID):
        """[str] ID of utility agent being used."""
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
        cls.heating_agents = [low_pressure_steam, 
                              medium_pressure_steam, 
                              high_pressure_steam]
        
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
        )
        chilled_water = UtilityAgent(
            'chilled_water',
            Water=1, T=280.372, P=101325,
            thermo=thermo_water,
            T_limit = 300.372,
            heat_transfer_price = 5e-6
        )
        chilled_brine = UtilityAgent(
            'chilled_brine',
            Water=1, T=255.372, P=101325,
            thermo=thermo_water,
            T_limit = 275.372,
            heat_transfer_price = 8.145e-6,
        )
        propane = UtilityAgent(
            'propane',
            Propane=1,
            thermo=cls.thermo_propane,
            T=273.15 - 42.11,
            P=101325,
            heat_transfer_price = 13.17e-6,
            phase='l',
        )
        cls.cooling_agents = [cooling_water, 
                              chilled_water, 
                              chilled_brine,
                              propane]

    def copy(self):
        hu = HeatUtility()
        hu.copy_like(self)
        return hu

    def copy_like(self, other):
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
        self.T_pinch = other.T_pinch
        self.iscooling = other.iscooling

    def scale(self, factor):
        """Scale utility data."""
        self.flow *= factor
        self.duty *= factor
        self.cost *= factor
        self.inlet_utility_stream.mol *= factor
        # No need to factor the outlet utility stream
        # because it shares the same flow rate data as the inlet

    def empty(self):
        """Remove utility requirements."""
        self.cost = self.flow = self.duty = self.unit_duty = 0
        self.iscooling = self.agent = self.T_pinch = None
        
    def set_utility_by_flow_rate(self, agent, F_mol):
        if F_mol == 0.: 
            self.empty()
            return
        self.load_agent(agent)
        heat_transfer_efficiency = self.heat_transfer_efficiency or agent.heat_transfer_efficiency
        if agent.T_limit: raise ValueError('agent must work by latent heat to set by flow rate')
        self.outlet_utility_stream.phase = 'g' if self.inlet_utility_stream.phase == 'l' else 'l'
        dH = self.outlet_utility_stream.Hvap
        self.duty = duty = dH * F_mol
        self.unit_duty = duty * heat_transfer_efficiency
        self.outlet_utility_stream.mol[:] = F_mol
        self.flow = F_mol
        self.cost = agent._heat_transfer_price * abs(duty) + agent._regeneration_price * F_mol
        
    def __call__(self, unit_duty, T_in, T_out=None, agent=None):
        """Calculate utility requirements given the essential parameters.
        
        Parameters
        ----------
        unit_duty : float
               Unit duty requirement (kJ/hr)
        T_in : float
               Inlet process stream temperature (K)
        T_out : float, optional
               Outlet process stream temperature (K)
        agent : UtilityAgent, optional
                Utility agent to use. Defaults to a suitable agent from 
                predefined heating/cooling utility agents.
        
        """
        if unit_duty == 0:
            self.empty()
            return
        T_out = T_out or T_in
        iscooling = unit_duty < 0.
        
        # Note: These are pinch temperatures at the utility inlet and outlet. 
        # Not to be confused with the inlet and outlet of the process stream.
        T_pinch_in, T_pinch_out = self.get_inlet_and_outlet_pinch_temperature(
            iscooling, T_in, T_out)

        ## Select heat transfer agent ##
        if agent:
            self.load_agent(agent)
        elif self.iscooling == iscooling and self.T_pinch == T_pinch_in:
            agent = self.agent
        else:
            self.iscooling = iscooling
            self.T_pinch = T_pinch_in
            if iscooling:
                agent = self.get_suitable_cooling_agent(T_pinch_in)
            else:
                agent = self.get_suitable_heating_agent(T_pinch_in)
            self.load_agent(agent)
            
        ## Calculate utility requirement ##
        heat_transfer_efficiency = self.heat_transfer_efficiency or agent.heat_transfer_efficiency
        duty = unit_duty/heat_transfer_efficiency
        if agent.T_limit:
            # Temperature change
            self.outlet_utility_stream.T = self.get_outlet_temperature(
                T_pinch_out, agent.T_limit, iscooling)
            dH = self.inlet_utility_stream.H - self.outlet_utility_stream.H
        else:
            # Phase change
            if self.inlet_utility_stream.phase == 'l' :
                self.outlet_utility_stream.phase = 'g'
                dH = -self.outlet_utility_stream.Hvap
            else:
                dH = self.outlet_utility_stream.Hvap
                self.outlet_utility_stream.phase = 'l'
        
        # Update utility flow
        self.outlet_utility_stream.mol[:] *= duty / dH
        
        # Update results
        self.unit_duty = unit_duty
        self.flow = F_mol = self.inlet_utility_stream.F_mol
        self.duty = duty
        self.cost = agent._heat_transfer_price * abs(duty) + agent._regeneration_price * F_mol

    @property
    def inlet_process_stream(self):
        """[:class`~thermosteam.Stream`] If a heat exchanger is available,
        this stream is the inlet process stream to the heat exchanger."""
        heat_exchanger = self.heat_exchanger
        if heat_exchanger:
            return heat_exchanger.ins[0]
        else:
            raise AttributeError('no heat exchanger available '
                                 'to retrieve process stream')

    @property
    def outlet_process_stream(self):
        """[:class`~thermosteam.Stream`] If a heat exchanger is available,
        this stream is the outlet process stream to the heat exchanger."""
        heat_exchanger = self.heat_exchanger
        if heat_exchanger:
            return heat_exchanger.outs[0]
        else:
            raise AttributeError('no heat exchanger available '
                                 'to retrieve process stream')
    
    @staticmethod
    def heat_utilities_by_agent(heat_utilities):
        """Return a dictionary of heat utilities sorted by agent."""
        heat_utilities = [i for i in heat_utilities if i.agent]
        heat_utilities_by_agent = {i.agent: [] for i in heat_utilities}
        for i in heat_utilities:
            heat_utilities_by_agent[i.agent].append(i)
        return heat_utilities_by_agent

    @classmethod
    def sum(cls, heat_utilities):
        """
        Return a HeatUtility object that reflects the sum of heat
        utilities.
        """
        heat_utility = cls()
        heat_utility.mix_from(heat_utilities)
        return heat_utility    

    @classmethod
    def sum_by_agent(cls, heat_utilities, agent=None):
        """
        Return a tuple of heat utilities that reflect the sum of heat utilities
        by agent. If an agent is specified, return a HeatUtility object that
        reflects the sum of all heat utilities with the given agent.
        
        """
        heat_utilities_by_agent = cls.heat_utilities_by_agent(heat_utilities)
        if agent:
            return cls.sum(heat_utilities_by_agent[agent])
        else:
            return tuple([cls.sum(i) for i in heat_utilities_by_agent.values()])

    @classmethod
    def get_agent(cls, ID):
        """Return utility agent with given ID."""
        for agent in cls.heating_agents + cls.cooling_agents:
            if agent.ID == ID: return agent
        raise LookupError(ID)

    @classmethod
    def get_heating_agent(cls, ID):
        """Return heating agent with given ID."""
        for agent in cls.heating_agents:
            if agent.ID == ID: return agent
        raise LookupError(ID)
  
    @classmethod
    def get_cooling_agent(cls, ID):
        """Return cooling agent with given ID."""
        for agent in cls.cooling_agents:
            if agent.ID == ID: return agent
        raise LookupError(ID)
    
    @classmethod
    def get_suitable_heating_agent(cls, T_pinch):
        """
        Return a heating agent that works at the pinch temperature.
        
        Parameters
        ----------
        T_pinch : float
            Pinch temperature [K].
        
        """
        for agent in cls.heating_agents:
            if T_pinch < agent.T: return agent
        raise RuntimeError(f'no heating agent that can heat over {T_pinch} K')    

    @classmethod
    def get_suitable_cooling_agent(cls, T_pinch):
        """Return a cooling agent that works at the pinch temperature.
        
        Parameters
        ----------
        T_pinch : float
            Pinch temperature [K].
        
        """
        for agent in cls.cooling_agents:
            if T_pinch > agent.T: return agent
        raise RuntimeError(f'no cooling agent that can cool under {T_pinch} K')    

    def load_agent(self, agent):
        """Initialize utility streams with given agent."""
        # Initialize streams
        self.inlet_utility_stream = agent.to_stream()
        self.outlet_utility_stream = self.inlet_utility_stream.flow_proxy()
        self.agent = agent

    def mix_from(self, heat_utilities):
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
            assert all([i.agent.ID == ID for i in other_heat_utilities]), (
                "utility agent must be the same to mix heat utilities")
            self.load_agent(agent)
            self.flow = self.inlet_utility_stream.F_mol = sum([i.flow for i in heat_utilities])
            self.duty = sum([i.duty for i in heat_utilities])
            self.unit_duty = sum([i.unit_duty for i in heat_utilities])
            self.cost = sum([i.cost for i in heat_utilities])
            self.heat_transfer_efficiency = None
            self.T_pinch = None
            self.iscooling = None

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
    def get_outlet_temperature(T_pinch, T_limit, iscooling):
        """
        Return outlet temperature of the utility in a counter current heat exchanger

        Parameters
        ----------
        T_pinch : float
                  Pinch temperature of utility stream [K].
        iscooling : bool
                  True if utility is loosing energy.

        """
        if iscooling:
            return T_limit if T_limit and T_limit < T_pinch else T_pinch
        else:
            return T_limit if T_limit and T_limit > T_pinch else T_pinch
        
    @classmethod
    def get_inlet_and_outlet_pinch_temperature(cls, iscooling, T_in, T_out):
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

    def _info_data(self, duty, flow, cost):
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
        
    # Representation
    def _info(self, duty, flow, cost):
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
                    +f' duty:{duty: .3g} {duty_units}\n'
                    +f' flow:{flow: .3g} {flow_units}\n'
                    +f' cost:{cost: .3g} {cost_units}')
            

    def show(self, duty=None, flow=None, cost=None):
        """Print all specifications"""
        print(self._info(duty, flow, cost))
    _ipython_display_ = show

HeatUtility.default_agents()