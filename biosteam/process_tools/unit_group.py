# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

import pandas as pd
import biosteam as bst
import numpy as np
from thermosteam.utils import repr_items
from . import utils
from .._unit import Unit
from ..utils import streams_from_units
from .._heat_utility import HeatUtility
from collections.abc import Mapping

__all__ = ('UnitGroup', 'GroupedUnit', 'create_connected_grouped_units')

INST_EQ_COST = 'Inst. eq. cost'
ELEC_CONS = 'Elec. cons.'
ELEC_PROD = 'Elec. prod.'
INSTALLED_EQUIPMENT_COST = 'Installed equipment cost'
COOLING = 'Cooling'
HEATING = 'Heating'
COOLING_DUTY = 'Cooling duty'
HEATING_DUTY = 'Heating duty'
ELECTRICITY_CONSUMPTION = 'Electricity consumption'
ELECTRICITY_PRODUCTION = 'Electricity production'
MATERIAL_COST = 'Material cost'
MAT_COST = 'Mat. cost'
CAPITAL_UNITS = 'MM$'
ELEC_UNITS = 'MW'
DUTY_UNITS = 'GJ/hr'
MAT_UNITS = 'USD/hr'

# %% Unit group for generating results

class UnitGroup:
    """
    Create a UnitGroup object for generating biorefinery results.
    
    Parameters
    ----------
    name : str, optional
        Name of group for bookkeeping.
    units : tuple[Unit], optional
        Unit operations.
    metrics=None : list[Metric], optional
        Metrics to generate results. These metrics are computed when 
        generating results as dictionaries, pandas series, and data frames.
    
    Examples
    --------
    Create a UnitGroup from BioSTEAM's example ethanol subsystem:
    
    >>> from biorefineries.sugarcane import chemicals
    >>> from biosteam import *
    >>> settings.set_thermo(chemicals)
    >>> water = Stream('water', Water=100., T=350.)
    >>> sucrose = Stream('sucrose', Sucrose=3.)
    >>> with System('example_sys') as example_sys:
    ...     P1 = Pump('P1', water)
    ...     T1 = MixTank('T1', (P1-0, sucrose))
    ...     H1 = HXutility('H1', T1-0, T=300.)
    >>> example_sys.simulate()
    >>> ugroup = UnitGroup('Example group', example_sys.units)
    
    We can autofill metrics to evaluate:
    
    >>> ugroup.autofill_metrics(electricity_production=True)
    >>> ugroup.metrics
    [<Indicator: Installed equipment cost (MM$)>, <Indicator: Cooling duty (GJ/hr)>, <Indicator: Heating duty (GJ/hr)>, <Indicator: Electricity consumption (MW)>, <Indicator: Electricity production (MW)>, <Indicator: Material cost (USD/hr)>]
    
    Get all metric results:
        
    >>> ugroup.to_dict()
    {'Installed equipment cost [MM$]': 0.05, 'Cooling duty [GJ/hr]': 0.37, 'Heating duty [GJ/hr]': 0.0, 'Electricity consumption [MW]': 0.0, 'Electricity production [MW]': 0.0, 'Material cost [USD/hr]': 0.0}
    
    Each result can be retrieved separately:
    
    >>> ugroup.get_installed_cost()
    0.05
    
    >>> ugroup.get_cooling_duty()
    0.37
    
    The `to_dict` method also returns user-defined metrics:
        
    >>> # First define metrics
    >>> @ugroup.metric # Name of metric defaults to function name
    ... def moisture_content():
    ...     product = H1.outs[0]
    ...     return product.imass['Water'] / product.F_mass
    
    >>> @ugroup.metric(units='kg/hr') # This helps for bookkeeping
    ... def sucrose_flow_rate():
    ...     return float(H1.outs[0].imass['Sucrose'])
    
    >>> ugroup.show()
    UnitGroup: Example group
     units: P1, T1, H1
     metrics: Installed equipment cost [MM$]
              Cooling duty [GJ/hr]
              Heating duty [GJ/hr]
              Electricity consumption [MW]
              Electricity production [MW]
              Material cost [USD/hr]
              Moisture content
              Sucrose flow rate [kg/hr]
    
    >>> ugroup.to_dict()
    {'Installed equipment cost [MM$]': 0.05, 'Cooling duty [GJ/hr]': 0.37, 'Heating duty [GJ/hr]': 0.0, 'Electricity consumption [MW]': 0.0, 'Electricity production [MW]': 0.0, 'Material cost [USD/hr]': 0.0, 'Moisture content': 0.63, 'Sucrose flow rate [kg/hr]': 1026.8}
    
    """
    __slots__ = ('name', 'units', 'metrics', 'filter_savings', 'extend_feed_ends')
    
    def __init__(self, name=None, units=(), metrics=None, 
                 filter_savings=True, extend_feed_ends=True):
        #: [str] Name of group for bookkeeping
        self.name = 'Unnamed' if name is None else str(name)
        
        #: list[Unit] Unit operations
        self.units = units if isinstance(units, list) else list(units) 
        
        if metrics is None: 
            metrics = []
        elif not isinstance(metrics, list):
            metrics = list(metrics)
        
        #: list[Metric] Metrics to generate results
        self.metrics = metrics
        
        #: [bool] Whether to only allow postive flows in utility results
        self.filter_savings = filter_savings
        
        #: [bool] Whether to consider feeds past external storage, pumps, and heat 
        #: exchangers for calculating material costs.
        self.extend_feed_ends = extend_feed_ends
    
    def autofill_metrics(self, shorthand=False, 
                         installed_cost=True,
                         cooling_duty=True,
                         heating_duty=True,
                         electricity_consumption=True,
                         electricity_production=False,
                         material_cost=True):
        if installed_cost:
            self.metric(self.get_installed_cost,
                        INST_EQ_COST if shorthand else INSTALLED_EQUIPMENT_COST,
                        CAPITAL_UNITS)
        if cooling_duty:
            self.metric(self.get_cooling_duty,
                        COOLING if shorthand else COOLING_DUTY,
                        DUTY_UNITS)
        if heating_duty:
            self.metric(self.get_heating_duty,
                        HEATING if shorthand else HEATING_DUTY,
                        DUTY_UNITS)
        if electricity_consumption:
            self.metric(self.get_electricity_consumption,
                        ELEC_CONS if shorthand else ELECTRICITY_CONSUMPTION,
                        ELEC_UNITS)
        if electricity_production:
            self.metric(self.get_electricity_production,
                        ELEC_PROD if shorthand else ELECTRICITY_PRODUCTION,
                        ELEC_UNITS)
        if material_cost:
            self.metric(self.get_material_cost,
                        MAT_COST if shorthand else MATERIAL_COST,
                        MAT_UNITS)
    
    def __iter__(self):
        return iter(self.units)
    
    def to_unit(self, ID=None, thermo=None):
        return GroupedUnit(ID, units=self.units)
    
    def to_system(self, ID=None):
        """Return a System object of all units."""
        return bst.System.from_units(ID, self.units)
    
    def split(self, stream, upstream_name=None, downstream_name=None):
        """
        Split unit group in two; upstream and downstream.
        
        Parameters
        ----------    
        stream : Iterable[:class:~thermosteam.Stream], optional
            Stream where unit group will be split.
        upstream_name : str, optional
            Name of upstream UnitGroup object.
        downstream_name : str, optional
            Name of downstream UnitGroup object.
        
        Examples
        --------
        >>> from biorefineries.cornstover import cornstover_sys, M201
        >>> from biosteam import default
        >>> ugroup = cornstover_sys.to_unit_group()
        >>> upstream, downstream = ugroup.split(M201-0)
        >>> upstream.show()
        UnitGroup: Unnamed
         units: U101, H2SO4_storage, T201, M201
        >>> for i in upstream: assert i not in downstream.units
        >>> assert set(upstream.units + downstream.units) == set(cornstover_sys.units)
        >>> default() # Reset to biosteam defaults
        
        """
        sys = self.to_system()
        units = sys.units
        index = units.index(stream.sink)
        return (UnitGroup(upstream_name, units[:index]),
                UnitGroup(downstream_name, units[index:]))
    
    def metric(self, getter=None, name=None, units=None, element=None):
        """
        Define and register metric.
        
        Parameters
        ----------    
        getter : function, optional
                 Should return metric.
        name : str, optional
               Name of parameter. If None, defaults to the name of the getter.
        units : str, optional
                Parameter units of measure
        
        Notes
        -----
        This method works as a decorator.
        
        """
        if element is None: element = ''
        if not getter: return lambda getter: self.metric(getter, name, units, element)
        metric = bst.Metric(name, getter, units, element)
        bst.Metric.check_index_unique(metric, self.metrics)
        self.metrics.append(metric)
        return metric 
    
    def register_utility_agent(self, agent, basis='duty'):
        """Register utility agent as a metric to UnitGroup."""
        if isinstance(agent, str): agent = HeatUtility.get_agent(agent)
        name = agent.ID.replace('_', ' ').capitalize()
        if basis == 'duty':
            self.metric(lambda: self.get_utility_duty(agent), name, 'GJ')
        elif basis == 'flow':
            self.metric(lambda: self.get_utility_flow(agent), name, 'MT')
        else:
            raise ValueError(f"basis must be either 'duty' or 'flow', not {repr(basis)}")
    
    @property
    def heat_utilities(self):
        """[tuple] All HeatUtility objects."""
        heat_utilities = utils.get_heat_utilities(self.units)
        if self.filter_savings:
            return utils.filter_out_heat_utility_savings(heat_utilities) 
        else: 
            return heat_utilities 
    
    @property
    def power_utilities(self):
        """[tuple] All PowerUtility objects."""
        return tuple(utils.get_power_utilities(self.units))
    
    @classmethod
    def filter_by_types(cls, name, units, types):
        """Create a UnitGroup object of given type(s)."""
        return cls(name, utils.filter_by_types(units, types))
    
    @classmethod
    def filter_by_lines(cls, name, units, lines):
        """Create a UnitGroup object of given line(s)."""
        return cls(name, utils.filter_by_lines(units, lines))
    
    @classmethod
    def group_by_types(cls, units, name_types=None):
        """Create a list of UnitGroup objects for each name-type pair."""
        if name_types:
            if isinstance(name_types, Mapping): name_types = name_types.items()
            return [cls.filter_by_types(name, units, types) for name, types in name_types]
        else:
            groups = utils.group_by_types(units)
            return [cls(i, j) for i, j in groups.items()]
    
    @classmethod
    def group_by_lines(cls, units, name_lines=None):
        """Create a list of UnitGroup objects for each name-line pair."""
        if name_lines:
            if isinstance(name_lines, Mapping): name_lines = name_lines.items()
            return [cls.filter_by_lines(name, units, lines) for name, lines in name_lines]
        else:
            groups = utils.group_by_lines(units)
            return [cls(i, j) for i, j in groups.items()]
    
    @classmethod
    def group_by_area(cls, units):
        """
        Create a list of UnitGroup objects for each area available.
        
        Examples
        --------
        >>> from biosteam import *
        >>> from biorefineries.cornstover import cornstover_sys
        >>> areas = UnitGroup.group_by_area(cornstover_sys.units)
        >>> areas[-1].show()
        UnitGroup: 700
         units: T701, P701, T702, P702, M701,
                T703
        
        >>> default() # Bring biosteam settings back to default
        """
        return [cls(i, j) for i, j in utils.group_by_area(units).items()]
    
    def get_inlet_flow(self, units, key=None):
        """
        Return total flow across all inlets.
        
        Parameters
        ----------
        units : str
            Units of measure.
        key : tuple[str] or str, optional
            Chemical identifiers. If none given, the sum of all chemicals returned
            
        Examples
        --------
        >>> from biosteam import Stream, Mixer, Splitter, UnitGroup, settings, main_flowsheet
        >>> settings.set_thermo(['Water', 'Ethanol'])
        >>> main_flowsheet.clear()
        >>> S1 = Splitter('S1', Stream(Ethanol=10, units='ton/hr'), split=0.1)
        >>> M1 = Mixer('M1', ins=[Stream(Water=10, units='ton/hr'), S1-0])
        >>> sys = main_flowsheet.create_system(operating_hours=330*24)
        >>> ugroup = UnitGroup('Example group', sys.units)
        >>> ugroup.get_inlet_flow('ton/hr') # Sum of all chemicals
        20.0
        >>> ugroup.get_inlet_flow('ton/hr', 'Water') # Just water
        10.0
        
        """
        if key:
            return sum([i.get_flow(units, key) for i in bst.utils.feeds_from_units(self.units)])
        else:
            return sum([i.get_total_flow(units) for i in bst.utils.feeds_from_units(self.units)])
    
    def get_outlet_flow(self, units, key=None):
        """
        Return total flow across all outlets.
        
        Parameters
        ----------
        units : str
            Units of measure.
        key : tuple[str] or str, optional
            Chemical identifiers. If none given, the sum of all chemicals returned
            
        Examples
        --------
        >>> from biosteam import Stream, Mixer, Splitter, UnitGroup, settings, main_flowsheet
        >>> settings.set_thermo(['Water', 'Ethanol'])
        >>> main_flowsheet.clear()
        >>> S1 = Splitter('S1', Stream(Ethanol=10, units='ton/hr'), split=0.1)
        >>> M1 = Mixer('M1', ins=[Stream(Water=10, units='ton/hr'), S1-0])
        >>> sys = main_flowsheet.create_system(operating_hours=330*24)
        >>> sys.simulate()
        >>> ugroup = UnitGroup('Example group', sys.units)
        >>> ugroup.get_inlet_flow('ton/hr') # Sum of all chemicals
        20.0
        >>> ugroup.get_inlet_flow('ton/hr', 'Water') # Just water
        10.0
        """
        if key:
            return sum([i.get_flow(units, key) for i in bst.utils.products_from_units(self.units)])
        else:
            return sum([i.get_total_flow(units) for i in bst.utils.products_from_units(self.units)])
    
    def get_material_cost(self):
        """Return the total material cost in USD/hr"""
        inlets = bst.utils.feeds_from_units(self.units)
        inlets = set(inlets)
        bst.utils.filter_out_missing_streams(inlets)
        if self.extend_feed_ends:
            get_inlet_origin = bst.utils.get_inlet_origin
            inlets = [get_inlet_origin(i) for i in inlets]
        feeds = bst.utils.feeds(inlets)
        return sum([i.cost for i in feeds])
    
    def get_utility_duty(self, agent):
        """Return the total utility duty for given agent in GJ/hr"""
        return utils.get_utility_duty(self.heat_utilities, agent)
    
    def get_utility_flow(self, agent):
        """Return the total utility flow for given agent in MT/hr"""
        return utils.get_utility_flow(self.heat_utilities, agent)
    
    def get_cooling_duty(self):
        """Return the total cooling duty in GJ/hr."""
        return utils.get_cooling_duty(self.heat_utilities)
    
    def get_heating_duty(self):
        """Return the total heating duty in GJ/hr."""
        return utils.get_heating_duty(self.heat_utilities)
    
    def get_installed_cost(self):
        """Return the total installed equipment cost in million USD."""
        return utils.get_installed_cost(self.units)
    
    def get_purchase_cost(self):
        """Return the total equipment purchase cost in million USD."""
        return utils.get_purchase_cost(self.units)
    
    def get_electricity_consumption(self):
        """Return the total electricity consumption in MW."""
        return utils.get_electricity_consumption(self.power_utilities)

    def get_electricity_production(self):
        """Return the total electricity production in MW."""
        return utils.get_electricity_production(self.power_utilities)
    
    def get_net_electricity_production(self):
        """Return the net electricity production in MW."""
        power_utilities = self.power_utilities
        return (utils.get_electricity_production(power_utilities)
                - utils.get_electricity_consumption(power_utilities))
    
    def to_dict(self, with_units=True):
        """Return dictionary of results."""
        metrics = self.metrics
        if not metrics: self.autofill_metrics()
        if with_units:
            return {i.name_with_units: i() for i in self.metrics}
        else:
            return {i.name: i() for i in self.metrics}
            
    def diagram(self, *args, **kwargs):
        return bst.System(None, self.units).diagram(*args, **kwargs)
    
    def to_series(self, with_units=True):
        """Return a pandas.Series object of metric results."""
        return pd.Series(self.to_dict(with_units), name=self.name)

    @classmethod
    def df_from_groups(cls, unit_groups, fraction=False, scale_fractions_to_positive_values=True):
        """
        Return a pandas DataFrame object of metric results from unit groups.
        
        Parameters
        ----------
        unit_groups : Sequence[UnitGroup]
            Metric results will be calculated from unit groups.
        fraction : bool, optional.
            Whether to divide metric results by the total sum across all groups. 
        scale_fractions_to_positive_values : bool, optional.
            Whether to compute fractions by dividing results by the sum of only 
            positive results.
        s
        Examples
        --------
        Create a pandas DataFrame of the net eletricity production across
        all areas in the sugarcane biorefinery:
        
        >>> import biosteam as bst
        >>> from biorefineries import sugarcane as sc
        >>> sc.load()
        >>> unit_groups = bst.UnitGroup.group_by_area(sc.sys.units)
        >>> for i in unit_groups: 
        ...     metric = i.metric(i.get_net_electricity_production, 
        ...                       'Net electricity production', 'kW')
        >>> bst.UnitGroup.df_from_groups(
        ...     unit_groups, fraction=True,
        ...     scale_fractions_to_positive_values=True,
        ... )
             Net electricity production
        0                           100
        100                       -2.99
        200                       -3.52
        300                      -0.913
        400                           0
        
        >>> bst.UnitGroup.df_from_groups(
        ...     unit_groups, fraction=True,
        ...     scale_fractions_to_positive_values=False,
        ... )
             Net electricity production
        0                           108
        100                       -3.23
        200                       -3.81
        300                      -0.986
        400                           0
        
        >>> bst.default() # Reset to biosteam defaults
        
        """
        with_units = not fraction
        data = [i.to_series(with_units) for i in unit_groups]
        df = pd.DataFrame(data)
        if fraction:
            values = df.values
            if scale_fractions_to_positive_values:
                postive_values = np.where(values > 0., values, 0.)
                values *= 100 / postive_values.sum(axis=0, keepdims=True)
            else:
                values *= 100 / values.sum(axis=0, keepdims=True)
        return df

    @classmethod
    def df_from_groups_across_coordinate(cls, unit_groups, f, xs, name=None):
        def get_df(x):
            f(x)
            return cls.df_from_groups(unit_groups)
        dfs = [get_df(x) for x in xs]
        df0 = dfs[0]
        columns = list(df0)
        data = sum([[df[i] for df in dfs] for i in columns], [])
        df = pd.DataFrame(np.array(data).transpose(),
            index=df0.index,
            columns=pd.MultiIndex.from_product([columns, xs], names=['Metric', name]),
        )
        return df

    def show(self):
        units = self.units
        if units:
            units = '\n' + repr_items(' units: ', units)
        else:
            units = "\n units: (No units)"
        metrics = self.metrics
        if metrics:
            metric_newline = "\n" + " " * len(' metrics: ')
            metrics = f"\n metrics: {metric_newline.join([i.describe() for i in self.metrics])}"
        else:
            metrics = ""
        print (
            f"{type(self).__name__}: {self.name}"
            + units
            + metrics
        )
        
    _ipython_display_ = show

    def __repr__(self):
        return f"{type(self).__name__}({repr(self.name)}, {self.units}, metrics={self.metrics})"

# %% UnitGroup as an actual unit operation

def create_connected_grouped_units(unit_groups):
    grouped_units = [GroupedUnit('.' + i.name, units=i) for i in unit_groups]
    utils.connect_by_ID(grouped_units)
    return grouped_units

class GroupedUnit(Unit):
    line = ''
    _ins_size_is_fixed = _outs_size_is_fixed = False
    _N_ins = _N_outs = 0
    
    def __init__(self, ID=None, thermo=None, *, units):
        ins = []
        outs = []
        self.units = set(units)
        for s in streams_from_units(units):
            source = s._source
            sink = s._sink
            if source in units and sink not in units:
                outs.append(s.copy('.' + s.ID))
            elif sink in units and source not in units:
                ins.append(s.copy('.' + s.ID))
        super().__init__(ID, ins, outs, thermo)
    
    def _assert_compatible_property_package(self): pass
    
    @property
    def auxiliary_unit_names(self):
        return [i.ID for i in self.units]
    
    @property
    def auxiliary_units(self):
        return self.units
    
    def _run(self): pass
