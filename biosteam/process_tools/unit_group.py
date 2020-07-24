# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

import pandas as pd
from matplotlib import pyplot as plt
from . import utils
from ..evaluation import Metric
from .._heat_utility import HeatUtility
from collections.abc import Mapping
from ..plots import style_axis

__all__ = ('UnitGroup',)

INST_COST = 'Inst. cost'
ELEC_CONSUMPTION = 'Elec. consumption'
ELEC_PRODUCTION = 'Elec. production'
INSTALLED_EQUIPMENT_COST = 'Installed equipment cost'
COOLING_DUTY = 'Cooling duty'
HEATING_DUTY = 'Heating duty'
ELECTRICITY_CONSUMPTION = 'Electricity consumption'
ELECTRICITY_PRODUCTION = 'Electricity production'
CAPITAL_UNITS = '[MM$]'
ELEC_UNITS = '[MW]'
DUTY_UNITS = '[GJ/hr]'

# %% Unit group for generating results

class UnitGroup:
    """
    Create a UnitGroup object for generating biorefinery results.
    
    Parameters
    ----------
    name : str
        Name of group for bookkeeping.
    units : tuple[Unit]
        Unit operations.
    metrics=None : list[Metric], optional
        Metrics to generate results. These metrics are computed when 
        generating results as dictionaries, pandas series and data frames,
        and plots.
    
    Examples
    --------
    Create a UnitGroup from BioSTEAM's example ethanol subsystem:
    
    >>> from biosteam.examples import ethanol_subsystem_example
    >>> ethanol_sys = ethanol_subsystem_example()
    >>> group = UnitGroup('Ethanol production', ethanol_sys.units)
    
    You can get main process results using UnitGroup methods:
        
    >>> group.to_dict(with_electricity_production=True)
    {'Installed equipment cost [MM$]': 13.57808163192647, 'Cooling duty [GJ/hr]': 104.4639430479312, 'Heating duty [GJ/hr]': 156.88073343258333, 'Electricity consumption [MW]': 0.40116044433570186, 'Electricity production [MW]': 0.0}
    
    Each result can be retrieved separately:
    
    >>> group.get_installed_cost()
    13.57808163192647
    
    >>> group.get_heating_duty()
    156.88073343258333
    
    """
    __slots__ = ('name', 'units', 'metrics')
    
    def __init__(self, name, units, metrics=None):
        self.name = str(name) #: [str] Name of group for bookkeeping
        self.units = tuple(units) #: tuple[Unit] Unit operations
        self.metrics = metrics or [] #: list[Metric] Metrics to generate results
    
    def metric(self, name, units=None, getter=None):
        """Add metric to UnitGroup."""
        if not getter: return lambda getter: self.metric(name, units, getter)
        metric = Metric(name, getter, units, self.name)
        self.metrics.append(metric)
        return metric 
    
    def register_utility_agent(self, agent, basis='duty'):
        """Register utility agent as a metric to UnitGroup."""
        if isinstance(agent, str): agent = HeatUtility.get_agent(agent)
        name = agent.ID.replace('_', ' ').capitalize()
        if basis == 'duty':
            self.metric(name, 'GJ', lambda: self.get_utility_duty(agent))
        elif basis == 'flow':
            self.metric(name, 'MT', lambda: self.get_utility_flow(agent))
        else:
            raise ValueError(f"basis must be either 'duty' or 'flow', not {repr(basis)}")
    
    @property
    def heat_utilities(self):
        """[tuple] All HeatUtility objects."""
        return utils.get_heat_utilities(self.units)
    
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
        """Create a list of UnitGroup object of given iterable of name-type(s) pairs."""
        if name_types:
            if isinstance(name_types, Mapping): name_types = name_types.items()
            return [cls.filter_by_types(name, units, types) for name, types in name_types]
        else:
            groups = utils.group_by_types(units)
            return [cls(i, j) for i, j in groups.items()]
    
    @classmethod
    def group_by_lines(cls, units, name_lines=None):
        """Create a list of UnitGroup object of given iterable of name-type(s) pairs."""
        if name_lines:
            if isinstance(name_lines, Mapping): name_lines = name_lines.items()
            return [cls.filter_by_lines(name, units, lines) for name, lines in name_lines]
        else:
            groups = utils.group_by_lines(units)
            return [cls(i, j) for i, j in groups.items()]
    
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
    
    def get_electricity_consumption(self):
        """Return the total electricity consumption in MW."""
        return utils.get_electricity_consumption(self.power_utilities)

    def get_electricity_production(self):
        """Return the total electricity production in MW."""
        return utils.get_electricity_production(self.power_utilities)
    
    def to_dict(self, with_electricity_production=False, shorthand=False, with_units=True):
        """Return dictionary of results."""
        if shorthand:
            inst_cost = INST_COST
            elec_consumption = ELEC_CONSUMPTION
            elec_production = ELEC_PRODUCTION
        else:
            inst_cost = INSTALLED_EQUIPMENT_COST
            
            elec_consumption = ELECTRICITY_CONSUMPTION
            elec_production = ELECTRICITY_PRODUCTION
        cooling_duty = COOLING_DUTY
        heating_duty = HEATING_DUTY
        if with_units:
            inst_cost += ' ' + CAPITAL_UNITS
            elec_consumption += ' ' + ELEC_UNITS
            elec_production += ' ' + ELEC_UNITS
            cooling_duty += ' ' + DUTY_UNITS
            heating_duty += ' ' + DUTY_UNITS
        dct = {inst_cost: self.get_installed_cost(),
               cooling_duty: self.get_cooling_duty(),
               heating_duty: self.get_heating_duty(),
               elec_consumption: self.get_electricity_consumption()}
        if with_electricity_production:
            dct[elec_production] = self.get_electricity_production()
        for i in self.metrics:
            dct[i.name_with_units] = i()
        return dct
                
    def to_series(self, with_electricity_production=False, shorthand=False, with_units=True):
        """Return a pandas.Series object of results."""
        return pd.Series(self.to_dict(with_electricity_production, shorthand, with_units), name=self.name)

    @classmethod
    def df_from_groups(cls, unit_groups,
                       with_electricity_production=False, 
                       shorthand=False, fraction=False):
        """Return a pandas.DataFrame object from unit groups."""
        with_units = not fraction
        data = [i.to_series(with_electricity_production, shorthand, with_units) for i in unit_groups]
        df = pd.DataFrame(data)
        if fraction:
            values = df.values
            values *= 100 / values.sum(axis=0, keepdims=True)
        return df

    @classmethod
    def plot_bars_from_groups(cls, unit_groups, with_electricity_production=False,
                              shorthand=True, fraction=True, horizontal_ticks=False, 
                              **kwargs):
        """Plot unit groups as a stacked bar chart."""
        df = cls.df_from_groups(unit_groups,
                                with_electricity_production,
                                shorthand,
                                fraction)
        df.T.plot(kind='bar', stacked=True, **kwargs)
        locs, labels = plt.xticks()
        plt.xticks(locs, ['\n['.join(i.get_text().split(' [')) for i in labels])
        plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
        if horizontal_ticks: plt.xticks(rotation=0)
        if fraction:
            plt.ylabel('[%]')
            plt.ylim(0, 100)
        style_axis(top=False)

    def __repr__(self):
        return f"{type(self).__name__}({repr(self.name)}, {self.units}, metrics={self.metrics})"