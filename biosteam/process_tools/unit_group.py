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

__all__ = ('UnitGroup',)


# %% Unit group for generating results

class UnitGroup:
    """Create a UnitGroup object for generating biorefinery results."""
    __slots__ = ('name', 'units', 'metrics')
    
    def __init__(self, name, units, metrics=None):
        self.name = str(name) #: [str] Name for group
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
            units = 'GJ'
            self.metric(name, units, self.name, lambda: self.get_utility_duty(agent))
        elif basis == 'flow':
            units = 'MT'
            self.metric(name, units, self.name, lambda: self.get_utility_flow(agent))
        else:
            raise ValueError(f"basis must be either 'duty' or 'flow', not {repr(basis)}")
    
    @property
    def heat_utilities(self):
        """[tuple] All HeatUtility objects."""
        return get_heat_utilities(self.units)
    
    @property
    def power_utilities(self):
        """[tuple] All PowerUtility objects."""
        return tuple(utils.get_power_utilities(self.units))
    
    @classmethod
    def group_by_type(cls, name, units, types):
        """Create a UnitGroup object of given type(s)."""
        isa = isinstance
        return cls(name, [i for i in units if isa(i, types)])
    
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
    
    def to_dict(self, with_electricity_production=False):
        """Return dictionary of results."""
        dct = {'Installed equipment cost [million USD]': self.get_installed_cost(),
               'Cooling duty [GJ/hr]': self.get_cooling_duty(),
               'Heating duty [GJ/hr]': self.get_heating_duty(),
               'Electricity consumption [MW]': self.get_electricity_consumption()}
        if with_electricity_production:
            dct['Electricity production [MW]'] = self.get_electricity_production()
        for i in self.metrics:
            dct[i.name_with_units] = i()
        return dct
                
    def to_series(self, with_electricity_production=False):
        """Return a pandas.Series object of results."""
        return pd.Series(self.to_dict(with_electricity_production), name=self.name)

    @classmethod
    def df_from_groups(cls, unit_groups, *, fraction=False):
        """Return a pandas.DataFrame object from unit groups."""
        data = [i.to_series() for i in unit_groups]
        df = pd.DataFrame(data)
        if fraction:
            values = df.values
            values *= 100 / values.sum(axis=0, keepdims=True)
        return df

    @classmethod
    def plot_bars_from_groups(cls, unit_groups, *, fraction=False, colormap='tab10', **kwargs):
        """Plot unit groups as a stacked bar chart."""
        df = cls.create_df(unit_groups, fraction)
        df.T.plot(kind='bar', stacked=True, 
                  colormap=colormap, 
                  **kwargs)
        plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
        plt.xticks(rotation=0)
        if fraction:
            plt.ylabel('[%]')
            plt.ylim(0, 100)

    def __repr__(self):
        return f"{type(self).__name__}({repr(self.name)}, {self.units}, metrics={self.metrics})"