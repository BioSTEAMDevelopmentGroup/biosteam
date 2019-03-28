# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 18:49:26 2019

@author: yoelr
"""
from biosteam.proxy_stream import ProxyStream
from biosteam.unit import Unit, metaUnit
from biosteam.utils import function
from biosteam.report.table import cost_table, heat_utilities_table, power_utilities_table, stream_table, save_system_results, results_table
from copy import copy


class metaProcess(metaUnit):
    
    def __new__(mcl, clsname, superclasses, new_definitions):
        return type.__new__(mcl, clsname, superclasses, new_definitions)

    def __repr__(cls):
        return type.__repr__(cls)
        
class Process(Unit, metaclass=metaProcess):
    """Create Process object as described by a system."""
    _utility_summary = []
    
    def __init__(self, system, ID=None, feeds=None, products=None):
        self._system = system
        self._units = system.units.union(system.offsite_units)
        self.kwargs = None
        super().__init__(ID or system.ID,
                        (ProxyStream(i.ID, i) for i in products),
                        (ProxyStream(i.ID, i) for i in feeds))
        for u in system._flattened_network:
            u._setup_linked_streams()
    
    @property
    def ins(self):
        return tuple(self._ins)
    
    @property
    def outs(self):
        return tuple(self._outs)
    
    def _install(self):
        self._utils = utils = []
        for u in self._units:
            utils.extend(u.heat_utilities)
            if u._has_power_utility:
                utils.append(u.power_utility)
    
    def _setup(self):
        self._check_kwargs = False
        sys = self._system
        for u in sys.units:
            if u._kwargs != u.kwargs:
                u._setup()
                u._kwargs = copy(u.kwargs)
    
    def _purchase_costs(self):
        s = sum
        for u in self._units:
            if u._has_cost: yield s(u._purchase_costs())
    
    def _run(self):
        if self._check_kwargs: self._setup()
        sys = self._system
        sys._reset_iter()
        sys._converge()
        
    def _simulate(self):
        self._system.simulate()
        
    def _summary(self):
        sys = self._system
        for u in sys.units: u._summary()
        for i in sys.facilities:
            if isinstance(i, function): i()
            else: i.simulate()
        Summary = self.results['Summary']
        Summary['Utility cost'] = sum(i.cost for i in self._utils)
        Summary['Purchase cost'] = sum(self._purchase_costs())
        self._check_kwargs = True
    
    def equipment_table(self):
        return results_table(self._units)
    
    # def save_excel(self, file='report.xlsx', **streamtable_properties):
    #     save_system_results(self._system, file)
    
    # def cost_table(self):
    #     return cost_table(self._system)
    
    # def stream_table(self, **properties):
    #     return stream_table(self._system.streams, **properties)
    
    def power_utilities_table(self):
        return power_utilities_table(self._units)
    
    def heat_utilities_table(self):
        return heat_utilities_table(self._units)