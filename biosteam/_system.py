# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>,
#                          Joy Zhang <joycheung1994@gmail.com>,
#                          Yalin Li <mailto.yalin.li@gmail.com>
#                          Sarang Bhagwat <sarangb2@illinois.edu>,
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from __future__ import annotations
from typing import Optional, Callable, Iterable, Sequence, Collection, TYPE_CHECKING
import flexsolve as flx
from .digraph import (digraph_from_units,
                      digraph_from_system,
                      minimal_digraph,
                      surface_digraph,
                      finalize_digraph)
from thermosteam import (
    AbstractStream, Stream, MultiStream, Chemical, 
    BipartitePhenomeGraph, PhenomeGraph, PhenomeNode, DotFile,
    all_subgraphs, 
)
from thermosteam.base import SparseArray
from . import HeatUtility, PowerUtility
from thermosteam.utils import registered
from scipy.optimize import root
from .exceptions import try_method_with_object_stamp, Converged, UnitInheritanceError
from thermosteam import Network, mark_disjunction, unmark_disjunction
from ._facility import Facility
from ._unit import Unit
from thermosteam.network import repr_ins_and_outs
from . import utils
from .utils import (
    repr_items, ignore_docking_warnings,
    piping, colors, list_available_names, dictionaries2array,
    Timer
)
from .process_tools import get_power_utilities, get_heat_utilities
from collections import abc
from warnings import warn
from inspect import signature
from thermosteam.utils import repr_kwargs
import biosteam as bst
import numpy as np
import pandas as pd
from numpy.linalg import solve
from scipy.integrate import solve_ivp
from . import report
from thermosteam.network import temporary_units_dump, TemporaryUnit
import os
import openpyxl
import thermosteam as tmo
from itertools import product
if TYPE_CHECKING: 
    from ._tea import TEA
    from .evaluation import Response

__all__ = ('System', 'AgileSystem', 'MockSystem',
           'AgileSystem', 'OperationModeResults',
           'mark_disjunction', 'unmark_disjunction')

# %% Miscillaneous

def _reformat(name):
    name = name.replace('_', ' ')
    if name.islower(): name= name.capitalize()
    return name

# %% System specification

class SystemSpecification:
    __slots__ = ('f', 'args')
    
    def __init__(self, f, args):
        self.f = f
        self.args = args
        
    def __call__(self): self.f(*self.args)


# %% Reconfiguration for phenomena oriented simulation and convergence

ObjectType = np.dtype('O')
TRACK_CONDITION_NUMBER = 2
TRACK_VARIABLES = 1

class Configuration:
    __slots__ = ('path', 'stages', 'nodes', 'streams',
                 'stream_ref', 'connections', 'aggregated',
                 'composition_sensitive_path', 
                 'composition_sensitive_nodes',
                 '_has_dynamic_coefficients')
    
    def __init__(self, path, stages, nodes, streams, stream_ref, connections, aggregated):
        self.path = path
        self.stages = stages
        self.nodes = nodes
        self.streams = streams
        self.stream_ref = stream_ref
        self.connections = connections
        self.aggregated = aggregated
        self.composition_sensitive_path = [i for i in path if getattr(i, 'composition_sensitive', False)]
        self.composition_sensitive_nodes =  [i for i in nodes if getattr(i, 'composition_sensitive', False) or isinstance(i, Stream)]
    
    def get_composition_sensitive_flows(self):
        streams = []
        past_streams = set()
        for i in self.composition_sensitive_nodes:
            for stream in (i.ins + i.outs):
                if stream in past_streams: continue
                past_streams.add(stream)
                streams.append(stream)
        phases = [i.phases for i in streams]
        flows = [stream.imol.data.copy() for stream in streams]
        return streams, phases, flows
    
    def get_composition_sensitive_data(self):
        streams = []
        past_streams = set()
        for i in self.composition_sensitive_nodes:
            for stream in (i.ins + i.outs):
                if stream in past_streams: continue
                past_streams.add(stream)
                streams.append(stream)
        data = [i.get_data() for i in streams]
        return streams, data
    
    def solve_nonlinearities(self):
        if self.aggregated:
            for i in self.path: 
                if hasattr(i, '_update_aggretated_nonlinearities'):
                    i._update_aggretated_nonlinearities()
                else:
                    i._update_nonlinearities()
        else:
            for i in self.path: 
                i._update_nonlinearities()
        
        if self.composition_sensitive_path:
            # for i in self.composition_sensitive_path: i._run()
            # breakpoint()
            containers, flows = self.solve_material_flows(composition_sensitive=True, full_output=True)
            def update_inner_material_balance_parameters(flows):
                for i in self.composition_sensitive_nodes: i._update_composition_parameters()
                return self.solve_material_flows(composition_sensitive=True)
            
            flx.fixed_point(
                update_inner_material_balance_parameters,
                flows, 
                xtol=tmo.LLE.pseudo_equilibrium_inner_loop_options['xtol'], 
                rtol=tmo.LLE.pseudo_equilibrium_inner_loop_options['rtol'], 
                maxiter=100, checkiter=False,
                checkconvergence=False,
            )
            for i in self.composition_sensitive_path: i._update_net_flow_parameters()
            self.solve_material_flows(composition_sensitive=True)
    
    def dynamic_coefficients(self, b):
        try:
            if not self._has_dynamic_coefficients: return None
        except:
            delayed = [(i, j) for i, j in enumerate(b) if j.dtype is ObjectType]
            self._has_dynamic_coefficients = bool(delayed)
        else:
            delayed = [(i, j) for i, j in enumerate(b) if j.dtype is ObjectType]
        return delayed
    
    def solve_material_flows(self, composition_sensitive=False, full_output=False):
        if composition_sensitive:
            nodes = self.composition_sensitive_nodes
        else:
            nodes = self.nodes
        A = []
        b = []
        for node in nodes:
            f = node._create_material_balance_equations
            try: eqs = f(composition_sensitive)
            except TypeError as e: 
                try: eqs = f()
                except: raise e from None
            for coefficients, value in eqs:
                coefficients = {
                    i.material_reference: j for i, j in coefficients.items()
                }
                A.append(coefficients)
                b.append(value)
        delayed = self.dynamic_coefficients(b)
        if delayed:
            raise NotImplementedError('delayed coefficients')
            delayed_index = []
            for _, t in delayed:
                delayed_index.extend([i for i, j in enumerate(t) if callable(j)])
            delayed_set = set(delayed_index)
            delayed_index = list(delayed_set)
            ready_index = [i for i in range(len(t)) if i not in delayed_set]
            A_ready = [{i: j[ready_index] for i, j in dct.items()} for dct in A]
            A_ready, objs = dictionaries2array(A_ready)
            b_ready = np.array([i[ready_index] for i in b], float)
            values = solve(A_ready, b_ready.T).T
            values[values < 0] = 0
            for obj, value in zip(objs, values): # update material flows
                indexer, phase = obj
                index = indexer._phase_indexer(phase)
                mol = indexer.data.rows[index]
                mol[ready_index] = value
            A_delayed = [{i: j[delayed_index] for i, j in dct.items()} for dct in A]
            A_delayed, objs = dictionaries2array(A_delayed)
            b_delayed = np.array([
                [(f(j) if callable(f:=bi[j]) else f) for j in delayed_index]
                for bi in b
            ], float)
            values = solve(A_delayed, b_delayed.T).T
            values[values < 0] = 0
            for obj, value in zip(objs, values): 
                obj._update_material_flows(value, delayed_index)
        else:
            A, objs = dictionaries2array(A)
            if np.isnan(A).any(): raise RuntimeError('invalid number encountered')
            values = solve(A, np.array(b).T).T
            if np.isnan(values).any(): raise RuntimeError('invalid number encountered')
            values[(values < 0) & (values > -1e-6)] = 0
            masks = values < 0
            if masks.any():
                containers = []
                for obj in objs:
                    indexer, phase = obj
                    index = indexer._phase_indexer(phase)
                    containers.append(indexer.data.rows[index])
                new_values = values[masks] 
                old_values = SparseArray.from_rows(containers).to_array()
                denominator = (old_values[masks] - new_values)
                denominator[denominator == 0] = 1
                weights = -new_values / denominator
                relaxation_factor = weights.max()
                values = relaxation_factor * old_values + (1 - relaxation_factor) * values
                values[values < 0] = 0
                for mol, value in zip(containers, values): # update material flows
                    mol[:] = value
            else:
                if full_output:
                    containers = []
                    for obj, value in zip(objs, values): # update material flows
                        indexer, phase = obj
                        index = indexer._phase_indexer(phase)
                        mol = indexer.data.rows[index]
                        containers.append(mol)
                        mol[:] = value
                else:
                    for obj, value in zip(objs, values): # update material flows
                        indexer, phase = obj
                        index = indexer._phase_indexer(phase)
                        mol = indexer.data.rows[index]
                        mol[:] = value
        for node in nodes: 
            if hasattr(node, '_update_auxiliaries'): node._update_auxiliaries()
        if full_output:
            return containers, values
        else:
            return values
        
    def solve_energy_flows(self):
        nodes = self.nodes
        A = []
        b = []
        for node in nodes:
            for coefficients, value in node._create_energy_balance_equations():
                A.append(coefficients)
                b.append(value)
        if not A: return
        for node in nodes:
            for coefficients, value in node._create_bulk_balance_equations():
                A.append(coefficients)
                b.append(value)
        A = [{(getattr(i, 'material_reference', i), v): j for (i, v), j in x.items()} for x in A]
        A, objs = dictionaries2array(A)
        values = solve(A, np.array(b).T).T
        # try: values = solve(A, np.array(b).T).T
        # except: 
        #     A = []
        #     b = []
        #     breakpoint()
        #     for node in nodes:
        #         for coefficients, value in node._create_energy_balance_equations():
        #             breakpoint()
        #             A.append(coefficients)
        #             b.append(value)
        #     if not A: return
        #     for node in nodes:
        #         for coefficients, value in node._create_bulk_balance_equations():
        #             breakpoint()
        #             A.append(coefficients)
        #             b.append(value)
        for (obj, var), value in zip(objs, values):
            if var == 'F_mol':
                indexer, phase = obj
                index = indexer._phase_indexer(phase)
                mol = indexer.data.rows[index]
                if value == 0:
                    mol[:] = 0
                else:
                    total = mol.sum()
                    if total != 0: mol[:] *= value / total
            else:
                obj._update_variable(var, value)
        for node in nodes:
            if hasattr(node, '_reset_bulk_variable'): node._reset_bulk_variable()
        
    def __enter__(self):
        units = self.stages
        streams = self.streams
        sinks = {}
        sources = {}
        for u in units:
            ins = u.flat_ins
            outs = u.flat_outs
            for i in ins: 
                i._sink = u
                sinks[i.imol] = u
            for i in outs: 
                i._source = u
                sources[i.imol] = u
        units_set = set(units)
        for i in streams:
            if i.source and i.source not in units_set and i.imol in sources: i._source = sources[i.imol]
            if i.sink and i.sink not in units_set and i.imol in sinks: i._sink = sinks[i.imol]
        return self
    
    def __exit__(self, type, exception, traceback):
        for s, source, sink in self.connections:
            s._source = source
            s._sink = sink
        if exception: raise exception
            
class MaterialFlows:
    __slots__ = ('containers', 'old_values', 'values')
    
    def __init__(self, containers, old_values, values):
        self.containers = containers
        self.old_values = old_values
        self.values = values


# %% Sequential modular convergence tools

class Methods(dict):
    
    def __repr__(self):
        return f"{type(self).__name__}({', '.join(self)})"


class RecycleData:
    __slots__ = (
        'recycle',
        'mol',
        'T',
        'P',
    )
    def __init__(self, recycle):
        self.recycle = recycle
        self.mol = recycle.mol.copy()
        self.T, self.P = recycle.thermal_condition
        
    def to_dict(self, dct=None):
        recycle = self.recycle
        IDs = recycle.chemicals.IDs
        if dct is None: dct = {}
        for i, j in self.mol.nonzero_items(): dct[recycle, IDs[i]] = j
        dct[recycle, 'T'] = self.T
        dct[recycle, 'P'] = self.P
        return dct
        
    def to_series(self):
        return pd.Series(self.to_list(), self.get_names())
    
    def get_keys(self, lst=None):
        recycle = self.recycle
        IDs = recycle.chemicals.IDs
        if lst is None: lst = []
        for i, j in self.mol.nonzero_items(): 
            lst.append(
                (recycle, IDs[i])
            )
        lst.append(
            (recycle, 'T')
        )
        lst.append(
            (recycle, 'P')
        )
        return lst
    
    def get_names(self, lst=None):
        recycle = self.recycle
        ID = recycle.ID
        IDs = recycle.chemicals.IDs
        if lst is None: lst = []
        for i, j in self.mol.nonzero_items(): 
            lst.append(
                f"{ID}.{IDs[i]}"
            )
        lst.append(
            f"{ID}.T"
        )
        lst.append(
            f"{ID}.P"
        )
        return lst
        
    def to_tuples(self, lst=None):
        recycle = self.recycle
        IDs = recycle.chemicals.IDs
        if lst is None: lst = []
        for i, j in self.mol.nonzero_items(): 
            lst.append(
                ((recycle, IDs[i]), j)
            )
        lst.append(
            ((recycle, 'T'), self.T)
        )
        lst.append(
            ((recycle, 'P'), self.P)
        )
        return lst
    
    def to_list(self, lst=None):
        if lst is None:
            lst = [*self.mol.nonzero_values(), self.T, self.P]
        else:
            lst.extend(self.mol.nonzero_values())
            lst.append(self.T)
            lst.append(self.P)
        return lst
    
    def to_array(self):
        return np.array(self.to_list())
    
    def reset(self):
        recycle = self.recycle
        recycle.mol.copy_like(self.mol)
        recycle.T = self.T
        recycle.P = self.P
        
    def update(self):
        recycle = self.recycle
        self.mol = recycle.mol.copy()
        self.T, self.P = recycle.thermal_condition
        
    def __repr__(self):
        return f"{type(self).__name__}({self.recycle})"
        
    
class JointRecycleData:
    __slots__ = ('recycle_data', 'responses')
    def __init__(self, recycles, responses):
        self.recycle_data = [RecycleData(i) for i in recycles]
        self.responses = {i: i.get() for i in responses}
    
    def get_keys(self):
        lst = []
        for i in self.recycle_data: i.get_keys(lst)
        lst.extend(self.responses)
        return lst
    
    def get_names(self):
        lst = []
        for i in self.recycle_data: i.get_names(lst)
        lst.extend([str(i) for i in self.responses])
        return lst
    
    def to_dict(self):
        dct = {}
        for i in self.recycle_data: i.to_dict(dct)
        for i, j in self.responses.items(): dct[i] = j
        return dct
    
    def to_series(self):
        return pd.Series(self.to_list(), self.get_names())
    
    def to_tuples(self):
        lst = []
        for i in self.recycle_data: i.to_tuples(lst)
        for i, j in self.responses.items(): lst.append((i, j))
        return lst
    
    def to_list(self):
        lst = []
        for i in self.recycle_data: i.to_list(lst)
        for i in self.responses.values(): lst.append(i)
        return lst
    
    def to_array(self):
        return np.array(self.to_list())
    
    def reset(self):
        for i in self.recycle_data: i.reset()
        
    def update(self):
        for i in self.recycle_data: i.update()

    def __repr__(self):
        recycles = ', '.join([str(i.recycle) for i in self.recycle_data])
        responses =  ', '.join([str(i) for i in self.responses])
        return f"{type(self).__name__}(recycles=[{recycles}], responses=[{responses}])"


def get_recycle_data(stream):
    """
    Return stream temperature, pressure, and molar flow rates as a
    1d array.
    
    """
    data = stream.imol.data
    size = data.size
    TP = stream._thermal_condition
    arr = np.zeros(size + 2)
    arr[0] = TP.T
    arr[1] = TP.P
    data.to_flat_array(arr[2:])
    return arr

def set_recycle_data(stream, arr):
    """
    Set stream temperature, pressure, and molar flow rates with a
    1d array.
    
    """
    data = stream.imol.data
    data.from_flat_array(arr[2:])
    TP = stream._thermal_condition
    T = float(arr[0]) # ndfloat objects are slow and parasitic (don't go away)
    P = float(arr[1])
    TP._T =  TP._T * 0.5 if T == 0 else T
    TP._P = TP._P * 0.5 if P == 0. else P
    
   
# %% System creation tools

def get_units(sys, units=None):
    if units is None: units = []
    isa = isinstance
    for i in sys._path:
        if isa(i, Unit):
            units.append(i)
        elif isa(i, System):
            get_units(i, units)
    return units

def get_missing_units(facility, units, path):
    isa = isinstance
    path.append(facility)
    for i in facility.outs:
        sink = i.sink
        if sink is None or isa(sink, Facility) or sink in units: continue
        units.add(sink)
        get_missing_units(sink, units, path)

def facilities_from_units(units):
    isa = isinstance
    return [i for i in units if isa(i, Facility)]

def find_blowdown_recycle(facilities):
    isa = isinstance
    for i in facilities:
        if isa(i, bst.BlowdownMixer): return i.outs[0]

def find_recycles(path):
    all_units = set(path)
    past_outlets = set()
    recycles = set()
    for i in path:
        for s in i.ins:
            if s in past_outlets or s.source not in all_units: continue
            recycles.add(s)
        past_outlets.update(i.outs)
    return list(recycles)

# %% LCA

class ProcessImpactItem:
    __slots__ = ('name', 'basis', 'inventory', 'CF')
    
    def __init__(self, name, basis, inventory, CF):
        self.name = name
        self.basis = basis
        self.inventory = inventory
        self.CF = CF
        
    def impact(self):
        return self.inventory() * self.CF
    
def set_impact_value(properties, name, value, groups):
    for group in groups:
        if group in name:
            group_name = _reformat(group)
            if group_name in properties:
                properties[group_name] += value
            elif value:
                properties[group_name] = value
            break
    else:
        if value: properties[_reformat(name)] = value
    

# %% Debugging and exception handling

def raise_recycle_type_error(recycle):
    raise ValueError(
       f"invalid recycle of type '{type(recycle).__name__}' encountered; "
        "recycle must be either a Stream object, a tuple of Stream objects, or None"
    )

def print_exception_in_debugger(self, func, e):
    print(f"{colors.exception(type(e).__name__+ ':')} {e}")
    try: self.show()
    except: pass

def update_locals_with_flowsheet(lcs):
    lcs.update(bst.main_flowsheet.to_dict())
    lcs.update(bst.__dict__)

def _method_debug(self, f):
    """Method decorator for debugging system."""
    def g(*args, **kwargs):
        try:
            f(*args, **kwargs)
        except Exception as e:
            print_exception_in_debugger(self, f, e)
            update_locals_with_flowsheet(locals())
            # All systems, units, streams, and flowsheets are available as
            # local variables. Although this debugging method is meant
            # for internal development, please feel free to give it a shot.
            breakpoint()
    g.__name__ = f.__name__
    g.__doc__ = f.__doc__
    g._original = f
    return g

def _method_profile(self, f):
    self._total_excecution_time_ = 0.
    t = utils.TicToc()
    def g():
        t.tic()
        f()
        self._total_excecution_time_ += t.elapsed_time
    g.__name__ = f.__name__
    g.__doc__ = f.__doc__
    g._original = f
    return g

# %% Converging recycle systems

class MockSystem:
    """
    Create a MockSystem object with inlets and outlets just like System
    objects, but without implementing any of the convergence methods nor
    path related attributes.

    Parameters
    ----------
    units : Iterable[:class:`~biosteam.Unit`], optional
        Unit operations in mock system.

    Notes
    -----
    This object is used to prevent the creation of unneeded systems for less
    computational effort.

    """
    __slots__ = ('units',
                 'flowsheet',
                 '_ins',
                 '_outs',
                 '_inlet_names',
                 '_outlet_names')

    def __init__(self, units=()):
        self.units = units or list(units)
        self._load_flowsheet()

    def _load_flowsheet(self):
        self.flowsheet = bst.main_flowsheet.get_flowsheet()

    @property
    def ins(self) -> piping.StreamPorts[piping.InletPort]:
        """All inlets to the system."""
        try:
            return self._ins
        except:
            inlets = utils.feeds_from_units(self.units)
            self._ins = ins = piping.StreamPorts.from_inlets(inlets, sort=True)
            return ins
    @property
    def outs(self) -> piping.StreamPorts[piping.InletPort]:
        """All outlets to the system."""
        try:
            return self._outs
        except:
            outlets = utils.products_from_units(self.units)
            self._outs = outs = piping.StreamPorts.from_outlets(outlets, sort=True)
            return outs

    def load_inlet_ports(self, inlets, names=None):
        """Load inlet ports to system."""
        all_inlets = utils.feeds_from_units(self.units)
        inlets = list(inlets)
        for i in inlets:
            if i not in all_inlets:
                raise ValueError(f'{i} is not an inlet')
        self._ins = piping.StreamPorts.from_inlets(inlets)
        self._inlet_names = {} if names is None else names
    
    def load_outlet_ports(self, outlets, names=None):
        """Load outlet ports to system."""
        all_outlets = utils.products_from_units(self.units)
        outlets = list(outlets)
        for i in outlets:
            if i not in all_outlets:
                raise ValueError(f'{i} is not an outlet')
        self._outs = piping.StreamPorts.from_outlets(outlets)
        self._outlet_names = {} if names is None else names

    def get_inlet(self, name):
        if name in self._inlet_names: 
            return self._ins[self._inlet_names[name]]
        raise ValueError(f"inlet named {repr(name)} does not exist")

    def get_outlet(self, name):
        if name in self._outlet_names: 
            return self._outs[self._outlet_names[name]]
        raise ValueError(f"outlet named {repr(name)} does not exist")

    def set_inlet(self, ID, inlet):
        stream = self.get_inlet(ID)
        sink = stream.sink
        if sink: sink.ins.replace(stream, inlet)
    
    def set_outlet(self, ID, outlet):
        stream = self.get_outlet(ID)
        source = stream.sink
        if source: source.outs.replace(stream, outlet)

    def __enter__(self):
        if self.units:
            raise RuntimeError("only empty mock systems can enter `with` statement")
        unit_registry = self.flowsheet.unit
        unit_registry.open_context_level()
        return self

    def __exit__(self, type, exception, traceback):
        if self.units:
            raise RuntimeError('mock system was modified before exiting `with` statement')
        unit_registry = self.flowsheet.unit
        dump = unit_registry.close_context_level()
        self.units = dump
        if exception: raise exception

    __sub__ = Unit.__sub__
    __rsub__ = Unit.__rsub__
    __pow__ = __sub__
    __rpow__ = __rsub__

    def show(self):
        ins = repr_items('    ins=', self.ins._ports, brackets='[]')
        outs = repr_items('    outs=', self.outs._ports, brackets='[]')
        units = repr_items('    units=', self.units, brackets='[]')
        args = ',\n'.join([ins, outs, units])
        print(f"{type(self).__name__}(\n{args}\n)")

    _ipython_display_ = show

    def __repr__(self):
        return f"{type(self).__name__}(ins={self.ins}, outs={self.outs})"


@registered('SYS')
class System:
    """
    Create a System object that can iteratively run each element in a path
    of BioSTREAM objects until the recycle stream is converged. A path can
    have :class:`~biosteam.Unit` and/or :class:`~biosteam.System` objects.
    When the path contains an inner System object, it converges/solves it in
    each loop/iteration.

    Parameters
    ----------
    ID :
        Unique identification. If ID is None, instance will not be
        registered in flowsheet.
    path :
        Path that is run element by element until the recycle converges.
    recycle : 
        Tear stream for the recycle loop.
    facilities : 
        Offsite facilities that are simulated only after
        completing the path simulation.
    N_runs : 
        Number of iterations to converge the system.
    operating_hours :
        Number of operating hours in a year. This parameter is used to
        compute annualized properties such as utility cost and material cost
        on a per year basis.
    lang_factor : 
        Lang factor for getting fixed capital investment from
        total purchase cost. If no lang factor, installed equipment costs are
        estimated using bare module factors.

    """
    __slots__ = (
        '_ID',
        '_path',
        '_facilities',
        '_recycle',
        '_N_runs',
        '_method',
        '_algorithm',
        '_mol_error',
        '_T_error',
        '_rmol_error',
        '_rT_error',
        '_iter',
        '_ins',
        '_outs',
        '_path_cache',
        '_prioritized_units',
        '_temporary_connections_log',
        '_integrated_facilities',
        '_shared_facilities',
        'maxtime',
        'maxiter',
        'molar_tolerance',
        'relative_molar_tolerance',
        'temperature_tolerance',
        'relative_temperature_tolerance',
        'operating_hours',
        'flowsheet',
        'lang_factor',
        'process_impact_items',
        'tracked_recycles',
        '_connections',
        '_TEA',
        '_LCA',
        '_subsystems',
        '_stages',
        '_aggregated_stages',
        '_units',
        '_unit_set',
        '_unit_path',
        '_cost_units',
        '_streams',
        '_feeds',
        '_products',
        '_inlet_names',
        '_outlet_names',
        # Specifications
        'simulate_after_specifications',
        '_specifications',
        '_running_specifications',
        '_simulation_default_arguments',
        '_simulation_outputs',
        # Convergence prediction
        '_responses',
        # Phenomena oriented simulation
        '_last_flows',
        '_last_error',
        '_diverged_count',
        '_aggregated_stage_configuration',
        'adaptive_phenomena_oriented_simulation',
        '_stage_configuration',
        'grouped_variables',
        'tracking_flag',
        'variable_profiles',
        'condition_number_profile',
        'equation_profiles',
        'edge_profiles',
        'variable_nodes',
        'equation_nodes',
        'timer',
        # Dynamic simulation
        '_isdynamic',
        '_state',
        '_state_idx',
        '_state_header',
        '_DAE',
        '_scope',
        'dynsim_kwargs',
    )

    take_place_of = Unit.take_place_of
    replace_with = Unit.replace_with

    ### Class attributes ###

    #: Default maximum number of iterations
    default_maxiter: int = 200

    #: Default molar tolerance for each component [kmol/hr]
    default_molar_tolerance: float = 1.

    #: Default relative molar tolerance for each component
    default_relative_molar_tolerance: float = 0.01

    #: Default temperature tolerance [K]
    default_temperature_tolerance: float = 0.10

    #: Default relative temperature tolerance
    default_relative_temperature_tolerance: float = 0.001

    #: Default convergence algorithm.
    default_algorithm: str = 'Sequential modular'

    #: Default method for convergence algorithm.
    default_methods: dict[str] = {
        'Sequential modular': 'Aitken',
        'Phenomena based': 'fixed-point',
        'Phenomena modular': 'fixed-point',
    }
    
    #: Whether to raise a RuntimeError when system doesn't converge
    strict_convergence: bool = True

    #: Method definitions for convergence
    available_methods: Methods[str, tuple[Callable, bool, dict]] = Methods()

    @classmethod
    def register_method(cls, name, solver, conditional=False, **kwargs):
        """
        Register new convergence method (root solver). Two solver signatures
        are supported:
        
        * If conditional is False, the signature must be solver(f, x, **kwargs) 
          where f(x) = 0 is the solution. This is common for scipy solvers.
        
        * If conditional is True, the signature must be solver(f, x, **kwargs) 
          where f(x) = (x, converged) is the solution and the solver stops 
          when converged is True. This method is preferred in BioSTEAM.
        
        """
        name = name.lower().replace('-', '').replace('_', '').replace(' ', '')
        cls.available_methods[name] = (solver, conditional, kwargs)

    @classmethod
    def from_feedstock(cls,
            ID: Optional[str]='', 
            feedstock: Stream=None, 
            feeds: Optional[Iterable[Stream]]=None, 
            facilities: Iterable[Facility]=(),
            ends: Iterable[Stream]=None,
            operating_hours: Optional[float]=None,
            **kwargs,
        ):
        """
        Create a System object from a feedstock.

        Parameters
        ----------
        ID : 
            Name of system.
        feedstock :
            Main feedstock of the process.
        feeds : 
            Additional feeds to the process.
        facilities : 
            Offsite facilities that are simulated only after
            completing the path simulation.
        ends : 
            Streams that not products, but are ultimately specified through
            process requirements and not by its unit source.
        operating_hours : 
            Number of operating hours in a year. This parameter is used to
            compute annualized properties such as utility cost and material cost
            on a per year basis.

        """
        if feedstock is None: raise ValueError('must pass feedstock stream')
        network = Network.from_feedstock(feedstock, feeds, ends)
        return cls._from_network(ID, network, facilities,
                                operating_hours,
                                **kwargs)

    @classmethod
    def from_units(cls, ID: Optional[str]="",
                   units: Optional[Iterable[Unit]]=None, 
                   ends: Optional[Iterable[Stream]]=None,
                   operating_hours: Optional[float]=None,
                   **kwargs):
        """
        Create a System object from all units given.

        Parameters
        ----------
        ID : 
            Name of system.
        units :
            Unit operations to be included.
        ends : 
            End streams of the system which are not products. Specify this
            argument if only a section of the complete system is wanted, or if
            recycle streams should be ignored.
        operating_hours : 
            Number of operating hours in a year. This parameter is used to
            compute annualized properties such as utility cost and material cost
            on a per year basis.

        """
        facilities = facilities_from_units(units)
        network = Network.from_units(units, ends)
        return cls._from_network(ID, network, facilities,
                                operating_hours,
                                **kwargs)

    @classmethod
    def from_segment(cls, ID: Optional[str]="", start: Optional[Unit]=None, 
                     end: Optional[Unit]=None, operating_hours: Optional[float]=None,
                     inclusive=False, **kwargs):
        """
        Create a System object from all units in between start and end.

        Parameters
        ----------
        ID : 
            Name of system.
        start :
            Only downstream units from start are included in the system.
        end : 
            Only upstream units from end are included in the system.
        operating_hours : 
            Number of operating hours in a year. This parameter is used to
            compute annualized properties such as utility cost and material cost
            on a per year basis.
        inclusive :
            Whether to include start and end units.

        """
        if start is None:
            if end is None: raise ValueError("must pass start and/or end")
            units = end.get_upstream_units(facilities=False)
        elif end is None:
            units = start.get_downstream_units(facilities=False)
        else:
            upstream_units = end.get_upstream_units(facilities=False)
            downstream_units = start.get_downstream_units(facilities=False)
            units = upstream_units.intersection(downstream_units)
        if inclusive:
            if start is not None: units.add(start)
            if end is not None: units.add(end)
        return bst.System.from_units(ID, units, operating_hours=operating_hours,
                                     **kwargs)
         
    @classmethod
    def _from_network(cls, ID, network, facilities=(),
                     operating_hours=None, **kwargs):
        """
        Create a System object from a network.

        Parameters
        ----------
        ID : str
            Name of system.
        network : Network
            Network that defines the simulation path.
        facilities : Iterable[Facility]
            Offsite facilities that are simulated only after
            completing the path simulation.
        operating_hours : float, optional
            Number of operating hours in a year. This parameter is used to
            compute annualized properties such as utility cost and material cost
            on a per year basis.

        """
        facilities = Facility.ordered_facilities(facilities)
        isa = isinstance
        ID_subsys = None if ID is None else ''
        path = [(cls._from_network(ID_subsys, i) if isa(i, Network) else i)
                for i in network.path]
        return cls(ID, path, network.recycle, facilities, None,
                   operating_hours, **kwargs)

    def __init__(self, 
            ID: Optional[str]= '', 
            path: Optional[Iterable[Unit|System]]=(), 
            recycle: Optional[Stream]=None, 
            facilities: Iterable[Facility]=(),
            N_runs: Optional[int]=None,
            operating_hours: Optional[float]=None,
            lang_factor: Optional[float]=None,
            responses: Optional[list[Response]]=None,
            algorithm: Optional[str]=None,
            method: Optional[str]=None,
            maxiter: Optional[int]=None,
            molar_tolerance: Optional[float]=None,
            relative_molar_tolerance: Optional[float]=None,
            temperature_tolerance: Optional[float]=None,
            relative_temperature_tolerance: Optional[float]=None,
            maxtime: Optional[float]=None,
            shared_facilities=None,
        ):
        self.N_runs = N_runs
        
        self.algorithm = self.default_algorithm if algorithm is None else algorithm
        
        self.method = self.default_methods[self.algorithm] if method is None else method
        
        #: Maximum number of iterations.
        self.maxiter: int = self.default_maxiter if maxiter is None else maxiter

        #: Maximum simulation time [s].
        self.maxtime: float|None = maxtime

        #: Molar tolerance [kmol/hr].
        self.molar_tolerance: float = self.default_molar_tolerance if molar_tolerance is None else molar_tolerance

        #: Relative molar tolerance.
        self.relative_molar_tolerance: float = self.default_relative_molar_tolerance if relative_molar_tolerance is None else relative_molar_tolerance

        #: Temperature tolerance [K].
        self.temperature_tolerance: float = self.default_temperature_tolerance if temperature_tolerance is None else temperature_tolerance

        #: Relative temperature tolerance.
        self.relative_temperature_tolerance: float = self.default_relative_temperature_tolerance if relative_temperature_tolerance is None else relative_temperature_tolerance

        #: Number of operating hours per year
        self.operating_hours: float|None = operating_hours if operating_hours is None else operating_hours
        
        #: Lang factor for computing fixed capital cost from purchase costs
        self.lang_factor: float|None = lang_factor

        #: Unit operations that have been integrated into the system configuration.
        self._prioritized_units = set()

        #: Cache for path segments and sections.
        self._path_cache = {}

        #: Log for all process specifications checked for temporary connections.
        self._temporary_connections_log = set()

        #: Whether subsystem facilities are shared with the entire system
        self._shared_facilities = True if shared_facilities is None else shared_facilities

        #: Whether to simulate system after running all specifications.
        self.simulate_after_specifications = False

        #: Whether to run sequential modular simulation when phenomena oriented simulation is not converging.
        self.adaptive_phenomena_oriented_simulation = True

        #: Whether facilities have a recycle loop.
        self._integrated_facilities = None

        #: Flag to track convergence at each iteration. 
        #: * 0 - no tracking
        #: * 1 - track variables
        #: * 2 - track variables and condition number
        self.tracking_flag = 0
        
        self._register(ID)
        self._set_path(path)
        self.recycle = recycle
        self._set_facilities(facilities)
        self._specifications = []
        self._running_specifications = False
        self._load_flowsheet()
        self._reset_errors()
        self._save_configuration()
        self._load_stream_links()
        self._state = None
        self._state_idx = None
        self._state_header = None
        self._DAE = None
        self.dynsim_kwargs = {}
        self.tracked_recycles = {}
        self._last_error = np.inf
        self._last_flows = 0
        self._diverged_count = 0
        subsystems = self.subsystems
        algorithm = self._algorithm
        method = self._method
        for i in subsystems: i._algorithm = algorithm
        for i in subsystems: i._method = method

    @property
    def responses(self):
        """Unit design decisions that need to converge to satisfy
        process specifications."""
        try:
            return self._responses
        except:
            self._responses = responses = []
            for unit in self.units: responses.extend(unit.responses)
            return responses

    def parameter(self, 
            setter: Optional[Callable]=None,
            element: Optional[Unit]=None, 
            coupled: Optional[bool]=None,
            name: Optional[str]=None, 
            distribution: Optional[str|object]=None, 
            units: Optional[str]=None, 
            baseline: Optional[float]=None, 
            bounds: Optional[tuple[float, float]]=None, 
            hook: Optional[Callable]=None, 
            description: Optional[str]=None, 
            optimized=False, 
            kind=None, # For backwards compatibility
        ):
        """
        Define and register parameter.
        
        Parameters
        ---------*    
        setter : 
            Should set parameter in the element.
        element : 
            Element in the system being altered.
        coupled : 
            Whether parameter is coupled to the system's mass and energy balances.
            This allows a ConvergenceModel to predict it's impact on recycle loops.
            Defaults to False.
        name : 
            Name of parameter. If None, default to argument name of setter.
        distribution : 
            Parameter distribution.
        units : 
            Parameter units of measure
        baseline : 
            Baseline value of parameter.
        bounds : 
            Lower and upper bounds of parameter.
        hook : 
            Should return the new parameter value given the sample.
        
        """
        if kind == 'coupled': coupled = True
        elif coupled is None: coupled = False
        if isinstance(setter, bst.Parameter):
            if element is None: element = setter.element
            if coupled is None: coupled = setter.coupled
            if name is None: name = setter.name
            if distribution is None: distribution = setter.distribution
            if units is None: units = setter.units
            if baseline is None: baseline = setter.baseline
            if bounds is None: bounds = setter.bounds
            if hook is None: hook = setter.hook
            if description is None: description = setter.description
            setter = setter.setter
        elif isinstance(setter, bst.MockFeature):
            if element is None: element = setter.element
            if name is None: name = setter.name
            if units is None: units = setter.units
        elif not setter:
            return lambda setter: self.parameter(setter, element, coupled, name,
                                                 distribution, units, baseline,
                                                 bounds, hook, description, optimized)
        p = bst.Parameter(name, setter, element,
                          self, distribution, units, 
                          baseline, bounds, coupled, hook, description)
        return p

    def track_recycle(self, recycle: Stream, collector: list[Stream]=None):
        if not isinstance(recycle, Stream):
            for i in recycle: self.track_recycle(recycle, collector)
        if collector is None: collector = []
        self.tracked_recycles[recycle] = collector
        unit = recycle.sink
        system = self.find_system(unit)
        if system is self: return
        system.track_recycle(recycle, collector)

    def update_configuration(self,
            units: Optional[Sequence[str]]=None,
        ):
        self._update_configuration(units)
        self._save_configuration()

    def _update_configuration(self,
            units: Optional[Sequence[str]]=None,
        ):
        old_IDs = [i.ID for i in self.subsystems]
        # Warning: This method does not save the configuration.
        if units is None: units = self.units
        self._delete_path_cache()
        isa = isinstance
        Facility = bst.Facility
        facilities = Facility.ordered_facilities([i for i in units if isa(i, Facility)])
        ID_subsys = None if '.' in self.ID else ''
        network = Network.from_units(units)
        path = [(type(self)._from_network(ID_subsys, i) if isa(i, Network) else i)
                for i in network.path]
        self._reset_errors()
        self._set_path(path)
        self.recycle = network.recycle
        self._set_facilities(facilities)
        self.set_tolerance(
            mol=self.molar_tolerance,
            rmol=self.relative_molar_tolerance,
            T=self.temperature_tolerance,
            rT=self.relative_temperature_tolerance,
            maxiter=self.maxiter,
            subsystems=True,
        )
        for i, j in zip(self.subsystems, old_IDs): i.ID = j

    def __enter__(self):
        if self._path or self._recycle or self._facilities:
            raise RuntimeError("only empty systems can enter `with` statement")
        del self._units
        unit_registry = self.flowsheet.unit
        unit_registry.open_context_level()
        return self

    def __exit__(self, type, exception, traceback):
        unit_registry = self.flowsheet.unit
        dump = unit_registry.close_context_level()
        if exception: raise exception
        if self._path or self._recycle or self._facilities:
            raise RuntimeError('system cannot be modified before exiting `with` statement')
        else:
            self.update_configuration(dump)
            self._load_stream_links()
            self.set_tolerance(
                algorithm=self._algorithm,
                method=self._method,
                mol=self.molar_tolerance,
                rmol=self.relative_molar_tolerance,
                T=self.temperature_tolerance,
                rT=self.relative_temperature_tolerance,
                subsystems=True,
            )

    def _save_configuration(self):
        self._connections = [i.get_connection() for i in self.streams]

    @ignore_docking_warnings
    def _load_configuration(self):
        for i in self._connections: i.reconnect()
        if self.recycle:
            for i in self.units: 
                i._system = i._recycle_system = self
                if hasattr(i, 'stages'):
                    for j in i.stages:
                        j._system = j._recycle_system = self
        else:
            for i in self.units: 
                i._system = self
                i._recycle_system = None
                if hasattr(i, 'stages'):
                    for j in i.stages:
                        j._system = self
                        j._recycle_system = None
            for sys in self.subsystems: 
                if sys.recycle: 
                    for i in sys.units: 
                        i._recycle_system = sys
                        if hasattr(i, 'stages'):
                            for j in i.stages: j._recycle_system = sys

    @ignore_docking_warnings
    def interface_property_packages(self):
        """Add junctions in between streams which have incompatible
        property packages."""
        path = self._path
        Stream = bst.Stream
        Interface = (bst.Junction, bst.Mixer, bst.MixTank)
        isa = isinstance
        new_path = []
        for obj in path:
            new_path.append(obj)
            outs = obj.outs
            for s in outs:
                source = s._source
                sink = s._sink
                if not sink or isa(sink, Interface): continue
                if sink.chemicals is not source.chemicals:
                    chemicals = s.chemicals
                    source_index = source._outs.index(s)
                    sink_index = sink._ins.index(s)
                    if chemicals is sink.chemicals:
                        s_sink = s
                        s_source = Stream(thermo=source.thermo)
                        s_source.copy_like(s)
                    else:
                        s_sink = Stream(thermo=sink.thermo)
                        s_sink.copy_like(s)
                        if chemicals is source.chemicals:
                            s_source = s
                        else:
                            s_source = Stream(thermo=source.thermo)
                            s_source.copy_like(s)
                    junction = bst.Junction(upstream=s_source, downstream=s_sink)
                    new_path.append(junction)
                    source._outs[source_index] = s_source
                    sink._ins[sink_index] = s_sink
        for obj in path:
            if isa(obj, System): obj.interface_property_packages()
        self._path = tuple(new_path)
        self._save_configuration()

    def _delete_path_cache(self):
        for i in ('_subsystems', '_units', '_unit_path', '_cost_units',
                  '_streams', '_feeds', '_products'):
            if hasattr(self, i): delattr(self, i)
        self._path_cache.clear()
        self._temporary_connections_log.clear()
        self._prioritized_units.clear()

    def copy(self, ID=None):
        """Copy system."""
        new = System(ID)
        new.copy_like(self)
        return new

    def copy_like(self, other: System):
        """Copy path, facilities and recycle from other system."""
        self._path = other._path
        self._facilities = other._facilities
        self._recycle = other._recycle
        self._connections = other._connections

    def set_tolerance(self, mol: Optional[float]=None, rmol: Optional[float]=None,
                      T: Optional[float]=None, rT: Optional[float]=None, 
                      subsystems: bool=False, maxiter: Optional[int]=None, 
                      subfactor: Optional[float]=None, method: Optional[str]=None,
                      algorithm: Optional[str]=None,
                      maxtime:Optional[float]=None):
        """
        Set the convergence tolerance and convergence method of the system.

        Parameters
        ----------
        mol :
            Molar tolerance [kmol/hr].
        rmol :
            Relative molar tolerance.
        T :
            Temperature tolerance [K].
        rT :
            Relative temperature tolerance.
        subsystems :
            Whether to set tolerance and method of subsystems as well.
        maxiter :
            Maximum number if iterations.
        subfactor :
            Factor to rescale tolerance in subsystems.
        method :
            Numerical method.
        algorithm :
            Convergence algorithm.
        maxtime :
            Maximum convergence time [s].
        
        """
        if mol is not None: self.molar_tolerance = float(mol)
        if rmol is not None: self.relative_molar_tolerance = float(rmol)
        if T is not None: self.temperature_tolerance = float(T)
        if rT is not None: self.temperature_tolerance = float(rT)
        if maxiter is not None: self.maxiter = int(maxiter)
        if method is not None: self.method = method
        if algorithm is not None: self.algorithm = algorithm
        if maxtime is not None: self.maxtime = maxtime
        if subsystems:
            if subfactor:
                for i in self.subsystems: i.set_tolerance(*[(i * subfactor if i else i) for i in (mol, rmol, T, rT)],
                                                          subsystems, maxiter, subfactor, method, algorithm, maxtime)
            else:
                for i in self.subsystems: i.set_tolerance(mol, rmol, T, rT, subsystems, maxiter, subfactor, method, algorithm, maxtime)

    ins = MockSystem.ins
    outs = MockSystem.outs
    load_inlet_ports = MockSystem.load_inlet_ports
    load_outlet_ports = MockSystem.load_outlet_ports
    get_inlet = MockSystem.get_inlet
    get_outlet = MockSystem.get_outlet
    set_inlet = MockSystem.set_inlet
    set_outlet = MockSystem.set_outlet
    _load_flowsheet  = MockSystem._load_flowsheet

    def _load_stream_links(self):
        for u in self.units: u._load_stream_links()

    @property
    def TEA(self) -> TEA:
        """TEA object linked to the system."""
        return getattr(self, '_TEA', None)

    @property
    def LCA(self):
        """QSDsan.LCA object linked to the system."""
        return getattr(self, '_LCA', None)

    specification = Unit.specification
    specifications = Unit.specifications
    add_bounded_numerical_specification = Unit.add_bounded_numerical_specification
    def add_specification(self, 
            specification: Optional[Callable]=None, 
            args: Optional[tuple]=(),
            simulate: Optional[bool]=None,
        ):
        """
        Add a specification.

        Parameters
        ----------
        specification : 
            Function runned for mass and energy balance.
        args : 
            Arguments to pass to the specification function.
        simulate :
            Whether to simulate after specification.

        Examples
        --------
        :doc:`../tutorial/Process_specifications`

        See Also
        --------
        add_bounded_numerical_specification
        specifications

        Notes
        -----
        This method also works as a decorator.

        """
        if not specification: return lambda specification: self.add_specification(specification, args, simulate)
        if not callable(specification): raise ValueError('specification must be callable')
        self._specifications.append(SystemSpecification(specification, args))
        if simulate is not None: self.simulate_after_specifications = simulate
        return specification

    def _extend_recycles(self, recycles):
        isa = isinstance
        recycle = self._recycle
        if recycle:
            if isa(recycle, Stream):
                recycles.append(recycle)
            elif isa(recycle, abc.Iterable):
                recycles.extend(recycle)
            else:
                raise_recycle_type_error(recycle)
        for i in self._path:
            if isa(i, System): i._extend_recycles(recycles)

    def get_all_recycles(self):
        recycles = []
        self._extend_recycles(recycles)
        return recycles

    def _extend_flattend_path_and_recycles(self, path, recycles, stacklevel):
        isa = isinstance
        recycle = self._recycle
        stacklevel += 1
        if recycle:
            if isa(recycle, Stream):
                recycles.append(recycle)
            elif isa(recycle, abc.Iterable):
                recycles.extend(recycle)
            else:
                raise_recycle_type_error(recycle)
        for i in self._path:
            if isa(i, System):
                if i.facilities:
                    warning = RuntimeWarning('subsystem with facilities could not be flattened')
                    warn(warning, stacklevel=stacklevel)
                    path.append(i)
                elif i.specifications:
                    warning = RuntimeWarning('subsystem with specification could not be flattened')
                    warn(warning, stacklevel=stacklevel)
                    path.append(i)
                else:
                    i._extend_flattend_path_and_recycles(path, recycles, stacklevel)
            else:
                path.append(i)

    def prioritize_unit(self, unit: Unit):
        """
        Prioritize unit operation to run first within it's recycle system,
        if there is one.

        Parameters
        ----------
        unit : 
            Unit operation to prioritize.

        Raises
        ------
        ValueError
            When unit is not in the system.
        RuntimeError
            When prioritization algorithm fails. This should never happen.

        Examples
        --------
        Create a simple recycle loop and prioritize a different unit operation:

        >>> from biosteam import main_flowsheet as f, Stream, settings, Mixer, Splitter
        >>> f.set_flowsheet('simple_recycle_loop')
        >>> settings.set_thermo(['Water'], cache=True)
        >>> feedstock = Stream('feedstock', Water=1000)
        >>> water = Stream('water', Water=10)
        >>> recycle = Stream('recycle')
        >>> product = Stream('product')
        >>> M1 = Mixer('M1', [feedstock, water, recycle])
        >>> S1 = Splitter('S1', M1-0, [product, recycle], split=0.5)
        >>> recycle_loop_sys = f.create_system('recycle_loop_sys')
        >>> recycle_loop_sys.print()
        System('recycle_loop_sys',
            [M1,
             S1],
            recycle=S1-1)
        >>> recycle_loop_sys.prioritize_unit(S1)
        >>> recycle_loop_sys.print()
        System('recycle_loop_sys',
            [S1,
             M1],
            recycle=S1-1)

        """
        isa = isinstance
        if unit not in self.unit_path:
            raise ValueError(f'unit {repr(unit)} not in system')
        path = self._path
        for index, other in enumerate(path):
            if unit is other:
                if (self._recycle or self.N_runs): 
                    self._path = path[index:] + path[:index]
                    self.method = 'fixed-point'
                del self._unit_path
                return
            elif isa(other, System) and unit in other.unit_path:
                other.prioritize_unit(unit)
                del self._unit_path
                return
        raise RuntimeError('problem in system algorithm')

    def find_system(self, unit: Unit):
        """
        Return system containing given unit within it's path.
        """
        isa = isinstance
        for i in self.path:
            if isa(i, System):
                if unit in i.units: return i.find_system(unit)
            elif i is unit:
                return self
        raise ValueError(f"unit {repr(unit)} not within system {repr(self)}")

    def path_section(self, starts, ends, inclusive=False):
        starts = tuple(starts)
        ends = tuple(ends)
        key = (starts, ends)
        if key in self._path_cache:
            path, end = self._path_cache[key]
        else:
            relevant_units = set(starts)
            for start in starts: relevant_units.update(start.get_downstream_units())
            unit_path = self.unit_path
            start_index = min([unit_path.index(start) for start in starts])
            end_index = max([unit_path.index(end) for end in ends])
            start = unit_path[start_index]
            end = unit_path[end_index]
            path = self.path_segment(start, end, False, relevant_units)
            self._path_cache[key] = (path, end)
        if inclusive: path = [*path, end]
        return path

    def path_segment(self, start, end, inclusive=False, 
                     relevant_units=None, critical_units=None):
        key = (start, end)
        if key in self._path_cache:
            segment = list(self._path_cache[key])
        else:
            if relevant_units is None: 
                relevant_units = start.get_downstream_units()
                relevant_units.add(start)
            if critical_units is None:
                critical_units = set(start.path_until(end))
            isa = isinstance
            if end not in relevant_units: return []
            path = self.path
            segment = []
            if start is None: # Need to make sure critical units are runned first in recycles
                for i, obj in enumerate(path):
                    if obj is end or isa(obj, System) and end in obj.units:
                        leftover_path = path[i + 1:]
                        break
                else:
                    raise ValueError(f"end unit {repr(end)} not in system")
                for obj in leftover_path:
                    if obj in critical_units:
                        critical_units.discard(obj)
                    elif isa(obj, System) and critical_units.intersection(obj.unit_set):
                        critical_units.difference_update(obj.unit_set)
                    else:
                        continue
                    segment.append(obj)
            else:
                for i, obj in enumerate(path):
                    if isa(obj, System):
                        if start in obj.unit_set:
                            if end in obj.unit_set:
                                return obj.path_segment(start, end, False, relevant_units, critical_units)
                            else:
                                path = path[i:]
                                break
                    elif obj is start:
                        if self.recycle:
                            path = path[i:] + path[:i] # recycle loop should start here
                        else:
                            path = path[i:] # start is appended in the next loop
                        break
                else:
                    raise ValueError(f"start unit {repr(start)} not in system")
            for i in path:
                if isa(i, System):
                    if end in i.unit_set:
                        segment.extend(i.path_segment(None, end, False, relevant_units, critical_units))
                        break
                elif i is end:
                    break  
                if isa(i, Unit) and i in relevant_units:
                    critical_units.discard(i)
                    segment.append(i)
                elif isa(i, System) and relevant_units.intersection(i.unit_set):
                    critical_units.difference_update(i.unit_set)
                    segment.append(i)
            else:
                raise ValueError(f"end unit {repr(end)} not in system")
            self._path_cache[key] = tuple(segment)
        if inclusive: segment = [*segment, end]
        return segment

    def split(self, 
              stream: Stream,
              ID_upstream: Optional[str]=None,
              ID_downstream: Optional[str]=None):
        """
        Split system in two; upstream and downstream.

        Parameters
        ----------
        stream : 
            Stream where unit group will be split.
        ID_upstream : 
            ID of upstream system.
        ID_downstream : 
            ID of downstream system.

        Examples
        --------
        >>> from biorefineries import cellulosic
        >>> from biosteam import default
        >>> cs = cellulosic.Biorefinery() # Create corn stover biorefinery
        >>> upstream_sys, downstream_sys = cs.cornstover_sys.split(cs.M201-0)
        >>> upstream_group = upstream_sys.to_unit_group()
        >>> upstream_group.show()
        UnitGroup: Unnamed
         units: U101, H2SO4_storage, T201, M201
        >>> downstream_group = downstream_sys.to_unit_group()
        >>> for i in upstream_group: assert i not in downstream_group.units
        >>> assert set(upstream_group.units + downstream_group.units) == set(cs.cornstover_sys.units)
        >>> default() # Reset to biosteam defaults

        """
        if self._recycle: raise RuntimeError('cannot split system with recycle')
        path = self._path
        streams = self.streams
        surface_units = {i for i in path if isinstance(i, Unit)}
        if stream.source in surface_units:
            index = path.index(stream.source) + 1
        elif stream.sink in surface_units:
            index = path.index(stream.sink)
        elif stream not in streams:
            raise ValueError('stream not in system')
        else:
            raise ValueError('stream cannot reside within a subsystem')
        return (System(ID_upstream, path[:index], None),
                System(ID_downstream, path[index:], None, self._facilities))

    def flatten(self):
        """Flatten system by removing subsystems."""
        recycles = []
        path = []
        self._extend_flattend_path_and_recycles(path, recycles, stacklevel=2)
        N = len(path)
        for i in range(N * N):
            stop = True
            for n, i in enumerate(path[:-1]):
                j = path[n+1]
                switched = [s.sink is i for s in j.outs if s.sink]
                if switched and all(switched):
                    path.remove(i)
                    path.insert(n+1, i)
                    stop = False
                    break
            if stop: break
        self._path = tuple(path)
        self._recycle = find_recycles(path)

    def to_unit_group(self, name: Optional[str]=None):
        """Return a UnitGroup object of all units within the system."""
        return bst.UnitGroup(name, self.units)

    def _set_path(self, path):
        #: tuple[Unit, function and/or System] A path that is run element
        #: by element until the recycle converges.
        self._path = path = tuple(path)

    def _set_facilities(self, facilities):
        #: tuple[Unit, function, and/or System] Offsite facilities that are simulated only after completing the path simulation.
        self._facilities = facilities = tuple(facilities)
        if self._shared_facilities: facilities = sum([i.facilities for i in self.subsystems], facilities)
        isa = isinstance
        units = self.cost_units
        for i in facilities:
            if isa(i, Facility):
                i._other_units = other_units = units.copy()
                other_units.discard(i)
        units = set(self.units)
        for i in facilities: units.discard(i)
        recycles = [
            outlet 
            for facility in facilities
            for outlet in facility.outs
            if outlet.sink in units
        ]
        if recycles:
            self._integrated_facilities = True
            system_recycles = self.recycle
            if system_recycles:
                if isinstance(system_recycles, abc.Iterable):
                    recycles.extend(i)
                else:
                    recycles.append(i)
            self.recycle = recycles
        else:
            self._integrated_facilities = False

    # Forward pipping
    __sub__ = Unit.__sub__
    __rsub__ = Unit.__rsub__

    # Backwards pipping
    __pow__ = __sub__
    __rpow__ = __rsub__

    def _get_source_stage(self, stream, stages_set):
        parent_source = stream.source
        if parent_source is None: return None
        ID = stream.ID
        while parent_source not in stages_set:
            IDs = [i.ID for i in parent_source.outs]
            index = IDs.index(ID)
            auxlet = parent_source.auxouts[index]
            parent_source = auxlet.source
            if parent_source is None: break
        return parent_source
    
    def _get_sink_stage(self, stream, stages_set):
        parent_sink = stream.sink
        if parent_sink in stages_set:
            return parent_sink
        else:
            index = parent_sink.ins.index(stream)
            auxlet = parent_sink.auxins[index]
            return auxlet.sink

    @property
    def subsystems(self) -> list[System]:
        """All subsystems in the system."""
        try:
            return self._subsystems
        except:
            self._subsystems = [i for i in self._path if isinstance(i, System)]
            return self._subsystems

    @property
    def stages(self):
        try:
            return self._stages
        except:
            self._stages = stages = []
            for i in self.units:
                stages.extend(getattr(i, 'stages', [i]))
            return stages

    @property
    def aggregated_stages(self):
        try:
            return self._aggregated_stages
        except:
            self._aggregated_stages = aggregated_stages = []
            for i in self.units:
                aggregated_stages.extend(getattr(i, 'aggregated_stages', [i]))
            return aggregated_stages

    @property
    def units(self) -> list[Unit]:
        """All unit operations as ordered in the path without repetitions."""
        try:
            return self._units
        except:
            self._units = units = []
            self._unit_set = past_units = set()
            isa = isinstance
            for i in self._path + self._facilities:
                if isa(i, Unit):
                    if i in past_units: continue
                    units.append(i)
                    past_units.add(i)
                elif isa(i, System):
                    sys_units = i.units
                    units.extend([i for i in sys_units if i not in past_units])
                    past_units.update(sys_units)
            return units

    @property
    def unit_set(self) -> set[Unit]:
        """Set of all unit operations."""
        try:
            return self._unit_set
        except:
            units = []
            isa = isinstance
            for i in self._path + self._facilities:
                if isa(i, Unit):
                    units.append(i)
                elif isa(i, System):
                    units.extend(i.units)
            self._unit_set = unit_set = set(units)
            return unit_set

    @property
    def unit_path(self) -> list[Unit]:
        """Unit operations as ordered in the path (some units may be repeated)."""
        try:
            return self._unit_path
        except:
            self._unit_path = units = []
            isa = isinstance
            for i in self._path + self._facilities:
                if isa(i, Unit):
                    units.append(i)
                elif isa(i, System):
                    units.extend(i.unit_path)
            return units

    @property
    def cost_units(self) -> set[Unit]:
        """All unit operations with costs."""
        try:
            return self._cost_units
        except:
            self._cost_units = units = set()
            isa = isinstance
            for i in self._path + self._facilities:
                if isa(i, Unit) and (i._design or i._cost):
                    units.add(i)
                elif isa(i, System):
                    units.update(i.cost_units)
            return units
      
    @property
    def streams(self) -> list[Stream]:
        """All streams within the system."""
        try:
            return self._streams
        except:
            self._streams = streams = []
            stream_set = set()
            for u in self.units:
                for s in u._ins + u._outs:
                    if not s: s = s.materialize_connection()
                    elif s in stream_set: continue
                    elif s.__class__ is AbstractStream: continue
                    streams.append(s)
                    stream_set.add(s)
            return streams 
    
    @property
    def feeds(self) -> list[Stream]:
        """All feeds to the system."""
        try:
            return self._feeds
        except:
            self._feeds = feeds = utils.feeds(self.streams)
            return feeds
    @property
    def products(self) -> list[Stream]:
        """All products of the system."""
        try:
            return self._products
        except:
            self._products = products = utils.products(self.streams)
            return products

    @property
    def facilities(self)  -> tuple[Facility, ...]:
        """All system facilities."""
        return self._facilities

    @property
    def recycle(self) -> Stream|None:
        """
        A tear stream for the recycle loop.
        """
        return self._recycle
    @recycle.setter
    def recycle(self, recycle):
        isa = isinstance
        if not recycle:
            self._recycle = None
        elif isa(recycle, Stream):
            self._recycle = recycle
        elif isa(recycle, abc.Iterable):
            real_recycles = [i for i in recycle if isa(i, Stream)]
            if len(real_recycles) == 0: 
                self.method = 'fixed-point'
                permanent = self.unit_set
                unit_path = self.unit_path
                for unit in unit_path:
                    if len(unit.outs) != 1: continue
                    stream = unit.outs[0]
                    if stream.sink in permanent:
                        self._recycle = stream
                        return
                for unit in unit_path:
                    self._recycle = unit.outs[0]
                    return
                return
            recycle = sorted(set(real_recycles), key=lambda x: x._ID)
            for i in recycle:
                if not isa(i, Stream):
                    raise AttributeError("recycle streams must be Stream objects; "
                                        f"not {type(i).__name__}")
            self._recycle = recycle
        elif recycle.__class__ is AbstractStream:
            self.method = 'fixed-point'
            permanent = self.unit_set
            unit_path = self.unit_path
            for unit in unit_path:
                if len(unit.outs) != 1: continue
                stream = unit.outs[0]
                if stream.sink in permanent:
                    self._recycle = stream
                    return
            for unit in unit_path:
                self._recycle = unit.outs[0]
                return
        else:
            raise_recycle_type_error(recycle)
    
    @property
    def depth(self):
        if self.recycle: 
            depth = [1]
        else:
            depth = [0]
        for i in self.subsystems:
            depth.append(i.depth + 1)
        return max(depth)

    @property
    def N_runs(self) -> int|None:
        """Number of times to converge the path."""
        return self._N_runs
    @N_runs.setter
    def N_runs(self, N_runs):
        self._N_runs = N_runs

    @property
    def path(self) -> list[Unit|System]:
        """A path that is run element by element until the recycle(s) converges (if any)."""
        return self._path

    @property
    def method(self) -> str:
        """Iterative convergence accelerator method ('wegstein', 'aitken', or 'fixedpoint')."""
        return self._method
    @method.setter
    def method(self, method):
        method = method.lower().replace('-', '').replace('_', '').replace(' ', '')
        if method in self.available_methods:
            self._method = method
        else:
            raise AttributeError(
                f"method '{method}' not available; only "
                f"{list_available_names(self.available_methods)} are available"
            )

    converge_method = method # For backwards compatibility

    @property
    def algorithm(self) -> str:
        """Iterative convergence algorithm ('Sequatial modular', or 'Phenomena based').
        
        Notes
        -----
        Timulation algorithms are available:
        
        Sequential modular - squentially runs each unit and subsystem until 
        the recycle (i.e. tear) stream convergences.
        
        Phenomena oriented - partitions and linearizes the
        underlying phenomenological equations for rapid and robust convergence.
        
        """
        return self._algorithm
    @algorithm.setter
    def algorithm(self, algorithm):
        algorithm = algorithm.capitalize().replace('-', ' ').replace('_', '')
        if algorithm in self.default_methods:
            self._algorithm = algorithm
        else:
            raise AttributeError(
                f"method '{algorithm}' not available; only "
                f"{list_available_names(self.default_methods)} are available"
            )

    @property
    def isdynamic(self) -> bool:
        """Whether the system contains any dynamic Unit."""
        if hasattr(self, '_isdynamic'):
            if self._isdynamic is not None:
                return self._isdynamic
        else:
            isdynamic = False
            for i in self.units:
                if i._isdynamic: 
                    isdynamic = True
                    break
            self._isdynamic = isdynamic
        return isdynamic
    @isdynamic.setter
    def isdynamic(self, isdynamic):
        self._isdynamic = isdynamic
    
    @property
    def _n_rotate(self):
        nr = 0
        for u in self.units:
            if not u._isdynamic: nr += 1
            else: break
        return nr

    def _downstream_path(self, unit):
        """Return a list composed of the `unit` and everything downstream."""
        if unit not in self.unit_path: return []
        elif self._recycle: return self._path
        isa = isinstance
        for index, i in enumerate(self._path):
            if unit is i:
                return self._path[index:]
            elif (isa(i, System) and unit in i.unit_path):
                return i._downstream_path(unit) + self._path[index+1:]
        return []

    def _downstream_facilities(self, unit):
        """Return a list of facilities composed of the `unit` and
        everything downstream."""
        isa = isinstance
        for index, i in enumerate(self._facilities):
            if unit is i or (isa(i, System) and unit in i.unit_path):
                return self._facilities[index:]
        return []

    def _downstream_system(self, unit):
        """Return a system with a path composed of the `unit` and
        everything downstream (facilities included)."""
        if self._recycle or unit is self._path[0]: return self
        path = self._downstream_path(unit)
        if path:
            facilities = self._facilities
        else:
            facilities = self._downstream_facilities(unit)
            if not facilities:
                raise RuntimeError(f'{unit} not found in system')
        system = System(None, path,
                        facilities=facilities)
        system._ID = f'{type(unit).__name__}-{unit} and downstream'
        return system

    def _minimal_digraph(self, facilities, graph_attrs):
        """Return digraph of the path as a box."""
        units = self.units
        if not facilities: units = [i for i in units if not isinstance(i, bst.Facility)]
        return minimal_digraph(self.ID, units, self.streams, **graph_attrs)

    def _surface_digraph(self, facilities, graph_attrs):
        units = self._path
        if not facilities: units = [i for i in units if not isinstance(i, bst.Facility)]
        return surface_digraph(units, **graph_attrs)

    def _thorough_digraph(self, facilities, graph_attrs):
        units = self.unit_path
        if not facilities: units = [i for i in units if not isinstance(i, bst.Facility)]
        return digraph_from_units(units, self.streams, **graph_attrs)

    def _cluster_digraph(self, facilities, graph_attrs):
        path = (*self.path, *self.facilities) if facilities else self.path
        return digraph_from_system(self, path, **graph_attrs)

    def _stage_digraph(self, graph_attrs):
        with self.stage_configuration(aggregated=False) as conf:
            return digraph_from_units(conf.stages, conf.streams, **graph_attrs)

    def diagram(self, kind: Optional[int|str]=None, file: Optional[str]=None, 
                format: Optional[str]=None, display: Optional[bool]=True,
                number: Optional[bool]=None, profile: Optional[bool]=None,
                label: Optional[bool]=None, title: Optional[str]=None,
                auxiliaries: Optional[bool]=None, facilities: Optional[bool]=True,
                **graph_attrs):
        """
        Display a `Graphviz <https://pypi.org/project/graphviz/>`__ diagram of
        the system.

        Parameters
        ----------
        kind :
            * 0 or 'cluster': Display all units clustered by system.
            * 1 or 'thorough': Display every unit within the system.
            * 2 or 'surface': Display only elements listed in the path.
            * 3 or 'minimal': Display a single box representing the system.
            * 4 or 'stage': Display every stage within the system.
        file : 
            File name to save diagram. 
        format:
            File format (e.g. "png", "svg"). Defaults to 'png'
        display : 
            Whether to display diagram in console or to return the graphviz
            object.
        number : 
            Whether to number unit operations according to their
            order in the system path.
        profile : 
            Whether to clock the simulation time of unit operations.
        label : 
            Whether to label the ID of streams with sources and sinks.
        auxiliaries:
            Whether to include auxiliary units.
        
        """
        self._load_configuration()
        if kind is None: kind = 1
        if title is None: title = ''
        graph_attrs['label'] = title
        if auxiliaries is not None: graph_attrs['auxiliaries'] = auxiliaries 
        preferences = bst.preferences
        with preferences.temporary():
            if number is not None: preferences.number_path = number
            if label is not None: preferences.label_streams = label
            if profile is not None: preferences.profile = profile
            if format is not None: preferences.graphviz_format = format
            if kind == 0 or kind == 'cluster':
                f = self._cluster_digraph(facilities, graph_attrs)
            elif kind == 1 or kind == 'thorough':
                f = self._thorough_digraph(facilities, graph_attrs)
            elif kind == 2 or kind == 'surface':
                f = self._surface_digraph(facilities, graph_attrs)
            elif kind == 3 or kind == 'minimal':
                f = self._minimal_digraph(facilities, graph_attrs)
            elif kind == 4 or kind == 'stage':
                f = self._stage_digraph(graph_attrs)
            else:
                raise ValueError("kind must be one of the following: "
                                 "0 or 'cluster', 1 or 'thorough', 2 or 'surface', "
                                 "3 or 'minimal'")
            if display or file:
                def size_key(units):
                    N = len(units)
                    for u in units: N += len(u.auxiliary_units)
                    if N < 2:
                        return 'unit'
                    elif N < 6:
                        return 'system'
                    elif N < 9:
                        return 'big-system'
                    else:
                        return 'huge-system'
                height = (
                    preferences.graphviz_html_height
                    [size_key(self.units)]
                    [preferences.tooltips_full_results]
                )
                finalize_digraph(f, file, format, height)
            else:
                return f

    # Methods for running one iteration of a loop
    def _iter_run_conditional(self, data):
        """
        Run the system at specified recycle molar flow rate.

        Parameters
        ----------
        data : numpy.ndarray
              Recycle molar flow rates, temperature, and pressure.

        Returns
        -------
        mol_new : numpy.ndarray
            New recycle molar flow rates.
        not_converged : bool
            True if recycle has not converged.

        """
        data[data < 0.] = 0.
        self._set_recycle_data(data)
        T = self._get_recycle_temperatures()
        mol = self._get_recycle_mol()
        self.run()
        recycle = self._recycle
        for i, j in self.tracked_recycles.items():
            if i is recycle: j.append(i.copy(None))
        mol_new = self._get_recycle_mol()
        T_new = self._get_recycle_temperatures()
        mol_errors = abs(mol - mol_new)
        if mol_errors.any():
            self._mol_error = mol_error = mol_errors.max()
            if mol_error > 1e-12:
                nonzero_index = mol_errors.nonzero_index()
                mol_errors = mol_errors[nonzero_index]
                max_errors = np.maximum.reduce([abs(mol[nonzero_index]), abs(mol_new[nonzero_index])])
                self._rmol_error = rmol_error = (mol_errors / max_errors).max()
            else:
                self._rmol_error = rmol_error = mol_error
        else:
            self._mol_error = mol_error = 0.
            self._rmol_error = rmol_error = 0.
        T_errors = abs(T - T_new)
        self._T_error = T_error = T_errors.max()
        self._rT_error = rT_error = (T_errors / T).max()
        self._iter += 1
        not_converged = not (
            (mol_error < self.molar_tolerance
             or rmol_error < self.relative_molar_tolerance)
            and
            (T_error < self.temperature_tolerance
             or rT_error < self.relative_temperature_tolerance)
        )
        if not_converged:
            if (self._iter >= self.maxiter
                or (self.maxtime and self.timer.elapsed_time > self.maxtime)): 
                if self.strict_convergence: 
                    raise RuntimeError(f'{repr(self)} could not converge' + self._error_info())
                else: 
                    not_converged = False
        return self._get_recycle_data(), not_converged
        
    def _iter_run(self, data):
        """
        Run the system at specified recycle molar flow rate.

        Parameters
        ----------
        data : numpy.ndarray
              Recycle temperature, pressure, and molar flow rates.

        Returns
        -------
        mol_new : numpy.ndarray
            New recycle molar flow rates.

        """
        data_new, not_converged = self._iter_run_conditional(data)
        if not_converged: return data_new - data
        else: raise Converged

    def _get_recycle_mol(self):
        recycles = self.get_all_recycles()
        N = len(recycles)
        if N == 1:
            return recycles[0].mol.copy()
        elif N > 1: 
            return SparseArray([i.mol.copy() for i in recycles])
        else:
            raise RuntimeError('no recycle available')

    def _get_recycle_data(self):
        recycles = self.get_all_recycles()
        N = len(recycles)
        if N == 1:
            return get_recycle_data(recycles[0])
        elif N > 1: 
            return np.hstack([get_recycle_data(i) for i in recycles])
        else:
            raise RuntimeError('no recycle available')

    def _set_recycle_data(self, data):
        recycles = self.get_all_recycles()
        N = len(recycles)
        if N == 1:
            set_recycle_data(recycles[0], data)
        elif N > 1:
            N = data.size
            M = sum([i.imol.data.size + 2 for i in recycles])
            if M != N: 
                raise IndexError(f'expected {N} elements; got {M} instead')
            index = 0
            for i in recycles:
                end = index + i.imol.data.size + 2
                set_recycle_data(i, data[index:end])
                index = end
        else:
            raise RuntimeError('no recycle available')

    def _get_recycle_temperatures(self):
        recycle = self._recycle
        if isinstance(recycle, Stream):
            T = self._recycle.T
        elif isinstance(recycle, abc.Iterable):
            T = [i.T for i in recycle]
        else:
            raise RuntimeError('no recycle available')
        return np.array(T, float)

    def _create_temporary_connections(self):
        """Create temporary connections based on process specifications."""
        temporary_connections_log = self._temporary_connections_log
        for u in self.units:
            for ps in u._specifications: 
                if ps in temporary_connections_log: continue
                ps.create_temporary_connections(u)
                temporary_connections_log.add(ps)

    def _setup_units(self):
        """Setup all unit operations."""
        prioritized_units = self._prioritized_units
        for u in self.units: 
            u._system = self
            u._setup()
            u._check_setup()
            if u not in prioritized_units:
                for ps in u._specifications: ps.compile_path(u)
                if u.prioritize: self.prioritize_unit(u)
                prioritized_units.add(u)
                
    def _setup(self, update_configuration=False, units=None, load_configuration=True):
        """Setup each element of the system."""
        if units is None: units = self.units
        if update_configuration:
            self._temporary_connections_log.clear()
            self._create_temporary_connections()
            self._update_configuration(units=[*units, *temporary_units_dump])
            temporary_units_dump.clear()
            self._setup_units()
            self._remove_temporary_units()
            self._save_configuration()
            self._load_stream_links()
        else:
            if load_configuration: self._load_configuration()
            self._create_temporary_connections()
            if temporary_units_dump:
                self._update_configuration(units=[*units, *temporary_units_dump])
                temporary_units_dump.clear() 
                self._setup_units()
                self._remove_temporary_units()
                self._save_configuration()
            else:
                self._setup_units()

    @piping.ignore_docking_warnings
    def _remove_temporary_units(self):
        isa = isinstance
        new_path = []
        for i in self.path:
            if isa(i, System): 
                i._remove_temporary_units()
            elif isa(i, TemporaryUnit):
                i.old_connection.reconnect()
                continue
            new_path.append(i)
        self._path = tuple(new_path)

    def run(self):
        algorithm = self.algorithm 
        if algorithm == 'Sequential modular':
            self.run_sequential_modular()
        else:
            if self._iter == 0:
                self.run_sequential_modular()
                with self.stage_configuration(aggregated=False) as conf:
                    conf.solve_material_flows()
                    if self.tracking_flag: self._collect_variables('material')
            elif algorithm == 'Phenomena based':
                self.run_phenomena()
            elif algorithm == 'Phenomena modular':
                self.run_phenomena_modular()
            else:
                raise RuntimeError(f'unknown algorithm {algorithm!r}')

    def run_sequential_modular(self):
        """Run mass and energy balances for each element in the path
        without costing unit operations."""
        isa = isinstance
        f = try_method_with_object_stamp
        flag = self.tracking_flag
        integrated = self._integrated_facilities
        for i in self._path:
            if isa(i, Unit):
                if integrated: f(i, i.simulate)
                else: f(i, i.run)
                if flag: self._collect_variables(i.ID)
            elif isa(i, System): 
                if integrated: f(i, i.simulate) 
                else: f(i, i.converge)
            else: raise RuntimeError('path elements must be either a unit or a system')
        if flag == TRACK_CONDITION_NUMBER:
            self._collect_recycle_condition_number()

    def run_phenomena(self):
        """Decouple and linearize material, equilibrium, summation, enthalpy,
        and reaction phenomena and iteratively solve them."""
        path = self.unit_path
        flag = self.tracking_flag
        with self.stage_configuration(aggregated=False) as conf:
            try:
                conf.solve_nonlinearities()
                if flag: self._collect_variables('phenomenode')
                conf.solve_energy_flows()
                if flag: self._collect_variables('energy')
                conf.solve_material_flows()
                if flag: self._collect_variables('material')
            except (NotImplementedError, UnboundLocalError, KeyError) as error:
                raise error
            except Exception as e:
                print('FAILED!')
                print(e)
                print('-------')
                for i in path: 
                    i.run()
                    if flag: self._collect_variables(i.ID)
                conf.solve_material_flows()
                if flag: self._collect_variables('material')

    def run_phenomena_modular(self):
        path = self.unit_path
        flag = self.tracking_flag
        with self.stage_configuration(aggregated=False) as conf:
            try:
                for n, i in enumerate(path): 
                    if isinstance(i, (bst.Mixer, bst.Splitter)) and not i.specifications: continue
                    i.run()
                    if flag: self._collect_variables((i.ID, 'phenomenode'))
                    conf.solve_material_flows()
                    if flag: self._collect_variables('material')
                    conf.solve_energy_flows()
                    if flag: self._collect_variables('energy')
            except (NotImplementedError, UnboundLocalError, KeyError) as error:
                raise error
            except Exception as e:
                print('-------')
                print('FAILED!')
                print(e)
                for i in path[n+1:]: 
                    i.run()
                    if flag: self._collect_variables(i.ID)

    def stage_configuration(self, aggregated=False):
        try:
            if aggregated:
                return self._aggregated_stage_configuration
            else:
                return self._stage_configuration
        except:
            stages = self.stages
            streams = list(set([i for s in stages for i in (*s.flat_ins, *s.flat_outs)]))
            connections = [(i, i.source, i.sink) for i in streams]
            if aggregated:
                stages = self.aggregated_stages
                streams = [i for u in stages for i in u.flat_ins + u.flat_outs]
                stream_ref = {i.material_reference: i for u in stages for i in (*u.flat_ins, *u.flat_outs)}
                nodes = [i for i in stages if not getattr(i, 'decoupled', False)]
                self._aggregated_stage_configuration = conf = Configuration(
                    self.path, stages, nodes, streams, stream_ref, connections, aggregated
                )
            else:
                stream_ref = {i.material_reference: i for i in streams}
                nodes = [i for i in stages if not getattr(i, 'decoupled', False)]
                self._stage_configuration = conf = Configuration(
                    self.path, stages, nodes, streams, stream_ref, connections, aggregated
                )
            return conf
        
    def track_convergence(
            self, 
            track=True, 
            condition_number=False,
        ):
        if track:
            self.tracking_flag = track + condition_number
            with self.stage_configuration(aggregated=False) as conf:
                for i in self.stages: 
                    i.create_equation_nodes(conf.stream_ref)
                    i.initialize_equation_nodes()
                self.equation_nodes = equations = tuple(set(sum([i.equation_nodes for i in self.stages], ())))
                self.variable_nodes = variables = tuple(set(sum([i.variables for i in equations], ())))
                self.variable_profiles = variable_profiles = {i: [] for i in variables}
                for i in conf.nodes:
                    specification_variables = [j for i in i.material_balance_specifications_nodes for j in i.variables]
                    for variable in specification_variables:
                        if variable not in variable_profiles: variable_profiles[variable] = []
            for i in self.units:
                i._system_collect_variables = self._collect_variables
                i.tracking = True
            self.equation_profiles = {}
            self.edge_profiles = {}
            self.grouped_variables = {}
            if condition_number:
                if self.algorithm == 'Sequential modular':
                    self.condition_number_profile = {
                        'recycle': [],
                    }
                else:
                    self.condition_number_profile = {
                        'material': [],
                        'energy': [],
                        'phenomenode': [],
                    }
            if not self.recycle:
                self.timer = Timer()
                self.timer.start()
                
        else:
            for i in self.variable_profiles: i.getter = None
            self.tracking_flag = False
            self.equation_nodes = None
            self.variable_nodes = None
            self.variable_profiles = None
            self.equation_profiles = None
            self.edge_profiles = None
            self.grouped_variables = None
            for i in self.units: i.tracking = False
    
    def get_profiles(self, equation_profiles, edge_profiles):
        time, variable_profiles = self.get_variable_profiles()
        if equation_profiles or edge_profiles:
            column_names = list(variable_profiles)
            namespace = self.flowsheet.to_dict()
            namespace['nan'] = np.nan
            chemicals = bst.settings.get_chemicals()
            directories = []
            stages = []
            for name in column_names:
                directory, chemical = name
                if directory[-1] == 'R':
                    directory = directory[:-1] + 'dmol'
                else:
                    if directory[-1] == 'F':
                        directory = directory[:-1] + 'mol'
                    if '.outs[' in directory:
                        dot_index = directory.rindex('.outs[')
                    elif '.ins[' in directory:
                        dot_index = directory.rindex('.ins[')
                    else:
                        dot_index = directory.rindex('.')
                stages.append(directory[:dot_index])
                if chemical != '-':
                    index = chemicals.index(chemical)
                    directory += f"[{index}]"
                directories.append(directory)
            stages = set(stages)
            if equation_profiles: equation_profiles = self.get_equation_profiles(variable_profiles, directories, stages, namespace) 
            if edge_profiles: edge_profiles = self.get_edge_profiles(variable_profiles, directories, stages, namespace) 
        return time, variable_profiles, equation_profiles, edge_profiles
    
    def get_phenomegraph(self, bipartite=False, **kwargs):
        if bipartite:
            return self.get_bipartite_phenomegraph(**kwargs)
        else:
            return self.get_monopartite_phenomegraph(**kwargs)
            
    def get_monopartite_phenomegraph(self, decomposition=None, dotfile=None):
        if decomposition is None: decomposition = self.algorithm.lower()
        subgraph_units = [i.ID for i in self.path]
        if 'phenomenode' in decomposition:
            if 'modular' in decomposition:
                criteria = [*all_subgraphs, *[(i, 'phenomenode') for i in subgraph_units]]
            else:
                criteria = all_subgraphs
        elif 'modular' in decomposition:
            criteria = list(product(subgraph_units, all_subgraphs))
        else:
            raise ValueError('unknown decompostion')
        stage_tags = [i.node_tag for i in self.stages]
        time, variable_profiles = self.get_variable_profiles()
        return PhenomeGraph(
            self.equation_nodes, criteria, stage_tags, variable_profiles,
            time,
        )
    
    def get_bipartite_phenomegraph(
            self,
            decomposition=None,
            get_equation_profiles=False,
            get_edge_profiles=False,
            dotfile=None, 
            **kwargs,
        ):
        if decomposition is None: decomposition = self.algorithm.lower()
        if not getattr(self, 'tracking_flag', False):
            with self.stage_configuration(aggregated=False) as conf:
                for i in self.stages: 
                    i.create_equation_nodes(conf.stream_ref)
                    i.initialize_equation_nodes()
                self.equation_nodes = equations = tuple(set(sum([i.equation_nodes for i in self.stages], ())))
                self.variable_nodes = variables = tuple(set(sum([i.variables for i in equations], ())))
                self.variable_profiles = variable_profiles = {i: [] for i in variables}
                for i in conf.nodes:
                    specification_variables = [j for i in i.material_balance_specifications_nodes for j in i.variables]
                    for variable in specification_variables:
                        if variable not in variable_profiles: variable_profiles[variable] = []
            self.equation_profiles = {}
            self.edge_profiles = {}
        if isinstance(dotfile, str): dotfile = DotFile(dotfile)
        variables = self.variable_nodes
        equations = self.equation_nodes
        edges = tuple(BipartitePhenomeGraph.get_node_edges(equations))
        name = self.ID + '_graph'
        subgraph_units = [i.ID for i in self.path]
        time, variable_profiles, equation_profiles, edge_profiles = self.get_profiles(
            get_equation_profiles, get_edge_profiles
        )
        graph = BipartitePhenomeGraph(
            name, equations, variables, edges, equation_profiles, 
            variable_profiles, edge_profiles, time,
            dotfile.file if dotfile else None, 
            subgraph_units=subgraph_units,
            **kwargs,
        )
        if 'phenomenode' in decomposition:
            if 'modular' in decomposition:
                names = [*all_subgraphs, *[(i, 'phenomenode') for i in subgraph_units]]
            else:
                names = all_subgraphs
        elif 'modular' in decomposition:
            names = subgraph_units
        else:
            raise ValueError('unknown decompostion')
        
        def equation_match(equation, subgraph):
            name = equation.name
            if isinstance(subgraph, tuple):
                return all([i in name for i in subgraph])
            else:
                return subgraph in name
            
        def variable_match(variable, subgraph):
            name = variable.name
            if subgraph in all_subgraphs:
                return True
            elif isinstance(subgraph, tuple):
                return subgraph[0] in name
            else:
                return subgraph in name
            
        for subgraph in names:
            subgraph_equations = [i for i in equations if equation_match(i, subgraph)]
            subgraph_variables = set(sum([i.outputs for i in subgraph_equations], []))
            subgraph_variables = set([i for i in subgraph_variables if variable_match(i, subgraph)])
            subgraph_edges = BipartitePhenomeGraph.get_node_edges(subgraph_equations)
            subgraph_edges = tuple([i for i in subgraph_edges if i.variable in subgraph_variables])
            subgraph_variables = tuple(subgraph_variables)
            subgraph_name = '.'.join([name, subgraph])
            subgraph_variable_profiles = variable_profiles[[i.name for i in subgraph_variables]]
            if get_edge_profiles:
                subgraph_edge_profiles = edge_profiles[[i.name for i in subgraph_edges]]
            else:
                subgraph_edge_profiles = None
            if get_equation_profiles:
                subgraph_equation_profiles = equation_profiles[[i.name for i in subgraph_equations]]
            else:
                subgraph_equation_profiles = None
            file = dotfile.subfile(subgraph) if dotfile else None
            graph.subgraph(
                subgraph_name, subgraph_equations, 
                subgraph_variables, subgraph_edges,
                subgraph_equation_profiles,
                subgraph_variable_profiles,
                subgraph_edge_profiles,
                file
            )
        return graph
    
    def get_variable_profiles(self):
        variable_profiles = self.variable_profiles
        variable_nodes = tuple(variable_profiles)
        names = [i.name for i in variable_nodes]
        chemicals = [i.ID for i in bst.settings.get_chemicals()]
        subnames = [('-' if np.ndim(i.last_value) == 0 else chemicals) for i in variable_nodes]
        columns = [(name, i) for name, subname in zip(names, subnames) for i in subname]
        M = len(next(iter(variable_profiles.values())))
        N = len(columns)
        values = np.zeros([M, N])
        for i in range(M):
            j = 0
            for lst in variable_profiles.values():
                values_i = lst[i]
                if np.ndim(values_i):
                    end = j + len(values_i)
                else:
                    end = j + 1
                values[i, j:end] = values_i
                j = end
        rows, = np.where((~np.isnan(values)).all(axis=1))
        df = pd.DataFrame(values[rows], columns=pd.MultiIndex.from_tuples(columns))
        record = np.array(self.timer.record)
        return record[rows], df
    
    def get_edge_profiles(self, variable_profiles_table, directories, stages, namespace):
        edge_profiles = self.edge_profiles
        M = variable_profiles_table.shape[0]
        steady_states = variable_profiles_table.iloc[-1]
        for i in range(M):
            row = variable_profiles_table.iloc[i]
            for directory, stage, value, steady_state in zip(directories, stages, row, steady_states):
                exec(f"{directory} = {value}", namespace)
                unit = eval(stage, namespace)
                for name, error in unit._collect_edge_errors():
                    if name in edge_profiles:
                        edge_profiles[name].append(error)
                    else:
                        edge_profiles[name] = [error]
                exec(f"{directory} = {steady_state}", namespace)
        columns = list(edge_profiles)
        N = len(columns)
        values = np.zeros([M, N])
        for i in range(M):
            for j, lst in enumerate(edge_profiles.values()):
                values_i = lst[i]
                values[i, j] = values_i
        return pd.DataFrame(values, columns=pd.MultiIndex.from_tuples(columns))
    
    def get_equation_profiles(self, variable_profiles_table, directories, stages, namespace):
        equation_profiles = self.equation_profiles
        M = variable_profiles_table.shape[0]
        for i in range(M):
            row = variable_profiles_table.iloc[i]
            for directory, value in zip(directories, row):
                exec(f"{directory} = {value}", namespace)
            for stage in stages:
                unit = eval(stage, namespace)
                for name, error in unit._collect_equation_errors():
                    if name in equation_profiles:
                        equation_profiles[name].append(error)
                    else:
                        equation_profiles[name] = [error]
                for name, error in unit._collect_specification_errors():
                    if name in equation_profiles:
                        equation_profiles[name].append(error)
                    else:
                        equation_profiles[name] = [error]
        columns = list(equation_profiles)
        N = len(columns)
        values = np.zeros([M, N])
        for j, lst in enumerate(equation_profiles.values()):
            for i in range(M):
                values_i = lst[i]
                values[i, j] = values_i
        return pd.DataFrame(values, columns=columns)
            
    def assert_tracking_ok(self):
        if not self.tracking_flag: return
        N, = set([len(i) for i in self.variable_profiles.values()])
        M = len(self.timer.record)
        assert M == N

    # TODO: Decide whether to remove condition number analysis.
    def _collect_recycle_condition_number(self):
        pass

    def _collect_material_condition_number(self):
        pass
    
    def _collect_energy_condition_number(self):
        pass
    
    def _collect_phenomena_condition_number(self):
        pass
            
    def _collect_variables(self, key):
        self.timer.measure()
        with self.timer.offset():
            if self.tracking_flag == TRACK_CONDITION_NUMBER:
                if key == 'material':
                    self._collect_material_condition_number()
                elif key == 'energy':
                    self._collect_energy_condition_number()
                elif key == 'phenomenode':
                    self._collect_phenomena_condition_number()
            grouped_variables = self.grouped_variables
            if isinstance(key, str): key = (key,)
            if key in grouped_variables:
                group = grouped_variables[key]
            else:
                equations = [
                    eq for eq in self.equation_nodes
                    if all([i in eq.name for i in key])
                ]
                group = [i for eq in equations for i in eq.tracked_outputs]
                grouped_variables[key] = group = set(group)
            variables = self.variable_profiles
            for variable, values in variables.items(): 
                if variable in group: 
                    value = variable.get_value()
                else:
                    value = variable.last_value
                values.append(value)
                
    def _mock_collect_variables(self, key, variable_errors, convergence_rate):
        grouped_variables = self.grouped_variables
        if isinstance(key, str): key = (key,)
        if key in grouped_variables:
            group = grouped_variables[key]
        else:
            equations = [
                eq for eq in self.equation_nodes
                if all([i in eq.name for i in key])
            ]
            group = [i for eq in equations for i in eq.tracked_outputs]
            grouped_variables[key] = group = set(group)
        variables = self.variable_profiles
        for variable, values in variables.items(): 
            if variable in group: 
                variable_errors[variable] = value = max(variable_errors[variable] - convergence_rate, 0)
            else:
                value = variable_errors[variable]
            values.append(value)
        self.timer.mock_measure()
        
    def _solve(self):
        """Solve the system recycle iteratively."""
        self._reset_iter()
        solver, conditional, kwargs = self.available_methods[self._method]
        data = self._get_recycle_data()
        f = self._iter_run_conditional if conditional else self._iter_run
        try: solver(f, data, **kwargs)
        except (IndexError, ValueError) as error:
            data = self._get_recycle_data()
            try: solver(f, data, **kwargs)
            except Converged: pass
            except: raise error
        except Converged: pass

    def get_recycle_data(self):
        """
        Return recycle data defining material flows and conditions of recycle streams.
        
        """
        return JointRecycleData(self.get_all_recycles(), self.responses)

    def converge(self, recycle_data: list[RecycleData]=None, update_recycle_data: bool=False):
        """
        Converge mass and energy balances. If material data was given, 
        return converged material flows at steady state. Shape will be M by N,
        where M is the number of recycle streams and N is the number of chemicals.
        
        Parameters
        ----------
        recycle_data : 
            Material data to set recycle streams.
            
        See Also
        --------
        :meth:`System~.get_recycle_data`
        
        Warning
        -------
        No design, costing, nor facility algorithms are run.
        To run full simulation algorithm, see :func:`~biosteam.System.simulate`.
        
        """
        if recycle_data is not None: recycle_data.reset()
        if self._recycle:
            method = self._solve
        else:
            method = self.run_sequential_modular
        if self._N_runs:
            for i in range(self._N_runs): method()
        else:
            method()
        if update_recycle_data:
            try: recycle_data.update()
            except AttributeError: raise ValueError('no recycle data to update')

    def _summary(self, force_path=False):
        simulated_units = set()
        isa = isinstance
        Unit = bst.Unit
        f = try_method_with_object_stamp
        if not self._integrated_facilities or force_path:
            for i in self._path:
                if isa(i, Unit):
                    if i in simulated_units: continue
                    simulated_units.add(i)
                f(i, i._summary)
        for i in self._facilities:
            if isa(i, Unit): f(i, i.simulate)
            elif isa(i, System):
                f(i, i.converge)
                i._summary()
            else: i() # Assume it is a function
        for i in self._facilities:
            if isa(i, bst.BoilerTurbogenerator): f(i, i.simulate)

    def _reset_iter(self):
        self._iter = 0
        if self.maxtime or self.tracking_flag: 
            self.timer = Timer()
            self.timer.start()
        for j in self.tracked_recycles.values(): j.clear()
        for system in self.subsystems: system._reset_iter()

    def _reset_errors(self):
        #: Molar flow rate error (kmol/hr)
        self._mol_error = 0

        #: Relative molar flow rate error
        self._rmol_error = 0

        #: Temperature error (K)
        self._T_error = 0

        #: Relative temperature error
        self._rT_error = 0

        #: Number of iterations
        self._iter = 0

    def empty_outlet_streams(self):
        """Reset all outlet streams to zero flow."""
        self._reset_errors()
        units = self.units
        streams = utils.streams_from_units(units)
        utils.filter_out_missing_streams(streams)
        streams_by_data = {}
        for i in streams:
            data = i.imol.data
            data_id = id(data)
            if data_id in streams_by_data:
                streams_by_data[data_id].append(i)
            else:
                streams_by_data[data_id] = [i]
        for streams in streams_by_data.values():
            if all([i.source in units for i in streams]):
                streams[0].empty()

    def empty_recycles(self):
        """Reset all recycle streams to zero flow."""
        self._reset_errors()
        recycle = self._recycle
        if recycle:
            if isinstance(recycle, Stream):
                recycle.empty()
            elif isinstance(recycle, abc.Iterable):
                for i in recycle: i.empty()
            else:
                raise_recycle_type_error(recycle)
        for system in self.subsystems:
            system.empty_recycles()

    def rescale(self, feedstock: Stream, ratio: float):
        """Rescale feedstock flow rate and update recycle stream flow rate guesses."""
        feedstock.rescale(ratio)
        for i in self.get_all_recycles(): i.rescale(ratio)

    def reset_cache(self):
        """Reset cache of all unit operations."""
        if self.isdynamic:
            self._DAE = None
            self._state = None
            for s in self.streams:
                s._state = None
                s._dstate = None
            self.dynsim_kwargs = {}
            self.scope.reset_cache()
        for unit in self.units:
            unit.reset_cache(self.isdynamic)
        for i in self.streams: i.reset_cache()
            
    def set_dynamic_tracker(self, *subjects, **kwargs):
        """
        Set up an :class:`SystemScope` object to track the dynamic data.

        Parameters
        ----------
        *subjects :
            Any subjects of the system to track, which must have an `.scope`
            attribute of type :class:`Scope`.
        """
        if self.isdynamic:
            self._scope = {'subjects':subjects, 'kwargs':kwargs}
        else:
            warn(f'{self.__repr__()} must have at least one dynamic unit to '
                 f'set up a dynamic tracker.')
        
    @property
    def scope(self) -> utils.SystemScope:
        """
        A tracker of dynamic data of the system, set up
        with :class:`System`.`set_dynamic_tracker()`
        """
        if not hasattr(self, '_scope'): self._scope = utils.SystemScope(self)
        elif isinstance(self._scope, dict): 
            # sys.converge() seems to break WasteStreamScope, so it's now 
            # set up to initiate the SystemScope object after converge() when 
            # the system is run the first time 
            subjects = self._scope['subjects']
            kwargs = self._scope['kwargs']
            self._scope = utils.SystemScope(self, *subjects, **kwargs)
        return self._scope
    
    # _hasode = lambda unit: hasattr(unit, '_compile_ODE')
    
    def _dstate_attr2arr(self, y):
        dy = y.copy()
        idx = self._state_idx
        for unit in self.units:
            if unit.hasode:
                start, stop = idx[unit._ID]
                dy[start: stop] = unit._dstate
        return dy

    def _update_state(self, arr):
        self._state[:] = arr
        for unit in self.units:
            if unit.hasode: unit._update_state()

    def _load_state(self):
        """Returns the initial state (a 1d-array) of the system for dynamic simulation."""
        nr = self._n_rotate
        units = self.units[nr:] + self.units[:nr]
        if self._state is None:
            for ws in self.feeds:
                if not ws.state.all(): ws._init_state()
            for inf in units[0].ins:
                if not inf.state.all(): inf._init_state()
            y = np.array([])
            idx = {}
            for unit in units: 
                unit._init_state()
                unit._update_state()
                unit._update_dstate()
                if unit.hasode:
                    start = len(y)
                    y = np.append(y, unit._state)
                    stop = len(y)
                    idx[unit._ID] = (start, stop)
            if len(y) == 0: y = np.array([0])
            self._state = y
            self._state_idx = idx
            for unit in units:
                if unit.hasode:
                    start, stop = idx[unit._ID]
                    unit._state = self._state[start: stop]
        else:
            y = self._state
            idx = self._state_idx
            for unit in units: unit._update_dstate()
        return y, idx, nr

    def _compile_DAE(self):
        nr = self._n_rotate
        units = self.units[nr:] + self.units[:nr]
        _update_state = self._update_state
        _dstate_attr2arr = self._dstate_attr2arr
        funcs = [u.ODE if u.hasode else u.AE for u in units]
        track = self.scope
        dk = self.dynsim_kwargs
        if dk.get('print_t'): # print integration time for debugging
            def dydt(t, y):
                _update_state(y)
                print(t)
                for unit, func in zip(units, funcs):
                    if unit.hasode:
                        QC_ins, QC, dQC_ins = unit._ins_QC, unit._state, unit._ins_dQC
                        func(t, QC_ins, QC, dQC_ins)   # updates dstate
                    else:
                        QC_ins, dQC_ins = unit._ins_QC, unit._ins_dQC
                        func(t, QC_ins, dQC_ins)   # updates both state and dstate
                track(t)
                return _dstate_attr2arr(y)
        else:
            def dydt(t, y):
                _update_state(y)
                for unit, func in zip(units, funcs):
                    if unit.hasode:
                        QC_ins, QC, dQC_ins = unit._ins_QC, unit._state, unit._ins_dQC
                        func(t, QC_ins, QC, dQC_ins)   # updates dstate
                    else:
                        QC_ins, dQC_ins = unit._ins_QC, unit._ins_dQC
                        func(t, QC_ins, dQC_ins)   # updates both state and dstate
                track(t)
                return _dstate_attr2arr(y)
        self._DAE = dydt            

    @property
    def DAE(self):
        """System-wide differential algebraic equations."""
        if self._DAE is None:
            try: self._compile_DAE()
            except AttributeError: return None
        return self._DAE

    def _write_state(self):
        for ws in [i for i in self.streams if i not in self.feeds]:
            ws._state2flows()

    def clear_state(self):
        """Clear all states and dstates (system, units, and streams)."""
        self._state = None
        for u in self.units:
            u._state = None
            u._dstate = None
        for ws in self.streams:
            y = ws.state.copy()
            ws._state = None
            ws._dstate = None
            ws.state = y*0.0
            ws.dstate = y*0.0

    def mock_convergence(
            self, 
            convergence_rate=10, 
        ):
        self.track_convergence()
        self._setup()
        self._reset_iter()
        algorithm = self.algorithm 
        if algorithm == 'Sequential modular':
            method = self._mock_run_sequential
        elif algorithm == 'Phenomena based':
            method = self._mock_run_phenomena
        elif algorithm == 'Phenomena modular':
            method = self._mock_run_phenomena_modular
        else:
            raise RuntimeError(f'unknown algorithm {algorithm!r}')
        variable_errors = {i: 100 for i in self.variable_nodes}
        while sum(variable_errors.values()): method(variable_errors, convergence_rate)
        
    def _mock_run_phenomena(self, variable_errors, convergence_rate):
        for i in ('phenomenode', 'energy', 'material'):
            self._mock_collect_variables(i, variable_errors, convergence_rate)
        
    def _mock_run_sequential(self, variable_errors, convergence_rate):
        subgraphs = ('phenomenode', 'energy', 'material')
        for i in self._path:
            if hasattr(i, 'stages'):
                for factor in (0.75, 0.25):
                    inner_rate = factor * convergence_rate
                    for subgraph in subgraphs:
                        self._mock_collect_variables(
                            (i.ID, subgraph), 
                            variable_errors, 
                            inner_rate
                        )
            else:
                for subgraph in subgraphs:
                    self._mock_collect_variables(
                        (i.ID, subgraph), 
                        variable_errors, 
                        convergence_rate
                    )
                    
    def _mock_run_phenomena_modular(self, variable_errors, convergence_rate):
        for i in self._path:
            if isinstance(i, (bst.Mixer, bst.Splitter)): continue
            if hasattr(i, 'stages'):
                for factor in (3 / 4, 3 / 16, 1 / 16):
                    inner_rate = factor * convergence_rate
                    self._mock_collect_variables(
                        (i.ID, 'phenomenode'), 
                        variable_errors, 
                        inner_rate
                    )
            else:
                self._mock_collect_variables(
                    (i.ID, 'phenomenode'), 
                    variable_errors, 
                    convergence_rate
                )
            for subgraph in ('energy', 'material'):
                self._mock_collect_variables(
                    subgraph, 
                    variable_errors, 
                    convergence_rate
                )
            
    def simulate(self, update_configuration: Optional[bool]=None, units=None, 
                 design_and_cost=None, **kwargs):
        """
        If system is dynamic, run the system dynamically. Otherwise, converge 
        the path of unit operations to steady state. After running/converging 
        the system, size and cost all unit operations.
        
        Parameters
        ----------
        **kwargs : 
            Additional parameters for :func:`dynamic_run` (if dynamic) or 
            :func:`converge` (if steady state).
        update_configuration :
            Whether to update system configuration if unit connections have
            changed. 
        units : 
            Unit operations of the system. If given, original unit operations of 
            the system will be replaced.
            
        """
        with self.flowsheet:
            specifications = self._specifications
            if specifications and not self._running_specifications:
                self._running_specifications = True
                try:
                    if self.simulate_after_specifications: 
                        for ss in specifications: ss()
                        outputs = self.simulate(
                            update_configuration, units, design_and_cost, **kwargs
                        )
                    else:
                        # Save simulation arguments and outputs when running specifications that simulate
                        self._simulation_default_arguments = dict(
                            update_configuration=update_configuration, 
                            design_and_cost=design_and_cost,
                            units=units,
                            kwargs=kwargs,
                        )
                        for ss in specifications: ss()
                        try: 
                            outputs = self._simulation_outputs
                        except AttributeError:
                            outputs = None
                        else:
                            del self._simulation_outputs
                        finally:
                            del self._simulation_default_arguments
                finally:
                    self._running_specifications = False
            else:
                if hasattr(self, '_simulation_default_arguments'):
                    sda = self._simulation_default_arguments
                    if update_configuration is None and 'update_configuration' in sda:
                        update_configuration = sda['update_configuration']
                    if units is None and 'units' in sda:
                        units = sda['units']
                    if design_and_cost is None and 'design_and_cost' in sda:
                        design_and_cost = sda['design_and_cost']
                    if kwargs is None and 'kwargs' in sda:
                        kwargs = sda['kwargs']
                if design_and_cost is None: design_and_cost = True
                if update_configuration is None: update_configuration = False
                self._setup(update_configuration, units)
                if self.isdynamic: 
                    outputs = self.dynamic_run(**kwargs)
                    if design_and_cost: self._summary()
                else:
                    try:
                        outputs = self.converge(**kwargs)
                        if design_and_cost: self._summary()
                    except Exception as error:
                        if update_configuration: raise error # Avoid infinite loop
                        new_connections = [i.get_connection() for i in self.streams]
                        if self._connections != new_connections:
                            # Connections has been updated within simulation.
                            outputs = self.simulate(
                                update_configuration=True, 
                                design_and_cost=design_and_cost,
                                **kwargs
                            )
                        else:
                            raise error
                    else:
                        if (not update_configuration # Avoid infinite loop
                            and self._connections != [i.get_connection() for i in self.streams]):
                            # Connections has been updated within simulation.
                            outputs = self.simulate(
                                update_configuration=True,
                                design_and_cost=design_and_cost,
                                **kwargs
                            )
                self._simulation_outputs = outputs
            return outputs

    def dynamic_run(self, **dynsim_kwargs):
        """
        Run system dynamically without accounting for
        the cost or environmental impacts of unit operations.
        
        Parameters
        ----------
        **dynsim_kwargs : 
            Dynamic simulation keyword arguments, could include:

                t_span : tuple[float, float]
                    Interval of integration (t0, tf).
                    The solver starts with t=t0 and integrates until it reaches t=tf.
                t_eval : iterable(float)
                    The time points where status will be saved.
                state_reset_hook: str|Callable
                    Hook function to reset the cache state between simulations
                    for dynamic systems).
                    Can be "reset_cache" or "clear_state" to call `System.reset_cache`
                    or `System.clear_state`, or None to avoiding resetting.
                export_state_to: str
                    If provided with a path, will save the simulated states over time to the given path,
                    supported extensions are ".xlsx", ".xls", "csv", and "tsv".
                sample_id : str
                    ID of the samples to run (for results exporting).
                print_msg : bool
                    Whether to print returned message from scipy.
                print_t : bool
                    Whether to print integration time in the console,
                    usually used for debugging.
                solve_ivp_kwargs
                    All remaining keyword arguments will be passed to ``solve_ivp``.
        
        See Also
        --------
        `scipy.integrate.solve_ivp <https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html>`_
        
        """
        dk = self.dynsim_kwargs
        dk.update(dynsim_kwargs)
        dk_cp = dk.copy()
        t_eval = dk_cp.pop('t_eval', None)
        state_reset_hook = dk_cp.pop('state_reset_hook', None)
        export_state_to = dk_cp.pop('export_state_to', '')
        sample_id = dk_cp.pop('sample_id', '')
        print_msg = dk_cp.pop('print_msg', False)
        print_t = dk_cp.pop('print_t', False)
        dk_cp.pop('y0', None) # will be updated later
        # Reset state, if needed
        if state_reset_hook:
            if isinstance(state_reset_hook, str):
                f = getattr(self, state_reset_hook)
            else:
                f = state_reset_hook # assume to be a function
            f()
        else: 
            for u in self.units:
                if not hasattr(u, '_state'): u._init_dynamic()
        # Load initial states
        self.converge()
        y0, idx, nr = self._load_state()
        dk['y0'] = y0
        # Integrate
        self.dynsim_kwargs['print_t'] = print_t # self.dynsim_kwargs might be reset by `state_reset_hook`
        self.scope.sol = sol = solve_ivp(fun=self.DAE, y0=y0, **dk_cp)
        if print_msg:
            if sol.status == 0:
                print('Simulation completed.')
            else: print(sol.message)
        self._write_state()          
        # Write states to file
        if export_state_to:
            try: file, ext = export_state_to.rsplit('.', 1)
            except ValueError: 
                file = export_state_to
                ext = 'npy'
            if sample_id != '':
                if ext != 'npy':
                    warn(f'file extension .{ext} ignored. Time-series data of '
                          f"{self.ID}'s tracked variables saved as a .npy file.")
                path = f'{file}_{sample_id}.npy'
            else: path = f'{file}.{ext}'
            self.scope.export(path=path, t_eval=t_eval)

    # User definitions
    
    def define_process_impact(self, key: str, name: str, basis: str, inventory: Callable, CF: float):
        """
        Define a process impact.

        Parameters
        ----------
        key :
            Impact indicator key.
        name :
            Name of process impact.
        basis :
            Functional unit for the characterization factor.
        inventory : 
            Should return the annualized (not hourly) inventory flow rate 
            given no parameters. Units should be in basis / yr
        CF : 
            Characterization factor in impact / basis.

        """
        item = ProcessImpactItem(name, basis, inventory, CF)
        try:
            if key in self.process_impact_items:
                items = self.process_impact_items[key]
                name = item.name
                for i in items:
                    if i.name == name:
                        raise ValueError(
                            f"item with name '{name}' already defined"
                        )
                items.append(item)
            else:
                self.process_impact_items[key] = [item]
        except AttributeError:
            self.process_impact_items = {key: [item]}

    # Convenience methods

    @property
    def heat_utilities(self) -> tuple[HeatUtility, ...]:
        """The sum of all heat utilities in the system by agent."""
        return HeatUtility.sum_by_agent(get_heat_utilities(self.cost_units))

    @property
    def power_utility(self) -> PowerUtility:
        """Sum of all power utilities in the system."""
        return PowerUtility.sum(get_power_utilities(self.cost_units))

    def get_inlet_cost_flows(self):
        """
        Return a dictionary with flow rates for inlet streams with fees/credits/utilities.
        """
        dct = {}
        for unit in self.units:
            for name, flow in unit.get_inlet_cost_flows().items():
                if flow:
                    if name in dct:
                        dct[name] += flow
                    else:
                        dct[name] = flow
        return dct

    def get_outlet_revenue_flows(self):
        """
        Return a dictionary with flow rates for outlet streams with fees/credits/utilities.
        """
        dct = {}
        for unit in self.units:
            for name, flow in unit.get_outlet_revenue_flows().items():
                if flow:
                    if name in dct:
                        dct[name] += flow
                    else:
                        dct[name] = flow
        return dct

    def get_inlet_flow(self, units: str, key: Optional[Sequence[str]|str]=None):
        """
        Return total flow across all inlets per year.

        Parameters
        ----------
        units : 
            Material units of measure (e.g., 'kg', 'gal', 'kmol').
        key : 
            Chemical identifiers. If none given, the sum of all chemicals returned

        Examples
        --------
        >>> from biosteam import Stream, Mixer, Splitter, settings, main_flowsheet
        >>> settings.set_thermo(['Water', 'Ethanol'])
        >>> main_flowsheet.clear()
        >>> S1 = Splitter('S1', Stream(Ethanol=10, units='ton/hr'), split=0.1)
        >>> M1 = Mixer('M1', ins=[Stream(Water=10, units='ton/hr'), S1-0])
        >>> sys = main_flowsheet.create_system(operating_hours=330*24)
        >>> sys.get_inlet_flow('Mton') # Sum of all chemicals
        0.1584
        >>> sys.get_inlet_flow('Mton', 'Water') # Just water
        0.0792

        """
        units += '/hr'
        if key:
            return self.operating_hours * sum([i.get_flow(units, key) for i in utils.feeds_from_units(self.units)])
        else:
            return self.operating_hours * sum([i.get_total_flow(units) for i in utils.feeds_from_units(self.units)])

    def get_outlet_flow(self, units: str, key: Optional[Sequence[str]|str]=None):
        """
        Return total flow across all outlets per year.

        Parameters
        ----------
        units : 
            Material units of measure (e.g., 'kg', 'gal', 'kmol').
        key : 
            Chemical identifiers. If none given, the sum of all chemicals returned

        Examples
        --------
        >>> from biosteam import Stream, Mixer, Splitter, settings, main_flowsheet
        >>> settings.set_thermo(['Water', 'Ethanol'])
        >>> main_flowsheet.clear()
        >>> S1 = Splitter('S1', Stream(Ethanol=10, units='ton/hr'), split=0.1)
        >>> M1 = Mixer('M1', ins=[Stream(Water=10, units='ton/hr'), S1-0])
        >>> sys = main_flowsheet.create_system(operating_hours=330*24)
        >>> sys.simulate()
        >>> sys.get_outlet_flow('Mton') # Sum of all chemicals
        0.1584
        >>> sys.get_outlet_flow('Mton', 'Water') # Just water
        0.0792

        """
        units += '/hr'
        if key:
            return self.operating_hours * sum([i.get_flow(units, key) for i in utils.products_from_units(self.units)])
        else:
            return self.operating_hours * sum([i.get_total_flow(units) for i in utils.products_from_units(self.units)])
    
    def get_mass_flow(self, stream: Stream):
        """Return the mass flow rate of a stream [kg/yr]."""
        return stream.F_mass * self.operating_hours
    
    def get_market_value(self, stream: Stream):
        """Return the market value of a stream [USD/yr]."""
        return stream.cost * self.operating_hours

    def has_market_value(self, stream: Stream):
        """Return whether stream has a market value."""
        return bool(stream.price) and not stream.isempty()
        
    def get_property(self, stream: Stream, name: str, units=None):
        """Return the annualized property of a stream."""
        return stream.get_property(name, units) * self.operating_hours

    def _price2cost(self, stream):
        """Get factor to convert stream price to cost."""
        F_mass = stream.F_mass
        if not F_mass: warn(RuntimeWarning(f"stream '{stream}' is empty"))
        price2cost = F_mass * self.operating_hours
        if stream.sink and not stream.source:
            return - price2cost
        elif stream.source:
            return price2cost
        else:
            raise ValueError("stream must be either a feed or a product")

    def get_net_heat_utility_impact(self, 
            agent: bst.UtilityAgent, key: str,
            heat_utilities: Optional[tuple[bst.HeatUtility]]=None,
            displace=True
        ):
        if isinstance(agent, str): 
            ID = agent
            agent = bst.HeatUtility.get_agent(ID)
        elif isinstance(agent, bst.UtilityAgent):
            ID = agent.ID
        else:
            raise TypeError(
                 "agent must be a UtilityAgent object or a string; "
                f"not a '{type(agent).__name__}' object"
            )
        CF, units = bst.HeatUtility.get_CF(ID, key)
        if CF == 0.: return 0.
        if heat_utilities is None: heat_utilities = self.heat_utilities
        for hu in heat_utilities:
            if hu.agent and hu.agent.ID == ID:
                if not displace and hu.flow < 0: continue
                if units == 'kg':
                    return hu.flow * CF * agent.MW * self.operating_hours
                elif units == 'kmol':
                    return hu.flow * CF * self.operating_hours
                elif units == 'kJ':
                    return hu.duty * CF * self.operating_hours
                else:
                    raise RuntimeError("unknown error")
        return 0.
    
    def get_net_electricity_impact(self, key, displace=True):
        power_utility = self.power_utility
        if not displace and power_utility.power < 0: return 0.
        try:
            return self.power_utility.get_impact(key) * self.operating_hours
        except KeyError:
            return 0.

    def get_net_utility_impact(self, key, displace=True):
        agents = (*bst.HeatUtility.cooling_agents,
                  *bst.HeatUtility.heating_agents)
        heat_utilities = self.heat_utilities
        return sum([self.get_net_heat_utility_impact(i, key, heat_utilities, displace) for i in agents]) + self.get_net_electricity_impact(key, displace)
    
    def get_total_feeds_impact(self, key):
        """
        Return the total annual impact of all feeds given 
        the impact indicator key.
        
        """
        return sum([s.F_mass * s.characterization_factors[key] for s in self.feeds
                    if key in s.characterization_factors]) * self.operating_hours
    
    def get_total_products_impact(self, key):
        """
        Return the total annual impact of all products given 
        the impact indicator key.
        
        """
        return sum([s.F_mass * s.characterization_factors[key] for s in self.products
                    if key in s.characterization_factors]) * self.operating_hours
    
    def get_material_impact(self, stream, key):
        """
        Return the annual material impact given the stream and the 
        impact indicator key.
        
        """
        return stream.get_impact(key) * self.operating_hours
    
    def get_total_input_impact(self, key):
        """
        Return total annual impact of inputs given the impact indicator key.
        
        """
        heat_utilities = self.heat_utilities
        power_utility = self.power_utility
        impact = self.get_total_feeds_impact(key) / self.operating_hours
        for hu in heat_utilities:
            if hu.flow > 0.: impact += hu.get_impact(key)
        if power_utility.rate > 0.:
            impact += power_utility.get_impact(key)
        return impact * self.operating_hours
    
    def get_process_impact(self, key, displace=True):
        """
        Return the annual process impact given the impact indicator key.
        
        """
        try:
            process_impact_items = self.process_impact_items
        except:
            return 0.
        else:
            if key in process_impact_items:
                if displace:
                    return sum([item.impact() for item in process_impact_items[key]])
                else:
                    return sum([flow * item.CF for item in process_impact_items[key] if (flow:=item.inventory()) > 0])
            else:
                return 0.
    
    def get_net_impact(self, key, displace=True):
        """
        Return net annual impact, given the impact indicator key. If displace 
        is True, impacts from coproducts will substracted. 
        
        """
        if displace:
            return (
                self.get_total_feeds_impact(key)
                + self.get_process_impact(key)
                + self.get_net_utility_impact(key)
                - self.get_total_products_impact(key)
            )
        else:
            return (
                self.get_total_feeds_impact(key)
                + self.get_process_impact(key, displace=displace)
                + self.get_net_utility_impact(key, displace=displace)
            )
    
    def get_property_allocated_impact(self, key, name, basis, ignored=None, products=None):
        if ignored is None: ignored = frozenset()
        total_property = 0.
        heat_utilities = self.heat_utilities
        power_utility = self.power_utility
        operating_hours = self.operating_hours
        units = None if basis is None else basis + '/hr'
        if name in bst.allocation_properties: name += '-allocation'
        if hasattr(bst.PowerUtility, name):
            if power_utility.rate < 0.:
                total_property += power_utility.get_property(name, units) * operating_hours
        if hasattr(bst.HeatUtility, name):
            for hu in heat_utilities:
                if hu.flow < 0.: total_property += hu.get_property(name, units) * operating_hours
        if hasattr(bst.Stream, name):
            if products is None: products = self.products
            for stream in products:
                if stream in ignored: continue
                total_property += self.get_property(stream, name, units)
        impact = self.get_total_feeds_impact(key)
        for hu in heat_utilities:
            if hu.flow > 0.: impact += hu.get_impact(key) * operating_hours
        if power_utility.rate > 0.:
            impact += power_utility.get_impact(key) * operating_hours
        impact += self.get_process_impact(key)
        return impact / total_property
    
    def get_property_allocation_factors(self, name, basis=None, groups=(), ignored=None, products=None):
        if ignored is None: ignored = frozenset()
        heat_utilities = self.heat_utilities
        power_utility = self.power_utility
        operating_hours = self.operating_hours
        units = None if basis is None else basis + '/hr'
        properties = {}
        if name in bst.allocation_properties: name += '-allocation'
        if isinstance(groups, str): groups = [groups]
        if hasattr(bst.PowerUtility, name):
            if power_utility.rate < 0.:
                value = power_utility.get_property(name, units)
                set_impact_value(properties, 'Electricity', value * operating_hours, groups)
        if hasattr(bst.HeatUtility, name):
            for hu in heat_utilities:
                if hu.flow < 0.: 
                    value = hu.get_property(name, units)
                    set_impact_value(properties, hu.agent.ID, value * operating_hours, groups)
        if hasattr(bst.Stream, name):
            if products is None: products = self.products
            for stream in products:
                if stream in ignored: continue
                value = self.get_property(stream, name, units)
                set_impact_value(properties, stream.ID, value, groups)
        total_property = sum(properties.values())
        allocation_factors = {i: j / total_property for i, j in properties.items()}
        return allocation_factors
    
    def get_displacement_allocation_factors(self, main_products, key, groups=()):
        heat_utilities = self.heat_utilities
        power_utility = self.power_utility
        allocation_factors = {}
        isa = isinstance
        if isa(main_products, bst.Stream): main_products = [main_products]
        CFs_original = []
        main_products = [i for i in main_products if not i.isempty()]
        for main_product in main_products:
            if isa(main_product, bst.Stream):
                CFs_original.append(main_product.get_CF(key))
            elif main_product == 'Electricity':
                main_product = power_utility
                CFs_original.append(main_product.get_CF(key, consumption=False))
            else:
                raise NotImplementedError(f"main product '{main_product}' is not yet an option for this method")
        items = main_products.copy()
        for stream in self.products:
            if key in stream.characterization_factors:
                items.append(stream)
        for hu in heat_utilities:
            if hu.flow > 0. and (hu.agent.ID, key) in hu.characterization_factors:
                items.append(hu)
        if power_utility.rate < 0. and key in power_utility.characterization_factors:
            items.append(power_utility)
        for item in items:
            total_input_impact = self.get_total_input_impact(key)
            net_impact = self.get_net_impact(key)
            if isa(item, bst.Stream):
                item_impact = self.get_material_impact(item, key)
            else:
                item_impact = -1 * item.get_impact(key) * self.operating_hours
            displaced_impact = net_impact + item_impact
            if isa(item, (bst.Stream, bst.HeatUtility)):
                set_impact_value(allocation_factors, item.ID, displaced_impact / total_input_impact, groups)
            elif isa(item, bst.PowerUtility):
                set_impact_value(allocation_factors, 'Electricity', displaced_impact / total_input_impact, groups)
            else:
                raise RuntimeError('unknown error')
            if item in main_products:
                if isa(item, bst.Stream):
                    item.set_CF(key, displaced_impact / self.get_mass_flow(item))
                elif item == 'Electricity':
                    item.set_CF(key, displaced_impact / (item.rate * self.operating_hours))
        for i, j in zip(main_products, CFs_original): i.set_CF(key, j)
        total = sum(allocation_factors.values())
        return {i: j / total for i, j in allocation_factors.items()}
    
    def displacement_allocation_table(self, key, product):
        return bst.report.lca_displacement_allocation_table(
            [self], key, 
            product if isinstance(product, list) else [product],
        )
    
    @property
    def sales(self) -> float:
        """Annual sales revenue [USD/yr]."""
        return self.operating_hours * (
            sum([s.cost for s in self.products if s.price]) 
            + sum([i._outlet_revenue for i in self.cost_units])
        )
    @property
    def material_cost(self) -> float:
        """Annual material cost [USD/yr]."""
        return self.operating_hours * (
            sum([s.cost for s in self.feeds if s.price])
            + sum([i._inlet_cost for i in self.cost_units])
        )
    @property
    def utility_cost(self) -> float:
        """Total utility cost [USD/yr]."""
        return sum([u.utility_cost for u in self.cost_units]) * self.operating_hours
    @property
    def purchase_cost(self) -> float:
        """Total purchase cost [USD]."""
        return sum([u.purchase_cost for u in self.cost_units])
    @property
    def installed_equipment_cost(self) -> float:
        """Total installed cost [USD]."""
        lang_factor = self.lang_factor
        if lang_factor:
            return sum([u.purchase_cost for u in self.cost_units]) * lang_factor
        else:
            return sum([u.installed_cost for u in self.cost_units])
    installed_cost = installed_equipment_cost
    def get_electricity_consumption(self):
        """Return the total electricity consumption [kWhr/yr]."""
        return self.operating_hours * self.power_utility.consumption

    def get_electricity_production(self):
        """Return the total electricity production [kWhr/yr]."""
        return self.operating_hours * self.power_utility.production
    
    def get_utility_duty(self, agent):
        """Return the total utility duty for given agent [kJ/yr]."""
        if not isinstance(agent, str): agent = agent.ID
        return self.operating_hours * sum([i.duty for i in self.heat_utilities if i.agent.ID == agent]) 
    
    def get_utility_flow(self, agent):
        """Return the total utility flow for given agent [kmol/yr]."""
        if not isinstance(agent, str): agent = agent.ID
        return self.operating_hours * sum([i.flow for i in self.heat_utilities if i.agent.ID == agent]) 
    
    def get_cooling_duty(self):
        """Return the total cooling duty [kJ/yr]."""
        return - self.operating_hours * sum([i.duty for i in self.heat_utilities if i.flow * i.duty < 0])
    
    def get_heating_duty(self):
        """Return the total heating duty [kJ/yr]."""
        return self.operating_hours * sum([i.duty for i in self.heat_utilities if i.flow * i.duty > 0])
    
    # Other
    def _to_network(self):
        """Return network that defines the system path."""
        isa = isinstance
        path = [(i._to_network() if isa(i, System) else i) for i in self._path]
        network = Network.__new__(Network)
        network.path = path
        recycle = self._recycle
        network.recycle = set(recycle) if isinstance(recycle, list) else recycle
        network.units = set(self.unit_path)
        return network

    def results(self, with_units=True):
        """
        Return a DataFrame of key results from simulation.
        """
        keys = []; addkey = keys.append
        vals = []; addval = vals.append
        stream_prices = bst.stream_prices
        all_utilities = self.heat_utilities
        power_utility = self.power_utility
        if with_units:
            if power_utility:
                addkey(('Electricity', 'Power'))
                addval(('kW', power_utility.power))
                addkey(('Electricity', 'Cost'))
                addval(('USD/hr', power_utility.cost))
            for heat_utility in HeatUtility.sum_by_agent(all_utilities):
                if heat_utility:
                    ID = heat_utility.ID.replace('_', ' ').capitalize()
                    addkey((ID, 'Duty'))
                    addval(('kJ/hr', heat_utility.duty))
                    addkey((ID, 'Flow'))
                    addval(('kmol/hr', heat_utility.flow))
                    addkey((ID, 'Cost'))
                    addval(('USD/hr', heat_utility.cost))
            for name, flow in self.get_inlet_cost_flows().items():
                ID = name + ' (inlet)'
                addkey((ID, 'Flow'))
                addval(('kg/hr', flow))
                addkey((ID, 'Cost'))
                addval(('USD/hr', flow * stream_prices[name]))
            for name, flow in self.get_outlet_revenue_flows().items():
                ID = name + ' (outlet)'
                addkey((ID, 'Flow'))
                addval(('kg/hr', flow))
                addkey((ID, 'Cost'))
                addval(('USD/hr', - flow * stream_prices[name]))
            addkey(('Total purchase cost', ''))
            addval(('USD', self.purchase_cost))
            addkey(('Installed equipment cost', ''))
            addval(('USD', self.installed_cost))
            addkey(('Utility cost', ''))
            addval(('USD/hr', sum([u.utility_cost for u in self.cost_units])))
            addkey(('Material cost', ''))
            addval(('USD/hr', sum([s.cost for s in self.feeds if s.price])))
            addkey(('Sales', ''))
            addval(('USD/hr', sum([s.cost for s in self.products if s.price])))
            if not keys: return None
            df = pd.DataFrame(vals,
                              pd.MultiIndex.from_tuples(keys),
                              ('Units', self.ID))
            df.columns.name = 'System'
            return df
        else:
            power_utility = self.power_utility
            if power_utility:
                addkey(('Electricity', 'Power'))
                addval(power_utility.power)
                addkey(('Electricity', 'Cost'))
                addval(power_utility.cost)
            for heat_utility in HeatUtility.sum_by_agent(all_utilities):
                if heat_utility:
                    ID = heat_utility.ID.replace('_', ' ').capitalize()
                    addkey((ID, 'Duty'))
                    addval(heat_utility.duty)
                    addkey((ID, 'Flow'))
                    addval(heat_utility.flow)
                    addkey((ID, 'Cost'))
                    addval(heat_utility.cost)
            for name, flow in self.get_inlet_utility_flows().items():
                ID = name + ' (inlet)'
                addkey((ID, 'Flow'))
                addval(flow)
                addkey((ID, 'Cost'))
                addval(flow * stream_prices[name])
            for name, flow in self.get_outlet_utility_flows().items():
                ID = name + ' (outlet)'
                addkey((ID, 'Flow'))
                addval(flow)
                addkey((ID, 'Cost'))
                addval(-flow * stream_prices[name])
            addkey(('Total purchase cost', ''))
            addval(self.purchase_cost)
            addkey(('Installed equipment cost', ''))
            addval(self.installed_cost)
            addkey(('Utility cost', ''))
            addval(sum([u.utility_cost for u in self.cost_units]))
            addkey(('Material cost', ''))
            addval(sum([s.cost for s in self.feeds if s.price]))
            addkey(('Sales', ''))
            addval(sum([s.cost for s in self.products if s.price]))
            if not keys: return None
            series = pd.Series(vals, pd.MultiIndex.from_tuples(keys))
            series.name = self.ID
            return series

    # Report summary
    def save_report(
            self, 
            file: Optional[str]='report.xlsx', 
            dpi: Optional[str]='900', 
            sheets=None,
            stage=False,
            **stream_properties
        ): 
        """
        Save a system report as an xlsx file.
        
        Parameters
        ----------
        file : 
            File name to save report
        dpi : 
            Resolution of the flowsheet. Defaults to '300'
        **stream_properties : str
            Additional stream properties and units as key-value pairs (e.g. T='degC', flow='gpm', H='kW', etc..)
            
        """
        if sheets is None:
            sheets = {
                'Flowsheet',
                'Itemized costs',
                'Stream table',
                'Utilities',
                'Design requirements',
                'Reactions',
                # 'Specifications'
            }
        writer = pd.ExcelWriter(file)
        units = sorted(self.units, key=lambda x: x.line)
        cost_units = [i for i in units if i._design or i._cost]
        if 'Flowsheet' in sheets:
            try:
                with bst.preferences.temporary() as p:
                    p.reset()
                    p.light_mode()
                    kind = 'stage' if stage else 'thorough'
                    self.diagram(kind, file='flowsheet', dpi=str(dpi), format='png')
            except:
                diagram_completed = False
                warn(RuntimeWarning('failed to generate diagram through graphviz'), stacklevel=2)
            else:
                import PIL.Image
                try:
                    # Assume openpyxl is used
                    worksheet = writer.book.create_sheet('Flowsheet')
                    flowsheet = openpyxl.drawing.image.Image('flowsheet.png')
                    worksheet.add_image(flowsheet, anchor='A1')
                except PIL.Image.DecompressionBombError:
                    PIL.Image.MAX_IMAGE_PIXELS = int(1e9)
                    flowsheet = openpyxl.drawing.image.Image('flowsheet.png')
                    worksheet.add_image(flowsheet, anchor='A1')
                except:
                    # Assume xlsx writer is used
                    try:
                        worksheet = writer.book.add_worksheet('Flowsheet')
                    except:
                        warn("problem in saving flowsheet; please submit issue to BioSTEAM with"
                             "your current version of openpyxl and xlsx writer", RuntimeWarning)
                    worksheet.insert_image('A1', 'flowsheet.png')
                diagram_completed = True
        
        if 'Itemized costs' in sheets:
            tea = self.TEA
            if tea:
                tea = self.TEA
                cost = report.cost_table(tea)
                cost.to_excel(writer, 'Itemized costs')
                tea.get_cashflow_table().to_excel(writer, 'Cash flow')
            else:
                warn(f'Cannot find TEA object in {repr(self)}. Ignoring TEA sheets.',
                     RuntimeWarning, stacklevel=2)
        
        if 'Stream table' in sheets:
            if stage:
                with self.stage_configuration() as conf:
                    stream_tables = report.stream_tables(conf.streams, **stream_properties)
            else:
                stream_tables = report.stream_tables(self.streams, **stream_properties)
            report.tables_to_excel(stream_tables, writer, 'Stream table')
        
        if 'Utilities' in sheets:
            # Heat utility tables
            heat_utilities = report.heat_utility_tables(cost_units)
            n_row = report.tables_to_excel(heat_utilities, writer, 'Utilities')
        
            # Power utility table
            power_utility = report.power_utility_table(cost_units)
            n_row = report.tables_to_excel([power_utility], writer, 'Utilities', n_row=n_row)
            
            # Fees table
            other_utilities = report.other_utilities_table(cost_units)
            n_row = report.tables_to_excel(other_utilities, writer, 'Utilities', n_row=n_row)
        
        if 'Design requirements' in sheets:
            # General desing requirements
            results = report.unit_result_tables(cost_units)
            report.tables_to_excel(results, writer, 'Design requirements')
        
        if 'Reactions' in sheets:
            # Reaction tables
            reactions = report.unit_reaction_tables(units)
            report.tables_to_excel(reactions, writer, 'Reactions')
        
        if 'Specifications' in sheets:
            if stage:
                with self.stage_configuration() as conf:
                    specifications = [
                        u._mass_and_energy_balance_specifications_table()
                        for u in conf.stages
                    ]
            else:
                specifications = [
                    u._mass_and_energy_balance_specifications_table()
                    for u in units
                ]
            report.tables_to_excel(
                specifications, writer, 
                'Specifications'
            )
        writer.close()
        if diagram_completed: os.remove("flowsheet.png")

    # Debugging
    def _turn_on(self, mode):
        """Turn on special simulation modes like `profile` or `debug`."""
        if not isinstance(mode, str):
            raise TypeError(f"mode must be a string; not a {type(mode).__name__} object")
        mode = mode.lower()
        if mode == 'debug':
            _wrap_method = _method_debug
        elif mode == 'profile':
            _wrap_method = _method_profile
        else:
            raise ValueError(f"mode must be either 'debug' or 'profile'; not '{mode}'")
        for u in self.units:
            if u._specifications:
                u._specifications = [_wrap_method(u, i) for i in u.specification]
            else:
                u.run = _wrap_method(u, u.run)
            u._design = _wrap_method(u, u._design)
            u._cost = _wrap_method(u, u._cost)
            u._lca = _wrap_method(u, u._lca)

    def _turn_off(self):
        """Turn off special simulation modes like `profile` or `debug`."""
        for u in self.units:
            if u.specification:
                u.specification = u.specification._original
            else:
                u.run = u.run._original
            u._design = u._design._original
            u._cost = u._cost._original
            u._lca = u._lca._original

    def debug(self):
        """Simulate in debug mode. If an exception is raised, it will
        automatically enter in a breakpoint."""
        self._turn_on('debug')
        try: self.simulate()
        finally: self._turn_off()

    def profile(self):
        """
        Simulate system in profile mode and return a DataFrame object of unit
        operation simulation times.

        """
        self._turn_on('profile')
        try: self.simulate()
        finally: self._turn_off()
        units = self.units
        units.sort(key=(lambda u: u._total_excecution_time_), reverse=True)
        data = [(u.line, 1000. * u._total_excecution_time_) for u in units]
        for u in units: del u._total_excecution_time_
        return pd.DataFrame(data, index=[u.ID for u in units],
                            columns=('Unit Operation', 'Time (ms)'))

    # Representation
    def print(self, spaces=''): # pragma: no cover
        """
        Print in a format that you can use recreate the system.
        """
        print(self._stacked_info())

    def _stacked_info(self, spaces=''): # pragma: no cover
        """
        Return info with inner layers of path and facilities stacked.
        """
        info = f"{type(self).__name__}({repr(self.ID)}"
        spaces += 4 * " "
        dlim = ',\n' + spaces
        update_info = lambda new_info: dlim.join([info, new_info])
        def get_path_info(path):
            isa = isinstance
            path_info = []
            for i in path:
                if isa(i, Unit):
                    path_info.append(str(i))
                elif isa(i, System):
                    path_info.append(i._stacked_info(spaces))
                else:
                    path_info.append(str(i))
            return '[' + (dlim + " ").join(path_info) + ']'
        path_info = get_path_info(self._path)
        info = update_info(path_info)
        facilities = self._facilities
        if facilities:
            facilities_info = get_path_info(facilities)
            facilities_info = f'facilities={facilities_info}'
            info = update_info(facilities_info)
        recycle = self._recycle
        if recycle:
            recycle = self._get_recycle_info()
            info = update_info(f"recycle={recycle}")
        if self.N_runs:
            info = update_info(f"N_runs={self.N_runs}")
        info += ')'
        return info

    def _get_recycle_info(self):
        recycle = self._recycle
        if isinstance(recycle, Stream):
            recycle = recycle._source_info()
        else:
            recycle = ", ".join([i._source_info() for i in recycle])
            recycle = '{' + recycle + '}'
        return recycle

    def _ipython_display_(self):
        if bst.preferences.autodisplay: self.diagram('minimal')
        self.show()

    def _error_info(self):
        """Return information on convergence."""
        recycle = self._recycle
        if recycle:
            if self._iter == 0: return None
            s = '' if isinstance(recycle, Stream) else 's'
            return (f"\nHighest convergence error among components in recycle"
                    f"\nstream{s} {self._get_recycle_info()} after {self._iter} loops:"
                    f"\n- flow rate   {self._mol_error:.2e} kmol/hr ({self._rmol_error*100.:.2g}%)"
                    f"\n- temperature {self._T_error:.2e} K ({self._rT_error*100.:.2g}%)")
        elif self.subsystems:
            sys = max(
                self.subsystems, key=lambda i: max(
                    [i._mol_error / i.molar_tolerance, i._rmol_error / i.relative_molar_tolerance,
                     i._T_error / i.temperature_tolerance, i._rT_error / i.relative_temperature_tolerance]
                )
            )
            return sys._error_info()

    def __str__(self):
        if self.ID: return self.ID
        else: return type(self).__name__

    def __repr__(self):
        if self.ID: return f'<{type(self).__name__}: {self.ID}>'
        else: return f'<{type(self).__name__}>'

    def show(self, layout=None, T=None, P=None, flow=None, composition=None, N=None,
             IDs=None, sort=None, data=True):
        """Prints information on system."""
        print(self._info(layout, T, P, flow, composition, N, IDs, sort, data))

    def _info(self, layout, T, P, flow, composition, N, IDs, sort, data):
        """Return string with all stream specifications."""
        ins_and_outs = repr_ins_and_outs(layout, self.ins, self.outs,
                                         T, P, flow, composition, N, IDs, sort, data)
        error = self._error_info()
        if error:
            return (f"System: {self.ID}"
                    + error + '\n'
                    + ins_and_outs)
        else:
            return f"System: {self.ID}\n{ins_and_outs}"

del ignore_docking_warnings

System.register_method('aitken', flx.conditional_aitken, conditional=True)
System.register_method('wegstein', flx.conditional_wegstein, conditional=True)
System.register_method('fixedpoint', flx.conditional_fixed_point, conditional=True)
options = dict(fatol=1e-24, xatol=1e-24, xtol=1e-24, ftol=1e-24, maxiter=int(1e6))
for name in ('anderson', 'diagbroyden', 'excitingmixing', 'linearmixing', 'broyden1', 'broyden2'):
    System.register_method(name, root, method=name, options=options)
del root, name, options


# %% Working with different operation modes

class OperationModeResults:
  
    __slots__ = ('unit_capital_costs', 
                 'utility_cost', 
                 'stream_properties',
                 'feeds', 'products',
                 'heat_utilities',
                 'power_utility')
    
    def __init__(self, unit_capital_costs, stream_properties, 
                 utility_cost, feeds, products, heat_utilities, power_utility):
        self.unit_capital_costs = unit_capital_costs
        self.stream_properties = stream_properties
        self.utility_cost = utility_cost
        self.feeds = feeds
        self.products = products
        self.heat_utilities = heat_utilities
        self.power_utility = power_utility

    @property
    def material_cost(self):
        flow_rates = self.stream_properties['F_mass']
        return sum([flow_rates[i] * i.price  for i in self.feeds])

    @property
    def sales(self):
        flow_rates = self.stream_properties['F_mass']
        return sum([flow_rates[i] * i.price for i in self.products])


class OperationMode:
    __slots__ = ('__dict__',)
    def __init__(self, **data):
        self.__dict__ = data

    def simulate(self):
        """
        Simulate operation mode and return an OperationModeResults object with
        data on variable operating costs (i.e. utility and material costs) and sales.

        """
        agile_system = self.agile_system
        operation_parameters = agile_system.operation_parameters
        mode_operation_parameters = agile_system.mode_operation_parameters
        stream_properties = agile_system.stream_properties
        for name, value in self.__dict__.items():
            if name in operation_parameters: operation_parameters[name](value)
            elif name in mode_operation_parameters: mode_operation_parameters[name](value, self)
        system = self.system
        system.simulate()
        feeds = system.feeds
        products = system.products
        cost_units = system.cost_units
        operating_hours = self.operating_hours
        streams = feeds + products
        return OperationModeResults(
            {i: i.get_design_and_capital() for i in cost_units},
            {name: {stream: getattr(stream, name) for stream in streams}
             for name in stream_properties},
            operating_hours * sum([i.utility_cost for i in cost_units]),
            feeds, products, system.heat_utilities, system.power_utility
        )

    def __repr__(self):
        return f"{type(self).__name__}({repr_kwargs(self.__dict__, start='')})"


class OperationMetric:
    __slots__ = ('getter', 'value')
    
    def __init__(self, getter):
        self.getter = getter
    
    def __call__(self):
        return self.value
    
    def __repr__(self):
        return f"{type(self).__name__}(getter={self.getter})"


class AgileSystem:
    """
    Class for creating objects which may serve to retrieve
    general results from multiple operation modes in such a way that it
    represents an agile production process. When simulated, an AgileSystem
    generates results from system operation modes and compile them to
    retrieve results later.

    Parameters
    ----------
    operation_modes : list[OperationMode], optional
        Defines each mode of operation with time steps and parameter values
    operation_parameters : dict[str: function], optional
        Defines all parameters available for all operation modes.
    lang_factor : float, optional
        Lang factor for getting fixed capital investment from
        total purchase cost. If no lang factor, installed equipment costs are
        estimated using bare module factors.

    """

    __slots__ = (
        'operation_modes', 'operation_parameters', 'active_operation_mode',
        'mode_operation_parameters', 'annual_operation_metrics',
        'operation_metrics', 'unit_capital_costs', 
        'net_electricity_consumption', 'utility_cost', 
        'stream_properties', 'flow_rates', 'feeds', 'products', 'purchase_cost', 
        'installed_equipment_cost', 'heat_utilities', 'power_utility',
        'process_impact_items', 'lang_factor', '_OperationMode', 
        '_TEA', '_LCA', '_units', '_streams',
    )
    isdynamic = False
    TEA = System.TEA
    LCA = System.LCA
    define_process_impact = System.define_process_impact
    get_net_heat_utility_impact = System.get_net_heat_utility_impact
    get_net_electricity_impact = System.get_net_electricity_impact
    get_net_utility_impact = System.get_net_utility_impact
    get_total_input_impact = System.get_total_input_impact
    get_process_impact = System.get_process_impact
    get_net_impact = System.get_net_impact
    get_property_allocated_impact = System.get_property_allocated_impact
    get_property_allocation_factors = System.get_property_allocation_factors
    get_displacement_allocation_factors = System.get_displacement_allocation_factors
    get_electricity_consumption = System.get_electricity_consumption
    get_electricity_production = System.get_electricity_production
    get_utility_duty = System.get_utility_duty
    get_utility_flow = System.get_utility_flow
    get_cooling_duty = System.get_cooling_duty
    get_heating_duty = System.get_heating_duty
    rescale = System.rescale

    def __init__(self, operation_modes=None, operation_parameters=None, 
                 mode_operation_parameters=None, annual_operation_metrics=None,
                 operation_metrics=None, lang_factor=None, 
                 stream_property_names=None):
        self.operation_modes = [] if operation_modes is None else operation_modes 
        self.operation_parameters = {} if operation_parameters  is None else operation_parameters
        self.mode_operation_parameters = {} if mode_operation_parameters is None else mode_operation_parameters
        self.annual_operation_metrics = [] if annual_operation_metrics is None else annual_operation_metrics
        self.operation_metrics = [] if operation_metrics is None else operation_metrics
        self.flow_rates = flow_rates = {}
        self.stream_properties = stream_properties = {'F_mass': flow_rates}
        if stream_property_names is not None:
            for i in stream_property_names: stream_properties[i] = {}
        self.lang_factor = lang_factor
        self.heat_utilities = None
        self.power_utility = None
        self.active_operation_mode = None
        self._OperationMode = type('OperationMode', (OperationMode,), {'agile_system': self})

    def _downstream_system(self, unit):
        return self

    def get_all_recycles(self):
        return set(sum([i.system.get_all_recycles() for i in self.operation_modes], []))

    def operation_mode(self, system, operating_hours, **data):
        """
        Define and register an operation mode.

        Parameters
        ----------
        operating_hours : function
            Length of operation in hours.
        **data : str
            Name and value-pairs of operation parameters.

        """
        for s in system.streams: s.unlink()
        om = self._OperationMode(system=system, operating_hours=operating_hours, **data)
        self.operation_modes.append(om)
        return om

    def operation_parameter(self, setter=None, name=None, mode_dependent=False):
        """
        Define and register operation parameter.

        Parameters
        ----------
        setter : function
            Should set parameter.
        name : str, optional
            Name of parameter. If None, default to argument name of setter.
        mode_dependent :
            Whether the setter accepts the OperationMode object as a second argument.

        """
        if not setter: return lambda setter: self.operation_parameter(setter, name, mode_dependent)
        if not name: name, *_ = signature(setter).parameters.keys()
        if mode_dependent:
            self.mode_operation_parameters[name] = setter
        else:
            self.operation_parameters[name] = setter
        return setter

    def operation_metric(self, getter=None, annualize=False):
        """
        Return an OperationMetric object.
        
        Parameters
        ----------    
        getter : Callable(<OperationMode>)
            Should return the metric value.
        annualize : bool, optional
            Whether to multiply by operating hours and sum across operation modes.
            If True, return value of the OperationMetric is a float. Otherwise, 
            the OperationMetric returns a dictionary of results with OperationMode
            objects as keys.
        
        """
        if not getter: return lambda getter: self.operation_metric(getter, annualize)
        operation_metric = OperationMetric(getter)
        if annualize:
            self.annual_operation_metrics.append(operation_metric)
        else:
            self.operation_metrics.append(operation_metric)
        return operation_metric

    def get_market_value(self, stream):
        """Return the market value of a stream [USD/yr]."""
        return self.flow_rates[stream] * stream.price
    
    def has_market_value(self, stream):
        """Return whether stream has a market value."""
        return bool(self.flow_rates[stream] and stream.price)
    
    def get_mass_flow(self, stream):
        """Return the mass flow rate of a stream [kg/yr]."""
        return self.flow_rates[stream]
    
    def get_property(self, stream, name, units=None):
        """Return the annualized property of a stream."""
        value = self.stream_properties[name][stream] * self.operating_hours
        if units is None:
            return value
        else:
            units_dct = bst.Stream._units_of_measure
            original_units = units_dct[name]
            return original_units.convert(value, units)
    
    def get_total_feeds_impact(self, key):
        """
        Return the total annual impact of all feeds given 
        the characterization factor key.
        
        """
        flow_rates = self.flow_rates
        return sum([flow_rates[s] * s.characterization_factors[key] for s in self.feeds
                    if key in s.characterization_factors])
    
    def get_total_products_impact(self, key):
        """
        Return the total annual impact of all products given 
        the characterization factor key.
        
        """
        flow_rates = self.flow_rates
        return sum([flow_rates[s] * s.characterization_factors[key] for s in self.products
                    if key in s.characterization_factors])
    
    def get_material_impact(self, stream, key):
        """
        Return the annual material impact given the stream and the 
        characterization factor key.
        
        """
        return self.flow_rates[stream] * stream.characterization_factor[key]
    
    def _price2cost(self, stream):
        """Get factor to convert stream price to cost for cash flow in solve_price method."""
        if stream in self.flow_rates:
            F_mass = self.flow_rates[stream]
        else:
            F_mass = 0.
        if not F_mass: warn(f"stream '{stream}' is empty", category=RuntimeWarning)
        if stream in self.products:
            return F_mass
        elif stream in self.feeds:
            return - F_mass
        else:
            raise ValueError("stream must be either a feed or a product")

    @property
    def material_cost(self):
        flow_rates = self.flow_rates
        return sum([flow_rates[i] * i.price  for i in self.feeds])

    @property
    def sales(self):
        flow_rates = self.flow_rates
        return sum([flow_rates[i] * i.price for i in self.products])

    streams = System.streams

    @property
    def units(self):
        try:
            return self._units
        except:
            self._units = units = []
            past_units = set()
            for i in self.operation_modes:
                for i in i.system.unit_path:
                    if i in past_units: continue
                    units.append(i)
            return units

    @property
    def cost_units(self):
        systems = set([i.system for i in self.operation_modes])
        if len(systems) == 1:
            return systems.pop().cost_units
        else:
            units = set()
            for i in systems: units.update(i.cost_units)
            return units

    def empty_recycles(self):
        for mode in self.operation_modes:
            mode.system.empty_recycles()

    def reset_cache(self):
        for mode in self.operation_modes:
            mode.system.reset_cache()

    @property
    def operating_hours(self):
        return sum([i.operating_hours for i in self.operation_modes])
    @operating_hours.setter
    def operating_hours(self, operating_hours):
        factor = operating_hours / self.operating_hours
        for i in self.operation_modes: i.operating_hours *= factor

    def simulate(self):
        operation_modes = self.operation_modes
        operation_metrics = self.operation_metrics
        annual_operation_metrics = self.annual_operation_metrics
        N_modes = len(operation_modes)
        N_metrics = len(operation_metrics)
        N_annual_metrics = len(annual_operation_metrics)
        annual_metric_range = range(N_annual_metrics)
        metric_range = range(N_metrics)
        mode_range = range(N_modes)
        operation_mode_results = N_modes * [None]
        annual_values = [N_modes * [None] for i in annual_metric_range]
        values = [{i: None for i in operation_modes} for i in metric_range]
        total_operating_hours = self.operating_hours
        for i in mode_range:
            self.active_operation_mode = mode = operation_modes[i]
            operation_mode_results[i] = results = mode.simulate()
            for j in annual_metric_range:
                metric = annual_operation_metrics[j]
                annual_values[j][i] = metric.getter(mode) * mode.operating_hours
            for j in metric_range:
                metric = operation_metrics[j]
                values[j][mode] = metric.getter(mode)
            scale = mode.operating_hours / total_operating_hours
            for hu in results.heat_utilities: hu.scale(scale)
            results.power_utility.scale(scale)
        self.active_operation_mode = None
        for i in annual_metric_range:
            metric = annual_operation_metrics[i]
            metric.value = sum(annual_values[i])
        for i in metric_range:
            metric = operation_metrics[i]
            metric.value = values[i]
        units = set(sum([list(i.unit_capital_costs) for i in operation_mode_results], []))
        unit_modes = {i: [] for i in units}
        for results in operation_mode_results:
            for i, j in results.unit_capital_costs.items(): unit_modes[i].append(j)
        self.heat_utilities = bst.HeatUtility.sum_by_agent(sum([r.heat_utilities for r in operation_mode_results], []))
        self.power_utility = bst.PowerUtility.sum([r.power_utility for r in operation_mode_results])
        self.unit_capital_costs = {i: i.get_agile_design_and_capital(j) for i, j in unit_modes.items()}
        self.utility_cost = sum([i.utility_cost for i in operation_mode_results])
        self.feeds = list(set(sum([i.feeds for i in operation_mode_results], [])))
        self.products = list(set(sum([i.products for i in operation_mode_results], [])))
        self.purchase_cost = sum([u.purchase_cost for u in self.unit_capital_costs])
        lang_factor = self.lang_factor
        if lang_factor:
            self.installed_equipment_cost = sum([u.purchase_cost * lang_factor for u in self.unit_capital_costs.values()])
        else:
            self.installed_equipment_cost = sum([u.installed_cost for u in self.unit_capital_costs.values()])
        stream_properties = self.stream_properties
        for dct in stream_properties.values(): dct.clear()
        for i in mode_range:
            results = operation_mode_results[i]
            operating_hours = operation_modes[i].operating_hours
            for name, dct in results.stream_properties.items():
                propdct = stream_properties[name]
                for stream, property in dct.items():
                    if stream in propdct: propdct[stream] += property * operating_hours
                    else: propdct[stream] = property * operating_hours

    def __repr__(self):
        return f"{type(self).__name__}(operation_modes={self.operation_modes}, operation_parameters={self.operation_parameters}, lang_factor={self.lang_factor})"
