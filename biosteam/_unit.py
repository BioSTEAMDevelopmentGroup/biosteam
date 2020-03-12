# -*- coding: utf-8 -*-
"""
Created on Wed Dec 18 07:18:50 2019

@author: yoelr
"""

import numpy as np
import pandas as pd
from graphviz import Digraph
from ._graphics import Graphics, default_graphics
from thermosteam import Stream
from ._heat_utility import HeatUtility
from .utils import Ins, Outs, NotImplementedMethod, \
                   format_unit_line, static
from ._power_utility import PowerUtility
from ._digraph import save_digraph
from thermosteam.utils import thermo_user, registered
from thermosteam.base import UnitsOfMeasure
import biosteam as bst

__all__ = ('Unit',)


lines_with_new_graphics = ['Unit', 'Mixer', 'Static',
                           'Splitter', 'Solids separator', 'Facility']

# %% Unit Operation

@thermo_user
@registered(ticket_name='U')
class Unit:
    """Abstract parent class for Unit objects. Child objects must contain `_run`, `_design` and `_cost` methods to estimate stream outputs of a Unit and find design and cost information.  

    **Abstract class methods**
    
    _setup()
        Set stream conditions and constant data.
    _run()
        Run simulation and update output streams.
    _design()
        Add design requirements to the `design_results` dictionary.
    _cost()
        Add itemized purchse costs to the `purchase_costs` dictionary.

    **Abstract class attributes**
    
    **line=None**
        [str] Name denoting the type of Unit class. Defaults to the class
        name of the first child class.
    **BM=None** 
        [float] Bare module factor (installation factor).
    **_units={}**
        [dict] Units of measure for `design_results` dictionary.
    **_N_ins=1**
        [int] Expected number of input streams.
    **_N_outs=2**
        [int] Expected number of output streams.
    **_ins_size_is_fixed=True**
        [bool] Whether the number of streams in ins is fixed.
    **_outs_size_is_fixed=True**
        [bool] Whether the number of streams in outs is fixed.
    **_N_heat_utilities=0**
        [int] Number of heat utilities created with each instance
    **_graphics**
        [biosteam.Graphics, abstract, optional] Settings for diagram
        representation. Defaults to a box with the same number of input
        and output edges as `_N_ins` and `_N_outs`

    Parameters
    ----------
    ID='' : str, defaults to a unique ID
        A unique identification. If ID is None, unit will not be
        registered in flowsheet.
    ins=None : Iterable[Stream, or str], Stream, or str
        Input streams or IDs to initialize input streams.
        If empty, default IDs will be given. If None, defaults to missing streams.
    outs=() : Iterable[Stream, or str], Stream, or str
        Output streams or IDs to initialize output streams.
        If empty, default IDs will be given.
        If None, leave streams missing.
    thermo=None : Thermo
        Thermo object to initialize input and output streams. Defaults to
        `biosteam.settings.thermo`.
    
    Attributes
    ----------
    ins : Ins[Stream]
        Input streams.
    outs : Outs[Stream]
        Output streams.
    power_utility : PowerUtility
        Electricity rate requirements are stored here.
    heat_utilities : tuple[HeatUtility]
        Cooling and heating requirements are stored here.
    design_results : dict
        All design requirements.
    purchase_costs : dict
        Itemized purchase costs.
    thermo : Thermo
        The thermodynamic property package used by the unit.
    
    
    Examples
    --------
    :doc:`tutorial/Creating a Unit`
    
    :doc:`tutorial/Using -pipe- notation`
    
    :doc:`tutorial/Inheriting from Unit`
    
    :doc:`tutorial/Unit decorators`
    
    """ 
    
    def __init_subclass__(cls, isabstract=False, new_graphics=True):
        dct = cls.__dict__
        if new_graphics:
            if 'line' not in dct and cls.line in lines_with_new_graphics:
                # Set new graphics for default line
                cls._graphics = Graphics.box(cls._N_ins, cls._N_outs)
                cls.line = format_unit_line(cls.__name__)
            elif '_graphics' not in dct:
                # Set new graphics for specified line
                cls._graphics = Graphics.box(cls._N_ins, cls._N_outs)
        if not isabstract and not cls._run: static(cls)
         
    _is_facility_ = False
        
    ### Abstract Attributes ###
    
    # [float] Bare module factor (installation factor).
    BM = 1.
    
    # [dict] Units for construction and design results
    _units = {}
    
    # [int] Expected number of input streams
    _N_ins = 1  
    
    # [int] Expected number of output streams
    _N_outs = 2  
    
    # [bool] Whether the number of streams in ins is fixed
    _ins_size_is_fixed = True
    
    # [bool] Whether the number of streams in outs is fixed
    _outs_size_is_fixed = True
    
    # [int] number of heat utilities
    _N_heat_utilities = 0
    
    # [StreamLinkOptions] Options of linking streams
    _stream_link_options = None
    
    # [biosteam Graphics] a Graphics object for diagram representation.
    _graphics = default_graphics
    
    # [str] The general type of unit, regardless of class
    line = 'Unit'

    ### Other defaults ###
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None):
        self._load_thermo(thermo)
        self._init_ins(ins)
        self._init_outs(outs)
        self._init_utils()
        self._init_results()
        self._assert_compatible_property_package()
        self._register(ID)
    
    def _init_ins(self, ins):
        # Ins[Stream] Input streams
        self._ins = Ins(self, self._N_ins, ins, self._thermo, self._ins_size_is_fixed)
    
    def _init_outs(self, outs):
        # Outs[Stream] Output streams
        self._outs = Outs(self, self._N_outs, outs, self._thermo, self._outs_size_is_fixed)
    
    def _init_utils(self):
        # tuple[HeatUtility] All heat utilities associated to unit
        self.heat_utilities = tuple([HeatUtility() for i in
                                     range(self._N_heat_utilities)])
        
        # [PowerUtility] Electric utility associated to unit
        self.power_utility = PowerUtility()
    
    def _init_results(self):
        # [dict] All purchase cost results in USD.
        self.purchase_costs = {}
        
        # [dict] All design results.
        self.design_results = {}
        
        # [dict] Greenhouse gas emissions
        self._GHGs = {}
    
    def _assert_compatible_property_package(self):
        chemical_IDs = self.chemicals.IDs
        streams = self.ins + self.outs
        assert all([s.chemicals.IDs == chemical_IDs for s in streams if s]), (
            "unit operation chemicals are incompatible with inlet and outlet streams; "
            "try using the `thermo` keyword argument to initialize unit operation "
            "with a compatible thermodynamic property package")
    
    def get_design_result(self, key, units):
        value = self.design_results[key]
        units_of_measure = UnitsOfMeasure(self._units[key])
        factor = units_of_measure.conversion_factor(units)
        return value * factor
    
    # Forward pipping
    def __sub__(self, other):
        """Source streams."""
        isa = isinstance
        int_types = (int, np.int)
        if isa(other, Unit):
            other._ins[:] = self._outs
            return other
        elif isa(other, int_types):
            return self.outs[other]
        elif isa(other, Stream):
            self._outs[:] = (other,)
            return self
        elif isa(other, (tuple, list, np.ndarray)):
            if all(isa(i, (int, np.int)) for i in other):
                return [self._outs[i] for i in other]
            else:
                if isa(other, Unit):
                    self._outs[:] = other._ins
                    return other
                else:
                    self._outs[:] = other
                    return self
        else:
            return other.__rsub__(self)

    def __rsub__(self, other):
        """Sink streams."""
        isa = isinstance
        int_types = (int, np.int)
        if isa(other, int_types):
            return self._ins[other]
        elif isa(other, Stream):
            self._ins[:] = (other,)
            return self
        elif isa(other, (tuple, list, np.ndarray)):
            if all(isa(i, int_types) for i in other):
                return [self._ins[i] for i in other]
            else:
                if isa(other, Unit):
                    self._ins[:] = other._outs
                    return other
                else:
                    self._ins[:] = other
                    return self
        else:
            raise ValueError(f"cannot pipe '{type(other).__name__}' object into unit")

    # Backwards pipping
    __pow__ = __sub__
    __rpow__ = __rsub__
    
    # Abstract methods
    _run = NotImplementedMethod
    _setup = NotImplementedMethod
    _design = NotImplementedMethod
    _cost = NotImplementedMethod
    def _get_design_specs(self):
        return ()
    
    def _load_stream_links(self):
        options = self._stream_link_options
        if options:
            s_in = self.ins[0]
            s_out = self.outs[0]
            s_out.link_with(s_in, options.flow, options.phase, options.TP)
    
    # Summary
    def _summary(self):
        """Calculate all results from unit run."""
        self._design()
        self._cost()
    
    @property
    def purchase_cost(self):
        """Total purchase cost [USD]."""
        return sum(self.purchase_costs.values())
    
    @property
    def installation_cost(self):
        """Installation cost [USD]."""
        return self.BM * sum(self.purchase_costs.values())
    
    @property
    def utility_cost(self):
        """Total utility cost [USD/hr]."""
        return sum([i.cost for i in self.heat_utilities]) + self.power_utility.cost

    def simulate(self):
        """Run rigourous simulation and determine all design requirements."""
        self._load_stream_links()
        self._setup()
        self._run()
        self._summary()

    def results(self, with_units=True, include_utilities=True,
                include_total_cost=True):
        """Return key results from simulation as a DataFrame if `with_units` is True or as a Series otherwise."""
        # TODO: Divide this into functions
        keys = []; addkey = keys.append
        vals = []; addval = vals.append
        if with_units:
            if include_utilities:
                power_utility = self.power_utility
                if power_utility:
                    addkey(('Power', 'Rate'))
                    addkey(('Power', 'Cost'))
                    addval(('kW', power_utility.rate))
                    addval(('USD/hr', power_utility.cost))
                for heat_utility in self.heat_utilities:
                    if heat_utility:
                        ID = heat_utility.ID.replace('_', ' ').capitalize()
                        addkey((ID, 'Duty'))
                        addkey((ID, 'Flow'))
                        addkey((ID, 'Cost'))
                        addval(('kJ/hr', heat_utility.duty))
                        addval(('kmol/hr', heat_utility.flow))
                        addval(('USD/hr', heat_utility.cost))
            units = self._units
            Cost = self.purchase_costs
            for ki, vi in self.design_results.items():
                addkey(('Design', ki))
                addval((units.get(ki, ''), vi))
            for ki, vi, ui in self._get_design_specs():
                addkey(('Design', ki))
                addval((ui, vi))
            for ki, vi in Cost.items():
                addkey(('Purchase cost', ki))
                addval(('USD', vi))
            if include_total_cost:
                addkey(('Total purchase cost', ''))
                addval(('USD', self.purchase_cost))
                addkey(('Utility cost', ''))
                addval(('USD/hr', self.utility_cost))
            if self._GHGs:
                a, b = self._totalGHG
                GHG_units =  self._GHG_units
                for ko, vo in self._GHGs.items():
                    for ki, vi in vo.items():
                        addkey((ko, ki))
                        addval((GHG_units.get(ko, ''), vi))
                a_key, b_key = GHG_units.keys()
                a_unit, b_unit = GHG_units.values()
                addkey(('Total ' + a_key, ''))
                addval((a_unit, a))
                addkey(('Total ' + b_key, ''))
                addval((b_unit, b))
            if not keys: return None
            df = pd.DataFrame(vals,
                              pd.MultiIndex.from_tuples(keys),
                              ('Units', self.ID))
            df.columns.name = self.line
            return df
        else:
            if include_utilities:
                power_utility = self.power_utility
                if power_utility:
                    addkey(('Power', 'Rate'))
                    addkey(('Power', 'Cost'))
                    addval(power_utility.rate)
                    addval(power_utility.cost)
                for heat_utility in self.heat_utilities:
                    if heat_utility:
                        if heat_utility:
                            ID = heat_utility.ID.replace('_', ' ').capitalize()
                            addkey((ID, 'Duty'))
                            addkey((ID, 'Flow'))
                            addkey((ID, 'Cost'))
                            addval(heat_utility.duty)
                            addval(heat_utility.flow)
                            addval(heat_utility.cost)
            for ki, vi in self.design_results.items():
                addkey(('Design', ki))
                addval(vi)
            for ki, vi, ui in self._get_design_specs():
                addkey(('Design', ki))
                addval((ui, vi))
            for ki, vi in self.purchase_costs.items():
                addkey(('Purchase cost', ki))
                addval(vi)    
            if self._GHGs:
                GHG_units = self._GHG_units
                for ko, vo in self._GHGs.items():
                    for ki, vi in vo.items():
                        addkey((ko, ki))
                        addval(vi)
                a, b = self._totalGHG
                a_key, b_key = GHG_units.keys()
                addkey(('Total ' + a_key, ''))
                addval(a)
                addkey(('Total ' + b_key, ''))
                addval(b)
            if include_total_cost:
                addkey(('Total purchase cost', ''))
                addval(self.purchase_cost)
                addkey(('Utility cost', ''))
                addval(self.utility_cost)
            if not keys: return None
            series = pd.Series(vals, pd.MultiIndex.from_tuples(keys))
            series.name = self.ID
            return series

    @property
    def thermo(self):
        """Thermodynamic property package"""
        return self._thermo
    @property
    def ins(self):
        """All input streams."""
        return self._ins    
    @property
    def outs(self):
        """All output streams."""
        return self._outs

    def _add_upstream_neighbors_to_set(self, set):
        """Add upsteam neighboring units to set."""
        for s in self._ins:
            u_source = s._source
            if u_source: set.add(u_source)

    def _add_downstream_neighbors_to_set(self, set):
        """Add downstream neighboring units to set."""
        for s in self._outs:
            u_sink = s._sink
            if u_sink: set.add(u_sink)

    @property
    def _downstream_units(self):
        """Return a set of all units downstream."""
        downstream_units = set()
        outer_periphery = set()
        self._add_downstream_neighbors_to_set(outer_periphery)
        inner_periphery = None
        old_length = -1
        new_length = 0
        while new_length != old_length:
            old_length = new_length
            inner_periphery = outer_periphery
            downstream_units.update(inner_periphery)
            outer_periphery = set()
            for unit in inner_periphery:
                unit._add_downstream_neighbors_to_set(outer_periphery)
            new_length = len(downstream_units)
        return downstream_units
        
    def _neighborhood(self, radius=1, upstream=True, downstream=True):
        """Return a set of all neighboring units within given radius.
        
        Parameters
        ----------
        radius : int
                 Maxium number streams between neighbors.
        downstream=True : bool, optional
            Whether to include downstream operations
        upstream=True : bool, optional
            Whether to include upstream operations
        
        """
        radius -= 1
        neighborhood = set()
        if radius < 0: return neighborhood
        if upstream:self._add_upstream_neighbors_to_set(neighborhood)
        if downstream: self._add_downstream_neighbors_to_set(neighborhood)
        direct_neighborhood = neighborhood
        for i in range(radius):
            neighbors = set()
            for unit in direct_neighborhood:
                if upstream: unit._add_upstream_neighbors_to_set(neighbors)
                if downstream: unit._add_downstream_neighbors_to_set(neighbors)
            if neighbors == direct_neighborhood: break
            direct_neighborhood = neighbors
            neighborhood.update(direct_neighborhood)
        return neighborhood

    def diagram(self, radius=0, upstream=True, downstream=True, 
                file=None, format='png', **graph_attrs):
        """Display a `Graphviz <https://pypi.org/project/graphviz/>`__ diagram of the unit and all neighboring units within given radius.
        
        Parameters
        ----------
        radius : int
                 Maxium number streams between neighbors.
        downstream=True : bool, optional
            Whether to show downstream operations
        upstream=True : bool, optional
            Whether to show upstream operations
        file : Must be one of the following:
            * [str] File name to save diagram.
            * [None] Display diagram in console.
        format : str
                 Format of file.
        
        """
        if radius > 0:
            neighborhood = self._neighborhood(radius, upstream, downstream)
            neighborhood.add(self)
            sys = bst.System('', neighborhood)
            return sys.diagram('thorough', file, format, **graph_attrs)
        
        graphics = self._graphics

        # Make a Digraph handle
        f = Digraph(name='unit', filename='unit', format=format)
        f.attr('graph', ratio='0.5', splines='normal', outputorder='edgesfirst',
               nodesep='1.1', ranksep='0.8', maxiter='1000')  # Specifications
        f.attr(rankdir='LR', **graph_attrs)  # Left to right

        # If many streams, keep streams close
        if (len(self.ins) >= 3) or (len(self.outs) >= 3):
            f.attr('graph', nodesep='0.4')

        # Initialize node arguments based on unit and make node
        type_ = graphics.node_function(self) or self.line
        name = self.ID + '\n' + type_
        f.attr('node', **self._graphics.node)
        f.node(name)

        # Set stream node attributes
        f.attr('node', shape='rarrow', fillcolor='#79dae8',
               style='filled', orientation='0', width='0.6',
               height='0.6', color='black', peripheries='1')

        # Make nodes and edges for input streams
        di = 0  # Destination position of stream
        for stream in self.ins:
            if not stream: continue
            f.node(stream.ID)  
            edge_in = self._graphics.edge_in
            if di >= len(edge_in): di = 0
            f.attr('edge', arrowtail='none', arrowhead='none',
                   tailport='e', **edge_in[di])
            f.edge(stream.ID, name)
            di += 1

        # Make nodes and edges for output streams
        oi = 0  # Origin position of stream
        for stream in self.outs:
            if not stream: continue
            f.node(stream.ID) 
            edge_out = self._graphics.edge_out  
            if oi >= len(edge_out): oi = 0
            f.attr('edge', arrowtail='none', arrowhead='none',
                   headport='w', **edge_out[oi])
            f.edge(name, stream.ID)
            oi += 1
        save_digraph(f, file, format)
    
    ### Net input and output flows ###
    
    # Molar flow rates
    @property
    def mol_in(self):
        """Molar flows going in [kmol/hr]."""
        return sum([s.mol for s in self._ins if s])
    @property
    def mol_out(self):
        """Molar flows going out [kmol/hr]."""
        return sum([s.mol for s in self._outs if s])

    @property
    def z_mol_in(self):
        """Molar fractions going in [kmol/hr]."""
        return self._mol_in/self.F_mol_in
    @property
    def z_mol_out(self):
        """Molar fractions going in."""
        return self._mol_out/self.F_mol_out

    @property
    def F_mol_in(self):
        """Net molar flow going in [kmol/hr]."""
        return sum([s.F_mol for s in self._ins if s])
    @property
    def F_mol_out(self):
        """Net molar flow going out [kmol/hr]."""
        return sum([s.F_mol for s in self._outs if s])

    # Mass flow rates
    @property
    def mass_in(self):
        """Mass flows going in [kg/hr]."""
        return sum([s.mol for s in self._ins if s]) * self._thermo.chemicals.MW
    @property
    def mass_out(self):
        """Mass flows going out [kg/hr]."""
        return sum([s.mol for s in self._outs if s]) * self._thermo.chemicals.MW

    @property
    def z_mass_in(self):
        """Mass fractions going in."""
        return self.mass_in/self.F_mass_in
    @property
    def z_mass_out(self):
        """Mass fractions going out."""
        return self.mass_out/self.F_mass_out

    @property
    def F_mass_in(self):
        """Net mass flow going in [kg/hr]."""
        return self.mass_in.sum()
    @property
    def F_mass_out(self):
        """Net mass flow going out [kg/hr]."""
        return self.mass_out.sum()

    # Volumetric flow rates
    @property
    def vol_in(self):
        """Volumetric flows going in [m3/hr]."""
        return sum([s.vol for s in self._ins if s])
    @property
    def F_vol_in(self):
        """Net volumetric flow going in [m3/hr]."""
        return sum(self.vol_in)

    @property
    def z_vol_in(self):
        """Volumetric fractions going in."""
        return self.vol_in/self.F_vol_in
    @property
    def vol_out(self):
        """Volumetric flows going out [m3/hr]."""
        return sum([s.vol for s in self._outs if s])

    @property
    def F_vol_out(self):
        """Net volumetric flow going out [m3/hr]."""
        return sum(self.vol_out)
    @property
    def z_vol_out(self):
        """Volumetric fractions going out."""
        return self.vol_out/self.F_vol_out

    # Enthalpy flow rates
    @property
    def H_in(self):
        """Enthalpy flow going in [kJ/hr]."""
        return sum([s.H for s in self._ins if s])

    @property
    def H_out(self):
        """Enthalpy flow going out [kJ/hr]."""
        return sum([s.H for s in self._outs if s])

    @property
    def Hf_in(self):
        """Enthalpy of formation flow going in [kJ/hr]."""
        return sum([s.Hf for s in self._ins if s])

    @property
    def Hf_out(self):
        """Enthalpy of formation flow going out [kJ/hr]."""
        return sum([s.Hf for s in self._outs if s])

    @property
    def Hnet(self):
        """Net enthalpy flow, including enthalpies of formation [kJ/hr]."""
        return self.H_out - self.H_in + self.Hf_out - self.Hf_in
    
    # Representation
    def _info(self, T, P, flow, composition, N):
        """Information on unit."""
        if self.ID:
            info = f'{type(self).__name__}: {self.ID}\n'
        else:
            info = f'{type(self).__name__}\n'
        info += f'ins...\n'
        i = 0
        for stream in self.ins:
            if not stream:
                info += f'[{i}] {stream}\n'
                i += 1
                continue
            stream_info = stream._info(T, P, flow, composition, N)
            unit = stream._source
            index = stream_info.index('\n')
            source_info = f'  from  {type(unit).__name__}-{unit}\n' if unit else '\n'
            info += f'[{i}] {stream.ID}' + source_info + stream_info[index+1:] + '\n'
            i += 1
        info += f'outs...\n'
        i = 0
        for stream in self.outs:
            if not stream:
                info += f'[{i}] {stream}\n'
                i += 1
                continue
            stream_info = stream._info(T, P, flow, composition, N)
            unit = stream._sink
            index = stream_info.index('\n')
            sink_info = f'  to  {type(unit).__name__}-{unit}\n' if unit else '\n'
            info += f'[{i}] {stream.ID}' + sink_info + stream_info[index+1:] + '\n'
            i += 1
        info = info.replace('\n ', '\n    ')
        return info[:-1]

    def show(self, T=None, P=None, flow=None, composition=False, N=None):
        """Prints information on unit."""
        print(self._info(T, P, flow, composition, N))
    
    def _ipython_display_(self):
        try: self.diagram()
        except: pass
        self.show()

    def __repr__(self):
        if self.ID:
            return f'<{type(self).__name__}: {self.ID}>'
        else:
            return f'<{type(self).__name__}>'

del thermo_user, registered