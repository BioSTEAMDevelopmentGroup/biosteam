# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2024, Yoel Cortes-Pena <yoelcortes@gmail.com>
#               2023-2024, Yalin Li <mailto.yalin.li@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import thermosteam as tmo
import biosteam as bst
from .._unit import streams
from biosteam.utils import as_stream, MissingStream
from biosteam.process_tools import utils
from biosteam import Unit
from typing import Optional
from inspect import signature

__all__ = ('SystemFactory', 'stream_kwargs')

def stream_kwargs(
        ID=None, flow=None, phase=None, T=None, P=None,
        units=None, price=None, total_flow=None, thermo=None, 
        characterization_factors=None, **kwargs
    ):
    """Return a dictionary of stream key word arguments as passed to 
    Stream objects."""
    if ID is not None:
        kwargs['ID'] = ID
    if flow is not None:
        kwargs['flow'] = flow
    if phase is not None:
        kwargs['phase'] = phase 
    if T is not None:
        kwargs['T'] = T
    if P is not None:
        kwargs['P'] = P
    if units is not None:
        kwargs['units'] = units
    if flow is not None:
        kwargs['flow'] = flow
    if price is not None:
        kwargs['price'] = price 
    if total_flow is not None:
        kwargs['total_flow'] = total_flow
    if thermo is not None:
        kwargs['thermo'] = thermo
    if characterization_factors is not None:
        kwargs['characterization_factors'] = characterization_factors
    return kwargs

not_chemical = {'characterization_factors', 'thermo', 'total_flow', 
                'price', 'flow', 'units', 'P', 'T', 'phase', 'flow', 'ID'}
def ignore_undefined_chemicals(kwargs, not_chemical=not_chemical):
    chemicals = tmo.settings.chemicals
    return {i: j for i, j in kwargs.items() if i in chemicals or i in not_chemical}

def get_stream(name, lst, ID):
    isa = isinstance
    isfunc = callable
    for i in lst:
        if isa(i, str):
            if ID == i: return i
        elif isa(i, dict):
           if ID == i.get('ID'): return i
        elif isfunc(i):
            if ID == i.__name__: return i
    raise ValueError(f"no {name} with ID {repr(ID)}")
    
def get_name(obj):
    if isinstance(obj, dict):
        return obj.get('ID')
    elif isinstance(obj, str):
        return obj
    elif callable(obj):
        return obj.__name__


# %% System factory

class SystemFactory:
    """
    Decorate a function to return a system from the unit operations it creates 
    when called, allowing it to default the ID, ins, and outs parameters.
    
    Parameters
    ----------
    f : Callable, optional
        Should create unit operations. `f` should have a signature of function(ins, outs, *args, **kwargs).
    ID : str, optional
        Default system name.
    ins: list[dict], optional
        List of key word arguments for initializing inlet streams.
    outs: list[dict], optional
        List of key word arguments for initializing outlet streams.
    fixed_ins_size : bool, optional
        Whether the number of inlets must match the number expected.
    fixed_outs_size : bool, optional
        Whether the number of outlets must match the number expected.
    fthermo : callable, optional
        Should return a :class:`~thermosteam.Thermo` object that may serve
        as a property package for the system. It may optionally accept an
        existing :class:`~thermosteam.Chemicals` object for compatibility with 
        the default property package (i.e., bst.settings.chemicals).
    
    Examples
    --------
    Create a heating system with just a pump and a heat exchanger:
    
    >>> from biosteam import *
    >>> @SystemFactory(
    ...     ID='heating_sys',
    ...     ins=[dict(ID='cold_stream', Water=100)],
    ...     outs=[dict(ID='hot_stream')]
    ... )
    ... def create_heating_system(ins, outs, T_out):
    ...     cold_stream, = ins
    ...     hot_stream, = outs
    ...     P1 = Pump('P1', ins=cold_stream)
    ...     H1 = HXutility('H1', ins=P1-0, outs=hot_stream, T=T_out)
    ...
    >>> create_heating_system.show()
    SystemFactory(
        f=<create_heating_system(ins, outs, T_out)>,
        ID='heating_sys',
        ins=[dict(ID='cold_stream',
                  Water=100)],
        outs=[dict(ID='hot_stream')]
    )
    >>> settings.set_thermo(['Water'], cache=True)
    >>> heating_sys = create_heating_system(T_out=350) 
    >>> heating_sys.simulate()
    >>> heating_sys.show()
    System: heating_sys
    ins...
    [0] cold_stream
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Water  100
    outs...
    [0] hot_stream
        phase: 'l', T: 350 K, P: 101325 Pa
        flow (kmol/hr): Water  100
    
    Create a mockup version, add a tank, then create the system:
    
    >>> main_flowsheet.clear() # Remove old unit operations
    >>> sys = create_heating_system(outs=[''], T_out=350, mockup=True) 
    >>> sys.show() # Mock systems have ins and outs, just like real systems
    MockSystem(
        ins=[0-P1],
        outs=[H1-0],
        units=[P1, H1]
    )
    >>> T1 = StorageTank('T1', sys-0, 'hot_stream_from_storage')
    >>> heating_sys = main_flowsheet.create_system('heating_sys')
    >>> heating_sys.simulate()
    >>> heating_sys.show() 
    System: heating_sys
    ins...
    [0] cold_stream
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Water  100
    outs...
    [0] hot_stream_from_storage
        phase: 'l', T: 350 K, P: 101325 Pa
        flow (kmol/hr): Water  100
    
    Create the system and assign unit operation IDs by area convention:
    
    >>> sys = create_heating_system(outs=[''], T_out=350, area=100, mockup=True) 
    >>> sorted(main_flowsheet.unit, key=lambda u: u.ID) # Note how previous unit operations still exist in registry
    [<HXutility: H1>, <HXutility: H101>, <Pump: P1>, <Pump: P101>, <StorageTank: T1>]
    
    To access unit operations by the original ID given in the system factory,
    you can request a unit dictionary as follows:
        
    >>> sys, udct = create_heating_system(outs=[''], T_out=350, mockup=True, area=200, udct=True)
    >>> udct['P1'] # Originally, this unit was named P1
    <Pump: P201>
    
    """
    
    def __new__(cls, f=None, ID=None, ins=None, outs=None,
                fixed_ins_size=True, fixed_outs_size=True,
                fthermo=None):
        if f:
            fsig = signature(f)
            params = list(fsig.parameters)
            if params[:2] != ['ins', 'outs']:
                raise ValueError('function must have a signature of function(ins, outs, *args, **kwargs)')
            other_params = params[2:]
            reserved_parameters = ('mockup', 'area', 'udct', 'autorename', 'operating_hours')
            for i in reserved_parameters:
                if i in other_params:
                    raise ValueError(f"function cannot accept '{i}' as an argument")
            isa = isinstance
            isfunc = callable
            self = super().__new__(cls)
            self.f = f
            self.ID = ID
            self.ins = ins = [] if ins is None else [i if isa(i, dict) or isfunc(i) else dict(ID=i) for i in ins] 
            self.outs = outs = [] if outs is None else [i if isa(i, dict) or isfunc(i) else dict(ID=i) for i in outs] 
            self.fixed_ins_size = fixed_ins_size
            self.fixed_outs_size = fixed_outs_size
            self.fthermo = fthermo
            self.__name__ = f.__name__
            self.__doc__ = f.__doc__
            self.__signature__ = fsig.replace(
                parameters=[*system_factory_parameters, *[fsig.parameters[i].replace(kind=3) for i in other_params]]
            )
            self.__annotations__ = annotations = f.__annotations__.copy()
            annotations['ins'] = annotations['outs'] = streams
            return self
        else:
            return lambda f: cls(f, ID, ins, outs, 
                                 fixed_ins_size, fixed_outs_size,
                                 fthermo)
    
    def __call__(self, ID=None, ins=None, outs=None,
            mockup=False, area=None, udct=None, 
            autorename=None, operating_hours=None,
            lang_factor=None, algorithm=None, 
            method=None, maxiter=None,
            molar_tolerance=None,
            relative_molar_tolerance=None,
            temperature_tolerance=None,
            relative_temperature_tolerance=None,
            box=False, network_priority=None,
            **kwargs
        ):
        fthermo = self.fthermo
        if fthermo: 
            fthermo_sig = signature(fthermo)
            if 'chemicals' in fthermo_sig.parameters:
                chemicals = getattr(bst.settings, 'chemicals', None)
                bst.settings.set_thermo(fthermo(chemicals=chemicals))
            elif not hasattr(bst.settings, '_thermo'):
                bst.settings.set_thermo(fthermo())
        if autorename is not None: 
            original_autorename = tmo.utils.Registry.AUTORENAME
            tmo.utils.Registry.AUTORENAME = autorename
        ins = create_streams(self.ins, ins, 'inlets', self.fixed_ins_size)
        outs = create_streams(self.outs, outs, 'outlets', self.fixed_outs_size)
        rename = area is not None
        options = dict(
            ID=ID or self.ID,
            operating_hours=operating_hours,
            lang_factor=lang_factor, algorithm=algorithm, 
            method=method, maxiter=maxiter,
            molar_tolerance=molar_tolerance,
            relative_molar_tolerance=relative_molar_tolerance,
            temperature_tolerance=temperature_tolerance,
            relative_temperature_tolerance=relative_temperature_tolerance,
        )
        if network_priority is not None: box = True
        if box:
            if mockup: raise ValueError('cannot box mockup system')
            if rename: 
                unit_registry = bst.main_flowsheet.unit
                irrelevant_units = set(unit_registry)
                unit_registry.untrack(irrelevant_units)
            elif udct:
                unit_registry = bst.main_flowsheet.unit
                irrelevant_units = set(unit_registry)
            if network_priority is None:
                module = bst.Module(ins=ins, outs=outs)
            else:
                module = bst.FacilityModule(ins=ins, outs=outs)
                module.network_priority = network_priority
            ins = tuple([module.auxin(i) for i in module.ins])
            outs = tuple([module.auxout(i) for i in module.outs])
            with bst.Flowsheet(ID), bst.System(**options) as system:
                self.f(ins, outs, **kwargs)
            module._init(system=system)
        else:        
            with (bst.MockSystem() if mockup else bst.System(**options)) as system:
                if rename: 
                    unit_registry = system.flowsheet.unit
                    irrelevant_units = set(unit_registry)
                    unit_registry.untrack(irrelevant_units)
                elif udct:
                    unit_registry = system.flowsheet.unit
                    irrelevant_units = set(unit_registry)
                self.f(ins, outs, **kwargs)
        system.load_inlet_ports(ins, {k: i for i, j in enumerate(self.ins) if (k:=get_name(j)) is not None})
        system.load_outlet_ports(outs, {k: i for i, j in enumerate(self.outs) if (k:=get_name(j)) is not None})
        if autorename is not None: tmo.utils.Registry.AUTORENAME = original_autorename
        if udct: 
            unit_dct = {}
            def add(key, unit):
                if key in unit_dct:
                    obj = unit_dct[key]
                    if isinstance(obj, list):
                        obj.append(unit)
                    else:
                        unit_dct[key] = [obj, unit]
                else:
                    unit_dct[key] = unit
            
            for ID, unit in system.flowsheet.unit.data.items():
                if unit in irrelevant_units: continue
                add(ID, unit)
                add(unit.line, unit)
        if rename: 
            unit_registry.track(irrelevant_units)
            utils.rename_units(system.units, area)
        return (system, unit_dct) if udct else system
    
    def get_inlet(self, ID):
        return get_stream('inlet', self._ins, ID)
    
    def get_outlet(self, ID):
        return get_stream('outlet', self._outs, ID)
    
    def __repr__(self):
        return f"{self.__name__}{self.__signature__}"
    
    def show(self):
        """Print decorator in nice format."""
        f = self.f
        ID = self.ID
        ins = self.ins
        outs = self.outs
        newline = '\n' + 9 * " "
        dlim = ',' + newline
        dct_dlim = "," + newline + 5 * " "
        repr_data = lambda kwargs: 'dict(' + dct_dlim.join([f"{i}={repr(j)}" for i,j in kwargs.items()]) + ')'
        repr_items = lambda items:'[' + dlim.join([repr_data(data) for data in items]) + ']'
        ins = repr_items(ins) if ins else str(ins)
        newline += " "
        dlim = ',' + newline
        dct_dlim = "," + newline + 5 * " "
        outs = repr_items(outs) if outs else str(outs)
        name = f"<{f.__name__}{signature(f)}>" if hasattr(f, '__name__') else str(f)
        print(
            f"SystemFactory(\n"
            f"    f={name},\n"
            f"    ID={repr(ID)},\n"
            f"    ins={ins},\n"
            f"    outs={outs}\n"
             ")"
        )    
        
    _ipython_display_ = show
        
def create_streams(defaults, user_streams, kind, fixed_size):
    Stream = tmo.Stream
    isfunc = callable
    isa = isinstance
    stream_types = (Stream, MissingStream)
    if user_streams is None:
        return [(kwargs() if isfunc(kwargs) else Stream(**ignore_undefined_chemicals(kwargs))) for kwargs in defaults]
    if isa(user_streams, stream_types):
        user_streams = [user_streams]
    N_defaults = len(defaults)
    N_streams = len(user_streams)
    if fixed_size and N_streams > N_defaults:
        raise ValueError(f'too many {kind} ({N_streams} given); '
                         f'number of {kind} must be {N_defaults} or less')
    streams = []
    index = 0
    for kwargs, stream in zip(defaults, user_streams):
        if not isa(stream, stream_types): 
            if isa(stream, str):
                if isfunc(kwargs): 
                    stream = kwargs()
                    stream.ID = stream
                else:
                    kwargs = kwargs.copy()
                    kwargs['ID'] = stream
                    stream = Stream(**ignore_undefined_chemicals(kwargs))
            elif stream:
                raise TypeError(
                    f"{kind} must be streams, strings, or None; "
                    f"invalid type '{type(stream).__name__}' at index {index}"
                )
            elif isfunc(kwargs):
                stream = kwargs()
            else:
                stream = Stream(**ignore_undefined_chemicals(kwargs))
        streams.append(stream)
        index += 1
    if N_streams < N_defaults:
        streams += [(kwargs() if isfunc(kwargs) else Stream(**ignore_undefined_chemicals(kwargs)))
                    for kwargs in defaults[index:]]
    elif N_streams > N_defaults:
        streams += [as_stream(i) for i in user_streams[N_defaults:]]
    return streams
    
system_factory_parameters = list(signature(SystemFactory.__call__).parameters.values())[1:-1]
