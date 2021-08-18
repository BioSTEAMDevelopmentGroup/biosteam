# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from thermosteam import Stream
import biosteam as bst
from biosteam.utils import as_stream
from biosteam.process_tools import utils
from inspect import signature

__all__ = ('SystemFactory', )


# %% System factory

class SystemFactory:
    """
    Decorate a function that returns a system when called, allowing it to
    default the ID, ins, and outs parameters.
    
    Parameters
    ----------
    f : Callable, optional
        Should return a System object given the ID, inlets, outlets, and 
        other parameters. `f` should have a signature of function(ID, ins, outs, *args, **kwargs).
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
    optional_ins_index : list[int]
        Indexes of inlets that need not be connected to any unit within the system.
    optional_outs_index : list[int]
        Indexes of inlets that need not be connected to any unit within the system.
    
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
    [<HXutility: H1>,
     <HXutility: H101>,
     <Pump: P1>,
     <Pump: P101>,
     <StorageTank: T1>]
    
    To access unit operations by the original ID given in the system factory,
    you can request a unit dictionary as follows:
        
    >>> sys, udct = create_heating_system(outs=[''], T_out=350, mockup=True, area=200, udct=True)
    >>> udct['P1'] # Originally, this unit was named P1
    <Pump: P201>
    
    """
    __slots__ = ('f', 'ID', 'ins', 'outs',
                 'fixed_ins_size',
                 'fixed_outs_size',
                 'optional_ins_index',
                 'optional_outs_index')
    
    def __new__(cls, f=None, ID=None, ins=None, outs=None,
                fixed_ins_size=True, fixed_outs_size=True,
                optional_ins_index=(), optional_outs_index=()):
        if f:
            params = list(signature(f).parameters)
            if params[:2] != ['ins', 'outs']:
                raise ValueError('function must have a signature of function(ins, outs, *args, **kwargs)')
            other_params = params[3:]
            reserved_parameters = ('mockup', 'area', 'udct')
            for i in reserved_parameters:
                if i in other_params:
                    raise ValueError(f"function cannot accept '{i}' as an argument")
            self = super().__new__(cls)
            self.f = f
            self.ID = ID
            self.ins = ins or []
            self.outs = outs or []
            self.fixed_ins_size = fixed_ins_size
            self.fixed_outs_size = fixed_outs_size
            self.optional_ins_index = optional_ins_index
            self.optional_outs_index = optional_outs_index
            return self
        else:
            return lambda f: cls(f, ID, ins, outs, 
                                 fixed_ins_size, fixed_outs_size,
                                 optional_ins_index, optional_outs_index)
    
    def __call__(self, ID=None, ins=None, outs=None, mockup=False, area=None, udct=None, 
                 operating_hours=None, **kwargs):
        ins = create_streams(self.ins, ins, 'inlets', self.fixed_ins_size)
        outs = create_streams(self.outs, outs, 'outlets', self.fixed_outs_size)
        rename = area is not None
        with (bst.MockSystem() if mockup else bst.System(ID or self.ID, operating_hours=operating_hours)) as system:
            if rename: 
                unit_registry = system.flowsheet.unit
                irrelevant_units = system._irrelevant_units
                unit_registry.untrack(irrelevant_units)
            self.f(ins, outs, **kwargs)
        system.load_inlet_ports(ins, optional=[ins[i] for i in self.optional_ins_index])
        system.load_outlet_ports(outs, optional=[outs[i] for i in self.optional_outs_index])
        if rename: 
            units = system.units
            if udct: unit_dct = {i.ID: i for i in units}
            unit_registry.track(irrelevant_units)
            utils.rename_units(units, area)
            if udct: return system, unit_dct
        elif udct:
            unit_dct = {i.ID: i for i in system.units}
            return system, unit_dct
        return system
    
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
    if user_streams is None:
        return [Stream(**kwargs) for kwargs in defaults]
    isa = isinstance
    if isa(user_streams, Stream):
        user_streams = [user_streams]
    N_defaults = len(defaults)
    N_streams = len(user_streams)
    if fixed_size and N_streams > N_defaults:
        raise ValueError(f'too many {kind} ({N_streams} given); '
                         f'number of {kind} must be {N_defaults} or less')
    streams = []
    index = 0
    for kwargs, stream in zip(defaults, user_streams):
        if not isa(stream, Stream): 
            if isa(stream, str):
                kwargs = kwargs.copy()
                kwargs['ID'] = stream
                stream = Stream(**kwargs)
            elif stream:
                raise TypeError(
                    f"{kind} must be streams, strings, or None; "
                    f"invalid type '{type(stream).__name__}' at index {index}"
                )
            else:
                stream = Stream(**kwargs)
        streams.append(stream)
        index += 1
    if N_streams < N_defaults:
        streams += [Stream(**kwargs) for kwargs in defaults[index:]]
    elif N_streams > N_defaults:
        streams += [as_stream(i) for i in user_streams[N_defaults:]]
    return streams
    