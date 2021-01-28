# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from thermosteam import Stream
from biosteam import System
from inspect import signature

__all__ = ('SystemFactory', )


# %% System factory

class SystemFactory:
    """
    Decorate a function that returns a system when called, allowing it to
    default the ID, ins, and outs parameters.
    
    Paramters
    ---------
    f : Callable, optional
        Should return a System object given the ID, inlets, outlets, and 
        other parameters. `f` should have a signature of function(ID, ins, outs, *args, **kwargs).
    ID : str, optional
        Default system name.
    ins: list[dict], optional
        List of kwargs for initializing inlet streams.
    outs: list[dict], optional
        List of kwargs for initializing outlet streams.
    
    Examples
    --------
    Create a heating system with just a pump and a heat exchanger:
    
    >>> from biosteam import *
    >>> @SystemFactory(
    ...     ID='heating_sys',
    ...     ins=[dict(ID='cold_stream', Water=100)],
    ...     outs=[dict(ID='hot_stream')]
    ... )
    ... def create_heating_system(ID, ins, outs, T_out):
    ...     cold_stream, = ins
    ...     hot_stream, = outs
    ...     P1 = Pump('P1', ins=cold_stream)
    ...     H1 = HXutility('H1', ins=P1-0, outs=hot_stream, T=T_out)
    ...
    >>> create_heating_system.show()
    SystemFactory(
        f=<create_heating_system(ID, ins, outs, T_out)>,
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
    
    """
    __slots__ = ('f', 'ID', 'ins', 'outs')
    
    def __new__(cls, f=None, ID=None, ins=None, outs=None):
        if f:
            ins = ins or []
            outs = outs or []
            self = super().__new__(cls)
            self.f = f
            self.ID = ID
            self.ins = ins
            self.outs = outs
            return self
        else:
            return lambda f: cls(f, ID, ins, outs)
    
    def __call__(self, ID=None, ins=None, outs=None, *args, **kwargs):
        ins = create_streams(self.ins, ins, 'inlets')
        outs = create_streams(self.outs, outs, 'outlets')
        if not ID: ID = self.ID
        with System(ID) as system:
            user_system = self.f(ID, ins, outs, *args, **kwargs)
            if user_system: system.copy_like(user_system)
        system.load_inlet_ports(ins)
        system.load_outlet_ports(outs)
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
        
def create_streams(defaults, user_streams, kind):
    if user_streams is None:
        return [Stream(**kwargs) for kwargs in defaults]
    isa = isinstance
    if isa(user_streams, Stream):
        user_streams = [user_streams]
    N_defaults = len(defaults)
    N_streams = len(user_streams)
    if N_streams > N_defaults:
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
    return streams