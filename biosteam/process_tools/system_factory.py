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

__all__ = ('SystemFactory', 'system_defaults')

class SystemFactory:
    """
    Create a SystemFactory object that serves as a wrapper around functions 
    that return system objects and defaults the ID, ins, and outs parameters.
    
    Paramters
    ---------
    f : function(ID, ins, outs, *args, **kwargs)
        Should return a System object given the ID, inlets, outlets, and other parameters.
    ID : str
        Default system name.
    ins: list[dict]
        List of kwargs for initializing inlet streams.
    outs: list[dict]
        List of kwargs for initializing outlet streams.
    
    See Also
    --------
    system_factory
    
    """
    
    __slots__ = ('f', 'ID', 'ins', 'outs')
    
    def __init__(self, f, ID, ins, outs):
        self.f = f
        self.ID = ID
        self.ins = ins
        self.outs = outs
    
    def __call__(self, ID=None, ins=None, outs=None, *args, **kwargs):
        ins = create_streams(self.ins, ins, 'inlets')
        outs = create_streams(self.outs, outs, 'outlets')
        system = self.f(ID or self.ID, ins, outs, *args, **kwargs)
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
        print(
            f"SystemFactory(\n"
            f"    f={f},\n"
            f"    ID={repr(ID)},\n"
            f"    ins={ins},\n"
            f"    outs={outs}\n"
             ")"
        )    
        
    _ipython_display_ = show

def system_defaults(f=None, ID=None, ins=None, outs=None):
    """
    Decorate a function that returns a system when called, allowing it to
    default the ID, ins, and outs parameters.
    
    Paramters
    ---------
    f : function(ID, ins, outs, *args, **kwargs), optional
        Should return a System object given the ID, inlets, outlets, and other parameters.
    ID : str, optional
        Default system name.
    ins: list[dict], optional
        List of kwargs for initializing inlet streams.
    outs: list[dict], optional
        List of kwargs for initializing outlet streams.
    
    Examples
    --------
    >>> from biosteam import *
    >>> @system_factory(
    ...     ID='heating_sys',
    ...     ins=[dict(ID='cold_stream', Water=100)],
    ...     outs=[dict(ID='hot_stream')]
    ... )
    >>> def create_heating_system(ID, ins, outs, T_out):
    ...     cold_stream, = ins
    ...     hot_stream, = outs
    ...     P1 = Pump('P1', ins=cold_stream)
    ...     H1 = HXutility('H1', ins=P1-0, outs=hot_stream, T=T_out)
    ...     return System(ID, [P1, H1])
    ...
    >>> create_heating_system.show()
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
    if f:
        ins = ins or []
        outs = outs or []
        return SystemFactory(f, ID, ins, outs)
    else:
        return lambda f: system_defaults(f, ID, ins, outs)
        
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