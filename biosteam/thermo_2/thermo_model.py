# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 22:23:44 2019

@author: yoelr
"""
from numba import njit
from math import log
from scipy.interpolate import interp1d
from numba.targets.registry import CPUDispatcher

def asTDependentModel(evaluate, Tmin=0, Tmax=1e6, *,
                      isaccurate=True, args=(), name=None):
    return QuickTDependentModel(evaluate, Tmin, Tmax,
                                isaccurate=isaccurate, args=args, name=name)

def asTPDependentModel(evaluate, Tmin=0, Tmax=1e6, Pmin=100000, Pmax=105000, *,
                       isaccurate=True, args=(), name=None):
    return QuickTDependentModel(evaluate, Tmin, Tmax, Pmin, Pmax,
                                isaccurate=isaccurate, args=args, name=name)

def TDependentThermoFunction(function, args):
    thermofunction = lambda T: function(T, *args)
    thermofunction.__name__ = function.__name__
    return thermofunction

def TPDependentThermoFunction(function, args):
    thermofunction = lambda T, P: function(T, P, *args)
    thermofunction.__name__ = function.__name__
    return thermofunction

class ThermoModel:
    __slots__ = ('name')
    
    def __repr__(self):
        return f"<{type(self).__name__}: {self.name}>"

asTfunc = lambda method, args: TDependentThermoFunction(method, args) if args else method
asTPfunc = lambda method, args: TPDependentThermoFunction(method, args) if args else method
asfast = lambda method: method if isinstance(method, CPUDispatcher) else njit(method)

class TDependentModel(ThermoModel):
    __slots__ = ('evaluate', 'integrate', 'integrate_over_T',
                 'differentiate', 'Tmin', 'Tmax', 'isaccurate', 'name')
    def __init__(self, evaluate, Tmin, Tmax,
                 integrate=None,
                 integrate_over_T=None,
                 differentiate=None,
                 *, compile=True,
                 isaccurate=True,
                 args=(),
                 name=None):
        function_names =  ('evaluate', 'integrate', 'integrate_over_T', 'differentiate')
        all_functions = (evaluate, integrate, integrate_over_T, differentiate)
        attrs, functions = zip(*[(attr,f) for attr,f in zip(function_names, all_functions) if f])
        if compile: functions = [asfast(f) for f in functions]
        if args: functions = [asTfunc(f, args) for f in functions]
        for attr, f in zip(attrs, functions): setattr(self, attr, f)
        self.name = name or evaluate.__name__
        self.Tmin = Tmin
        self.Tmax = Tmax
        self.isaccurate = isaccurate
    
    def __call__(self, T, P=None):
        return self.evaluate(T)
    
    def indomain(self, T):
        return self.Tmin < T < self.Tmax
    
    def show(self):
        print(f"{type(self).__name__}: {self.name}\n"
              f" Tmin: {self.Tmin:.2f}\n"
              f" Tmax: {self.Tmax:.2f}")
        
    _ipython_display_ =  show


class TPDependentModel(ThermoModel):
    __slots__ = ('evaluate', 'Tmin', 'Tmax', 'Pmin', 'Pmax',
                 'integrate_T', 'integrate_P',
                 'differentiate_T', 'differentiate_P',
                 'isaccurate', 'name')
    def __init__(self, evaluate, Tmin, Tmax, Pmin, Pmax,
                 integrate_T=None, integrate_P=None,
                 differentiate_T=None, differentiate_P=None,
                 *, compile=True, isaccurate=True, args=(), name=None):        
        function_names =  ('evaluate', 'integrate_T', 'integrate_P', 'differentiate_T', 'differentiate_P')
        all_functions = (evaluate, integrate_T, integrate_P, differentiate_T, differentiate_P)
        attrs, functions = zip(*[(attr,f) for attr,f in zip(function_names, all_functions) if f])
        if compile: functions = [asfast(f) for f in functions]
        if args: functions = [asTfunc(f, args) for f in functions]
        for attr, f in zip(attrs, functions): setattr(self, attr, f)
        self.Tmin = Tmin
        self.Tmax = Tmax
        self.Pmin = Pmin
        self.Pmax = Pmax
        self.isaccurate = isaccurate
        self.name = name or evaluate.__name__
    
    def __call__(self, T, P=101325.):
        return self.evaluate(T, P)
    
    def indomain(self, T, P):
        return (self.Tmin < T < self.Tmax) and (self.Pmin < P < self.Pmax)

    def show(self):
        print(f"{type(self).__name__}: {self.name}\n"
              f" Tmin: {self.Tmin:.2f} K\n"
              f" Tmax: {self.Tmax:.2f} K\n"
              f" Pmin: {self.Pmin:.5g} Pa\n"
              f" Pmax: {self.Pmax:.5g} Pa")
        
    _ipython_display_ = show

class NotImplementedModel(ThermoModel):
    __slots__ = ('handler',)
    def __init__(self, handler):
        self.handler = handler
    
    def get_active_model(self):
        handler = self.handler
        models = handler.models
        if models:
            active_model = models[0]
            handler.active_model = active_model
            return active_model
    
    def __getattr__(self, attr):
        active_model = self.get_active_model()
        if active_model:
            return active_model
        else: 
            raise NotImplementedError("no thermo model available")
    
    def __repr__(self):
        try:
            return super().__repr__()
        except:
            return f"<{type(self).__name__}>"

    
class QuickTDependentModel(TDependentModel):
    __slots__ = ()

    def __init__(self, evaluate, Tmin=0, Tmax=1e6, *, isaccurate=True, args=(), name=None):
        self.evaluate = asTfunc(asfast(evaluate), args)
        self.Tmin = Tmin
        self.Tmax = Tmax
        self.isaccurate = isaccurate
        self.name = name or evaluate.__name__
    
    def integrate(self, Ta, Tb):
        return self.evaluate((Tb+Ta)/2)*(Tb - Ta)
    
    def integrate_over_T(self, Ta, Tb): 
        return self.evaluate((Tb+Ta)/2)*log(Tb/Ta)

    def differentiate(self, T, dT=1e-12):
        return (self.evaluate(T+dT) - self.evaluate(T))/dT


class QuickTPDependentModel(TPDependentModel):
    __slots__ = ()

    def __init__(self, evaluate, Tmin=0, Tmax=1e6, Pmin=100000, Pmax=105000,
                 *, isaccurate=True, args=(), name=None):
        self.evaluate = asTfunc(asfast(evaluate), args)
        self.Tmin = Tmin
        self.Tmax = Tmax
        self.Pmin = Pmin
        self.Pmax = Pmax
        self.isaccurate = isaccurate
        self.name = name or evaluate.__name__
    
    def integrate_T(self, Ta, Tb, P):
        return self.evaluate((Tb+Ta)/2, P)*(Tb - Ta)
    
    def integrate_P(self, Pa, Pb, T):
        return self.evaluate((Pb+Pa)/2, T)*(Pb - Pa)

    def differentiate_T(self, T, P, dT=1e-12):
        return (self.evaluate(T+dT, P) - self.evaluate(T, P))/dT
    
    def differentiate_P(self, T, P, dP=1e-12):
        return (self.evaluate(T, P+dP) - self.evaluate(T, P))/dP
    

class ConstantTDependentModel(TDependentModel):
    __slots__ = ('X',)
    def __init__(self, X, Tmin, Tmax, isaccurate=False, name="Constant"):
        self.X = X
        self.Tmin = Tmin
        self.Tmax = Tmax
        self.isaccurate = isaccurate
        self.name = name
        
    def evaluate(self, T):
        return self.X
    
    def integrate(self, Ta, Tb):
        return self.X*(Tb - Ta)
    
    def differentiate(self, T, P):
        return 0
    
    def integrate_over_T(self, Ta, Tb): 
        return self.X*log(Tb/Ta)    


class ConstantTPDependentModel(TPDependentModel):
    __slots__ = ('X',)
    def __init__(self, X, Tmin, Tmax, Pmin=100000, Pmax=105000,
                 isaccurate=False, name='Constant'):
        self.X = X
        self.Tmin = Tmin
        self.Tmax = Tmax
        self.Pmin = Pmin
        self.Pmax = Pmax
        self.isaccurate = isaccurate
        self.name = name
    
    def evaluate(self, T):
        return self.X
    
    def integrate_T(self, Ta, Tb, P):
        return self.X*(Tb - Ta)
    
    def integrate_P(self, Pa, Pb, T): 
        return self.X*(Pb - Pa)
    
    def differentiate_T(self, T, P):
        return 0
    
    def differentiate_P(self, T, P):
        return 0
    

class InterpolatedTDependentModel(TDependentModel):
    __slots__ = ('T_lb', 'T_ub', 'extrapolator', 'spline')
    def __init__(self, Ts, Ys, Tmin, Tmax, kind='cubic', name='Interpolated'):
        # Only allow linear extrapolation, but with whatever transforms are specified
        self.extrapolator = interp1d(Ts, Ys, fill_value='extrapolate')
        # If more than 5 property points, create a spline interpolation
        self.spline = interp1d(Ts, Ys, kind=kind) if len(Ts)>5 else self.extrapolator
        self.Tmin = Tmin
        self.Tmax = Tmax
        self.T_lb = Ts[0]
        self.T_ub = Ts[-1]
        self.name = name
    
    def evaluate(self, T):
        return self.spline(T) if (self.T_lb <= T <= self.T_ub) else self.extrapolator(T)

    def integrate(self, T):
        return NotImplemented
        
    def differentiate(self, T):
        return NotImplemented
    
    
class TInterpolatedTPDependentModel(TPDependentModel):
    __slots__ = ('T_lb', 'T_ub', 'extrapolator', 'spline')
    def __init__(self, Ts, Ys, Tmin, Tmax, Pmin, Pmax, kind='cubic', name='Interpolated'):
        # Only allow linear extrapolation, but with whatever transforms are specified
        self.extrapolator = interp1d(Ts, Ys, fill_value='extrapolate')
        # If more than 5 property points, create a spline interpolation
        self.spline = interp1d(Ts, Ys, kind=kind) if len(Ts)>5 else self.extrapolator
        self.T_lb = Ts[0]
        self.T_ub = Ts[-1]
        self.Tmin = Tmin
        self.Tmax = Tmax
        self.Pmin = Pmin
        self.Pmax = Pmax
        self.name = name
    
    def evaluate(self, T):
        return self.spline(T) if (self.T_lb <= T <= self.T_ub) else self.extrapolator(T)

    def integrate(self, T):
        return NotImplemented
        
    def differentiate(self, T):
        return NotImplemented