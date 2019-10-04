# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 22:23:44 2019

@author: yoelr
"""
from numba import njit
from math import log
from scipy.interpolate import interp1d
from numba.targets.registry import CPUDispatcher

class ThermoModel:
    __slots__ = ()
    
    @property
    def name(self):
        return self.evaluate.__name__
    
    def __repr__(self):
        return f"<{type(self).__name__}: {self.name}>"

asfast = lambda method: method if isinstance(method, CPUDispatcher) else njit(method)

class TDependentModel(ThermoModel):
    __slots__ = ('evaluate', 'integrate', 'integrate_over_T',
                 'differentiate', 'Tmin', 'Tmax', 'isaccurate')
    def __init__(self, evaluate, Tmin, Tmax,
                 integrate=None,
                 integrate_over_T=None,
                 differentiate=None,
                 compile=True,
                 isaccurate=True):
        if compile:
            self.evaluate = asfast(evaluate)
            if integrate: self.integrate = asfast(integrate)
            if integrate_over_T: self.integrate_over_T = asfast(integrate_over_T)
            if differentiate: self.differentiate = asfast(differentiate)
        else:
            self.evaluate = evaluate
            if integrate: self.integrate = integrate
            if integrate_over_T: self.integrate_over_T = integrate_over_T
            if differentiate: self.differentiate = differentiate
        self.Tmin = Tmin
        self.Tmax = Tmax
        self.isaccurate = isaccurate
    
    def __call__(self, T):
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
                 'isaccurate')
    def __init__(self, evaluate, Tmin, Tmax, Pmin, Pmax,
                 integrate_T=None, integrate_P=None,
                 differentiate_T=None, differentiate_P=None,
                 compile=True, isaccurate=True):        
        if compile:
            self.evaluate = asfast(evaluate)
            if integrate_T: self.integrate_T = asfast(integrate_T)
            if integrate_P: self.integrate_P = asfast(integrate_P)
            if differentiate_T: self.differentiate_T = asfast(differentiate_T)
            if differentiate_P: self.differentiate_P = asfast(differentiate_P)
        else:
            self.evaluate = evaluate
            if integrate_T: self.integrate_T = integrate_T
            if integrate_P: self.integrate_P = integrate_P
            if differentiate_T: self.differentiate_T = differentiate_T
            if differentiate_P: self.differentiate_P = differentiate_P
        self.Tmin = Tmin
        self.Tmax = Tmax
        self.Pmin = Pmin
        self.Pmax = Pmax
        self.isaccurate = isaccurate
    
    def __call__(self, T, P):
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

    
def SimpleTDependentModel(evaluate, Tmin=0, Tmax=1e6):
    def integrate(Ta, Tb): return evaluate((Tb+Ta)/2)*(Tb - Ta)
    def integrate_over_T(Ta, Tb): return evaluate((Tb+Ta)/2)*log(Tb/Ta)    
    return TDependentModel(evaluate, Tmin, Tmax, integrate=integrate,
                           integrate_over_T=integrate_over_T,
                           isaccurate=False)    

def SimpleTPDependentModel(evaluate, Tmin=0, Tmax=1e6):
    def integrate_T(Ta, Tb, P): return evaluate((Tb+Ta)/2, P)*(Tb - Ta)
    def integrate_P(Pa, Pb, T): return evaluate(T, (Pb+Pa)/2)*(Pb - Pa)
    return TDependentModel(evaluate, Tmin, Tmax,
                           integrate_T=integrate_T,
                           integrate_P=integrate_P,
                           isaccurate=False)

def ConstantTDependentModel(X, Tmin, Tmax):
    def Constant(T): return X
    def Constant_integrate(Ta, Tb): return X*(Tb - Ta)
    def Constant_integrate_over_T(Ta, Tb): return X*log(Tb/Ta)    
    return TDependentModel(Constant, Tmin, Tmax, integrate=Constant_integrate,
                           integrate_over_T=Constant_integrate_over_T,
                           isaccurate=False)

def ConstantTPDependentModel(X, Tmin, Tmax, Pmin=101000, Pmax=102000):
    def Constant(T, P): return X
    def Constant_integrate(Ta, Tb, P): return X*(Tb - Ta)   
    return TPDependentModel(Constant, Tmin, Tmax, Pmin, Pmax,
                            integrate_T=Constant_integrate,
                            isaccurate=False)

def InterpolatedTDependentModel(Ts, Ys, Tmin, Tmax, kind='cubic'):
    # Only allow linear extrapolation, but with whatever transforms are specified
    extrapolator = interp1d(Ts, Ys, fill_value='extrapolate')
    # If more than 5 property points, create a spline interpolation
    spline = interp1d(Ts, Ys, kind=kind) if len(Ts)>5 else extrapolator
    Tstart = Ts[0]
    Tend = Ts[-1]
    def interpolate(T): return extrapolator(T) if (T < Tstart or T > Tend) else spline(T)
    return TDependentModel(interpolate, Tmin, Tmax, compile=False, isaccurate=False)
    
def InterpolatedTPDependentModel(Ts, Ys, Tmin, Tmax, Pmin, Pmax, kind='cubic'):
    # Only allow linear extrapolation, but with whatever transforms are specified
    extrapolator = interp1d(Ts, Ys, fill_value='extrapolate')
    # If more than 5 property points, create a spline interpolation
    spline = interp1d(Ts, Ys, kind=kind) if len(Ts)>5 else extrapolator
    Tstart = Ts[0]
    Tend = Ts[-1]
    def interpolate(T, P): return extrapolator(T) if (T < Tstart or T > Tend) else spline(T)
    return TPDependentModel(interpolate, Tmin, Tmax, Pmin, Pmax,
                            compile=False, isaccurate=False)