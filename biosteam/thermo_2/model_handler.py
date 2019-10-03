# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 23:02:53 2019

@author: yoelr
"""

class ThermoModelHandler:
    __slots__ = ('models', 'active_model')
    def __init__(self, models, active_model):
        self.models = models
        self.active_model = active_model
    
    def __repr__(self):
        return f"<{type(self).__name__}: {', '.join([i.name for i in self.models])}>"
       
    def show(self):
        active = self.active_model
        name = lambda model:(f"{model.name} -> active") if model is active else model.name
        models = ("\n ").join([name(i) for i in self.models])
        print(f"{type(self).__name__}: \n"
              f" {models}")
        
    _ipython_display_ = show
        
    
class TDependentModelHandler(ThermoModelHandler):
    __slots__ = ()
    
    @property
    def Tmin(self):
        return min([i.Tmin for i in self.models])
    @property
    def Tmax(self):
        return max([i.Tmax for i in self.models])
    
    def evaluate(self, T):
        active_model = self.active_model
        if active_model.indomain(T):
            return active_model.evaluate(T)
        else:
            other_models = [i for i in self.models if i is not active_model]
            for active_model in other_models:
                if active_model.indomain(T):
                    if active_model.isaccurate: self.active_model = active_model
                    return active_model.evaluate(T)
            raise ValueError(f"no valid model at T={T:.2f} K")
          
    __call__ = evaluate
            
    def differentiate(self, T):
        active_model = self.active_model
        if active_model.indomain(T):
            return active_model.differentiate(T)
        else:
            other_models = [i for i in self.models if i is not active_model]
            for active_model in other_models:
                if active_model.indomain(T):
                    if active_model.isaccurate: self.active_model = active_model
                    self.active_model = active_model
                    return active_model.differentiate(T)
            raise ValueError(f"no valid model at T={T:.2f} K")
            
    def integrate(self, Ta, Tb):
        integral = 0
        for model in self.models:
            if not hasattr(model, 'integrate'): continue
            lb_satisfied = Ta > model.Tmin
            ub_satisfied = Tb < model.Tmax
            if lb_satisfied:
                if ub_satisfied:
                    return integral + model.integrate(Ta, Tb)
                else:
                    Ti = model.Tmax
                    integral += model.integrate(Ta, Ti)
                    Ta = Ti
            elif ub_satisfied:
                Ti = model.Tmin
                integral += model.integrate(Ti, Tb)
                Tb = Ti
        raise ValueError(f"no valid model between T={Ta:.2f} to {Tb:.2f} K")
    
    def integrate_over_T(self, Ta, Tb):
        integral = 0
        for model in self.models:
            if not hasattr(model, 'integrate_over_T'): continue
            lb_satisfied = Ta > model.Tmin
            ub_satisfied = Tb < model.Tmax
            if lb_satisfied:
                if ub_satisfied:
                    return integral + model.integrate_over_T(Ta, Tb)
                else:
                    Ti = model.Tmax
                    integral += model.integrate_over_T(Ta, Ti)
                    Ta = Ti
            elif ub_satisfied:
                Ti = model.Tmin
                integral += model.integrate_over_T(Ti, Tb)
                Tb = Ti
        raise ValueError(f"no valid model between T={Ta:.2f} to {Tb:.2f} K")
    
    
class TPDependentModelHandler(ThermoModelHandler):
    __slots__ = ()
    
    Tmin = TDependentModelHandler.Tmin
    Tmax = TDependentModelHandler.Tmax
    @property
    def Pmin(self):
        return min([i.Tmin for i in self.models])
    @property
    def Pmax(self):
        return max([i.Tmax for i in self.models])
    
    def evaluate(self, T, P=101325.):
        active_model = self.active_model
        if active_model.indomain(T, P):
            return active_model.evaluate(T, P)
        else:
            other_models = [i for i in self.models if i is not active_model]
            for active_model in other_models:
                if active_model.indomain(T, P):
                    if active_model.isaccurate: self.active_model = active_model
                    return active_model.evaluate(T, P)
            raise ValueError(f"no valid model at T={T:.2f} K and P={P:5g} Pa")

    __call__ = evaluate

    def differentiate_T(self, T, P=101325.):
        active_model = self.active_model
        if active_model.indomain(T, P):
            return active_model.differentiate_T(T, P)
        else:
            other_models = [i for i in self.models if i is not active_model]
            for active_model in other_models:
                if active_model.indomain(T, P):
                    if active_model.isaccurate: self.active_model = active_model
                    return active_model.differentiate_T(T, P)
            raise ValueError(f"no valid model at T={T:.2f} K and P={P:5g} Pa")
            
    def differentiate_P(self, T, P=101325.):
        active_model = self.active_model
        if active_model.indomain(T, P):
            return active_model.differentiate_P(T, P)
        else:
            other_models = [i for i in self.models if i is not active_model]
            for active_model in other_models:
                if active_model.indomain(T, P):
                    if active_model.isaccurate: self.active_model = active_model
                    return active_model.differentiate_P(T, P)
            raise ValueError(f"no valid model at T={T:.2f} K and P={P:5g} Pa")

    # def integrate_T(self, Ta, Tb, P):
    #     active_model = self.active_model
    #     if (active_model.Tmin < Ta and active_model.Tmax > Tb
    #         and active_model.Pmin < P < active_model.Pmax):
    #         return active_model.integrate_T(Ta, Tb, P)
    #     else:
    #         other_models = [i for i in self.models if i is not active_model]
    #         for active_model in other_models:
    #             if (active_model.Tmin < Ta and active_model.Tmax > Tb
    #                 and active_model.Pmin < P < active_model.Pmax):
    #                 self.active_model = active_model
    #                 try: return active_model.integrate_T(Ta, Tb, P)
    #                 except: pass
    #         raise ValueError(f"no valid model between T={Ta:.2f} to {Tb:.2f} K")
            
    # def integrate_P(self, Pa, Pb, T):
    #     active_model = self.active_model
    #     if (active_model.Pmin < Pa and active_model.Pmax > Pb) and (active_model.Tmin < T < active_model.Tmax):
    #         return active_model.integrate_T(Pa, Pb, T)
    #     else:
    #         other_models = [i for i in self.models if i is not active_model]
    #         for active_model in other_models:
    #             if (active_model.Pmin < Pa and active_model.Pmax > Pb
    #                 and active_model.Tmin < T < active_model.Tmax):
    #                 self.active_model = active_model
    #                 try: return active_model.integrate_T(Pa, Pb, T)
    #                 except: pass
    #         raise ValueError(f"no valid model between P={Pa:5g} to {Pb:5g} Pa and T={T:.2f} K")

    