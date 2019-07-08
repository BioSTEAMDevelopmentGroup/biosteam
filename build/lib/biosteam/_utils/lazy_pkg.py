# -*- coding: utf-8 -*-
"""
Created on Sun Jul  7 20:26:38 2019

@author: yoelr
"""
import sys
import importlib
from types import ModuleType

__all__ = ('LazyPkg',)

class LazyPkg(ModuleType):
    
    def __init__(self, __name__, modules):
        super().__init__(__name__)
        dct = self.__dict__
        dct.update(sys.modules[__name__].__dict__)
        sys.modules[__name__] = self
        try:
            __all__ = self.__all__
            for i in modules:
                if i in __all__:
                    raise ValueError(f"module '{i}' already in {__name__}.__all__")
            __all__.extend(modules)
        except:
            raise ValueError("__all__ must be present in the module as a list")
    
    def __dir__(self):
        return self.__all__
    
    def __getattr__(self, name):
        if name in self.__all__:
            module = importlib.import_module('.'+name, self.__package__)
            setattr(self, name, module)
            return module
        else:
            for i in self.__all__:
                module = getattr(self, i)
                if isinstance(module, ModuleType) and hasattr(module, name):
                    attr = getattr(module, name)
                    setattr(self, name, attr)
                    return attr
        raise AttributeError(f"module {self.__name__} has not attribute '{name}'")
        
        