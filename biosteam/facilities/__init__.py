# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from ._chemical_capital_investment import *
from ._boiler_turbogenerator import *
from ._cooling_tower import *
from ._chilled_water_package import *
from ._process_water_center import *
from ._air_distribution_package import *
from ._blowdown_mixer import *
from ._cleaning_in_place import *
from ._refrigeration_package import *
from ._fire_water_tank import *
from .systems import *

from . import _chemical_capital_investment
from . import _blowdown_mixer
from . import _boiler_turbogenerator
from . import _cooling_tower
from . import _chilled_water_package 
from . import _process_water_center
from . import _air_distribution_package
from . import _cleaning_in_place
from . import _refrigeration_package
from . import _fire_water_tank
from . import systems

__all__ = (
    *_chemical_capital_investment.__all__,
    *_blowdown_mixer.__all__,
    *_air_distribution_package.__all__,
    *_cooling_tower.__all__,
    *_boiler_turbogenerator.__all__,
    *_chilled_water_package.__all__,
    *_process_water_center.__all__,
    *_cleaning_in_place.__all__,
    *_refrigeration_package.__all__,
    *_fire_water_tank.__all__,
    *systems.__all__,
)

# %% Lazy re-exports from the hensmith package (PEP 562)
#
# HeatExchangerNetwork and the network synthesis helpers now live in the
# hensmith package (github.com/BioSTEAMDevelopmentGroup/hensmith). They are
# re-exported lazily so that the biosteam <-> hensmith circular dependency
# is import-safe: importing biosteam never initializes hensmith, and
# hensmith can import biosteam eagerly while defining its classes. The
# names must stay out of __all__ because biosteam/__init__ star-imports
# this module during initialization, which would resolve them eagerly and
# re-create the import cycle (guarded by tests/test_hensmith_integration.py).
# biosteam/__init__ delegates its own __getattr__/__dir__ here so the two
# shims cannot drift apart.

#: Names re-exported lazily from hensmith and still fully supported
#: under the biosteam namespace.
_HENSMITH_LAZY_NAMES = ('HeatExchangerNetwork',)

#: Names re-exported lazily from hensmith as deprecated aliases, kept for
#: one release cycle after the 2026-08 migration.
_HENSMITH_DEPRECATED_NAMES = (
    'StreamLifeCycle', 'ProblemTable', 'problem_table',
    'synthesize_network', 'plot_pinch_diagram',
)

def _import_from_hensmith(name, module_globals):
    """
    Resolve `name` from hensmith on behalf of the module owning
    `module_globals`, raising AttributeError (not ImportError, which would
    leak through `hasattr` and default-valued `getattr`) when hensmith is
    unavailable. Supported names are cached into the requesting module so
    later accesses bypass __getattr__; deprecated names are not cached so
    each access warns.
    """
    try:
        import hensmith
    except ImportError as e:
        raise AttributeError(
            f"module {module_globals['__name__']!r} has no attribute {name!r}: "
            f"{name!r} moved to the hensmith package "
            "(https://github.com/BioSTEAMDevelopmentGroup/hensmith); "
            "install it with `pip install hensmith`"
        ) from e
    value = getattr(hensmith, name)
    if name in _HENSMITH_DEPRECATED_NAMES:
        import warnings
        warnings.warn(
            f"{module_globals['__name__']}.{name} is deprecated and will be "
            f"removed in a future release; import it from hensmith instead "
            f"(`from hensmith import {name}`)",
            DeprecationWarning, stacklevel=3,
        )
    else:
        module_globals[name] = value
    return value

def __getattr__(name):
    if name in _HENSMITH_LAZY_NAMES or name in _HENSMITH_DEPRECATED_NAMES:
        return _import_from_hensmith(name, globals())
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

def __dir__():
    return sorted({*globals(), *_HENSMITH_LAZY_NAMES})
