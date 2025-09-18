# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from __future__ import annotations
from thermosteam import settings
from thermosteam.utils import define_units_of_measure
from thermosteam.units_of_measure import (
    DisplayUnits, convert, power_utility_units_of_measure, UnitsOfMeasure
)
from typing import Optional

__all__ = ('PowerUtility',)


impact_indicator_basis = UnitsOfMeasure('kWhr')
default_price = 0.0782

@define_units_of_measure(power_utility_units_of_measure)
class PowerUtility:
    """
    Create an PowerUtility object that stores data on consumption and production
    of electricity.
            
    Notes
    -----
    The default price is 0.0782 USD/kWhr as suggested in [1]_.
    
    References
    ----------
    .. [1] Seider, W. D., Lewin,  D. R., Seader, J. D., Widagdo, S., Gani, R.,
        & Ng, M. K. (2017). Product and Process Design Principles. Wiley.
    
    Examples
    --------
    Create a PowerUtility object:
    
    >>> pu = PowerUtility()
    >>> pu
    PowerUtility(consumption=0.0, production=0.0)
    
    PowerUtility objects have `consumption` and `production` attributes
    which are updated when setting the power with the assumption that
    a positive power means no production (only consumption) and a negative
    power means no consumption (only production).
    
    >>> pu(power=-500)
    >>> pu.consumption, pu.production
    (0.0, 500.0)
    
    >>> pu(power=500.)
    >>> pu.consumption, pu.production
    (500.0, 0.0)
    
    It is possible to have both consumption and production by setting these
    attributes individually (instead of setting power)
    
    >>> pu.production = 100. 
    >>> pu.power
    400.0
    
    Notice how the power is equal to the consumption minus the production.
    
    The cost is available as a property:
    
    >>> pu.cost # USD/hr
    31.28
    
    It may be useful to print results in different units of measure:
    
    >>> pu.show(power='BTU/s')
    PowerUtility:
     consumption: 474 BTU/s
     production: 94.8 BTU/s
     power: 379 BTU/s
     cost: 31.3 USD/hr
    
    """
    __slots__ = ('consumption', 'production')
    
    #: Characterization factors for life cycle assessment [impact/kWhr] by impact key and kind (None, 'consumption', or 'production').
    characterization_factors: dict[tuple[str, str], float] = {}
    
    #: Units of measure for IPython display
    display_units: DisplayUnits = DisplayUnits(power='kW', cost='USD/hr')
    
    def __init__(self, consumption: float=0., production: float=0.):
        #: Electricity consumption [kW]
        self.consumption: float = consumption
        
        #: Electricity production [kW]
        self.production: float = production
    
    def empty(self):
        """Set consumption and production to zero."""
        self.consumption = self.production = 0.
    
    @classmethod
    def get_CF(cls, key: str, consumption: Optional[bool]=True, 
               production: Optional[bool]=True, basis: Optional[str]=None, 
               units: Optional[str]=None):
        """
        Return the life-cycle characterization factor for consumption and 
        production on a kWhr basis given the impact key.
        
        Parameters
        ----------
        key :
            Name of impact indicator.
        consumption :
            Whether to return impact indicator for electricity consumption.
        production :
            Whether to return impact indicator for electricity production.
        basis :
            Basis of characterization factor. Energy is the only valid dimension. 
            Defaults to 'kWhr'.
        units :
            Units of impact indicator. Before using this argument, the default units 
            of the impact indicator should be defined with 
            :meth:`settings.define_impact_indicator <thermosteam._settings.ProcessSettings.define_impact_indicator>`.
            Units must also be dimensionally consistent with the default units.
        
        """
        try:
            value = cls.characterization_factors[key]
        except KeyError:
            value = (0., 0.)
        if consumption:
            if not production:
                value = value[0]
        elif production:
            value = value[1]
        else:
            return None
        if units is not None:
            original_units = settings.get_impact_indicator_units(key)
            f = original_units.conversion_factor(units)
            if consumption and production:
                value = (consumption * f, 
                         production * f)
            else:
                value *= f
        if basis is not None:
            f = impact_indicator_basis.conversion_factor(basis)
            consumption /= f
            production /= f
        return value

    @classmethod
    def set_CF(cls, key: str, consumption: Optional[float]=None,
               production: Optional[float]=None, basis: Optional[str]=None, 
               units: Optional[str]=None):
        """
        Set the life-cycle characterization factors for consumption and production
        on a kWhr basis given the impact key.
        
        Parameters
        ----------
        key :
            Name of impact indicator.
        consumption :
            Impact indicator for electricity consumption.
        production :
            Impact indicator for electricity production.
        basis :
            Basis of characterization factor. Energy is the only valid dimension. 
            Defaults to 'kWhr'.
        units :
            Units of impact indicator. Before using this argument, the default units 
            of the impact indicator should be defined with 
            :meth:`settings.define_impact_indicator <thermosteam._settings.ProcessSettings.define_impact_indicator>`.
            Units must also be dimensionally consistent with the default units.
        
        """
        if consumption is None:
            if production is None:
                raise ValueError("must pass at least one of either 'consumption' or 'production'")
            consumption = production
        elif production is None:
            production = consumption
        if units is not None:
            original_units = settings.get_impact_indicator_units(key)
            f = original_units.conversion_factor(units)
            consumption /= f 
            production /= f
        if basis is not None:
            f = impact_indicator_basis.conversion_factor(basis)
            consumption *= f
            production *= f
        cls.characterization_factors[key] = (consumption, production)
    
    @classmethod
    def default_price(cls):
        """Reset price back to BioSTEAM's default."""
        cls.price: float = default_price #: Electricity price [USD/kWhr]
    
    @property
    def power(self) -> float:
        """Power requirement [kW]."""
        return self.consumption - self.production
    @power.setter
    def power(self, power: float):
        power = float(power)
        if power >= 0.:
            self.consumption = power
            self.production = 0.
        else:
            self.consumption = 0.
            self.production = -power
    rate = power # For backwards compatibility
    
    @property
    def cost(self) -> float:
        """Cost [USD/hr]"""
        return self.price * self.power
    
    def get_impact(self, key: str):
        """Return the impact [impact/hr] given characterization factor keys 
        for consumption and production. If no production key given, it defaults
        to the consumption key."""
        power = self.consumption - self.production
        try: 
            cf = self.characterization_factors[key]
        except:
            return 0.
        return (cf[0] if power > 0. else cf[1]) * power
    
    def __bool__(self):
        return bool(self.consumption or self.production)
    
    def __call__(self, power: float):
        """Set power [kW]."""
        self.power = power
    
    def copy(self):
        return self.__class__(self.consumption, self.production)
    
    def mix_from(self, power_utilities: list[PowerUtility]):
        """
        Mix in requirements of power utilities.
        
        Examples
        --------
        >>> pus = [PowerUtility(production=100),
        ...        PowerUtility(consumption=50),
        ...        PowerUtility(production=20)]
        >>> pu = PowerUtility()
        >>> pu.mix_from(pus)
        >>> print(pu)
        PowerUtility(consumption=50.0, production=120.0)
        
        """
        self.consumption = sum([i.consumption for i in power_utilities])
        self.production = sum([i.production for i in power_utilities])
    
    def copy_like(self, power_utility: PowerUtility):
        """Copy consumption and production from another power utility."""
        self.consumption = power_utility.consumption
        self.production = power_utility.production
    
    def scale(self, scale: int):
        """Scale consumption and production accordingly."""
        self.consumption *= scale
        self.production *= scale
    
    @classmethod
    def sum(cls, power_utilities: list[PowerUtility]):
        """
        Return a PowerUtility object that represents the sum of power utilities.
        
        Examples
        --------
        >>> pus = [PowerUtility(production=100),
        ...        PowerUtility(consumption=50),
        ...        PowerUtility(production=20)]
        >>> pu = PowerUtility.sum(pus)
        >>> print(pu)
        PowerUtility(consumption=50.0, production=120.0)
        
        """
        power_utility = cls()
        power_utility.mix_from(power_utilities)
        return power_utility
    
    def show(self, power: Optional[str]=None, cost: Optional[str]=None):
        # Get units of measure
        display_units = self.display_units
        power_units = power or display_units.power
        cost_units = cost or display_units.cost
        production = convert(self.production, 'kW', power_units)
        consumption = convert(self.consumption, 'kW', power_units)
        power = consumption - production
        cost = convert(self.cost, 'USD/hr', cost_units)
        print(f'{type(self).__name__}:\n'
              f'consumption: {consumption:.3g} {power_units}\n'
              f'production: {production:.3g} {power_units}\n'
              f'power: {power:.3g} {power_units}\n'
              f'cost: {cost:.3g} {cost_units}')
    _ipython_display_ = show    
    
    def __repr__(self) -> str:
        return f'{type(self).__name__}(consumption={self.consumption}, production={self.production})'

    def __add__(self, other: PowerUtility) -> PowerUtility:
        if other == 0: return self # Special case to get Python built-in sum to work
        return PowerUtility.sum([self, other])

    def __radd__(self, other: PowerUtility) -> PowerUtility:
        return self.__add__(other)
    
PowerUtility.default_price()
settings.__class__.set_electricity_CF = PowerUtility.set_CF
settings.__class__.get_electricity_CF = PowerUtility.get_CF
del define_units_of_measure, UnitsOfMeasure
