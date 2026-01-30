# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
#                          Yalin Li <mailto.yalin.li@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from __future__ import annotations
import pandas as pd
import numpy as np
import flexsolve as flx
from copy import copy as copy_
from numba import njit
from math import ceil
from warnings import warn
import biosteam as bst
from typing import Optional, Sequence, Collection, TYPE_CHECKING
from ._unit import Unit
from numpy.typing import NDArray
if TYPE_CHECKING: from ._system import System

__all__ = ('TEA',)


cashflow_columns = ('Depreciable capital [MM$]',
                    'Fixed capital investment [MM$]',
                    'Working capital [MM$]',
                    'Depreciation [MM$]',
                    'Loan [MM$]',
                    'Loan interest payment [MM$]',
                    'Loan payment [MM$]',
                    'Loan principal [MM$]',
                    'Annual operating cost (excluding depreciation) [MM$]',
                    'Sales [MM$]',
                    'Tax [MM$]',
                    'Incentives [MM$]',
                    'Taxed earnings [MM$]',
                    'Forwarded losses [MM$]',
                    'Net earnings [MM$]',
                    'Cash flow [MM$]',
                    'Discount factor',
                    'Net present value (NPV) [MM$]',
                    'Cumulative NPV [MM$]')

# %% Depreciation utilities

@njit(cache=True)
def generate_DDB_schedule(years):
    val = 1.
    arr = np.ones(years)
    factor = 2. / years
    for i in range(years):
        depreciation = val * factor
        arr[i] = depreciation
        val -= depreciation
    return arr
    
@njit(cache=True)
def generate_SYD_schedule(years):
    digit_sum = years * (years + 1.) * 0.5
    arr = np.ones(years)
    for i in range(years):
        arr[i] = (years - i) / digit_sum
    return arr

# %% Utilities for TEA calculations

@njit(cache=True)
def NPV_at_IRR(IRR, cashflow_array, duration_array):
    """Return NPV at given IRR and cashflow data."""
    return (cashflow_array/(1.+IRR)**duration_array).sum()

@njit(cache=True)
def loan_principal_with_interest(loan, interest):
    principal = 0
    k = 1. + interest
    for i in loan:
        principal *= k
        principal += i
    return principal

@njit(cache=True)
def solve_payment(loan_principal, interest, years):
    f = 1 + interest
    fn = f ** years
    return loan_principal * interest * fn / (fn - 1)

@njit(cache=True)
def taxable_earnings_with_fowarded_losses(taxable_cashflow): # Forwards losses to later years to reduce future taxes
    taxed_earnings = taxable_cashflow.copy()
    for i in range(taxed_earnings.size - 1):
        x = taxed_earnings[i]
        if x < 0:
            taxed_earnings[i] = 0
            taxed_earnings[i + 1] += x
    if taxed_earnings[-1] < 0: taxed_earnings[-1] = 0
    return taxed_earnings

@njit(cache=True)
def add_replacement_cost_to_cashflow_array(equipment_installed_cost, 
                                           equipment_lifetime, 
                                           cashflow_array, 
                                           venture_years,
                                           start):
    N_purchases = ceil(venture_years / equipment_lifetime) 
    for i in range(1, N_purchases):
        cashflow_array[start + i * equipment_lifetime] += equipment_installed_cost

def add_all_replacement_costs_to_cashflow_array(unit_capital_cost, cashflow_array, 
                                                venture_years, start,
                                                lang_factor):
    equipment_lifetime = unit_capital_cost.equipment_lifetime
    if equipment_lifetime:
        if lang_factor:
            installed_costs =  {i: j*lang_factor for i, j in unit_capital_cost.purchase_costs.items()}
        else:
            installed_costs = unit_capital_cost.installed_costs
        if isinstance(equipment_lifetime, int):
             add_replacement_cost_to_cashflow_array(sum(installed_costs.values()), 
                                                    equipment_lifetime,
                                                    cashflow_array,
                                                    venture_years,
                                                    start)
        elif isinstance(equipment_lifetime, dict):
            for name, installed_cost in installed_costs.items():
                lifetime = equipment_lifetime.get(name)
                if lifetime:
                    add_replacement_cost_to_cashflow_array(installed_cost, 
                                                           lifetime,
                                                           cashflow_array,
                                                           venture_years,
                                                           start)

@njit(cache=True)
def fill_taxable_and_nontaxable_cashflows_without_loans(
        D, C, S, C_FC, C_WC, FCI, WC, TDC, VOC, FOC, sales,
        startup_time,
        startup_VOCfrac,
        startup_FOCfrac,
        startup_salesfrac,
        construction_schedule,
        start):
    # Cash flow data and parameters
    # C_FC: Fixed capital
    # C_WC: Working capital
    # D: Depreciation
    # C: Annual operating cost (excluding depreciation)
    # S: Sales
    w0 = startup_time
    w1 = 1. - w0
    C[start] = (w0 * startup_VOCfrac * VOC + w1 * VOC
                + w0 * startup_FOCfrac * FOC + w1 * FOC)
    S[start] = w0 * startup_salesfrac * sales + w1 * sales
    start1 = start + 1
    C[start1:] = VOC + FOC
    S[start1:] = sales
    C_FC[:start] = FCI * construction_schedule
    C_WC[start-1] = WC
    C_WC[-1] = -WC

def taxable_and_nontaxable_cashflows(
        unit_capital_costs,
        D, C, S, C_FC, C_WC, Loan, LP,
        FCI, WC, TDC, VOC, FOC, sales,
        startup_time,
        startup_VOCfrac,
        startup_FOCfrac,
        startup_salesfrac,
        construction_schedule,
        finance_interest,
        finance_years,
        finance_fraction,
        start, years,
        lang_factor,
        accumulate_interest_during_construction,
    ):
    # Cash flow data and parameters
    # C_FC: Fixed capital
    # C_WC: Working capital
    # Loan: Money gained from loan
    # LP: Loan payment
    # D: Depreciation
    # C: Annual operating cost (excluding depreciation)
    # S: Sales
    fill_taxable_and_nontaxable_cashflows_without_loans(
        D, C, S, C_FC, C_WC, FCI, WC, TDC, VOC, FOC, sales,
        startup_time,
        startup_VOCfrac,
        startup_FOCfrac,
        startup_salesfrac,
        construction_schedule,
        start,
    )
    for i in unit_capital_costs:
        add_all_replacement_costs_to_cashflow_array(i, C_FC, years, start, lang_factor)
    if finance_interest:
        interest = finance_interest
        years = finance_years
        Loan[:start] = loan = finance_fraction*(C_FC[:start])
        if accumulate_interest_during_construction:
            loan_principal = loan_principal_with_interest(loan, interest)
        else:
            loan_principal = loan.sum()
        LP[start:start + years] = solve_payment(loan_principal, interest, years)
        taxable_cashflow = S - C - D - LP
        nontaxable_cashflow = D + Loan - C_FC - C_WC 
        if not accumulate_interest_during_construction: 
            nontaxable_cashflow[:start] -= loan * interest 
    else:
        taxable_cashflow = S - C - D
        nontaxable_cashflow = D - C_FC - C_WC
    return taxable_cashflow, nontaxable_cashflow

def NPV_with_sales(
        sales, 
        taxable_cashflow, 
        nontaxable_cashflow,
        depreciation,
        sales_coefficients,
        discount_factors,
        fill_tax_and_incentives,
    ):
    """Return NPV with an additional annualized sales."""
    taxable_cashflow = taxable_cashflow + sales * sales_coefficients
    tax = np.zeros_like(taxable_cashflow, dtype=float)
    incentives = tax.copy()
    fill_tax_and_incentives(
        incentives, 
        taxable_earnings_with_fowarded_losses(taxable_cashflow),
        nontaxable_cashflow, tax, depreciation
    )
    cashflow = nontaxable_cashflow + taxable_cashflow + incentives - tax
    return (cashflow/discount_factors).sum()

# %% Techno-Economic Analysis

_duration_array_cache = {}

class TEA:
    """
    Abstract TEA class for cash flow analysis.
    
    **Abstract methods**
    
    _DPI(installed_equipment_cost) -> float
        Should return the direct permanent investment (DPI) given the
        installed equipment cost.
    _TDC(DPI) -> float
        Should take direct permanent investment (DPI) as an argument
        and return total depreciable capital (TDC).
    _FCI(TDC) -> float
        Should take total depreciable capital (TDC) as an argument and return
        fixed capital investment (FCI).
    _FOC(FCI) -> float
        Should take fixed capital investment (FCI) as an arguments and return
        fixed operating cost without depreciation (FOC). 
    _fill_tax_and_incentives(incentives, taxable_cashflow, nontaxable_cashflow, tax)
        Should take two empty 1d arrays and fill them with incentive and tax cash flows.
        Additional parameters include taxable_cashflow (sales - costs - 
        depreciation - payments), and nontaxable_cashflow (depreciation - capital 
        cost - working capital).
    
    Parameters
    ----------
    system : 
        Should contain feed and product streams.
    IRR : 
        Internal rate of return (fraction).
    duration : 
        Start and end year of venture (e.g. (2018, 2038)).
    depreciation : 
        Depreciation schedule array or a string with format '{schedule}{years}', 
        where years is the number of years until the property value is zero 
        and schedule is one of the following: 'MACRS' (Modified Accelerated Cost Recovery System), 
        'SL' (straight line), 'DDB' (double declining balance), or 
        'SYD' (sum-of-the-years' digits). If years is not given, it defaults 
        to the venture years at run time. 
    operating_days :  
        Number of operating days per year.
    income_tax : 
        Combined federal and state income tax rate (fraction).
    lang_factor : 
        Lang factor for getting fixed capital investment
        from total purchase cost. If no lang factor, estimate
        capital investment using bare module factors.
    construction_schedule : 
        Construction investment fractions per year (e.g. (0.5, 0.5) for 50%
        capital investment in the first year and 50% investment in the second).
    startup_months : 
        Startup time in months.
    startup_FOCfrac : 
        Fraction of fixed operating costs incurred during startup.
    startup_VOCfrac : 
        Fraction of variable operating costs incurred during startup.
    startup_salesfrac : 
        Fraction of sales achieved during startup.
    WC_over_FCI : 
        Working capital as a fraction of fixed capital investment.
    finance_interest : 
        Yearly interest of capital cost financing as a fraction.
    finance_years : 
        Number of years the loan is paid for.
    finance_fraction :
        Fraction of capital cost that needs to be financed.
    
    Warning
    -------
    When using a Lang factor, a short-cut to get the FCI, we cannot work 
    backwards to get installed equipment cost. For practical purposes, the
    code assumes that installed equipment cost, total depreciable capital 
    (TDC), and fixed capital investment (FCI) are the same when the Lang
    factor is in use. In actuality, the installed equipment cost should be 
    less than the fixed capital investment. 
    
    Examples
    --------
    :doc:`../tutorial/Techno-economic_analysis` 

    """
    __slots__ = ('system', 'income_tax', 'WC_over_FCI',
                 'finance_interest', 'finance_years', 'finance_fraction',
                 '_construction_schedule', '_startup_time',
                 'startup_FOCfrac', 'startup_VOCfrac', 'startup_salesfrac',
                 '_startup_schedule', '_operating_days',
                 '_duration', '_depreciation_key', '_depreciation',
                 '_years', '_duration', '_start',  'IRR', '_IRR', '_sales',
                 '_duration_array_cache', 'accumulate_interest_during_construction')
    
    #: Available depreciation schedules. Defaults include modified 
    #: accelerated cost recovery system from U.S. IRS publication 946 (MACRS),
    #: half-year convention.
    depreciation_schedules: dict[tuple(str, int), NDArray[float]] = {
        ('MACRS', 3): np.array([.3333, .4445, .1481, .0741]),

        ('MACRS', 5): np.array([.2000, .3200, .1920,
                                .1152, .1152, .0576]),
              
        ('MACRS', 7):  np.array([.1429, .2449, .1749,
                                 .1249, .0893, .0892,
                                 .0893, .0446]),
        
        ('MACRS', 10): np.array([.1000, .1800, .1440,
                                 .1152, .0922, .0737,
                                 .0655, .0655, .0656,
                                 .0655, .0328]),
      
        ('MACRS', 15): np.array([.0500, .0950, .0855,
                                 .0770, .0693, .0623,
                                 .0590, .0590, .0591,
                                 .0590, .0591, .0590,
                                 .0591, .0590, .0591,
                                 .0295]),
      
        ('MACRS', 20): np.array([0.03750, 0.07219, 0.06677,
                                 0.06177, 0.05713, 0.05285,
                                 0.04888, 0.04522, 0.04462,
                                 0.04461, 0.04462, 0.04461,
                                 0.04462, 0.04461, 0.04462,
                                 0.04461, 0.04462, 0.04461,
                                 0.04462, 0.04461, 0.02231])
    }
    
    #: Investment site factors used to multiply the total permanent 
    #: investment (TPI), also known as total fixed capital (FCI), to 
    #: account for locality cost differences based on labor availability,
    #: workforce efficiency, local rules, etc.
    investment_site_factors: dict[str, float] = {
        'U.S. Gulf Coast': 1.0,
        'U.S. Southwest': 0.95,
        'U.S. Northwest': 1.10,
        'U.S. Midwest': 1.15,
        'U.S. West Coast': 1.25,
        'Western Europe': 1.2,
        'Mexico': 0.95,
        'Japan': 1.15,
        'Pacific Rim': 1.0,
        'India': 0.85,
    }

    class Accounting:
        __slots__ = ('index', 'data', 'names', 'units')
        
        def __init__(self, units, names=None):
            self.units = units
            self.names = names
            self.index = []
            self.data =[]
        
        def entry(self, index: str, cost: list|float, notes: str = ''):
            self.index.append(index)
            if getattr(cost, 'ndim', 0) == 0: cost = float(cost)
            if isinstance(cost, (float, int, str)):
                self.data.append([notes, cost])
            else:
                self.data.append([notes, *cost])

        @property
        def total_costs(self):
            if self.names is None:
                return sum([i[1] for i in self.data])
            else:
                N = len(self.data[0])
                return np.array([sum([i[index] for i in self.data]) for index in range(1, N)])
        
        def table(self):
            names = self.names
            units = self.units
            index = self.index
            return pd.DataFrame(
                self.data, 
                index=pd.MultiIndex.from_tuples(index) if isinstance(index[0], tuple) else index,
                columns=('Notes', f'Cost [{units}]') if names is None else ('Notes', *[f'{i}\n[{units}]' for i in names]),
            )

    def __init_subclass__(cls, isabstract=False):
        if isabstract: return
        for method in ('_DPI', '_TDC', '_FCI', '_FOC'):
            if not hasattr(cls, method):
                breakpoint()
                raise NotImplementedError(
                    f"subclass must implement a '{method}' method unless the "
                     "'isabstract' keyword argument is True"
                )

    def copy(self, system=None):
        """Create a copy."""
        new = copy_(self)
        if system is not None:
            new.system = system
            system._TEA = new
        return new

    def __init__(self, system: System, IRR: float, duration: tuple[int, int], 
                 depreciation: str|NDArray[float], income_tax: float,
                operating_days: float, lang_factor: float|None, 
                construction_schedule: Sequence[float],
                startup_months: float, startup_FOCfrac: float, startup_VOCfrac: float,
                startup_salesfrac: float, WC_over_FCI: float,  finance_interest: float,
                finance_years: int, finance_fraction: float,
                accumulate_interest_during_construction: bool=False):
        #: System being evaluated.
        self.system: System = system
        
        self.duration = duration
        self.depreciation = depreciation
        self.construction_schedule = construction_schedule
        self.startup_months = startup_months
        self.operating_days = operating_days
        
        #: Internal rate of return (fraction).
        self.IRR: float = IRR
        
        #: Combined federal and state income tax rate (fraction).
        self.income_tax: float = income_tax
        
        self.lang_factor = lang_factor
        
        #: Fraction of fixed operating costs incurred during startup.
        self.startup_FOCfrac: float = startup_FOCfrac
        
        #: Fraction of variable operating costs incurred during startup.
        self.startup_VOCfrac: float = startup_VOCfrac
        
        #: Fraction of sales achieved during startup.
        self.startup_salesfrac: float = startup_salesfrac
        
        #: Working capital as a fraction of fixed capital investment.
        self.WC_over_FCI: float = WC_over_FCI
        
        #: Yearly interest of capital cost financing as a fraction.
        self.finance_interest: float = finance_interest
        
        #: Number of years the loan is paid for.
        self.finance_years: int = finance_years
        
        #: Fraction of capital cost that needs to be financed.
        self.finance_fraction: float = finance_fraction
        
        #: Guess IRR for solve_IRR method
        self._IRR: float = IRR
        
        #: Guess cost for solve_price method
        self._sales: float = 0
        
        #: Whether to immediately pay interest before operation or to accumulate interest during construction
        self.accumulate_interest_during_construction = accumulate_interest_during_construction
        
        #: For convenience, set a TEA attribute for the system
        system._TEA = self

    def _get_duration(self):
        return (self._start, self._years)

    def _DPI(self, installed_equipment_cost):
        return installed_equipment_cost # For compatibility with Lang factors

    def _TDC(self, DPI):
        return DPI # For compatibility with Lang factors
    
    def _FCI(self, TDC):
        return TDC # For compatibility with Lang factors

    @property
    def save_report(self):
        return self.system.save_report

    @property
    def units(self) -> set[Unit]:
        """All unit operations with costs."""
        return self.system.cost_units  

    @property
    def feeds(self) -> list[Unit]:
        """All feed streams."""
        return self.system.feeds  
      
    @property
    def products(self) -> list[Unit]:
        """All product streams."""
        return self.system.products

    @property
    def operating_days(self) -> float:
        """Number of operating days per year."""
        return self.system.operating_hours / 24
    @operating_days.setter
    def operating_days(self, days):
        """Number of operating days per year."""
        self.operating_hours = 24 * days
    
    @property
    def operating_hours(self) -> float:
        """Number of operating hours per year."""
        return self.system.operating_hours
    @operating_hours.setter
    def operating_hours(self, hours):
        self.system.operating_hours = hours
    
    @property
    def lang_factor(self) -> float|None:
        """
        Lang factor for getting fixed capital investment from 
        total purchase cost. If no lang factor, estimate capital investment 
        using bare module factors.
        
        """
        return self.system.lang_factor
    @lang_factor.setter
    def lang_factor(self, lang_factor):
        self.system.lang_factor = lang_factor
    
    @property
    def duration(self) -> tuple[int, int]:
        """Start and end year of venture."""
        return self._duration
    @duration.setter
    def duration(self, duration):
        start, end = [int(i) for i in duration]
        self._duration = (start, end)
        self._years = end - start

    @property
    def depreciation(self) -> str|NDArray[float]:
        """
        Depreciation schedule array or a string with format '{schedule}{years}', 
        where years is the number of years until the property value is zero 
        and schedule is one of the following: 'MACRS' (Modified Accelerated Cost Recovery System), 
        'SL' (straight line), 'DDB' (double declining balance), or 
        'SYD' (sum-of-the-years' digits). If years is not given, it defaults 
        to the venture years at run time.
        """
        return self._depreciation
    @depreciation.setter
    def depreciation(self, depreciation):
        if isinstance(depreciation, str):
            self._depreciation_key = self._depreciation_key_from_name(depreciation)
            self._depreciation = depreciation
        else:
            try:
                self._depreciation = np.array(depreciation, dtype=float)
            except:
                raise TypeError(
                   f"invalid depreciation type '{type(depreciation).__name__}'; "
                    "depreciation must be either an array or a string"
                ) from None
            else:
                self._depreciation_key = None
    
    @classmethod
    def _depreciation_key_from_name(cls, name):
        for prefix in ('MACRS', 'SL', 'DDB', 'SYD'):
            if name.startswith(prefix):
                years = name[len(prefix):]
                key = (prefix, int(years) if years else None)
                if prefix == 'MACRS' and key not in cls.depreciation_schedules:
                    raise ValueError(
                        f"depreciation name {repr(name)} has a valid "
                         "format, but is not yet implemented in BioSTEAM"
                    )
                return key
        raise ValueError(
               f"invalid depreciation name {repr(name)}; "
                "name must have format '{schedule}{years}', "
                "where years is the number of years until the property value is zero "
                "and schedule is one of the following: 'MACRS' (Modified Accelerated Cost Recovery System), "
                "'SL' (straight line), 'DDB' (double declining balance), or "
                "'SYD' (sum-of-the-years' digits)"
            )

    @classmethod
    def _depreciation_array_from_key(cls, key):
        depreciation_schedules = cls.depreciation_schedules
        if key in depreciation_schedules:
            return depreciation_schedules[key]
        else:
            schedule, years = key
            if schedule == 'SL':
                arr = np.full(years, 1./years)
            elif schedule == 'DDB':
                arr = generate_DDB_schedule(years)
            elif schedule == 'SYD':
                arr = generate_SYD_schedule(years)
            else: # pragma: no cover
                raise RuntimeError(f'unknown depreciation schedule {repr(schedule)}')
            depreciation_schedules[key] = arr
            return arr
    
    @property
    def construction_schedule(self) -> Sequence[float]:
        """Construction investment fractions per year, starting from year 0.
        For example, for 50% capital investment in year 0 and 50% investment 
        in year 1, use (0.5, 0.5)."""
        return self._construction_schedule
    @construction_schedule.setter
    def construction_schedule(self, schedule):
        self._construction_schedule = np.array(schedule, dtype=float)
        self._start = len(schedule)
    
    @property
    def startup_months(self) -> float:
        return self._startup_time * 12.
    @startup_months.setter
    def startup_months(self, months):
        assert months <= 12., "startup time must be less than a year"
        self._startup_time = months/12.
    
    @property
    def sales(self) -> float:
        """Total sales [USD/yr]."""
        return self.system.sales
    @property
    def material_cost(self) -> float:
        """Total material cost [USD/yr]."""
        return self.system.material_cost
    @property
    def utility_cost(self) -> float:
        """Total utility cost [USD/yr]."""
        return self.system.utility_cost
    @property
    def purchase_cost(self):
        """Total purchase cost [USD]."""
        return self.system.purchase_cost
    @property
    def installed_equipment_cost(self) -> float:
        """Total installed cost [USD]."""
        return self.system.installed_equipment_cost
    @property
    def DPI(self) -> float:
        """Direct permanent investment [USD]."""
        return self._DPI(self.installed_equipment_cost)
    @property
    def TDC(self) -> float:
        """Total depreciable capital [USD]."""
        return self._TDC(self.DPI)
    @property
    def FCI(self) -> float:
        """Fixed capital investment [USD]."""
        return self._FCI(self.TDC)
    @property
    def TCI(self) -> float:
        """Total capital investment [USD]."""
        return (1. + self.WC_over_FCI)*self.FCI
    @property
    def FOC(self) -> float:
        """Fixed operating costs [USD/yr]."""
        return self._FOC(self.FCI)
    @property
    def VOC(self) -> float:
        """Variable operating costs [USD/yr]."""
        return self.material_cost + self.utility_cost
    @property
    def AOC(self) -> float:
        """Annual operating cost excluding depreciation [USD/yr]."""
        return self.FOC + self.VOC
    @property
    def working_capital(self) -> float:
        return self.WC_over_FCI * self.FCI
    
    @property
    def annual_depreciation(self) -> float:
        """Depreciation [USD/yr] equivalent to TDC divided by the duration of the venture."""
        return self.TDC/(self.duration[1]-self.duration[0])

    @property
    def ROI(self) -> float:
        """Return on investment [1/yr] without accounting for annualized depreciation."""
        return self.net_earnings / self.TCI
    @property
    def net_earnings(self) -> float:
        """Net earnings without accounting for annualized depreciation."""
        net_earnings = self.sales - self.AOC
        if net_earnings < 0:
            return net_earnings
        else:
            return (1 - self.income_tax) * net_earnings
    @property
    def PBP(self) -> float:
        """Pay back period [yr] without accounting for annualized depreciation."""
        FCI = self.FCI
        net_earnings = self.net_earnings
        return FCI/net_earnings

    def _get_duration_array(self):
        key = start, years = (self._start, self._years)
        if key in _duration_array_cache:
            duration_array = _duration_array_cache[key]
        else:
            if len(_duration_array_cache) > 100: _duration_array_cache.clear()
            _duration_array_cache[key] = duration_array = np.arange(-start+1, years+1, dtype=float)
        return duration_array

    def _get_depreciation_array(self):
        key = self._depreciation_key
        if key is None: 
            return self._depreciation
        else:
            schedule, years = self._depreciation_key
            if years is None:
                years = self._years
                key = (schedule, years)
            return self._depreciation_array_from_key(key)
            
    def _fill_depreciation_array(self, D, start, years, TDC):
        depreciation_array = self._get_depreciation_array()
        N_depreciation_years = depreciation_array.size
        if N_depreciation_years > years:
            raise RuntimeError('depreciation schedule is longer than plant lifetime')
        D[start:start + N_depreciation_years] = TDC * depreciation_array

    def get_cashflow_table(self):
        """Return DataFrame of the cash flow analysis."""
        # Cash flow data and parameters
        # index: Year since construction until end of venture
        # C_D: Depreciable capital
        # C_FC: Fixed capital
        # C_WC: Working capital
        # D: Depreciation
        # L: Loan revenue
        # LI: Loan interest payment
        # LP: Loan payment
        # LPl: Loan principal
        # C: Annual operating cost (excluding depreciation)
        # S: Sales
        # T: Tax
        # I: Incentives
        # TE: Taxed earnings
        # FL: Forwarded losses
        # NE: Net earnings
        # CF: Cash flow
        # DF: Discount factor
        # NPV: Net present value
        # CNPV: Cumulative NPV
        TDC = self.TDC
        FCI = self._FCI(TDC)
        start = self._start
        years = self._years
        FOC = self._FOC(FCI)
        VOC = self.VOC
        sales = self.sales
        length = start + years
        C_D, C_FC, C_WC, D, L, LI, LP, LPl, C, S, T, I, TE, FL, NE, CF, DF, NPV, CNPV = data = np.zeros((19, length))
        self._fill_depreciation_array(D, start, years, TDC)
        w0 = self._startup_time
        w1 = 1. - w0
        C[start] = (w0*self.startup_VOCfrac*VOC + w1*VOC
                    + w0*self.startup_FOCfrac*FOC + w1*FOC)
        S[start] = w0*self.startup_salesfrac*sales + w1*sales
        start1 = start + 1
        C[start1:] = VOC + FOC
        S[start1:] = sales
        WC = self.WC_over_FCI * FCI
        C_D[:start] = TDC*self._construction_schedule
        C_FC[:start] = FCI*self._construction_schedule
        C_WC[start-1] = WC
        C_WC[-1] = -WC
        system = self.system
        lang_factor = system.lang_factor
        unit_capital_costs = system.unit_capital_costs.values() if isinstance(system, bst.AgileSystem) else system.cost_units
        for i in unit_capital_costs: add_all_replacement_costs_to_cashflow_array(i, C_FC, years, start, lang_factor)
        if self.finance_interest:
            interest = self.finance_interest
            years = self.finance_years
            end = start + years
            L[:start] = loan = self.finance_fraction*(C_FC[:start])
            accumulate_interest_during_construction = self.accumulate_interest_during_construction
            if accumulate_interest_during_construction:
                initial_loan_principal = loan_principal_with_interest(loan, interest)
            else:
                initial_loan_principal = loan.sum()
            LP[start:end] = solve_payment(initial_loan_principal, interest, years)
            loan_principal = 0
            if accumulate_interest_during_construction:
                for i in range(end):
                    LI[i] = li = (loan_principal + L[i]) * interest 
                    LPl[i] = loan_principal = loan_principal - LP[i] + li + L[i]
            else:
                for i in range(end):
                    if i < start: 
                        li = 0
                    else:
                        li = (loan_principal + L[i]) * interest 
                    LI[i] = li
                    LPl[i] = loan_principal = loan_principal - LP[i] + li + L[i]
                LI[:start] = L[:start] * interest # Interest still needs to be payed
            taxable_cashflow = S - C - D - LP
            nontaxable_cashflow = D + L - C_FC - C_WC
            if not accumulate_interest_during_construction:
                nontaxable_cashflow[:start] -= LI[:start] # Subtract the interest during construction from NPV
        else:
            taxable_cashflow = S - C - D
            nontaxable_cashflow = D - C_FC - C_WC
        TE[:] = taxable_earnings_with_fowarded_losses(taxable_cashflow)
        FL[1:] = (taxable_cashflow - TE).cumsum()[:-1]
        self._fill_tax_and_incentives(
            I, TE, nontaxable_cashflow, T, D
        )
        NE[:] = taxable_cashflow + I - T
        CF[:] = NE + nontaxable_cashflow
        DF[:] = 1/(1.+self.IRR)**self._get_duration_array()
        NPV[:] = CF * DF
        CNPV[:] = NPV.cumsum()
        DF *= 1e6
        data /= 1e6
        return pd.DataFrame(data.transpose(),
                            index=np.arange(self._duration[0]-start, self._duration[1]),
                            columns=cashflow_columns)
    @property
    def NPV(self) -> float:
        """Net present value."""
        taxable_cashflow, nontaxable_cashflow, depreciation = self._taxable_nontaxable_depreciation_cashflows()
        tax = np.zeros_like(taxable_cashflow)
        incentives = tax.copy()
        self._fill_tax_and_incentives(
            incentives,
            taxable_earnings_with_fowarded_losses(taxable_cashflow),
            nontaxable_cashflow, tax, depreciation
        )
        cashflow = nontaxable_cashflow + taxable_cashflow + incentives - tax
        return NPV_at_IRR(self.IRR, cashflow, self._get_duration_array())
    
    def _AOC(self, FCI):
        """Return AOC at given FCI"""
        return self._FOC(FCI) + self.VOC
    
    def _taxable_nontaxable_depreciation_cashflows(self):
        """Return taxable, nontaxable and depreciation cash flows by year as a tuple[1d array, 1d array, 1d array]."""
        # Cash flow data and parameters
        # C_FC: Fixed capital
        # C_WC: Working capital
        # Loan: Money gained from loan
        # LP: Loan payment
        # D: Depreciation
        # C: Annual operating cost (excluding depreciation)
        # S: Sales
        # NE: Net earnings
        # CF: Cash flow
        TDC = self.TDC
        FCI = self._FCI(TDC)
        start = self._start
        years = self._years
        FOC = self._FOC(FCI)
        VOC = self.VOC
        D, C_FC, C_WC, Loan, LP, C, S = np.zeros((7, start + years))
        self._fill_depreciation_array(D, start, years, TDC)
        WC = self.WC_over_FCI * FCI
        system = self.system
        return (
            *taxable_and_nontaxable_cashflows(
                system.unit_capital_costs if isinstance(system, bst.AgileSystem) else system.cost_units,
                D, C, S, C_FC, C_WC, Loan, LP,
                FCI, WC, TDC, VOC, FOC, self.sales,
                self._startup_time,
                self.startup_VOCfrac,
                self.startup_FOCfrac,
                self.startup_salesfrac,
                self._construction_schedule,
                self.finance_interest,
                self.finance_years,
                self.finance_fraction,
                start, years,
                self.lang_factor,
                self.accumulate_interest_during_construction,
            ),
            D
        )
    
    def _fill_tax_and_incentives(self, incentives, taxable_cashflow, nontaxable_cashflow, tax, depreciation):
        tax[:] = self.income_tax * taxable_cashflow
    
    def _net_earnings_and_nontaxable_cashflow_arrays(self):
        taxable_cashflow, nontaxable_cashflow, depreciation = self._taxable_nontaxable_depreciation_cashflows()
        size = taxable_cashflow.size
        tax = np.zeros(size)
        incentives = tax.copy()
        self._fill_tax_and_incentives(
            incentives, 
            taxable_earnings_with_fowarded_losses(taxable_cashflow), 
            nontaxable_cashflow, tax, depreciation
        )
        net_earnings = taxable_cashflow + incentives - tax
        return net_earnings, nontaxable_cashflow
    
    @property
    def cashflow_array(self) -> NDArray[float]:
        """Cash flows by year."""
        return sum(self._net_earnings_and_nontaxable_cashflow_arrays())
    
    @property
    def net_earnings_array(self) -> NDArray[float]:
        """Net earnings by year."""
        return self._net_earnings_and_nontaxable_cashflow_arrays()[0]
    
    def production_costs(self, products: Sequence[bst.Stream], with_annual_depreciation: Optional[bool]=True):
        """
        Return production cost of products [USD/yr].
        
        Parameters
        ----------
        products : 
            Main products of the system.
        with_annual_depreciation: 
            Whether to add annualized depreciation to the production costs.
        
        Notes
        -----
        If there is more than one main product, The production cost is
        proportionally allocated to each of the main products with respect to
        their marketing values. The marketing value of each product is
        determined by the annual production multiplied by its selling price.
        
        """
        system = self.system
        market_values = np.array([system.get_market_value(i) for i in products])
        total_market_value = market_values.sum()
        weights = market_values/total_market_value
        return weights * self.total_production_cost(products, with_annual_depreciation)
        
    def total_production_cost(self, products: Collection[bst.Stream], with_annual_depreciation: Optional[bool]=True):
        """
        Return total production cost of products [USD/yr].
        
        Parameters
        ----------
        products : 
            Main products of the system.
        with_annual_depreciation: 
            Whether to add annualized depreciation to the production costs.
        
        """
        system = self.system
        coproduct_sales = self.sales - sum([system.get_market_value(i) for i in products])
        if with_annual_depreciation:
            TDC = self.TDC
            annual_depreciation = TDC/(self.duration[1]-self.duration[0])
            AOC = self._AOC(self._FCI(TDC))
            return AOC - coproduct_sales + annual_depreciation
        else:
            return self.AOC - coproduct_sales
    
    def solve_IRR(self, financing=True, bounds=None):
        """Return the IRR at the break even point (NPV = 0) through cash flow analysis."""
        IRR = self._IRR
        if not IRR or np.isnan(IRR) or IRR < 0.: IRR = 0.01
        if financing:
            args = (self.cashflow_array, self._get_duration_array())
            if bounds:
                IRR = flx.IQ_interpolation(
                    NPV_at_IRR, *bounds, x=IRR, xtol=1e-6, ytol=10.,
                    maxiter=200, args=args, checkiter=False
                )
            else:
                IRR = flx.aitken_secant(
                    NPV_at_IRR, IRR, 1.0001 * IRR + 1e-3, xtol=1e-6, ytol=10.,
                    maxiter=200, args=args, checkiter=False
                )
        else:
            financing_values = self.finance_fraction, self.finance_interest
            self.finance_fraction = self.finance_interest = None
            try:
                args = (self.cashflow_array, self._get_duration_array())
                if bounds:
                    IRR = flx.IQ_interpolation(
                        NPV_at_IRR, *bounds, x=IRR, xtol=1e-6, ytol=10.,
                        maxiter=200, args=args, checkiter=False
                    )
                else:
                    IRR = flx.aitken_secant(
                        NPV_at_IRR, IRR, 1.0001 * IRR + 1e-3, xtol=1e-6, ytol=10.,
                        maxiter=200, args=args, checkiter=False
                    )
            finally:
                self.finance_fraction, self.finance_interest = financing_values
        self._IRR = IRR
        return IRR
        
    def solve_price(self, streams: bst.Stream|Collection[bst.Stream]):
        """
        Return the price [USD/kg] of a stream(s) at the break even point (NPV = 0)
        through cash flow analysis. 
        
        Parameters
        ----------
        streams :
            Streams with variable selling price.
            
        """
        if isinstance(streams, bst.Stream): streams = [streams]
        system = self.system
        price2cost = sum([system._price2cost(i) for i in streams])
        if price2cost == 0.: raise ValueError('cannot solve price of empty streams')
        try:
            sales = self.solve_sales()
        except:
            original_prices = [i.price for i in streams]
            for i in streams: i.price = 0.
            sales = self.solve_sales()
            current_price = 0.
            for i, j in zip(streams, original_prices): i.price = j 
        else:
            current_price = sum([system.get_market_value(i) for i in streams]) / abs(price2cost)
        return current_price + sales / price2cost 
        
    def VOC_table(
            self, products, functional_unit='MT', with_products=False, 
        ):
        return bst.report.voc_table(
            [self.system], products, 
            system_names=None, functional_unit=functional_unit,
            with_products=with_products, dataframe=True
        )
        
    def CAPEX_table(self):
        return NotImplemented
    
    def FOC_table(self):
        return NotImplemented
    
    def solve_sales(self):
        """
        Return the required additional sales [USD] to reach the breakeven 
        point (NPV = 0) through cash flow analysis. 
        
        """
        discount_factors = (1 + self.IRR)**self._get_duration_array()
        sales_coefficients = np.ones_like(discount_factors, dtype=float)
        start = self._start
        sales_coefficients[:start] = 0
        w0 = self._startup_time
        sales_coefficients[start] =  w0*self.startup_salesfrac + (1.-w0)
        taxable_cashflow, nontaxable_cashflow, depreciation = self._taxable_nontaxable_depreciation_cashflows()
        if np.isnan(taxable_cashflow).any():
            warn('nan encountered in cashflow array; resimulating system', category=RuntimeWarning)
            self.system.simulate()
            taxable_cashflow, nontaxable_cashflow, depreciation = self._taxable_nontaxable_depreciation_cashflows()
            if np.isnan(taxable_cashflow).any():
                raise RuntimeError('nan encountered in cashflow array')
        args = (taxable_cashflow, 
                nontaxable_cashflow, 
                depreciation,
                sales_coefficients,
                discount_factors,
                self._fill_tax_and_incentives)
        x0 = self._sales if np.isfinite(self._sales) else 0
        f = NPV_with_sales
        y0 = f(x0, *args)
        x1 = x0 - y0 / self._years # First estimate
        try:
            sales = flx.aitken_secant(f, x0, x1, xtol=10, ytol=100.,
                                      maxiter=1000, args=args, checkiter=True)
        except:
            bracket = flx.find_bracket(f, x0, x1, args=args)
            sales = flx.IQ_interpolation(f, *bracket, args=args, xtol=10, ytol=100, maxiter=1000, checkiter=False)
        self._sales = sales
        return sales
    
    def __repr__(self):
        return f'{type(self).__name__}({self.system.ID}, ...)'
    
    def _info(self):
        return (f'{type(self).__name__}: {self.system}\n'
                f'NPV: {self.NPV:,.0f} USD at {self.IRR:.1%} IRR')
    
    def show(self):
        """Prints information on unit."""
        print(self._info())
    _ipython_display_ = show