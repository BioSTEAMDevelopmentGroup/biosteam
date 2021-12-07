# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
#                          Yalin Li <zoe.yalin.li@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import pandas as pd
import numpy as np
import flexsolve as flx
from copy import copy as copy_
from numba import njit
from math import ceil
from warnings import warn
import biosteam as bst

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
def initial_loan_principal(loan, interest):
    principal = 0
    k = 1. + interest
    for i in loan:
        principal += i
        principal *= k
    return principal

@njit(cache=True)
def final_loan_principal(payment, principal, interest, years):
    for iter in range(years):
        principal += principal * interest - payment
    return principal

def solve_payment(payment, loan, interest, years):
    principal = initial_loan_principal(loan, interest)
    payment = flx.aitken_secant(final_loan_principal,
                                payment, payment+10., 1., 10.,
                                args=(principal, interest, years),
                                maxiter=200, checkiter=False)
    return payment

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
        start
    )
    for i in unit_capital_costs:
        add_all_replacement_costs_to_cashflow_array(i, C_FC, years, start, lang_factor)
    if finance_interest:
        interest = finance_interest
        years = finance_years
        Loan[:start] = loan = finance_fraction*(C_FC[:start]+C_WC[:start])
        LP[start:start + years] = solve_payment(loan.sum()/years * (1. + interest),
                                                loan, interest, years)
        taxable_cashflow = S - C - D - LP
        nontaxable_cashflow = D + Loan - C_FC - C_WC 
    else:
        taxable_cashflow = S - C - D
        nontaxable_cashflow = D - C_FC - C_WC
    return taxable_cashflow, nontaxable_cashflow

def NPV_with_sales(
        sales, 
        taxable_cashflow, 
        nontaxable_cashflow,
        sales_coefficients,
        discount_factors,
        fill_tax_and_incentives,
    ):
    """Return NPV with an additional annualized sales."""
    taxable_cashflow = taxable_cashflow + sales * sales_coefficients
    tax = np.zeros_like(taxable_cashflow)
    incentives = tax.copy()
    fill_tax_and_incentives(incentives, taxable_cashflow, nontaxable_cashflow, tax)
    cashflow = nontaxable_cashflow + taxable_cashflow + incentives - tax
    return (cashflow/discount_factors).sum()

# %% Techno-Economic Analysis

_duration_array_cache = {}

class TEA:
    """
    Abstract TEA class for cash flow analysis.
    
    **Abstract methods**
    
    _DPI(installed_equipment_cost) -> DPI
        Should return the direct permanent investment given the
        installed equipment cost.
    _TDC(DPI) -> TDC
        Should take direct permanent investment as an argument
        and return total depreciable capital.
    _FCI(TDC) -> FCI
        Should take total depreciable capital as an argument and return
        fixed capital investment.
    _FOC(FCI) -> FOC
        Should take fixed capital investment as an arguments and return
        fixed operating cost without depreciation. 
    _fill_tax_and_incentives(incentives, taxable_cashflow, nontaxable_cashflow, tax)
        Should take two empty 1d arrays and fill them with incentive and tax cash flows.
        Additional parameters include taxable_cashflow (sales - costs - 
        depreciation - payments), and nontaxable_cashflow (depreciation - capital 
        cost - working capital).
    
    Parameters
    ----------
    system : System
        Should contain feed and product streams.
    IRR : float
        Internal rate of return (fraction).
    duration : tuple[int, int]
        Start and end year of venture (e.g. (2018, 2038)).
    depreciation : array or str
        Depreciation schedule array or a string with format '{schedule}{years}', 
        where years is the number of years until the property value is zero 
        and schedule is one of the following: 'MACRS' (Modified Accelerated Cost Recovery System), 
        'SL' (straight line), 'DDB' (double declining balance), or 
        'SYD' (sum-of-the-years' digits). If years is not given, it defaults 
        to the venture years at run time. 
    operating_days : float 
        Number of operating days per year.
    income_tax : float
        Combined federal and state income tax rate (fraction).
    lang_factor : float
        Lang factor for getting fixed capital investment
        from total purchase cost. If no lang factor, estimate
        capital investment using bare module factors.
    construction_schedule : 1d array [float]
        Construction investment fractions per year (e.g. (0.5, 0.5) for 50%
        capital investment in the first year and 50% investment in the second).
    startup_months : float
        Startup time in months.
    startup_FOCfrac : float
        Fraction of fixed operating costs incurred during startup.
    startup_VOCfrac : float
        Fraction of variable operating costs incurred during startup.
    startup_salesfrac : float
        Fraction of sales achieved during startup.
    WC_over_FCI : float
        Working capital as a fraction of fixed capital investment.
    finanace_interest : float
        Yearly interest of capital cost financing as a fraction.
    finance_years : int
                    Number of years the loan is paid for.
    finance_fraction : float
                       Fraction of capital cost that needs to be financed.
    
    Warning
    -------
    When using a Lang factor, a short-cut to get the FCI, we cannot work 
    backwards to get installed equipment cost. For practical purposes, the
    code assumes that installed equipment cost, total depreciable capital 
    (TDC), and fixed capital investment (FCI) are the same when the Lang
    factor is in use. In actuallity, the installed equipment cost should be 
    less than the fixed capital investment. 
    
    Examples
    --------
    :doc:`tutorial/Techno-economic_analysis` 

    """
    __slots__ = ('system', 'income_tax', 'WC_over_FCI',
                 'finance_interest', 'finance_years', 'finance_fraction',
                 '_construction_schedule', '_startup_time',
                 'startup_FOCfrac', 'startup_VOCfrac', 'startup_salesfrac',
                 '_startup_schedule', '_operating_days',
                 '_duration', '_depreciation_key', '_depreciation',
                 '_years', '_duration', '_start',  'IRR', '_IRR', '_sales',
                 '_duration_array_cache')
    
    # Defaults include modified accelerated cost recovery system from
    # U.S. IRS publicaiton 946 (MACRS), half-year convention
    #: dict[tuple(str, int), 1d-array] Available depreciation schedules.
    depreciation_schedules = {
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
    
    #: dict[str, float] Investment site factors used to multiply the total permanent 
    #: investment (TPI), also known as total fixed capital (FCI), to 
    #: account for locality cost differences based on labor availability,
    #: workforce efficiency, local rules, etc.
    investment_site_factors = {
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

    def __init__(self, system, IRR, duration, depreciation, income_tax,
                operating_days, lang_factor, construction_schedule,
                startup_months, startup_FOCfrac, startup_VOCfrac,
                startup_salesfrac, WC_over_FCI,  finance_interest,
                finance_years, finance_fraction):
        #: [System] System being evaluated.
        self.system = system
        
        self.duration = duration
        self.depreciation = depreciation
        self.construction_schedule = construction_schedule
        self.startup_months = startup_months
        self.operating_days = operating_days
        
        #: [float]  Internal rate of return (fraction).
        self.IRR = IRR
        
        #: [float] Combined federal and state income tax rate (fraction).
        self.income_tax = income_tax
        
        self.lang_factor = lang_factor
        
        #: [float] Fraction of fixed operating costs incurred during startup.
        self.startup_FOCfrac = startup_FOCfrac
        
        #: [float] Fraction of variable operating costs incurred during startup.
        self.startup_VOCfrac = startup_VOCfrac
        
        #: [float] Fraction of sales achieved during startup.
        self.startup_salesfrac = startup_salesfrac
        
        #: [float] Working capital as a fraction of fixed capital investment.
        self.WC_over_FCI = WC_over_FCI
        
        #: [float] Yearly interest of capital cost financing as a fraction.
        self.finance_interest = finance_interest
        
        #: [int] Number of years the loan is paid for.
        self.finance_years = finance_years
        
        #: [float] Fraction of capital cost that needs to be financed.
        self.finance_fraction = finance_fraction
        
        #: Guess IRR for solve_IRR method
        self._IRR = IRR
        
        #: Guess cost for solve_price method
        self._sales = 0
        
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
    def units(self):
        """set[Unit] All unit operations with costs."""
        return self.system.cost_units  

    @property
    def feeds(self):
        """list[Unit] All feed streams."""
        return self.system.feeds  
      
    @property
    def products(self):
        """list[Unit] All product streams."""
        return self.system.products

    @property
    def operating_days(self):
        """[float] Number of operating days per year."""
        return self.system.operating_hours / 24
    @operating_days.setter
    def operating_days(self, days):
        """[float] Number of operating days per year."""
        self.operating_hours = 24 * days
    
    @property
    def operating_hours(self):
        """[float] Number of operating hours per year."""
        return self.system.operating_hours
    @operating_hours.setter
    def operating_hours(self, hours):
        self.system.operating_hours = hours
    
    @property
    def lang_factor(self):
        """
        [float] Lang factor for getting fixed capital investment from 
        total purchase cost. If no lang factor, estimate capital investment 
        using bare module factors.
        
        """
        return self.system.lang_factor
    @lang_factor.setter
    def lang_factor(self, lang_factor):
        self.system.lang_factor = lang_factor
    
    @property
    def duration(self):
        """tuple[int, int] Start and end year of venture."""
        return self._duration
    @duration.setter
    def duration(self, duration):
        start, end = [int(i) for i in duration]
        self._duration = (start, end)
        self._years = end - start

    @property
    def depreciation(self):
        """
        [array or str] Depreciation schedule array or a string with format '{schedule}{years}', 
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
    def construction_schedule(self):
        """tuple[float] Construction investment fractions per year, starting from year 0. For example, for 50% capital investment in year 0 and 50% investment in year 1: (0.5, 0.5)."""
        return self._construction_schedule
    @construction_schedule.setter
    def construction_schedule(self, schedule):
        self._construction_schedule = np.array(schedule, dtype=float)
        self._start = len(schedule)
    
    @property
    def startup_months(self):
        return self._startup_time * 12.
    @startup_months.setter
    def startup_months(self, months):
        assert months <= 12., "startup time must be less than a year"
        self._startup_time = months/12.
    
    @property
    def sales(self):
        """Total sales (USD/yr)."""
        return self.system.sales
    @property
    def material_cost(self):
        """Total material cost (USD/yr)."""
        return self.system.material_cost
    @property
    def utility_cost(self):
        """Total utility cost (USD/yr)."""
        return self.system.utility_cost
    @property
    def purchase_cost(self):
        """Total purchase cost (USD)."""
        return self.system.purchase_cost
    @property
    def installed_equipment_cost(self):
        """Total installed cost (USD)."""
        return self.system.installed_equipment_cost
    @property
    def DPI(self):
        """Direct permanent investment."""
        return self._DPI(self.installed_equipment_cost)
    @property
    def TDC(self):
        """Total depreciable capital."""
        return self._TDC(self.DPI)
    @property
    def FCI(self):
        """Fixed capital investment."""
        return self._FCI(self.TDC)
    @property
    def TCI(self):
        """Total capital investment."""
        return (1. + self.WC_over_FCI)*self.FCI
    @property
    def FOC(self):
        """Fixed operating costs (USD/yr)."""
        return self._FOC(self.FCI)
    @property
    def VOC(self):
        """Variable operating costs (USD/yr)."""
        return self.material_cost + self.utility_cost
    @property
    def AOC(self):
        """Annual operating cost excluding depreciation (USD/yr)."""
        return self.FOC + self.VOC
    @property
    def working_capital(self):
        return self.WC_over_FCI * self.FCI
    
    @property
    def annual_depreciation(self):
        """Depreciation (USD/yr) equivalent to TDC dived by the the duration of the venture."""
        return self.TDC/(self.duration[1]-self.duration[0])

    @property
    def ROI(self):
        """Return on investment (1/yr) without accounting for annualized depreciation."""
        FCI = self.FCI
        net_earnings = (1-self.income_tax)*(self.sales-self._AOC(FCI))
        TCI = FCI*(1.+self.WC_over_FCI)
        return net_earnings/TCI
    @property
    def net_earnings(self):
        """Net earnings without accounting for annualized depreciation."""
        return (1-self.income_tax)*(self.sales-self.AOC)
    @property
    def PBP(self):
        """Pay back period (yr) without accounting for annualized depreciation."""
        FCI = self.FCI
        net_earnings = (1-self.income_tax)*(self.sales-self._AOC(FCI))
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
        C_D, C_FC, C_WC, D, L, LI, LP, LPl, C, S, T, I, NE, CF, DF, NPV, CNPV = data = np.zeros((17, length))
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
            L[:start] = loan = self.finance_fraction*(C_FC[:start]+C_WC[:start])
            f_interest = (1. + interest)
            LP[start:end] = solve_payment(loan.sum()/years * f_interest,
                                          loan, interest, years)
            loan_principal = 0
            for i in range(end):
                LI[i] = li = (loan_principal + L[i]) * interest 
                LPl[i] = loan_principal = loan_principal - LP[i] + li + L[i]
            taxable_cashflow = S - C - D - LP
            nontaxable_cashflow = D + L - C_FC - C_WC
        else:
            taxable_cashflow = S - C - D
            nontaxable_cashflow = D - C_FC - C_WC
        self._fill_tax_and_incentives(I, taxable_cashflow, nontaxable_cashflow, T)
        NE[:] = taxable_cashflow + I - T
        CF[:] = NE + nontaxable_cashflow
        DF[:] = 1/(1.+self.IRR)**self._get_duration_array()
        NPV[:] = CF*DF
        CNPV[:] = NPV.cumsum()
        DF *= 1e6
        data /= 1e6
        return pd.DataFrame(data.transpose(),
                            index=np.arange(self._duration[0]-start, self._duration[1]),
                            columns=cashflow_columns)
    @property
    def NPV(self):
        """Net present value."""
        taxable_cashflow, nontaxable_cashflow = self._taxable_and_nontaxable_cashflow_arrays()
        tax = np.zeros_like(taxable_cashflow)
        incentives = tax.copy()
        self._fill_tax_and_incentives(incentives, taxable_cashflow, nontaxable_cashflow, tax)
        cashflow = nontaxable_cashflow + taxable_cashflow + incentives - tax
        return NPV_at_IRR(self.IRR, cashflow, self._get_duration_array())
    
    def _AOC(self, FCI):
        """Return AOC at given FCI"""
        return self._FOC(FCI) + self.VOC
    
    def _taxable_and_nontaxable_cashflow_arrays(self):
        """Return taxable and nontaxable cash flows by year as a tuple[1d array, 1d array]."""
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
        return taxable_and_nontaxable_cashflows(system.unit_capital_costs if isinstance(system, bst.AgileSystem) else system.cost_units,
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
                                                self.lang_factor)
    
    def _fill_tax_and_incentives(self, incentives, taxable_cashflow, nontaxable_cashflow, tax):
        index = taxable_cashflow > 0.
        tax[index] = self.income_tax * taxable_cashflow[index]
    
    def _net_earnings_and_nontaxable_cashflow_arrays(self):
        taxable_cashflow, nontaxable_cashflow = self._taxable_and_nontaxable_cashflow_arrays()
        size = taxable_cashflow.size
        tax = np.zeros(size)
        incentives = tax.copy()
        self._fill_tax_and_incentives(incentives, taxable_cashflow, nontaxable_cashflow, tax)
        net_earnings = taxable_cashflow + incentives - tax
        return net_earnings, nontaxable_cashflow
    
    @property
    def cashflow_array(self):
        """[1d array] Cash flows by year."""
        return sum(self._net_earnings_and_nontaxable_cashflow_arrays())
    
    @property
    def net_earnings_array(self):
        """[1d array] Net earnings by year."""
        return self._net_earnings_and_nontaxable_cashflow_arrays()[0]
    
    def production_costs(self, products, with_annual_depreciation=True):
        """
        Return production cost of products [USD/yr].
        
        Parameters
        ----------
        products : Iterable[:class:`~thermosteam.Stream`]
            Main products of the system
        with_annual_depreciation=True : bool, optional
        
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
        
    def total_production_cost(self, products, with_annual_depreciation):
        """
        Return total production cost of products [USD/yr].
        
        Parameters
        ----------
        products : Iterable[:class:`~thermosteam.Stream`]
                    Main products of the system
        with_annual_depreciation=True : bool, optional
        
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
    
    def solve_IRR(self):
        """Return the IRR at the break even point (NPV = 0) through cash flow analysis."""
        IRR = self._IRR
        if not IRR or np.isnan(IRR) or IRR < 0.: IRR = self.IRR
        if not IRR or np.isnan(IRR) or IRR < 0.: IRR = 0.10
        args = (self.cashflow_array, self._get_duration_array())
        IRR = flx.aitken_secant(NPV_at_IRR,
                                IRR, 1.0001 * IRR + 1e-3, xtol=1e-6, ytol=10.,
                                maxiter=200, args=args, checkiter=False)
        self._IRR = IRR
        return IRR
        
    def solve_price(self, streams):
        """
        Return the price (USD/kg) of a stream(s) at the break even point (NPV = 0)
        through cash flow analysis. 
        
        Parameters
        ----------
        streams : tuple[:class:`~thermosteam.Stream`] or :class:`~thermosteam.Stream`
            Streams with variable selling price.
            
        """
        if isinstance(streams, bst.Stream): streams = [streams]
        system = self.system
        price2costs = [system._price2cost(i) for i in streams]
        price2cost = sum(price2costs)
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
        
    def solve_sales(self):
        """
        Return the required additional sales (USD) to reach the breakeven 
        point (NPV = 0) through cash flow analysis. 
        
        """
        discount_factors = (1 + self.IRR)**self._get_duration_array()
        sales_coefficients = np.ones_like(discount_factors)
        start = self._start
        sales_coefficients[:start] = 0
        w0 = self._startup_time
        sales_coefficients[self._start] =  w0*self.startup_VOCfrac + (1-w0)
        sales = self._sales
        if not np.isfinite(sales): sales = 0.
        taxable_cashflow, nontaxable_cashflow = self._taxable_and_nontaxable_cashflow_arrays()
        args = (taxable_cashflow, 
                nontaxable_cashflow, 
                sales_coefficients,
                discount_factors,
                self._fill_tax_and_incentives)
        x0 = sales
        x1 = 1.01 * sales + 10
        f = NPV_with_sales
        try:
            sales = flx.aitken_secant(f, x0, x1, xtol=10, ytol=1000.,
                                      maxiter=1000, args=args, checkiter=True)
        except Exception as e:
            y0 = f(x0, *args)
            y1 = f(x1, *args)
            if y0 == y1 and y0 < 0.:
                x0 = -y1 / self._years
            else:
                raise e
            sales = flx.aitken_secant(f, x0, x1, xtol=10, ytol=1000.,
                                      maxiter=1000, args=args, checkiter=True)
        self._sales = sales
        return sales
    
    def __repr__(self):
        return f'{type(self).__name__}({self.system.ID}, ...)'
    
    def _info(self):
        return (f'{type(self).__name__}: {self.system}\n'
                f' NPV: {self.NPV:,.0f} USD at {self.IRR:.1%} IRR')
    
    def show(self):
        """Prints information on unit."""
        print(self._info())
    _ipython_display_ = show