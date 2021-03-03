# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
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
from flexsolve import njitable
from math import floor
from warnings import warn

__all__ = ('TEA', 'CombinedTEA')


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


# %% Depreciation data

_MACRS = {'MACRS5': np.array([.2000, .3200, .1920,
                              .1152, .1152, .0576]),
          
          'MACRS7':  np.array([.1429, .2449, .1749,
                               .1249, .0893, .0892,
                               .0893, .0446]),
          
          'MACRS10': np.array([.1000, .1800, .1440,
                               .1152, .0922, .0737,
                               .0655, .0655, .0656,
                               .0655, .0328]),

          'MACRS15': np.array([.0500, .0950, .0855,
                               .0770, .0693, .0623,
                               .0590, .0590, .0591,
                               .0590, .0591, .0590,
                               .0591, .0590, .0591,
                               .0295]),

          'MACRS20': np.array([0.03750, 0.07219, 0.06677,
                               0.06177, 0.05713, 0.05285,
                               0.04888, 0.04522, 0.4462,
                               0.04461, 0.04462, 0.04461,
                               0.04462, 0.04461, 0.04462,
                               0.04461, 0.04462, 0.04461,
                               0.04462, 0.04461, 0.02231])}

# TODO: Add 'SL', 'DB', 'DDB', 'SYD', 'ACRS' and 'MACRS' functions to generate depreciation data

# %% Utilities

@njitable(cache=True)
def NPV_at_IRR(IRR, cashflow_array, duration_array):
    """Return NPV at given IRR and cashflow data."""
    return (cashflow_array/(1.+IRR)**duration_array).sum()

@njitable(cache=True)
def initial_loan_principal(loan, interest):
    principal = 0
    k = 1. + interest
    for i in loan:
        principal += i
        principal *= k
    return principal

@njitable(cache=True)
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

@njitable(cache=True)
def add_replacement_cost_to_cashflow_array(equipment_installed_cost, 
                                           equipment_lifetime, 
                                           cashflow_array, 
                                           venture_years,
                                           start):
    N_purchases = floor(venture_years / equipment_lifetime) 
    for i in range(1, N_purchases):
        cashflow_array[start + i * equipment_lifetime] += equipment_installed_cost

def add_all_replacement_costs_to_cashflow_array(unit, cashflow_array, 
                                                venture_years, start,
                                                lang_factor):
    equipment_lifetime = unit._equipment_lifetime
    if equipment_lifetime:
        if isinstance(equipment_lifetime, int):
             installed_cost = unit.purchase_cost * lang_factor if lang_factor else unit.installed_cost
             add_replacement_cost_to_cashflow_array(installed_cost, 
                                                    equipment_lifetime,
                                                    cashflow_array,
                                                    venture_years,
                                                    start)
        elif isinstance(equipment_lifetime, dict):
            if lang_factor:
                installed_costs = unit.purchase_costs.copy()
                for i in installed_costs:
                    installed_costs[i] *= lang_factor
            else:
                installed_costs = unit.installed_costs
            for name, installed_cost in installed_costs.items():
                lifetime = equipment_lifetime.get(name)
                if lifetime:
                    add_replacement_cost_to_cashflow_array(installed_cost, 
                                                           lifetime,
                                                           cashflow_array,
                                                           venture_years,
                                                           start)

@njitable(cache=True)
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
        unit_operations,
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
    for i in unit_operations:
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
    depreciation : str
        'MACRS' + number of years (e.g. 'MACRS7').
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
    
    Examples
    --------
    :doc:`tutorial/Techno-economic_analysis` 

    """
    __slots__ = ('system', 'income_tax', 'lang_factor', 'WC_over_FCI',
                 'finance_interest', 'finance_years', 'finance_fraction',
                 'feeds', 'products', '_construction_schedule', '_startup_time',
                 'startup_FOCfrac', 'startup_VOCfrac', 'startup_salesfrac',
                 'units', '_startup_schedule', '_operating_days',
                 '_operating_hours', '_duration', 
                 '_depreciation_array', '_depreciation', '_years',
                 '_duration', '_start',  'IRR', '_IRR', '_sales',
                 '_duration_array_cache')
    
    def __init_subclass__(self, isabstract=False):
        if isabstract: return
        for method in ('_DPI', '_TDC', '_FCI', '_FOC'):
            if not hasattr(self, method):
                raise NotImplementedError(f"subclass must implement a '{method}' method unless the 'isabstract' keyword argument is True")

    @staticmethod
    def like(system, other):
        """Create a TEA object from `system` with the same settings as `other`."""
        self = copy_(other)
        self.units = sorted([i for i in system.units if i._design or i._cost], key=lambda x: x.line)
        self.system = system
        self.feeds = system.feeds
        self.products = system.products
        system._TEA = self
        return self

    def __init__(self, system, IRR, duration, depreciation, income_tax,
                 operating_days, lang_factor, construction_schedule,
                 startup_months, startup_FOCfrac, startup_VOCfrac,
                 startup_salesfrac, WC_over_FCI, finance_interest,
                 finance_years, finance_fraction):
        self.duration = duration
        self.depreciation = depreciation
        self.construction_schedule = construction_schedule
        self.startup_months = startup_months
        self.operating_days = operating_days
        
        #: [float]  Internal rate of return (fraction).
        self.IRR = IRR
        
        #: [float] Combined federal and state income tax rate (fraction).
        self.income_tax = income_tax
        
        #: [float] Lang factor for getting fixed capital investment from total purchase cost. If no lang factor, estimate capital investment using bare module factors.
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
        
        #: list[Unit] All unit operations considered
        self.units = sorted([i for i in system.units if i._design or i._cost], key=lambda x: x.line)
        
        #: [System] System being evaluated.
        self.system = system
        
        #: set[Stream] All product streams.
        self.products = system.products
        
        #: set[Stream] All feed streams.
        self.feeds = system.feeds
        
        system._TEA = self

    def _DPI(self, installed_equipment_cost):
        return installed_equipment_cost # Default for backwards compatibility

    def _TDC(self, DPI):
        return DPI # Default for backwards compatibility

    @property
    def operating_days(self):
        """[float] Number of operating days per year."""
        return self._operating_days
    @operating_days.setter
    def operating_days(self, days):
        """[float] Number of operating days per year."""
        self._operating_days = days
        self._operating_hours = days*24
    
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
        """[str] 'MACRS' + number of years (e.g. 'MACRS7')."""
        return self._depreciation
    @depreciation.setter
    def depreciation(self, depreciation):
        try:
            self._depreciation_array = _MACRS[depreciation]
        except KeyError:
            raise ValueError(f"depreciation must be either 'MACRS5', 'MACRS7', 'MACRS10' or 'MACRS15 (not {repr(depreciation)})")
        self._depreciation = depreciation
    
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
    def utility_cost(self):
        """Total utility cost (USD/yr)."""
        return sum([u.utility_cost for u in self.units]) * self._operating_hours
    @property
    def purchase_cost(self):
        """Total purchase cost (USD)."""
        return sum([u.purchase_cost for u in self.units])
    @property
    def installed_equipment_cost(self):
        """Total installed cost (USD)."""
        return self.purchase_cost * self.lang_factor if self.lang_factor else sum([u.installed_cost for u in self.units])
    @property
    def DPI(self):
        """Direct permanent investment."""
        return self._DPI(self.purchase_cost * self.lang_factor if self.lang_factor else self.installed_equipment_cost)
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
        return self.WC_over_FCI * self.TDC
    @property
    def material_cost(self):
        """Annual material cost."""
        return sum([s.cost for s in self.feeds if s.price]) * self._operating_hours
    @property
    def annual_depreciation(self):
        """Depreciation (USD/yr) equivalent to FCI dived by the the duration of the venture."""
        return self.TDC/(self.duration[1]-self.duration[0])
    @property
    def sales(self):
        """Annual sales revenue."""
        return sum([s.cost for s in self.products if s.price]) * self._operating_hours
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
        depreciation = self._depreciation_array
        D[start:start+len(depreciation)] = TDC * depreciation
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
        lang_factor = self.lang_factor
        for i in self.units: add_all_replacement_costs_to_cashflow_array(i, C_FC, years, start, lang_factor)
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
        return NPV_at_IRR(self.IRR, self.cashflow_array, self._get_duration_array())
    
    def _AOC(self, FCI):
        """Return AOC at given FCI"""
        return self._FOC(FCI) + self.VOC
    
    def production_cost(self, products, with_annual_depreciation=True):
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
        market_values = np.array([i.cost for i in products])
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
        coproducts = self.products.copy()
        for i in products: coproducts.remove(i)
        coproduct_sales = sum([s.cost for s in coproducts if s.price]) * self._operating_hours
        if with_annual_depreciation:
            TDC = self.TDC
            annual_depreciation = TDC/(self.duration[1]-self.duration[0])
            AOC = self._AOC(self._FCI(TDC))
            return AOC - coproduct_sales + annual_depreciation
        else:
            return self.AOC - coproduct_sales
    
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
        D, C_FC, C_WC, Loan, LP, C, S = np.zeros((7, start+years))
        depreciation_array = self._depreciation_array
        D[start:start + depreciation_array.size] = TDC * depreciation_array
        WC = self.WC_over_FCI * FCI
        return taxable_and_nontaxable_cashflows(self.units,
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
    
    def _price2cost(self, stream):
        """Get factor to convert stream price to cost for cash flow in solve_price method."""
        F_mass = stream.F_mass
        if not F_mass: warn(RuntimeWarning(f"stream '{stream}' is empty"))
        return F_mass * self._operating_hours
    
    def _NPV_with_sales(self, sales, 
                       taxable_cashflow, 
                       nontaxable_cashflow,
                       sales_coefficients,
                       discount_factors):
        """Return NPV with an additional annualized sales."""
        taxable_cashflow = taxable_cashflow + sales * sales_coefficients
        tax = np.zeros_like(taxable_cashflow)
        incentives = tax.copy()
        self._fill_tax_and_incentives(incentives, taxable_cashflow, nontaxable_cashflow, tax)
        cashflow = nontaxable_cashflow + taxable_cashflow + incentives - tax
        return (cashflow/discount_factors).sum()
    
    def solve_price(self, stream):
        """
        Return the price (USD/kg) of stream at the break even point (NPV = 0)
        through cash flow analysis. 
        
        Parameters
        ----------
        stream : :class:`~thermosteam.Stream`
            Stream with variable selling price.
            
        """
        sales = self.solve_incentive()
        price2cost = self._price2cost(stream)
        if price2cost == 0.:
            return np.inf
        elif stream.sink:
            return stream.price - sales / price2cost
        elif stream.source:
            return stream.price + sales / price2cost
        else:
            raise ValueError("stream must be either a feed or a product")
    
    def solve_incentive(self):
        """
        Return the required incentive (USD) to reach the break even point (NPV = 0)
        through cash flow analysis. 
        
        """
        discount_factors = (1 + self.IRR)**self._get_duration_array()
        sales_coefficients = np.ones_like(discount_factors)
        start = self._start
        sales_coefficients[:start] = 0
        w0 = self._startup_time
        sales_coefficients[self._start] =  w0*self.startup_VOCfrac + (1-w0)
        sales = self._sales
        if not sales or np.isnan(sales): sales = 0.
        taxable_cashflow, nontaxable_cashflow = self._taxable_and_nontaxable_cashflow_arrays()
        args = (taxable_cashflow, 
                nontaxable_cashflow, 
                sales_coefficients,
                discount_factors)
        sales = flx.aitken_secant(self._NPV_with_sales,
                                  sales, 1.0001 * sales + 1e-4, xtol=1e-6, ytol=10.,
                                  maxiter=300, args=args, checkiter=False)
        self._sales = sales
        return sales
    
    def __repr__(self):
        return f'<{type(self).__name__}: {self.system.ID}>'
    
    def _info(self):
        return (f'{type(self).__name__}: {self.system.ID}\n'
                f' NPV: {self.NPV:,.0f} USD at {self.IRR:.1%} IRR')
    
    def show(self):
        """Prints information on unit."""
        print(self._info())
    _ipython_display_ = show
                
    
class CombinedTEA(TEA):
    """
    Create a CombinedTEA object that performs techno-economic analysis by 
    using data from many TEA objects that correspond to different areas
    of a biorefinery. A CombinedTEA object serves to accomodate for areas of
    a biorefinery with different assumptions. For example, an area of a
    biorefinery may not operate the same days a year, may have different 
    depreciation schedules, and may even be taxed differently than other 
    areas. Ultimately, a CombinedTEA object is an aggregation of TEA objects
    that work together to conduct techno-economic analysis of a whole biorefinery.
    
    Parameters
    ----------
    TEAs : Iterable[TEA]
        TEA objects used to conduct techno-economic analysis.
    IRR : float
        Internal rate of return.
    
    """
    _TDC = _FCI = _FOC = NotImplemented
    
    __slots__ = ('TEAs',)
    
    def __init__(self, TEAs, IRR, system=None):
        #: Iterable[TEA] All TEA objects for cashflow calculation
        self.TEAs = TEAs
        
        #: [float] Internal rate of return (fraction)
        self.IRR = IRR
        
        #: Guess IRR for solve_IRR method
        self._IRR = IRR
        
        #: Guess sales for solve_price method
        self._sales = 0
        
        if system: system._TEA = self
    
    def _get_duration_array(self):
        TEAs = self.TEAs
        if not TEAs: raise RuntimeError('no TEA objects available')
        first, *others = [(i._start, i._years) for i in TEAs]
        for i in others:
            if first != i: 
                raise RuntimeError('each TEA object must have the same number '
                                   'of operating years and the same number of '
                                   'construction years (as determined by the '
                                   'length of the contruction schedule)')
        return TEAs[0]._get_duration_array()
    
    @property
    def units(self):
        """All unit operations used for TEA."""
        return tuple(sum([i.units for i in self.TEAs], []))
    
    @property
    def operating_days(self):
        v_all = [i.operating_days for i in self.TEAs]
        v0, *vs = v_all
        if all([v0 == v for v in vs]): return v0
        else: raise AttributeError("TEAs don't have the same operating days")
    @operating_days.setter
    def operating_days(self, operating_days):
        vector = np.zeros(len(self.TEAs))
        vector[:] = operating_days
        for i, j in zip(self.TEAs, vector): i.operating_days = j
    
    @property
    def duration(self):
        v_all = [i.duration for i in self.TEAs]
        v0, *vs = v_all
        if all([np.all(v0 == v) for v in vs]): return v0
        else: raise AttributeError("TEAs don't have the same duration")
    @duration.setter
    def duration(self, duration):
        vector = np.zeros([len(self.TEAs), 2], dtype=int)
        vector[:] = duration
        for i, j in zip(self.TEAs, vector): i.duration = j
    
    @property
    def startup_months(self):
        v_all = [i.startup_months for i in self.TEAs]
        v0, *vs = v_all
        if all([v0 == v for v in vs]): return v0
        else: raise AttributeError("TEAs don't have the same startup months")
    @startup_months.setter
    def startup_months(self, startup_months):
        vector = np.zeros(len(self.TEAs))
        vector[:] = startup_months
        for i, j in zip(self.TEAs, vector): i.startup_months = j
    
    @property
    def income_tax(self):
        v_all = [i.income_tax for i in self.TEAs]
        v0, *vs = v_all
        if all([v0 == v for v in vs]): return v0
        else: raise AttributeError("TEAs don't have the same income tax")
    @income_tax.setter
    def income_tax(self, income_tax):
        vector = np.zeros(len(self.TEAs))
        vector[:] = income_tax
        for i, j in zip(self.TEAs, vector): i.income_tax = j
    
    def _taxable_and_nontaxable_cashflow_arrays(self):
        taxable_cashflow_arrays, nontaxable_cashflow_arrays = zip(*[i._taxable_and_nontaxable_cashflow_arrays() for i in self.TEAs])
        return sum(taxable_cashflow_arrays), sum(nontaxable_cashflow_arrays)
    
    @property
    def utility_cost(self):
        """Total utility cost (USD/yr)."""
        return sum([i.utility_cost for i in self.TEAs])
    
    @property
    def purchase_cost(self):
        """Total purchase cost (USD)."""
        return sum([i.purchase_cost for i in self.TEAs])
    
    @property
    def installed_equipment_cost(self):
        """Total installation cost (USD)."""
        return sum([i.installed_equipment_cost for i in self.TEAs])
    
    @property
    def DPI(self):
        """Direct permanent investment."""
        return sum([i.DPI for i in self.TEAs])
    
    @property
    def TDC(self):
        """Total depreciable capital."""
        return sum([i.TDC for i in self.TEAs])
    
    @property
    def FCI(self):
        """Fixed capital investment."""
        return sum([i.FCI for i in self.TEAs])
    
    @property
    def TCI(self):
        """Total capital investment."""
        return sum([i.TCI for i in self.TEAs])
    
    @property
    def FOC(self):
        """Fixed operating costs (USD/yr)."""
        return sum([i.FOC for i in self.TEAs])
    
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
        """Working capital."""
        return sum([i.working_capital for i in self.TEAs])
    
    @property
    def material_cost(self):
        """Annual material cost."""
        return sum([i.material_cost for i in self.TEAs])
    
    @property
    def annual_depreciation(self):
        """Depreciation (USD/yr) equivalent to FCI dived by the the duration of the venture."""
        return sum([i.annual_depreciation for i in self.TEAs])
    
    @property
    def sales(self):
        """Annual sales revenue."""
        return sum([i.sales for i in self.TEAs])
    
    @property
    def net_earnings(self):
        """Net earnings without accounting for annualized depreciation."""
        return sum([i.net_earnings for i in self.TEAs])
    
    @property
    def ROI(self):
        """Return on investment (1/yr) without accounting for annualized depreciation."""
        return sum([i.ROI for i in self.TEAs])
    
    @property
    def PBP(self):
        """Pay back period (yr) without accounting for annualized depreciation."""
        return self.FCI/self.net_earnings
    
    def get_cashflow_table(self):
        """Return DataFrame of the cash flow analysis."""
        IRR = self.IRR
        TEAs = self.TEAs
        IRRs = [i.IRR for i in TEAs]
        for i in TEAs: i.IRR = IRR
        TEA, *other_TEAs = TEAs
        table = TEA.get_cashflow_table()
        for i in other_TEAs:
            i_table = i.get_cashflow_table()
            if (i_table.index != table.index).any():
                raise NotImplementedError('cannot yet create cashflow table from TEAs with different venture years')
            table[:] += np.asarray(i_table)
        table['Net earnings'] = self.net_earnings_array / 1e6
        table['Cash flow'] = CF = self.cashflow_array  / 1e6
        table['Net present value (NPV)'] = NPV = CF * table['Discount factor']
        table['Cumulative NPV'] = NPV.cumsum()
        for IRR, TEA in zip(IRRs, TEAs): TEA.IRR = IRR
        return table    
    
    def production_cost(self, products, with_annual_depreciation=True):
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
        market_values = np.array([i.cost for i in products])
        total_market_value = market_values.sum()
        weights = market_values/total_market_value
        total_production_cost = 0
        for TEA in self.TEAs:
            total_production_cost += TEA.total_production_cost(products, with_annual_depreciation)
        return weights * total_production_cost
    
    def solve_IRR(self):
        """Return the IRR at the break even point (NPV = 0) through cash flow analysis."""
        IRR = self._IRR
        if not IRR or np.isnan(IRR) or IRR < 0.: IRR = self.IRR
        if not IRR or np.isnan(IRR) or IRR < 0.: IRR = 0.10
        args = (self.cashflow_array, self._get_duration_array())
        self._IRR = flx.aitken_secant(NPV_at_IRR,
                                      IRR, 1.0001 * IRR + 1e-3, xtol=1e-6, ytol=10.,
                                      maxiter=200, args=args, checkiter=False)
        return self._IRR
    
    def solve_price(self, stream, TEA=None):
        """
        Return the price (USD/kg) of stream at the break even point (NPV = 0)
        through cash flow analysis. 
        
        Parameters
        ----------
        stream : :class:`~thermosteam.Stream`
                 Stream with variable selling price.
        TEA : TEA, optional
              Stream should belong here.
            
        """
        if not TEA: TEA = self.TEAs[0]
        sales = self.solve_incentive(TEA)
        price2cost = TEA._price2cost(stream)
        if price2cost == 0:
            return np.inf
        elif stream.sink:
            return stream.price - sales / price2cost
        elif stream.source:
            return stream.price + sales / price2cost
        else:
            raise ValueError("stream must be either a feed or a product")
    
    def solve_incentive(self, TEA=None):
        """
        Return the required incentive (USD) to reach the break even point (NPV = 0)
        through cash flow analysis. 
        
        Parameters
        ----------
        TEA=None : TEA, optional
              Incentive will be added using settings from given TEA. Defaults to
              first TEA object in the `TEAs` attribute.
        
        """
        IRR = self.IRR
        if not TEA: TEA = self.TEAs[0]
        discount_factors = (1. + IRR)**TEA._get_duration_array()
        sales_coefficients = np.ones_like(discount_factors)
        start = TEA._start
        sales_coefficients[:start] = 0
        w0 = TEA._startup_time
        sales_coefficients[TEA._start] =  w0 * TEA.startup_VOCfrac + (1 - w0)
        sales = self._sales
        if not sales or np.isnan(sales): sales = 0.
        taxable_cashflow, nontaxable_cashflow = self._taxable_and_nontaxable_cashflow_arrays()
        args = (taxable_cashflow, 
                nontaxable_cashflow, 
                sales_coefficients,
                discount_factors)
        sales = flx.aitken_secant(TEA._NPV_with_sales,
                                  sales, 1.0001 * sales + 1e-4, xtol=1e-6, ytol=10.,
                                  maxiter=300, args=args, checkiter=False)
        self._sales = sales
        return sales
    
    def __repr__(self):
        return f'<{type(self).__name__}: {", ".join([i.system.ID for i in self.TEAs])}>'
    
    def _info(self):
        return (f'{type(self).__name__}: {", ".join([i.system.ID for i in self.TEAs])}\n'
                f' NPV: {self.NPV:,.0f} USD at {self.IRR:.1%} IRR')
    
    def show(self):
        """Prints information on unit."""
        print(self._info())
    _ipython_display_ = show
    
    
    
# def update_loan_principal(loan_principal, loan, loan_payment, interest):
#     principal = 0
#     for i, loan_i in enumerate(loan):
#         loan_principal[i] = principal = loan_i + principal * interest - loan_payment[i]