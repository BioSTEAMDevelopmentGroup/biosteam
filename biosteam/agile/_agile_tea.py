# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from .._tea import TEA
from warnings import warn
from ._agile_scenario import AgileScenario
import biosteam as bst

__all__ = ('AgileTEA',)


# %% Agile TEA


class AgileTEA(AgileScenario):
    """
    Abstract AgileTEA class for cash flow analysis.
    
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
    IRR : float
        Internal rate of return (fraction).
    duration : tuple[int, int]
        Start and end year of venture (e.g. (2018, 2038)).
    depreciation : str
        'MACRS' + number of years (e.g. 'MACRS7').
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

    """
    __slots__ = ('scenarios', 'income_tax', 'lang_factor', 'WC_over_FCI',
                 'finance_interest', 'finance_years', 'finance_fraction',
                 'startup_FOCfrac', 'startup_VOCfrac', 'startup_salesfrac',
                 '_construction_schedule', '_startup_time',
                 '_startup_schedule', '_duration', 
                 '_depreciation_array', '_depreciation', '_years',
                 '_duration', '_start',  'IRR', '_IRR', '_sales',
                 '_duration_array_cache', 'units', 'operating_hours',
                 'utility_cost', 'material_cost', 'sales', 'flow_rates',
                 'feeds', 'products')
    
    def create_scenario(self, system):
        return system.get_scenario_costs()
    
    def compile_scenarios(self, scenarios):
        units = set(sum([list(i.unit_capital_costs) for i in scenarios], []))
        unit_scenarios = {i: [] for i in units}
        for scenario in scenarios:
            unit_capital_costs = scenario.unit_capital_costs
            for i, j in unit_capital_costs.items(): unit_scenarios[i].append(j)
        self.units = [i.get_agile_capital_costs(j) for i, j in unit_scenarios.items()]
        self.operating_hours = sum([i.operating_hours for i in scenarios])
        self.utility_cost = sum([i.utility_cost for i in scenarios])
        self.material_cost = sum([i.material_cost for i in scenarios])
        self.sales = sum([i.sales for i in scenarios])
        self.flow_rates = flow_rates = {}
        self.feeds = set(sum([i.feeds for i in scenarios], []))
        self.products = set(sum([i.products for i in scenarios], []))
        for scenario in scenarios:
            for stream, F_mass in scenario.flow_rates.items():
                if stream in flow_rates: flow_rates[stream] += F_mass
                else: flow_rates[stream] = F_mass
    
    # TODO: Add 'SL', 'DB', 'DDB', 'SYD', 'ACRS' and 'MACRS' functions to generate depreciation data
    #: dict[str, 1d-array] Available depreciation schedules.
    depreciation_schedules = TEA.depreciation_schedules
    
    #: dict[str, float] Investment site factors used to multiply the total permanent 
    #: investment (TPI), also known as total fixed capital (FCI), to 
    #: account for locality cost differences based on labor availability,
    #: workforce efficiency, local rules, etc.
    investment_site_factors = TEA.investment_site_factors

    def __init_subclass__(cls, isabstract=False):
        if isabstract: return
        for method in ('_DPI', '_TDC', '_FCI', '_FOC'):
            if not hasattr(cls, method):
                breakpoint()
                raise NotImplementedError(
                    f"subclass must implement a '{method}' method unless the "
                     "'isabstract' keyword argument is True"
                )

    def __init__(self, IRR, duration, depreciation, income_tax,
                 lang_factor, construction_schedule, startup_months, 
                 startup_FOCfrac, startup_VOCfrac,
                 startup_salesfrac, WC_over_FCI, finance_interest,
                 finance_years, finance_fraction):
        self.duration = duration
        self.depreciation = depreciation
        self.construction_schedule = construction_schedule
        self.startup_months = startup_months
        
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
        
    _get_duration_array = TEA._get_duration_array
    _fill_depreciation_array = TEA._fill_depreciation_array
    _get_duration = TEA._get_duration
    _taxable_and_nontaxable_cashflow_arrays = TEA._taxable_and_nontaxable_cashflow_arrays
    _net_earnings_and_nontaxable_cashflow_arrays = TEA._net_earnings_and_nontaxable_cashflow_arrays
    _fill_tax_and_incentives = TEA._fill_tax_and_incentives
    _DPI = TEA._DPI
    _TDC = TEA._TDC
    _AOC = TEA._AOC
    duration = TEA.duration
    depreciation = TEA.depreciation
    construction_schedule = TEA.construction_schedule
    startup_months = TEA.startup_months
    working_capital = TEA.working_capital
    annual_depreciation = TEA.annual_depreciation
    net_earnings = TEA.net_earnings
    purchase_cost = TEA.purchase_cost
    installed_equipment_cost = TEA.installed_equipment_cost
    cashflow_array = TEA.cashflow_array
    net_earnings_array = TEA.net_earnings_array
    production_costs = TEA.production_costs
    total_production_cost = TEA.total_production_cost
    solve_IRR = TEA.solve_IRR
    solve_price = TEA.solve_price
    solve_sales = TEA.solve_sales
    DPI = TEA.DPI
    TDC = TEA.TDC
    FCI = TEA.FCI
    TCI = TEA.TCI
    FOC = TEA.FOC
    VOC = TEA.VOC
    AOC = TEA.AOC
    ROI = TEA.ROI
    PBP = TEA.PBP
    NPV = TEA.NPV
    
    def market_value(self, stream):
        """Return the market value of a stream [USD/yr]."""
        return self.flow_rates[stream] * stream.price
    
    def _price2cost(self, stream):
        """Get factor to convert stream price to cost for cash flow in solve_price method."""
        if stream in self.flow_rates:
            F_mass = self.flow_rates[stream] 
        else:
            F_mass = 0.
        if not F_mass: warn(RuntimeWarning(f"stream '{stream}' is empty"))
        if stream in self.products:
            return F_mass
        elif stream in self.feeds:
            return - F_mass
        else:
            raise ValueError("stream must be either a feed or a product")
    
    def _info(self):
        return (f'{type(self).__name__}: \n'
                f' NPV: {self.NPV:,.0f} USD at {self.IRR:.1%} IRR')
    
    def show(self):
        """Prints information on unit."""
        print(self._info())
    _ipython_display_ = show
    