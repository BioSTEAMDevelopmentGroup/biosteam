# -*- coding: utf-8 -*-
"""
Created on Mon Feb  4 19:38:37 2019

@author: Guest Group
"""

import pandas as pd
import numpy as np
from ._utils import wegstein_secant, aitken_secant
from copy import copy as copy_

__all__ = ('TEA',)


# TODO: Add 'SL', 'DB', 'DDB', 'SYD', 'ACRS' and 'MACRS' functions to generate depreciation data

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
                               .0295])}


# %% Cash flow and results info

_cashflow_columns = ('Depreciable capital',
                     'Depreciation',
                     'Fixed capital',
                     'Working capital',
                     'Annual operating cost (excl. depr.)',
                     'Sales',
                     'Net earnings',
                     'Cash flow',
                     'Discounted cash flow',
                     'Cumulative cash flow')


# %% Techno-Economic Analysis

class TEA:
    """Abstract TEA class for cash flow analysis.
    
        **Abstract methods**
        
            **_TDC:** [function] Should take direct permanent investment as an argument and return total depreciable capital (e.g. _TDC(self, DPI) -> TDC).
            
            **_FCI:** [function] Should take total depreciable capital as an argument and return fixed capital investment (e.g. _FCI(self, TDC) -> FCI).
            
            **_AOC:** [function] Should take fixed capital investment as an arguments and return annual operation cost without depreciation (e.g. _AOC(self, FCI) -> AOC).
        
        **Parameters**
        
            **system:** [System] Should contain feed and product streams.
            
            **IRR:** [float]  Internal rate of return (fraction).
            
            **duration:** tuple[int, int] Start and end year of venture (e.g. (2018, 2038)).
            
            **depreciation:** [str] 'MACRS' + number of years (e.g. 'MACRS7').
            
            **operating_days:** [float] Number of operating days per year.
            
            **income_tax:** [float] Combined federal and state income tax rate (fraction).
            
            **lang_factor:** [float] Lang factor for getting fixed capital investment from total purchase cost. If no lang factor, estimate capital investment using bare module factors.
            
            **startup_schedule:** tuple[float] Startup investment fractions per year (e.g. (0.5, 0.5) for 50% capital investment in the first year and 50% investment in the second).
            
            **WC_over_FCI**: [float] Working capital as a fraction of fixed capital investment.
                        
        **Examples**
        
            :doc:`Techno-economic analysis of a biorefinery` 
    
    """
    
    __slots__ = ('system', 'income_tax', 'lang_factor', 'WC_over_FCI',
                 '_startup_schedule', 'IRR', '_IRR', '_cost', '_units',
                 '_operating_days', '_annual_factor', '_duration',
                 '_duration_array', '_cashflow_table',
                 '_cashflow_data', '_cashflow',
                 '_depreciation_array', '_depreciation',
                 '_duration', '_start')
    
    def __init_subclass__(self):
        if not hasattr(self, '_TDC'):
            raise NotImplementedError(f"subclass must implement '_TDC' method")
        if not hasattr(self, '_FCI'):
            raise NotImplementedError(f"subclass must implement '_FCI' method")
        if not hasattr(self, '_AOC'):
            raise NotImplementedError(f"subclass must implement '_AOC' method")

    def __init__(self, system, IRR, duration, depreciation, income_tax,
                 operating_days, lang_factor, startup_schedule, WC_over_FCI):
        self.IRR = IRR
        self.duration = duration
        self.depreciation = depreciation
        self.income_tax = income_tax
        self.operating_days = operating_days
        self.lang_factor = lang_factor
        self.startup_schedule = startup_schedule
        self.WC_over_FCI = WC_over_FCI
                
        #: Guess IRR for solve_IRR method
        self._IRR = 0.15
        
        #: Guess stream cost for solve_price method
        self._cost = 0
        
        #: list[Unit] All unit operations considered
        self._units = sorted(system._costunits, key=lambda x: x.line)
        
        self.system = system
        system._TEA = self

    @property
    def operating_days(self):
        """[float] Number of operating days per year."""
        return self._operating_days
    @operating_days.setter
    def operating_days(self, days):
        """[float] Number of operating days per year."""
        self._operating_days = days
        self._annual_factor = days*24
    
    @property
    def duration(self):
        """tuple[int, int] Start and end year of venture."""
        return self._duration
    @duration.setter
    def duration(self, duration):
        self._duration = duration
        start, end = duration
        index = tuple(range(start, end))
        years = end-start
        data = np.zeros((years, 10))
        self._cashflow_table = pd.DataFrame(data, index,
                                            _cashflow_columns,
                                            dtype=float)
        self._cashflow_data = data.transpose()
        self._cashflow = self._cashflow_data[-3]
        self._duration_array = np.array(range(years))
        
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
    def startup_schedule(self):
        """tuple[float] Startup investment fractions per year, starting from year 0. For example, for 50% capital investment in year 0 and 50% investment in year 1: (0.5, 0.5)."""
        return tuple(self._startup_schedule)
    @startup_schedule.setter
    def startup_schedule(self, schedule):
        self._startup_schedule = np.array(schedule, dtype=float)
        self._start = len(schedule)
    
    @property
    def utility_cost(self):
        """Total utility cost (USD/yr)."""
        return sum([u.utility_cost for u in self._units]) * self._annual_factor
    @property
    def purchase_cost(self):
        """Total purchase cost (USD)."""
        return sum([u.purchase_cost for u in self._units])
    @property
    def installation_cost(self):
        """Total installation cost (USD)."""
        return sum([u.installation_cost for u in self._units])
    @property
    def DPI(self):
        """Direct permanent investment."""
        return self.purchase_cost * self.lang_factor if self.lang_factor else self.installation_cost
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
    def AOC(self):
        """Annual operating cost excluding depreciation."""
        return self._AOC(self.FCI)
    @property
    def working_capital(self):
        return self.WC_over_FCI * self.TDC
    @property
    def material_cost(self):
        """Annual material cost."""
        return sum([s.cost for s in self.system.feeds if s.price]) * self._annual_factor
    @property
    def annualized_depreciation(self):
        """Depreciation (USD/yr) equivalent to FCI dived by the the duration of the venture."""
        return self.FCI * (self.duration[1]-self.duration[0])
    @property
    def sales(self):
        """Annual sales revenue."""
        return sum([s.cost for s in self.system.products if s.price]) * self._annual_factor
    @property
    def ROI(self):
        """Return on investment (1/yr) accounting for annualized depreciation."""
        FCI = self.FCI
        depreciation = FCI*(self.duration[1]-self.duration[0])
        net_earnings = (1-self.income_tax)*(self.sales-self._AOC(FCI)-depreciation)
        TCI = FCI*(1.+self.WC_over_FCI)
        return net_earnings/TCI
    @property
    def PBP(self):
        """Pay back period (yr) accounting for annualized depreciation."""
        FCI = self.FCI
        depreciation = FCI*(self.duration[1]-self.duration[0])
        net_earnings = (1-self.income_tax)*(self.sales-self._AOC(FCI)-depreciation)
        return FCI/(net_earnings + depreciation)

    def get_cashflow(self, copy=True):
        """Return DataFrame of the cash flow analysis."""
        self._update_cashflow()
        CF, DCF, CPV = self._cashflow_data[-3:]
        DCF[:] = CF/(1.+self.IRR)**self._duration_array
        CPV[:] = DCF.cumsum()
        return copy_(self._cashflow_table) if copy else self._cashflow_table

    @property
    def NPV(self):
        """Net present value."""
        self._update_cashflow()
        return self._NPV_at_IRR(self.IRR)
    
    def production_cost(self, *products):
        """Return production cost of products.
        
        **Parameters**
        
            ***products:** [Stream] Main products of the system
        
        .. Note::
           If there is more than one main product, The production cost is proportionally allocated to each of the main products with respect to their marketing values. The marketing value of each product is determined by the annual production multiplied by its selling price.
        """
        market_values = np.array([i.cost for i in products])
        weights = market_values/market_values.sum()
        return weights*self.AOC
    
    def _update_cashflow(self):
        """Perform cash flow analysis and return net present value."""
        # Cash flow data and parameters
        # C_DC: Depreciable capital
        # C_FC: Fixed capital
        # C_WC: Working capital
        # D: Depreciation
        # C: Annual operating cost (excluding depreciation)
        # S: Sales
        # NE: Net earnings
        # CF: Cash flow
        # DCF: Discounted cash flow
        # CPV: Cumulative present value
        TDC = self.TDC
        FCI = self._FCI(TDC)
        WC = self.WC_over_FCI * FCI
        depreciation = self._depreciation_array
        schedule = self._startup_schedule
        C_DC, D, C_FC, C_WC, C, S, NE, CF, DCF, CPV = self._cashflow_data
        D[:] = 0.0
        start = self._start
        C_DC[:start] = TDC*schedule
        C_FC[:start] = FCI*schedule
        D[start+1:start+len(depreciation)+1] = TDC*depreciation
        C_WC[start] = WC
        C_WC[-1] = -WC
        C[start:] = self._AOC(FCI)
        S[start:] = self.sales
        NE[:] = (S - C - D)*(1 - self.income_tax)
        CF[:] = (NE + D) - C_FC - C_WC
    
    def _NPV_at_IRR(self, IRR):
        """Return NPV at given IRR and cashflow data."""
        return (self._cashflow/(1+IRR)**self._duration_array).sum()

    def _NPV_with_cost(self, cost):
        """Return NPV with an additional anualized cost."""
        cashflow = self._cashflow.copy()
        cashflow[self._start:] -= cost
        return (cashflow/(1+self.IRR)**self._duration_array).sum()

    def solve_IRR(self):
        """Return the IRR at the break even point (NPV = 0) through cash flow analysis."""
        self._update_cashflow()
        try:
            self._IRR = aitken_secant(self._NPV_at_IRR,
                                      self._IRR, self._IRR+1e-6,
                                      xtol=1e-6, maxiter=200)
        except:
            self._IRR = aitken_secant(self._NPV_at_IRR,
                                      0.15, 0.15001,
                                      xtol=1e-6, maxiter=200)
        return self._IRR
    
    def solve_price(self, stream):
        """Return the price (USD/kg) of stream at the break even point (NPV = 0) through cash flow analysis. 
        
        **Parameters**
        
            **stream:** [Stream] Stream with variable selling price.
            
        """
        self._update_cashflow()
        try:
            self._cost = aitken_secant(self._NPV_with_cost,
                                       self._cost, self._cost+1e-6,
                                       xtol=1e-6, maxiter=200)
        except:
            self._cost = aitken_secant(self._NPV_with_cost, 0, 1e-6,
                                       xtol=1e-6, maxiter=200)        
        price2cost = stream.massnet*self._annual_factor*(1-self.income_tax)
        if stream.sink:
            return stream.price + self._cost/price2cost
        elif stream.source:
            print(self._cost/price2cost)
            return stream.price - self._cost/price2cost
        else:
            raise ValueError(f"stream must be either a feed or a product")
    
    def __repr__(self):
        return f'<{type(self).__name__}: {self.system.ID}>'
    
    def _info(self):
        return (f'{type(self).__name__}: {self.system.ID}\n'
                f' NPV: {self.NPV:.3g} USD at {self.IRR:.1%} IRR\n'
                f' ROI: {self.ROI:.3g} 1/yr\n'
                f' PBP: {self.PBP:.3g} yr')
    
    def show(self):
        """Prints information on unit."""
        print(self._info())
    _ipython_display_ = show
                