# -*- coding: utf-8 -*-
"""
Created on Mon Feb  4 19:38:37 2019

@author: Guest Group
"""

import pandas as pd
import numpy as np
from ._utils import wegstein_secant, aitken_secant
from copy import copy

__all__ = ('TEA',)

_DataFrame = pd.DataFrame
_array = np.array

# TODO: Add 'SL', 'DB', 'DDB', 'SYD', 'ACRS' and 'MACRS' functions to generate depreciation data

# %% Depreciation data

_MACRS = {'MACRS5':  _array([.2000, .3200, .1920,
                             .1152, .1152, .0576]),
          
          'MACRS7':  _array([.1429, .2449, .1749,
                             .1249, .0893, .0892,
                             .0893, .0446]),
          
          'MACRS10': _array([.1000, .1800, .1440,
                             .1152, .0922, .0737,
                             .0655, .0655, .0656,
                             .0655, .0328]),

          'MACRS15': _array([.0500, .0950, .0855,
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
    """Create a TEA object that can perform cash flow analysis on a System object.
    
        **Parameters**
        
            **system:** [System] Should contain feed and product streams.
            
        **Examples**
        
            :doc:`Techno-economic analysis of a biorefinery` 
    
    """
    __slots__ = ('system', '_cached', '_cashflow',
                 '_options', '_IRR_guess', '_cost_guess',
                 '_costarr', '_units')
    
    #: Default cash flow options
    _default = {'Lang factor': 3,
                'Operating days': 330,
                'IRR': 0.15,
                'Wage': 5e4,
                'Year': None,
                'Employees': 50,
                'Fringe benefits': 0.40,
                'Income tax': 0.35,
                'Property tax': 0.001,
                'Property insurance': 0.005,
                'Duration': 20,
                'Supplies': 0.20,
                'Maintenance': 0.01,
                'Administration': 0.005,
                'Working capital': 0.05,
                'Startup': 0,
                'Land': 0,
                'Royalties': 0,
                'Contingency': 0.18,
                'Depreciation': 'MACRS7',
                'Startup schedule': (0.4, 0.6),
                'Other recurring costs': 0,
                'Other fixed capital': 0}
    
    @property
    def options(self):
        """
        [dict] Options for cash flow analysis [-]:
            * **Lang factor:** Lang factor for getting fixed capital investment from total purchase cost. If no lang factor, estimate capital investment using bare module factors.
            * **Operating days:** (day).
            * **IRR:** Internal rate of return (fraction).
            * **Wage:** Wage per employee (USD/yr).
            * **Year:** Start year of venture.
            * **Employees:** Number of employees.
            * **Fringe benefits:** Cost of fringe benefits as a fraction of labor cost.
            * **Income tax:** Combined federal and state income tax rate (fraction).
            * **Property tax:** Fee as a fraction of fixed capital investment.
            * **Property insurance:** Fee as a fraction of fixed capital investment.
            * **Duration:** Duration of venture (years).
            * **Supplies:** Yearly fee as a fraction of labor costs.
            * **Maintenance:** Yearly fee as a fraction of fixed capital investment.
            * **Administration:** Yearly fee as a fraction of fixed capital investment.
            * **Working capital:** Fee as a fraction of fixed capital investment.
            * **Startup:** Cost of start up as a fraction of depreciable capital.
            * **Land:** Cost of land as a fraction of depreciable capital.
            * **Royalties:** Cost of royalties as a fraction of depreciable capital.
            * **Depreciation:** 'MACRS' + number of years (e.g. 'MACRS7').
            * **Startup schedule:** tuple of startup investment fractions per year, starting from year 0. For example, for 50% capital investment in year 0 and 50% investment in year 1: (0.5, 0.5).
            * **Contingency:** Cost as a fraction of direct permanent invesment.
            * **Other recurring costs:** Any additional recurring costs per year of operation.
            * **Other fixed capital:** Any additional investment not accounted for.
        """
        return self._options
    
    def __init__(self, system):
        self.system = system
        self._options = dict(**self._default)
        
        # [dict] Cached values from '_init' methods:
        #  {operating_days: flowrate_factor,
        #   year_duration: (duration_array, cashflow_data),
        #   depreciation_schedule: (start, Depreciation, D_len, schedule)}
        self._cached = {}
        
        #: Guess IRR for solve_IRR method
        self._IRR_guess = 0.15
        
        #: Guess stream cost for solve_price method
        self._cost_guess = 0
        
        units = system._costunits
        
        #: list[Unit] All unit operations considered
        self._units = units = sorted(units, key=lambda x: x.line)
        self._costarr = np.zeros([len(units), 2])
        
        system._TEA = self

    @property
    def utility_cost(self):
        """Total utility cost (USD/hr)."""
        return sum([u.utility_cost for u in self._units])
    @property
    def purchase_cost(self):
        """Total purchase cost (USD)."""
        return sum([u.purchase_cost for u in self._units])
    @property
    def intallation_cost(self):
        """Total installation cost (USD)."""
        return sum([u.installation_cost for u in self._units])

    def get_cashflow(self):
        """Return DataFrame of the cash flow analysis."""
        flow_factor, cashflow_info, depreciation_data = self._get_cached_data()
        cashflow_data, duration_array = cashflow_info
        parameters = self._calc_parameters(flow_factor)
        self._calc_cashflow(cashflow_data,
                            parameters[:-3],
                            depreciation_data)
        CF, DCF, CPV = cashflow_data[-3:]
        DCF[:] = CF/(1+self.options['IRR'])**duration_array
        CPV[:] = DCF.cumsum()
        return copy(self._cashflow)

    def _results(self):
        """Return a dictionary of summarized results of cash flow analysis."""
        flow_factor, cashflow_info, depreciation_data = self._get_cached_data()
        cashflow_data, duration_array = cashflow_info
        parameters = self._calc_parameters(flow_factor)
        self._calc_cashflow(cashflow_data,
                            parameters[:-3],
                            depreciation_data)
        NPV = self._calc_NPV(self.options['IRR'],
                             cashflow_data[-3],
                             duration_array)
        DC, FCI, WC, S, C, tax, UC, MC, L = parameters
        D = FCI/self._options['Duration']
        AOC = C + D
        TCI = FCI + WC
        net_earnings = (1-tax)*(S-AOC)
        ROI = net_earnings/TCI
        PBP = FCI/(net_earnings + D)
        r = {}
        r['Depreciable capital'] = DC
        r['Fixed capital investment'] = FCI
        r['Working capital'] = WC
        r['Total capital investment'] = TCI
        r['Depreciation'] = D
        r['Utility cost'] = UC
        r['Material cost'] = MC
        r['Sales'] = S
        r['Labor'] = L
        r['Annual operating cost'] = AOC
        r['Net present value'] = NPV
        r['Return on investment'] = ROI
        r['Pay back period'] = PBP
        return r

    def results(self, with_units=True):
        """Return results of techno-economic analysis as a DataFrame object if `with_units` is True or as a Series otherwise."""
        r = self._results()
        keys = []; addkey = keys.append
        vals = []; addval = vals.append
        if with_units:
            results_units = {'Pay back period': 'yr',
                             'Return on investment': '1/yr',
                             'Net present value': 'USD',
                             'Depreciable capital': 'USD',
                             'Fixed capital investment': 'USD',
                             'Total capital investment': 'USD',
                             'Depreciation': 'USD/yr',
                             'Annual operating cost': 'USD/yr',
                             'Working capital': 'USD',
                             'Utility cost': 'USD/yr',
                             'Material cost': 'USD/yr',
                             'Sales': 'USD/yr',
                             'Labor': 'USD/yr'}
            for ki, vi in r.items():
                addkey(ki)
                addval((results_units.get(ki, ''), vi))
            return pd.DataFrame(vals, keys, ('Units', 'Value'))
        else:
            return pd.Series(r)

    def NPV(self):
        """Calculate NPV by cash flow analysis."""
        flow_factor, cashflow_info, depreciation_data = self._get_cached_data()
        cashflow_data, duration_array = cashflow_info
        parameters = self._calc_parameters(flow_factor)
        self._calc_cashflow(cashflow_data,
                            parameters[:-3],
                            depreciation_data)
        return self._calc_NPV(self.options['IRR'],
                              cashflow_data[-3],
                              duration_array)
    
    def production_cost(self, *products):
        """Return production cost of products.
        
        **Parameters**
        
            ***products:** [Stream] Main products of the system
        
        .. Note::
           If there is more than one main product, The production cost is proportionally allocated to each of the main products with respect to their marketing values. The marketing value of each product is determined by the annual production multiplied by its selling price.
        """
        system = self.system
        sysfeeds = system.feeds
        sysproducts = system.products
        o = self.options
        flow_factor = 24*o['Operating days']
        # DPI_: Direct permanent investment (USD)
        # DC_: Depreciable capital (USD)
        # UC_: utility cost (USD/hr)
        costarr = self._costarr
        F_lang = o['Lang factor']
        if F_lang:
            costarr[:] = [(u.purchase_cost, u.utility_cost) for u in self._units]
            DC_, UC_ = costarr.sum(0) 
            DC_ *=  F_lang
        else:
            costarr[:] = [(u.installation_cost, u.utility_cost) for u in self._units]
            DPI_, UC_ = costarr.sum(0) 
            DC_ = DPI_ * (1 + o['Contingency'])
        FC_ = DC_ * (1 + o['Startup'] + o['Land'] + o['Royalties']) + o['Other fixed capital']
        MC_ = 0 # Material cost USD/hr
        CP_ = 0 # Coproducts USD/hr
        for s in sysfeeds:
            price = s.price
            if price: MC_ += price*s.massnet
        for s in sysproducts:
            if s not in products:
                price = s.price
                if price: CP_ += price*s.massnet
        # Multiply by flow_factor for USD/yr
        UC_ *= flow_factor
        MC_ *= flow_factor
        CP_ *= flow_factor
        fb = o['Fringe benefits'] + o['Supplies']
        f =  (o['Maintenance']
            + o['Administration']
            + o['Property tax']
            + o['Property insurance'])
        d = DC_/o['Duration'] # Depreciation
        L_ = o['Wage']*o['Employees']
        total_operating_cost = UC_ + (1+fb)*L_ + MC_ + f*FC_ + d + o['Other recurring costs']
        market_values = np.array([i.cost for i in products])
        weights = market_values/market_values.sum()
        production_cost = weights*total_operating_cost
        return production_cost
    
    def _get_cached_data(self):
        """Return cached data.
        
        **Return:** [tuple] including:
        
            **flow_factor:** [float] Factor to convert flow rate from kg/hr to kg/yr of operation.
            
            **cashflow_info:** [tuple] including:
                * cashflow_data: Cash flow data table
                * duration_array: Range from 1 to the end of the project length
        
            **depreciation_data:** [tuple] including:
                * start: Index of year when operation starts
                * Depreciation: Array of depreciation as a fraction of fixed cost
                * D_len: Lenght of depreciation
                * index_startup: Year index and fixed capital investment fractions
            
        """
        # Keys for cached data
        o = self._options
        cached = self._cached
        
        # Get and update cached data
        operating_days = o['Operating days']
        flow_factor = cached.get(operating_days)
        if not flow_factor:
            cached[operating_days] = flow_factor = 24*operating_days
        
        year_duration = (o['Year'], o['Duration'])
        cashflow_info = cached.get(year_duration)
        if not cashflow_info:
            year, duration = year_duration
            index = tuple(range(year, year + duration)) if year else None
            data = np.zeros((duration, 10))
            cashflow_data = data.transpose()
            self._cashflow = _DataFrame(data, index, _cashflow_columns, dtype=float)
            duration_array = _array(range(duration))
            cached[year_duration] = cashflow_info = (cashflow_data, duration_array)
        
        depreciation_schedule = (o['Depreciation'], o['Startup schedule'])
        depreciation_data = cached.get(depreciation_schedule)
        if not depreciation_data:
            depreciation, schedule = depreciation_schedule
            Depreciation = _MACRS[depreciation]
            start = len(schedule)
            end = start + len(Depreciation)
            schedule = np.array(schedule)
            depreciation_data = (start,
                                 Depreciation,
                                 end,
                                 schedule)
            cached[depreciation_schedule] = depreciation_data
        
        return flow_factor, cashflow_info, depreciation_data
    
    def _calc_parameters(self, flow_factor):
        """Return elementary cash flow parameters."""
        # Cash flow parameters (USD or USD/yr)
        # DC_: Depreciable capital
        # FC_: Fixed capital
        # WC_: Working capital
        # UC_: Utility cost
        # MC_: Material cost
        # S_: Sales
        # L_: Labor cost
        # C_: Annual operating cost (excluding depreciation)
        
        system = self.system
        feeds = system.feeds
        products = system.products
        o = self.options
        # DC_: Depreciable capital (USD)
        # UC_: utility cost (USD/hr)
        costarr = self._costarr
        F_lang = o['Lang factor']
        if F_lang:
            costarr[:] = [(u.purchase_cost, u.utility_cost) for u in self._units]
            DC_, UC_ = costarr.sum(0) 
            DC_ *=  F_lang
        else:
            costarr[:] = [(u.installation_cost, u.utility_cost) for u in self._units]
            DPI_, UC_ = costarr.sum(0) 
            DC_ = DPI_ * (1 + o['Contingency'])
        FC_ = DC_ * (1 + o['Startup'] + o['Land'] + o['Royalties']) + o['Other fixed capital']
        MC_ = 0 # Material cost USD/hr
        S_  = 0 # Sales USD/hr
        for s in feeds:
            price = s.price
            if price: MC_ += price*s.massnet
        for s in products:
            price = s.price
            if price: S_ += price*s.massnet
        # Multiply by flow_factor for USD/yr
        UC_ *= flow_factor
        MC_ *= flow_factor
        S_ *= flow_factor
        FC_ = DC_*(1 + o['Startup'] + o['Land'] + o['Royalties']) + o['Other fixed capital']
        fb = o['Fringe benefits'] + o['Supplies']
        f =  (o['Maintenance']
            + o['Administration']
            + o['Property tax']
            + o['Property insurance'])
        WC_ = o['Working capital']*FC_
        L_ = o['Wage']*o['Employees']
        C_ = UC_ + (1+fb)*L_ + MC_ + f*FC_ + o['Other recurring costs']
        return DC_, FC_, WC_, S_, C_, o['Income tax'], UC_, MC_, L_
    
    @staticmethod
    def _calc_cashflow(cashflow_data,
                       parameters,
                       depreciation_data):
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
        C_DC, D, C_FC, C_WC, C, S, NE, CF, DCF, CPV = cashflow_data
        DC_, FC_, WC_, S_, C_, tax = parameters
        start, Depreciation, end, schedule = depreciation_data
        
        # Calculate
        D[:] = 0.0
        C_DC[:start] = DC_*schedule
        C_FC[:start] = FC_*schedule
        D[start:end] = DC_*Depreciation
        C_WC[start] = WC_
        C_WC[-1] = -WC_
        C[start:] = C_
        S[start:] = S_
        NE[:] = (S - C - D)*(1 - tax)
        CF[:] = (NE + D) - C_FC - C_WC
    
    @staticmethod
    def _calc_NPV(IRR, CF, duration_array):
        """Return NPV at given IRR and cashflow data."""
        return (CF/(1+IRR)**duration_array).sum()

    def solve_IRR(self):
        """Return the IRR at the break even point (NPV = 0) through cash flow analysis."""
        # Calculate cashflow table
        flow_factor, cashflow_info, depreciation_data = self._get_cached_data()
        cashflow_data, duration_array = cashflow_info
        parameters = self._calc_parameters(flow_factor)
        self._calc_cashflow(cashflow_data,
                            parameters[:-3],
                            depreciation_data)
        # Solve
        IRR = aitken_secant(self._calc_NPV, self._IRR_guess, self._IRR_guess+1e-6,
                            args=(cashflow_data[-3], duration_array),
                            xtol=1e-6,
                            maxiter=200)
        self._IRR_guess = IRR if (0 < IRR < 1) else 0.15
        return IRR
    
    def solve_price(self, stream):
        """Return the price (USD/kg) of stream at the break even point (NPV = 0) through cash flow analysis. 
        
        **Parameters**
        
            **stream:** [Stream] Stream with variable selling price.
            
        """
        # Calculate cashflow table
        flow_factor, cashflow_info, depreciation_data = self._get_cached_data()
        cashflow_data, duration_array = cashflow_info
        parameters = self._calc_parameters(flow_factor)
        self._calc_cashflow(cashflow_data,
                            parameters[:-3],
                            depreciation_data)
        
        # Create function that adjusts cash flow with new stream price
        start = depreciation_data[0]
        tax = parameters[-4]
        IRR = self.options['IRR']
        cost_factor = stream.massnet*flow_factor*(1-tax)
        CF = cashflow_data[-3][start:]
        CF_copy = _array(CF)
        system = self.system
        if stream in system.feeds:
            adjust = CF_copy.__sub__
        elif stream in system.products:
            adjust = CF_copy.__add__
        else:
            raise ValueError(f"stream must be either a feed or a product of the system")
        
        # Solve
        calc_NPV = self._calc_NPV
        data_subset = cashflow_data[-3]
        def break_even_point(cost):
            CF[:] = adjust(cost)
            return calc_NPV(IRR, data_subset, duration_array)
        self._cost_guess = cost = aitken_secant(break_even_point,
                                                self._cost_guess,
                                                self._cost_guess+1e-6,
                                                xtol=1e-6,
                                                maxiter=200)
        return stream.price + cost/cost_factor
    
    def __repr__(self):
        return f'<{type(self).__name__}: {self.system.ID}>'
    
    def _info(self):
        r = self._results()
        out = f'{type(self).__name__}: {self.system.ID}\n'
        IRR = self.options['IRR']*100
        if r:
            NPV = r['Net present value']
            ROI = r['Return on investment']
            PBP = r['Pay back period']
            out += f' NPV: {NPV:.3g} USD at {IRR:.1f}% IRR\n'
            out += f' ROI: {ROI:.3g} 1/yr\n'
            out += f' PBP: {PBP:.3g} yr' 
        else:
            out += f' NPV: None\n'
            out += f' ROI: None\n'
            out += f' PBP: None'
        return out
    
    def show(self):
        """Prints information on unit."""
        print(self._info())
    _ipython_display_ = show
        
        
    