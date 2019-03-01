# -*- coding: utf-8 -*-
"""
Created on Mon Feb  4 19:38:37 2019

@author: Guest Group
"""

from . import pd, np
from bookkeep import SmartBook
from scipy.optimize import newton
from biosteam.utils import CS

_DataFrame = pd.DataFrame
_array = np.array
_asarray = np.asarray

# TODO: Add 'SL', 'DB', 'DDB', 'SYD', 'ACRS' and 'MACRS' functions to generate depreciation data

# %% Depreciation data

_MACRS = {'MACRS5':  _array([.2000, .3200, .1920,
                             .1152, .1152, .0576]),
          
          'MACRS7':  _array([.1429, .2449, .1749,
                             .1249, .0893, .0892,
                             .0883, .0446]),
          
          'MACRS10': _array([.1000, .1800, .1440,
                             .1152, .0922, .0737,
                             .0655, .0655, .0328]),

          'MACRS15': _array([.0500, .0950, .0855,
                             .0770, .0693, .0623,
                             .0590, .0590, .0591,
                             .0590, .0591, .0590,
                             .0591, .0590, .0591,
                             .0295])}


# %% Cash flow and results info

_cashflow_columns = ('Depreciable capital',
                     'Working capital',
                     'Depreciation',
                     'Annual operating cost (excl. depr.)',
                     'Sales',
                     'Net earnings',
                     'Cash flow',
                     'Discounted cash flow',
                     'Cumulative cash flow')

_results_units = {'Pay back period': 'yr',
                  'Return on investment': '1/yr',
                  'Net present value': 'USD',
                  'Fixed capital investment': 'USD',
                  'Total capital investment': 'USD',
                  'Depreciation': 'USD/yr',
                  'Annual operating cost': 'USD/yr',
                  'Working capital': 'USD',
                  'Utility cost': 'USD/yr',
                  'Material cost': 'USD/yr',
                  'Sales': 'USD/yr',
                  'Labor': 'USD/yr'}


# %% Techno-Economic Analysis

class TEA:
    """Create a TEA object that can perform cash flow analysis on a System object.
    
        **Parameters**
        
            **system:** [System] Should contain feed and product streams.
            
        **Examples**
        
            :doc:`TEA Example` 
    
    """
    __slots__ = ('results', 'system', 'cashflow', '_cached',
                 '_options', '_IRR_guess', '_cost_guess')
    
    #: Default cash flow options
    _default = {'Lang factor': 4.37,
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
                'Mantainance': 0.01,
                'Administration': 0.005,
                'Working capital': 0.05,
                'Startup cost': 0,
                'Delivery': 0,
                'Land': 0,
                'Depreciation': 'MACRS7',
                'Startup schedule': (0.4, 0.6)}
    
    @property
    def options(self):
        """
        [dict] Options for cash flow analysis [-]:
            * **Lang factor:** Used to get fixed capital investment from total purchase cost
            * **Operating days:** (day)
            * **IRR:** Internal rate of return (fraction)
            * **Wage:** Wage per employee (USD/yr)
            * **Year:** Start year of venture
            * **Employees:** Number of employees
            * **Fringe benefits:** Cost of fringe benefits as a fraction of labor cost
            * **Income tax:** Combined federal and state income tax rate (fraction)
            * **Property tax:** Fee as a fraction of fixed capital investment
            * **Property insurance:** Fee as a fraction of fixed capital investment
            * **Duration:** Duration of venture (years)
            * **Supplies:** Yearly fee as a fraction of labor costs
            * **Mantainance:** Yearly fee as a fraction of fixed capital investment
            * **Administration:** Yearly fee as a fraction of fixed capital investment
            * **Working capital:** Fee as a fraction of fixed capital investment
            * **Startup cost:** Cost of start up as a fraction of fixed capital investment
            * **Delivery:** Delibery of equipment fee as a fraction of fixed capital investment 
            * **Land:** Cost of land as a fraction of fixed capital investment
            * **Depreciation:** 'MACRS' + number of years (e.g. 'MACRS7').
            * **Startup schedule:** tuple of startup investment fractions per year, starting from year 0. For example, for 50% capital investment in year 0 and 50% investment in year 1: array([0.5, 0.5]).
        """
        return self._options
    
    def __init__(self, system):
        self.system = system
        self._options = dict(**self._default)
        
        # [dict] Cached values from '_init' methods:
        #  {operating_days: flowrate_factor,
        #   year_duration: (duration_array, cashflow_data),
        #   depreciation_schedule: (start, Depreciation, D_len, index_schedule)}
        self._cached = {}
        
        #: [dict] Summarized results of cash flow analysis
        self.results = SmartBook(_results_units)
        
        #: [DataFrame] Cash flow table
        self.cashflow = None
        
        #: Guess IRR for solve_IRR method
        self._IRR_guess = None
        
        #: Guess stream cost for solve_price method
        self._cost_guess = None
        system.TEA = self

    def __call__(self):
        """Perform cash flow analysis and update the "results" and "cashflow" attributes."""
        flow_factor, cashflow_info, depreciation_data = self._get_cached_data()
        cashflow_data, duration_array = cashflow_info
        parameters = self._calc_parameters(flow_factor)
        self._calc_cashflow(cashflow_data,
                            parameters[:-3],
                            depreciation_data)
        self._calc_NPV(self.options['IRR'],
                       cashflow_data[-3:],
                       duration_array)
        self._update_results(parameters, cashflow_data[-1, -1])
    
    def _get_cached_data(self):
        """Return cached data.
        
        **Return**
        
            **flow_factor:** [float] Factor to convert flow rate from kg/hr to kg/yr of operation.
            
            **cashflow_info:** tuple[array] including:
                * cashflow_data: Cash flow data table
                * duration_array: Range from 1 to the end of the project length
        
            **depreciation_data:** [tuple] including:
                * start: Index of year when operation starts
                * Depreciation: Array of depreciation as a fraction of fixed cost
                * D_len: Lenght of depreciation
                * index_startup: Year index and fixed capital investment fractions
            
        """
        # Keys for cached data
        o = self.options
        cached = self._cached
        year_duration = (o['Year'], o['Duration'])
        depreciation_schedule = (o['Depreciation'], o['Startup schedule'])
        operating_days = o['Operating days']
        
        # Get and update cached data
        flow_factor = cached.get(operating_days)
        cashflow_info = cached.get(year_duration)
        depreciation_data = cached.get(depreciation_schedule)
        if not flow_factor:
            cached[operating_days] = flow_factor = 24*operating_days
        if not cashflow_info:
            year, duration = year_duration
            index = tuple(range(year, year + duration)) if year else None
            data = np.zeros((duration, 9))
            self.cashflow = _DataFrame(data, index, _cashflow_columns, dtype=float)
            cashflow_data = _asarray(self.cashflow).transpose()
            duration_array = _array(range(duration))
            cached[year_duration] = cashflow_info = (cashflow_data, duration_array)
        if not depreciation_data:
            depreciation, schedule = depreciation_schedule
            Depreciation = _MACRS[depreciation]
            start = len(schedule)
            end = start + len(Depreciation)
            index_schedule = tuple(enumerate(schedule))
            depreciation_data = (start,
                                 Depreciation,
                                 end,
                                 index_schedule)
            cached[depreciation_schedule] = depreciation_data
        
        return flow_factor, cashflow_info, depreciation_data
    
    def _calc_parameters(self, flow_factor):
        """Return elementary cash flow parameters."""
        # Cash flow parameters
        # FC_: Fixed capital cost
        # WC_: Working capital cost
        # UC_: Utility cost
        # MC_: Material cost
        # S_: Sales
        # L_: Labor cost
        # C_: Annual operating cost (excluding depreciation)
        
        system = self.system
        units = system.units.union(system.offsite_units)
        feeds = system.feeds
        products = system.products
        FC_ = 0 # Fixed capital USD
        UC_ = 0 # Utility cost USD/hr
        MC_ = 0 # Material cost USD/hr
        S_ = 0  # Sales USD/hr
        for u in units:
            c = iter(u.results['Summary'].values())
            FC_ += next(c)
            UC_ += next(c)
        for s in feeds:
            price = s.price
            if price: MC_ += price*(s._mol*s._MW).sum()
        for s in products:
            price = s.price
            if price: S_ += price*(s._mol*s._MW).sum()
        o = self.options
        UC_ *= flow_factor
        MC_ *= flow_factor
        S_ *= flow_factor
        FC_ *= (1 + o['Startup cost'] + o['Delivery'] + o['Land'])*o['Lang factor']
        wcf = o['Working capital']
        fb = o['Fringe benefits'] + o['Supplies']
        f =  (o['Mantainance']
            + o['Administration']
            + o['Property tax']
            + o['Property insurance'])
        WC_ = wcf*FC_
        L_ = o['Wage']*o['Employees']
        C_ = UC_ + (1+fb)*L_ + MC_ + f*FC_
        return FC_, WC_, S_, C_, o['Income tax'], UC_, MC_, L_
    
    @staticmethod
    def _calc_cashflow(cashflow_data,
                       parameters,
                       depreciation_data):
        """Perform cash flow analysis and return net present value."""
        # Cash flow data and parameters
        # C_DC: Depreciable capital cost
        # C_WC: Working capital cost
        # D: Depreciation
        # C: Annual operating cost (excluding depreciation)
        # S: Sales
        # NE: Net earnings
        # CF: Cash flow
        # DCF: Discounted cash flow
        # CPV: Cumulative present value
        start, Depreciation, end, index_schedule = depreciation_data
        C_DC, C_WC, D, C, S, NE, CF, DCF, CPV = cashflow_data
        FC_, WC_, S_, C_, tax = parameters
        
        # Calculate
        D[:] = 0.0
        for i, f in index_schedule:
            C_DC[i] = f*FC_
        D[start:end] = FC_*Depreciation
        C_WC[start] = WC_
        C_WC[-1] = -WC_
        C[start:] = C_
        S[start:] = S_
        NE[:] = (S - C - D)*(1 - tax)
        CF[:] = (NE + D) - C_DC - C_WC
    
    @staticmethod
    def _calc_NPV(IRR, CF_subset, duration_array):
        """Return NPV at given IRR and cashflow data."""
        CF, DCF, CPV = CF_subset
        DCF[:] = CF/(1+IRR)**duration_array
        CPV[:] = DCF.cumsum()
        return CPV[-1]
        
    def _update_results(self, parameters, NPV):
        """Update results attribute."""
        FCI, WC, S, C, tax, UC, MC, L = parameters
        D = FCI/self.options['Duration']
        AOC = C + D
        TCI = FCI + WC
        net_earnings = (1-tax)*(S-AOC)
        ROI = net_earnings/TCI
        PBP = FCI/(net_earnings + D)
        
        r = self.results
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

    def solve_IRR(self, update=True):
        """Return the IRR at the break even point (NPV = 0) through cash flow analysis.
        
        **Parameters**
        
            **update:** [bool] If True, update IRR, cashflow, and results.
           
        """
        # Calculate cashflow table
        flow_factor, cashflow_info, depreciation_data = self._get_cached_data()
        cashflow_data, duration_array = cashflow_info
        parameters = self._calc_parameters(flow_factor)
        self._calc_cashflow(cashflow_data,
                            parameters[:-3],
                            depreciation_data)
        
        # Setup arguments for solver
        data_subset = cashflow_data[-3:]
        old_data_subset = _array(data_subset)
        guess = self._IRR_guess
        IRR_guess = guess if guess else self.options['IRR']
        args = (data_subset, duration_array)
        
        # Solve
        self._IRR_guess = IRR = newton(self._calc_NPV, IRR_guess, args=args)
        if update:
            self.options['IRR'] = IRR
            self._calc_cashflow(cashflow_data,
                                parameters[:-3],
                                depreciation_data)
            self._update_results(parameters, data_subset[-1, -1])
        else:
            data_subset[:] = old_data_subset
        return IRR
    
    def solve_price(self, stream, update=True):
        """Return the price (USD/kg) of stream at the break even point (NPV = 0) through cash flow analysis. 
        
        **Parameters**
        
            **stream:** [Stream] Stream with variable selling price.
            
            **update:** [bool] If True, update stream price, cashflow, and results.
            
        """
        # Exclude stream cost when calculating cash flow
        data = _asarray(self.cashflow)
        data_old = _array(data)
        price_old = stream.price
        stream.price = 0
        
        # Calculate cashflow table
        flow_factor, cashflow_info, depreciation_data = self._get_cached_data()
        cashflow_data, duration_array = cashflow_info
        parameters = self._calc_parameters(flow_factor)
        self._calc_cashflow(cashflow_data,
                            parameters[:-3],
                            depreciation_data)
        
        # Create function that adjusts cash flow with new stream price
        data_subset = cashflow_data[-3:]
        args = (data_subset, duration_array)
        start = depreciation_data[0]
        tax = parameters[-4]
        IRR = self.options['IRR']
        calc_NPV = self._calc_NPV
        cost_factor = stream.massnet*flow_factor*(1-tax)
        CF = data_subset[0][start:]
        CF_copy = _array(CF)
        system = self.system
        if stream in system.feeds:
            adjust = CF_copy.__sub__
        elif stream in system.products:
            adjust = CF_copy.__add__
        else:
            raise ValueError(f"Stream '{stream.ID}' is not a feed nor a product of the system")
        
        def break_even_point(cost):
            CF[:] = adjust(cost)
            return calc_NPV(IRR, *args)
        
        # Solve
        guess = self._cost_guess
        cost_guess = guess if guess else price_old*cost_factor
        self._cost_guess = cost = newton(break_even_point, cost_guess)
        price = cost/cost_factor
        if update:
            parameters = self._calc_parameters(flow_factor)
            self._calc_cashflow(cashflow_data,
                                parameters[:-3],
                                depreciation_data)
            self._update_results(parameters, data_subset[-1, -1])
            stream.price = price
        else:
            data[:] = data_old
            stream.price = price_old
        return price
    
    def __repr__(self):
        return f'<{type(self).__name__}: {self.system.ID}>'
    
    def _info(self):
        r = self.results
        out = f'{type(self).__name__}: {self.system.ID}\n'
        IRR = self.options['IRR']*100
        if r:
            NPV = r['Net present value']
            ROI = r['Return on investment']
            PBP = r['Pay back period']
            out += f' NPV: {NPV:.3g} ' + CS.dim(f'USD at %{IRR:.1f} IRR') + '\n'
            out += f' ROI: {ROI:.3g} ' + CS.dim('1/yr') + '\n'
            out += f' PBP: {PBP:.3g} ' + CS.dim('yr') 
        else:
            out += f' NPV: ' + CS.dim('None') + '\n'
            out += f' ROI: ' + CS.dim('None') + '\n'
            out += f' PBP: ' + CS.dim('None')
        return out
    
    def show(self):
        """Print information on TEA."""
        print(self._info())
        
        
        
    