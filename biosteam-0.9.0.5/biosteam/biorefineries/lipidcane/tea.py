# -*- coding: utf-8 -*-
"""
Created on Thu Aug  1 11:48:19 2019

@author: yoelr
"""
from biosteam import TEA

__all__ = ('LipidcaneTEA',)

class LipidcaneTEA(TEA):
    """Create a LipidcaneTEA object for techno-economic analysis of a biorefinery [1]
    
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
    
            **labor_cost:** [float] Total labor cost (USD/yr).

            **fringe_benefits:** [float] Cost of fringe benefits as a fraction of labor cost.
            
            **property_tax:** [float] Fee as a fraction of fixed capital investment.

            **property_insurance:** [float] Fee as a fraction of fixed capital investment.    
    
            **supplies:** [float] Yearly fee as a fraction of labor cost.

            **maintenance:** [float] Yearly fee as a fraction of fixed capital investment.

            **administration:** [float] Yearly fee as a fraction of fixed capital investment.
        
    **References**
    
        [1] Huang, H., Long, S., & Singh, V. (2016). Techno-economic analysis of biodiesel and ethanol co-production from lipid-producing sugarcane. Biofuels, Bioproducts and Biorefining, 10(3), 299–315. https://doi.org/10.1002/bbb.1640
    
    """
    __slots__ = ('labor_cost', 'fringe_benefits', 'maintenance',
                 'property_tax', 'property_insurance', '_FCI_cached',
                 'supplies', 'maintanance', 'administration')
    
    def __init__(self, system, IRR, duration, depreciation, income_tax,
                 operating_days, lang_factor, construction_schedule, WC_over_FCI,
                 labor_cost, fringe_benefits, property_tax,
                 property_insurance, supplies, maintenance, administration):
        super().__init__(system, IRR, duration, depreciation, income_tax,
                         operating_days, lang_factor, construction_schedule,
                         startup_months=0, startup_FOCfrac=0, startup_VOCfrac=0,
                         startup_salesfrac=0, finance_interest=0, finance_years=0, 
                         finance_fraction=0, WC_over_FCI=WC_over_FCI)
        self.labor_cost = labor_cost
        self.fringe_benefits = fringe_benefits
        self.property_tax = property_tax
        self.property_insurance = property_insurance
        self.supplies= supplies
        self.maintenance = maintenance
        self.administration = administration
        
    def _TDC(self, DPI):
        return DPI
    
    def _FCI(self, TDC):
        self._FCI_cached = TDC
        return TDC
    
    def _FOC(self, FCI):
        return (FCI*(self.property_tax + self.property_insurance
                     + self.maintenance + self.administration)
                + self.labor_cost*(1+self.fringe_benefits+self.supplies))