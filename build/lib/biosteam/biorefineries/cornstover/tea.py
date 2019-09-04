# -*- coding: utf-8 -*-
"""
Created on Thu Aug  1 21:48:12 2019

@author: yoelr
"""

from biosteam import TEA

__slots__ = ('CornstoverTEA',)

class CornstoverTEA(TEA):
    
    __slots__ = ('OSBL_units', 'warehouse', 'site_development',
                 'additional_piping', 'proratable_costs', 'field_expenses',
                 'construction', 'contingency', 'other_indirect_costs', 
                 'labor_cost', 'labor_burden', 'property_insurance',
                 'maintenance', '_ISBL_DPI_cached', '_FCI_cached',
                 '_utility_cost_cached')
    
    def __init__(self, system, IRR, duration, depreciation, income_tax,
                 operating_days, lang_factor, construction_schedule,
                 startup_months, startup_FOCfrac, startup_VOCfrac,
                 startup_salesfrac, WC_over_FCI,  finance_interest,
                 finance_years, finance_fraction, OSBL_units, warehouse,
                 site_development, additional_piping, proratable_costs,
                 field_expenses, construction, contingency,
                 other_indirect_costs, labor_cost, labor_burden,
                 property_insurance, maintenance):
        super().__init__(system, IRR, duration, depreciation, income_tax,
                         operating_days, lang_factor, construction_schedule,
                         startup_months, startup_FOCfrac, startup_VOCfrac,
                         startup_salesfrac, WC_over_FCI,  finance_interest,
                         finance_years, finance_fraction)
        self.OSBL_units = OSBL_units
        self.warehouse = warehouse
        self.site_development = site_development
        self.additional_piping = additional_piping
        self.proratable_costs = proratable_costs
        self.field_expenses = field_expenses
        self.construction = construction
        self.contingency = contingency
        self.other_indirect_costs = other_indirect_costs
        self.labor_cost = labor_cost
        self.labor_burden = labor_burden
        self.property_insurance = property_insurance
        self.maintenance = maintenance
    
    @property
    def utility_cost(self):
        self._utility_cost_cached = utility_cost = super().utility_cost
        return utility_cost
    
    def _ISBL_DPI(self, DPI):
        """Direct permanent investment of units inside battery limits."""
        if self.lang_factor:
            self._ISBL_DPI_cached = DPI - sum([i.purchase_cost for i in self.OSBL_units]) * self.lang_factor
        else:
            self._ISBL_DPI_cached = DPI - sum([i.installation_cost for i in self.OSBL_units])
        return self._ISBL_DPI_cached
        
    def _TDC(self, DPI):
        return DPI + self._ISBL_DPI(DPI)*(self.warehouse
                                          + self.site_development
                                          + self.additional_piping)
    
    def _indirect_costs(self, TDC):
        return TDC*(self.proratable_costs + self.field_expenses
                    + self.construction + self.contingency
                    + self.other_indirect_costs)
    
    def _FCI(self, TDC):
        self._FCI_cached = FCI = TDC + self._indirect_costs(TDC)
        return FCI
    
    def _FOC(self, FCI):
        return (FCI * self.property_insurance
                + self._ISBL_DPI_cached * self.maintenance
                + self.labor_cost*(1+self.labor_burden))