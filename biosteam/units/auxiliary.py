# -*- coding: utf-8 -*-
"""
This module contains functions for adding auxliary unit operations.
"""
import biosteam as bst

__all__ = ('Auxiliary',)
    

class Auxiliary:
    """Abstract class for light-weight auxiliary unit. The class should 
    compute all results during initialization."""
    __slots__ = (
        'owner', 
        'auxname',
        'auxiliary_units',
        'power_utility',
        'heat_utilities', 
        'design_results',
        'baseline_purchase_costs',
        'purchase_costs',
        'installed_costs',
        'F_M', 'F_D', 'F_P', 'F_BM',
        'auxiliary_unit_names',
        'parallel',
    )
    add_power_utility = bst.Unit.add_power_utility
    add_heat_utility = bst.Unit.add_heat_utility 
    create_heat_utility = bst.Unit.create_heat_utility 
    get_design_result = bst.Unit.get_design_result
    set_design_result = bst.Unit.set_design_result
    convert_design_result = bst.Unit.convert_design_result
    
    def __init_subclass__(cls):
        if '_units' not in cls.__dict__: cls._units = {}
    
    def __init__(self):
        self.power_utility = bst.PowerUtility()
        self.heat_utilities = []
        self.baseline_purchase_costs = {}
        self.purchase_costs = {}
        self.installed_costs = {}
        self.design_results = {}
        self.F_M = {}
        self.F_D = {} 
        self.F_P = {}
        self.F_BM = {}
        self.parallel = {}
        self.auxiliary_unit_names = ()
        
    def _design(self): pass
    def _cost(self): pass
        
    def _setup(self):
        results = (self.baseline_purchase_costs, self.purchase_costs, 
                   self.installed_costs, self.F_M, self.F_D, self.F_P,
                   self.F_BM, self.design_results)
        for i in results: i.clear()
        for i in self.heat_utilities: i.empty()
        self.heat_utilities.clear()
        self.power_utility.empty()
        
    def _load_costs(self):
        r"""
        Calculate and save free on board (f.o.b.) purchase costs and
        installed equipment costs (i.e. bare-module cost) for each item in the 
        :attr:`~Auxiliary.baseline_purchase_costs` dictionary.
        
        Notes
        -----
        As explained in [1]_, the f.o.b. purchase cost is given by:
        
        .. math::
           
           C_{P} = C_{Pb}F_{D}F_{P}F_{M}
        
        And the installed equipment cost is given by:
        
        .. math::
           
           C_{BM} = C_{Pb} (F_{BM} + F_{D}F_{P}F_{M} - 1)
        
        Where:
            * :math:`C_{Pb}`: Baseline purchase cost.
            * :math:`F_{BM}`: Bare module factor.
            * :math:`F_{D}`: Design factor.
            * :math:`F_{P}`: Pressure factor.
            * :math:`F_{M}`: Material factor.
        
        Values for the bare-module, design, pressure, and material factors of 
        each equipment should be stored in the :attr:`~Auxiliary.F_BM`, :attr:`~Auxiliary.F_D`, 
        :attr:`~Auxiliary.F_P`, and :attr:`~Auxiliary.F_M` dictionaries.
        
        Warning
        -------
        If an item is listed in the :attr:`~Auxiliary.purchase_costs` dictionary but not in the
        :attr:`~Auxiliary.baseline_purchase_costs` dictionary, the baseline purchase cost is 
        assumed to be the same as the purchase cost.
        
        References
        ----------
        .. [1] Seider, W. D., Lewin,  D. R., Seader, J. D., Widagdo, S., Gani, R., & Ng, M. K. (2017). Product and Process Design Principles. Wiley. Cost Accounting and Capital Cost Estimation (Chapter 16)
        
        """
        F_BM = self.F_BM
        F_D = self.F_D
        F_P = self.F_P
        F_M = self.F_M
        baseline_purchase_costs = self.baseline_purchase_costs
        purchase_costs = self.purchase_costs
        installed_costs = self.installed_costs
        
        # Load main costs
        for i in purchase_costs:
            if i not in baseline_purchase_costs:
                baseline_purchase_costs[i] = purchase_costs[i]
        for name, Cpb in baseline_purchase_costs.items(): 
            if name in installed_costs and name in purchase_costs:
                continue # Assume costs already added elsewhere using another method
            F = F_D.get(name, 1.) * F_P.get(name, 1.) * F_M.get(name, 1.)
            try:
                installed_costs[name] = Cpb * (F_BM[name] + F - 1.)
            except KeyError:
                F_BM[name] = 1.
                installed_costs[name] = purchase_costs[name] = Cpb * F
            else:
                purchase_costs[name] = Cpb * F
