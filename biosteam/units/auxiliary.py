# -*- coding: utf-8 -*-
"""
This module contains functions for adding auxliary unit operations.
"""
from . import design_tools as design
import biosteam as bst

__all__ = ('Auxiliary', 'VacuumSystem', 'Agitator')
    

class Auxiliary:
    """Abstract class for light-weight auxiliary unit. The class should 
    compute all results during initialization."""
    __slots__ = (
        'power_utility',
        'heat_utilities', 
        'baseline_purchase_costs',
        'purchase_costs',
        'installed_costs',
        'F_M', 'F_D', 'F_P', 'F_BM',
        'owner', 
    )
    
    add_power_utility = bst.Unit.add_power_utility
    add_heat_utility = bst.Unit.add_heat_utility 
    create_heat_utility = bst.Unit.create_heat_utility 
    
    def __init__(self):
        self.power_utility = bst.PowerUtility()
        self.heat_utilities = []
        self.baseline_purchase_costs = {}
        self.purchase_costs = {}
        self.installed_costs = {}
        self.F_M = {}
        self.F_D = {} 
        self.F_P = {}
        self.F_BM = {}
        
    def _setup(self):
        results = (self.baseline_purchase_costs, self.purchase_costs, 
                   self.installed_costs, self.F_M, self.F_D, self.F_P,
                   self.F_BM)
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

        
class Agitator(Auxiliary):
    __slots__ = ()
    
    def __init__(self, kW):
        super().__init__()
        self.add_power_utility(kW) # kW
        hp = kW * 1.34102
        self.baseline_purchase_costs['Agitator'] = design.compute_closed_vessel_turbine_purchase_cost(hp)


class VacuumSystem(Auxiliary): 
    __slots__ = ()
    
    def __init__(self, unit=None, vacuum_system_preference=None, 
                 F_mass=None, F_vol=None, P_suction=None, vessel_volume=None):
        super().__init__()
        if unit:
            # Deduce arguments if unit given
            if F_mass is None or F_vol is None:
                # If vapor is condensed, assume vacuum system is after condenser
                vents = [i for i in unit.outs if 'g' in i.phase]
                F_mass = 0.
                F_vol = 0.
                for vapor in vents:
                    hx = vapor.sink
                    if isinstance(hx, bst.HX):
                        index = hx.ins.index(vapor)
                        vapor = hx.outs[index]
                        if 'g' not in vapor.phase: continue
                    if isinstance(vapor, bst.MultiStream):
                        F_mass += vapor['g'].F_mass
                        F_vol += vapor['g'].F_vol
                    else:
                        F_mass += vapor.F_mass
                        F_vol += vapor.F_vol
            if P_suction is None:
                P_suction = getattr(unit, 'P', None) or getattr(unit, 'P_suction', None)
                if P_suction is None: P_suction = min([i.P for i in unit.outs])
            if vessel_volume is None:
                if 'Total volume' in unit.design_results:
                    vessel_volume = unit.get_design_result('Total volume', 'm3')
                else:
                    raise ValueError("'Total volume' was not found in design results; "
                                     "'vesse_volume' parameter could not be deduced")
        else:
            # In case user does not supply all arguments
            params = [('F_mass', F_mass), 
                      ('F_vol', F_vol),
                      ('P_suction', P_suction), 
                      ('vessel_volume', vessel_volume)]
            for name, value in params:
                if value is None: 
                    raise ValueError(
                        f"missing argument '{name}'; cannot deduce '{name}'"
                         "when no unit is specified"
                    )
        capex = self.baseline_purchase_costs # Assume all costs the same
        vacuum_results = design.compute_vacuum_system_power_and_cost(
            F_mass, F_vol, P_suction, vessel_volume, vacuum_system_preference
        )
        name = vacuum_results['Name']
        capex[name] = vacuum_results['Cost']
        heating_agent = vacuum_results['Heating agent']
        if heating_agent: # Steam turbine
            vacuum_steam = self.create_heat_utility()
            vacuum_steam.set_utility_by_flow_rate(heating_agent, vacuum_results['Steam flow rate'])
            if vacuum_results['Condenser']: 
                self.add_heat_utility(-vacuum_steam.unit_duty, 373.15) # Vacuum cooling water
        self.add_power_utility(vacuum_results['Work'])