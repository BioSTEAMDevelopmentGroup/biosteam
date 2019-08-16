# -*- coding: utf-8 -*-
"""
Created on Sat Aug 18 14:25:34 2018

@author: yoelr
"""
from ._exceptions import DimensionError
from ._utils import DisplayUnits
from ._species import Species
from ._stream import Stream, mol_flow_dim, mass_flow_dim, vol_flow_dim
from . import _Q


__all__ = ('HeatUtility',)

# Costs from Table 17.1 in Warren D.  Seider et al.-Product and Process Design Principles_ Synthesis, Analysis and Evaluation-Wiley (2016)
# ^This table was made using data from Busche, 1995
# Entry temperature conditions of coolants taken from Table 12.1 in Warren, 2016

# %% Default Heat Transfer species

_Water = Species('Water')


# %% Utility classes


class HeatUtility:
    """Create an HeatUtility object that can choose a utility stream and calculate utility requirements. It can calculate required flow rate, temperature, or phase change of utility. Calculations assume counter current flow rate.
    
    **Parameters**
    
        **efficiency:** [float] Fraction of heat transfered after accounting for heat loss.
    
    **Class Attributes**

        **cooling_agents:** dict[str: UtilityAgent] All cooling utilities available
        
        **heating_agents:** dict[str: UtilityAgent] All heating utilities available
    
    **Examples**
    
        Create a heat utility:
            
        .. code-block:: python
        
           >>> hu = HeatUtility()
           >>> hu
           HeatUtility: None
            duty: 0
            flow: 0
            cost: 0
        
        Calculate utility requirement by calling it with a duty (kJ/hr), and enter and exit temperature (K):
            
        .. code-block:: python
        
           >>> hu(1000, 300, 350)
           >>> hu
           HeatUtility: Low pressure steam
            duty: 1.11e+03 kJ/hr
            flow: 0.0284 kmol/hr,
            cost: 0.00674 USD/hr
       
        All results are accessible:
            
        .. code-block:: python
        
           >>> hu.ID, hu.duty, hu.flow, hu.cost
           ('Low pressure steam',
            1111.111111111111,
            0.028351477551759364,
            0.006741981361808377)
           
    """
    __slots__ = ('_fresh', '_used', 'ID', 'duty',
                 'flow', 'cost', 'efficiency',
                 '_args')
    dT = 5  #: [float] Pinch temperature difference
    
    #: Units of measure for results dictionary
    _units = dict(duty='kJ/hr', flow='kmol/hr', cost='USD/hr')
    
    #: [DisplayUnits] Units of measure for IPython display
    display_units = DisplayUnits(**_units)

    class UtilityAgent:
        """Create a UtilityAgent object that defines a utility option.
        
            **Parameters**
        
                species: [Species] 
                
                molfrac: [tuple] Molar fraction of species
                
                phase: [str] {'g', 'l'}
                
                Hvap: [float] Latent heat of vaporization (kJ/kmol)
                
                T_limit: [float] Temperature limit (K)
                
                price_kJ: [float] Price (USD/kJ)
                
                price_kmol: [float] Price (USD/kmol)
                
                efficiency: [float] Heat transfer efficiency (between 0 to 1)
                
        """
        __slots__ = ('species', 'molfrac', 'T', 'P', 'phase', 'Hvap', 'T_limit',
                     'price_kJ', 'price_kmol', 'efficiency')
        def __init__(self, species, molfrac, T, P, phase, Hvap, T_limit,
                     price_kJ, price_kmol, efficiency):
            self.species = species
            self.molfrac = molfrac
            self.T = T
            self.P = P
            self.phase = phase
            self.Hvap = Hvap
            self.T_limit = T_limit
            self.price_kJ = price_kJ
            self.price_kmol = price_kmol
            self.efficiency = efficiency
        
        def __repr__(self):
            return f"<{type(self).__name__}: T={self.T} K>"
        
        def _ipython_display_(self):
            print(f"{type(self).__name__}:\n"
                 +f" species     .         {repr(self.species)}\n"
                 +f" molfrac     .         ({', '.join([format(i,'.2f') for i in self.molfrac])},)\n"
                 +f" T           K         {self.T:,.2f}\n"
                 +f" P           Pa        {self.P:,.0f}\n"
                 +f" phase       .         '{self.phase}'\n"
                 +f" Hvap        kJ/kmol   {self.Hvap and format(self.Hvap, ',.4g')}\n"
                 +f" T_limit     K         {self.T_limit and format(self.T_limit, ',.2f')}\n"
                 +f" price_kJ    USD/kJ    {self.price_kJ:.4g}\n"
                 +f" price_kmol  USD/kmol  {self.price_kmol:.4g}\n"
                 +f" efficiency  .         {self.efficiency:.2f}")
            
        show = _ipython_display_
            

    cooling_agents = {
        'Cooling water': UtilityAgent(_Water, (1,), 305.372, 101325, 'l',
                                      None, 324.817, 0, 4.8785e-4, 1),
        'Chilled water': UtilityAgent(_Water, (1,), 280.372, 101325, 'l',
                                      None, 300.372, -5e-6, 0, 1),
        'Chilled Brine': UtilityAgent(_Water, (1,), 255.372, 101325, 'l',
                                      None, 275.372, -8.145e-6, 0, 1)}

    heating_agents = {
        'Low pressure steam': UtilityAgent(_Water, (1,), 411.494, 344738.0, 'g',
                                           3.89e+04, None, 0, 0.2378, 0.95),
        'Medium pressure steam': UtilityAgent(_Water, (1,), 454.484, 1034214.0,
                                              'g', 3.63e+04, None, 0, 0.2756, 0.90),
        'High pressure steam': UtilityAgent(_Water, (1,), 508.858, 3102642.0, 'g',
                                            3.21e+04, None, 0, 0.3171, 0.85)}

    def __init__(self, efficiency=None):
        self.ID = ''
        self.cost = self.flow = self.duty = 0
        self.efficiency = efficiency
        
        #: tuple[bool, float] Cached arguments:
        #:     * [0]: [bool] True if duty is negative 
        #:     * [1]: [float] Temperature of operation
        self._args = None

    def _init_streams(self, flow, species, T, P, phase):
        """Initialize utility streams."""
        self._fresh = Stream(None, flow, species, T=T, P=P, phase=phase)
        self._used = Stream.proxy(None, self._fresh)

    def __call__(self, duty, T_in, T_out=None):
        """Calculate utility requirements given the essential parameters.
        
        **Parameters**
        
            **duty:** [float] Unit duty requirement (kJ/hr)
            
            **T_in:** [float] Entering process stream temperature (K)
            
            **T_out:** [float] Exit process stream temperature (K)
        
        """
        # Set pinch and operating temperature
        if T_out is None:
            T_pinch = T_op = T_in
        else:
            T_pinch, T_op = self._get_pinch(duty, T_in, T_out)
        
        if duty == 0:
            self.ID = ''
            self.flow = self.duty = self.cost = 0
            self._args = None
            return
        
        negduty = duty < 0
        args = (negduty, T_op)
        if self._args == args:
            if negduty:
                agent = self.cooling_agents[self.ID]
            else:
                agent = self.heating_agents[self.ID]
        else:
            self._args = args
            # Select heat transfer agent
            if negduty:
                agent = self._select_cooling_agent(T_op)
            else:
                agent = self._select_heating_agent(T_op)
            
        # Calculate utility flow rate requirement
        efficiency = self.efficiency if self.efficiency else agent.efficiency
        duty = duty/efficiency
        if agent.Hvap:
            self._update_flow_wt_phase_change(duty, agent.Hvap)
        else:
            self._update_flow_wt_pinch_T(duty, T_pinch, agent.T_limit, negduty)
        
        # Update and return results
        self.flow = mol = self._fresh.molnet
        self.duty = duty
        self.cost = agent.price_kJ*duty + agent.price_kmol*mol

    @staticmethod
    def _get_pinch(duty, T_in, T_out):
        """Return pinch temperature and operating temperature."""
        if duty < 0:
            return (T_in, T_out) if T_in > T_out else (T_out, T_in)
        else:
            return (T_in, T_out) if T_in < T_out else (T_out, T_in)
    
    # Selection of a heat transfer agent
    def _select_cooling_agent(self, T_pinch):
        """Select a cooling agent that works at the pinch temperature and return relevant information.
        
        **Parameters**

             **T_pinch:**  [float] Pinch temperature of process stream (K)
        
        **Returns:** [UtilityAgent]
        
        """
        dt = 2*self.dT
        T_max = T_pinch - dt
        for ID, agent in self.cooling_agents.items():
            if T_max > agent.T:
                if self.ID != ID:
                    self._init_streams(agent.molfrac, agent.species,
                                       agent.T, agent.P, agent.phase)
                    self.ID = ID
                return agent
        raise RuntimeError(f'no cooling agent that can cool under {T_pinch} K')
            
    def _select_heating_agent(self, T_pinch):
        """Return a heating agent that works at the pinch temperature and return relevant information.
        
        **Parameters**

             **T_pinch:**  [float] Pinch temperature of process stream (K)
        
        **Returns:** [UtilityAgent]
        
        """
        dt = 2*self.dT
        T_min = T_pinch + dt
        for ID, agent in self.heating_agents.items():
            if T_min < agent.T:
                if self.ID != ID:
                    self._init_streams(agent.molfrac, agent.species,
                                       agent.T, agent.P, agent.phase)
                    self.ID = ID
                return agent
        raise RuntimeError(f'no heating agent that can heat over {T_pinch} K')

    # Main Calculations
    def _update_flow_wt_pinch_T(self, duty, T_pinch, T_limit, negduty):
        """Set utility Temperature at the pinch, calculate and set minimum net flowrate of the utility to satisfy duty and update."""
        self._used.T = self._T_exit(T_pinch, self.dT, T_limit, negduty)
        self._update_utility_flow(self._fresh, self._used, duty)

    def _update_flow_wt_phase_change(self, duty, latent_heat):
        """Change phase of utility, calculate and set minimum net flowrate of the utility to satisfy duty and update."""
        f = self._fresh
        u = self._used
        u._phase = 'g' if f._phase=='l' else 'l'
        u._mol[:] = duty/latent_heat

    # Subcalculations
    @staticmethod
    def _update_utility_flow(fresh, utility, duty):
        """Changes flow rate of utility such that it can meet the duty requirement"""
        utility._mol *= duty/(fresh.H - utility.H)

    @staticmethod
    def _T_exit(T_pinch, dT, T_limit, negduty):
        """Return exit temperature of the utility in a counter current heat exchanger

        **Parameters**

             **T_pinch:** [float] Pinch temperature of process stream (K)

             **dT:** [float] Pinch temperature difference (K)

             **negduty:** [bool] True if exit temperature should be lower (process stream is gaining energy)

        """
        if negduty:
            T_exit = T_pinch - dT
            if T_limit and T_limit < T_exit: T_exit = T_limit
        else:
            T_exit = T_pinch + dT
            if T_limit and T_limit > T_exit: T_exit = T_limit
        return T_exit

    def _info_data(self, duty, flow, cost):
        # Get units of measure
        su = self.display_units
        units = self._units
        duty_units = duty or su.duty
        flow_units = flow or su.flow
        cost_units = cost or su.cost
        
        # Select flow dimensionality
        flow_dim = _Q(0, flow_units).dimensionality
        if flow_dim == mol_flow_dim:
            flowattr = 'molnet'
        elif flow_dim == mass_flow_dim:
            flowattr = 'massnet'
        elif flow_dim == vol_flow_dim:
            flowattr = 'volnet'
        else:
            raise DimensionError(f"dimensions for flow units must be in molar, mass or volumetric flow rates, not '{flow_dim}'")
        
        # Change units and return info string
        try:
            u_in = self._fresh
            flownet = getattr(u_in, flowattr)
            flowunits = u_in.units[flowattr]
            flow = _Q(flownet, flowunits).to(flow_units).magnitude
        except:
            flow = _Q(self.flow, 'kmol/hr').to(flow_units).magnitude
        
        duty = _Q(self.duty, units['duty']).to(duty_units).magnitude
        cost = _Q(self.cost, units['cost']).to(cost_units).magnitude
        return duty, flow, cost, duty_units, flow_units, cost_units
        
    def __repr__(self):
        if self.ID:
            duty, flow, cost, duty_units, flow_units, cost_units = self._info_data(None, None, None)
            return f'<{self.ID}: {self.duty:.3g} {duty_units}, {self.flow:.3g} {flow_units}, {self.cost:.3g} {cost_units}>'
        else:
            return f'<{type(self).__name__}: None>'
        
    # Representation
    def _info(self, duty, flow, cost):
        """Return string related to specifications"""
        if not self.ID:
            return (f'{type(self).__name__}: None\n'
                    +f' duty: 0\n'
                    +f' flow: 0\n'
                    +f' cost: 0')
        else:
            (duty, flow, cost, duty_units,
             flow_units, cost_units) = self._info_data(duty, flow, cost)
            return (f'{type(self).__name__}: {self.ID}\n'
                    +f' duty:{duty: .3g} {duty_units}\n'
                    +f' flow:{flow: .3g} {flow_units}\n'
                    +f' cost:{cost: .3g} {cost_units}')
            

    def show(self, duty=None, flow=None, cost=None):
        """Print all specifications"""
        print(self._info(duty, flow, cost))
    _ipython_display_ = show


del _Water