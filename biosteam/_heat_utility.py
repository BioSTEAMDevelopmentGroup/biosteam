# -*- coding: utf-8 -*-
"""
Created on Sat Aug 18 14:25:34 2018

@author: yoelr
"""
from thermosteam._stream import Stream, mol_units, mass_units, vol_units
from thermosteam.base.units_of_measure import convert, get_dimensionality, DisplayUnits
from thermosteam.exceptions import DimensionError
from thermosteam import Thermo

__all__ = ('HeatUtility',)

# Costs from Table 17.1 in Warren D.  Seider et al.-Product and Process Design Principles_ Synthesis, Analysis and Evaluation-Wiley (2016)
# ^This table was made using data from Busche, 1995
# Entry temperature conditions of coolants taken from Table 12.1 in Warren, 2016

# %% Default Heat Transfer chemicals

_Water = Thermo(['Water'])

# %% Utility classes

class UtilityAgent:
    """Create a UtilityAgent object that defines a utility option.
    
    Parameters
    ----------
    thermo: Thermo
    z_mol: tuple[float]
            Molar fraction of chemicals.
    phase : {'g', 'l'}
            Either gas ("g") or liquid ("l").
    Hvap : float
           Latent heat of vaporization (kJ/kmol)
    T_limit : float
              Temperature limit (K)
    price_kJ : float
               Price (USD/kJ)
    price_kmol : float
                 Price (USD/kmol)
    efficiency : float
                 Heat transfer efficiency (between 0 to 1)
            
    """
    __slots__ = ('thermo', 'z_mol', 'T', 'P', 'phase', 'Hvap', 'T_limit',
                 'price_kJ', 'price_kmol', 'efficiency')
    def __init__(self, thermo, z_mol, T, P, phase, Hvap, T_limit,
                 price_kJ, price_kmol, efficiency):
        self.thermo = thermo
        self.z_mol = z_mol
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
             +f" thermo                Thermo({self.thermo.chemicals}, ...)\n"
             +f" z_mol                 ({', '.join([format(i,'.2f') for i in self.z_mol])},)\n"
             +f" T           K         {self.T:,.2f}\n"
             +f" P           Pa        {self.P:,.0f}\n"
             +f" phase                 '{self.phase}'\n"
             +f" Hvap        kJ/kmol   {self.Hvap and format(self.Hvap, ',.4g')}\n"
             +f" T_limit     K         {self.T_limit and format(self.T_limit, ',.2f')}\n"
             +f" price_kJ    USD/kJ    {self.price_kJ:.4g}\n"
             +f" price_kmol  USD/kmol  {self.price_kmol:.4g}\n"
             +f" efficiency            {self.efficiency:.2f}")
        
    show = _ipython_display_

            
class UtilityAgents:
    __slots__ = ('_agents',)
    
    def __init__(self, agents):
        self._agents = dict(agents)
        
    def __getitem__(self, key):
        return self._agents[key]
    
    def __setitem__(self, ID, agent):
        assert isinstance(agent, UtilityAgent)
        agents = sorted([*self._agents.items(), (ID, agent)], key=self._sorting_key)
        self._agents = dict(agents)
    
    def __iter__(self):
        yield from self._agents.items()
        
    def __repr__(self):
        return f"{type(self).__name__}({self._agents})"


class HeatingAgents(UtilityAgents):
    __slots__ = ()
    
    @staticmethod
    def _sorting_key(ID_agent):
        _, agent = ID_agent
        return agent.T
    
    def select_agent(self, T_pinch):
        """Return a heating agent that works at the pinch temperature and return UtilityAgent.
        
        Parameters
        ----------
        T_pinch : float
                  Pinch temperature of process stream (K)
        
        """
        for ID, agent in self:
            if T_pinch < agent.T: return ID, agent
        raise RuntimeError(f'no heating agent that can heat over {T_pinch} K')
    

class CoolingAgents(UtilityAgents):
    __slots__ = ()
    
    @staticmethod
    def _sorting_key(ID_agent):
        _, agent = ID_agent
        return 1 / agent.T

    # Selection of a heat transfer agent
    def select_agent(self, T_pinch):
        """Select a cooling agent that works at the pinch temperature and return UtilityAgent.
        
        Parameters
        ----------
        T_pinch : float
                  Pinch temperature of process stream (K)
        
        """
        for ID, agent in self:
            if T_pinch > agent.T: return ID, agent
        raise RuntimeError(f'no cooling agent that can cool under {T_pinch} K')


class HeatUtility:
    """Create an HeatUtility object that can choose a utility stream and calculate utility requirements. It can calculate required flow rate, temperature change, or phase change of utility. Calculations assume counter current flow rate.
    
    Parameters
    ----------
    efficiency=None : float, optional
        Enforced fraction of heat transfered from utility (due
        to losses to environment).
    
    Attributes
    ----------
    cooling_agents : dict[str: UtilityAgent]
                     All cooling utilities available
    heating_agents : dict[str: UtilityAgent]
                     All heating utilities available
    
    Examples
    --------
    Create a heat utility:
        
    .. code-block:: python
    
       >>> from biosteam import HeatUtility
       >>> hu = HeatUtility()
       >>> hu.show()
       HeatUtility: None
        duty: 0
        flow: 0
        cost: 0
    
    Calculate utility requirement by calling it with a duty (kJ/hr), and entrance and exit temperature (K):
        
    .. code-block:: python
    
       >>> hu(1000, 300, 350)
       >>> hu.show()
       HeatUtility: Low pressure steam
        duty: 1.05e+03 kJ/hr
        flow: 0.0271 kmol/hr
        cost: 0.00644 USD/hr
   
    All results are accessible:
        
    .. code-block:: python
    
       >>> hu.ID, hu.duty, hu.flow, hu.cost
       ('Low pressure steam', 1052.6315789473686, 0.027090580063500323, 0.006442139939100377)
           
    """
    __slots__ = ('inlet_utility_stream', 'outlet_utility_stream', 'ID', 'duty',
                 'flow', 'cost', 'efficiency', '_args')
    dT = 5  #: [float] Pinch temperature difference
    
    #: [DisplayUnits] Units of measure for IPython display
    display_units = DisplayUnits(duty='kJ/hr', flow='kmol/hr', cost='USD/hr')

    cooling_agents = CoolingAgents({
        'Cooling water': UtilityAgent(_Water, (1,), 305.372, 101325, 'l',
                                      None, 324.817, 0, 4.8785e-4, 1),
        'Chilled water': UtilityAgent(_Water, (1,), 280.372, 101325, 'l',
                                      None, 300.372, -5e-6, 0, 1),
        'Chilled Brine': UtilityAgent(_Water, (1,), 255.372, 101325, 'l',
                                      None, 275.372, -8.145e-6, 0, 1)})

    heating_agents = HeatingAgents({
        'Low pressure steam': UtilityAgent(_Water, (1,), 411.494, 344738.0, 'g',
                                           38856., None, 0, 0.2378, 0.95),
        'Medium pressure steam': UtilityAgent(_Water, (1,), 454.484, 1034214.0,
                                              'g', 36264., None, 0, 0.2756, 0.90),
        'High pressure steam': UtilityAgent(_Water, (1,), 508.858, 3102642.0, 'g',
                                            32124., None, 0, 0.3171, 0.85)})

    def __init__(self, efficiency=None):
        self.efficiency = efficiency
        self.empty()

    def _init_streams(self, flow, thermo, T, P, phase):
        """Initialize utility streams."""
        self.inlet_utility_stream = Stream(None, flow, thermo=thermo, T=T, P=P, phase=phase)
        self.outlet_utility_stream = self.inlet_utility_stream.flow_proxy()

    def empty(self):
        self.ID = ''
        self.cost = self.flow = self.duty = 0
        self._args = None
        #: tuple[bool, float] Cached arguments:
        #:     * [0]: [bool] True if duty is negative 
        #:     * [1]: [float] Pinch temperature of fresh utility
            

    def __call__(self, duty, T_in, T_out=None):
        """Calculate utility requirements given the essential parameters.
        
        Parameters
        ----------
        duty : float
               Unit duty requirement (kJ/hr)
        T_in : float
               Entering process stream temperature (K)
        T_out : float, optional
                Exit process stream temperature (K)
        
        """
        if duty == 0:
            self.empty()
            return
        T_out = T_out or T_in
        iscooling = duty < 0
        T_pinch_in, T_pinch_out = self._get_pinch(iscooling, T_in, T_out, self.dT)

        ## Select heat transfer agent ##
        
        args = (iscooling, T_pinch_in)
        agents = self.cooling_agents if iscooling else self.heating_agents
        if self._args == args:            
            # Use cached agent
            agent = agents[self.ID]
        else:
            # Select agent
            self._args = args
            ID, agent = agents.select_agent(T_pinch_in)
            self.load_agent(ID, agent)
            
        ## Calculate utility requirement ##
        
        efficiency = self.efficiency or agent.efficiency
        duty = duty/efficiency
        if agent.Hvap:
            # Phase change
            self.outlet_utility_stream.phase = 'g' if self.inlet_utility_stream.phase == 'l' else 'l'
            dH = agent.Hvap * self.outlet_utility_stream.F_mol
        else:
            # Temperature change
            self.outlet_utility_stream.T = self._T_exit(T_pinch_out, agent.T_limit, iscooling)
            dH = self.inlet_utility_stream.H - self.outlet_utility_stream.H
        
        # Update utility flow
        self.outlet_utility_stream.mol[:] *= duty / dH
        
        # Update and return results
        self.flow = F_mol = self.inlet_utility_stream.F_mol
        self.duty = duty
        self.cost = agent.price_kJ*duty + agent.price_kmol*F_mol

    @staticmethod
    def _get_pinch(iscooling, T_in, T_out, dT):
        """Return pinch entrance and exit temperature of utility."""
        if iscooling:
            assert T_in + 1e-6 >= T_out, "entrance temperature must be larger than exit temperature if cooling"
            T_pinch_in = T_out - dT
            T_pinch_out = T_in - dT
        else:
            assert T_in <= T_out + 1e-6, "entrance temperature must be smaller than exit temperature if heating"
            T_pinch_in = T_out + dT
            T_pinch_out = T_in + dT
        return T_pinch_in, T_pinch_out
    
    def load_agent(self, ID, agent):
        if self.ID != ID:
            self._init_streams(agent.z_mol, agent.thermo,
                               agent.T, agent.P, agent.phase)
            self.ID = ID

    # Subcalculations

    @staticmethod
    def _T_exit(T_pinch, T_limit, iscooling):
        """Return exit temperature of the utility in a counter current heat exchanger

        Parameters
        ----------
        T_pinch : float
                  Pinch temperature of utility stream (K)
        iscooling : bool
                  True if exit temperature should be lower (process stream is gaining energy)

        """
        if iscooling:
            return T_limit if T_limit and T_limit < T_pinch else T_pinch
        else:
            return T_limit if T_limit and T_limit > T_pinch else T_pinch
        

    def _info_data(self, duty, flow, cost):
        # Get units of measure
        su = self.display_units
        duty_units = duty or su.duty
        flow_units = flow or su.flow
        cost_units = cost or su.cost
        
        # Change units and return info string
        flow = self.inlet_utility_stream.get_total_flow(flow_units)
        duty = convert(self.duty, 'kJ/hr', duty_units)
        cost = convert(self.cost, 'USD/hr', cost_units)
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