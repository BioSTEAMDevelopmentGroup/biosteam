# -*- coding: utf-8 -*-
"""
Created on Sat Aug 18 14:25:34 2018

@author: yoelr
"""
from biosteam.exceptions import DimensionError, biosteamError
from biosteam.species import Species
from biosteam.stream import Stream, mol_flow_dim, mass_flow_dim, vol_flow_dim
from bookkeep import SmartBook, UnitManager
from biosteam import Q_, pd, np


__all__ = ('HeatUtility',)

# Costs from Table 17.1 in Warren D.  Seider et al.-Product and Process Design Principles_ Synthesis, Analysis and Evaluation-Wiley (2016)
# ^This table was made using data from Busche, 1995
# Entry temperature conditions of coolants taken from Table 12.1 in Warren, 2016

# %% Default Heat Transfer species

_Water = Species('Water')

# %% Data for cooling utilities 

_columns = ('Type', 'species', 'molar fraction', 'T (K)',
            'P (Pa)', 'phase', 'T limit (K)',
            'price (USD/kJ)', 'price (USD/kmol)', 'efficiency')

_cooling_index = ('Cooling water',
                  'Chilled water',
                  'Chilled Brine')

_cooling_water = ('sensible',      # Type [str]: {'sensible', 'latent'}
                  _Water,          # species [Species]
                  (1,),            # flow: [tuple]
                  305.372,         # T (K)
                  101325,          # P (Pa)
                  'liq',           # phase: {'g', 'l'}
                  324.817,         # T limit (K)
                  0,               # price (USD/kJ)
                  4.8785e-4,       # price (USD/kmol)
                  1)               # heat transfer efficiency

_chilled_water = ('sensible',
                  _Water,
                  (1,),
                  280.372,
                  101325,
                  'liq',
                  300.372,
                  -5e-6,
                  0,
                  1)

_chilled_brine = ('sensible',
                  _Water,
                  (1,),
                  255.372,
                  101325,
                  'liq',
                  275.372,
                  -8.145e-6,
                  0,
                  1)


# %% Data for heating utilities

_heating_index = ('Low pressure steam',
                  'Medium pressure steam',
                  'High pressure steam')

_low_pressure_steam = ('latent',
                       _Water,
                       (1,),
                       411.494,
                       344738.0,
                       'gas',
                       None,
                       0,
                       0.2378,
                       0.90)

_medium_pressure_steam = ('latent',
                          _Water,
                          (1,),
                          454.484,
                          1034214.0,
                          'gas',
                          None,
                          0,
                          0.2756,
                          0.80)

_high_pressure_steam = ('latent',
                        _Water,
                        (1,),
                        508.858,
                        3102642.0,
                        'gas',
                        None,
                        0,
                        0.3171,
                        0.70)


# %% Utility classes


class HeatUtility:
    """Create an HeatUtility object that can choose a utility stream and calculate utility requirements. It can calculate required flow rate, temperature, or phase change of utility. Calculations assume counter current flow rate.
    
    **Parameters**
        
        **efficiency:** [float] Fraction of heat transfered after accounting for heat loss.
    
    **Class Attributes**
    
        **cooling_agents:** [DataFrame] All cooling utilities available
        
        **heating_agents:** [DataFrame] All heating utilities available
    
    **Examples**
    
        Create a heat utility with a source:
            
        .. code-block:: python
        
           >>> hu = HeatUtility('Flash')
           >>> hu.show()
           HeatUtility: None
            duty: 0
            flow: 0
            cost: 0
        
        Calculate utility requirement by calling it with a duty (kJ/hr) and temperature (K):
            
        .. code-block:: python
        
           >>> hu(1000, 300, 350)
           >>> hu.show()
           HeatUtility: Low pressure steam
            duty: 0.00607 USD/hr
            flow: 0.0255 kmol/hr
            cost: 1e+03 kJ/hr
       
        All results are accessible:
            
        .. code-block:: python
        
           >>> hu.ID, hu.duty, hu.flow, hu.cost
           ('Low pressure steam', 0.00607, 0.0255, 1e+03)
           
    """
    __slots__ = ('_fresh', '_vap', '_liq', 'ID', 'duty', 'flow', 'cost', 'efficiency')
    dT = 5  #: [float] Pinch temperature difference
    
    #: Units of measure for results dictionary
    _units = dict(duty='kJ/hr', flow='kmol/hr', cost='USD/hr')

    # All cooling utilities available
    cooling_agents = pd.DataFrame([_cooling_water,
                                   _chilled_water,
                                   _chilled_brine],
                                  columns=_columns,
                                  index=_cooling_index).transpose()

    # All heating utilities available
    heating_agents = pd.DataFrame([_low_pressure_steam,
                                   _medium_pressure_steam,
                                   _high_pressure_steam],
                                  columns=_columns,
                                  index=_heating_index).transpose()

    def __init__(self, efficiency=None):
        self.ID = ''
        self.duty = 0
        self.flow = 0
        self.cost = 0
        self.efficiency = efficiency

    def _init_streams(self, flow, species, T, P, phase):
        """Initialize utility streams."""
        self._fresh = Stream('*', flow, species, T=T, P=P, phase=phase)
        self._liq = Stream('*', species=species, phase='l')
        self._vap = Stream('*', species=species, phase='g')

    def __call__(self, duty:'kJ/hr', T_in:'K', T_out:'K'=None):
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
        
        
        # Select heat transfer agent
        if duty < 0:
            (ID, Type, price_duty,
             price_mol, T_limit,
             efficiency) = self._select_cooling_agent(T_op)
        elif duty > 0:
            (ID, Type, price_duty,
            price_mol, T_limit,
            efficiency) = self._select_heating_agent(T_op)
        else:
            self.ID = ''
            self.flow = self.duty = self.cost = 0
            return
        
        # Calculate utility flow rate requirement
        efficiency = self.efficiency if self.efficiency else efficiency
        duty = duty/efficiency
        if Type == 'sensible':
            self._update_flow_wt_pinch_T(duty, T_pinch, T_limit)
        elif Type == 'latent':
            self._update_flow_wt_phase_change(duty)
        else:
            raise biosteamError(f"Invalid heat transfer agent '{ID}'.")
        
        # Update and return results
        self.ID = ID
        self.flow = mol = self._fresh.molnet
        self.duty = duty
        self.cost = price_duty*duty + price_mol*mol

    @staticmethod
    def _get_pinch(duty, T_in, T_out) -> '(T_pinch, T_op)':
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
        
        **Returns**
        
            **Type:** [str] {'latent', 'sensible'}
            
            **price_duty:** [float] (USD/kJ)
            
            **price_mass:** [float] (USD/kg)
            
            **T_limit:** [float] Maximum or minimum temperature of agent (K)
            
            **efficiency:** [float] Heat transfer efficiency
        
        """
        dt = 2*self.dT
        T_max = T_pinch - dt
        cooling_agents = self.cooling_agents
        for ID in cooling_agents:
            (Type, species, flow, T, P,
             phase, T_limit, price_duty,
             price_mol, efficiency) = cooling_agents[ID]
            if T_max > T:
                self._init_streams(flow, species, T, P, phase[0])
                return ID, Type, price_duty, price_mol, T_limit, efficiency
        raise biosteamError(f'No cooling agent that can cool under {T_pinch} K')
            
    def _select_heating_agent(self, T_pinch):
        """Select a heating agent that works at the pinch temperature and return relevant information.
        
        **Parameters**

             **T_pinch:**  [float] Pinch temperature of process stream (K)
        
        **Returns**
        
            **Type:** [str] {'latent', 'sensible'}
            
            **price_duty:** [float] (USD/kJ)
            
            **price_mass:** [float] (USD/kg)
            
            **T_limit:** [float] Maximum or minimum temperature of agent (K)
            
            **efficiency:** [float] Heat transfer efficiency
        
        """
        dt = 2*self.dT
        T_min = T_pinch + dt
        heating_agents = self.heating_agents
        for ID in heating_agents:
            (Type, species, flow, T, P,
             phase, T_limit, price_duty,
             price_mol, efficiency) =  heating_agents[ID]
            if T_min < T:
                self._init_streams(flow, species, T, P, phase[0])
                return ID, Type, price_duty, price_mol, T_limit, efficiency
        raise biosteamError(f'No heating agent that can heat over {T_pinch} K')

    # Main Calculations
    def _update_flow_wt_pinch_T(self, duty, T_pinch, T_limit):
        """Set utility Temperature at the pinch, calculate and set minimum net flowrate of the utility to satisfy duty and update."""
        f, v, l = self._fresh, self._vap, self._liq
        v.P = l.P = f.P
        if f.phase == 'l':
            l.copylike(f)
            u = l
        else:
            v.copylike(f)
            u = v
        u.T = self._T_exit(T_pinch, self.dT, T_limit, duty > 0)
        self._update_utility_flow(f, u, duty)
        f.mol = u.mol

    def _update_flow_wt_phase_change(self, duty):
        """Change phase of utility, calculate and set minimum net flowrate of the utility to satisfy duty and update."""
        f, v, l = self._fresh, self._vap, self._liq
        v.P = l.P = f.P
        v.T = l.T = f.T
        f_is_liq = f.phase == 'l'
        duty_positive = duty > 0
        if duty_positive and not f_is_liq:
            u = l
        elif duty_positive < 0 and f_is_liq:
            u = v
        else:
            sign = 'positive' if duty_positive else 'negative'
            raise ValueError(
                f"Fresh utility stream phase is '{f.phase}' while duty is {sign}, consider switching phase")
        self._update_utility_flow(f, u, duty)
        f.mol = u.mol

    # Subcalculations
    @staticmethod
    def _update_utility_flow(fresh, utility, duty):
        """Changes flow rate of utility such that it can meet the duty requirement"""
        utility.mol = fresh.mol
        dH = fresh.H - utility.H
        utility.mol *= abs(duty/dH)

    @staticmethod
    def _T_exit(T_pinch, dT, T_limit, duty_positive):
        """Return exit temperature of a stream in a counter current heat exchanger

        **Parameters**

             **T_pinch:** [float] Pinch temperature of process stream (K)

             **dT:** [float] Pinch temperature difference (K)

             **duty_positve:** [bool] True if exit temperature should be higher (stream is loosing energy)

        """
        if duty_positive:
            T_exit = T_pinch + dT
            if T_limit and T_limit > T_exit:
                T_exit = T_limit
        else:
            T_exit = T_pinch - dT
            if T_limit and T_limit < T_exit:
                T_exit = T_limit
        return T_exit

    # Representation
    def _info(self, **show_units):
        """Return string related to specifications"""
        # Get units of measure
        su = show_units
        units = self._units
        Duty = su.get('duty') or units['duty']
        Flow = su.get('flow') or units['flow']
        Cost = su.get('cost') or units['cost']
        
        # Select flow dimensionality
        flow_dim = Q_(0, Flow).dimensionality
        if flow_dim == mol_flow_dim:
            flowattr = 'molnet'
        elif flow_dim == mass_flow_dim:
            flowattr = 'massnet'
        elif flow_dim == vol_flow_dim:
            flowattr = 'volnet'
        else:
            raise DimensionError(f"Dimensions for flow units must be in molar, mass or volumetric flow rates, not '{flow_dim}'.")
        
        # Change units and return info string
        if self.ID:
            u_in = self._fresh
            flownet = getattr(u_in, flowattr)
            flowunits = u_in.units[flowattr]
            duty = Q_(self.duty, units['duty']).to(Duty).magnitude
            flow = Q_(flownet, flowunits).to(Flow).magnitude
            cost = Q_(self.cost, units['cost']).to(Cost).magnitude
            return (f'{type(self).__name__}: {self.ID}\n'
                   +f' duty:{duty: .3g} {Duty}\n'
                   +f' flow:{flow: .3g} {Flow}\n'
                   +f' cost:{cost: .3g} {Cost}')
        else:
            return (f'{type(self).__name__}: None\n'
                   +f' duty: 0\n'
                   +f' flow: 0\n'
                   +f' cost: 0')

    def show(self, **show_units):
        """Print all specifications"""
        print(self._info(**show_units))
        
    def __repr__(self):
        if self.ID:
            units = self._units
            ud = units['duty']
            uf = units['flow']
            uc = units['cost']
            return f'<{self.ID}: {self.duty:.3g} {ud}, {self.flow:.3g} {uf}, {self.cost:.3g} {uc}>'
        else:
            return f'<{type(self).__name__}: None>'
