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
            'price (USD/kJ)', 'price (USD/kg)', 'efficiency')

_cooling_index = ('Cooling water',
                  'Chilled water',
                  'Chilled Brine')

_cooling_water = ('sensible',      # Type [str]: {'sensible', 'latent'}
                  _Water,           # species [Species]
                  (1,),            # flow: [tuple]
                  305.372,         # T (K)
                  101325,          # P (Pa)
                  'liq',           # phase: {'g', 'l'}
                  324.817,         # T limit (K)
                  0,               # price (USD/kJ)
                  2.708e-5,        # price (USD/kg)
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
                       0.01320,
                       0.90)

_medium_pressure_steam = ('latent',
                          _Water,
                          (1,),
                          454.484,
                          1034214.0,
                          'gas',
                          None,
                          0,
                          0.01530,
                          0.80)

_high_pressure_steam = ('latent',
                        _Water,
                        (1,),
                        508.858,
                        3102642.0,
                        'gas',
                        None,
                        0,
                        0.01760,
                        0.70)


# %% Utility classes


class HeatUtility:
    """Create an HeatUtility object that can choose a utility stream and calculate utility requirements. It can calculate required flow rate, temperature, or phase change of utility. Calculations assume counter current flow rate.
    
    **Parameters**
    
        **source:** [Object or str] source of the heat utility
        
        **efficiency:** [float] Fraction of heat transfered after accounting for heat loss.
    
    **Class Attributes**
    
        **cooling_agents:** [DataFrame] All cooling utilities available
        
        **heating_agents:** [DataFrame] All heating utilities available
    
    **Examples**
    
        Create a heat utility with a source:
            
        .. code-block:: python
        
           >>> heat_util = HeatUtility('Flash')
           >>> heat_util.show()
           HeatUtility: None
            Duty: 0
            Flow: 0
            Cost: 0
        
        Calculate utility requirement by calling it with a duty (kJ/hr) and temperature (K):
            
        .. code-block:: python
        
           >>> heat_util(1000, 300, 350)
           {'Cost': 0.00607 (USD/hr),
            'Flow': 0.46 (kg/hr),
            'Duty': 1e+03 (kJ/hr),
            'ID': Low pressure steam}
       
        All results are cached:
            
        .. code-block:: python
        
           >>> heat_util.results
           {'Cost': 0.00607 (USD/hr),
            'Flow': 0.46 (kg/hr),
            'Duty': 1e+03 (kJ/hr),
            'ID': Low pressure steam}
           
    """
    __slots__ = ('_fresh', '_vap', '_liq', 'results', 'efficiency')
    dT = 5  #: [float] Pinch temperature difference
    
    #: Units of measure for results dictionary
    _units = UnitManager([], Duty='kJ/hr', Flow='kg/hr', Cost='USD/hr')

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

    def __init__(self, source, efficiency=None):
        self.results = SmartBook(self._units, ID='', Cost=0,
                                 Flow=0, Duty=0, source=source)
        self.efficiency = efficiency

    def _init_streams(self, ID, flow, species, T, P, phase):
        """Initialize utility streams."""
        ID = '*'+ID
        self._fresh = Stream(ID, species=species,
                             flow=flow,
                             T=T,
                             P=P,
                             phase=phase)
        self._liq = Stream(ID, species=species, phase='l')
        self._vap = Stream(ID, species=species, phase='g')

    def __call__(self, duty:'kJ/hr', T_in:'K', T_out:'K'=None) -> 'results [dict]':
        """Return dictionary of utility requirements given the essential parameters.
        
        **Parameters**
        
            **duty:** [float] Unit duty requirement (kJ/hr)
            
            **T_in:** [float] Entering process stream temperature (K)
            
            **T_out:** [float] Exit process stream temperature (K)
        
        **Returns**
            * 'ID': ID of utility
            * 'Duty': Unit duty requirement, including heat transfer losses (kJ/hr)
            * 'Flow': HeatUtility flow rate (kg/hr)
            * 'Cost': Cost of utility (USD/hr)
        
        """
        r = self.results
        
        # Set pinch and operating temperature
        if T_out is None:
            T_pinch = T_op = T_in
        else:
            T_pinch, T_op = self._get_pinch(duty, T_in, T_out)
        
        
        # Select heat transfer agent
        if duty < 0:
            (Type, price_duty,
             price_mass, T_limit,
             efficiency) = self._select_cooling_agent(T_op)
        elif duty > 0:
            (Type, price_duty,
            price_mass, T_limit,
            efficiency) = self._select_heating_agent(T_op)
        else:
            return r
        
        # Calculate utility flow rate requirement
        efficiency = self.efficiency if self.efficiency else efficiency
        duty = duty/efficiency
        if Type == 'sensible':
            self._update_flow_wt_pinch_T(duty, T_pinch, T_limit)
        elif Type == 'latent':
            self._update_flow_wt_phase_change(duty)
        else:
            raise biosteamError(f"Invalid heat transfer agent '{self._fresh.ID[1:]}'.")
        
        # Update and return results
        ID = self._fresh.ID[1:]
        r['ID'] = ID
        r['Flow'] = mass = self._fresh.massnet
        r['Duty'] = duty
        r['Cost'] = price_duty*duty + price_mass*mass
        return r

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
             price_mass, efficiency) = cooling_agents[ID]
            if T_max > T:
                self._init_streams(ID, flow, species, T, P, phase[0])
                return Type, price_duty, price_mass, T_limit, efficiency
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
             price_mass, efficiency) =  heating_agents[ID]
            if T_min < T:
                self._init_streams(ID, flow, species, T, P, phase[0])
                return Type, price_duty, price_mass, T_limit, efficiency
        raise biosteamError(f'No heating agent that can heat over {T_pinch} K')

    # Main Calculations
    def _update_flow_wt_pinch_T(self, duty, T_pinch, T_limit):
        """Set utility Temperature at the pinch, calculate and set minimum net flowrate of the utility to satisfy duty and update."""
        f, v, l = self._fresh, self._vap, self._liq
        v.P = l.P = f.P
        if f.phase == 'l':
            l.copy_like(f)
            u = l
        else:
            v.copy_like(f)
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
        r = self.results
        Duty = su.get('Duty') or su.get('duty') or r.units['Duty']
        Flow = su.get('Flow') or su.get('flow') or r.units['Flow']
        Cost = su.get('Cost') or su.get('cost') or r.units['Cost']
        
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
        if hasattr(self,'_fresh') and self.results:
            u_in = self._fresh
            flownet = getattr(u_in, flowattr)
            flowunits = u_in.units[flowattr]
            duty = Q_(r['Duty'], r.units['Duty']).to(Duty).magnitude
            flow = Q_(flownet, flowunits).to(Flow).magnitude
            cost = Q_(r['Cost'], r.units['Cost']).to(Cost).magnitude
            ID = r['ID']
            return (f'{type(self).__name__}: {ID}\n'
                   +f' Duty:{duty: .3g} {Duty}\n'
                   +f' Flow:{flow: .3g} {Flow}\n'
                   +f' Cost:{cost: .3g} {Cost}')
        else:
            return (f'{type(self).__name__}: None\n'
                   +f' Duty: 0\n'
                   +f' Flow: 0\n'
                   +f' Cost: 0')

    def show(self, **show_units):
        """Print all specifications"""
        print(self._info(**show_units))

    def __repr__(self):
        if hasattr(self, '_fresh'):
            name = self._fresh.ID.replace('*','')
            return f'<{type(self).__name__}: {name}>'
        else:
            return f'<{type(self).__name__}: None>'
