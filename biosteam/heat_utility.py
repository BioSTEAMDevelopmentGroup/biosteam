# -*- coding: utf-8 -*-
"""
Created on Sat Aug 18 14:25:34 2018

@author: yoelr
"""
from biosteam.exceptions import DimensionError
from biosteam.species import Species
from biosteam.stream import Stream, mol_flow_dim, mass_flow_dim, vol_flow_dim
from biosteam.mixed_stream import MixedStream
from bookkeep import SmartBook, UnitManager
from biosteam import Q_

# %% Default Heat Transfer Streams

CoolingAir = Species('Nitrogen', 'Oxygen',
                     'Water', 'CO2', 'Argon')

Water = Species('Water')

# %% Utilities


class HeatUtility:
    """Create an HeatUtility object that can choose a utility stream and calculate utility requirements. It can calculate required flow rate, temperature, or phase change of utility. Calculations assume counter current flow rate."""
    __slots__ = ('_fresh', '_vap', '_liq', '_phase_ch', '_cost', 
                 '_T_limit', 'results', '_user_specified')
    dT = 5  #: [float] Pinch temperature difference
    
    #: Units of measure for results dictionary
    _units = UnitManager([], Duty='kJ/hr', Flow='kg/hr', Cost='USD/hr')

    def __init__(self):
        self.results = SmartBook(self._units, Cost=0, Flow=0, Duty=0)
        self._T_limit = None
        self._user_specified = False

    def define_utility(self, fresh, vap, liq, phase_ch, cost, T_limit=None):
        """Define utility streams and operation method.
        
        **Parameters**
            
            **fresh:** [Stream] entering utility stream

            **vap:** [Stream] exiting utility vapor

            **liq:** [Stream] exiting utility liquid

            **phase_ch:** [bool] True if utility changes phase

            **cost:** [function] Calculates net cost in $/hr

            **T_limit:** [float] Minimum/maximum exit temperature for sensible heat transfer agents
        
        """
        self._fresh, self._vap, self._liq, self._phase_ch, self._T_limit = fresh, vap, liq, phase_ch, T_limit
        self._user_specified = True

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
            if T_limit:
                if T_limit > T_exit:
                    T_exit = T_limit
        else:
            T_exit = T_pinch - dT
            if T_limit:
                if T_limit < T_exit:
                    T_exit = T_limit
        return T_exit

    def _set_CoolingAir(self):
        """Set cooling air as the utility."""
        ID = '*Cooling Air'
        self._fresh = Stream(ID, species=CoolingAir,
                                flow=[0.7809, 0.2095, 0.004, 0.0004, 0.0093],
                                T=305.372,
                                P=101325,
                                phase='g')
        self._liq = Stream(ID, species=CoolingAir, phase='l')
        self._vap = Stream(ID, species=CoolingAir,phase='g')


    def _set_Water(self, ID, T, P, phase):
        """Set water as the utility."""
        ID = '*'+ID
        self._fresh = Stream(ID, species=Water,
                             flow=[1],
                             T=T,
                             P=P,
                             phase=phase)
        self._liq = Stream(ID, species=Water, phase='l')
        self._vap = Stream(ID, species=Water, phase='g')

    # Default Heat Transfer agents
    def _default_cooling_agent(self, duty, T_pinch):
        """Select a cooling agent that works at set temperature"""
        # Costs from Table 17.1 in Warren D.  Seider et al.-Product and Process Design Principles_ Synthesis, Analysis and Evaluation-Wiley (2016)
        # ^This table was made using data from Busche, 1995
        # Entry temperature conditions of coolants taken from Table 12.1 in Warren, 2016
        dt = 2*self.dT
        if T_pinch > 305.372 + dt:
            self._set_Water('Cooling Water', 305.372, 101325, 'l')
            def cost(): return self._fresh.massnet*0.027/997
            self._cost = cost
            self._T_limit = 324.817
        elif T_pinch > 280.372 + dt:
            self._set_Water('Chilled Water', 280.372, 101325, 'l')
            def cost(): return -duty/10**6 * 5  # in $ per gigajoule
            self._cost = cost
        elif T_pinch > 255.372 + dt:
            self._set_Water('Chilled Brine', 255.372, 101325, 'l')
            def cost(): return -duty/10**6 * 8.145  # in $ per gigajoule
            self._cost = cost
        else:
            raise Exception(f'No deault cooling agent that can cool under {T_pinch} K')

    def _default_heating_agent(self, T_pinch):
        """Select a heating agent that works at set temperature"""
        # Costs from Table 17.1 in Warren D.  Seider et al.-Product and Process Design Principles_ Synthesis, Analysis and Evaluation-Wiley (2016)
        # ^This table was made using data from Busche, 1995
        # Entry temperature conditions of coolants taken from Table 12.1 in Warren, 2016
        dt = 2*self.dT
        if T_pinch < 411.494 - dt:
            self._set_Water('Low Pressure Steam', 411.494, 344738.0, 'g')
            self._cost = lambda: self._fresh.massnet*0.01320
        elif T_pinch < 454.484 - dt:
            self._set_Water('Medium Pressure Steam', 454.484, 1034214.0, 'g')
            self._cost = lambda: self._fresh.massnet*0.01530
        elif T_pinch < 508.858 - dt:
            self._set_Water('High Pressure Steam', 508.859, 3102642.0, 'g')
            self._cost = lambda: self._fresh.massnet*0.01760
        else:
            raise Exception(f'No default heating agent that can heat over {T_pinch} K')

    # Main Calculations
    def _update_flow_wt_pinch_T(self, duty, T_pinch):
        """Set utility Temperature at the pinch, calculate and set minimum net flowrate of the utility to satisfy duty and update."""
        f, v, l = self._fresh, self._vap, self._liq
        v.P = l.P = f.P
        if f.phase == 'l':
            l.copy_like(f)
            u = l
        else:
            v.copy_like(f)
            u = v
        u.T = self._T_exit(T_pinch, self.dT, self._T_limit, duty > 0)
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

    @staticmethod
    def _get_pinch(duty, T_in, T_out) -> '(T_pinch, T_op)':
        """Return pinch temperature and operating temperature."""
        if duty < 0:
            return (T_in, T_out) if T_in > T_out else (T_out, T_in)
        else:
            return (T_in, T_out) if T_in < T_out else (T_out, T_in)

    def _default_run(self, duty, T_pinch, T_op):
        """Select a heat transfer agent, calculate and set required utility."""
        # Select heating agent
        if duty < 0:
            self._default_cooling_agent(duty, T_op)
        else:
            self._default_heating_agent(T_op)
        
        # Calculate utility flow rate requirement
        if self._fresh.phase == 'l':
            self._update_flow_wt_pinch_T(duty, T_pinch)
        else:
            self._update_flow_wt_phase_change(duty)

    def __call__(self, duty:'kJ/hr', T_in:'K', T_out:'K'=None) -> 'results [dict]':
        """Return dictionary of utility requirements given the essential parameters.
        
        **Parameters**
        
            **duty:** [float] Unit duty requirement (kJ/hr)
            
            **T_in:** [float] Entering process stream temperature (K)
            
            **T_out:** [float] Exit process stream temperature (K)
        
        **Returns**
            * 'HeatUtility': Name of utility ()
            * 'Duty': Unit duty requirement (kJ/hr)
            * 'Flow': HeatUtility flow rate (kg/hr)
            * 'Cost': Cost of utility (USD/hr)
        
        """
        # Set pinch and operating temperature
        if T_out is None:
            T_pinch = T_op = T_in
        else:
            T_pinch, T_op = self._get_pinch(duty, T_in, T_out)
        
        if self._user_specified:
            # Run utility as set by user
            if self._phase_ch:
                self._update_flow_wt_phase_change(duty)
            else:
                self._update_flow_wt_pinch_T(duty, T_pinch)
        else:
            # Choose a utility stream and find flow rate and exit temperature/phase
            self._default_run(duty, T_pinch, T_op)
        
        # Update and return results
        r = self.results
        u_T = self._fresh.T
        name = self._fresh.ID.replace('*','')
        r['HeatUtility'] = f'{name} at {u_T:.0f} K'
        r['Flow'] = self._fresh.massnet
        r['Duty'] = duty
        r['Cost'] = self._cost()
        return r

    # Representation
    def _info(self, **show_units):
        # Get units of measure
        su = show_units
        r = self.results
        Duty = su.get('Duty') or su.get('duty') or r.units['Duty']
        Flow = su.get('Flow') or su.get('flow') or r.units['Flow']
        Cost = su.get('Cost') or su.get('cost') or r.units['Cost']
        T = su.get('T') or 'K'
        
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
            u_T = Q_(self._fresh.T, 'K').to(T).magnitude
            name = self._fresh.ID.replace('*','')
            return (f'{type(self).__name__}: {name} at T={u_T:.0f} {T}\n'
                   +f' Duty:{duty: .3g} {Duty}\n'
                   +f' Flow:{flow: .3g} {Flow}\n'
                   +f' Cost:{cost: .3g} {Cost}')
        else:
            return (f'{type(self).__name__}: None\n'
                   +f' Duty: None\n'
                   +f' Flow: None\n'
                   +f' Cost: None')

    def show(self, **show_units):
        print(self._info(**show_units))

    def __repr__(self):
        if hasattr(self, '_fresh'):
            name = self._fresh.ID.replace('*','')
            return f'<{type(self).__name__}: {name}>'
        else:
            return f'<{type(self).__name__}: None>'
