# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import numpy as np
from .design_tools.mechanical import (
    calculate_NPSH,
    pump_efficiency, 
    nearest_NEMA_motor_size,
    nema_sizes_hp,
    electric_motor_cost
)
from .design_tools.specification_factors import (
    pump_material_factors,
    pump_centrifugal_factors)
from .._unit import Unit
from ..utils import static_flow_and_phase
from math import ceil
from warnings import warn
import biosteam as bst

__all__ = ('Pump',)

ln = np.log
exp = np.exp

# %% Data

max_hp = nema_sizes_hp[-1]
pump_types = ('Default', 'CentrifugalSingle', 'CentrifugalDouble', 'Gear')
    

# %% Classes

# TODO: Fix pump selection to include NPSH available and required.
@static_flow_and_phase
class Pump(Unit):
    """
    Create a pump that sets the pressure of the outlet.

    Parameters
    ----------
    ins : 
        Inlet.
    outs : 
        Outlet.
    P=101325 : float, optional
        Pressure of output stream (Pa). If the pressure of the outlet is 
        the same as the inlet, pump is designed to increase of pressure
        to `dP_design`.
    pump_type='Default' : str
        If 'Default', a pump type is selected based on heuristics.
    material : str, optional
        Contruction material of pump. Defaults to 'Cast iron'.
    dP_design=101325 : float
        Design pressure change for pump.
    ignore_NPSH=True : bool
        Whether to take into consideration NPSH in the selection of a
        pump type.
    
    Notes
    -----
    Default pump selection and design and cost algorithms are based on [1]_.
    
    Examples
    --------
    Simulate Pump for pressure increase:
    
    >>> from biosteam import Stream, settings
    >>> from biosteam.units import Pump
    >>> settings.set_thermo(['Water', 'Ethanol'], cache=True)
    >>> feed = Stream('feed', Water=200, T=350)
    >>> P1 = Pump('P1', ins=feed, outs='out', P=2e5)
    >>> P1.simulate()
    >>> P1.show()
    Pump: P1
    ins...
    [0] feed
        phase: 'l', T: 350 K, P: 101325 Pa
        flow (kmol/hr): Water  200
    outs...
    [0] out
        phase: 'l', T: 350 K, P: 200000 Pa
        flow (kmol/hr): Water  200
    >>> P1.results()
    Pump                              Units           P1
    Electricity         Power            kW        0.288
                        Cost         USD/hr       0.0225
    Design              Type                 Centrifugal
                        Ideal power      hp        0.136
                        Flow rate       gpm         16.3
                        Efficiency                 0.352
                        Power            hp        0.386
                        Head             ft         96.3
                        Motor size       hp          0.5
    Purchase cost       Pump            USD     4.37e+03
                        Motor           USD          273
    Total purchase cost                 USD     4.64e+03
    Utility cost                     USD/hr       0.0225
    
    References
    ----------
    .. [1] Seider, Warren D., et al. (2017). "Cost Accounting and Capital Cost
        Estimation". In Product and Process Design Principles: Synthesis,
        Analysis, and Evaluation (pp. 450-455). New York: Wiley.
    
    """
    _units = {'Ideal power': 'hp',
              'Motor size': 'hp',
              'Power': 'hp',
              'Head': 'ft',
              'NPSH': 'ft',
              'Flow rate': 'gpm'}
    _F_BM_default = {'Pump': 3.3,
                     'Motor': 3.3}

    @property
    def pump_type(self):
        """Pump type"""
        return self._pump_type
    @pump_type.setter
    def pump_type(self, pump_type):
        if pump_type not in pump_types:
            raise ValueError('Type must be one of the following: {", ".join(pump_types)}')
        self._pump_type = pump_type
    
    @property
    def material(self):
        """Pump material"""
        return self._material
    @material.setter    
    def material(self, material):
        try:
            self.F_M['Pump'] = pump_material_factors[material]
        except KeyError:
            raise ValueError("material must be one of the following: "
                            f"{', '.join(pump_material_factors)}")
        self._material = material  
        
    def _init(self, 
             P=None,
             pump_type='Default',
             material='Cast iron',
             dP_design=101325,
             ignore_NPSH=True
        ):
        self.P = P
        self.pump_type = pump_type
        self.material = material
        self.dP_design = dP_design
        self.ignore_NPSH = ignore_NPSH
    
    def _run(self):
        s_in, = self.ins
        s_out, = self.outs
        s_out.copy_like(s_in)
        if self.P: s_out.P = self.P 
    
    def _design(self):
        Design = self.design_results
        si, = self.ins
        so, = self.outs
        if si.isempty(): 
            self.design_results.clear()
            return
        Pi = si.P
        Po = so.P
        Qi = si.F_vol
        mass = si.F_mass
        nu = si.nu
        dP = Po - Pi
        if dP < 1: dP = self.dP_design
        power_ideal = Qi*dP*3.725e-7 # hp
        q = Qi*4.403 # gpm
        N = 1
        
        # Note:
        # A pump motor can have maximally `max_hp` power, so need to first solve
        # for how many pumps with maximally `max_hp` are needed. 
        # Assume 100% efficiency to estimate the number of pumps.
        # Additionally, pump types have their own maximum size.
        N_guess = -1
        while N != N_guess:
            N_guess = N
            power_ideal_i = power_ideal / N_guess
            q_i = q / N_guess
            efficiency = pump_efficiency(q_i, power_ideal_i) 
            power_i = power_ideal_i / efficiency
            head_i = N * power_i / mass * 897806 # ft
            # Note that:
            # head = power / (mass * gravity)
            # head [ft] = power[hp]/mass[kg/hr]/9.81[m/s^2] * conversion_factor
            # and 897806 = (conversion_factor/9.81)
            
            if self.ignore_NPSH:
                NPSH_satisfied = True
            else:
                Design['NPSH'] = NPSH = calculate_NPSH(Pi, si.P_vapor, si.rho)
                NPSH_satisfied = NPSH > 1.52
            
            pump_type = self.pump_type
            if pump_type == 'Default':
                if (15.24 <= head_i <= 3200
                    and nu <= 0.00002
                    and NPSH_satisfied):
                    pump_type = 'Centrifugal'
                elif (head_i <= 914.4
                      and 0.00001 <= nu <= 0.252):
                    N = max(ceil(q_i / 1500), N)
                    pump_type = 'Gear'
                elif (head_i <= 20000
                      # and power_i <= 200
                      and 1 <= power_i <= 200
                      and nu <= 0.01):
                    N = max(ceil(q_i / 500), N)
                    pump_type = 'MeteringPlunger'
                else:
                    NPSH = calculate_NPSH(Pi, si.P_vapor, si.rho)
                    warn(f'{repr(self)} no pump type available at current power '
                         f'({power_i:.3g} hp), head ({head_i:.3g} ft), kinematic '
                         f'viscosity ({nu:.3g} m2/s), and NPSH ({NPSH:.3g} ft); '
                          'assuming centrigugal pump',
                         RuntimeWarning)
                    pump_type = 'Centrifugal'
        Design['Type'] = pump_type
        Design['Ideal power'] = power_ideal_i = power_ideal / N
        Design['Flow rate'] = q_i = q / N_guess
        Design['Efficiency'] = pump_efficiency(q_i, power_ideal_i) 
        Design['Power'] = power_i = power_ideal_i / efficiency
        Design['Head'] = N * power_i / mass * 897806 # ft
        Design['Motor size'] = motor_rating = nearest_NEMA_motor_size(power_i)
        self.add_power_utility(power_i / 1.341) # Add power for individual pump in kW
        self.parallel['self'] = N # BioSTEAM will multiply all costs (both capital and utility) by this number
        self.parallel['Motor'] = N * ceil(power_i / motor_rating) # Multiple motors per pump is posible
    
    def _cost(self):
        Design = self.design_results
        if not Design: return
        Cost = self.baseline_purchase_costs
        pump_type = Design['Type']
        q = Design['Flow rate']
        h = Design['Head']
        p = Design['Power']
        I = bst.CE/567
        # TODO: Add cost equation for small pumps
        # Head and flow rate is too small, so make conservative estimate on cost
        if q < 50: q = 50
        if h < 50: h = 50
        F_T = 1 # Assumption
        
        # Cost pump
        if 'Centrifugal' in pump_type:
            # Find pump factor
            F_Tdict = pump_centrifugal_factors
            if p < 75 and 50 <= q <= 900 and 50 <= h <= 400:
                F_T = F_Tdict['VSC3600']
            elif p < 200 and 50 <= q <= 3500 and 50 <= h <= 2000:
                F_T = F_Tdict['VSC1800']
            elif p < 150 and 100 <= q <= 1500 and 100 <= h <= 450:
                F_T = F_Tdict['HSC3600']
            elif p < 250 and 250 <= q <= 5000 and 50 <= h <= 500:
                F_T = F_Tdict['HSC1800']
            elif p < 250 and 50 <= q <= 1100 and 300 <= h <= 1100:
                F_T = F_Tdict['2HSC3600']
            elif p < 1450 and 100 <= q <= 1500 and 650 <= h <= 3200:
                F_T = F_Tdict['2+HSC3600']
            elif p < 1450 and 1500 <= q <= 5000 and 650 <= h <= 3200:
                # TODO: This flow rate is allowable in the design section
                # but not in the estimation of purchase cost.
                # For now, we ignore this problem, but additional code on
                # selecting number of pumps should be added.
                F_T = F_Tdict['2+HSC3600']
            else:
                # TODO: This flow rate is allowable in the design section
                # but not in the estimation of purchase cost.
                # For now, we ignore this problem, but additional code on
                # selecting number of pumps should be added.
                F_T = F_Tdict['2+HSC3600']
                # raise NotImplementedError(f'no centrifugal pump available at current power ({p:.3g} hp), flow rate ({q:.3g} gpm), and head ({h:.3g} ft)')
            S = q*h**0.5 # Size factor
            S_new = S if S > 400 else 400
            lnS = ln(S_new)
            Cb = exp(12.1656-1.1448*lnS+0.0862*lnS**2)
            Cb *= S/S_new
            Cost['Pump'] = Cb*I
        elif pump_type == 'Gear':
            q_new = q if q > 50 else 50
            lnq = ln(q_new)
            Cb = exp(8.2816 - 0.2918*lnq + 0.0743*lnq**2)
            Cb *= q/q_new
            Cost['Pump'] = Cb*I
        elif pump_type == 'MeteringPlunger':
            lnp = ln(p)
            lnp2 = lnp * lnp
            Cb = exp(7.9361 + 0.26986*lnp + 0.06718*lnp2)
            Cost['Pump'] = Cb*I
        
        self.F_D['Pump'] = F_T
        
        # Cost electric motor
        Cost['Motor'] = electric_motor_cost(p)
    
        
        
        
