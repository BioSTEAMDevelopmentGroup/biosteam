# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 15:53:14 2018

@author: yoelr
"""
import numpy as np
from fluids.pump import nema_sizes_hp
from .design_tools.mechanical import calculate_NPSH, pump_efficiency, nearest_NEMA_motor_size
from .design_tools.specification_factors import (
    pump_material_factors,
    pump_centrifugal_factors)
from .._unit import Unit
from ..utils import static_flow_and_phase
import biosteam as bst

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
        [0] Input stream
    outs :
        [1] Output stream
    P=101325 : float, optional
        Pressure of output stream (Pa). If the pressure of the outlet is 
        the same as the inlet, pump is design to increase of pressure
        to `P_design`.
    pump_type='Default' : str
        If 'Default', a pump type is selected based on heuristics.
    material='Cast iron' : str
        Contruction material of pump.
    P_design=405300 : float
        Design pressure for pump.
    ignore_NPSH=True : bool
        Whether to take into consideration NPSH in the selection of a
        pump type.
    
    Notes
    -----
    Default pump selection and design and cost algorithms are based on [0]_.
    
    Examples
    --------
    Simulate Pump for pressure increase:
    
    >>> from thermosteam import Chemicals, Stream, settings
    >>> from biosteam.units import Pump
    >>> chemicals = Chemicals(['Water', 'Ethanol'])
    >>> settings.set_thermo(chemicals)
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
    
    References
    ----------
    .. [0] Seider, Warren D., et al. (2017). "Cost Accounting and Capital Cost
        Estimation". In Product and Process Design Principles: Synthesis,
        Analysis, and Evaluation (pp. 450-455). New York: Wiley.
    
    """
    _units = {'Ideal power': 'hp',
              'Power': 'hp',
              'Head': 'ft',
              'NPSH': 'ft',
              'Flow rate': 'gpm'}
    BM = 3.3

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
            self._F_M = pump_material_factors[material]
        except KeyError:
            raise ValueError("material must be one of the following: "
                            f"{', '.join(pump_material_factors)}")
        self._F_Mstr = material  
        
    def __init__(self, ID='', ins=None, outs=(), P=101325,
                 pump_type='Default',
                 material='Cast iron',
                 P_design=405300,
                 ignore_NPSH=True):
        Unit.__init__(self, ID, ins, outs)
        self.P = P
        self.pump_type = pump_type
        self.material = material
        self.P_design = P_design
        self.ignore_NPSH = ignore_NPSH
    
    def _setup(self):
        s_in, = self.ins
        s_out, = self.outs
        s_out.P = self.P
    
    def _run(self):
        s_in, = self.ins
        s_out, = self.outs
        s_out.T = s_in.T
    
    def _design(self):
        Design = self.design_results
        si, = self.ins
        so, = self.outs
        Pi = si.P
        Po = so.P
        Qi = si.F_vol
        mass = si.F_mass
        nu = si.nu
        dP = Po - Pi
        if dP < 1: dP = self.P_design - Pi
        Design['Ideal power'] = power_ideal = Qi*dP*3.725e-7 # hp
        Design['Flow rate'] = q = Qi*4.403 # gpm
        if power_ideal <= max_hp:
            Design['Efficiency'] = efficiency = pump_efficiency(q, power_ideal)
            Design['Actual power'] = power = power_ideal/efficiency
            Design['Pump power'] = nearest_NEMA_motor_size(power)
            Design['N'] = N = 1
            Design['Head'] = head = power/mass*897806 # ft
            # Note that:
            # head = power / (mass * gravity)
            # head [ft] = power[hp]/mass[kg/hr]/9.81[m/s^2] * conversion_factor
            # and 897806 = (conversion_factor/9.81)
        else:
            power_ideal /= 2
            q /= 2
            if power_ideal <= max_hp:
                Design['Efficiency'] = efficiency = pump_efficiency(q, power_ideal)
                Design['Actual power'] = power = power_ideal/efficiency
                Design['Pump power'] = nearest_NEMA_motor_size(power)
                Design['N'] = N = 2
                Design['Head'] = head = power/mass*897806 # ft
            else:
                raise NotImplementedError('more than 2 pump required, but not yet implemented')
        
        if self.ignore_NPSH:
            NPSH_satisfied = True
        else:
            Design['NPSH'] = NPSH = calculate_NPSH(Pi, si.P_vapor, si.rho)
            NPSH_satisfied = NPSH > 1.52
        
        # Get type
        pump_type = self.pump_type
        if pump_type == 'Default':
            if (0.00278 < q < 5000
                and 15.24 < head < 3200
                and nu < 0.00002
                and NPSH_satisfied):
                pump_type = 'Centrifugal'
            elif (q < 1500
                  and head < 914.4
                  and 0.00001 < nu < 0.252):
                pump_type = 'Gear'
            elif (head < 20000
                  and q < 500
                  and power < 200
                  and nu < 0.01):
                pump_type = 'MeteringPlunger'
            else:
                NPSH = calculate_NPSH(Pi, si.P_vapor, si.rho)
                raise NotImplementedError(f'no pump type available at current power ({power:.3g} hp), flow rate ({q:.3g} gpm), and head ({head:.3g} ft), kinematic viscosity ({nu:.3g} m2/s), and NPSH ({NPSH:.3g} ft)')
                
        Design['Type'] = pump_type
        self.power_utility(power/N/1.341) # Set power in kW
    
    def _cost(self):
        # Parameters
        Design = self.design_results
        Cost = self.purchase_costs
        pump_type = Design['Type']
        q = Design['Flow rate']
        h = Design['Head']
        p = Design['Pump power']
        F_M = self._F_M
        I = bst.CE/567
        lnp = ln(p)
        
        # TODO: Add cost equation for small pumps
        # Head and flow rate is too small, so make conservative estimate on cost
        if q < 50: q = 50
        if h < 50: h = 50
        
        # Cost pump
        if 'Centrifugal' in pump_type:
            # Find pump factor
            F_Tdict = pump_centrifugal_factors
            F_T = 1 # Assumption
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
            else:
                raise NotImplementedError(f'no centrifugal pump available at current power ({p:.3g} hp), flow rate ({q:.3g} gpm), and head ({h:.3g} ft)')
            S = q*h**0.5 # Size factor
            S_new = S if S > 400 else 400
            lnS = ln(S_new)
            Cb = exp(12.1656-1.1448*lnS+0.0862*lnS**2)
            Cb *= S/S_new
            Cost['Pump'] = F_M*F_T*Cb*I
        elif pump_type == 'Gear':
            q_new = q if q > 50 else 50
            lnq = ln(q_new)
            Cb = exp(8.2816 - 0.2918*lnq + 0.0743*lnq**2)
            Cb *= q/q_new
            Cost['Pump'] = F_M*Cb*I
        elif pump_type == 'MeteringPlunger':
            Cb = exp(7.9361 + 0.26986*lnp + 0.06718*lnp**2)
            Cost['Pump'] = F_M*Cb*I
        
        # Cost electric motor
        lnp2 = lnp**2
        lnp3 = lnp2*lnp
        lnp4 = lnp3*lnp
        Cost['Motor'] = exp(5.9332 + 0.16829*lnp
                            - 0.110056*lnp2 + 0.071413*lnp3
                            - 0.0063788*lnp4)*I
    
        
        
        