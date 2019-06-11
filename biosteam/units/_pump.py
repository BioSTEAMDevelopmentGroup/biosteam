# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 15:53:14 2018

@author: yoelr
"""
import numpy as np
from .._unit import Unit, metaUnit
from fluids.pump import nema_sizes_hp
from .designtools._vacuum import _calc_MotorEfficiency, _calc_BreakEfficiency

ln = np.log
exp = np.exp

# %% Data

max_hp = nema_sizes_hp[-1]

# Material factors
F_Mdict = {'Cast iron':       1,
           'Ductile iron':    1.15,
           'Cast steel':      1.35,
           'Bronze':          1.9,
           'Stainless steel': 2,
           'Hastelloy C':     2.95,
           'Monel':           3.3,
           'Nickel':          3.5,
           'Titanium':        9.7}

# Gear factor 
F_Tgear = {'OpenDripProof':           1,
           'EnclosedFanCooled':       1.4,
           'ExplosionProofEnclosure': 1.8}

# Centrifugal factor
F_Tcentrifugal = {'VSC3600':   1,
                  'VSC1800':   1.5,
                  'HSC3600':   1.7,
                  'HSC1800':   2,
                  '2HSC3600':  2.7,
                  '2+HSC3600': 8.9}

# %% Pump selection

# Pump types
Types = ('CentrifugalSingle', 'CentrifugalDouble', 'Gear')


def calc_NPSH(P_suction, P_vapor, rho_liq):
    """Return NPSH in ft given suction and vapor pressure in Pa and density in kg/m^3."""
    # Note: NPSH = (P_suction - P_vapor)/(rho_liq*gravity)
    # Taking into account units, NPSH will be equal to return value
    return 0.334438*(P_suction - P_vapor)/rho_liq
    

# %% Classes

class metaPump(metaUnit):
    @property
    def material(cls):
        """Pump material"""
        return Pump._material
    @material.setter    
    def material(cls, material):
        try:
            Pump._F_M = F_Mdict[material]
        except KeyError:
            dummy = str(F_Mdict.keys())[11:-2]
            raise ValueError(f"material must be one of the following: {dummy}")
        Pump._F_Mstr = material   

# TODO: Fix pump selection to include NPSH available and required.
class Pump(Unit, metaclass=metaPump):
    """Create a pump that sets the pressure of the 0th output stream.

    **Parameters**

        **P:** [float] Pressure of output stream (Pa). If None, cost pump as increase of pressure to P_startup.
    
    **ins**
    
        [0] Input stream
        
    **outs**
    
        [1] Output stream
    
    **Examples**
    
        :doc:`Pump Example`
    
    **References**
    
        [0] Seider, Warren D., et al. (2017). "Cost Accounting and Capital Cost Estimation". In Product and Process Design Principles: Synthesis, Analysis, and Evaluation (pp. 450-455). New York: Wiley.
    
    """
    _units = {'Ideal power': 'hp',
              'Power': 'hp',
              'Head': 'ft',
              'NPSH': 'ft',
              'Flow rate': 'gpm'}
    _N_ins = 1
    _N_outs = 1
    _has_power_utility = True
    _linkedstreams = True
    _kwargs = {'P': None}
    
    # Pump type
    _Type = 'Default'
    
    # Material factor
    _F_Mstr = 'Cast iron'
    _F_M = 1
    
    # Gear factor 
    _F_Tstr = 'VSC3600'
    _F_T = 1
    
    #: Pressure for costing purposes (Pa).
    P_startup = 405300 
    
    #: [bool] Ignore minimum NPSH requirement if True.
    ignore_NPSH = True

    @property
    def Type(self):
        """Pump type"""
        return self._Type
    @Type.setter
    def Type(self, Type):
        if Type not in Types:
            Types_str = str(Types)[1:-1]
            raise ValueError('Type must be one of the following: ' + Types_str)
        self._Type = Type
    
    @property
    def material(self):
        """Pump material"""
        return self._material
    @material.setter    
    def material(self, material):
        try:
            self._F_M = F_Mdict[material]
        except KeyError:
            dummy = str(F_Mdict.keys())[11:-2]
            raise ValueError(f"material must be one of the following: {dummy}")
        self._F_Mstr = material  
        
    def _run(self):
        out = self._outs[0]
        feed = self._ins[0]
        out.P = self._kwargs['P'] or feed.P
        out.T = feed.T
        out._phase = feed._phase
    
    def _design(self):
        Design = self._results['Design']
        si = self.ins[0]
        so = self.outs[0]
        Pi = si.P
        Po = so.P
        Qi = si.volnet
        mass = si.massnet
        nu = si.nu
        
        if abs(Po - Pi) < 1: Po = self.P_startup
        Design['Ideal power'] = power_ideal = Qi*(Po - Pi)*3.725e-7 # hp
        Design['Flow rate'] = q = Qi*4.403 # gpm
        if power_ideal <= max_hp:
            Design['Efficiency'] = e = self._calc_Efficiency(q, power_ideal)
            Design['Actual power'] = power = power_ideal/e
            Design['Pump power'] = self._nearest_PumpPower(power)
            Design['N'] = 1
            Design['Head'] = head = power/mass*897806 # ft
            # Note that:
            # head = power / (mass * gravity)
            # head [ft] = power[hp]/mass[kg/hr]/9.81[m/s^2] * conversion_factor
            # and 897806 = (conversion_factor/9.81)
        else:
            raise NotImplementedError('more than 1 pump required, but not yet implemented')
        
        if self.ignore_NPSH:
            NPSH_satisfied = True
        else:
            Design['NPSH'] = NPSH = calc_NPSH(Pi, si.P_vapor, si.rho)
            NPSH_satisfied = NPSH > 1.52
        
        # Get type
        Type = self.Type
        if Type == 'Default':
            if (0.00278 < q < 5000
                and 15.24 < head < 3200
                and nu < 0.00001
                and NPSH_satisfied):
                Type = 'Centrifugal'
            elif (q < 1500
                  and head < 914.4
                  and 0.00001 < nu < 0.252):
                Type = 'Gear'
            elif (head < 20000
                  and q < 500
                  and power < 200
                  and nu < 0.01):
                Type = 'MeteringPlunger'
            else:
                raise NotImplementedError(f'no pump type available at current power ({power:.3g} hp), flow rate ({q:.3g} gpm), and head ({head:.3g} ft), kinematic viscosity ({nu:.3g} m2/s), and NPSH ({NPSH:.3g} ft)')
                
        Design['Type'] = Type
        self._power_utility(power/1.341) # Set power in kW
    
    def _cost(self):
        # Parameters
        results = self._results
        Cost = results['Cost']
        Design = results['Design']
        Type = Design['Type']
        q = Design['Flow rate']
        h = Design['Head']
        p = Design['Pump power']
        F_M = self._F_M
        I = self.CEPCI/567
        lnp = ln(p)
        
        # TODO: Add cost equation for small pumps
        # Head and flow rate is too small, so make conservative estimate on cost
        if q < 50: q = 50
        if h < 50: h = 50
        
        # Cost pump
        if 'Centrifugal' in Type:
            # Find pump factor
            F_Tdict = F_Tcentrifugal
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
        elif Type == 'Gear':
            q_new = q if q > 50 else 50
            lnq = ln(q_new)
            Cb = exp(8.2816 - 0.2918*lnq + 0.0743*lnq**2)
            Cb *= q/q_new
            Cost['Pump'] = F_M*Cb*I
        elif Type == 'MeteringPlunger':
            Cb = exp(7.9361 + 0.26986*lnp + 0.06718*lnp**2)
            Cost['Pump'] = F_M*Cb*I
        
        # Cost electric motor
        lnp2 = lnp**2
        lnp3 = lnp2*lnp
        lnp4 = lnp3*lnp
        Cost['Motor'] = exp(5.9332 + 0.16829*lnp - 0.110056*lnp2 + 0.071413*lnp3 - 0.0063788*lnp4)*I
        return Cost    
    
    @staticmethod
    def _calc_PowerFlow(Qi, Qo, Dpipe, mass=None):
        """Return ideal power due to flow rate change.
        
        **Parameters**
            
            Qi: [Stream] Input flow rate
            
            Qo: [Stream] Output flow rate
            
            Dpipe: [float] Pipe inside diameter
            
            mass: [float] Mass flow rate
        
        """
        if mass:
            A = np.pi*((Dpipe/2)**2)
            # Kinetic energy change term
            vi = Qi/A # velocity in
            vo = Qo/A # velocity out
            K_term = (mass*vo**2 - mass*vi**2)/2
        else:
            K_term = 0
        
        return K_term
    
    @staticmethod
    def _nearest_PumpPower(p:'hp'):
        for power in nema_sizes_hp:
            if power >= p: return power
        return power
        
    @staticmethod
    def _calc_Efficiency(q:'gpm', p:'hp'):
        mup = _calc_BreakEfficiency(q)
        mum = _calc_MotorEfficiency(p/mup)
        return mup*mum
        
        