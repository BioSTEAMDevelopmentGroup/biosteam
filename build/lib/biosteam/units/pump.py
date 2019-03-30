# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 15:53:14 2018

@author: yoelr
"""
from biosteam import Unit
from biosteam import np
from fluids.pump import nema_sizes_hp
from .designtools.vacuum import _calc_MotorEfficiency, _calc_BreakEfficiency

ln = np.log
exp = np.exp

max_hp = nema_sizes_hp[-1]

# Ranges of flow rate (m3/hr) and working pressure where pumps work/
# Assume lower minimum flow rate and pressure for costing in Centrifugal single and Regenerative
pump_ranges = {'CentrifugalSingle': ((0.0000, 2271  ), (      0, 689500  )), 
               'CentrifugalDouble': ((0.2271, 2271  ), (      0, 6895000 )),
               'Gear':              ((0.2271, 681.3 ), (  68950, 13790000)),
               'MeteringPlunger':   ((0.2271, 6.813 ), (  68950, 20685000)),
               'Screw':             ((0.2271, 13.626), (  68950, 5516000 )),
               'MeteringDiaphragm': ((0.2271, 1.1355), (  68950, 5516000 )),
               'DirectactingSteam': ((0.2271, 454.2 ), (  68950, 13790000)),
               'AxialFlow':         (( 227.1, 2271  ), (  68950, 137900  ))}

# Regenerative pump (Cost not yet implemented)
regenerative_ranges = ((0.0100, 34.065), (1379000, 6895000 )) 

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

# Pump types
Types = ('CentrifugalSingle', 'CentrifugalDouble', 'Gear', 'MeteringPlunger', 'Default')

class Pump(Unit):
    """Create a pump that sets the pressure of the 0th output stream.

    **Parameters**

        **P:** *[float]* Pressure of output stream (Pa). If None, cost pump as increase of pressure to P_startup.
    
    **ins**
    
        [0] Input stream
        
    **outs**
    
        [1] Output stream
    
    **Examples**
    
        :doc:`Pump Example`
    
    """
    _N_ins = 1
    _N_outs = 1
    _has_power_utility = True
    _has_linked_streams = True
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

    @property
    def Type(self):
        """Pump type"""
        return self._Type
    @Type.setter
    def Type(self, Type):
        if Type not in Types:
            Types_str = str(Types)[1:-1]
            raise ValueError('Pump type must be one of the following: ' + Types_str)
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
            raise ValueError(f"Material must be one of the following: {dummy}")
        self._F_Mstr = material  
        
    def _run(self):
        P = self._kwargs['P']
        for i, o in zip(self.ins, self.outs):
            o.copylike(i)
            if P: o.P = P
    
    def _operation(self):
        """
        * 'Type': [str] Pump type
        * 'Ideal power': (hp)
        * 'Efficiency': ()
        * 'Power': (hp)
        * 'Head': (ft)
        * 'Flow rate': (gpm)
        """
        Oper = self.results['Operation']
        si = self.ins[0]
        so = self.outs[0]
        Pi = si.P
        Po = so.P
        Qi = si.volnet
        mass = si.massnet
        
        # Get type
        Type = self.Type
        if Type == 'Default':
            Type = self._default_type(Pi, Po, Qi)
        Oper['Type'] = Type
        
        # Get ideal power
        if abs(Po - Pi) < 1:    
            Po = self.P_startup
        power_ideal = self._calc_PowerPressure(Pi, Po, Qi)*3.725e-7
        Oper['Ideal power'] = power_ideal # hp
        
        Oper['Flow rate'] = q = Qi*4.403
        if power_ideal <= max_hp:
            Oper['Efficiency'] = e = self._calc_Efficiency(q, power_ideal)
            Oper['Power'] = power =  self._nearest_PumpPower(power_ideal/e)
            Oper['N_pumps'] = 1
            Oper['Head'] = self._calc_Head(power, mass)
        else:
            raise Exception('More than 1 pump required, but not yet implemented.')
        
        self.power_utility(power/1.341)
        return Oper
    
    def _cost(self):
        """
        * 'Pump': (USD)
        * 'Motor': (USD)
        """
        # Parameters
        results = self.results
        Cost = results['Cost']
        Oper = results['Operation']
        Type = Oper['Type']
        q = Oper['Flow rate']
        h = Oper['Head']
        p = Oper['Power']
        F_M = self._F_M
        I = self.CEPCI/567
        lnp = ln(p)
        
        # Cost pump
        if 'Centrifugal' in Type:
            # Find pump factor
            F_Tdict = F_Tcentrifugal
            F_T = 1 # Assumption
            if Type == 'CentrifugalSingle':
                if 50 <= q <= 900 and 50 <= h <= 400:
                    F_T = F_Tdict['VSC3600']
                elif 50 <= q <= 3500 and 50 <= h <= 2000:
                    F_T = F_Tdict['VSC1800']
                elif 100 <= q <= 1500 and 100 <= h <= 450:
                    F_T = F_Tdict['HSC3600']
                elif 250 <= q <= 5000 and 50 <= h <= 500:
                    F_T = F_Tdict['HSC1800']
            elif Type == 'CentrifugalDouble':
                if 50 <= q <= 1100 and 300 <= h <= 1100:
                    F_T = F_Tdict['2HSC3600']
                elif 100 <= q <= 1500 and 650 <= h <= 3200:
                    F_T = F_Tdict['2+HSC3600']
            # Calculate cost
            self._S = S = self._calc_SizeFactor(q, h)
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
    def _calc_SizeFactor(q:'gal/min', h:'ft') -> 'S':
        return q*h**0.5
    
    @classmethod
    def _default_type(cls, Pi, Po, Qi):
        """Return default selection of pump type."""
        if Pi <= Po:
            for key, ranges in pump_ranges.items():
                flow_range, pressure_range = ranges
                flow_min, flow_max = flow_range
                pressure_min, pressure_max = pressure_range
                if flow_min <  Qi < flow_max and pressure_min < Pi and pressure_max > Po:
                    return key
        elif Pi > Po:
            flow_range, pressure_range = regenerative_ranges
            flow_min, flow_max = flow_range
            pressure_min, pressure_max = pressure_range
            if flow_min <  Qi < flow_max and pressure_min < Pi and pressure_max > Po:
                return 'Regenerative'
        raise Exception(f'No valid pump option for pressure, {Pi} Pa, and flow rate, {Qi} m3/hr.')
    
    @staticmethod
    def _available_PumpTypes(Pi, Po, Qi):
        """Return tuple of available pump types.
        
        **Parameters**
        
            Pi: [Stream] Input pressure (Pa)
            
            Po: [Stream] Output pressure (Pa)
            
            Qi: [Stream] Input flow rate (m^3/hr)
            
        """
        # Check which pumps work for flow rates
        pumps = [] # Pumps that work
        for key, ranges in pump_ranges.items():
            flow_range, pressure_range = ranges
            flow_min, flow_max = flow_range
            pressure_min, pressure_max = pressure_range
            if flow_min <  Qi < flow_max and pressure_min < Pi and pressure_max > Po:
                pumps.append(key)
        return pumps

    @staticmethod
    def _calc_PowerPressure(Pi, Po, Qi):
        """Return ideal power due to pressure change.
        
        **Parameters**
        
            Pi: [Stream] Input pressure
            
            Po: [Stream] Output pressure
            
            Qi: [Stream] Input flow rate
        
        """
        return Qi*(Po - Pi)
    
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
    def _calc_Head(p:'hp', mass:'kg/hr') -> '(ft)':
        return p/mass*897806 # 897806 = 1/g * unit_conversion_factor
    
    @staticmethod
    def _nearest_PumpPower(p:'hp'):
        for power in nema_sizes_hp:
            if power >= p:
                return power
        return power
        
    @staticmethod
    def _calc_Efficiency(q:'gpm', p:'hp'):
        mup = _calc_BreakEfficiency(q)
        mum = _calc_MotorEfficiency(p/mup)
        return mup*mum
        
        