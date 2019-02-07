# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 15:53:14 2018

@author: yoelr
"""
from biosteam import Unit
from biosteam import np
from fluids.pump import Corripio_motor_efficiency, \
    Corripio_pump_efficiency, motor_round_size, nema_sizes_hp

ln = np.log
exp = np.exp

max_hp = nema_sizes_hp[-1]

# Ranges of flow rate where pumps work (m3/hr)
Caps = {'Screw':             (0.2271, 13.626),
        'MeteringPlunger':   (0.2271, 6.813),
        'MeteringDiaphragm': (0.2271, 1.1355),
        'Regenerative':      (0.0100, 34.065), # Assume lower minimum for costing
        'CentrifugalSingle': (0.0100, 2271), # Assume lower minimum for costing
        'CentrifugalDouble': (0.2271, 2271),
        'DirectactingSteam': (0.2271, 454.2),
        'Gear':              (0.2271, 681.3),
        'AxialFlow':         (227.1, 2271)}

# Working pressures of said pumps
Pressures = {'Screw':             (68950, 5516000),
             'MeteringPlunger':   (68950, 20685000),
             'MeteringDiaphragm': (68950, 5516000),
             'Regenerative':      (1379000, 6895000),
             'CentrifugalSingle': (0, 689500), # Assume lower minimum for costing
             'CentrifugalDouble': (0, 6895000), # Assume lower minimum for costing
             'DirectactingSteam': (68950, 13790000),
             'Gear':              (68950, 13790000),
             'AxialFlow':         (68950, 137900)}

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
    _power_util = True
    kwargs = {'P': None}
    
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
        ins = self.ins
        outs = self.outs
        P = self.kwargs['P']
        for i, o in zip(ins, outs):
            o.mol = i.mol
            o.T = i.T
            o.phase = i.phase
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
        Qo = so.volnet
        mass = si.massnet
        
        # Get type
        Type = self.Type
        if Type == 'Default':
            Type = self._default_type(Pi, Po, Qi, Qo)
        Oper['Type'] = Type
        
        # Get ideal power
        if abs(Po - Pi) < 1:    
            Po = self.P_startup
        power_ideal = self._calc_PowerPressure(Pi, Po, Qi, Qo)*3.725e-7
        Oper['Ideal power'] = power_ideal # hp
        
        Oper['Flow rate'] = q = Qi*4.403
        if power_ideal <= max_hp:
            Oper['Efficiency'] = e = self._calc_Efficiency(q, power_ideal)
            Oper['Power'] = power =  self._nearest_PumpPower(power_ideal*e)
            Oper['N_pumps'] = 1
            Oper['Head'] = self._calc_Head(power, mass)
        else:
            raise Exception('Do this soon')
        
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
    def _default_type(cls, Pi, Po, Qi, Qo):
        """Return default selection of pump type."""
        options = cls._options(Pi, Po, Qi, Qo)
        Type = None
        if Pi <= Po:
            if 'CentrifugalSingle' in options:
                Type = 'CentrifugalSingle'
            elif 'CentrifugalDouble' in options:
                Type = 'CentrifugalDouble'
        elif Pi > Po:
            if 'Regenerative' in options:
                Type = 'Regenerative'
        else:
            if 'Gear' in options:
                Type = 'Gear'
        if not Type:
            try:
                Type = options[0]
            except IndexError:
                raise Exception(f'No valid pump option for pressure, {Pi} Pa, and flow rate, {Qi} m3/hr.')
        
        return Type
    
    @staticmethod
    def _options(Pi, Po, Qi, Qo):
        """Return tuple of pump type options.
        
        **Parameters**
        
            Pi: [Stream] Input pressure
            
            Po: [Stream] Output pressure
            
            Qi: [Stream] Input flow rate (m^3/hr)
            
            Qo: [Stream] Output flow rate (m^3/hr)
        
        """
        # Check which pumps work for flow rates
        pumpsV = [] # Pumps that work
        for key, range_ in Caps.items():
            min_, max_ = range_
            if min_ <  Qi < max_ and min_ < Qo  < max_:
                pumpsV.append(key)
                
        # Check which pumps work for pressures
        pumpsP = []
        for key, range_ in Pressures.items():
            min_, max_ = range_
            if min_ < Pi < max_ and min_ < Po < max_:
                pumpsP.append(key)
                
        # Pumps that work
        return tuple(set(pumpsV) & set(pumpsP))

    @staticmethod
    def _calc_PowerPressure(Pi, Po, Qi, Qo):
        """Return ideal power due to pressure change.
        
        **Parameters**
        
            Pi: [Stream] Input pressure
            
            Po: [Stream] Output pressure
            
            Qi: [Stream] Input flow rate
            
            Qo: [Stream] Output flow rate
        
        """
        return Qo*Po - Qi*Pi
    
    @staticmethod
    def _calc_PowerFlow(Qi, Qo, Dpipe, mass=None):
        """Return ideal power due to flow rate change.
        
        **Parameters**
            
            Qi: [Stream] Input flow rate
            
            Qo: [Stream] Output flow rate
            
            Dpipe: [float] Pipe inside diameter
            
            mass: [float] Mass flow rate
        
        """
        A = np.pi*((Dpipe/2)**2)
        
        if mass:
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
    def _calc_Efficiency(q:'gal/min', p:'hp'):
        if 50 < q < 5000 and 1 < p < 1500:
            mup = -0.316 + 0.24015*ln(q) - 0.01199*ln(q)**2
            mum = 0.8 + 0.0319*ln(p) - 0.00182*ln(p)**2
            e = 1/mup*mum
            if e > 1:
                e = 1
        else:
            e = 0.9
        return e
        
        