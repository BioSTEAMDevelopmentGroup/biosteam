# -*- coding: utf-8 -*-
"""
Created on Sat Aug 18 13:26:29 2018

@author: yoelr
"""
from scipy import integrate
from math import log
from .. import _Q
from .._utils import DisplayUnits

__all__ = ('Compound',)

_R = 8.3144598  # Universal gas constant (J/mol-K)

# %% Integrals methods

# Methods taking the integral of a function
def _integral_rigorous(func, Xi, Xf):
    """Take the full integral."""
    return integrate.quad(func, Xi, Xf)[0]

def _integral_initial(func, Xi, Xf):
    """Take the integral with a constant derivative at the initial value."""
    return func(Xf)*(Xf-Xi)

def _integral_final(func, Xi, Xf):
    """Take the integral with a constant derivative at the final value."""
    return func(Xf)*(Xf-Xi)

def _integral_average(func, Xi, Xf):
    """Take the integral with a constant derivative at the average value."""
    return func((Xi+Xf)/2)*(Xf-Xi)


# %% Compound class

class Compound:
    """Abstract class, Compound, for making objects that contain pure component thermodynamic and transport properties. The enthalpy property, H, ignores excess ethalpies and calculates based on latent heats and the heat capacity at a constant pressure of 101325 Pa. It is suitable for solids and liquids which are weak functions of pressure

    **Parameters**

         **IDs:** [str] A unique identification

         **T:** [float] Temperature (K)

         **P:** [float] Pressure (Pa)

         **phase:** [str] 'l'(iquid) or 'g'(as)           

    **Class Attributes**

         **phase_ref** = 'l': Reference phase

         **H_ref** = 0: Reference enthalpy (kJ/hr)

         **T_ref** = 298.15: Reference temperature (K)
         
         **P_ref** = 101325: Reference pressure (Pa)

    """
    #: [DisplayUnits] Units of measure for IPython display
    display_units = DisplayUnits(T='K', P='Pa')
    
    #: [dict] Units of measure for material properties (class attribute). 
    units = dict(MW='g/mol',
                 T='K',
                 P='Pa',
                 H='J/mol',
                 S='J/mol',
                 G='J/mol',
                 U='J/mol',
                 A='J/mol',
                 Hf='J/mol',
                 Vm='m^3/mol',
                 Cpm='J/mol',
                 Cp='J/g/K',
                 rho='kg/m^3',
                 rhom='mol/m^3',
                 nu='m^2/s',
                 mu='Pa*s',
                 sigma='N/m',
                 k='W/m/K',
                 alpha='m^2/s')
    
    phase_ref = 'l' # Reference phase
    P_ref = 101325  # Reference pressure (Pa)
    T_ref = 298.15  # Reference temperature (K)
    S_ref = 0       # Reference entropy (kJ/kmol)
    H_ref = 0       # Reference enthalpy (kJ/kmol)
    MW = 1          # Arbitrary molecular weight (g/mol)

    #: Method for taking the integral of Cp
    _integral = staticmethod(_integral_average)

    def quantity(self, prop_ID):
        """Return a material property as a Quantity object as described in the `pint package <https://pint.readthedocs.io/en/latest/>`__ 

        **Parameters**

             **prop_ID:** *[str]* name of the property (e.g. 'mol', 'H', 'k' ...)
             
        """
        return _Q(getattr(self, prop_ID), self.units[prop_ID])

    ### Abstract methods ###
    
    def HeatCapacitySolid(self, T): return self.Cpsm
    def HeatCapacityLiquid(self, T): return self.Cplm
    def HeatCapacityGas(self, T): return self.Cpgm
    def calc_H_excess(self, T, P): return 0    
    def calc_S_excess(self, T, P): return 0

    ### Entropy methods ###
    
    def _delSs_delT(self, T0, T):
        """Change in Entropy of a solid due to change in temperature at constant pressure."""
        return self.HeatCapacitySolid((T+T0)/2)*log(T/T0)

    def _delSl_delT(self, T0, T):
        """Change in Entropy of a liquid due to change in temperature at constant pressure."""
        return self.HeatCapacityLiquid((T+T0)/2)*log(T/T0)

    def _delSg_delT(self, T0, T):
        """Change in Entropy of a gas due to change in temperature at constant pressure."""
        return self.HeatCapacityGas((T+T0)/2)*log(T/T0)

    ### Energies ###

    @property
    def H(self):
        """Enthapy (kJ/kmol) disregarding pressure and assuming the specified phase."""
        phase_ref = self.phase_ref
        phase = self.phase

        # Perform enthalpy calculations at 101325 Pa
        if phase_ref == phase:
            if phase == 'l':
                H = self._integral(self.HeatCapacityLiquid, self.T_ref, self.T)
            elif phase == 'g':
                H = self._integral(self.HeatCapacityGas, self.T_ref, self.T)
            elif phase == 's':
                H = self._integral(self.HeatCapacitySolid, self.T_ref, self.T)
        elif phase_ref == 'l' and phase == 'g':
            H = self.H_int_l_T_ref_l_to_Tb + self.Hvap_Tbm + \
                self._integral(self.HeatCapacityGas, self.Tb, self.T)
        elif phase_ref == 'g' and phase == 'l':
            H = -self.H_int_Tb_to_T_ref_g - self.Hvap_Tbm + \
                self._integral(self.HeatCapacityLiquid, self.Tb, self.T)
        elif phase_ref == 's' and phase == 'l':
            H = self.H_int_T_ref_s_to_Tm + self.Hfusm + \
                self._integral(self.HeatCapacityLiquid, self.Tm, self.T)
        elif phase_ref == 'l' and phase == 's':
            H = -self.H_int_l_Tm_to_T_ref_l - self.Hfusm + \
                self._integral(self.HeatCapacitySolid, self.Tm, self.T)
        elif phase_ref == 's' and phase == 'g':
            H = self.H_int_T_ref_s_to_Tm + self.Hfusm + self.H_int_l_Tm_to_Tb + \
                self.Hvap_Tbm + self._integral(self.HeatCapacityGas, self.Tb, self.T)
        elif phase_ref == 'g' and phase == 's':
            H = -self.H_int_Tb_to_T_ref_g - self.Hvap_Tbm - self.H_int_l_Tm_to_Tb - \
                self.Hfusm + self._integral(self.HeatCapacitySolid, self.Tm, self.T)
        return self.H_ref + H

    @property
    def S(self):
        """Entropy (kJ/kmol) assuming the specified phase."""
        S = self.S_ref
        phase = self.phase
        phase_ref = self.phase_ref

        # Add Pressure Entropy
        if self.phase == 'g':
            S += -_R*log(self.P/self.P_ref)

        # Add Temperature Entropy
        if phase == phase_ref:
            if phase == 'l':
                S += self._delSl_delT(self.T_ref, self.T)
            elif phase == 'g':
                S += self._delSg_delT(self.T_ref, self.T)
            elif phase == 's':
                S += self._delSs_delT(self.T_ref, self.T)
        elif phase_ref == 'l' and phase == 'g':
            S += self.S_int_l_T_ref_l_to_Tb + self.Hvap_Tbm / \
                self.Tb + self._delSg_delT(self.Tb, self.T)
        elif phase_ref == 'g' and phase == 'l':
            S += - self.S_int_Tb_to_T_ref_g - self.Hvapm / \
                self.Tb + self._delSl_delT(self.Tb, self.T)
        elif phase_ref == 's' and phase == 'l':
            S += self.S_int_T_ref_s_to_Tm + self.Hfusm / \
                self.Tm + self._delSl_delT(self.Tm, self.T)
        elif self.phase_ref == 'l' and phase == 's':
            S += - self.S_int_l_Tm_to_T_ref_l - self.Hfusm / \
                self.Tm + self._delSs_delT(self.Tm, self.T)
        elif phase_ref == 's' and phase == 'g':
            S += self.S_int_T_ref_s_to_Tm + self.Hfusm/self.Tm + \
                self.S_int_l_Tm_to_Tb + self.Hvap_Tbm / \
                self.Tb + self._delSg_delT(self.Tb, self.T)
        elif phase_ref == 'g' and phase == 's':
            S += - self.S_int_Tb_to_T_ref_g - self.Hvap_Tbm/self.Tb - \
                self.S_int_l_Tm_to_Tb - self.Hfusm / \
                self.Tm + self._delSs_delT(self.Tm, self.T)
        else:
            raise Exception(f'Error in Compound object "{self.ID}" with phase "{phase}" and reference "{phase_ref}.')
        return S
    
    ### Class factory ###

    @classmethod
    def Factory(cls, ID, const_prop_molar_dct={}, state_dep_prop_molar_dct={}):
        """Make new child class with constant and temperature dependent properties.

        **Parameters**

             **ID:** *[str]* ID of new Compound object

             **const_prop_molar_dct:** *[dict]* constant properties to set to new object

             **state_dep_prop_molar_dct:** *[dict]* methods that will turn into Python properties for the new object

        """
        child = type(ID, (cls,), {})
        
        # Set contant attributes
        for attr, const in const_prop_molar_dct.items():
            setattr(child, attr, const)
        
        # Set state dependent property
        for attr, fun in state_dep_prop_molar_dct.items():            
            setattr(child, attr, property(fun))

        # Replace the heat capacity objects for constants
        Cpsm = const_prop_molar_dct.get('Cpsm')
        Cplm = const_prop_molar_dct.get('Cplm')
        Cpgm = const_prop_molar_dct.get('Cpgm')
        if Cpsm:
            setattr(child, 'HeatCapacitySolid', lambda self, T: Cpsm)
        if Cplm:
            setattr(child, 'HeatCapacityLiquid', lambda self, T: Cplm)
        if Cpgm:
            setattr(child, 'HeatCapacityGas', lambda self, T: Cpgm)

        # Replace the heat capacity objects for functions
        Cpsm = const_prop_molar_dct.get('Cpsm')
        Cplm = const_prop_molar_dct.get('Cplm')
        Cpgm = const_prop_molar_dct.get('Cpgm')
        if Cpsm:
            setattr(child, 'HeatCapacitySolid', property(Cpsm))
        if Cplm:
            setattr(child, 'HeatCapacityLiquid', property(Cplm))
        if Cpgm:
            setattr(child, 'HeatCapacityGas', property(Cpgm))
        
        # Replace the heat capacity objects if just one value
        Cpm = const_prop_molar_dct.get('Cpm')
        if Cpm:
            # Replace the heat capacity objects
            setattr(child, 'HeatCapacitySolid', lambda T: const)
            setattr(child, 'HeatCapacityLiquid', lambda T: const)
            setattr(child, 'HeatCapacityGas', lambda T: const)
        
        return child

    @classmethod
    def set_integral_type(cls, integral_type):
        """Change the enthalpy and entropy property calculator to different levels of rigour.

        **Parameters**

             int_type: [str] One of the following:
                      * 'rigorous': Compute integral by full integration
                      * 'average': Compute integral at an average function value
                      * 'constant': Compute integral using current function value
        """
        if integral_type == 'rigorous':
            cls._integral = staticmethod(_integral_rigorous)
        elif integral_type == 'average':
            cls._integral = staticmethod(_integral_average)
        elif integral_type == 'final':
            cls._integral = staticmethod(_integral_final)
        elif integral_type == 'initial':
            cls._integral = staticmethod(_integral_final)
        else:
            raise ValueError("Must pass one of the following types of Cp integration: 'rigorous', 'average', 'initial', or 'final'")
        cls._integral_type = integral_type

    # Representation
    def _info(self, T, P):
        """Return string with all specifications."""
        units = self.units
        T_units, P_units = self.display_units
        T_units = T or T_units
        P_units = P or P_units
        
        # First line
        info = f"{type(self).__name__}: {self.ID}\n"
        
        # Second line (thermo)
        T = _Q(self.T, units['T']).to(T_units).magnitude
        P = _Q(self.P, units['P']).to(P_units).magnitude
        info += f" phase: '{self.phase}', T: {T:.5g} {T_units}, P: {P:.6g} {P_units}"
        
        return info

    def __str__(self):
        return self.ID

    def show(self, T=None, P=None):
        """print information on self"""
        print(self._info(T, P))
    _ipython_display_ = show

    def __repr__(self):
        return f'<{type(self).__name__}: {self.ID}>'

