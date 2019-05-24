# -*- coding: utf-8 -*-
"""
Created on Tue Sep 18 10:28:29 2018

@author: Guest Group
"""
from scipy.optimize import least_squares, brentq, minimize_scalar, newton
from ._species import Species
from ._stream import Stream, nonzero_species, MassFlow, \
                            mol_flow_dim, mass_flow_dim, vol_flow_dim, \
                            phases, phase_index, flow
from ._utils import fraction, tuple_array, \
                           PropertyFactory, property_array
from ._exceptions import EquilibriumError, DimensionError
from . import _Q
import numpy as np

__all__ = ('MixedStream',)


# %% Dictionaries for finding the right attribute given the phase

# Used for the show method
letter_phases = {'s': 'solid',
                 'l': 'liquid',
                 'L': 'LIQUID',
                 'g': 'vapor'}

# %% Volumetric flow property

@PropertyFactory    
def VolumetricFlow(self):
    """Volumetric flow (m^3/hr)."""
    stream, mol, phase = self.data
    if mol:
        c = self.name # c: compound
        c.T = stream.T
        c.P = stream.P
        c.phase = phase
        return (c.Vm * mol[0] * 1000.)
    else:
        return 0.

@VolumetricFlow.setter
def VolumetricFlow(self, value):
    if value:
        stream, mol, phase = self.data
        c = self.name # c: compound
        c.T = stream.T
        c.P = stream.P
        c.phase = phase
        mol[0] = value/(c.Vm * 1000)
    else:
        mol[0] = 0.


# %% MixedStream (4 phase flow)

# TODO: Add solids handling:
# -Viscosity with class attribute J, the coefficient in Mc = Mw(1 - Vs/Vm)^(1/J).
# A generalized mixture rule for estimating the viscosity of solid-liquid suspensions and mechanical properties of polyphase rocks and composite materials, https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/2004JB003124, suggests that J = -0.5 for real suspensions including carbon fibers

# TODO: Make MixedStream capable of handling equilibrium, and mixing with Stream


class MixedStream(Stream):
    """An extension of the :doc:`Stream` class. All inherited properties and methods are supported. Create a MixedStream object that can manage mixtures of solid ('s'), liquid ('l'), LIQUID ('L'), and vapor ('g') phases. A MixedStream object is capable of vapor-liquid equilibrium (VLE) and liquid-LIQUID equilibrium (LLE).

    **Parameters**

         **ID:** [str] ID of the stream. If ID is '', a default ID will be chosen.

         **species:** [Species] Species object for the stream. If None, assume same as class' Species object.

         **T:** [float] Temperature (K)

         **P:** [float] Pressure (Pa)

    **Examples**

         :doc:`MixedStream objects and thermodynamic equilibrium`
    """
    __slots__ = ()
    def __init__(self, ID='', species=None, T=298.15, P=101325):
        # Get species and set species information
        if species is None: self._species = self._cls_species
        elif isinstance(species, Species):
            self._species = self._species
            species._read_only()
        else:
            raise TypeError(f"species must be a Species object, not '{type(species).__name__}'")
        self._T_QinP = self._T_VP = self._link = self._source = self._sink = self._ID = None
        self.ID = ID
        self.T = T  #: [float] Temperature (K)
        self.P = P  #: [float] Pressure (Pa)
        self._lL_split_cached = self._y_cached = self._bubble_cached = self._dew_cached = (None,)
        self._setflows(np.zeros((4, self._species._Nspecies), dtype=np.float64))

    def getflow(self, phase, *species, units='kmol/hr'):
        """Get flow rates of species in given units."""
        index = self.indices(*species) if species else ...
        p = phase_index[phase]
        if units == 'kmol/hr': return self._mol[p, index]
        q = _Q(1, units)
        dim = q.dimensionality
        if dim == mol_flow_dim:
            return self._mol[p, index]*q.to('kmol/hr').magnitude
        elif dim == mass_flow_dim:
            return self._mass[p, index]*q.to('kg/hr').magnitude
        elif dim == vol_flow_dim:
            return self._vol[p, index]*q.to('m3/hr').magnitude
        else:
            raise DimensionError(f"dimensions for flow units must be in molar, mass or volumetric flow rates, not '{dim}'")
        
    def setflow(self, phase, flow=(), species=(), units='kmol/hr',
                 inplace='', **flow_pairs):
        """Set flow rates according to the species order and flow_pairs. `inplace` can be any operation that can be performed in place (e.g. +, -, *, /, |, **, etc.)."""
        if species and isinstance(species[0], str):
            species = [self._species_dict[ID] for ID in species]
        species = (*species, *flow_pairs.keys())
        flow = (*flow, *flow_pairs.values())
        p = phase_index[phase]
        index = self.indices(*species) if species else ... 
        if units == 'kmol/hr':
            exec(f"self._mol[p, index] {inplace}= flow", locals())
            return
        q = _Q(flow, units)
        dim = q.dimensionality
        if dim == mol_flow_dim:
            exec(f"self._mol[p, index] {inplace}= q.to('kmol/hr').magnitude", locals())
        elif dim == mass_flow_dim:
            exec(f"self._mass[p, index] {inplace}= q.to('kg/hr').magnitude", locals())
        elif dim == vol_flow_dim:
            exec(f"self._vol[p, index] {inplace}= q.to('m3/hr').magnitude", locals())
        else:
            raise DimensionError(f"dimensions for flow units must be in molar, mass or volumetric flow rates, not '{dim}'")

    def _setflows(self, mol):
        """ Set molar, mass and volumetric flow rates for all phases.
        
        **Parameters**
        
            **mol_array:** Array of molar flow rates (kmol/hr) with species by column and phases ('s', 'l', 'L', 'g') by row.
        
        """
        species = self._species
        MW = species._MW
        num_compounds = species._num_compounds
        massflows = [] # Mass flow rates    
        volflows = [] # Volumetric flow rates    
        for phase, m in zip(phases, mol):
            vol = []
            mass = []
            for i, s in num_compounds:
                mi = m[i:i+1]
                mass.append(MassFlow(s.ID, (mi, MW[i])))
                vol.append(VolumetricFlow(s, (self, mi, phase.lower())))
            massflows.append(mass)
            volflows.append(vol)
        self._mol = mol
        self._mass = property_array(massflows)
        self._vol = property_array(volflows)

    @classmethod
    def proxy(cls, ID, link=None, T=298.15, P=101325):
        """Create a Stream object that serves as a proxy for the `link` stream."""
        self = object.__new__(cls)
        self._link = self._source = self._sink = self._ID = None
        self.T = T
        self.P = P
        self.ID = ID
        if link: self.link = link
        return self

    @property
    def phase(self):
        """String denoting all phases present"""
        phase = ''
        arr = self._mol != 0
        if arr[0].any(): phase += 's'
        if arr[1].any(): phase += 'l'
        if arr[2].any(): phase += 'L'
        if arr[3].any(): phase += 'g'
        return phase
    @phase.setter
    def phase(self, phase):
        if phase not in phases:
            raise ValueError("phase must be one of the following: 's', 'l', 'L' or 'g'")
        mol = self._mol
        val = mol.sum(0)
        mol[:] = 0
        mol[phase_index[phase]] = val

    ### Solids ###

    # Molar
    @flow
    def solid_mol(self):
        """solid phase molar flow rates (kmol/hr)."""
        return self._mol[0]
    @property
    def solid_molfrac(self):
        """solid phase molar fractions."""
        return fraction(self._mol[0]).view(tuple_array)
    @property
    def solid_molnet(self):
        """solid phase net molar flow rate (kg/hr)."""
        return self._mol[0].sum()

    # Mass
    @flow
    def solid_mass(self):
        """solid phase mass flow rates (kg/hr)."""
        return self._mass[0]
    @property
    def solid_massfrac(self):
        """solid phase mass fractions."""
        return fraction(self._mass[0]).view(tuple_array)
    @property
    def solid_massnet(self):
        """solid phase net mass flow rate (kg/hr)."""
        return self._mass[0].sum()

    # Volumetric
    @flow
    def solid_vol(self):
        """solid phase volumetric flow rates (m3/hr)."""
        return self._vol[0]
    @property
    def solid_volfrac(self):
        """solid phase volumetric fractions."""
        return fraction(self._vol[0]).view(tuple_array)
    @property
    def solid_volnet(self):
        """solid phase net volumetric flow rate (m3/hr)."""
        return self._vol[0].sum()

    ### Liquids ###

    # Molar
    @flow
    def liquid_mol(self):
        """liquid phase molar flow rates (kmol/hr)."""
        return self._mol[1]
    @property
    def liquid_molfrac(self):
        """liquid phase molar fractions."""
        return fraction(self._mol[1]).view(tuple_array)
    @property
    def liquid_molnet(self):
        """liquid phase net molar flow rate (kg/hr)."""
        return self._mol[1].sum()

    # Mass
    @flow
    def liquid_mass(self):
        """liquid phase mass flow rates (kg/hr)."""
        return self._mass[1]
    @property
    def liquid_massfrac(self):
        """liquid phase mass fractions."""
        return fraction(self._mass[1]).view(tuple_array)
    @property
    def liquid_massnet(self):
        """liquid phase net mass flow rate (kg/hr)."""
        return self._mass[1].sum()

    # Volumetric
    @flow
    def liquid_vol(self):
        """liquid phase volumetric flow rates (m3/hr)."""
        return self._vol[1]
    @property
    def liquid_volfrac(self):
        """liquid phase volumetric fractions."""
        return fraction(self._vol[1]).view(tuple_array)
    @property
    def liquid_volnet(self):
        """liquid phase net volumetric flow rate (m3/hr)."""
        return self._vol[1].sum()

    ### LIQUID ###

    # Molar
    @flow
    def LIQUID_mol(self):
        """LIQUID phase molar flow rates (kmol/hr)."""
        return self._mol[2]
    @property
    def LIQUID_molfrac(self):
        """LIQUID phase molar fractions."""
        return fraction(self._mol[2]).view(tuple_array)
    @property
    def LIQUID_molnet(self):
        """LIQUID phase net molar flow rate (kg/hr)."""
        return self._mol[2].sum()

    # Mass
    @flow
    def LIQUID_mass(self):
        """LIQUID phase mass flow rates (kg/hr)."""
        return self._mass[2]
    @property
    def LIQUID_massfrac(self):
        """LIQUID phase mass fractions."""
        return fraction(self._mass[2]).view(tuple_array)
    @property
    def LIQUID_massnet(self):
        """LIQUID phase net mass flow rate (kg/hr)."""
        return self._mass[2].sum()

    # Volumetric
    @flow
    def LIQUID_vol(self):
        """LIQUID phase volumetric flow rates (m3/hr)."""
        return self._vol[2]
    @property
    def LIQUID_volfrac(self):
        """LIQUID phase volumetric fractions."""
        return fraction(self._vol[2]).view(tuple_array)
    @property
    def LIQUID_volnet(self):
        """LIQUID phase net volumetric flow rate (m3/hr)."""
        return self._vol[2].sum()

    ### Vapor ###

    # Molar
    @flow
    def vapor_mol(self):
        """vapor phase molar flow rates (kmol/hr)."""
        return self._mol[3]
    @property
    def vapor_molfrac(self):
        """vapor phase molar fractions."""
        return fraction(self._mol[3]).view(tuple_array)
    @property
    def vapor_molnet(self):
        """vapor phase net molar flow rate (kg/hr)."""
        return self._mol[3].sum()

    # Mass
    @flow
    def vapor_mass(self):
        """vapor phase mass flow rates (kg/hr)."""
        return self._mass[3]
    @property
    def vapor_massfrac(self):
        """vapor phase mass fractions."""
        return fraction(self._mass[3]).view(tuple_array)
    @property
    def vapor_massnet(self):
        """vapor phase net mass flow rate (kg/hr)."""
        return self._mass[3].sum()

    # Volumetric
    @flow
    def vapor_vol(self):
        """vapor phase volumetric flow rates (m3/hr)."""
        return self._vol[3]
    @property
    def vapor_volfrac(self):
        """vapor phase volumetric fractions."""
        return fraction(self._vol[3]).view(tuple_array)
    @property
    def vapor_volnet(self):
        """vapor phase net volumetric flow rate (m3/hr)."""
        return self._vol[3].sum()
    
    ### Overall ###

    # Flow rates
    @property
    def mol(self):
        """Molar flow rates (kmol/hr)."""
        return self._mol.sum(0).view(tuple_array)
    @property
    def mass(self):
        """Mass flow rates (kmol/hr)."""
        return (self._mol.sum(0) * self._species._MW).view(tuple_array)
    @property
    def vol(self):
        """Volumetric flow rates (kmol/hr)."""
        return self._vol.sum(0).view(tuple_array)

    # Fractions
    @property
    def molfrac(self):
        """Molar fractions."""
        return fraction(self._mol.sum(0)).view(tuple_array)
    @property
    def massfrac(self):
        """Mass fractions."""
        return fraction(self._mol.sum(0) * self._species._MW).view(tuple_array)
    @property
    def volfrac(self):
        """Volumetric fractions."""
        return fraction(self._vol.sum(0)).view(tuple_array)

    # Net flow rates
    @property
    def molnet(self):
        """Mol net flow rate (kmol/hr)."""
        return self._mol.sum()
    @property
    def massnet(self):
        """Mass net flow rate (kmol/hr)."""
        return (self._mol.sum(0)*self._species._MW).sum()
    @property
    def volnet(self):
        """Volumetric net flow rate (kmol/hr)."""
        return self._vol.sum()

    # Other properties
    @property
    def mu(self):
        # Katti, P.K.; Chaudhri, M.M. (1964). "Viscosities of Binary Mixtures of Benzyl Acetate with Dioxane, Aniline, and m-Cresol". Journal of Chemical and Engineering Data. 9 (1964): 442–443.
        molnet = self.molnet
        if self.molnet == 0: return 0
        N = self._Nspecies
        mol = self._mol
        mus = np.zeros(4, N)
        Vms = np.zeros(4, N)
        props = self._species._props
        T = self._T
        P = self._P
        mus[0] = props('mu', mol[0], T, P, 's')
        Vms[0] = props('Vm', mol[0], T, P, 's')
        mus[1] = props('mu', mol[1], T, P, 'l')
        Vms[1] = props('Vm', mol[1], T, P, 'l')
        mus[2] = props('mu', mol[2], T, P, 'l')
        Vms[2] = props('Vm', mol[2], T, P, 'l')
        mus[3] = props('mu', mol[3], T, P, 'g')
        Vms[3] = props('Vm', mol[3], T, P, 'g')
        pos = np.where(mol != 0)
        mol = mol[pos]
        molfrac = mol/molnet
        return np.exp((molfrac*np.log(mus[pos]*Vms[pos])).sum())/self.Vm
    mu.__doc__ = Stream.mu.__doc__
    
    ### Energy flows ###

    @property
    def H(self):
        """Enthalpy flow rate as a function of T, excluding formation energies (kJ/hr).

        >>> s1.H # The stream is at the reference state
        0.0
        >>> s1.H = 1000
        >>> s1.T
        304.80600901692236
        >>> s1.H
        1002.1604769775776

        .. Note:: The solver reaches a level of tolerance and the new temperature is an approximation.
        """
        return self._species._mixedpropflow('H', self._mol, self.T, self.P)
    H = H.setter(Stream.H.fset)
        
    @property
    def S(self):
        """Entropy flow rate as a function of T and P, excluding formation energies (kJ/hr).

        >>> s1.S # The stream is at the reference state
        0.0
        >>> s1.T = 320
        >>> s1.S
        10.642980522582695
        """
        return self._species._mixedpropflow('S', self._mol, self.T, self.P)

    @property
    def C(self):
        """Heat capacity flow rate as a function of T (kJ/K/hr).

        >>> s2.C
        262.74816631655267
        """
        return self._species._mixedpropflow('Cpm', self._mol, self.T, self.P)

    @property
    def Vm(self):
        """Molar volume as a function of T and P (m^3/mol).

        >>> s1.Vm
        1.8069039870122814e-05
        """
        return self._species._mixedprop('Vm', self._mol, self.T, self.P)

    @property
    def k(self):
        """Thermal conductivity as a function of T and P (W/m/k).

        >>> s1.k
        0.5942044328004411
        """
        return self._species._mixedprop('k', self._mol, self.T, self.P)

    @property
    def alpha(self):
        """Thermal diffusivity as a function of T and P (m^2/s).

        >>> s1.alpha
        1.4255801521655763e-07
        """
        return self._species._mixedprop('alpha', self._mol, self.T, self.P)

    @property
    def sigma(self):
        """Surface tension as a function of T (N/m).

        >>> s1.sigma
        0.07205503890847455
        """
        return self._species._mixedprop('sigma', self._mol, self.T, self.P)

    @property
    def Pr(self):
        """Prandtl number as a function of T and P (non-dimensional).

        >>> s1.Pr
        6.421553249380879
        """
        return self._species._mixedprop('Pr', self._mol, self.T, self.P)
    
    ### Specifications ###

    def copylike(self, stream):
        """Copy mol, T, P, and phase of stream to self."""
        if self._species is not stream._species:
            raise ValueError('species must be the same to copy stream specifications')
        if isinstance(stream, MixedStream):
            self._mol[:] = stream._mol
        elif isinstance(stream, Stream):
            self.empty()
            self._mol[phase_index[stream._phase]] = stream.mol
        else:
            raise TypeError('argument must be a Stream object')
        self.P = stream.P
        self.T = stream.T

    ### Equilibrium ###

    def disable_phases(self, phase):
        """Cast stream into a Stream object.
        
        **Parameters**
        
            **phase:** {'s', 'l', 'g'} desired phase of stream
            
        """
        self.__class__ = Stream
        self.__init__(ID=self._ID, flow=self.mol, T=self.T, P=self.P, phase=phase)

    def enable_phases(self):
        """Cast stream into a MixedStream object."""

    def VLE(self, species_IDs=None, LNK=None, HNK=None, P=None,
            Qin=None, T=None, V=None, x=None, y=None):
        """Perform vapor-liquid equilibrium. Pressure defaults to current pressure.

        **Optional parameters**

            **P:** [float] Operating pressure (Pa)
                
            **species_IDs:** [tuple] IDs of equilibrium species
                 
            **LNK:** tuple[str] Light non-keys
    
            **HNK:** tuple[str] Heavy non-keys

            **User may define one:**
                * **Qin:** [float] Energy input (kJ/hr)
                * **T:** [float] Operating temperature (K)
                * **V:** [float] Overall molar vapor fraction
                * **x:** [numpy array] Molar composition of liquid (for binary mixture)
                * **y:** [numpy array] Molar composition of vapor (for binary mixture)     

        .. Note:
           LNK and HNK are not taken into account for equilibrium. Assumes constant pressure if no pressure is provided.

        """
        ### Decide what kind of equilibrium to run ###
        
        Nspecs = 0
        for i in (P, T, x, y, V, Qin):
            Nspecs += (i is not None)
            
        raise_error = False
        if Nspecs == 0:
            Qin = 0
        elif Nspecs == 2 and not P:
            raise_error = True
        elif Nspecs > 2:
            raise_error = True
        
        if T is not None:
            self.T = T
            eq = 'TP'
        elif x is not None:
            x = np.asarray(x)
            eq = 'xP'
        elif y is not None:
            y = np.asarray(y)
            eq = 'yP'
        elif V is not None:
            eq = 'VP'
        elif Qin is not None:
            eq = 'QinP'
            Hin = self.H + Qin
        else:
            raise_error = True
        
        if raise_error:
            raise ValueError("invalid specification. Must pass one and only one of the following specifications: T, x, y, V, or Qin")
        
        if P:
            self.P = P
        else:
            P = self.P
            
        ### Set Up Arguments ###
        
        # Reused attributes
        _species = self._species
        sp_dict = _species.__dict__
        sp_index = _species._IDs.index

        # Get flow rates
        liquid_mol = self._mol[1]
        vapor_mol = self._mol[3]
        all_mol = liquid_mol + vapor_mol

        # Set up indices for both equilibrium and non-equilibrium species
        if species_IDs is None:
            species, index = _species._equilibrium_species(all_mol)
            species = tuple(species)
        else:
            index = [sp_index(specie) for specie in species_IDs]
            species = tuple([sp_dict[ID] for ID in species_IDs])
        if LNK:
            LNK_index = [sp_index(specie) for specie in LNK]
        else:
            LNK, LNK_index = _species._light_species(all_mol)
        if HNK:
            HNK_index = [sp_index(specie) for specie in HNK]
        else:
            HNK, HNK_index = _species._heavy_species(all_mol)
        mol = all_mol[index]

        # Set light and heavy keys
        vapor_mol[LNK_index] = all_mol[LNK_index]
        vapor_mol[HNK_index] = 0
        liquid_mol[HNK_index] = all_mol[HNK_index]
        liquid_mol[LNK_index] = 0

        ### Single component equilibrium case ###
        
        Nspecies = len(species)
        if Nspecies == 1:
            # Equilibrium based on one specie
            s = species[0]
            
            if eq == 'TP':
                # Either liquid or gas
                Psat = s.VaporPressure(T)
                if P < Psat:
                    liquid_mol[index] = 0
                    vapor_mol[index] = mol
                else:
                    liquid_mol[index] = mol
                    vapor_mol[index] = 0
                return
            
            elif eq == 'VP':
                # Set vapor fraction
                self.T = s.Tsat(P)
                vapor_mol[index] = mol*V
                liquid_mol[index] = mol - vapor_mol[index]
                return
            
            elif eq == 'QinP': 
                # Set temperature in equilibrium
                self.T = s.Tsat(P)
                
                # Check if super heated vapor
                vapor_mol[index] = mol
                liquid_mol[index] = 0
                H_dew = self.H
                if Hin >= H_dew:
                    self.H = Hin
                    return
    
                # Check if subcooled liquid
                vapor_mol[index] = 0
                liquid_mol[index] = mol
                H_bubble = self.H
                if Hin <= H_bubble:
                    self.H = Hin
                    return
                
                # Adjust vapor fraction accordingly
                V = (Hin - H_bubble)/(H_dew - H_bubble)
                vapor_mol[index] = mol*V
                liquid_mol[index] = mol - vapor_mol[index]
                return
                
            else:
                ValueError("Invalid specification. Pure component equilibrium requires one and only one of the following specifications: T, V, or Qin.")
                
        ### Begin multi-component equilibrium ###
        
        # Get overall composition
        molnet = mol.sum()
        if molnet == 0:
            return  # No equilibrium necessary
        else:
            zf = mol/molnet
        
        min_ = np.zeros(Nspecies)

        ### Get Equilibrium ###

        if eq == 'xP' or eq == 'yP':
            # Get temperature based on phase composition
            if eq == 'xP':
                self.T, y = self._bubble_point(species, x, P)
            else:
                self.T, x = self._dew_point(species, y, P)
                
            if Nspecies > 2:
                raise ValueError(f'more than two components in equilibrium. Only binary component equilibrium can be solved for specification, {eq}')
            # Get flow rates based on lever rule
            split_frac = (zf[0]-x[0])/(y[0]-x[0])
            if split_frac > 1.0001 or split_frac < -0.0001:
                raise EquilibriumError('desired composition is not feasible')
            elif split_frac > 1:
                split_frac = 1
            elif split_frac < 0:
                split_frac = 0
            v_net = molnet * split_frac
            vapor_mol[index] = v_net*np.array([y[0], 1 - y[0]])
            liquid_mol[index] = mol - vapor_mol[index]
            return

        # Need to find these for remainding equilibrium types
        T_dew, x_dew = self._dew_point(species, zf, P)
        T_bubble, y_bubble = self._bubble_point(species, zf, P)

        # Guestimate vapor composition for 'VP', 'TP', 'QinP' equilibrium
        y = self._y_cached[1] if self._y_cached[0] == species else None

        def v_given_TPmol(v_guess, T, P, mol):
            """Return equilibrium vapor flow rates given T, P and molar flow rates."""
            Psat_P = np.array([s.VaporPressure(T) for s in species])/P
            actcoef = self.activity_coefficients
            def v_error(v):
                """Error function for constant T and P, where v represents vapor flow rates."""
                l = mol - v
                x = l/l.sum()
                y = v/v.sum()
                err = x*Psat_P * actcoef(species, x, T) - y
                return abs(err)
            
            return least_squares(v_error, v_guess,
                                 bounds=(min_, mol),
                                 ftol=0.001).x

        if eq == 'VP':
            def V_error(T):
                nonlocal v
                if T >= T_dew:
                    return V - 1 - (T - T_dew)
                elif T <= T_bubble:
                    return V + (T - T_bubble)
                v = v_given_TPmol(v, T, P, mol)
                return (V - sum(v)/molnet)
            
            # Guess vapor flow rates
            pos = np.where(mol <= 0)
            mol[pos] = 10**-9
            if y is None: y = V*zf + (1-V)*y_bubble
            v = y*V*molnet
            
            # Solve
            T_VP = self._T_VP or (1-V)*T_bubble + V*T_dew
            try:
                self._T_VP = self.T = newton(V_error, T_VP)
            except:
                self._T_VP = self.T = brentq(V_error, T_bubble, T_dew)
            
            v[pos] = 0
            vapor_mol[index] = v
            liquid_mol[index] = mol - vapor_mol[index]
            return
        
        elif eq == 'TP':
            if y is None:
                # Guess composition in the vapor is a weighted average of bubble/dew points
                f = (T - T_dew)/(T_bubble - T_dew)
                y = f*zf + (1-f)*y_bubble
            
            # Check if there is equilibrium
            if T > T_dew:
                vapor_mol[index] = v = mol
            elif T < T_bubble:
                vapor_mol[index] = v = min_
            else:
                # Guess overall vapor fraction
                V_guess = (T - T_bubble)/(T_dew - T_bubble)

                # Guess vapor flow rates
                v_guess = y * V_guess * mol

                # Solve
                vapor_mol[index] = v = v_given_TPmol(v_guess, T, P, mol)
                self._y = v/sum(v)
            liquid_mol[index] = mol - v
            return
        
        elif eq == 'QinP':
            # Check if super heated vapor
            vapor_mol[index] = mol
            liquid_mol[index] = 0
            self.T = T_dew
            H_dew = self.H
            if Hin >= 1.01*H_dew:
                self.H = Hin
                return

            # Check if subcooled liquid
            vapor_mol[index] = 0
            liquid_mol[index] = mol
            self.T = T_bubble
            H_bubble = self.H
            if Hin <= 0.99*H_bubble:
                self.H = Hin
                return

            def H_error_T(T):
                nonlocal v
                self.T = T
                vapor_mol[index] = v = v_given_TPmol(v, T, P, mol)
                liquid_mol[index] = mol - v
                return abs(Hin - (self.H))

            # Guess T, overall vapor fraction, and vapor flow rates
            T_QinP = self._T_QinP or T_bubble + (Hin - H_bubble)/(H_dew - H_bubble) * (T_dew - T_bubble)
            V_guess = (T_QinP - T_bubble)/(T_dew - T_bubble)
            if y is None:
                # Guess composition in the vapor is a weighted average of boiling points
                f = (T_QinP - T_dew)/(T_bubble - T_dew)
                y = f*zf + (1-f)*y_bubble
            v = y * V_guess * mol
            
            # Solve
            try:
                try:
                    self._T_QinP = self.T = newton(H_error_T, T_QinP)
                except:
                    T_QinP = T_bubble + (Hin - H_bubble)/(H_dew - H_bubble) * (T_dew - T_bubble)
                    V_guess = (T_QinP - T_bubble)/(T_dew - T_bubble)
                    f = (T_QinP - T_dew)/(T_bubble - T_dew)
                    y = f*zf + (1-f)*y_bubble
                    v = y * V_guess * mol
                    self._T_QinP = self.T = newton(H_error_T, T_QinP)
            except:
                self._T_QinP = self.T = minimize_scalar(H_error_T,
                                                        bounds=(T_bubble, T_dew),
                                                        method='bounded').x
            self._y_cached = (species, v/v.sum())

    def LLE(self, species_IDs=None, split=None, lNK=(), LNK=(),
            solvents=(), solvent_split=(),
            P=None, T=None, Qin=0):
        """Perform LIQUID-liquid equilibrium.

        **Optional Parameters**

            **species_IDs:** *tuple[str]* IDs of equilibrium species
            
            **split:** *tuple[float]* Initial guess split fractions of each specie to the 'liquid'

            **lNK:** *tuple[str]* Species assumed to completely remain in the 'liquid' phase.

            **LNK:** *tuple[str]* Species assumed to completely remain in the 'LIQUID' phase.

            **solvents:** *tuple[str]* Species corresponding to specified solvent_split

            **solvent_split:** *tuple[float]* Split fractions of each specie to the 'liquid' phase.                
            
            **Qin:** *[float]* Energy input (kJ/hr)
            
            **T:** *[float]* Operating temperature (K)
            
            **P:** *[float]* Operating pressure (Pa)    

        .. Note:
           lNK and LNK are not taken into account for equilibrium. Assumes constant pressure and temperatrue if none are provided.

        """
        ### Set Up ###

        # Get molar flow rates
        liquid_mol = self._liquid_mol
        LIQUID_mol = self._LIQUID_mol
        all_mol = liquid_mol + LIQUID_mol

        # Set up indices for both equilibrium and non-equilibrium species
        _species = self._species
        sp_dict = _species.__dict__
        sp_index = _species._IDs.index
        if species_IDs is None:
            species, index = _species._equilibrium_species(all_mol)
        else:
            species = [sp_dict[ID] for ID in species_IDs]
            index = [sp_index(ID) for ID in species_IDs]
        
        for s, i in zip(species, index):
            ID = s.ID
            if (ID in solvents) or (ID in lNK) or (ID in LNK):
                species.remove(s)
                index.remove(i)

        species = tuple(species)
        lNK_index = [sp_index(s) for s in lNK]
        LNK_index = [sp_index(s) for s in LNK]
        solvent_index = [sp_index(s) for s in solvents]

        # Get min and max splits
        Nspecies = len(species)
        min_ = np.zeros(Nspecies)

        # Set up activity coefficients
        solvent_species = tuple([sp_dict[ID] for ID in solvents])
        activity_species = species + solvent_species
        act_coef = self.activity_coefficients

        # Make sure splits are given as arrays
        if split is None:
            split = self._lL_split_cached[1] if self._lL_split_cached[0] == species else None
            if split is None:
                split_dipoles = [s.dipole for s in species]
                solvent_dipoles = [s.dipole for s in solvent_species]
                dipoles = split_dipoles + solvent_dipoles
                split_dipoles = np.array(split_dipoles)
                try:
                    split = split_dipoles/max(dipoles)
                except TypeError as TE:
                    missing_dipoles = []
                    for i, is_missing in enumerate(split_dipoles==None):
                        if is_missing:
                            missing_dipoles.append(activity_species[i])
                    raise EquilibriumError(f'cannot make estimate for split due to missing dipole values for species: {missing_dipoles}')
                split = np.array(split)
            solvent_split = np.array(solvent_split)

        # Equilibrium by constant temperature or with duty
        if T:
            self.T = T
        elif Qin:
            self.H += Qin
            T = self.T
        else:
            T = self.T
        
        if P:
            self.P = P
        else:
            P = self.P

        # Equilibrium species
        mol = all_mol[index]  # Solute species flow rates
        Kmol = all_mol[solvent_index]  # Solvent species flow rates
        Kguess = Kmol * solvent_split  # Solvent species flow rates to 'l' phase
        molKmol = np.concatenate((mol, Kmol))  # Both

        # Non-equilibrium species
        lNK_mol = all_mol[lNK_index]
        LNK_mol = all_mol[LNK_index]

        ### Solve ###
        
        # Error function
        def guess_error(l_guess):
            # Error function for constant T and P, where l represents 'l' liquid flow rates
            l_guess = np.concatenate((l_guess, Kguess))
            L_guess = molKmol - l_guess
            x1_guess = l_guess/sum(l_guess)
            x2_guess = L_guess/sum(L_guess)
            x1_gamma = act_coef(activity_species, x1_guess, T)
            x2_gamma = act_coef(activity_species, x2_guess, T)
            err = (x1_guess * x1_gamma - x2_guess * x2_gamma)
            return abs(err)[0:Nspecies]

        # Guestimates
        pos = mol == 0
        mol[pos] = 10**-20  # Give it a little to not hit bounds for least squares
        l_guess = split * mol  # Initial guess

        # Solve
        sol = least_squares(guess_error, l_guess, bounds=(min_, mol))
        l_guess = sol.x
        l_guess[pos] = 0
        split = l_guess/mol
        
        # Make sure two phases are given
        xl = l_guess/sum(l_guess)
        L_guess = mol-l_guess
        xL = L_guess/sum(L_guess)
        if (xl - xL < 0.1).all():
            raise EquilibriumError('could not solve equilibrium, please input better split guesses')

        # Set results
        self._lL_split_cached = (species, split)
        liquid_mol[index] = l_guess
        liquid_mol[lNK_index] = lNK_mol
        liquid_mol[LNK_index] = 0
        liquid_mol[solvent_index] = Kmol * solvent_split
        LIQUID_mol[index] = mol-l_guess
        LIQUID_mol[lNK_index] = 0
        LIQUID_mol[LNK_index] = LNK_mol
        LIQUID_mol[solvent_index] = Kmol - liquid_mol[solvent_index]

    def _info(self, T, P, flow, fraction):
        """Return string with all specifications."""
        units = self.units
        basic_info = self._info_header() + '\n'
        T_units, P_units, flow_units, fraction = [(i if i is not None else j) for i, j in
                                                  zip((T, P, flow, fraction), self.display_units)]
        phases = self.phase
        basic_info += self._info_phaseTP(phases, T_units, P_units)
        
        # Get flow attribute name
        flow_dim = _Q(0, flow_units).dimensionality
        if flow_dim == mol_flow_dim:
            flow_type = 'mol'
        elif flow_dim == mass_flow_dim:
            flow_type = 'mass'
        elif flow_dim == vol_flow_dim:
            flow_type = 'vol'
        else:
            raise DimensionError(f"dimensions for flow units must be in molar, mass or volumetric flow rates, not '{flow_dim}'")
        
        def flow_rate_skeleton(phase, flow_type, flow_units, fraction, maxlen):
            phase_full_name = letter_phases[phase]
            beginning = f' {phase_full_name}: '
            new_line_spaces = len(beginning) * ' '

            if fraction:
                flowname = flow_type + 'net'
                flownet = getattr(self, phase_full_name + '_' + flowname)
                end = ('\n' + new_line_spaces + 'net' +(maxlen-3)*' ' + '  ' +
                       f'{flownet:.3g} {flow_units}')
            else:
                end = ''

            return beginning, new_line_spaces, end

        # Length of species column
        nonzero, species = self.nonzero_species
        if len(nonzero) == 0:
            return basic_info + ' flow: 0'
        all_lengths = [len(sp) for sp in species]
        maxlen = max(all_lengths + [7])  # include length of the word 'species'

        # Set up flow rate information for all phases
        first_phase = True
        phases_flowrates_info = ''
        num_IDs = self._species._num_IDs
        for phase in phases:
            # Get flow rates of phase in nice structure
            flow_attr = letter_phases[phase] + '_' + flow_type
            if fraction: flow_attr += 'frac'
            flow = getattr(self, flow_attr)
            nonzero, species = nonzero_species(num_IDs, flow)
            if fraction:
                flow_nonzero = flow[nonzero]
            else:
                flow_nonzero = _Q(flow[nonzero], units[flow_type]).to(flow_units).magnitude
            
            # Get basic structure for phase flow rates
            beginning, new_line_spaces, end = flow_rate_skeleton(phase, flow_type, flow_units, fraction, maxlen)

            # Make up for soild and vapor phase having 5 letters
            # (both liquid and LIQUID have 6 letters)
            if phase in 'sg':
                beginning += ' '
                new_line_spaces += ' '
                if fraction: end = '\n ' + end[1:]

            # Set flow rates info
            flowrates = ''
            l = len(nonzero)
            lengths = [len(sp) for sp in species]
            for i in range(l-1):
                spaces = ' ' * (maxlen - lengths[i])
                flowrates += f'{species[i]} ' + spaces + \
                    f' {flow_nonzero[i]:.3g}\n' + new_line_spaces
            spaces = ' ' * (maxlen - lengths[l-1])
            flowrates += (f'{species[l-1]} ' + spaces
                          + f' {flow_nonzero[l-1]:.3g}')

            # Add header to flow rates
            if first_phase:
                spaces = ' ' * (maxlen - 7)
                beginning = (f'{new_line_spaces}species{spaces}  {flow_units}\n'
                             + beginning)

            # Put it together
            phases_flowrates_info += beginning + flowrates + end + '\n\n'
            first_phase = False
            
        return basic_info + phases_flowrates_info[:-2]


Stream.VLE.__doc__ = MixedStream.VLE.__doc__
Stream.LLE.__doc__ = MixedStream.LLE.__doc__
MixedStream.mu.__doc__ = Stream.mu.__doc__


