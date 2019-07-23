# -*- coding: utf-8 -*-
"""
Created on Tue Sep 18 10:28:29 2018

@author: Guest Group
"""
from scipy.optimize import least_squares
from ._equilibrium import VLEsolver, solve_v, BubblePoint, DewPoint, \
                          DortmundActivityCoefficients
from ._species import Species, WorkingSpecies
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
        if species is None:
            assert self._cls_species, 'must define Stream.species first'
            self._species = self._cls_species
        elif isinstance(species, Species):
            self._species = WorkingSpecies(species)
        elif isinstance(species, WorkingSpecies):
            self._species = species
        else:
            raise TypeError(f"species must be a Species object, not '{type(species).__name__}'")
        self._source_link = self
        self._link = self._source = self._sink = self._ID = None
        self.ID = ID
        self.T = T  #: [float] Temperature (K)
        self.P = P  #: [float] Pressure (Pa)
        self._lL_split_cached = (None,)
        self._gamma = gamma = DortmundActivityCoefficients()
        self._bp = BubblePoint(gamma)
        self._dp = DewPoint(gamma)
        self._VLEsolver = VLEsolver()
        self._setflows(np.zeros((4, self._species._N), dtype=np.float64))

    def getflow(self, phase, *species, units='kmol/hr'):
        """Get flow rates of species in given units."""
        index = self.indices(species) if species else ...
        p = phase_index[phase]
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
        index = self.indices(species) if species else ... 
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
        cmps = species._compounds
        index = species._index
        massflows = [] # Mass flow rates    
        volflows = [] # Volumetric flow rates    
        for phase, m in zip(phases, mol):
            vol = []
            mass = []
            for i in index:
                s = cmps[i]
                mi = m[i:i+1]
                mass.append(MassFlow(s.ID, (mi, MW[i])))
                vol.append(VolumetricFlow(s, (self, mi, phase.lower())))
            massflows.append(mass)
            volflows.append(vol)
        self._mol = mol
        self._mass = property_array(massflows)
        self._vol = property_array(volflows)

    @property
    def link(self):
        return None
    @link.setter
    def link(self, stream):
        raise TypeError(f'{type(self).__name__} objects do not support linking')

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
        # molnet = self.molnet
        # if self.molnet == 0: return 0
        # N = self._N
        # mol = self._mol
        # mus = np.zeros(4, N)
        # Vms = np.zeros(4, N)
        # props = self._species._props
        # T = self._T
        # P = self._P
        # mus[0] = props('mu', mol[0], T, P, 's')
        # Vms[0] = props('Vm', mol[0], T, P, 's')
        # mus[1] = props('mu', mol[1], T, P, 'l')
        # Vms[1] = props('Vm', mol[1], T, P, 'l')
        # mus[2] = props('mu', mol[2], T, P, 'l')
        # Vms[2] = props('Vm', mol[2], T, P, 'l')
        # mus[3] = props('mu', mol[3], T, P, 'g')
        # Vms[3] = props('Vm', mol[3], T, P, 'g')
        # pos = np.where(mol != 0)
        # mol = mol[pos]
        # molfrac = mol/molnet
        # return np.exp((molfrac*np.log(mus[pos]*Vms[pos])).sum())/self.Vm
        return self._species._mixedprop('mu', self._mol, self.T, self.P)
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
        """Copy flow rates, T, P, and phase of stream to self."""
        self.copyflow(stream)
        self.P = stream.P
        self.T = stream.T

    def copyflow(self, stream, species=None, remove=False, exclude=False):
        """Copy flow rates of stream to self.
        
        **Parameters**
        
            **stream:** [Stream] Flow rates will be copied from here.
            
            **species:** iterable[str] Species IDs. Defaults to all species.
            
            **remove:** [bool] If True, copied species will be removed from `stream`.
            
            **exclude:** [bool] If True, exclude `species` when copying.
            
        .. Note::
            
            Species flow rates that are not copied will be set to zero.
        
        """
        assert self._species is stream._species, ('species must be the same to '
                                                  'copy stream specifications')
        if isinstance(stream, MixedStream):
            if species is None:
                self._mol[:] = stream._mol
                if remove: stream._mol[:] = 0
            else:
                indices = self.indices(species)
                if exclude:
                    self._mol[:] = stream._mol
                    self._mol[indices] = 0
                    if remove:
                        mol = stream._mol
                        mol[:], mol[indices] = 0, mol[indices]
                else:
                    mol[:] = stream._mol[:, indices]
                    if remove: stream.mol[:, indices] = 0
        elif isinstance(stream, Stream):
            if species is None:
                self._mol[:] = 0
                self._mol[phase_index[stream._phase]] = stream.mol
                if remove: stream.mol[:] = 0
            else:
                indices = self.indices(species)
                if exclude:
                    p = phase_index[stream._phase]
                    self._mol[p, :] = stream._mol
                    self._mol[p, indices] = 0
                    if remove:
                        mol = stream._mol
                        mol[:], mol[indices] = 0, mol[indices]
                else:
                    self._mol[:] = 0
                    self._mol[phase_index[stream._phase], indices] = stream._mol[indices]
                    if remove: stream.mol[indices] = 0
        else:
            raise TypeError('argument must be a Stream object')
            
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

    def VLE(self, species_IDs=(), LNK=(), HNK=(), P=None,
            Q=None, T=None, V=None, x=None, y=None):
        """Perform vapor-liquid equilibrium.

        **Parameters**
        
            **Specify two:**
                * **P:** [float] Operating pressure (Pa)
                * **Q:** [float] Energy input (kJ/hr)
                * **T:** [float] Operating temperature (K)
                * **V:** [float] Molar vapor fraction
                * **x:** [numpy array] Molar composition of liquid (for binary mixture)
                * **y:** [numpy array] Molar composition of vapor (for binary mixture)

            **species_IDs:** [tuple] IDs of species in equilibrium.
                 
            **LNK:** tuple[str] Light non-keys that remain as a vapor.
    
            **HNK:** tuple[str] Heavy non-keys that remain as a liquid.

        .. Note:
           LNK and HNK are not taken into account for equilibrium.

        """
        ### Decide what kind of equilibrium to run ###
        
        Nspecs = ((P is not None) + (T is not None) + (x is not None)
                 +(y is not None) + (V is not None) + (Q is not None))
        assert Nspecs == 2, ("must pass two of the following specifications: "
                             "P, T, x, y, V, or Q")
        
        if P:
            self.P = P
            if T:
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
            elif Q is not None:
                eq = 'QP'
                Hin = self.H + Q
        elif T:
            self.T = T
            if x is not None:
                x = np.asarray(x)
                eq = 'xT'
            elif y is not None:
                y = np.asarray(y)
                eq = 'yT'
            elif V is not None:
                eq = 'VT'
            elif Q is not None:
                eq = 'QT'
                Hin = self.H + Q
        elif V is not None:
            if x is not None:
                x = np.asarray(x)
                eq = 'xV'
            elif y is not None:
                y = np.asarray(y)
                eq = 'yV'
            elif Q is not None:
                eq = 'VQ'
                raise NotImplementedError('specification V and Q not implemented yet')
        elif Q:
            if x is not None:
                x = np.asarray(x)
                eq = 'xQ'
            elif y is not None:
                y = np.asarray(y)
                eq = 'yQ'
            raise NotImplementedError('specification {eq[1]} and Q not implemented yet')
        else:
            raise ValueError("can only pass either 'x' or 'y' arguments, not both")
            
        if eq[0] in ('x', 'y') and len(species_IDs) != 2:
            raise ValueError("species_IDs must be of length 2 to specify '{eq[0]}'")    
        
        ### Set Up Arguments ###

        # Get flow rates
        liquid_mol = self._mol[1]
        vapor_mol = self._mol[3]
        all_mol = liquid_mol + vapor_mol
        notzero = all_mol > 0

        # Reused attributes
        _species = self._species
        sp_dict = _species.__dict__
        sp_index = _species._indexdct

        # Set up indices for both equilibrium and non-equilibrium species
        if species_IDs:
            index = [sp_index[specie] for specie in species_IDs]
            species = [sp_dict[ID] for ID in species_IDs]
        else:
            index = _species._equilibrium_indices(notzero)
            cmps = _species._compounds
            species = [cmps[i] for i in index]
        if LNK:
            LNK_index = [sp_index[i] for i in LNK]
        else:
            LNK_index = _species._light_indices(notzero)
        if HNK:
            HNK_index = [sp_index[i] for i in HNK]
        else:
            HNK_index = _species._heavy_indices(notzero)
        mol = all_mol[index]

        # Set light and heavy keys
        vapor_mol[LNK_index] = all_mol[LNK_index]
        vapor_mol[HNK_index] = 0
        liquid_mol[HNK_index] = all_mol[HNK_index]
        liquid_mol[LNK_index] = 0

        ### Single component equilibrium case ###
        
        N = len(species)
        if N == 1:
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
            
            elif eq == 'QP': 
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
            
            elif eq == 'VT':
                # Set vapor fraction
                self.P = s.VaporPressure(T)
                vapor_mol[index] = mol*V
                liquid_mol[index] = mol - vapor_mol[index]
                return
            
            elif eq == 'QT': 
                # Set Pressure in equilibrium
                self.P = s.VaporPressure(P)
                
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
                
        ### Begin multi-component equilibrium ###
        
        # Get overall composition
        molnet = mol.sum()
        if molnet == 0:
            return  # No equilibrium necessary
        else:
            zf = mol/molnet

        ### Get Equilibrium ###
        gamma = self._gamma
        gamma.species = species

        if eq[0] in ('x', 'y'):
            # Get temperature based on phase composition
            if eq[1] == 'P':
                if eq[0] == 'x':
                    self.T, y = self._bp.solve_Ty(x, P)
                else:
                    self.T, x = self._dp.solve_Tx(y, P)
            elif eq[1] == 'T':
                if eq[0] == 'x':
                    self.P, y = self._bp.solve_Py(x, T)
                else:
                    self.P, x = self._dp.solve_Px(y, T)
            
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

        # Bounderies for remainding equilibrium types
        if eq[1] == 'P':
            T_dew, x_dew = self._dp.solve_Tx(zf, P)
            T_bubble, y_bubble = self._bp.solve_Ty(zf, P)
        elif eq[1] == 'T':
            P_dew, x_dew = self._dp.solve_Px(zf, T)
            P_bubble, y_bubble = self._bp.solve_Py(zf, T)

        solver = self._VLEsolver
        zs = mol/molnet
        if eq == 'VP':
            if V == 1:
                self.T = T_dew
                vapor_mol[index] = mol
                liquid_mol[index] = 0
                return
            elif V == 0:
                self.T = T_bubble
                vapor_mol[index] = 0
                liquid_mol[index] = mol
                return
            else:
                y = y or V*zf + (1-V)*y_bubble
                v = y*V*molnet
                
                def V_error(T):
                    v[:] = solve_v(v, T, P, mol, molnet, zs, N, gamma)
                    return v.sum()/molnet
                
                # Guess
                self.T = solver('T', 'V', V_error, T_bubble, T_dew, 0, 1, V)
                vapor_mol[index] = v
                liquid_mol[index] = mol - v
                solver.P = self.P
                solver.Q = self.H/self.massnet
        elif eq == 'VT':
            if V == 1:
                self.P = P_dew
                vapor_mol[index] = mol
                liquid_mol[index] = 0
                return
            elif V == 0:
                self.P = P_bubble
                vapor_mol[index] = 0
                liquid_mol[index] = mol
                return
            else:
                def V_error(P):
                    v[:] = solve_v(v, T, P, mol, molnet, zs, N, gamma)
                    return v.sum()/molnet
                
                v = (V*zf + (1-V)*y_bubble)*V*molnet
                self.P = solver('P', 'V', V_error, P_bubble, P_dew, 0, 1, V)
                vapor_mol[index] = v
                liquid_mol[index] = mol - v
                solver.T = self.T
                solver.Q = self.H/self.massnet
        elif eq == 'TP':
            # Guess composition in the vapor is a
            # weighted average of bubble/dew points
            f = (T - T_dew)/(T_bubble - T_dew)
            y = f*zf + (1-f)*y_bubble
            
            # Check if there is equilibrium
            if T >= T_dew:
                vapor_mol[index] = mol
                liquid_mol[index] = 0
                return
            elif T <= T_bubble:
                vapor_mol[index] = 0
                liquid_mol[index] = mol
                return
            else:
                # Guess overall vapor fraction
                V_guess = solver.V or (T - T_bubble)/(T_dew - T_bubble)
                
                # Guess vapor flow rates
                v = y * V_guess * mol

                # Solve
                vapor_mol[index] = v = solve_v(v, T, P, mol, molnet, zs, N, gamma)
            liquid_mol[index] = mol - v
            solver.T = T
            solver.P = P
            solver.Q = self.H/self.massnet
            solver.V = v.sum()/molnet
        elif eq == 'QP':
            # Check if super heated vapor
            vapor_mol[index] = mol
            liquid_mol[index] = 0
            self.T = T_dew
            H_dew = self.H
            dH_dew = Hin - H_dew
            if dH_dew >= 0:
                self.H = Hin
                return

            # Check if subcooled liquid
            vapor_mol[index] = 0
            liquid_mol[index] = mol
            self.T = T_bubble
            H_bubble = self.H
            dH_bubble = Hin - H_bubble
            if dH_bubble <= 0:
                self.H = Hin
                return
            
            # Guess T, overall vapor fraction, and vapor flow rates
            V = dH_bubble/(H_dew - H_bubble)
            v = V*zf + (1-V)*y_bubble*V*molnet
            mass = self.massnet
            
            def Q_error(T):
                self.T = T
                v[:] = solve_v(v, T, P, mol, molnet, zs, N, gamma)
                vapor_mol[index] = v
                liquid_mol[index] = mol - v
                return self.H/mass
            
            self.T = solver('T', 'Q', Q_error, T_bubble, T_dew, 
                            H_bubble/mass, H_dew/mass, Hin/mass)
            solver.P = self.P
            solver.V = v.sum()/molnet

        elif eq == 'QT':
            # Check if super heated vapor
            vapor_mol[index] = mol
            liquid_mol[index] = 0
            self.P = P_dew
            H_dew = self.H
            dH_dew = (Hin - H_dew)
            if dH_dew >= 0:
                self.H = Hin
                return

            # Check if subcooled liquid
            vapor_mol[index] = 0
            liquid_mol[index] = mol
            self.P = P_bubble
            H_bubble = self.H
            dH_bubble = (Hin - H_bubble)
            if dH_bubble <= 0:
                self.H = Hin
                return

            # Guess overall vapor fraction, and vapor flow rates
            V = dH_bubble/(H_dew - H_bubble)
            # Guess composition in the vapor is a weighted average of boiling points
            v = V*zf + (1-V)*y_bubble*V*molnet
            
            # Solve
            mass = self.massnet
            def Q_error(P):
                self.P = P
                v[:] = solve_v(v, T, P, mol, molnet, zs, N, gamma)
                vapor_mol[index] = v
                liquid_mol[index] = mol - v
                return self.H/mass
            
            self.P = solver('P', 'Q', Q_error, P_bubble, P_dew,
                            H_bubble/mass, H_dew/mass, Hin/mass)
            solver.T = self.T
            solver.V = v.sum()/molnet

    def LLE(self, species_IDs=None, split=None, lNK=(), LNK=(),
            solvents=(), solvent_split=(),
            P=None, T=None, Q=0):
        """Perform LIQUID-liquid equilibrium.

        **Optional Parameters**

            **species_IDs:** *tuple[str]* IDs of equilibrium species
            
            **split:** *tuple[float]* Initial guess split fractions of each specie to the 'liquid'

            **lNK:** *tuple[str]* Species assumed to completely remain in the 'liquid' phase.

            **LNK:** *tuple[str]* Species assumed to completely remain in the 'LIQUID' phase.

            **solvents:** *tuple[str]* Species corresponding to specified solvent_split

            **solvent_split:** *tuple[float]* Split fractions of each specie to the 'liquid' phase.                
            
            **Q:** *[float]* Energy input (kJ/hr)
            
            **T:** *[float]* Operating temperature (K)
            
            **P:** *[float]* Operating pressure (Pa)    

        .. Note:
           lNK and LNK are not taken into account for equilibrium. Assumes constant pressure and temperatrue if none are provided.

        """
        ### Set Up ###

        # Get molar flow rates
        liquid_mol = self.liquid_mol
        LIQUID_mol = self.LIQUID_mol
        all_mol = liquid_mol + LIQUID_mol
        notzero = all_mol > 0

        # Reused attributes
        _species = self._species
        sp_dict = _species.__dict__
        sp_index = _species._indexdct

        # Set up indices for both equilibrium and non-equilibrium species
        if species_IDs:
            index = [sp_index[specie] for specie in species_IDs]
            species = [sp_dict[ID] for ID in species_IDs]
        else:
            index = _species._equilibrium_indices(notzero)
            cmps = _species._compounds
            species = [cmps[i] for i in index]
        
        for s, i in zip(species, index):
            ID = s.ID
            if (ID in solvents) or (ID in lNK) or (ID in LNK):
                species.remove(s)
                index.remove(i)

        lNK_index = [sp_index[s] for s in lNK]
        LNK_index = [sp_index[s] for s in LNK]
        solvent_index = [sp_index[s] for s in solvents]

        # Get min and max splits
        Nspecies = len(species)
        min_ = np.zeros(Nspecies)

        # Set up activity coefficients
        solvent_species = [sp_dict[ID] for ID in solvents]
        activity_species = species + solvent_species
        self._gamma.species = species
        gamma = self._gamma

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
        elif Q:
            self.H += Q
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
            x1_gamma = gamma(x1_guess, T)
            x2_gamma = gamma(x2_guess, T)
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
        if hasattr(self, '_mol'):
            nonzero, species = self.nonzero_species
        else:
            return basic_info + f' link: {self._link}'
        T_units, P_units, flow_units, fraction = [(i if i is not None else j) for i, j in
                                                  zip((T, P, flow, fraction), self.display_units)]
        phases = self.phase
        basic_info += self._info_phaseTP(phases, T_units, P_units)
        len_ = len(nonzero)
        if len_ == 0: return basic_info + ' flow: 0'
        
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
        all_lengths = [len(sp) for sp in species]
        maxlen = max(all_lengths + [7])  # include length of the word 'species'

        # Set up flow rate information for all phases
        first_phase = True
        phases_flowrates_info = ''
        for phase in phases:
            # Get flow rates of phase in nice structure
            flow_attr = letter_phases[phase] + '_' + flow_type
            if fraction: flow_attr += 'frac'
            flow = getattr(self, flow_attr)
            nonzero, species = nonzero_species(self._species, flow)
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

#         def v_given_TPmol(v_guess, T, P, mol):
#             # Return equilibrium vapor flow rates
#             # given T, P and molar flow rates.
#             Psat_P = np.array([s.VaporPressure(T) for s in species])/P
#             actcoef = self.activity_coefficients
#             def v_error(v):
#                 # Error function for constant T and P,
#                 # where v represents vapor flow rates.
#                 l = mol - v
#                 x = l/l.sum()
#                 return abs(x*Psat_P * actcoef(species, x, T) - v/v.sum())
            
#             return least_squares(v_error, v_guess,
#                                  bounds=(np.zeros(Nspecies), mol),
#                                  ftol=0.001).x

#         if eq == 'VP':
#             if V == 1:
#                 self.T = T_dew
#                 vapor_mol[index] = mol
#                 liquid_mol[index] = 0
#                 return
#             elif V == 0:
#                 self.T = T_bubble
#                 vapor_mol[index] = 0
#                 liquid_mol[index] = mol
#                 return
#             else:
#                 def V_error(T):
#                     nonlocal v
#                     if T > T_dew:
#                         return V - (1 + T - T_dew)
#                     elif T < T_bubble:
#                         return V - (T - T_bubble)
#                     v = v_given_TPmol(v, T, P, mol)
#                     return V - v.sum()/molnet
                
#                 # Guess vapor flow rates and solve
#                 T = self.T
#                 if (T_dew < T < T_dew):
#                     v = vapor_mol[index]/molnet
#                     self.T = golden_ratio(V_error, T_bubble-0.1, T_dew+0.1, 0.0001, T)
#                 else:
#                     y = V*zf + (1-V)*y_bubble
#                     v = y*V*molnet
#                     P = V*T_dew + (1-V)*T_bubble
#                     self.T = golden_ratio(V_error, T_bubble-0.1, T_dew+0.1, 0.0001, T)
                
#                 vapor_mol[index] = v
#                 liquid_mol[index] = mol - vapor_mol[index]
#         elif eq == 'VT':
#             if V == 1:
#                 self.P = P_dew
#                 vapor_mol[index] = mol
#                 liquid_mol[index] = 0
#                 return
#             elif V == 0:
#                 self.P = P_bubble
#                 vapor_mol[index] = 0
#                 liquid_mol[index] = mol
#                 return
#             else:
#                 def V_error(P):
#                     nonlocal v
#                     if P < P_dew:
#                         return V - (1 + P_dew - P)
#                     elif P > P_bubble:
#                         return V - (P_bubble - P)
#                     v = v_given_TPmol(v, T, P, mol)
#                     return V - v.sum()/molnet
                
#                 # Guess vapor flow rates and solve
#                 P = self.P
#                 if (P_dew < P < P_bubble):
#                     v = vapor_mol[index]/molnet
#                     self.P = golden_ratio(V_error, P_dew-10, P_bubble+10, 10, P)
#                 else:
#                     y = V*zf + (1-V)*y_bubble
#                     v = y*V*molnet
#                     P = V*P_dew + (1-V)*P_bubble
#                     self.P = golden_ratio(V_error, P_dew-10, P_bubble+10, 10, P)
                    
#                 vapor_mol[index] = v
#                 liquid_mol[index] = mol - vapor_mol[index]        
#         elif eq == 'TP':
#             y = self._y_cached[1] if self._y_cached[0] == species else None
#             if y is None:
#                 # Guess composition in the vapor is a
#                 # weighted average of bubble/dew points
#                 f = (T - T_dew)/(T_bubble - T_dew)
#                 y = f*zf + (1-f)*y_bubble
            
#             # Check if there is equilibrium
#             if T >= T_dew:
#                 vapor_mol[index] = mol
#                 liquid_mol[index] = 0
#                 return
#             elif T <= T_bubble:
#                 vapor_mol[index] = 0
#                 liquid_mol[index] = mol
#                 return
#             else:
#                 # Guess overall vapor fraction
#                 V_guess = (T - T_bubble)/(T_dew - T_bubble)
                
#                 # Guess vapor flow rates
#                 v_guess = y * V_guess * mol

#                 # Solve
#                 vapor_mol[index] = v = v_given_TPmol(v_guess, T, P, mol)
#             self._y_cached = (species, v/v.sum())
#             liquid_mol[index] = mol - v
#         elif eq == 'QP':
#             # Check if super heated vapor
#             vapor_mol[index] = mol
#             liquid_mol[index] = 0
#             self.T = T_dew
#             H_dew = self.H
#             if Hin >= H_dew:
#                 self.H = Hin
#                 return

#             # Check if subcooled liquid
#             vapor_mol[index] = 0
#             liquid_mol[index] = mol
#             self.T = T_bubble
#             H_bubble = self.H
#             if Hin <= H_bubble:
#                 self.H = Hin
#                 return

#             def H_error(T):
#                 nonlocal v
#                 self.T = T
#                 vapor_mol[index] = v = v_given_TPmol(v, T, P, mol)
#                 liquid_mol[index] = mol - v
#                 return Hin - self.H

#             # Guess T, overall vapor fraction, and vapor flow rates
#             T = self.T
#             if T_bubble < T < T_dew:
#                 v = vapor_mol[index]/molnet
#             else:
#                 T = T_bubble + (Hin - H_bubble)/(H_dew - H_bubble) * (T_dew - T_bubble)
#                 V_guess = (T - T_bubble)/(T_dew - T_bubble)
#                 # Guess composition in the vapor is a weighted average of boiling points
#                 f = (T - T_dew)/(T_bubble - T_dew)
#                 y = f*zf + (1-f)*y_bubble
#                 v = y * V_guess * mol
            
#             # Solve
#             v = y * V_guess * mol
#             self.T = golden_ratio(H_error, T_bubble, T_dew, 0.0001, T)
#         elif eq == 'QT':
#             # Check if super heated vapor
#             vapor_mol[index] = mol
#             liquid_mol[index] = 0
#             self.P = P_dew
#             H_dew = self.H
#             if Hin >= H_dew:
#                 self.H = Hin
#                 return

#             # Check if subcooled liquid
#             vapor_mol[index] = 0
#             liquid_mol[index] = mol
#             self.P = P_bubble
#             H_bubble = self.H
#             if Hin <= H_bubble:
#                 self.H = Hin
#                 return
 
#             def H_error(P):
#                 nonlocal v
#                 self.P = P
#                 vapor_mol[index] = v = v_given_TPmol(v, T, P, mol)
#                 liquid_mol[index] = mol - v
#                 return Hin - self.H

#             # Guess P, overall vapor fraction, and vapor flow rates
#             P = self.P
#             if P_dew < P < P_bubble:
#                 v = vapor_mol[index]/molnet
#             else:
#                 P = P_bubble + (Hin - H_bubble)/(H_dew - H_bubble) * (P_dew - P_bubble)
#                 V_guess = (P - P_bubble)/(P_dew - P_bubble)
#                 f = (P - P_dew)/(P_bubble - P_dew)
#                 y = f*zf + (1-f)*y_bubble
#                 v = y * V_guess * mol
            
#             # Solve
#             self.P = golden_ratio(H_error, P_dew, P_bubble, 10, P)