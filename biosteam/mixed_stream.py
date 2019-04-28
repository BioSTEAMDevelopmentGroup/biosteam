# -*- coding: utf-8 -*-
"""
Created on Tue Sep 18 10:28:29 2018

@author: Guest Group
"""
from scipy.optimize import least_squares, brentq, minimize_scalar, newton
from .species import Species
from .stream import Stream, nonzero_species, MassFlow, \
                            mol_flow_dim, mass_flow_dim, vol_flow_dim
from .utils import fraction, tuple_array, material_array, \
                           PropertyFactory, property_array
from .exceptions import EquilibriumError, DimensionError
from . import np, Q_

ln = np.log
exp = np.exp

__all__ = ('MixedStream',)


# %% Dictionaries for finding the right attribute given the phase

phases = ('s', 'l', 'L', 'g')
range4 = (0, 1, 2, 3)
phase_index = dict(zip(phases, range4))

# Used for the show method
letter_phases = {'s': 'solid',
                 'l': 'liquid',
                 'L': 'LIQUID',
                 'g': 'vapor'}

# Used during MixedStream __init__
molflows_ID = ('_solid_mol',
               '_liquid_mol',
               '_LIQUID_mol',
               '_vapor_mol')

volflows_ID = ('_solid_vol',
               '_liquid_vol',
               '_LIQUID_vol',
               '_vapor_vol')

massflows_ID = ('_solid_mass',
                '_liquid_mass',
                '_LIQUID_mass',
                '_vapor_mass')

phases_volflow = {'s': 'solid_vol',
                  'l': 'liquid_vol',
                  'L': 'LIQUID_vol',
                  'g': 'vapor_vol'}


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
        return (c.Vm * mol.item(0) * 1000)
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
    #: The default phase when flow rates are zero.
    _default_phase = 'l'
    
    #: Unit source
    _source = None
    
    #: Unit sink
    _sink = None
    
    def __init__(self, ID='', species=None, T=298.15, P=101325):
        # Get species and set species information
        if species is None:
            self._copy_species(self)
        elif isinstance(species, Species):
            self._set_species(species)
        else:
            raise TypeError(f"species must be a Species object, not '{type(species).__name__}'")
        
        self.ID = ID
        
        #: [float] Temperature (K)
        self.T = T
        
        #: [float] Pressure (Pa)
        self.P = P
        self._dew_cached = {}
        self._bubble_cached = {}
        self._y_cached = {}
        self._lL_split_cached = {}
        self._setflows(np.zeros((4, self._Nspecies), dtype=np.float64))
        

    def getflow(self, phase, *species, units='kmol/hr'):
        """Get flow rates of species in given units."""
        index = self.indices(*species) if species else ...
        p = phase_index[phase]
        if units == 'kmol/hr': return self._molarray[p, index]
        flow_wt_units = Q_(1, units)
        dim = flow_wt_units.dimensionality
        if dim == mol_flow_dim:
            return self._molarray[p, index]*flow_wt_units.to(self.units['mol']).magnitude
        elif dim == mass_flow_dim:
            return self._massarray[p, index]*flow_wt_units.to(self.units['mass']).magnitude
        elif dim == vol_flow_dim:
            return self._volarray[p, index]*flow_wt_units.to(self.units['vol']).magnitude
        else:
            raise DimensionError(f"Dimensions for flow units must be in molar, mass or volumetric flow rates, not '{dim}'.")
        
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
            exec(f"self._molarray[p, index] {inplace}= flow", locals())
            return
        flow_wt_units = Q_(flow, units)
        dim = flow_wt_units.dimensionality
        if dim == mol_flow_dim:
            exec(f"self._molarray[p, index] {inplace}= flow_wt_units.to(self.units['mol']).magnitude", locals())
        elif dim == mass_flow_dim:
            exec(f"self._massarray[p, index] {inplace}= flow_wt_units.to(mass_units).magnitude", locals())
        elif dim == vol_flow_dim:
            exec(f"self._volarray[p, index] {inplace}= flow_wt_units.to(vol_units).magnitude", locals())
        else:
            raise DimensionError(f"Dimensions for flow units must be in molar, mass or volumetric flow rates, not '{dim}'.")

    def _setflows(self, molarray):
        """ Set molar, mass and volumetric flow rates for all phases.
        
        **Parameters**
        
            **mol_array:** Array of molar flow rates (kmol/hr) with species by column and phases ('s', 'l', 'L', 'g') by row.
        
        """
        MW = self._MW
        num_compounds = self._num_compounds
        
        # Molar flow rates
        self._molarray = molarray
        self._phases_molflow = dict(zip(phases, molarray))
        molarray = molarray.view(material_array)
        for mol, ID in zip(molarray, molflows_ID):
            setattr(self, ID, mol)
        
        massflows = [] # Mass flow rates    
        volflows = [] # Volumetric flow rates    
        for phase, mol in zip(phases, molarray):
            vol = []
            mass = []
            for i, s in num_compounds:
                mol_i = mol[i:i+1]
                m = MassFlow(s.ID, (mol_i, MW[i]))
                v = VolumetricFlow(s, (self, mol_i, phase.lower()))
                mass.append(m)
                vol.append(v)
            massflows.append(mass)
            volflows.append(vol)
        
        self._massarray = property_array(massflows)
        self._volarray = property_array(volflows)
        for ID_mass, ID_vol, mass, vol in zip(massflows_ID, volflows_ID, self._massarray, self._volarray):
            setattr(self, ID_mass, mass)
            setattr(self, ID_vol, vol)

    @property
    def phase(self):
        """String denoting all phases present"""
        phase = ''
        if (self._solid_mol != 0).any(): phase += 's'
        if (self._liquid_mol != 0).any(): phase += 'l'
        if (self._LIQUID_mol != 0).any(): phase += 'L'
        if (self._vapor_mol != 0).any(): phase += 'g'
        if phase == '': phase = self._default_phase
        return phase

    @phase.setter
    def phase(self, phase):
        if phase not in phases:
            raise ValueError("Phase must be one of the following: 's', 'l', 'L' or 'g'")
        mol = self._molarray.sum(0)
        self._default_phase = phase
        self.empty()
        self._phases_molflow[phase][:] = mol

    ### Solids ###

    # Molar
    @property
    def solid_mol(self):
        """solid phase molar flow rates (kmol/hr)."""
        return self._solid_mol

    @solid_mol.setter
    def solid_mol(self, val):
        self._solid_mol[:] = val

    @property
    def solid_molfrac(self):
        """solid phase molar fractions."""
        return fraction(self._solid_mol).view(tuple_array)

    @property
    def solid_molnet(self):
        """solid phase net molar flow rate (kg/hr)."""
        return self._solid_mol.sum()

    # Mass
    @property
    def solid_mass(self):
        """solid phase mass flow rates (kg/hr)."""
        return self._solid_mass

    @solid_mass.setter
    def solid_mass(self, val):
        self._solid_mass[:] = val

    @property
    def solid_massfrac(self):
        """solid phase mass fractions."""
        return fraction(self._solid_mol * self._MW).view(tuple_array)

    @property
    def solid_massnet(self):
        """solid phase net mass flow rate (kg/hr)."""
        return self._solid_mass.sum()

    # Volumetric
    @property
    def solid_vol(self):
        """solid phase volumetric flow rates (m3/hr)."""
        return self._solid_vol

    @solid_vol.setter
    def solid_vol(self, val):
        self._solid_vol[:] = val

    @property
    def solid_volfrac(self):
        """solid phase volumetric fractions."""
        return fraction(self._solid_vol).view(tuple_array)

    @property
    def solid_volnet(self):
        """solid phase net volumetric flow rate (m3/hr)."""
        return self._solid_vol.sum()

    ### Liquids ###

    # Molar
    @property
    def liquid_mol(self):
        """liquid phase molar flow rates (kmol/hr)."""
        return self._liquid_mol

    @liquid_mol.setter
    def liquid_mol(self, val):
        self._liquid_mol[:] = val

    @property
    def liquid_molfrac(self):
        """liquid phase molar fractions."""
        return fraction(self._liquid_mol).view(tuple_array)

    @property
    def liquid_molnet(self):
        """liquid phase net molar flow rate (kg/hr)."""
        return self._liquid_mol.sum()

    # Mass
    @property
    def liquid_mass(self):
        """liquid phase mass flow rates (kg/hr)."""
        return self._liquid_mass

    @liquid_mass.setter
    def liquid_mass(self, val):
        self._liquid_mass[:] = val

    @property
    def liquid_massfrac(self):
        """liquid phase mass fractions."""
        return fraction(self._liquid_mass).view(tuple_array)

    @property
    def liquid_massnet(self):
        """liquid phase net mass flow rate (kg/hr)."""
        return self._liquid_mass.sum()

    # Volumetric
    @property
    def liquid_vol(self):
        """liquid phase volumetric flow rates (m3/hr)."""
        return self._liquid_vol

    @liquid_vol.setter
    def liquid_vol(self, val):
        self._liquid_vol[:] = val

    @property
    def liquid_volfrac(self):
        """liquid phase volumetric fractions."""
        return fraction(self._liquid_vol).view(tuple_array)

    @property
    def liquid_volnet(self):
        """liquid phase net volumetric flow rate (m3/hr)."""
        return self._liquid_vol.sum()

    ### LIQUID ###

    # Molar
    @property
    def LIQUID_mol(self):
        """LIQUID phase molar flow rates (kmol/hr)."""
        return self._LIQUID_mol

    @LIQUID_mol.setter
    def LIQUID_mol(self, val):
        self._LIQUID_mol[:] = val

    @property
    def LIQUID_molfrac(self):
        """LIQUID phase molar fractions."""
        return fraction(self._LIQUID_mol).view(tuple_array)

    @property
    def LIQUID_molnet(self):
        """LIQUID phase net molar flow rate (kg/hr)."""
        return self._LIQUID_mol.sum()

    # Mass
    @property
    def LIQUID_mass(self):
        """LIQUID phase mass flow rates (kg/hr)."""
        return self._LIQUID_mass

    @LIQUID_mass.setter
    def LIQUID_mass(self, val):
        self._LIQUID_mass[:] = val

    @property
    def LIQUID_massfrac(self):
        """LIQUID phase mass fractions."""
        return fraction(self._LIQUID_mol * self._MW).view(tuple_array)

    @property
    def LIQUID_massnet(self):
        """LIQUID phase net mass flow rate (kg/hr)."""
        return self._LIQUID_mass.sum()

    # Volumetric
    @property
    def LIQUID_vol(self):
        """LIQUID phase volumetric flow rates (m3/hr)."""
        return self._LIQUID_vol

    @LIQUID_vol.setter
    def LIQUID_vol(self, val):
        self._LIQUID_vol[:] = val

    @property
    def LIQUID_volfrac(self):
        """LIQUID phase volumetric fractions."""
        return fraction(self._LIQUID_vol).view(tuple_array)

    @property
    def LIQUID_volnet(self):
        """LIQUID phase net volumetric flow rate (m3/hr)."""
        return self._LIQUID_vol.sum()

    ### Vapor ###

    # Molar
    @property
    def vapor_mol(self):
        """vapor phase molar flow rates (kmol/hr)."""
        return self._vapor_mol

    @vapor_mol.setter
    def vapor_mol(self, val):
        self._vapor_mol[:] = val

    @property
    def vapor_molfrac(self):
        """vapor phase molar fractions."""
        return fraction(self._vapor_mol).view(tuple_array)

    @property
    def vapor_molnet(self):
        """vapor phase net molar flow rate (kg/hr)."""
        return self._vapor_mol.sum()

    # Mass
    @property
    def vapor_mass(self):
        """vapor phase mass flow rates (kg/hr)."""
        return self._vapor_mass

    @vapor_mass.setter
    def vapor_mass(self, mass):
        self._vapor_mass[:] = mass

    @property
    def vapor_massfrac(self):
        """vapor phase mass fractions."""
        return fraction(self._vapor_mol * self._MW).view(tuple_array)

    @property
    def vapor_massnet(self):
        """vapor phase net mass flow rate (kg/hr)."""
        return self._vapor_mass.sum()

    # Volumetric
    @property
    def vapor_vol(self):
        """vapor phase volumetric flow rates (m3/hr)."""
        return self._vapor_vol

    @vapor_vol.setter
    def vapor_vol(self, vol):
        self._vapor_vol[:] = vol

    @property
    def vapor_volfrac(self):
        """vapor phase volumetric fractions."""
        return fraction(self._vapor_vol).view(tuple_array)

    @property
    def vapor_volnet(self):
        """vapor phase net volumetric flow rate (m3/hr)."""
        return self._vapor_vol.sum()
    
    ### Overall ###

    # Flow rates
    @property
    def mol(self):
        """Molar flow rates (kmol/hr)."""
        return self._molarray.sum(0).view(tuple_array)

    @property
    def mass(self):
        """Mass flow rates (kmol/hr)."""
        return (self._molarray.sum(0) * self._MW).view(tuple_array)
    
    @property
    def vol(self):
        """Volumetric flow rates (kmol/hr)."""
        return self._volarray.sum(0).view(tuple_array)

    # Fractions
    @property
    def molfrac(self):
        """Molar fractions."""
        return fraction(self._molarray.sum(0)).view(tuple_array)

    @property
    def massfrac(self):
        """Mass fractions."""
        return fraction(self._molarray.sum(0) * self._MW).view(tuple_array)

    @property
    def volfrac(self):
        """Volumetric fractions."""
        return fraction(self._volarray.sum(0)).view(tuple_array)

    # Net flow rates
    @property
    def molnet(self):
        """Mol net flow rate (kmol/hr)."""
        return self._molarray.sum()

    @property
    def massnet(self):
        """Mass net flow rate (kmol/hr)."""
        return (self._molarray.sum(0)*self._MW).sum()

    @property
    def volnet(self):
        """Volumetric net flow rate (kmol/hr)."""
        return self.solid_volnet + self.liquid_volnet + self.LIQUID_volnet + self.vapor_volnet

    # Other properties
    @property
    def mu(self):
        # Katti, P.K.; Chaudhri, M.M. (1964). "Viscosities of Binary Mixtures of Benzyl Acetate with Dioxane, Aniline, and m-Cresol". Journal of Chemical and Engineering Data. 9 (1964): 442–443.
        molnet = self.molnet
        if self.molnet == 0:
            return 0
        N = self._Nspecies
        mol = self._molarray.flatten()
        mus = np.zeros(N*4)
        Vms = np.zeros(N*4)
        start = 0
        end = N
        for phase in phases:
            mus[start:end] = self._phaseprop_list('mu', phase)
            Vms[start:end] = self._phaseprop_list('Vm', phase)
            start += N
            end += N
        pos = np.where(mol != 0)
        mol = mol[pos]
        molfrac = mol/molnet
        return exp(sum(molfrac*ln(mus[pos]*Vms[pos])))/self.Vm
    mu.__doc__ = Stream.mu.__doc__
    
    ### Material Properties ###

    # General for a given phase
    def _phaseprop_list(self, prop_ID, phase):
        """Return component property lists for given phase."""
        getattr_ = getattr
        out = np.zeros(self._Nspecies)
        phase = phase.lower()
        P = self.P
        T = self.T
        ic = self._num_compounds.__iter__().__next__
        for m in self._phases_molflow[phase]:
            if m:
                i, c = ic()
                c.P = P
                c.T = T
                c.phase = phase
                out[i] = getattr_(c, prop_ID)
            else: ic()
        return out

    def _phaseprop_molar_flow(self, prop_ID, phase):
        """Return array of component properties * kmol/hr for given phase."""
        return self._phaseprop_list(prop_ID, phase) * self._phases_molflow[phase]

    def _phaseprop_molar_flownet(self, prop_ID, phase):
        """Return sum of component properties * kmol/hr for given phase."""
        return self._phaseprop_molar_flow(prop_ID, phase).sum()

    def _phaseprop_molar_mean(self, prop_ID, phase):
        """Return molar weighted average property for given phase."""
        ppfn = self._phaseprop_list(prop_ID, phase) * self._phases_molflow[phase]
        return ppfn / ppfn.sum()

    # General for all phases
    def _prop_molar_flow(self, prop_ID):
        """Return array of component properties * kmol/hr."""
        ppf = self._phaseprop_molar_flow
        return ppf(prop_ID, 's') + ppf(prop_ID, 'l') + ppf(prop_ID, 'L') + ppf(prop_ID, 'g')

    def _prop_molar_flownet(self, prop_ID):
        """Return sum of component properties * kmol/hr"""
        return sum(self._prop_molar_flow(prop_ID))

    def _prop_molar_mean(self, prop_ID):
        """Return molar weighted average property"""
        return self._prop_molar_flownet(prop_ID) / self.molnet

    ### Specifications ###

    def copylike(self, stream):
        """Copy mol, T, P, and phase of stream to self."""
        if self._species is not stream._species:
            raise ValueError('species must be the same to copy stream specifications.')
        if isinstance(stream, MixedStream):
            self._molarray[:] = stream._molarray
        elif isinstance(stream, Stream):
            self.empty()
            self._phases_molflow[stream.phase][:] = stream.mol
        else:
            raise TypeError('Must pass a valid Stream instance')
        self.P = stream.P
        self.T = stream.T

    ### Equilibrium ###

    def disable_phases(self, phase):
        """Cast stream into a Stream object.
        
        **Parameters**
        
            **phase:** {'s', 'l', 'g'} desired phase of stream
            
        """
        mol = self.mol
        self.__class__ = Stream
        self.__init__(ID=self.ID, flow=mol, T=self.T, P=self.P, phase=phase)

    def enable_phases(self):
        """Cast stream into a MixedStream object."""

    # Vapor-liquid equilibrium
    def VLE(self, species_IDs=None, LNK=None, HNK=None, P=None,
            Qin=None, T=None, V=None, x=None, y=None):
        """Partition flow rates into vapor and liquid phases by equilibrium. Pressure defaults to current pressure.

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
            raise ValueError("Invalid specification. Must pass one and only one of the following specifications: T, x, y, V, or Qin.")
        
        if P:
            self.P = P
        else:
            P = self.P
            
        ### Set Up Arguments ###
        
        # Reused attributes
        sp_dict = self._species_dict 
        sp_index = self._IDs.index

        # Set up indices for both equilibrium and non-equilibrium species
        if species_IDs is None:
            species, index = self._equilibrium_species()
            species = tuple(species)
        else:
            index = [sp_index(specie) for specie in species_IDs]
            species = tuple(sp_dict[ID] for ID in species_IDs)
        if LNK:
            LNK_index = [sp_index(specie) for specie in LNK]
        else:
            LNK, LNK_index = self._light_species()
        if HNK:
            HNK_index = [sp_index(specie) for specie in HNK]
        else:
            HNK, HNK_index = self._heavy_species()

        # Get flow rates
        liquid_mol = self._liquid_mol
        vapor_mol = self._vapor_mol
        all_mol = liquid_mol + vapor_mol
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
                raise ValueError(f'More than two components in equilibrium. Only binary component equilibrium can be solved for specification, {eq}.')
            # Get flow rates based on lever rule
            split_frac = (zf[0]-x[0])/(y[0]-x[0])
            if split_frac > 1.0001 or split_frac < -0.0001:
                raise EquilibriumError('Desired composition is not feasible')
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
        y = self._y_cached.get(species)

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
            if not hasattr(self, '_T_VP'):
                self._T_VP = (1-V)*T_bubble + V*T_dew
            try:
                self._T_VP = self.T = newton(V_error, self._T_VP)
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
            if not hasattr(self, '_T_QinP'):
                self._T_QinP = T_bubble + (Hin - H_bubble)/(H_dew - H_bubble) * (T_dew - T_bubble)
            V_guess = (self._T_QinP - T_bubble)/(T_dew - T_bubble)
            if y is None:
                # Guess composition in the vapor is a weighted average of boiling points
                f = (self._T_QinP - T_dew)/(T_bubble - T_dew)
                y = f*zf + (1-f)*y_bubble
            v = y * V_guess * mol
            
            # Solve
            try:
                try:
                    self._T_QinP = self.T = newton(H_error_T, self._T_QinP)
                except:
                    self._T_QinP = T_bubble + (Hin - H_bubble)/(H_dew - H_bubble) * (T_dew - T_bubble)
                    V_guess = (self._T_QinP - T_bubble)/(T_dew - T_bubble)
                    f = (self._T_QinP - T_dew)/(T_bubble - T_dew)
                    y = f*zf + (1-f)*y_bubble
                    v = y * V_guess * mol
                    self._T_QinP = self.T = newton(H_error_T, self._T_QinP)
            except:
                self._T_QinP = self.T = minimize_scalar(H_error_T,
                                                        bounds=(T_bubble, T_dew),
                                                        method='bounded').x
                
            self._y_cached[species] = v/v.sum()

    # LIQUID-liquid equilibrium
    def LLE(self, species_IDs=None, split=None, lNK=(), LNK=(),
            solvents=(), solvent_split=(),
            P=None, T=None, Qin=0):
        """Partition flow rates into liquid and LIQUID phases by equilibrium.

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

        # Set up indices for both equilibrium and non-equilibrium species
        sp_dict = self._species_dict
        sp_index = self._IDs.index
        if species_IDs is None:
            species, index = self._equilibrium_species()
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
        solvent_species = tuple(sp_dict[ID] for ID in solvents)
        activity_species = species + solvent_species
        act_coef = self.activity_coefficients

        # Make sure splits are given as arrays
        if split is None:
            split = self._lL_split_cached.get(species)
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
                    raise EquilibriumError(f'Cannot make estimate for split. Missing dipole values for species: {missing_dipoles}')
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

        # Get molar flow rates
        liquid_mol = self._liquid_mol
        LIQUID_mol = self._LIQUID_mol
        all_mol = liquid_mol + LIQUID_mol

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
            raise EquilibriumError('Could not solve equilibrium, please input better split guesses')

        # Set results
        self._lL_split_cached[species] = split
        liquid_mol[index] = l_guess
        liquid_mol[lNK_index] = lNK_mol
        liquid_mol[LNK_index] = 0
        liquid_mol[solvent_index] = Kmol * solvent_split
        LIQUID_mol[index] = mol-l_guess
        LIQUID_mol[lNK_index] = 0
        LIQUID_mol[LNK_index] = LNK_mol
        LIQUID_mol[solvent_index] = Kmol - liquid_mol[solvent_index]

    def _info(self, **show_units):
        """Return string with all specifications."""
        units = self.units
        basic_info = self._info_header() + '\n'
        T_units, P_units, flow_units, fraction = self._info_units(show_units)
        phases = self.phase
        basic_info += self._info_phaseTP(phases, T_units, P_units)
        
        # Get flow attribute name
        flow_dim = Q_(0, flow_units).dimensionality
        if flow_dim == mol_flow_dim:
            flow_type = 'mol'
        elif flow_dim == mass_flow_dim:
            flow_type = 'mass'
        elif flow_dim == vol_flow_dim:
            flow_type = 'vol'
        else:
            raise DimensionError(f"Dimensions for flow units must be in molar, mass or volumetric flow rates, not '{flow_dim}'.")
        
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
        for phase in phases:
            # Get flow rates of phase in nice structure
            flow_attr = letter_phases[phase] + '_' + flow_type
            if fraction: flow_attr += 'frac'
            flow = getattr(self, flow_attr)
            nonzero, species = nonzero_species(self._num_IDs, flow)
            if fraction:
                flow_nonzero = flow[nonzero]
            else:
                flow_nonzero = Q_(flow[nonzero], units[flow_type]).to(flow_units).magnitude
            
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


