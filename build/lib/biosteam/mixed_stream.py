# -*- coding: utf-8 -*-
"""
Created on Tue Sep 18 10:28:29 2018

@author: Guest Group
"""
from scipy.optimize import least_squares, brentq, minimize_scalar
from biosteam.species import Species
from biosteam.stream import Stream, nonzero_species, MassFlow, vol_units, \
                            mol_flow_dim, mass_flow_dim, vol_flow_dim
from biosteam.utils import get_frac, tuple_array, material_array, \
                           PropertyFactory, property_array, reorder, CS
from biosteam.exceptions import EquilibriumError, DimensionError
from bookkeep import Q_
from biosteam import np
import copy

ln = np.log
exp = np.exp


# %% Dictionaries for finding the right attribute given the phase

phases = ('s', 'l', 'L', 'g')
range4 = (0, 1, 2, 3)

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

# Used throughout
phases_molflow = {'s': 'solid_mol',
                  'l': 'liquid_mol',
                  'L': 'LIQUID_mol',
                  'g': 'vapor_mol'}

phases_massflow = {'s': 'solid_mass',
                   'l': 'liquid_mass',
                   'L': 'LIQUID_mass',
                   'g': 'vapor_mass'}

phases_volflow = {'s': 'solid_vol',
                  'l': 'liquid_vol',
                  'L': 'LIQUID_vol',
                  'g': 'vapor_vol'}

phases_molfrac = {'s': 'solid_molfrac',
                  'l': 'liquid_molfrac',
                  'L': 'LIQUID_molfrac',
                  'g': 'vapor_molfrac'}

phases_massfrac = {'s': 'solid_massfrac',
                   'l': 'liquid_massfrac',
                   'L': 'LIQUID_massfrac',
                   'g': 'vapor_massfrac'}

phases_volfrac = {'s': 'solid_volfrac',
                  'l': 'liquid_volfrac',
                  'L': 'LIQUID_volfrac',
                  'g': 'vapor_volfrac'}


# %% Volumetric flow property

@PropertyFactory    
def VolumetricFlow(self):
    """Volumetric flow (m^3/hr)."""
    specie = self.name
    stream, mol, phase = self.data
    if mol:
        specie.T = stream.T
        specie.P = stream.P
        specie.phase = phase
        return float(getattr(specie, 'Vm') * mol[0] * 1000) 
    else:
        return 0.

@VolumetricFlow.setter
def VolumetricFlow(self, value):
    specie = self.name
    stream, mol, phase = self.data
    if value:
        specie.T = stream.T
        specie.P = stream.P
        specie.phase = phase
        mol[0] = value/(getattr(specie, 'Vm') * 1000)
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

         **solid_flow:** [array_like] solid flow rates by specie

         **liquid_flow:** [array_like] liquid flow rates by specie

         **LIQUID_flow:** [array_like] LIQUID flow rates by specie

         **vapor_flow:** [array_like] vapor flow rates by specie

         **species:** tuple[str] or [Species] Species corresponding to `flow` parameters. If empty, assume same species as class.

         **units:** [str] The units of the flow rates (only mass, molar, and volumetric flow rates are valid)

         **T:** [float] Temperature (K)

         **P:** [float] Pressure (Pa)

    **Examples**

         :doc:`MixedStream Example`
    """
    #: Species in LLE
    _species_eq = None
    
    #: The default phase when flow rates are zero.
    _default_phase = 'l'
    
    #: tuple with source unit ID and outs position
    _source = (None, None)
    
    #: tuple with sink unit ID and ins position
    _sink = (None, None)
    
    def __init__(self, ID='', solid_flow=(), liquid_flow=(), LIQUID_flow=(), vapor_flow=(), species=(), units='kmol/hr', T=298.15, P=101325):
        # Get species ID and set species information
        if isinstance(species, Species):
            # Set species as provided
            self._set_species(species)
            specie_IDs = ()
        elif hasattr(species, '__iter__'):
            specie_IDs = species
            species = self._species
            # Copy class species data (in case it changes)
            self._copy_species(self)
        else:
            raise TypeError(f'species must be a Species object, not {type(species).__name__}')
        
        self.ID = ID
        self.T = T
        self.P = P
        self._dew_cached = {}
        self._bubble_cached = {}
        
        # Match species and flow rates
        l_species = len(specie_IDs)
        flows = [solid_flow, liquid_flow, LIQUID_flow, vapor_flow]
        for i in range4:
            flow = flows[i]
            l_flow = len(flow)
            if l_flow == 0:
                flows[i] = np.zeros(self._Nspecies)
            elif l_species == l_flow:
                flows[i] = reorder(flow, specie_IDs, self._specie_IDs)
            elif l_flow != self._Nspecies:
                raise ValueError('length of flow rates must be equal to length of species')

        # Set flow rates to the right units
        mol_array = np.array(flows, dtype=np.float64)
        flow_wt_units = Q_(mol_array, units)
        dim = flow_wt_units.dimensionality
        if dim == mol_flow_dim:
            mol_array = flow_wt_units.to(self.units['mol']).magnitude
            self._set_flows(mol_array)
        elif dim == mass_flow_dim:
            mol_array = flow_wt_units.to(self.units['mass']).magnitude / self._MW
            self._set_flows(mol_array)
        elif dim == vol_flow_dim:
            self._set_flows(mol_array)
            for i, phase in zip(range4, phases):
                volflow = phases_volflow[phase]
                setattr(self, volflow, mol_array[i, :])
        else:
            raise DimensionError(f"Dimensions for flow rate units must be in molar, mass or volumetric flow rates, not '{dim}'.")

    def _set_flows(self, mol_array):
        """ Set molar, mass and volumetric flow rates for all phases.
        
        **Parameters**
        
            **mol_array:** Array of molar flow rates (kmol/hr) with species by column and phases ('s', 'l', 'L', 'g') by row.
        
        """
        species = self._species
        MW = self._MW
        
        # Molar flow rates
        self._molflows = molflows = []
        mol_array = mol_array.view(material_array)
        for mol, ID in zip(mol_array, molflows_ID):
            setattr(self, ID, mol)
            molflows.append(mol)
            
        # Mass flow rates    
        self._massflows = massflows = []
        for mol, ID in zip(molflows, massflows_ID):
            mass = []
            for i, spID in enumerate(species.ID):
                m = MassFlow(spID, (mol[i:i+1], MW[i]))
                mass.append(m)
            mass = property_array(mass)
            setattr(self, ID, mass)
            massflows.append(mass)
        
        # Volumetric flow rates    
        self._volflows = volflows = []   
        for phase, mol, ID in zip(phases, molflows, volflows_ID):
            vol = []
            for i, specie in enumerate(species):
                v = VolumetricFlow(specie, (self, mol[i:i+1], phase.lower()))
                vol.append(v)
            vol = property_array(vol)
            setattr(self, ID, vol)
            volflows.append(vol)

    @property
    def phase(self):
        """String denoting all phases present"""
        phase = ''
        if any(self._solid_mol != 0): phase += 's'
        if any(self._liquid_mol != 0): phase += 'l'
        if any(self._LIQUID_mol != 0): phase += 'L'
        if any(self._vapor_mol != 0): phase += 'g'
        if phase == '': phase = self._default_phase
        return phase

    @phase.setter
    def phase(self, phase):
        if phase not in phases:
            raise ValueError("Phase must be one of the following: 's', 'l', 'L' or 'g'")
        mol = sum(self._molflows)
        self._default_phase = phase
        self.empty()
        setattr(self, phases_molflow[phase], mol)

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
        return get_frac(self._solid_mol).view(tuple_array)

    @property
    def solid_molnet(self):
        """solid phase net molar flow rate (kg/hr)."""
        return sum(self._solid_mol)

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
        return get_frac(self.solid_mass).view(tuple_array)

    @property
    def solid_massnet(self):
        """solid phase net mass flow rate (kg/hr)."""
        return sum(self.solid_mass)

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
        return get_frac(self.solid_vol).view(tuple_array)

    @property
    def solid_volnet(self):
        """solid phase net volumetric flow rate (m3/hr)."""
        return self._phaseprop_molar_flownet('Vm', 's')*1000

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
        return get_frac(self._liquid_mol).view(tuple_array)

    @property
    def liquid_molnet(self):
        """liquid phase net molar flow rate (kg/hr)."""
        return sum(self.liquid_mol)

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
        return get_frac(self.liquid_mass).view(tuple_array)

    @property
    def liquid_massnet(self):
        """liquid phase net mass flow rate (kg/hr)."""
        return sum(self.liquid_mass)

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
        return get_frac(self.liquid_vol).view(tuple_array)

    @property
    def liquid_volnet(self):
        """liquid phase net volumetric flow rate (m3/hr)."""
        return self._phaseprop_molar_flownet('Vm', 'l')*1000

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
        return get_frac(self._LIQUID_mol).view(tuple_array)

    @property
    def LIQUID_molnet(self):
        """LIQUID phase net molar flow rate (kg/hr)."""
        return sum(self._LIQUID_mol)

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
        return get_frac(self.LIQUID_mass).view(tuple_array)

    @property
    def LIQUID_massnet(self):
        """LIQUID phase net mass flow rate (kg/hr)."""
        return sum(self.LIQUID_mass)

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
        return get_frac(self.LIQUID_vol).view(tuple_array)

    @property
    def LIQUID_volnet(self):
        """LIQUID phase net volumetric flow rate (m3/hr)."""
        return self._phaseprop_molar_flownet('Vm', 'L')*1000

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
        return get_frac(self._vapor_mol).view(tuple_array)

    @property
    def vapor_molnet(self):
        """vapor phase net molar flow rate (kg/hr)."""
        return sum(self._vapor_mol)

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
        return get_frac(self.vapor_mass).view(tuple_array)

    @property
    def vapor_massnet(self):
        """vapor phase net mass flow rate (kg/hr)."""
        return sum(self.vapor_mass)

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
        return get_frac(self.vapor_vol).view(tuple_array)

    @property
    def vapor_volnet(self):
        """vapor phase net volumetric flow rate (m3/hr)."""
        return self._phaseprop_molar_flownet('Vm', 'g')*1000

    ### Overall ###

    # Flow rates
    @property
    def mol(self):
        """Molar flow rates (kmol/hr)."""
        molflow = phases_molflow.get(self.phase)
        if molflow:
            return getattr(self, molflow)
        return sum(self._molflows).view(tuple_array)
        # return mixed_mol( self._flows )

    @mol.setter
    def mol(self, val):
        molflow = phases_molflow.get(self.phase)
        if molflow:
            return setattr(self, molflow, val)
        raise AttributeError('Cannot set flow rates if more than one phase is present')

    @property
    def mass(self):
        """Mass flow rates (kmol/hr)."""
        return (sum(self._molflows, 0) * self._MW).view(tuple_array)
    
    @property
    def vol(self):
        """Volumetric flow rates (kmol/hr)."""
        return sum(self._volflows).view(tuple_array)

    # Fractions
    @property
    def molfrac(self):
        """Molar fractions."""
        return get_frac(sum(self._molflows, 0)).view(tuple_array)

    @property
    def massfrac(self):
        """Mass fractions."""
        return get_frac(self.mass).view(tuple_array)

    @property
    def volfrac(self):
        """Volumetric fractions."""
        return get_frac(self.vol).view(tuple_array)

    # Net flow rates
    @property
    def molnet(self):
        """Mol net flow rate (kmol/hr)."""
        return np.array(self._molflows).sum()

    @property
    def massnet(self):
        """Mass net flow rate (kmol/hr)."""
        return self.solid_massnet + self.liquid_massnet + self.LIQUID_massnet + self.vapor_massnet

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
        mol = np.zeros(N*4)
        mus = np.zeros(N*4)
        Vms = np.zeros(N*4)
        start = 0
        end = N
        for phase, flow in phases_molflow.items():
            mol[start:end] = getattr(self, flow)
            mus[start:end] = self._phaseprop_list('mu', phase)
            Vms[start:end] = self._phaseprop_list('Vm', phase)
            start += N
            end += N
        Vm = self.Vm
        pos = np.where(mol != 0)
        mol = mol[pos]
        molfrac = mol/molnet
        return exp(sum(molfrac*ln(mus[pos]*Vms[pos])))/Vm 
    mu.__doc__ = Stream.mu.__doc__
    
    ### Material Properties ###

    # General for a given phase
    def _phaseprop_list(self, prop_ID, phase):
        """Return component property lists for given phase."""
        phase = phase.lower()
        out = np.zeros(self._Nspecies)
        mol = getattr(self, phases_molflow[phase])
        _species = self._species
        P = self.P
        T = self.T
        for i, sp_ID in self._index_ID.items():
            if mol[i] != 0:
                specie = getattr(_species, sp_ID)
                specie.P = P
                specie.T = T
                specie.phase = phase
                out[i] = getattr(specie, prop_ID)
        return out

    def _phaseprop_molar_flow(self, prop_ID, phase):
        """Return array of component properties * kmol/hr for given phase."""
        return self._phaseprop_list(prop_ID, phase) * getattr(self, phases_molflow[phase])

    def _phaseprop_molar_flownet(self, prop_ID, phase):
        """Return sum of component properties * kmol/hr for given phase."""
        return sum(self._phaseprop_molar_flow(prop_ID, phase))

    def _phaseprop_molar_mean(self, prop_ID, phase):
        """Return molar weighted average property for given phase."""
        mol = getattr(self, phases_molflow[phase])
        ppfn = self._phaseprop_list(prop_ID, phase) * mol
        return ppfn / sum(ppfn)

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

    def copy_like(self, stream):
        """Copy mol, T, P, and phase of stream to self."""
        if self._species is not stream._species:
            self._copy_species(stream)
        if isinstance(stream, MixedStream):
            self.solid_mol = copy.copy(stream.solid_mol)
            self.liquid_mol = copy.copy(stream.liquid_mol)
            self.LIQUID_mol = copy.copy(stream.LIQUID_mol)
            self.vapor_mol = copy.copy(stream.vapor_mol)
        elif isinstance(stream, Stream):
            self.empty()
            setattr(self, phases_molflow[stream.phase], stream.mol)
        else:
            raise TypeError('Must pass a valid Stream instance')
        self.P = stream.P
        self.T = stream.T

    def empty(self):
        """Set flow rates to zero."""
        self.solid_mol = 0
        self.liquid_mol = 0
        self.LIQUID_mol = 0
        self.vapor_mol = 0

    ### Equilibrium ###

    # Vapor-liquid equilibrium
    def VLE(self, specie_IDs=None, LNK=None, HNK=None, P=None,
            Qin=None, T=None, V=None, x=None, y=None):
        """Partition flow rates into vapor and liquid phases by equilibrium. Pressure defaults to current pressure.

        **Optional parameters**

            **P:** [float] Operating pressure (Pa)
                
            **specie_IDs:** [list] IDs of equilibrium species
                 
            **LNK:** list[str] Light non-keys
    
            **HNK:** list[str] Heavy non-keys

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
        sp = self._species
        _get_Psat = sp._get_Psat
        sp_index = self._ID_index

        # Set up indeces for both equilibrium and non-equilibrium species
        if specie_IDs == None:
            specie_IDs, index = self._equilibrium_species()
        else:
            index = [sp_index[specie] for specie in specie_IDs]
        if LNK:
            LNK_index = [sp_index[specie] for specie in LNK]
        else:
            LNK, LNK_index = self._light_species()
        if HNK:
            HNK_index = [sp_index[specie] for specie in HNK]
        else:
            HNK, HNK_index = self._heavy_species()

        # Get flow rates
        liquid_mol = self.liquid_mol
        vapor_mol = self.vapor_mol
        all_mol = self.liquid_mol + self.vapor_mol
        mol = all_mol[index]

        # Set light and heavy keys
        vapor_mol[LNK_index] = all_mol[LNK_index]
        vapor_mol[HNK_index] = 0
        liquid_mol[HNK_index] = all_mol[HNK_index]
        liquid_mol[LNK_index] = 0

        ### Single component equilibrium case ###
        
        Nspecies = len(specie_IDs)
        if Nspecies == 1:
            # Equilibrium based on one specie
            specie = getattr(sp, specie_IDs[0])
            
            if eq == 'TP':
                # Either liquid or gas
                Psat = specie.VaporPressure(T)
                if P < Psat:
                    liquid_mol[index] = 0
                    vapor_mol[index] = mol
                else:
                    liquid_mol[index] = mol
                    vapor_mol[index] = 0
                return
            
            elif eq == 'VP':
                # Set vapor fraction
                self.T = specie.Tsat(P)
                vapor_mol[index] = mol*V
                liquid_mol[index] = mol - vapor_mol[index]
                return
            
            elif eq == 'QinP': 
                # Set temperature in equilibrium
                self.T = specie.Tsat(P)
                
                # Check if super heated vapor
                vapor_mol[index] = mol
                liquid_mol[index] = 0
                H_dew = self.H
                if Hin > H_dew:
                    self.H = Hin
                    return
    
                # Check if subcooled liquid
                vapor_mol[index] = 0
                liquid_mol[index] = mol
                H_bubble = self.H
                if Hin < H_bubble:
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
        molnet = sum(mol)
        if molnet == 0:
            return  # No equilibrium necessary
        else:
            zf = mol/molnet
        
        min_ = np.zeros(Nspecies)
        act_coef = self.activity_coefficients
        
        ### Functions ###

        def v_given_TPmol(v_guess, T, P, mol):
            """Return equilibrium vapor flow rates given T, P and molar flow rates."""
            Psat = _get_Psat(specie_IDs, T)

            def v_error(v):
                """Error function for constant T and P, where v represents vapor flow rates."""
                l = mol - v
                x = l/sum(l)
                y = v/sum(v)
                err = x*Psat * act_coef(specie_IDs, x, T)/P - y
                return abs(err)
            
            return least_squares(v_error, v_guess,
                                 bounds=(min_, mol),
                                 ftol=0.001).x

        ### Get Equilibrium ###

        if eq == 'xP' or eq == 'yP':
            # Get temperature based on phase composition
            if eq == 'xP':
                self.T, y = self._bubble_point(specie_IDs, x, P)
            else:
                self.T, x = self._dew_point(specie_IDs, y, P)
                
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
        T_dew, x = self._dew_point(specie_IDs, zf, P)
        T_bubble, y = self._bubble_point(specie_IDs, zf, P)

        # Guestimate vapor composition for 'VP', 'TP', 'QinP' equilibrium
        if hasattr(self, '_y') and self._species_eq == specie_IDs:
            # Guess composition in the vapor is same as last time equlibrium was run
            y = self._y
        else:
            self._species_eq = specie_IDs
            # Guess composition in the vapor is a weighted average of boiling points
            Tb = sp._get_Tb(specie_IDs)
            try:
                T_min = min(Tb)
                T_max = max(Tb)
                split = (0.5 + (Tb - T_min)/(T_max - T_min))/1.5
                y = split * zf
                y = y/sum(y)
            except:
                y = zf

        if eq == 'VP':
            # Guess vapor flow rates
            pos = np.where(mol <= 0)
            mol[pos] = 10**-9
            v = y*V*molnet

            def V_error(T):
                nonlocal v
                if T >= T_dew:
                    return V - 1 - (T - T_dew)
                elif T <= T_bubble:
                    return V + (T - T_bubble)
                v = v_given_TPmol(v, T, P, mol)
                return (V - sum(v)/molnet)

            # Solve
            self.T = brentq(V_error, T_bubble, T_dew)
            v[pos] = 0
            vapor_mol[index] = v
            liquid_mol[index] = mol - vapor_mol[index]
            return
        
        elif eq == 'TP':
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
            if Hin >= H_dew:
                self.H = Hin
                return

            # Check if subcooled liquid
            vapor_mol[index] = 0
            liquid_mol[index] = mol
            self.T = T_bubble
            H_bubble = self.H
            if Hin <= H_bubble:
                self.H = Hin
                return

            # Guess T, overall vapor fraction, and vapor flow rates
            T_guess = T_bubble + (Hin - H_bubble)/(H_dew - H_bubble) * (T_dew - T_bubble)
            V_guess = (T_guess - T_bubble)/(T_dew - T_bubble)
            v = y * V_guess * mol

            def H_error_T(T):
                nonlocal v
                self.T = T
                vapor_mol[index] = v = v_given_TPmol(v, T, P, mol)
                liquid_mol[index] = mol - v
                return abs(Hin - (self.H))

            # Solve
            self.T = minimize_scalar(H_error_T, bounds=(T_bubble, T_dew),
                                     method='bounded').x
            self._y = v/sum(v)
            return 

    # LIQUID-liquid equilibrium
    def LLE(self, specie_IDs=None, split=None, lNK=(), LNK=(),
            solvents=(), solvent_split=(),
            P=None, T=None, Qin=0):
        """Partition flow rates into liquid and LIQUID phases by equilibrium.

        **Optional Parameters**

            **specie_IDs:** *tuple[str]* IDs of equilibrium species
            
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

        # Set up indeces for both equilibrium and non-equilibrium species
        sp_index = self._ID_index
        if specie_IDs is None:
            specie_IDs, _ = self._equilibrium_species()
        
        specie_IDs = list(specie_IDs)
        for specie in specie_IDs:
            if (specie in solvents) or (specie in lNK) or (specie in LNK):
                specie_IDs.remove(specie)
        specie_IDs = tuple(specie_IDs)
        
        index = [sp_index[specie] for specie in specie_IDs]
        lNK_index = [sp_index[specie] for specie in lNK]
        LNK_index = [sp_index[specie] for specie in LNK]
        solvent_index = [sp_index[specie] for specie in solvents]

        # Get min and max splits
        Nspecies = len(specie_IDs)
        min_ = np.zeros(Nspecies)

        # Set up activity coefficients
        activity_IDs = specie_IDs + solvents
        act_coef = self.activity_coefficients

        # Make sure splits are given as arrays
        if hasattr(self, '_lL_split') and self._species_eq == specie_IDs:
            split = self._lL_split  # cached split
        else:
            if split is None:
                species = self.species
                split_dipoles = [getattr(species, ID).dipole for ID in specie_IDs]
                solvent_dipoles = [getattr(species, ID).dipole for ID in solvents]
                dipoles = split_dipoles + solvent_dipoles
                split_dipoles = np.array(split_dipoles)
                try:
                    split = split_dipoles/max(dipoles)
                except TypeError as TE:
                    missing_dipoles = []
                    for i, is_missing in enumerate(split_dipoles==None):
                        if is_missing:
                            missing_dipoles.append(specie_IDs[i])
                    raise EquilibriumError(f'Cannot make estimate for split. Missing dipole values for species: {missing_dipoles}')
            self._species_eq = specie_IDs
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
            x1_gamma = act_coef(activity_IDs, x1_guess, T)
            x2_gamma = act_coef(activity_IDs, x2_guess, T)
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
        self._lL_split = split
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
        T_units = show_units.get('T')
        P_units = show_units.get('P')
        flow_units = show_units.get('flow')
        fraction = show_units.get('fraction')
        phases = self.phase
        index_ID = self._index_ID
        
        # First line
        unit_ID = self._source[0]
        if unit_ID is None:
            source = ''
        else:
            source = f'  from  {unit_ID}'
        unit_ID = self._sink[0]
        if unit_ID is None:
            sink = ''
        else:
            sink = f'  to  {unit_ID}'
        basic_info = f"{type(self).__name__}: {self.ID}{source}{sink}\n"

        # Default units
        show_format = self.show_format
        if not T_units:
            T_units = show_format.T
        if not P_units:
            P_units = show_format.P
        if not flow_units:
            flow_units = show_format.flow
        if fraction is None:
            fraction = show_format.fraction

        # Second line (thermo)
        T = Q_(self.T, units['T']).to(T_units).magnitude
        P = Q_(self.P, units['P']).to(P_units).magnitude
        basic_info += f" phase: '{self.phase}', T: {T:.5g} {T_units}, P: {P:.6g} {P_units}\n"
        
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
            return basic_info + ' flow:  0'
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
            nonzero, species = nonzero_species(index_ID, flow)
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
                beginning = (new_line_spaces[:-1] +
                             CS.dim(' species' + 
                             spaces + '  ' + flow_units) + 
                             '\n' + beginning)

            # Put it together
            phases_flowrates_info += beginning + flowrates + end + '\n\n'
            first_phase = False

        return basic_info + phases_flowrates_info[:-2]


Stream.VLE.__doc__ = MixedStream.VLE.__doc__
Stream.LLE.__doc__ = MixedStream.LLE.__doc__
MixedStream.mu.__doc__ = Stream.mu.__doc__


