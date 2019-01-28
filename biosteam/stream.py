#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 18 14:05:10 2018

@author: yoelr
"""
from biosteam import Q_, units_of_measure
import copy
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import newton, least_squares
from biosteam.utils import material_array, property_array, PropertyFactory, \
                           tuple_array, reorder, get_frac, Sink, Source
from biosteam.find import find
from biosteam.species import Species
from biosteam.exceptions import SolverError, DimensionError
from biosteam.UNIFAC import DORTMUND

ln = np.log
exp = np.exp
ndarray = np.ndarray
#plt.style.use('dark_background')
#CS = color_scheme

# %% TODOs

# TODO: add material property interphase when using Cape Open package


# %% Functions

def nonzero_species(index_ID, flow):
    index = []
    species = []
    for i, s in index_ID.items():
        if flow[i] != 0:
            index.append(i)
            species.append(s)
    return index, species


# %% Units of measure

mol_units = units_of_measure['mol']
mass_units = units_of_measure['mass']
vol_units = units_of_measure['vol']
T_units = units_of_measure['T']

mol_flow_dim = Q_(0, mol_units).dimensionality
mass_flow_dim = Q_(0, mass_units).dimensionality
vol_flow_dim = Q_(0, vol_units).dimensionality

class ShowFormat:
    __slots__ = ['_T', '_P', '_flow', '_fraction']
    def __init__(self):
        self._T = 'K'
        self._P = 'Pa'
        self._flow = 'kmol/hr'
        self._fraction = False
    
    @property
    def T(self):
        return self._T
    @T.setter
    def T(self, T):
        # Test to make sure it has temperature units
        Q_(0, T).ito('K')
        self._T = T
    
    @property
    def P(self):
        return self._P
    @P.setter
    def P(self, P):
        # Test to make sure it has pressure units
        Q_(0, P).ito('Pa')
        self._P = P
    
    @property
    def flow(self):
        return self._flow
    @flow.setter
    def flow(self, flow):
        # Test to make sure it has flow units
        flow_dim = Q_(0, flow).dimensionality
        invalid_dimension = True
        for dim in (mol_flow_dim, mass_flow_dim, vol_flow_dim):
            if flow_dim == dim:
                invalid_dimension = False
                break
        if invalid_dimension:
            raise DimensionError(f"Dimensions for flow units must be in molar, mass or volumetric flow rates, not '{flow_dim}'.")
        self._flow = flow
    
    @property
    def fraction(self):
        return self._fraction
    @fraction.setter
    def fraction(self, fraction):
        self._fraction = bool(fraction)
    
    def __repr__(self):
        dummy = '-in fractions-' if self.fraction else ''
        return f"<{type(self).__name__}: T='{self.T}', P='{self.P}', flow='{self.flow}'{dummy}>"

# %% Flow properties

@PropertyFactory
def MassFlow(self):
    """Mass flow (kg/hr)."""
    mol, MW = self.data
    if mol:
        return float(mol * MW)
    else:
        return 0.

@MassFlow.setter
def MassFlow(self, value):
    mol, MW = self.data
    if value:
        mol[0] = value/MW
    else:
        mol[0] = 0.

@PropertyFactory    
def VolumetricFlow(self):
    """Volumetric flow (m^3/hr)."""
    specie = self.name
    stream, mol = self.data
    if mol:
        specie.T = stream.T
        specie.P = stream.P
        specie.phase = stream.phase.lower()
        return float(getattr(specie, 'Vm') * mol * 1000)
    else:
        return 0.

@VolumetricFlow.setter
def VolumetricFlow(self, value):
    specie = self.name
    stream, mol = self.data
    if value:
        specie.T = stream.T
        specie.P = stream.P
        specie.phase = stream.phase.lower()
        mol[0] = value/(getattr(specie, 'Vm') * 1000)
    else:
        mol[0] = 0.


# %% Stream classes


class metaStream(type):
    """Metaclass for Stream."""
    __instances = []
    
    def __init__(cls, clsname, superclasses, new_definitions):
        # Add a Species object if not inherited
        if not hasattr(cls, '_species'):
            cls._species = Species()
            cls._show_format = ShowFormat()
        type(cls).__instances.append(cls)

    @property
    def species(cls):
        """[Species] Contains pure component thermodynamic properties for computing overall properties of Stream instances."""
        return cls._species

    @species.setter
    def species(cls, species):
        # Set Species object and working species
        if isinstance(species, Species):
            for cl in type(cls).__instances:
                # Set species
                cl._species = species
                
                # ID of all species
                cl._specie_IDs = IDs = species.ID
                
                # CAS of all species
                cl._CAS = CAS = species.CAS
                
                # Number of species
                cl._Nspecies = len(IDs)  
                
                # Indexes of the species (for flow rates)
                cl._index = index = np.array(range(cls._Nspecies))
                
                # Dictionary of ID to index
                cl._ID_index = dict(zip(IDs, index))
                
                # Dictionary of index to ID
                cl._index_ID = dict(zip(index, IDs))
                
                # Dictionary of CAS to index
                cl._CAS_index = dict(zip(CAS, index))
                
                # Dictionary of index to CAS 
                cl._index_CAS = dict(zip(index, CAS))
                
                # Set molecular weights for all species
                cl._MW = tuple_array(species.get_props(IDs, 'MW'))
        else:
            raise ValueError('Must pass a Species object')
            
    @property
    def show_format(cls):
        """[ShowFormat] Contains attributes with the units of measure for the show method."""
        return cls._show_format

    def __repr__(cls):
        if cls is Stream:
            return f'biosteam.{cls.__name__}'
        else:
            return f'Stream.{cls.__name__}'


class Stream(metaclass=metaStream):
    """Create a Stream object that defines a thermodynamic state along with material flow rates. Stream objects calculate thermodynamic and transport properties of the mixture and its components. The Species object define the chemical species avalilable for the Stream class. Ideal mixture is assumed for stream properties and excess energies [enthalpy and entropy] are ignored as a simplifying assumption for low pressure processes.

**Parameters**

    **ID:** [str] ID of the stream. If ID is '', a default ID will be chosen.

    **flow:** [tuple] All flow rates corresponding to `species`.

    **species:** tuple[str] or [Species] Species corresponding to `flow`. If empty, assume same species as class.

    **units:** [str] The units of the flow rates (only mass, molar, and volumetric flow rates are valid)

    **phase:** [str] 'l' for liquid, 'g' for gas

    **T:** [float] Temperature (K)

    **P:** [float] Pressure (Pa)

    ****flow_pairs:** Specie-flow pairs

**Class Properties**

    **species** = Species(): [Species] Contains pure component thermodynamic properties for computing overall properties of Stream instances

    **show_format**: [ShowFormat] Object with attributes that describe how the show method presents units

**Class Attributes**

    **_default_ID** = ['d', 1]: list[str, int] Default starting letter and current number for ID.

    **units** = {'C': 'kJ/K/hr', 'Cp': 'J/g/K', 'Cpm': 'J/mol/K', 'H': 'kJ/hr', 'S': 'kJ/hr', 'G': 'kJ/hr', 'U': 'kJ/hr', 'A': 'kJ/hr', 'Hf': 'kJ/hr', 'MW': 'g/mol', 'P': 'Pa', 'T': 'K', 'Vm': 'm^3/mol', 'alpha': 'm^2/s', 'k': 'W/m/K', 'mass': 'kg/hr', 'massfrac': 'kg/kg', 'massnet': 'kg/hr', 'mol': 'kmol/hr', 'molfrac': 'kmol/kmol', 'molnet': 'kmol/hr', 'mu': 'Pa*s', 'nu': 'm^2/s', 'rho': 'kg/m^3', 'rhom': 'mol/m^3', 'sigma': 'N/m', 'vol': 'm^3/hr', 'volfrac': 'm^3/m^3', 'volnet': 'm^3/hr'}
    : *[dict]* Units of measure for material properties. Items in this dictionary cannot be changed.

**Examples**

    Before making a stream, set the species using a Species object:

    .. code-block:: python

       >>> # Set Species object
       >>> Stream.species = Species('Ethanol', 'Water') 

    Stream objects may be created a variety of ways:

    .. code-block:: python

       >>> # Create a stream specifying specie and flow rate pairs:
       >>> s1 = Stream(ID='s1', Water=2)
       >>> s1.show()
       Stream: s1
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Water  2

       >>> # Create a stream assuming same specie order as given in specie_IDs:
       >>> s2 = Stream(ID='s2', flow=(1, 2))
       >>> s2.show()
       Stream: s2
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Ethanol  1
                        Water    2

       >>> # Create a stream passing flow rate units, phase, temperature and pressure:
       >>> s3 = Stream(ID='s3', flow=(0.278, 0.556), units='mol/s', phase='g', T=400, P=101325)
       >>> s3.show()
       Stream: s3
        phase: 'g', T: 400 K, P: 101325 Pa
        flow (kmol/hr): Ethanol  1.
                        Water    2.

       >>> # The working units do not change
       >>> s3.mol 
       material_array([1., 2.]) (kmol/hr)

    .. Warning:: Stream objects do not automatically calculate thermodynamic equilibrium. They simply assume the given phase, temperature and pressure are correct. To find equilibrium, use the VLE or LLE method.

    .. Note:: All following examples will use the Stream objects created here (s1, s2, and s3) starting from its initial values.

    Use the `show` method to print all specifications with desired units:

    .. code-block:: python

       >>> # Temperature in degree Celsius
       >>> s2.show(T='degC')
       Stream: s2 
        phase: 'l', T: 25. degC, P: 101325 Pa
        flow (kg/hr): Ethanol  46.1
                      Water    36
       
       >>> # Flow in kg/hr
       >>> s2.show(flow='kg/hr')
       Stream: s2 
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kg/hr): Ethanol  46.1
                      Water    36

       >>> # Flow in fractions
       >>> s2.show(fraction=True)
       Stream: s2
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow: Ethanol  0.333
              Water    0.667
                       3 kmol/hr

    Flow rates are stored internally as a material_array in the ‘mol’ attribute. A material_array object is a numpy ndarray which issues a RuntimeWarning when a negative or non-finite number is encountered:

    .. code-block:: python

       >>> # Set Water flow rate
       >>> s2.mol[1] = 18
       >>> s2.mol
       material_array([1, 18]) (kmol/hr)
       
       >>> # A negative flow issues a RuntimeWarning
       >>> s2.mol[1] = -1
       __main__:1: RuntimeWarning:
       Encountered negative or non-finite value in 'material_array' object.
       
       >>> # Setting the property sets its values in place
       >>> s2.mol = [1, 2]
       >>> s2.mol
       material_array([1, 2]) (kmol/hr)

    .. Note::

       material_array objects are numpy ndarrays which issue a RuntimeWarning when a negative or non-finite number is encountered.

    Mass and volumetric flow rates are also available as property_arrays of the molar flow rate. As such, they are always up to date with the molar flow rate and altering them also alters the molar flow rate:

    .. code-block:: python

       >>> # Altering mass or volumetric flows alters the molar flow
       >>> s2.vol
       property_array([0.059, 0.036]) (m^3/hr)
       >>> s2.vol[:] = [1, 1]
       >>> s2.mol
       material_array([17.06 , 55.343]) (kmol/hr)
       
       >>> # Values are always up to date with the molar flow
       >>> s2.mol[:] = [1, 2]
       >>> s2.vol
       property_array([0.059, 0.036]) (m^3/hr)

    .. Note::

       property_array objects are significantly slower than material_array objects. This is because flow rate data is internally stored as molar flow rates. Also, property_array objects are arrays of python objects, which add overhead over the C implementation of numpy. Whenever possible, use the material_array to manage flow rates. 

    Some thermodynamic/material properties, including enthalpy and heat capacity, are dependent only on temperature:

    .. code-block:: python
    
       >>> # Enthalpy (kJ/hr) and heat capacity (kJ/(kg-K))
       >>> s2.T = 298.15
       >>> s2.Cp
       3.200382054794244
       >>> s2.H 
       0.0
       
       >>> # Change of temperature
       >>> s2.T = 310
       >>> s2.H
       3140.9389548625936
       >>> s2.Cp
       3.2590872180092956
       
    .. Note:: Thermodynamic energies are relative to 25 degC and 1 atm.
    
    Some thermodynamic/material properties, including volume and density, are dependent on both temperature and pressure:

    .. code-block:: python
    
       >>> # Volumetric flow rate (m^3/hr) and density (kg/m^3)
       >>> s1.volnet
       0.036138079740245625
       >>> s1.rho
       997.0247522552814
      
       >>> # Change of pressure
       >>> s1.P *= 4
       >>> s1.volnet
       0.036136231155141196
       >>> s1.rho
       997.0757560552587
      
       >>> # Change of temperature
       >>> s1.T += 30
       >>> s1.volnet
       0.03663597700908213
       >>> s1.rho
       983.4747955832584

    A dictionary of available material properties is available in `Stream.units`. You may also find it useful to use the `help` method to search for a property:
        
    .. code-block:: python

       >>> Stream.help('conductivity')
       Thermal conductivity, k, is dependent on T and P with units, W/m/K, and type, float.

    """
    activity_coefficients = DORTMUND

    # [dict] Units of measure for material properties (class attribute). Items in this dictionary cannot be changed.
    units = units_of_measure

    # Utils information regarding properties
    _prop_molar_info = [
        # ID         # Description               # Dependency # Units      # Type
        ['T',        'temperature',              '',          'K',         'float'],
        ['H',        'enthalpy',                 'T',         'kJ/hr',     'float'],
        ['S',        'entropy',                  'TP',        'kJ/hr',     'float'],
        ['G',        'Gibbs free energy',        'TP',        'kJ/hr',     'float'],
        ['U',        'interal energy',           'TP',        'kJ/hr',     'float'],
        ['A',        'Helmholtz free energy',    'TP',        'kJ/hr',     'float'],
        ['Hf',       'enthalpy of formation',    '',          'kJ/hr',     'float'],
        ['P',        'pressure',                 '',          'Pa',        'float'],
        ['Cpm',      'molar heat capacity',      'T',         'J/mol/K',   'float'],
        ['Cp',       'specific heat capacity',   'T',         'J/kg/K',    'float'],
        ['Vm',       'molar volume',             'TP',        'm^3/mol',   'float'],
        ['rho',      'density',                  'TP',        'kg/m^3',    'float'],
        ['rhom',     'molar density',            'TP',        'mol/m^3',   'float'],
        ['nu',       'kinematic viscosity',      'TP',        'm^2/s',     'float'],
        ['mu',       'hydraulic viscosity',      'TP',        'Pa*s',      'float'],
        ['sigma',    'surface tension',          'T',         'N/m',       'float'],
        ['k',        'thermal conductivity',     'TP',        'W/m/K',     'float'],
        ['alpha',    'thermal diffusivity',      'TP',        'm^2/s',     'float'],
        ['Pr',       'Prantl number',            'TP',        "''",        'float'],
        ['mass',     'mass flow rates',          '',          'kg/hr',     'np.array'],
        ['mol',      'molar flow rates',         '',          'kmol/hr',   'np.array'],
        ['vol',      'volumetric flow rates',    'TP',        'm^3/hr',    'np.array'],
        ['massnet',  'net mass flow rate',       '',          'kg/hr',     'float'],
        ['molnet',   'net molar flow rate',      '',          'kmol/hr',   'float'],
        ['volnet',   'net volumetric flow rate', 'TP',        'm^3/hr',    'float'],
        ['massfrac', 'mass fractions',           '',          'kg/kg',     'np.array'],
        ['molfrac',  'molar fractions',          '',          'kmol/kmol', 'np.array'],
        ['volfrac',  'volumetric fractions',     'TP',        'm^3/m^3',   'np.array']]

    ### Class attributes for working species ###
    
    #: UNIFAC parameters cached
    # {specie_IDs: (chemgroups, rs, qs, group_counts)}
    _cached = {}
    
    #: Number of species
    _Nspecies = 0  
    
    #: Mol indexes
    _index = ()
    
    #: Dictionary to look up index by specie
    _ID_index = {}  
    
    #: Dictionary to look up specie by index
    _index_ID = {}  
    
    #: Array of molecular weights
    _MW = ()
    
    #: Species order for mass balance
    _specie_IDs = ()

    # [list] Default starting letter and current number for ID (class attribute)
    _default_ID = ['d', 1]

    #: Default ID
    _ID = ''

    def __init__(self, ID='', flow=(), species=(), units='kmol/hr',
                 phase='l', T=298.15, P=101325, **flow_pairs):
        # Get species ID and set species information
        if isinstance(species, Species):
            # Set species as provided
            self.species = species
            specie_IDs = ()
        elif hasattr(species, '__iter__'):
            specie_IDs = species
            # Copy class species data (in case it changes)
            self._species = species = self._species
            self._specie_IDs = self._specie_IDs
            self._CAS = self._CAS
            self._Nspecies = self._Nspecies
            self._index = self._index
            self._ID_index = self._ID_index
            self._index_ID = self._index_ID
            self._CAS_index = self._CAS_index 
            self._index_CAS = self._index_CAS
            self._MW = self._MW
        else:
            raise TypeError(f'species must be a Species object, not {type(species).__name__}')
        
        self.ID = ID
        
        #: Dew point cached
        # {specie_IDs: (pressure,
        #               temperature,
        #               vapor composition,
        #               liquid composition)}
        self._dew_cached = {}
        
        
        #: Bubble point cached
        # {specie_IDs: (pressure,
        #               temperature,
        #               vapor composition,
        #               liquid composition)}
        self._bubble_cached = {}
        
        #: tuple with source unit ID and outs position
        self._source = (None, None)
        
        #: tuple with sink unit ID and ins position
        self._sink = (None, None)
        
        #: [str] 'l'(iquid), or 'g'(as)
        self.phase = phase
        
        #: [float] Temperature (K)
        self.T = T
        
        #: [str] 'l' (liquid), or 'g' (gas)
        self.P = P  
        
        # Add specie flow pairs
        specie_IDs = tuple(specie_IDs) + tuple(flow_pairs.keys())
        flow = tuple(flow) + tuple(flow_pairs.values())
        
        # Match species and flow rates
        l_flow = len(flow)
        l_species = len(specie_IDs)
        if l_species == 0 and l_flow == 0:
            flow = np.zeros(self._Nspecies)
        elif l_species == l_flow:
            flow = reorder(flow, specie_IDs, self._specie_IDs)
        elif l_flow != self._Nspecies:
            raise ValueError('Length of flowrates must be equal to length of species')

        # Set molar flow rates
        flow = np.array(flow, dtype=np.float64)
        flow_wt_units = Q_(flow, units)
        dim = flow_wt_units.dimensionality
        if dim == mol_flow_dim:
            arr = flow_wt_units.to(self.units['mol']).magnitude
            self._mol = material_array(arr)
        elif dim == mass_flow_dim:
            arr = flow_wt_units.to(mass_units).magnitude / self._MW
            self._mol = material_array(arr)
        elif dim == vol_flow_dim:
            self._mol = material_array(flow)
            self.set_vol(flow_wt_units.to(vol_units).magnitude)
        else:
            raise DimensionError(f"Dimensions for flow units must be in molar, mass or volumetric flow rates, not '{dim}'.")

        mol = self._mol
        MW = self._MW
        
        # Mass flow rates
        mass = []
        for i, ID in enumerate(species.ID):
            m = MassFlow(ID, (mol[i:i+1], MW[i]))
            mass.append(m)
        self._mass = property_array(mass)
        
        # Volumetric flow rates    
        vol = []
        for i, specie in enumerate(species):
            v = VolumetricFlow(specie, (self, mol[i:i+1]))
            vol.append(v)
        self._vol = property_array(vol)

    # Forward pipping
    def __sub__(self, index):
        if isinstance(index, int):
            return Sink(self, index)
        elif isinstance(index, Stream):
            raise TypeError(f"Unsupported operand type(s) for -: '{type(self)}' and '{type(index)}'")
        return index.__rsub__(self)
        
    def __rsub__(self, index):
        if isinstance(index, int):
            return Source(self, index)
        elif isinstance(index, Stream):
            raise TypeError(f"Unsupported operand type(s) for -: '{type(self)}' and '{type(index)}'")
        return index.__sub__(self)

    # Backward pipping    
    __pow__ = __sub__
    __rpow__ = __rsub__
    
    @property
    def species(self):
        """[Species] Contains pure component thermodynamic properties for computing overall properties of Stream instances."""
        return self._species

    @species.setter
    def species(self, species):
        if self._ID != '':
            raise AttributeError("can't set attribute")
        # Set Species object and working species
        elif isinstance(species, Species):
            # Set species
            self._species = species
            
            # ID of all species
            self._specie_IDs = IDs = species.ID
            
            # CAS of all species
            self._CAS = CAS = species.CAS
            
            # Number of species
            self._Nspecies = len(IDs)  
            
            # Indexes of the species (for flow rates)
            self._index = index = np.array(range(self._Nspecies))
            
            # Dictionary of ID to index
            self._ID_index = dict(zip(IDs, index))
            
            # Dictionary of index to ID
            self._index_ID = dict(zip(index, IDs))
            
            # Dictionary of CAS to index
            self._CAS_index = dict(zip(CAS, index))
            
            # Dictionary of index to CAS 
            self._index_CAS = dict(zip(index, CAS))
            
            # Set molecular weights for all species
            self._MW = tuple_array(species.get_props(IDs, 'MW'))
        else:
            raise ValueError('Must pass an instance of Species')
            
    show_format = metaStream.show_format
    
    @property
    def ID(self):
        """Unique identification (str). If set as '', it will choose a default ID."""
        return self._ID

    @ID.setter
    def ID(self, ID):
        # Remove old reference to this object
        if self._ID != '' and ID != '':
            del find.stream[self._ID]

        # Get current default ID
        Stream = type(self)
        letter, number = Stream._default_ID

        # Make sure given ID is not a default ID
        if ID.startswith(letter):
            if ID[1:].isdigit():
                raise ValueError(f"IDs starting with '{letter}' and followed by only digits are reserved for defaults. To use the default ID, set ID as 'Default'.")

        # Select a default ID if requested
        if ID == '':
            Stream._default_ID[1] += 1
            ID = letter + str(number)

        # Add ID to find dictionary and set it
        find.stream[ID] = self
        self._ID = ID

    @property
    def MW(self):
        """Molecular weight of all species (array g/mol):     

        >>> s2.MW 
        [46.06844, 18.01528]
        """
        return self._MW

    ### Flow properties ###
    
    @classmethod
    def get_index(cls, IDs, CAS=False):
        """Get flow indeces of specified species.

        **Parameters**

             IDs: list[str] or [str] Specie names

        Example

        Get indeces by user defined ID:
        
        >>> s1.get_index('Water')
        1
        >>> s1.get_index(['Water', 'Ethanol'])
        [1, 0]

        Get indeces by CAS number:
        
        >>> s1,get_index('64-17-5'):
        0
        >>> s1,get_index(['64-17-5', '7732-18-5']):
        [0, 1]

        """
        if CAS:
            index_dict = cls._CAS_index
        else:
            index_dict = cls._ID_index
        
        if type(IDs) is str:
            indices = index_dict[IDs]
        else:
            indices = [index_dict[ID] for ID in IDs]
        
        return indices
    get_indx = get_index #: Alias of get_index.
    
    # Molar flow
    @property
    def mol(self):
        """Array of molar flow rates by specie (array kmol/hr):

        >>> s1.mol
        material_array([0, 2]) (kmol/hr)
        >>> s1.mol = [1, 2]
        >>> s1.mol
        material_array([1, 2]) (kmol/hr)
        """
        return self._mol

    @mol.setter
    def mol(self, val):
        try:
            self._mol[:] = val
        except ValueError:
            raise ValueError(f'Length of flow property ({len(val)} given) must be the same as the number of species ({self._Nspecies})')

    @property
    def molfrac(self):
        """Array of molar fractions (molar fractions).

        >>> s1.molfrac
        tuple_array([0., 1.])
        """
        return get_frac(self._mol).view(tuple_array)

    @property
    def molnet(self):
        """Net molar flow rate (kmol/hr)

        >>> s2.molnet
        3
        """
        return sum(self._mol)

    # Mass flow
    @property
    def mass(self):
        """Array of mass flow rates by specie (array kg/hr)

        >>> s1.mass
        property_array([ 0.   , 36.031])
        """
        return self._mass

    @mass.setter
    def mass(self, val):
        self._mass[:] = val

    @property
    def massfrac(self):
        """Array of mass fractions (mass fractions).

        >>> s1.massfrac
        tuple_array([0, 1])
        """
        return get_frac(self._mol * self._MW).view(tuple_array)

    @property
    def massnet(self):
        """Net mass flow rate (kg/hr)

        >>> s2.massnet
        82.099
        """
        return sum(self.mass)

    # Volumetric flow
    @property
    def vol(self):
        """Array of volumetric flow rates by specie (array m^3/hr). T and P dependent.

        >>> s2.vol
        property_array([0.059, 0.036])
        
        """
        return self._vol

    @vol.setter
    def vol(self, val):
        self._vol[:] = val

    @property
    def volfrac(self):
        """Array of volumetric fractions (volumetric fractions). T and P dependent.

        >>> s2.volfrac
        tuple_array([0.619, 0.381])
        """

    @property
    def volnet(self):
        """Net volumetric flow rate (m^3/hr). T and P dependent.

        >>> s2.volnet # a liquid stream
        0.09475552896916632
        >>> s3.volnet # a gas stream
        97.32699378482434
        """
        return self._prop_molar_flownet('Vm')*1000

    # Set non-molar flows
    def set_mass(self, mass_flow, index=None):
        """Set flow rates by mass in kg/hr.

        **Parameters**

             mass_flow: [array_like] Mass flow rates.

             index: array_like[int] Indeces of species corresponding to mass flow rates.

        .. code-block:: python

           >>> s1.set_mass(20 , index = 1)
           >>> s1.mass
           tuple_array([0., 20])
           >>> # No need to specify indeces when all flow rates are specified
           >>> s1.set_mass([0, 36.031])
           >>> s1.mass
           tuple_array([0, 36.031])

        """
        if index is None:
            self.mol = mass_flow / self._MW
        else:
            self.mol[index] = mass_flow / self._MW[index]
    
    def set_vol(self, vol_flow, index=None):
        """Set flow rates by volume in m^3/hr.

        **Parameters**

             vol_flow: [array_like] Volumetric flow rates.

             index: array_like[int] Indeces of species corresponding to volumetric flow rates.

        .. code-block:: python

           >>> s1.set_vol(20 , index = 1)
           >>> s1.vol
           tuple_array([0., 20])
           >>> # No need to specify indeces when all flow rates are specified
           >>> s1.set_vol([0, 1])
           >>> s1.vol
           tuple_array([0, 1])

        """
        vol_flow = np.array(vol_flow)
        species = self.species
        if index is None:
            Vms = species.get_props(self._specie_IDs, 'Vm', self.T, self.P, self.phase)
            self.mol = vol_flow / Vms / 1000
        else:
            index_sp = self._index_ID
            if hasattr(index, '__len__'):
                specie_IDs = (index_sp[i] for i in index)
            else:
                specie_IDs = index_sp[index]
            Vms = species.get_props(specie_IDs, 'Vm', self.T,
                                    self.P, self.phase)
            self.mol[index] = vol_flow / Vms / 1000

    ### Energy flows ###

    # Enthalpy
    @property
    def H(self):
        """Enthalpy flow rate (kJ/hr) without formation energies. Only T dependent.

        >>> s1.H # The stream is at the reference state
        0.0
        >>> s1.H = 1000
        >>> s1.T
        304.80600901692236
        >>> s1.H
        1002.1604769775776

        .. Note:: The solver reaches a level of tolerance and the new temperature is an approximation.
        """
        return self._prop_molar_flownet('H')

    @H.setter
    def H(self, H):
        if all(self.mol == 0):
            if H != 0:
                raise ValueError(f"Cannot set nonzero enthalpy "
                                +f"to empty stream '{self.ID}'.")
            else:
                return

        # First approximation
        T_old = self.T
        C = self.C
        self.T += (H - self.H)/C

        # Solve enthalpy by iteration
        it = 1
        while abs(self.T - T_old) > 0.01:
            T_old = self.T
            self.T += (H - self.H)/C
            it += 1
            if it > 40:
                raise SolverError(f"Could not solve temperature "
                                 +f"for given enthalpy in stream '{self.ID}'")

    @property
    def Hf(self):
        """Heat of formation flow rate (kJ/hr).

        >>> s1.Hf
        -483640.0
        """
        return self._prop_molar_flownet('Hfm')

    # Entropy
    @property
    def S(self):
        """Entropy flow rate (kJ/hr) without formation energies. T and P dependent.

        >>> s1.S # The stream is at the reference state
        0.0
        >>> s1.T = 320
        >>> s1.S
        10.642980522582695
        """
        return self._prop_molar_flownet('S')

    # Gibbs free energy
    @property
    def G(self):
        """Gibbs free energy flow rate (kJ/hr) without formation energies. T and P dependent.

        >>> s1.G # The stream is at the reference state
        0.0
        >>> s1.T = 320
        >>> s1.G
        -117.6450204179555
        """
        return self.H - self.S*self.T

    # Internal energy
    @property
    def U(self):
        """Internal energy flow rate (kJ/hr) without formation energies. T and P dependent.

        >>> s1.U # The stream is at the reference state
        0.0
        >>> s1.T = 320
        >>> s1.U
        -409.16014979602505
        """
        return self.H - self.P*self.volnet

    # Helmholtz
    @property
    def A(self):
        """Helmholtz energy flow rate (kJ/hr) without formation energies. T and P dependent.

        >>> s1.A # The stream is at the reference state
        0.0
        >>> s1.T = 320
        >>> s1.A
        -3814.9139170224876
        """
        return self.U - self.T*self.S

    # Capacity flow rate
    @property
    def C(self):
        """Heat capacity flow rate (kJ/K/hr). T dependent.

        >>> s2.C
        262.74816631655267
        """
        return self._prop_molar_flownet('Cpm')

    # Material properties
    @property
    def Cp(self):
        """Specific heat capacity (J/g/K). T dependent.

        >>> s1.Cp
        4.180597021827335
        """
        return self.Cpm*self.molnet/self.massnet

    @property
    def Cpm(self):
        """Molar heat capacity (J/mol/K). T dependent.

        >>> s1.Cpm # (J/mol/K)
        75.31462591538556
        """
        return self.C/self.molnet

    @property
    def Vm(self):
        """Molar volume (m^3/mol). T and P dependent.

        >>> s1.Vm
        1.8069039870122814e-05
        """
        return self._prop_molar_mean('Vm')

    @property
    def rho(self):
        """Density (kg/m^3). T and P dependent.

        >>> s1.rho
        997.0247522552814
        """
        return self.massnet/self.volnet

    @property
    def rhom(self):
        """Molar density (mol/m^3). T and P dependent.

        >>> s1.rhom
        55343.283715561534
        """
        return self.molnet/self._prop_molar_flownet('Vm')

    @property
    def nu(self):
        """Kinematic viscosity (m^2/s). T and P dependent.

        >>> s1.nu
        9.154438858391748e-07
        """
        return self.mu/self.rho

    @property
    def mu(self):
        """Hydraulic viscosity (Pa*s). T and P dependent.

        >>> s1.mu
        0.0009127202134824155
        """
        # Katti, P.K.; Chaudhri, M.M. (1964). "Viscosities of Binary Mixtures of Benzyl Acetate with Dioxane, Aniline, and m-Cresol". Journal of Chemical and Engineering Data. 9 (1964): 442–443.
        molfrac = self.molfrac
        mus = np.array(self._prop_list('mu'))
        Vms = np.array(self._prop_list('Vm'))
        Vm = self.Vm
        pos = np.where(molfrac != 0)
        return exp(sum(molfrac[pos]*ln(mus[pos]*Vms[pos])))/Vm 

    @property
    def k(self):
        """Thermal conductivity (W/m/k). T and P dependent.

        >>> s1.k
        0.5942044328004411
        """
        return self._prop_molar_mean('k')

    @property
    def alpha(self):
        """Thermal diffusivity (m^2/s). T and P dependent.

        >>> s1.alpha
        1.4255801521655763e-07
        """
        return self._prop_molar_mean('alpha')

    @property
    def sigma(self):
        """Surface tension (N/m). T dependent.

        >>> s1.sigma
        0.07205503890847455
        """
        return self._prop_molar_mean('sigma')

    @property
    def Pr(self):
        """Prandtl number (non-dimensional). T and P dependent.

        >>> s1.Pr
        6.421553249380879
        """
        return self._prop_molar_mean('Pr')

    # Other properties
    @property
    def source(self):
        """tuple with the ID of the unit it exits and the position"""
        return self._source

    @property
    def sink(self):
        """tuple with the ID of the unit it enters and the position"""
        return self._sink

    @property
    def nonzero_species(self):
        """Return flow indeces and species that have a non-zero flow rate.

        **Return**

             **index:** list[float] flow indexes that are not zero

             **species:** list[str] IDs for species corresponding to index

        >>> s1.nonzero_species
        [1], ['Water']
        """
        return nonzero_species(self._index_ID, self.mol)
    
    def prop_quantity(self, prop_ID):
        """Return a property as a Quantity object as described in the `pint package <https://pint.readthedocs.io/en/latest/>`__ 

        **Parameters**

             **prop_ID:** [str] name of the property (e.g. 'mol', 'H', 'k' ...)

        Example:

        >>> s1.prop_quantity('Cp')
        <Quantity(4.180597021827335, 'joule / gram / kelvin')>

        """
        return Q_(getattr(self, prop_ID), self.units[prop_ID])

    # Sensitivity analysis
    def _prop_table(self, prop_ID, var_ID, var_vals, index=None):
        """Return an array of property values at 'var_vals'.

        **Parameters**

             **prop_ID:** [str] Name of the dependent property (e.g. 'rho', 'Cp', 'H', 'volnet' ...)

             **var_ID:** [str] Name of the independent variable (e.g 'T', 'P', 'molnet' ...)

             **var_vals:** [numpy array] Independent variable values at which to compute the property

             **index:** [array_like or None] Optional indices if the variable is an array and the given values are a subset of the array.

        Example

        >>> s1._prop_table( 'H', 'T', [298.15, 299.15, 300.15, 301.15, 302.15])
        array([  0.   ,  65.041, 130.231, 195.569, 261.055])

        >>> s1._prop_table( 'rho', 'mol' , [0, 1, 2, 3, 4], index = 0) # vary the amount of Ethanol
        array([997.025, 866.43 , 835.659, 821.904, 814.109])
        """
        if index is not None:
            def get_var():
                return getattr(self, var_ID)[index]

            def set_var(xf):
                getattr(self, var_ID)[index] = xf
        else:
            def get_var():
                return getattr(self, var_ID)

            def set_var(xf):
                setattr(self, var_ID, xf)
        x0 = get_var()
        l = len(var_vals)
        y = np.zeros(l)
        for i in range(l):
            val = var_vals[i]
            set_var(val)
            y[i] = getattr(self, prop_ID)
        set_var(x0)
        return y

    def prop_plot(self, prop_ID, var_ID, var_vals, index=None):
        """Make a plot of a property vs a variable.

        **Parameters**

             **prop_ID:** [str] Name of the property (e.g. 'rho', 'Cp', 'H', 'volnet' ...)

             **var_ID:** [str] Name of the variable (e.g 'T', 'P', 'molnet' ...)

             **var_vals:** [numpy array] Independent variable values at which to compute the property

             **index:** [array_like or None] Optional indices if the variable is an array

        Example

        >>> s1.prop_plot( 'H', 'T', [298.15, 299.15, 300.15, 301.15, 302.15])

        .. image:: prop_plot.png
           :scale: 70 %
           :align: center
        """
        y = self._prop_table(prop_ID, var_ID, var_vals, index)
        plt.plot(var_vals, y)
        units = self.units
        if index is not None:
            var_ID = self._index_ID[index] + ' ' + var_ID
        plt.xlabel(var_ID + ' (' + units[var_ID] + ')', fontsize=16)
        plt.ylabel(prop_ID + ' (' + units[prop_ID] + ')', fontsize=16)

    # Derivative of a property and sensivity analysis
    def prop_derivative(self, prop_ID, var_ID, index=None):
        """Return the derivative of property 'prop_ID' and variable 'var_ID'. If the property given by var_ID is an array, use index to specify which element.

         **Parameters**

             **prop_ID:** [str] Name of the property (e.g. 'rho', 'Cp', 'H', 'volnet' ...)

             **var_ID:** [str] Name of the variable (e.g 'T', 'P', 'molnet' ...)

             **index:** [array_like or None] Optional indices if the variable is an array

        Example

        >>> s1.prop_derivative('rho', 'T')
        -0.4010588564718449

        """
        if index is not None:
            def get_var():
                return getattr(self, var_ID)[index]

            def set_var(xf):
                getattr(self, var_ID)[index] = xf
        else:
            def get_var():
                return getattr(self, var_ID)

            def set_var(xf):
                setattr(self, var_ID, xf)
        x0 = get_var()
        y0 = getattr(self, prop_ID)
        xf = x0 + 10**-6
        set_var(xf)
        yf = getattr(self, prop_ID)
        set_var(x0)
        return (yf-y0)/(xf-x0)

    def _prop_derivative_table(self, prop_ID, var_ID, var_vals, index=None):
        """Return an array of property derivatives at 'var_vals' 

        **Parameters**

             **prop_ID:** [str] Name of the dependent property (e.g. 'rho', 'Cp', 'H', 'volnet' ...)

             **var_ID:** [str] Name of the independent variable (e.g 'T', 'P', 'molnet' ...)

             **var_vals:** [numpy array] Independent variable values at which to compute the property

             **index:** [array_like or None] Optional indices if the variable is an array

        Example

        >>> s1._prop_derivative_table('rho', 'T', [298.15, 300.15, 302.15, 304.15, 306.15, 308.15])
        array([-0.401, -0.408, -0.415, -0.422, -0.429, -0.436])

        """
        if index is not None:
            def get_var():
                return getattr(self, var_ID)[index]

            def set_var(xf):
                getattr(self, var_ID)[index] = xf
        else:
            def get_var():
                return getattr(self, var_ID)

            def set_var(xf):
                setattr(self, var_ID, xf)

        def Dvar(x0):
            set_var(x0)
            y0 = getattr(self, prop_ID)
            xf = x0 + 10**-6
            set_var(xf)
            yf = getattr(self, prop_ID)
            return (yf-y0)/(xf-x0)
        
        x0 = get_var()
        length = len(var_vals)
        y = np.zeros(length)
        for i in range(length):
            val = var_vals[i]
            y[i] = Dvar(val)
        set_var(x0)
        return y

    def prop_derivative_plot(self, prop_ID, var_ID, var_vals, index=None):
        """Return a property derivative plot.

        **Parameters**

             **prop_ID:** [str] Name of the dependent property (e.g. 'rho', 'Cp', 'H', 'volnet')

             **var_ID:** [str] Name of the independent variable (e.g 'T', 'P', 'molnet')

             **var_vals:** [numpy array] Independent variable values at which to compute the property

             **index:** [array_like or None] Optional indices if the variable is an array

        Example

        >>> s1.prop_derivative_plot('rho', 'T', [298.15, 300.15, 302.15, 304.15, 306.15, 308.15])

        .. image:: prop_derivative_plot.png
           :scale: 70 %
           :align: center
        """
        y = self._prop_derivative_table(prop_ID, var_ID, var_vals, index)
        plt.plot(var_vals, y)
        units = self.units
        if index is not None:
            var_ID = self._index_ID[index] + ' ' + var_ID
        plt.xlabel(var_ID + ' (' + units[var_ID] + ')', fontsize=16)
        plt.ylabel(f'd{prop_ID}/d{var_ID}', fontsize=16)

    # General methods for getting ideal mixture properties
    def _prop_list(self, prop_ID):
        """Return component property list."""
        out = np.zeros(self._Nspecies)
        _species = self._species
        _specie_IDs = self._specie_IDs
        mol = self.mol
        P = self.P
        T = self.T
        phase = self.phase.lower()
        for i in self._index:
            if mol[i] != 0:
                specie = getattr(_species, _specie_IDs[i])
                specie.P = P
                specie.T = T
                specie.phase = phase
                out[i] = getattr(specie, prop_ID)
        return out

    def _prop_molar_flow(self, prop_ID):
        """Return array of component properties * kmol/hr."""
        return self._prop_list(prop_ID) * self.mol

    def _prop_molar_flownet(self, prop_ID):
        """Return sum of component properties * kmol/hr."""
        return sum(self._prop_molar_flow(prop_ID))

    def _prop_molar_mean(self, prop_ID):
        """Return molar weighted average property."""
        return self._prop_molar_flownet(prop_ID) / self.molnet

    # Specifications
    @staticmethod
    def like(stream, ID=''):
        """Create either a Stream or MixedStream object just like the given stream, depending on whether multiple phases are present."""
        s = stream
        phase = s.phase
        if len(phase) == 1:
            out = Stream(ID, flow=s.mol, phase=phase, T=s.T, P=s.P)
        else:
            MS = mixed_stream.MixedStream
            out = MS(ID, solid_flow=s.solid_mol, liquid_flow=s.liquid_mol,
                     LIQUID_flow=s.LIQUID_mol, vapor_flow=s.vapor_mol,
                     T=s.T, P=s.P)
        return out
    
    def copy_like(self, stream):
        """Copy mol, T, P, and phase of stream to self.

        >>> s1.copy_like(s3)
        >>> s1.show()
        Stream: s1
         phase: 'g', T: 400.00 K, P: 101325 Pa
         mol (kmol/hr): [Ethanol, 1]
                        [Water,   2]
        """
        phase = stream.phase
        if len(phase) == 1:    
            self.mol = copy.copy(stream.mol)
            self.P = stream.P
            self.phase = phase
            self.T = stream.T
        else:
            self.enable_phases()
            self.copy_like(stream)

    def empty(self):
        """Set net flow rate to zero

        >>> s1.empty()
        >>> s1.mol
        mol([0, 0]) kmol/hr
        """
        self.mol = 0
        
    def _equilibrium_species(self):
        """Return IDs and indexes of species in equilibrium."""
        specie_IDs, index = self.species._equilibrium_species
        mol = self.mol[index]
        index = mol!=0
        specie_IDs = tuple(sp for sp, i in zip(specie_IDs, index) if i)
        index = self.get_index(specie_IDs)
        return specie_IDs, index

    def _heavy_species(self):
        """Return IDs and indexes of heavy species not in equilibrium."""
        specie_IDs, index = self.species._heavy_species
        mol = self.mol[index]
        index = mol!=0
        specie_IDs = tuple(sp for sp, i in zip(specie_IDs, index) if i)
        index = self.get_index(specie_IDs)
        return specie_IDs, index

    def _light_species(self):
        """Return IDs and indexes of light species not in equilibrium."""
        specie_IDs, index = self.species._light_species
        mol = self.mol[index]
        index = mol!=0
        specie_IDs = tuple(sp for sp, i in zip(specie_IDs, index) if i)
        index = self.get_index(specie_IDs)
        return specie_IDs, index

    # Equilibrium
    def _default_eq_args(self, specie_IDs, z, P):
        """Return default arguments for equilibrium calculations."""
        # Make sure species and compositions are ok
        bol1 = z is None; bol2 = specie_IDs is None
        if (bol1 and not bol2):
            index = self.get_index(specie_IDs)
            mol = self.mol[index]
            z = mol/sum(mol)
        elif (not bol1 and bol2):
            raise ValueError('Composition was given with out specie IDs.')
        
        # Choose species that have UNIFAC groups
        if specie_IDs is None:
            # Get the molar fraction for those species that are not zero
            specie_IDs, index = self._equilibrium_species()
            mol = self.mol[index]
            z = mol/sum(mol)
        else:
            z = np.asarray(z)
            index = z!=0
            specie_IDs = tuple(sp for sp, i in zip(specie_IDs, index) if i)
            z = z[index]
        
        # Default to self pressure
        if P is None:
            P = self.P
        return specie_IDs, z, P
    
    def bubble_point(self, specie_IDs=None, x=None, P=None):
        """Bubble point given x and P. If no parameters are given, assume same composition as self.

        **Parameters**

            **specie_IDs:** iterable[str] Species corresponding to x.

            **x:** [array_like] The composition of the liquid phase.

            **P:** [float] The pressure in Pascals.

        **Returns**

            **T_bubble:** [float] The bubble point temperature

            **y:** [numpy ndarray] The composition of the vapor phase considering only the equilibrium species.
             
            **eq_species:** tuple[str] Species in equilibrium.

        >>> s1.bubble_point(specie_IDs = ('Ethanol', 'Water'),
        ...                 x = (0.6, 0.4), P = 101325)
        (352.2820850833474, array([0.703, 0.297]), ('Ethanol', 'Water'))

        """
        # If given just one specie, return Tsat
        if specie_IDs is not None:
            if isinstance(specie_IDs, str):
                return getattr(self.species, specie_IDs).Tsat(self.P), 1, specie_IDs
            if len(specie_IDs) == 1:
                return getattr(self.species, specie_IDs[0]).Tsat(self.P), np.array((1,)), specie_IDs
        
        # If just one specie in equilibrium, return Tsat
        specie_IDs, x, P = self._default_eq_args(specie_IDs, x, P)
        if len(specie_IDs) == 1:
            return getattr(self.species, specie_IDs[0]).Tsat(self.P), 1, specie_IDs
        
        # Solve and return bubble point
        x, P = self._bubble_point(specie_IDs, x, P)
        return x, P, specie_IDs
    
    def _bubble_point(self, specie_IDs, x, P):
        """Same as bubble_point but does not return equilibrium species."""
        # Setup functions to calculate vapor pressure and activity coefficients
        sp = self._species
        _get_Psat = sp._get_Psat
        act_coef = self.activity_coefficients
        x = np.array(x)
        
        # Retrive cached info
        # Even if different composition, use previous bubble point as guess
        cached = self._bubble_cached.get(specie_IDs) # c means cached
        if cached:
            cP, T_bubble, y_bubble, cx = cached
            if sum(abs(x - cx)) < 1e-6 and abs(cP - P) < 1:
                # Return cached data
                return T_bubble, y_bubble 
        else:
            # Initial temperature guess
            T_bubble = sum(np.array(x) * sp._get_Tb(specie_IDs))

        # Bubble point given x and P
        gamma = None
        Psat = None
        def bubble_point_error(T):
            nonlocal gamma, Psat
            gamma = act_coef(specie_IDs, x, T)
            Psat = _get_Psat(specie_IDs, T)
            return 1 - sum(x * Psat * gamma)/P

        # Solve and return
        T_bubble = newton(bubble_point_error, T_bubble, tol= 1e-7)
        x /= sum(x)
        y_bubble = x * Psat * gamma / P
        y_bubble /= sum(y_bubble)
        self._bubble_cached[specie_IDs] = (P, T_bubble, y_bubble, x)
        return T_bubble, y_bubble

    def dew_point(self, specie_IDs=None, y=None, P=None):
        """Dew point given y and P.

        **Parameters**

            **specie_IDs:** iterable[str] Species corresponding to y.

            **y:** [array_like] The composition of the vapor phase.

            **P:** [float] The pressure in Pascals.

        **Returns**

            **T_dew:** [float] The dew point temperature

            **x:** [numpy array] The composition of the liquid phase considering only the equilibrium species.
             
            **eq_species:** tuple[str] Species in equilibrium.


        >>> s1.dew_point(specie_IDs = ('Ethanol', 'Water'),
        ...              y = (0.5, 0.5, P = 101325)
        (353.7420742149825, array([0.556, 0.444]), ('Ethanol', 'Water'))

        """
        # If given just one specie, get Tsat
        if specie_IDs is not None:
            if isinstance(specie_IDs, str):
                return getattr(self.species, specie_IDs).Tsat(self.P), 1, specie_IDs
            if len(specie_IDs) == 1:
                return getattr(self.species, specie_IDs[0]).Tsat(self.P), np.array((1,)), specie_IDs
        
        # If just one specie in equilibrium, return Tsat
        specie_IDs, y, P = self._default_eq_args(specie_IDs, y, P)
        if len(specie_IDs) == 1:
            return getattr(self.species, specie_IDs[0]).Tsat(self.P), 1, specie_IDs
        
        # Solve and return dew point
        y, P = self._dew_point(specie_IDs, y, P)
        return y, P, specie_IDs
    
    def _dew_point(self, specie_IDs, y, P):
        """Same as bubble_point but does not return equilibrium species."""
        # Setup functions to calculate vapor pressure and activity coefficients
        sp = self._species
        _get_Psat = sp._get_Psat
        act_coef = self.activity_coefficients
        y = np.array(y)
        
        # Retrive cached info
        # Even if different composition, use previous bubble point as guess
        cached = self._dew_cached.get(specie_IDs)
        if cached:
            cP, T_dew, cy, x_dew = cached # c means cached
            if sum(abs(y - cy)) < 1e-6 and abs(cP - P) < 1:
                # Return cached data
                return T_dew, x_dew 
        else:
            # Initial temperature guess
            T_dew = sum(y * np.asarray(sp._get_Tb(specie_IDs)) )
            x_dew = y

        # Error function with constant y and P
        def dew_point_error(x):
            nonlocal T_dew, x_dew
            gamma = None
            Psat = None

            # Error for dew point temperature, constant x and P
            def dew_point_error_T(T):
                nonlocal gamma, Psat
                Psat = np.asarray(_get_Psat(specie_IDs, T))
                gamma = act_coef(specie_IDs, x, T)
                return 1 - sum(y/(Psat*gamma)) * P

            x_dew = x/sum(x)
            T_dew = newton(dew_point_error_T, T_dew, tol= 1e-7)
            return abs(x_dew - y*P/(Psat*gamma))

        # Get min and max splits
        Nspecies = len(specie_IDs)
        min_ = np.zeros(Nspecies)
        max_ = np.ones(Nspecies)

        # Find x by least squares
        try:
            x_dew = least_squares(dew_point_error, x_dew, bounds=(min_, max_), xtol= 1e-6).x
            y /= sum(y)
        except:
            Tmax = [getattr(sp, s).VaporPressure.Tmax for s in specie_IDs]
            T_dew = min(Tmax) - 1
            x_dew = None
        
        self._dew_cached[specie_IDs] = (P, T_dew, y, x_dew)
        return T_dew, x_dew

    # Dimensionless number methods
    def Re(self, L, A=None):
        """Return Reynolds number.

        **Parameters**

             **L:** [float] Characteristic lenght (m)

             **A:** [float] Cross-sectional area (m^2). If A is None, assume flow through a cylindrical pipe.

        Example

        >>> s1.Re(0.2, 0.031416)
        373.43769007101764
        """
        volnet = Q_(self.volnet, self.units['vol']).to('m^3/s').magnitude
        if A is None:  # Assume pipe
            return 4*volnet/(self.nu*np.pi*L)
        return volnet*L/(self.nu*A)

    # Class methods
    @classmethod
    def _helplist(cls, prop):
        """Return information related to a property.

        **Parameters**

             **prop:** [str] name or description of a property

        Return

             **out:** list[str] where elements are:
                  * **ID:** [str] name of the property

                  * **Description:** [str] description of the property

                  * **Dependency:** [str] TP dependency

                  * **Units:** [str] units of measure

                  * **Datatype:** [str] type returned by property

        Example

        >>> Stream.help('k')
        ['k', 'thermal conductivity', 'TP', 'W/m/K', 'float']

        >>> Stream.help('viscosity')
        [['nu', 'kinematic viscosity', 'TP', 'm^2/s', 'float'],
         ['mu', 'hydraulic viscosity', 'TP', 'Pa*s', 'float']]

        >>> Stream.help('kmol/hr')
        [['mol', 'molar flow rates', '', 'kmol/hr', 'np.array'],
         ['molnet', 'net molar flow rate', '', 'kmol/hr', 'float']]
        """
        out = []
        # Compare by both ID and description
        for index in (0, 1):
            # Search for exact matches
            for l in cls._prop_molar_info:
                if prop == l[index]:
                    return l
            # If no matches found, search harder
            if len(out) == 0:
                for l in cls._prop_molar_info:
                    if prop.lower() in l[index].lower():
                        out.append(l)
        return out

    @classmethod
    def help(cls, prop):
        """Print information related to a property.

        **Parameters**

             **prop:** [str] name or description of a property     

        Example

        >>> Stream.help('rho')
        Density, rho, is dependent on T and P with units, kg/m^3, and type, float.
        >>> Stream.help('density')
        Density, rho, is dependent on T and P with units, kg/m^3, and type, float.
        """
        def make_nice(helplist):
            """Make the helplist into a nice string"""

            # No helplist
            if helplist == []:
                pass

            # Only one helplist, print a nice string
            elif type(helplist[0]) is str:
                propID, description, dependency, units, datatype = helplist
                if dependency == 'TP':
                    dependency = 'is dependent on T and P '
                elif dependency == 'T':
                    dependency = 'is dependent on T '
                elif dependency == 'P':
                    dependency = 'is dependent on P '
                print(f"{description.capitalize()}, {propID}, {dependency}with units, {units}, and type, {datatype}.")

            # Many helplists, print all the nice strings
            elif type(helplist[0]) is list:
                for i in helplist:
                    make_nice(i)
                    print()
            else:
                raise Exception('Unknown error with ' + str(helplist))

        make_nice(cls._helplist(prop))

    @staticmethod
    def eqT(streams, T_guess=None, Q_in=0, approximate=True):
        """Bring all streams to temperature equilibrium.

        **Parameters**

            **streams:** [Stream] All streams
        
            **T_guess:** [float] Equilibrium temperature guess (K)

            **Q_in:** [float] Heat addition to streams (kJ/hr)     
            
            **approximate:** [bool] If True, approximate energy balance with constant heat capacity

        **Return**

             **T:** [float] New temperature (K)

        >>> Stream.eqT((s1, s2, s3))
        >>> (s1.T, s2.T, s3.T)
        (346.06735249999997,  346.06735249999997,  346.06735249999997)
        """
        # Set up problem
        C = sum(s.C for s in streams)
        H_ = sum(s.H for s in streams) + Q_in
        
        # Set initial temperature
        if not T_guess:
            T_guess = sum(s.T for s in streams)/len(streams)
        for s in streams:
            s.T = T_guess
        
        # Find new T
        H = sum(s.H for s in streams)
        T = T_guess + (H_ - H)/C
        
        # Update
        T_guess = T
        for s in streams:
            s.T = T_guess
        if approximate:
            return
        
        # Solve enthalpy by iteration
        it = 1
        while abs(T - T_guess) > 0.01:
            # Calculate
            H_ = sum(s.H for s in streams)
            C = sum(s.C for s in streams)
            T = T_guess + (H_- H)/C
            
            # Check iteration
            if it > 40:
                raise SolverError(f"Could not solve temperature.")
            
            # Update
            T_guess = T
            for s in streams:
                s.T = T_guess
            it += 1

    @staticmethod
    def sum(s_sum, streams, approximate=False):
        """Mixes streams and sets resulting mass and energy in 's_sum'. Assumes the same pressure as streams[0], and no phase equilibrium. If 's_sum' is a Stream object, it assumes the phase does not change.

        **Parameters** 

             **s_sum:** [Stream] container for the resulting mixture of streams

             **streams:** [tuple] Stream objects to be mixed

             **approximate:** [bool] If True, approximate energy balance with constant heat capacity

        **Return**

             **s_sum:** same object as in arguments

        .. Note:: Ignores contents of s_sum as if it was empty

        Example

        >>> s2.T = 340 # Set higher temperature to show energy balance
        >>> s_sum = Stream.sum(Stream('s_sum'), s1, s2)
        >>> s_sum.show()
        Stream: s_sum
         phase: 'l', T: 325.32 K, P: 101325 Pa
         flow (kmol/hr): Ethanol  1
                         Water    4
        """
        # Check if energy balance is required
        T_init = streams[0].T
        Ts = np.array([s.T for s in streams[1:]])
        energy_balance = any(T_init != Ts)
        if energy_balance:
            H = sum([s.H for s in streams])
        else:
            s_sum.T = T_init

        # Copy starting values
        s_sum.copy_like(streams[0])
        streams = streams[1:]

        # Mass balance
        MStream = mixed_stream.MixedStream
        if isinstance(s_sum, MStream):
            # For MixedStream objects
            phases_molflow = mixed_stream.phases_molflow
            solid_mol = s_sum.solid_mol
            liquid_mol = s_sum.liquid_mol
            LIQUID_mol = s_sum.LIQUID_mol
            vapor_mol = s_sum.vapor_mol
            for s in streams:
                if isinstance(s, MStream):
                    solid_mol += s.solid_mol
                    liquid_mol += s.liquid_mol
                    LIQUID_mol += s.LIQUID_mol
                    vapor_mol += s.vapor_mol
                elif isinstance(s, Stream):
                    phase = s.phase
                    attr = phases_molflow[phase]
                    phase_mol = getattr(s_sum, attr)
                    phase_mol += s.mol
        elif isinstance(s_sum, Stream):
            # For Stream objects
            for s in streams:
                s_sum.mol[:] += s.mol
        else:
            raise TypeError('Must pass a valid Stream object')

        # Energy Balance
        if energy_balance:
            if approximate:
                s_sum.T += (H - s_sum.H)/s_sum.C
            else:
                s_sum.H = H

        return s_sum

    # MixedStream compatibility
    def enable_phases(self):
        """Turn stream into a MixedStream object."""
        mS = mixed_stream.MixedStream
        if type(self) is not mS:
            default_phase = self.phase
            self.__class__ = mS
            self.__init__(T=self.T, P=self.P)
            self._default_phase = default_phase
            self.mol = self._mol

    def disable_phases(self, phase):
        """Turn stream into a Stream object."""
        if type(self) is not Stream:
            mol = self.mol
            self.__class__ = Stream
            self.phase = phase
            self._mol = mol

    def VLE(self, specie_IDs=None, LNK=None, HNK=None, P=None,
            T=None, V=None, x=None, y=None, Qin=None):
        self.enable_phases()
        self.VLE(specie_IDs, LNK, HNK, P, T, V, x, y, Qin)
        
    def LLE(self, specie_IDs, split=None, lNK=(), LNK=(),
            P=None, T=None, Qin=None,
            solvents=(), solvent_split=()):
        self.enable_phases()
        self.LLE(specie_IDs, split, lNK, LNK, P, T, Qin,
                 solvents, solvent_split)

    # Representation
    def _info(self, **show_units):
        """Return string with all specifications."""
        units = self.units
        T_units = show_units.get('T')
        P_units = show_units.get('P')
        flow_units = show_units.get('flow')
        fraction = show_units.get('fraction')
        #dim = CS.dim
        
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
        
        # Start of third line (flow rates)
        flow_dim = Q_(0, flow_units).dimensionality
        if fraction:
            if flow_dim == mol_flow_dim:
                flownet = Q_(self.molnet, units['molnet']).to(flow_units).magnitude
                flow = 'molfrac'
            elif flow_dim == mass_flow_dim:
                flownet = Q_(self.massnet, units['massnet']).to(flow_units).magnitude
                flow = 'massfrac'
            elif flow_dim == vol_flow_dim:
                flownet = Q_(self.volnet, units['volnet']).to(flow_units).magnitude
                flow = 'volfrac'
            else:
                raise DimensionError(f"Dimensions for flow units must be in molar, mass or volumetric flow rates, not '{flow_dim}'.")
            beginning = ' flow: '
            end = f'\n*{flownet:.3g} {flow_units} '
        else:
            beginning = f' flow ({flow_units}): '
            end = ''
            if flow_dim == mol_flow_dim:
                flow = 'mol'
            elif flow_dim == mass_flow_dim:
                flow = 'mass'
            elif flow_dim == vol_flow_dim:
                flow = 'vol'
            else:
                raise DimensionError(f"Dimensions for flow units must be in molar, mass or volumetric flow rates, not '{flow_dim}'.")
        # Remaining lines (all flow rates)
        new_line_spaces = len(beginning) * ' '
        nonzero, species = self.nonzero_species
        if 'frac' in flow:
            flow = getattr(self, flow)[nonzero]
        else:
            flow = Q_(getattr(self, flow)[nonzero], units[flow]).to(flow_units).magnitude
        len_ = len(nonzero)
        if len_ == 0:
            return basic_info + ' flows:  0'
        else:
            flowrates = ''
            lengths = [len(sp) for sp in species]
            maxlen = max(lengths) + 1
            for i in range(len_-1):
                spaces = ' ' * (maxlen - lengths[i])
                flowrates += species[i] + spaces + \
                    f' {flow[i]:.3g}\n' + new_line_spaces
            spaces = ' ' * (maxlen - lengths[len_-1])
            flowrates += species[len_-1] + spaces + f' {flow[len_-1]:.3g}'
        return basic_info + beginning + flowrates + end.replace('*', new_line_spaces + 'net' + (maxlen-4)*' ' + '  ')

    def show(self, **show_units):
        """Print all specifications."""
        print(self._info(**show_units))

    def __str__(self):
        return self.ID

    def __repr__(self):
        return '<' + type(self).__name__ + ': ' + self.ID + '>'


from . import mixed_stream


