#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 18 14:05:10 2018

@author: yoelr
"""
from . import _Q
import numpy as np
from scipy.optimize import newton, least_squares
from ._utils import property_array, PropertyFactory, DisplayUnits, \
                    tuple_array, fraction, Sink, Source, missing_stream
from ._flowsheet import find
from ._species import Species
from ._exceptions import SolverError, EquilibriumError, DimensionError
from ._equilibrium import DORTMUND

__all__ = ('Stream',)


# %% TODOs

# TODO: add material property interphase when using Cape Open package

# %% Functions

def nonzero_species(num_IDs, flow):
    index = []
    species = []
    for i, ID in num_IDs:
        if flow[i] != 0:
            index.append(i)
            species.append(ID)
    return index, species

def _print_helpdata(helpdata):
    """Print help data."""
    # Only one helplist, print a nice string
    if isinstance(helpdata[0], str):
        propID, description, dependency, units, datatype = helpdata
        if dependency == 'TP':
            dependency = 'as a function of T and P '
        elif dependency == 'T':
            dependency = 'as a function of T '
        elif dependency == 'P':
            dependency = 'as a function of P '
        print(f"{propID}: [{datatype}] {description.capitalize()} {dependency}({units}).")

    # Many helpdata, print all the nice strings
    else:
        for i in helpdata:
            _print_helpdata(i)


# %% Units of measure

# Biosteam units of measure
units_of_measure = dict(MW='g/mol',
                        mass='kg/hr',
                        mol='kmol/hr',
                        vol='m^3/hr',
                        massnet='kg/hr',
                        molnet='kmol/hr',
                        volnet='m^3/hr',
                        massfrac='kg/kg',
                        molfrac='kmol/kmol',
                        volfrac='m^3/m^3',
                        T='K',
                        P='Pa',
                        H='kJ/hr',
                        S='kJ/hr',
                        G='kJ/hr',
                        U='kJ/hr',
                        A='kJ/hr',
                        Hf='kJ/hr',
                        C='kJ/K/hr',
                        Vm='m^3/mol',
                        Cpm='J/mol/K',
                        Cp='J/g/K',
                        rho='kg/m^3',
                        rhom='mol/m^3',
                        nu='m^2/s',
                        mu='Pa*s',
                        sigma='N/m',
                        k='W/m/K',
                        alpha='m^2/s')

mol_flow_dim = _Q(0, units_of_measure['mol']).dimensionality
mass_flow_dim = _Q(0, units_of_measure['mass']).dimensionality
vol_flow_dim = _Q(0, units_of_measure['vol']).dimensionality

# %% Flow properties

@PropertyFactory
def MassFlow(self):
    """Mass flow (kg/hr)."""
    return self.data[0][0] * self.data[1] # mol[0] * MW

@MassFlow.setter
def MassFlow(self, value):
    self.data[0][0] = value/self.data[1] # mol[0] = value/MW

@PropertyFactory    
def VolumetricFlow(self):
    """Volumetric flow (m^3/hr)."""
    stream, mol = self.data
    m = mol[0]
    if m:
        c = self.name # c = compound
        c.T = stream.T
        c.P = stream.P
        c.phase = stream._phase
        return c.Vm * m * 1000
    else:
        return 0.

@VolumetricFlow.setter
def VolumetricFlow(self, value):
    stream, mol = self.data
    if value:
        c = self.name # c = compound
        c.T = stream.T
        c.P = stream.P
        c.phase = stream._phase
        mol[0] = value/(c.Vm * 1000)
    else:
        mol[0] = 0.

phases = ('s', 'l', 'L', 'g')
phase_index = dict(zip(phases, (0, 1, 2, 3)))

def flow(fget):
    def fset(self, value):
        if fget(self) is not value: raise AttributeError(f"can't set attribute")
    return property(fget, fset)


# %% Stream classes

class metaStream(type):
    """Metaclass for Stream."""
    @property
    def species(cls):
        """[Species] Contains pure component thermodynamic properties for computing overall properties of Stream instances."""
        return cls._cls_species
    @species.setter
    def species(cls, species):
        # Set Species object and related parameters
        if species is cls._cls_species: pass
        elif isinstance(species, Species):
            for cl in (Stream, MS.MixedStream):
                cl._cls_species = species
            species._read_only()
        else: raise ValueError('must pass a Species object')
    _species = species
        

class Stream(metaclass=metaStream):
    """Create a Stream object that defines material flow rates along its thermodynamic state. Thermodynamic and transport properties of a stream are readily available. Ideal mixture is assumed for stream properties and excess thermodynamic energies are neglected as a simplifying assumption for low pressure processes.

**Parameters**

    **ID:** [str] A unique identification. If ID is an empty string (i.e. '' ), a default ID will be chosen. If ID is None, stream will not be registered in flowsheet.

    **flow:** [tuple] All flow rates corresponding to `species`.

    **species:** tuple[str] or [Species] Species corresponding to `flow`. If empty, assume same species as class.

    **units:** [str] The units of the flow rates (only mass, molar, and volumetric flow rates are valid)

    **phase:** [str] 'l' for liquid, 'g' for gas, or 's' for solid.

    **T:** [float] Temperature (K).

    **P:** [float] Pressure (Pa).

    **price:** [float] (USD/kg).
    
    ****flow_pairs:** Compound-flow pairs

**Class Properties**

    **species** = Species(): [Species] Contains pure component thermodynamic properties for computing overall properties of Stream instances.

**Examples**

    Before making a stream, set the species using a Species object:

    .. code-block:: python

       >>> # Set Species object
       >>> Stream.species = Species('Ethanol', 'Water') 

    Stream objects may be created a variety of ways:

    .. code-block:: python

       >>> # Create a stream specifying compound and flow rate pairs:
       >>> s1 = Stream(ID='s1', Water=2)
       >>> s1.show()
       Stream: s1
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Water  2

       >>> # Create a stream assuming same order as given in species:
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
       array([1., 2.])

    .. Warning:: Stream objects do not automatically calculate thermodynamic equilibrium. They simply assume the given phase, temperature and pressure are correct. To find equilibrium, use the VLE or LLE method.

    .. Note:: All following examples will use the Stream objects created here (s1, s2, and s3) starting from its initial values.

    Use the `show` method to print all specifications with desired units:

    .. code-block:: python

       >>> # Temperature in degree Celsius
       >>> s2.show(T='degC')
       Stream: s2 
        phase: 'l', T: 25. degC, P: 101325 Pa
        flow (kmol/hr): Ethanol  1
                        Water    2
       
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
              net      3 kmol/hr

    Flow rates are stored internally as an array in the ‘mol’ attribute.

    .. code-block:: python

       >>> # Set Water flow rate
       >>> s2.mol[1] = 18
       >>> s2.mol # kmol/hr
       array([1, 18])

    Mass and volumetric flow rates are also available as property_arrays of the molar flow rate. As such, they are always up to date with the molar flow rate and altering them also alters the molar flow rate:

    .. code-block:: python

       >>> # Altering mass or volumetric flows alters the molar flow
       >>> s2.vol # m^3/hr
       property_array([0.059, 0.036])
       >>> s2.vol[:] = [1, 1]
       >>> s2.mol
       array([17.06 , 55.343])
       >>> # Values are always up to date with the molar flow
       >>> s2.mol[:] = [1, 2]
       >>> s2.vol
       property_array([0.059, 0.036])

    .. Note::

       property_array objects are significantly slower than array objects. This is because flow rate data is internally stored as molar flow rates. Also, property_array objects are arrays of python objects, which add overhead over the C implementation of numpy. Whenever possible, use the array to manage flow rates. 

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

    A dictionary of available stream properties is available in `Stream.units`. You may also find it useful to use the `help` method to search for a property:
        
    .. code-block:: python

       >>> Stream.help('conductivity')
       k: [float] Thermal conductivity as a function of T and P (W/m/K).

    """
    activity_coefficients = staticmethod(DORTMUND)

    # [dict] Units of measure for material properties (class attribute). 
    units = units_of_measure

    # Information regarding properties
    _prop_info = (
        # ID         # Description               # Dependency # Units      # Type
        ('T',        'temperature',              '',          'K',         'float'),
        ('H',        'enthalpy',                 'T',         'kJ/hr',     'float'),
        ('S',        'entropy',                  'TP',        'kJ/hr',     'float'),
        ('G',        'Gibbs free energy',        'TP',        'kJ/hr',     'float'),
        ('U',        'interal energy',           'TP',        'kJ/hr',     'float'),
        ('A',        'Helmholtz free energy',    'TP',        'kJ/hr',     'float'),
        ('Hf',       'enthalpy of formation',    '',          'kJ/hr',     'float'),
        ('P',        'pressure',                 '',          'Pa',        'float'),
        ('Cpm',      'molar heat capacity',      'T',         'J/mol/K',   'float'),
        ('Cp',       'specific heat capacity',   'T',         'J/kg/K',    'float'),
        ('Vm',       'molar volume',             'TP',        'm^3/mol',   'float'),
        ('rho',      'density',                  'TP',        'kg/m^3',    'float'),
        ('rhom',     'molar density',            'TP',        'mol/m^3',   'float'),
        ('nu',       'kinematic viscosity',      'TP',        'm^2/s',     'float'),
        ('mu',       'hydraulic viscosity',      'TP',        'Pa*s',      'float'),
        ('sigma',    'surface tension',          'T',         'N/m',       'float'),
        ('k',        'thermal conductivity',     'TP',        'W/m/K',     'float'),
        ('alpha',    'thermal diffusivity',      'TP',        'm^2/s',     'float'),
        ('Pr',       'Prantl number',            'TP',        "''",        'float'),
        ('mass',     'mass flow rates',          '',          'kg/hr',     'ndarray'),
        ('mol',      'molar flow rates',         '',          'kmol/hr',   'ndarray'),
        ('vol',      'volumetric flow rates',    'TP',        'm^3/hr',    'ndarray'),
        ('massnet',  'net mass flow rate',       '',          'kg/hr',     'float'),
        ('molnet',   'net molar flow rate',      '',          'kmol/hr',   'float'),
        ('volnet',   'net volumetric flow rate', 'TP',        'm^3/hr',    'float'),
        ('massfrac', 'mass fractions',           '',          'kg/kg',     'ndarray'),
        ('molfrac',  'molar fractions',          '',          'kmol/kmol', 'ndarray'),
        ('volfrac',  'volumetric fractions',     'TP',        'm^3/m^3',   'ndarray'))

    __slots__ = ('T', 'P', '_mol', '_mass', '_vol', 'price', '_ID', '_link',
                 '_species', '_sink', '_source', '_dew_cached', '_bubble_cached',
                 '_phase', '_y_cached', '_lL_split_cached', '_y', '_T_VP', '_T_QinP',
                 '_yP', '_T_dew', '_upstream_connection', '_downstream_connection',
                 '_y_over_Psat_gamma', '__weakref__')

    ### Class attributes for working species ###
    _cls_species = None
    
    # [list] Default starting letter and current number for ID (class attribute)
    _default_ID = ['d', 1]
    
    #: [bool] If True, approximate energy balance. False otherwise.
    lazy_energy_balance = True

    #: [DisplayUnits] Units of measure for IPython display
    display_units = DisplayUnits(T='K', P='Pa',
                                 flow=('kmol/hr', 'kg/hr', 'm3/hr'),
                                 fraction=False)
    line = 'Stream'

    def __init__(self, ID='', flow=(), species=(), units='kmol/hr',
                 phase='l', T=298.15, P=101325, *, price=0, **flow_pairs):
        # Get species and set species information
        if isinstance(species, Species):
            self._species = species
            species._read_only()
        else: 
            self._species = self._cls_species
        self._ID = self._link = self._sink = self._source = None
        self.phase = phase
        self.T = T  #: [float] Temperature (K)
        self.P = P  #: [float] Pressure (Pa)
        self.price = price  #: Price of stream (USD/kg)
        # Dew point cached:
        # (species, pressure, temperature, vapor composition, liquid composition)
        # Bubble point cached:
        # (species, pressure, temperature, vapor composition, liquid composition)
        self._dew_cached = self._bubble_cached = (None,)
        
        # Initialize flows
        try:
            self._setflows(flow, species, flow_pairs)
        except Exception as Error:
            if not self._species: raise RuntimeError(f'must specify Stream.species')
            else: raise Error
        mol = self._mol
        MW = self._species._MW
        mass = [] # Mass flow rates
        vol = [] # Volumetric flow rates    
        for i, s in self._species._num_compounds:
            mol_i = mol[i:i+1]
            mass.append(MassFlow(s.ID, (mol_i, MW[i])))
            vol.append(VolumetricFlow(s, (self, mol_i)))
        self._mass = property_array(mass)
        self._vol = property_array(vol)
        if units == 'kmol/hr': pass
        if units == 'kg/hr': self._mass[:] = mol
        elif units == 'm3/hr': self._vol[:] = mol
        else:
            q = _Q(mol, units)
            dim = q.dimensionality
            if dim == mol_flow_dim:
                self._mol[:] = q.to('kmol/hr').magnitude
            elif dim == mass_flow_dim:
                self._mass[:] = q.to('kg/hr').magnitude
            elif dim == vol_flow_dim:
                self._vol[:] = q.to('m3/hr').magnitude
            else:
                raise DimensionError(f"dimensions for flow units must be in molar, mass or volumetric flow rates, not '{dim}'")
        self.ID = ID


    def setflow(self, flow=(), species=(), units='kmol/hr', inplace='', **flow_pairs):
        """Set `flow` rates according to the `species` order and `flow_pairs`. `inplace` can be any operation that can be performed in place (e.g. +, -, *, /, |, **, etc.)."""
        species = (*species, *flow_pairs.keys())
        flow = (*flow, *flow_pairs.values())
        index = self.indices(*species) if species else ... 
        if units == 'kmol/hr': exec(f"self._mol[index] {inplace}= flow", locals())
        elif units == 'kg/hr': exec(f"self._mass[index] {inplace}= flow", locals())
        elif units == 'm3/hr': exec(f"self._vol[index] {inplace}= flow", locals())
        else:
            q = _Q(flow, units)
            dim = q.dimensionality
            if dim == mol_flow_dim:
                exec(f"self._mol[index] {inplace}= q.to('kmol/hr').magnitude", locals())
            elif dim == mass_flow_dim:
                exec(f"self._mass[index] {inplace}= q.to('kg/hr').magnitude", locals())
            elif dim == vol_flow_dim:
                exec(f"self._vol[index] {inplace}= q.to('m3/hr').magnitude", locals())
            else:
                raise DimensionError(f"dimensions for flow units must be in molar, mass or volumetric flow rates, not '{dim}'")
    
    def getflow(self, *species, units='kmol/hr'):
        """Get flow rates of species in given units."""
        index = self.indices(*species) if species else ...
        if units=='kmol/hr': return self._mol[index]    
        q = _Q(1, units)
        dim = q.dimensionality
        if dim == mol_flow_dim:
            return self._mol[index]*q.to('kmol/hr').magnitude
        elif dim == mass_flow_dim:
            return self._mass[index]*q.to('kg/hr').magnitude
        elif dim == vol_flow_dim:
            return self._vol[index]*q.to('m3/hr').magnitude
        else:
            raise DimensionError(f"dimensions for flow units must be in molar, mass or volumetric flow rates, not '{dim}'")
    
    def _setflows(self, flow, species, flow_pairs):
        """Initialize molar flow rates accoring the species order, and flow_pairs. Instance species do not change."""
        lenflow = len(flow)
        empty_flow = lenflow == 0
        if species is self._species:
            if flow_pairs:
                raise ValueError('cannot specify flow pairs when a Species object is passed')
            elif empty_flow:
                self._mol = np.zeros(species._Nspecies, float)
            elif lenflow == self._species._Nspecies:
                self._mol = np.array(flow, float)
            else:
                raise ValueError('if flow is not empty, length of flow must be equal to length of Species object')
        else:
            lenspecies = len(species)
            empty_species = lenspecies == 0
            N = self._species._Nspecies
            self._mol = np.zeros(N, float)
            if flow_pairs and empty_flow and empty_species:
                self._mol[self.indices(*flow_pairs.keys())] = (*flow_pairs.values(),)
            elif lenspecies == lenflow:
                if flow_pairs:
                    self._mol[self.indices(*species, *flow_pairs.keys())] = (*flow, *flow_pairs.values())
                elif not empty_flow:
                    self._mol[self.indices(*species)] = flow
            elif empty_species and lenflow == N:
                self._mol = np.array(flow, float)
            else:
                raise ValueError('length of flow rates must be equal to length of species')

    # Forward pipping
    def __sub__(self, index):
        if isinstance(index, int):
            return Sink(self, index)
        elif isinstance(index, Stream):
            raise TypeError(f"unsupported operand type(s) for -: '{type(self)}' and '{type(index)}'")
        return index.__rsub__(self)
        
    def __rsub__(self, index):
        if isinstance(index, int):
            return Source(self, index)
        elif isinstance(index, Stream):
            raise TypeError(f"unsupported operand type(s) for -: '{type(self)}' and '{type(index)}'")
        return index.__sub__(self)

    # Backward pipping    
    __pow__ = __sub__
    __rpow__ = __rsub__
    
    @property
    def phase(self):
        """[str] 'l' for liquid, 'g' for gas, or 's' for solid."""
        return self._phase
    @phase.setter
    def phase(self, phase):
        if phase not in ('l', 'g', 's'):
            raise ValueError(f"phase must be either 's', 'l', or 'g'")
        self._phase = phase
    
    @classmethod
    def proxy(cls, ID, link=None):
        """Create a Stream object that serves as a proxy for the `link` stream."""
        self = object.__new__(cls)
        self._link = self._source = self._sink = self._ID = None
        self.ID = ID
        if link: self.link = link
        return self
    
    @property
    def link(self):
        """When another Stream object is set as a link, it will share data with that object."""
        return self._link
    
    @link.setter
    def link(self, stream):
        if not stream: return
        try:
            if self._link is not stream:
                self._species = stream._species
                self._mass = stream._mass
                self._mol = stream._mol
                self._vol = stream._vol
                self._dew_cached = stream._dew_cached
                self._bubble_cached = stream._bubble_cached
                self._link = stream
                self.P = stream.P
                self.T = stream.T
                self._phase = stream.phase
        except Exception as Error:
            if isinstance(stream, Stream):
                raise Error
            else:
                raise TypeError(f"link must be a Stream object or None, not a '{type(stream).__name__}' object.")
              
    
    @property
    def species(self):
        """[Species] Contains pure component thermodynamic properties for computing overall properties of Stream instances."""
        return self._species
    
    @property
    def ID(self):
        """Unique identification (str). If set as '', it will choose a default ID."""
        return self._ID

    @ID.setter
    def ID(self, ID):
        if ID == '':
            # Select a default ID if requested
            letter, number = self._default_ID
            self._default_ID[1] += 1
            num = str(number)
            ID = letter + num
            self._ID = ID
            find.stream[ID] = self
        elif ID and ID != self._ID:
            ID = ID.replace(' ', '_')
            ID_words = ID.split('_')
            if not all([word.isalnum() for word in ID_words]):
                raise ValueError('ID cannot have any special characters')
            self._ID = ID
            find.stream[ID] = self

    @property
    def cost(self):
        """Cost of stream (USD/hr)"""
        return self.price*self.massnet

    @property
    def MW(self):
        """Molecular weight of all species (array g/mol):     

        >>> s2.MW 
        tuple_array([46.06844, 18.01528])
        """
        return self._species._MW.view(tuple_array)

    ### Flow properties ###
    
    def indices(self, *species, CAS=False):
        """Return flow indices of specified species.

        **Parameters**

             ***species:** [iterable] species IDs or CAS numbers.
             
             **CAS:** [bool] True if species are CAS numbers.

        Example

        Indices by ID:
        
        >>> s1.indices('Water', 'Ethanol')
        [1, 0]

        Indices by CAS number:
        
        >>> s1.indices('7732-18-5', '64-17-5', CAS=True):
        [1, 0]

        """
        try:
            lookup = self._species._CAS.index if CAS else self._species._IDs.index
            return [lookup(i) for i in species]
        except ValueError as VE:
            if not self._species:
                raise RuntimeError(f'no species defined in stream')
            if all([isinstance(i, str) for i in species]):
                if CAS:
                    stream_species = self._species._CAS
                    start = 'CAS numbers '
                else:
                    stream_species = self._species._IDs
                    start = 'IDs '
            else:
                raise TypeError(f'species must all be str objects')
            not_included = tuple([i for i in species if i not in stream_species])
            if not_included: raise ValueError(f'{start}{not_included} not defined in stream species')
            else: raise VE
                
    # Molar flow
    @flow
    def mol(self):
        """Array of molar flow rates (kmol/hr):

        >>> s1.mol
        array([0, 2])
        >>> s1.mol = [1, 2]
        >>> s1.mol
        array([1, 2])
        """
        return self._mol

    @property
    def molfrac(self):
        """Array of molar fractions.

        >>> s1.molfrac
        tuple_array([0., 1.])
        """
        return fraction(self._mol).view(tuple_array)

    @property
    def molnet(self):
        """Net molar flow rate (kmol/hr)

        >>> s2.molnet
        3
        """
        return self._mol.sum()

    # Mass flow
    @flow
    def mass(self):
        """Array of mass flow rates (kg/hr)

        >>> s1.mass
        property_array([ 0.   , 36.031])
        """
        return self._mass

    @property
    def massfrac(self):
        """Array of mass fractions.

        >>> s1.massfrac
        tuple_array([0, 1])
        """
        return fraction(self._mol * self._species._MW).view(tuple_array)

    @property
    def massnet(self):
        """Net mass flow rate (kg/hr)

        >>> s2.massnet
        82.099
        """
        return (self._species._MW * self._mol).sum()

    # Volumetric flow
    @flow
    def vol(self):
        """Array of volumetric flow rates as a function of T and P (m^3/hr).

        >>> s2.vol
        property_array([0.059, 0.036])
        
        """
        return self._vol

    @property
    def volfrac(self):
        """Array of volumetric fractions as a function of T and P.

        >>> s2.volfrac
        tuple_array([0.619, 0.381])
        """
        return fraction(self._vol).view(tuple_array)
        
    @property
    def volnet(self):
        """Net volumetric flow rate as a function of T and P (m^3/hr).

        >>> s2.volnet # a liquid stream
        0.09475552896916632
        >>> s3.volnet # a gas stream
        97.32699378482434
        """
        return self._species._propflow('Vm', self._mol, self.T, self.P, self._phase)*1000

    ### Energy flows ###

    # Enthalpy
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
        return self._species._propflow('H', self._mol, self.T, self.P, self._phase)

    @H.setter
    def H(self, H):
        try:
            if self.lazy_energy_balance:
                self.T += (H - self.H)/self.C    
            else:
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
                        raise SolverError(f"could not solve temperature given enthalpy")
        except Exception as Error:
            if (self._mol == 0).all():
                raise ValueError(f"cannot set enthalpy to empty stream")
            else:
                raise Error

    @property
    def Hf(self):
        """Heat of formation flow rate (kJ/hr).

        >>> s1.Hf
        -483640.0
        """
        return (self.mol * [i.Hfm or 0 for i in self._species._compounds]).sum()

    @property
    def Hc(self):
        """Heat of combustion flow rate (kJ/hr).

        >>> s2.Hc
        -1409983.0
        """
        return (self.mol * [i.Hc or 0 for i in self._species._compounds]).sum()
        
    # Entropy
    @property
    def S(self):
        """Entropy flow rate as a function of T and P, excluding formation energies (kJ/hr).

        >>> s1.S # The stream is at the reference state
        0.0
        >>> s1.T = 320
        >>> s1.S
        10.642980522582695
        """
        return self._species._propflow('S', self._mol, self.T, self.P, self._phase)

    # Gibbs free energy
    @property
    def G(self):
        """Gibbs free energy flow rate as a function of T and P, excluding formation energies (kJ/hr).

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
        """Internal energy flow rate as a function of T and P, excluding formation energies (kJ/hr).

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
        """Helmholtz energy flow rate as a function of T and P, excluding formation energies (kJ/hr).

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
        """Heat capacity flow rate as a function of T (kJ/K/hr).

        >>> s2.C
        262.74816631655267
        """
        return self._species._propflow('Cpm', self._mol, self.T, self.P, self._phase)

    # Material properties
    @property
    def Cp(self):
        """Specific heat capacity as a function of T (kJ/kg/K).

        >>> s1.Cp
        4.180597021827335
        """
        return self.Cpm*self.molnet/self.massnet

    @property
    def Cpm(self):
        """Molar heat capacity as a function of T (kJ/kmol/K).

        >>> s1.Cpm # (J/mol/K)
        75.31462591538556
        """
        return self.C/self.molnet

    @property
    def Vm(self):
        """Molar volume as a function of T and P (m^3/mol).

        >>> s1.Vm
        1.8069039870122814e-05
        """
        return self._species._prop('Vm', self._mol, self.T, self.P, self._phase)

    @property
    def rho(self):
        """Density as a function of T (kg/m^3).

        >>> s1.rho
        997.0247522552814
        """
        return self.massnet/self.volnet

    @property
    def rhom(self):
        """Molar density as a function of T and P (mol/m^3).

        >>> s1.rhom
        55343.283715561534
        """
        return self.molnet/self.volnet

    @property
    def nu(self):
        """Kinematic viscosity as a function of T and P (m^2/s).

        >>> s1.nu
        9.154438858391748e-07
        """
        return self.mu/self.rho

    @property
    def mu(self):
        """Hydraulic viscosity as a function of T and P (Pa*s).

        >>> s1.mu
        0.0009127202134824155
        """
        # Katti, P.K.; Chaudhri, M.M. (1964). "Viscosities of Binary Mixtures of Benzyl Acetate with Dioxane, Aniline, and m-Cresol". Journal of Chemical and Engineering Data. 9 (1964): 442–443.
        molfrac = self.molfrac
        props = self._species._props
        mus = np.array(props('mu', self._mol, self.T, self.P, self._phase))
        Vms = np.array(props('mu', self._mol, self.T, self.P, self._phase))
        pos = np.where(molfrac != 0)
        return np.exp((molfrac[pos]*np.log(mus[pos]*Vms[pos])).sum())/self.Vm

    @property
    def k(self):
        """Thermal conductivity as a function of T and P (W/m/k).

        >>> s1.k
        0.5942044328004411
        """
        return self._species._prop('k', self._mol, self.T, self.P, self._phase)

    @property
    def alpha(self):
        """Thermal diffusivity as a function of T and P (m^2/s).

        >>> s1.alpha
        1.4255801521655763e-07
        """
        return self._species._prop('alpha', self._mol, self.T, self.P, self._phase)

    @property
    def sigma(self):
        """Surface tension as a function of T (N/m).

        >>> s1.sigma
        0.07205503890847455
        """
        return self._species._prop('sigma', self._mol, self.T, self.P, self._phase)

    @property
    def Pr(self):
        """Prandtl number as a function of T and P (non-dimensional).

        >>> s1.Pr
        6.421553249380879
        """
        return self._species._prop('Pr', self._mol, self.T, self.P, self._phase)

    # Other properties
    @property
    def source(self):
        """Unit source."""
        return self._source

    @property
    def sink(self):
        """Unit sink."""
        return self._sink

    @property
    def nonzero_species(self):
        """Flow indices and species that have a non-zero flow rate.

        >>> s1.nonzero_species
        [1], ['Water']
        """
        return nonzero_species(self._species._num_IDs, self.mol)
    
    def quantity(self, prop_ID):
        """Return a property as a Quantity object as described in the `pint package <https://pint.readthedocs.io/en/latest/>`__ 

        **Parameters**

             **prop_ID:** [str] name of the property (e.g. 'mol', 'H', 'k' ...)

        Example:

        >>> s1.quantity('Cp')
        <Quantity(4.180597021827335, 'joule / gram / kelvin')>

        """
        attr = getattr(self, prop_ID)
        if isinstance(attr, np.ndarray):
            attr = np.array(attr.astype(float))
        return _Q(attr, self.units[prop_ID])

    # Derivative of a property and sensivity analysis
    def derivative(self, prop_ID, var_ID, index=None):
        """Return the derivative of property 'prop_ID' and variable 'var_ID'. If the property given by var_ID is an array, use index to specify which element.

         **Parameters**

             **prop_ID:** [str] Name of the property (e.g. 'rho', 'Cp', 'H', 'volnet' ...)

             **var_ID:** [str] Name of the variable (e.g 'T', 'P', 'molnet' ...)

             **index:** [array_like or None] Optional indices if the variable is an array

        Example

        >>> s1.derivative('rho', 'T')
        -0.4010588564718449

        """
        if index is not None:
            def getvar():
                return getattr(self, var_ID)[index]

            def setvar(xf):
                getattr(self, var_ID)[index] = xf
        else:
            def getvar():
                return getattr(self, var_ID)

            def setvar(xf):
                setattr(self, var_ID, xf)
        x0 = getvar()
        y0 = getattr(self, prop_ID)
        xf = x0 + 10**-6
        setvar(xf)
        yf = getattr(self, prop_ID)
        setvar(x0)
        return (yf-y0)/(xf-x0)

    # Specifications
    @staticmethod
    def like(stream, ID=''):
        """Create either a Stream or MixedStream object just like the given stream, depending on whether multiple phases are present."""
        s = stream
        if isinstance(stream, MS.MixedStream):
            out = MS.MixedStream(ID,
                                 species=stream._species,
                                 T=s.T, P=s.P)
            out._molarray[:] = stream._molarray
        else:
            out = Stream(ID, flow=s.mol, species=s._species,
                         phase=s.phase, T=s.T, P=s.P)
        return out
    
    def copylike(self, stream):
        """Copy mol, T, P, and phase of stream to self.

        >>> s1.copylike(s3)
        >>> s1.show()
        Stream: s1
         phase: 'g', T: 400.00 K, P: 101325 Pa
         flow (kmol/hr): Ethanol  1
                         Water    2
        """
        if isinstance(stream, MS.MixedStream):
            self.enable_phases()
            self.copylike(stream)
        elif not (self._species is stream._species):
            raise ValueError('species must be the same to copy stream specifications')
        else:
            self._mol[:] = stream.mol
            self.P = stream.P
            self._phase = stream.phase
            self.T = stream.T

    def empty(self):
        """Set flow rates to zero

        >>> s1.empty()
        >>> s1.mol
        array([0, 0])
        """
        self._mol[:] = 0
    
    # Equilibrium
    def _equilibrium_args(self):
        """Return default arguments for equilibrium calculations."""        
        # Choose species that have UNIFAC groups
        species, index = self._species._equilibrium_species(self.mol)
        
        # Get the molar fraction for those species that are not zero
        mol = self._mol[index]
        z = mol/mol.sum()

        return tuple(species), z, self.P
    
    def bubble_point(self):
        """Bubble point at current composition and Pressure.

        **Returns**

            **T_bubble:** [float] The bubble point temperature

            **y:** [numpy ndarray] The composition of the vapor phase considering only the equilibrium species.
             
            **eq_species:** tuple[str] Species in equilibrium.

        >>> stream = Stream(flow=(0.6, 0.4),
        ...                 species=Species('Ethanol', 'Water'))
        >>> stream.bubble_point()
        (352.2820850833474, array([0.703, 0.297]), (<Chemical: Ethanol>, <Chemical: Water>))

        """
        # If just one specie in equilibrium, return Tsat
        species, x, P = self._equilibrium_args()
        N_species = len(species)
        if N_species == 1:
            return species[0].Tsat(P), 1, species
        elif N_species == 0:
            raise EquilibriumError('no species available for phase equilibrium')
        
        # Solve and return bubble point
        x, P = self._bubble_point(species, x, P)
        return x, P, species
    
    def _bubble_point_error(self, T, species, x, P):
        """Bubble point given x and P"""
        gamma = self.activity_coefficients(species, x, T)
        Psat = [s.VaporPressure(T) for s in species]
        self._yP =  x * Psat * gamma
        return 1 - self._yP.sum()/P
    
    def _bubble_point(self, species, x, P):
        """Bubble point at given composition and Pressure

        **Parameters**

            **species:** tuple[Compound] Species corresponding to x.

            **x:** [array_like] The composition of the liquid phase.

            **P:** [float] The pressure in Pascals.
        
        **Returns**

            **T_bubble:** [float] The bubble point temperature

            **y:** [numpy ndarray] The composition of the vapor phase considering only the equilibrium species.
             
            **eq_species:** tuple[str] Species in equilibrium.

        >>> s1._bubble_point(species=Species('Ethanol', 'Water'),
        ...                 x=(0.6, 0.4), P=101325)
        (352.2820850833474, array([0.703, 0.297]))
        
        """
        # Setup functions to calculate vapor pressure and activity coefficients
        x = np.array(x)
        x /= x.sum()
        
        # Retrive cached info
        # Even if different composition, use previous bubble point as guess
        if self._bubble_cached[0] == species:
            _, cP, T_bubble, y_bubble, cx = self._bubble_cached # c means cached
            if abs(x - cx).sum() < 1e-6 and abs(cP - P) < 1:
                # Return cached data
                return T_bubble, y_bubble
        else:
            # Initial temperature guess
            T_bubble = (x * [s.Tb for s in species]).sum()

        # Solve and return
        T_bubble = newton(self._bubble_point_error, T_bubble,
                          args=(species, x, P), tol= 1e-7)
        y_bubble = self._yP / P
        y_bubble /= y_bubble.sum()
        self._bubble_cached = (species, P, T_bubble, y_bubble, x)
        return T_bubble, y_bubble

    def dew_point(self):
        """Dew point at current composition and Pressure.
        
        **Returns**

            **T_dew:** [float] The dew point temperature

            **x:** [numpy array] The composition of the liquid phase considering only the equilibrium species.
             
            **eq_species:** tuple[str] Species in equilibrium.

        >>> stream = Stream(flow=(0.5, 0.5),
        ...                 species=Species('Ethanol', 'Water'))
        >>> stream.dew_point()
        (353.7420742149825, array([0.556, 0.444]), (<Chemical: Ethanol>, <Chemical: Water>))
        
        """
        # If just one specie in equilibrium, return Tsat
        species, y, P = self._equilibrium_args()
        N_species = len(species)
        if N_species == 1:
            return species[0].Tsat(P), 1, species
        elif N_species == 0:
            raise EquilibriumError('no species available for phase equilibrium')
        
        # Solve and return dew point
        y, P = self._dew_point(species, y, P)
        return y, P, species
    
    def _dew_point_error_T(self, T, x, y, P, species):
            """Error for dew point temperature, constant x and P."""
            Psat = np.array([s.VaporPressure(T) for s in species])
            gamma = self.activity_coefficients(species, x, T)
            self._y_over_Psat_gamma = y/(Psat*gamma)
            return 1 - self._y_over_Psat_gamma.sum() * P
    
    def _dew_point_error_x(self, x, y, P, species):
        """Error function for dew point composition with constant y and P."""
        x_dew = x/x.sum()
        self._T_dew = newton(self._dew_point_error_T, self._T_dew,
                             args=(x, y, P, species), tol= 1e-7)
        return abs(x_dew - self._y_over_Psat_gamma*P)
    
    def _dew_point(self, species, y, P):
        """Dew point given y and P.

        **Parameters**

            **species:** tuple[Compound] Species corresponding to x.

            **y:** [array_like] The composition of the vapor phase.

            **P:** [float] The pressure in Pascals.

        **Returns**

            **T_dew:** [float] The dew point temperature

            **x:** [numpy array] The composition of the liquid phase considering only the equilibrium species.
             
            **eq_species:** tuple[str] Species in equilibrium.


        >>> s1.dew_point(species=Species('Ethanol', 'Water'),
        ...              y=(0.5, 0.5), P=101325)
        (353.7420742149825, array([0.556, 0.444]))
        """
        y = np.array(y)
        y /= y.sum()
        
        # Retrive cached info
        # Even if different composition, use previous bubble point as guess
        if self._dew_cached[0] == species:
            _, cP, self._T_dew, cy, x_dew = self._dew_cached # c means cached
            if abs(y - cy).sum() < 1e-6 and abs(cP - P) < 1:
                # Return cached data
                return self._T_dew, x_dew 
        else:
            # Initial temperature guess
            self._T_dew = (y * [s.Tb for s in species]).sum()
            x_dew = y

        # Get min and max splits
        Nspecies = len(species)
        min_ = np.zeros(Nspecies)
        max_ = np.ones(Nspecies)

        # Find x by least squares
        try: x_dew = least_squares(self._dew_point_error_x, x_dew,
                                   bounds=(min_, max_),
                                   args=(y, P, species),
                                   xtol= 1e-6).x
        except: pass
        self._dew_cached = (species, P, self._T_dew, y, x_dew)
        return self._T_dew, x_dew

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
        volnet = self.volnet/3600 # m3/s
        if A is None:  # Assume pipe
            return 4*volnet/(self.nu*np.pi*L)
        return volnet*L/(self.nu*A)

    # Class methods
    @classmethod
    def _helpdata(cls, prop):
        """Return information related to a property.

        **Parameters**

             **prop:** [str] name or description of a property

        Return

             **out:** tupple[str] where elements are:
                  * **ID:** [str] name of the property

                  * **Description:** [str] description of the property

                  * **Dependency:** [str] TP dependency

                  * **Units:** [str] units of measure

                  * **Datatype:** [str] type returned by property

        Example

        >>> Stream._helplist('k')
        ('k', 'thermal conductivity', 'TP', 'W/m/K', 'float')

        >>> Stream._helplist('viscosity')
        (('nu', 'kinematic viscosity', 'TP', 'm^2/s', 'float'),
         ('mu', 'hydraulic viscosity', 'TP', 'Pa*s', 'float'))

        >>> Stream._helplist('kmol/hr')
        (('mol', 'molar flow rates', '', 'kmol/hr', 'np.array'),
         ('molnet', 'net molar flow rate', '', 'kmol/hr', 'float'))
        """
        # Compare by both ID and description
        for index in (0, 1):
            # Search for exact matches
            for l in cls._prop_info:
                if prop == l[index]: return l
        
        out = []
        for index in (0, 1):    
            # If no matches found, search harder
            for l in cls._prop_info:
                if prop.lower() in l[index].lower(): out.append(l)
        return out

    @classmethod
    def help(cls, prop):
        """Print information related to a property.

        **Parameters**

             **prop:** [str] name or description of a property     

        Example

        >>> Stream.help('rho')
        rho: [float] Density as a function of T and P (kg/m^3).
        >>> Stream.help('density')
        rho: [float] Density as a function of T and P (kg/m^3).
        """
        data = cls._helpdata(prop)
        if data: _print_helpdata(data)
        else: print(f"No matching property '{prop}'.")

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
        sum_ = sum
        # Set up problem
        C = sum_([s.C for s in streams])
        H_ = sum_([s.H for s in streams]) + Q_in
        
        # Set initial temperature
        if not T_guess:
            T_guess = sum_([s.T for s in streams])/len(streams)
        for s in streams:
            s.T = T_guess
        
        # Find new T
        H = sum_([s.H for s in streams])
        T = T_guess + (H_ - H)/C
        
        # Update
        T_guess = T
        for s in streams:
            s.T = T_guess
        if approximate: return
        
        # Solve enthalpy by iteration
        it = 1
        while abs(T - T_guess) > 0.01:
            # Calculate
            H_ = sum_([s.H for s in streams])
            T = T_guess + (H_- H)/C
            
            # Check iteration
            if it > 40:
                raise SolverError(f"could not solve temperature")
            
            # Update
            T_guess = T
            for s in streams:
                s.T = T_guess
            it += 1

    @staticmethod
    def sum(s_sum, streams):
        """Mixes streams and sets resulting mass and energy in 's_sum'. Assumes the same pressure as streams[0], and no phase equilibrium. If 's_sum' is a Stream object, it assumes the phase does not change.

        **Parameters** 

             **s_sum:** [Stream] Container for the resulting mixture of streams

             **streams:** [tuple] Stream objects to be mixed
             
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
        other_streams = streams[1:]
        energy_balance = any([T_init!=s.T for s in other_streams])
        if energy_balance: H = sum([s.H for s in streams])

        # Copy starting values
        s_sum.copylike(streams[0])

        # Mass balance
        inst = isinstance
        MStream = MS.MixedStream
        try:
            if inst(s_sum, MStream):
                # For MixedStream objects
                for s in other_streams:
                    if inst(s, MStream):
                        s_sum._mol += s._mol
                    else:
                        s_sum._mol[phase_index[s._phase]] += s._mol
            else:
                # For Stream objects
                for s in other_streams:
                    s_sum._mol[:] += s.mol
        except Exception as Error:
            if not inst(s_sum, Stream):
                raise TypeError(f's_sum must be a Stream object, not {type(s_sum).__name__}')
            for s in streams:
                if not inst(s, Stream):
                    raise TypeError('streams must only contain Stream objects, not {type(s).__name__} objects')
            raise Error

        # Energy Balance
        if energy_balance: s_sum.H = H
        return s_sum

    # MixedStream compatibility
    def enable_phases(self):
        """Cast stream into a MixedStream object."""
        mol = self._mol
        self.__class__ = MS.MixedStream
        self._setflows(np.zeros((4, self._species._Nspecies)))
        self._mol[phase_index[self._phase]] = mol
        self._T_QinP = self._T_VP = None
        self._lL_split_cached = self._y_cached = (None,)

    def disable_phases(self, phase):
        """Cast stream into a Stream object.
        
        **Parameters**
        
            **phase:** {'s', 'l', 'g'} desired phase of stream
            
        """
        self._phase = phase
            
    def VLE(self, species_IDs=None, LNK=None, HNK=None, P=None,
            Qin=None, T=None, V=None, x=None, y=None):
        self.enable_phases()
        self.VLE(species_IDs, LNK, HNK, P, Qin, T, V, x, y)
        
    def LLE(self, species_IDs=None, split=None, lNK=(), LNK=(),
            solvents=(), solvent_split=(),
            P=None, T=None, Qin=None):
        self.enable_phases()
        self.LLE(species_IDs, split, lNK, LNK,
                 solvents, solvent_split, P, T, Qin)

    def _info_header(self):
        """Return stream information header."""
        # First line
        unit = self._source
        if unit is None:
            source = ''
        else:
            source = f'  from  {type(unit).__name__}-{unit}'
        unit = self._sink
        if unit is None:
            sink = ''
        else:
            sink = f'  to  {type(unit).__name__}-{unit}'
        if self.ID:
            return f"{type(self).__name__}: {self.ID}{source}{sink}"
        else:
            return f"{type(self).__name__}{source}{sink}"
        
    def _info_phaseTP(self, phases, T_units, P_units):
        T = _Q(self.T, self.units['T']).to(T_units).magnitude
        P = _Q(self.P, self.units['P']).to(P_units).magnitude
        return f" phase: '{phases}', T: {T:.5g} {T_units}, P: {P:.6g} {P_units}\n"

    # Representation
    def _info(self, T, P, flow, fraction):
        """Return string with all specifications."""
        units = self.units
        basic_info = self._info_header() + '\n'
        T_units, P_units, flow_units, fraction = [(i if i is not None else j) for i, j in
                                                  zip((T, P, flow, fraction), self.display_units)]
        basic_info += self._info_phaseTP(self._phase, T_units, P_units)
        
        # Start of third line (flow rates)
        flow_dim = _Q(0, flow_units).dimensionality
        if fraction:
            if flow_dim == mol_flow_dim:
                flownet = _Q(self.molnet, units['molnet']).to(flow_units).magnitude
                flow = 'molfrac'
            elif flow_dim == mass_flow_dim:
                flownet = _Q(self.massnet, units['massnet']).to(flow_units).magnitude
                flow = 'massfrac'
            elif flow_dim == vol_flow_dim:
                flownet = _Q(self.volnet, units['volnet']).to(flow_units).magnitude
                flow = 'volfrac'
            else:
                raise DimensionError(f"dimensions for flow units must be in molar, mass or volumetric flow rates, not '{flow_dim}'")
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
                raise DimensionError(f"dimensions for flow units must be in molar, mass or volumetric flow rates, not '{flow_dim}'")
        # Remaining lines (all flow rates)
        new_line_spaces = len(beginning) * ' '
        nonzero, species = self.nonzero_species
        if 'frac' in flow:
            flow = getattr(self, flow)[nonzero]
        else:
            flow = _Q(getattr(self, flow)[nonzero], units[flow]).to(flow_units).magnitude
        len_ = len(nonzero)
        if len_ == 0:
            return basic_info + ' flow: 0'
        else:
            flowrates = ''
            lengths = [len(sp) for sp in species]
            maxlen = max(lengths) + 1
            for i in range(len_-1):
                spaces = ' ' * (maxlen - lengths[i])
                flowrates += species[i] + spaces + f' {flow[i]:.3g}\n' + new_line_spaces
            spaces = ' ' * (maxlen - lengths[len_-1])
            flowrates += species[len_-1] + spaces + f' {flow[len_-1]:.3g}'
        
        return (basic_info 
              + beginning
              + flowrates
              + end.replace('*', new_line_spaces + 'net' + (maxlen-4)*' ' + '  ')
              + (f'\n link: {self._link}' if self._link else ''))

    def _ipython_display_(self, T=None, P=None, flow=None, fraction=None):
        """Print all specifications."""
        print(self._info(T, P, flow, fraction))
    show = _ipython_display_

    def _delete(self):
        source = self._source and self._source._outs
        if source: source[source.index(self)] = missing_stream
        sink = self._sink and self._sink._ins
        if sink: sink[sink.index(self)] = missing_stream

    def delete(self):
        delattr(find.stream, self._ID)

    def __str__(self):
        if self.ID: return self.ID
        else: return type(self).__name__

    def __repr__(self):
        if self.ID: return f'<{type(self).__name__}: {self.ID}>'
        else: return f'<{type(self).__name__}>'
        
        
from . import _mixed_stream as MS