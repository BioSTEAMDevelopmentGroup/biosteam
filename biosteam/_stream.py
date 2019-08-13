#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 18 14:05:10 2018

@author: yoelr
"""
from . import _Q
import numpy as np
from ._utils import property_array, PropertyFactory, DisplayUnits, \
                    tuple_array, fraction, Sink, Source, MissingStream
from ._flowsheet import find
from ._species import Species, WorkingSpecies
from ._exceptions import SolverError, EquilibriumError, DimensionError
from ._equilibrium import Dortmund, VLE, BubblePoint, DewPoint


__all__ = ('Stream',)

# %% TODOs

# TODO: add material property interphase when using Cape Open package

# %% Functions

def nonzero_species(species, flow):
    index_ = []
    IDs_ = []
    IDs = species._IDs
    for i in species._index:
        if flow[i] != 0:
            index_.append(i)
            IDs_.append(IDs[i])
    return index_, IDs_

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
        print(f"{propID}: [{datatype}] {description.capitalize()} "
              "{dependency}({units}).")

    # Many helpdata, print all the nice strings
    else:
        for i in helpdata:
            _print_helpdata(i)


# %% Units of measure

# Biosteam units of measure
units_of_measure = dict(cost='USD/hr',
                        MW='g/mol',
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
        if isinstance(species, Species):
            Stream._cls_species = WorkingSpecies(species)
        elif isinstance(species, WorkingSpecies):
            Stream._cls_species = species
        else: raise ValueError('must pass a Species object')
    _species = species
    @property
    def indices(self):
        return self._cls_species.indices
    @property
    def index(self):
        return self._cls_species.index
    @property
    def MW(cls):
        return cls._MW

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
    
       >>> s2.T = 298.15
       >>> s2.Cp # Heat capacity (kJ/(kg-K))
       3.200382054794244
       >>> s2.H  # Enthalpy (kJ/hr)
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
    
       >>> s1.volnet # Volumetric flow rate (m^3/hr)
       0.036138079740245625
       >>> s1.rho # Density (kg/m^3)
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

    A dictionary of available stream properties and respective units of measure is available in `Stream.units`. You may also find it useful to use the `help` method to search for a property:
        
    .. code-block:: python

       >>> Stream.help('conductivity')
       k: [float] Thermal conductivity as a function of T and P (W/m/K).

    """
    
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
                 '_species', '_sink', '_source', '_dew_point',
                 '_bubble_point', '_gamma', '_phase', '_lL_split_cached',
                 '__weakref__', '_source_link', '_VLE')

    line = 'Stream'

    ### Class attributes for working species ###    
    _cls_species = _MW = None
    
    # [list] Default starting letter and current number for ID (class attribute)
    _default_ID = ['d', 1]
    
    #: [bool] If True, approximate energy balance. False otherwise.
    lazy_energy_balance = True

    #: [DisplayUnits] Units of measure for IPython display
    display_units = DisplayUnits(T='K', P='Pa',
                                 flow=('kmol/hr', 'kg/hr', 'm3/hr'),
                                 fraction=False)

    def __init__(self, ID='', flow=(), species=(), units='kmol/hr',
                 phase='l', T=298.15, P=101325, *, price=0, **flow_pairs):
        # Get species and set species information
        if isinstance(species, Species):
            self._species = WorkingSpecies(species)
            species = ()
        elif isinstance(species, WorkingSpecies):
            self._species = species
            species = ()
        else: 
            assert self._cls_species, 'must define Stream.species first'
            self._species = self._cls_species
        self._link = self._ID = self._sink = self._source = None
        self._source_link = self
        self.phase = phase
        self.T = T  #: [float] Temperature (K)
        self.P = P  #: [float] Pressure (Pa)
        self.price = price  #: Price of stream (USD/kg)
        
        # Initialize flows
        self._setflows(flow, species, flow_pairs)
        mol = self._mol
        MW = self._species._MW
        mass = [] # Mass flow rates
        vol = [] # Volumetric flow rates    
        cmps = self._species._compounds
        for i in self._species._index:
            mol_i = mol[i:i+1]
            s = cmps[i]
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
        self._gamma = gamma = Dortmund()
        self._bubble_point = BubblePoint(gamma)
        self._dew_point = DewPoint(gamma)

    def setflow(self, flow=(), species=(), units='kmol/hr', inplace='', **flow_pairs):
        """Set `flow` rates according to the `species` order and `flow_pairs`. `inplace` can be any operation that can be performed in place (e.g. +, -, *, /, |, **, etc.)."""
        species = (*species, *flow_pairs.keys())
        flow = (*flow, *flow_pairs.values())
        index = self.indices(species) if species else ... 
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
        index = self.indices(species) if species else ...
        q = _Q(1, units)
        dim = q.dimensionality
        if dim == mol_flow_dim:
            return self._mol[index]*q.to('kmol/hr').magnitude
        elif dim == mass_flow_dim:
            return self._mass[index]*q.to('kg/hr').magnitude
        elif dim == vol_flow_dim:
            return self._vol[index]*q.to('m3/hr').magnitude
        else:
            raise DimensionError(f"dimensions for flow units must be in molar, "
                                 f"mass or volumetric flow rates, not '{dim}'")
    
    def _setflows(self, flow, species, flow_pairs):
        """Initialize molar flow rates according to the species order, and flow_pairs. Instance species do not change."""
        flowlen = len(flow)
        specieslen = len(species)
        if flowlen:
            if flow_pairs:
                raise ValueError('cannot specify flow pairs when species is passed')
            elif flowlen == specieslen:
                self._mol = self._species.array(species, flow)
            elif (not specieslen) and (flowlen == self._species._N):
                self._mol = np.array(flow, float)
            else:
                raise ValueError('length of flow rates must be equal to length of species')
        elif flow_pairs:
            self._mol = self._species.array(flow_pairs, [*flow_pairs.values()])
        else:
            self._mol = np.zeros(self.species._N, float)
                
    # Forward pipping
    def __sub__(self, index):
        if isinstance(index, int):
            return Sink(self, index)
        elif isinstance(index, Stream):
            raise TypeError("unsupported operand type(s) for -: "
                            f"'{type(self)}' and '{type(index)}'")
        return index.__rsub__(self)
        
    def __rsub__(self, index):
        if isinstance(index, int):
            return Source(self, index)
        elif isinstance(index, Stream):
            raise TypeError("unsupported operand type(s) for -: "
                            "'{type(self)}' and '{type(index)}'")
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
    
    @staticmethod
    def proxy(ID, link=MissingStream):
        """Create a Stream object that serves as a proxy for its `link` stream."""
        self = object.__new__(Stream)
        self._source = self._sink = self._ID = None
        self._source_link = MissingStream
        self.ID = ID
        self.price = 0
        if link: self.link = link
        else: self._link = MissingStream
        return self
    
    @property
    def link(self):
        """When another Stream object is set as a link, it will share data with that object."""
        return self._link
    
    @link.setter
    def link(self, stream):
        try:
            if self._source_link is stream._source_link:
                self._link = stream
            elif stream is None:
                self.__init__(ID=self._ID, flow=self._mol,
                              species=self._species,
                              T=self.T, P=self.P, phase=self._phase)
            else:
                self._species = stream._species
                self._mass = stream._mass
                self._mol = stream._mol
                self._vol = stream._vol
                self._dew_point = stream._dew_point
                self._bubble_point = stream._bubble_point
                self._gamma = stream._gamma
                self._source_link = stream._source_link
                self._link = stream
                self.P = stream.P
                self.T = stream.T
                self._phase = stream.phase
        except Exception as Error:
            if isinstance(stream, Stream): raise Error
            else: raise TypeError(f"link must be a Stream object, not a "
                                  "'{type(stream).__name__}' object.")
    
    # def __getattr__(self, key):
    #     if self._link is MissingStream:
    #         if self.source and self.source._ins[0] is not MissingStream:
    #             self.source._link_streams()
    #             try: return object.__getattribute__(self, key)
    #             except: raise AttributeError(f'{repr(self.source)}.ins[0] is missing stream')
    #         raise AttributeError(f'{self}.link is missing stream')
    #     elif not self._source_link:
    #         raise AttributeError(f'{repr(self.link.source)} must be simulated first')
    #     else:
    #         raise AttributeError(f"'{type(self).__name__}' has no attribute '{key}'")
    
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
            # Select a default ID
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

        from biosteam import *
        >>> Stream.species = Species('Water', 'Ethanol')
        >>> Stream.MW
        tuple_array([18.01528, 46.06844])
        """
        return self._species._MW.view(tuple_array)

    ### Flow properties ###
    
    @property
    def indices(self):
        """Return flow indices of specified species.

        **Parameters**

             **IDs:** [iterable] species IDs or CAS numbers.

        Example

        Indices by ID:
        
        >>> from biosteam import *
        >>> Stream.species = Species(['Ethanol', 'Water'])
        >>> s1 = Stream()
        >>> s1.indices(['Water', 'Ethanol'])
        [1, 0]

        Indices by CAS number:
        
        >>> s1.indices(['7732-18-5', '64-17-5']):
        [1, 0]

        """
        return self._species.indices
    
    @property
    def index(self):
        """Return index of specified compound.

        **Parameters**

             **ID:** [str] compound ID

        Example

        Index by ID:
        
        >>> from biosteam import *
        >>> Stream.species = Species(['Ethanol', 'Water'])
        >>> s1 = Stream()
        >>> s1.index('Water')
        1

        Indices by CAS number:
        
        >>> s1.index('7732-18-5'):
        1

        """
        return self._species.index
    
    # Molar flow
    @flow
    def mol(self):
        """Array of molar flow rates (kmol/hr):

        >>> from biosteam import *
        >>> Stream.species = Species('Water', 'Ethanol')
        >>> s1 = Stream(Ethanol=2)
        >>> s1.mol
        array([0, 2])
        >>> s1.mol[:] = [1, 2]
        >>> s1.mol
        array([1, 2])
        """
        return self._mol

    @property
    def molfrac(self):
        """Array of molar fractions.

        >>> from biosteam import *
        >>> Stream.species = Species('Water', 'Ethanol')
        >>> s1 = Stream(Ethanol=2)
        >>> s1.molfrac
        tuple_array([0., 1.])
        """
        return fraction(self._mol).view(tuple_array)

    @property
    def molnet(self):
        """Net molar flow rate (kmol/hr)

        >>> from biosteam import *
        >>> Stream.species = Species('Water', 'Ethanol')
        >>> s1 = Stream(Ethanol=2, Water=1)
        >>> s1.molnet
        3
        """
        return self._mol.sum()

    # Mass flow
    @flow
    def mass(self):
        """Array of mass flow rates (kg/hr)

        >>> from biosteam import *
        >>> Stream.species = Species('Ethanol', 'Water')
        >>> s1 = Stream(Water=2)
        >>> s1.mass
        property_array([ 0.   , 36.031])
        """
        return self._mass

    @property
    def massfrac(self):
        """Array of mass fractions.

        >>> from biosteam import *
        >>> Stream.species = Species('Ethanol', 'Water')
        >>> s1 = Stream(Water=2)
        >>> s1.massfrac
        tuple_array([0, 1])
        """
        return fraction(self._mol * self._species._MW).view(tuple_array)

    @property
    def massnet(self):
        """Net mass flow rate (kg/hr)

        >>> from biosteam import *
        >>> Stream.species = Species('Ethanol', 'Water')
        >>> s1 = Stream(Water=2)
        >>> s1.mass
        36.031
        """
        return (self._species._MW * self._mol).sum()

    # Volumetric flow
    @flow
    def vol(self):
        """Array of volumetric flow rates as a function of T and P (m^3/hr).

        >>> from biosteam import *
        >>> Stream.species = Species('Ethanol', 'Water')
        >>> s1 = Stream(Water=2)
        >>> s1.vol
        property_array([0.   , 0.036])
        
        """
        return self._vol

    @property
    def volfrac(self):
        """Array of volumetric fractions as a function of T and P.

        >>> from biosteam import *
        >>> Stream.species = Species('Ethanol', 'Water')
        >>> s1 = Stream(Water=2)
        >>> s1.volfrac
        tuple_array([0.0, 1.0], dtype=object)
                """
        return fraction(self._vol).view(tuple_array)
        
    @property
    def volnet(self):
        """Net volumetric flow rate as a function of T and P (m^3/hr).

        from biosteam import *
        >>> Stream.species = Species('Ethanol', 'Water')
        >>> s1 = Stream(Water=2)
        >>> s1.volnet
        0.036138079740245625
        """
        return self._species._propflow('Vm', self._mol, self.T, self.P, self._phase)*1000

    ### Energy flows ###

    def H_at(self, mol=None, T=None, P=None, phase=None):
        """Return enthalpy flow rate at given arguments, excluding formation energies (kJ/hr). Arguments with None values default to stream specifications."""
        return self._species._propflow('H', mol or self._mol, T or self.T,
                                       P or self.P, phase or self._phase)

    # Enthalpy
    @property
    def H(self):
        """Enthalpy flow rate as a function of T, excluding formation energies (kJ/hr).

        >>> from biosteam import *
        >>> Stream.species = Species('Ethanol', 'Water')
        >>> s1 = Stream(Water=2)
        >>> s1.H # The stream is at the reference state
        0.0
        >>> s1.H = 1000
        >>> s1.T
        304.7888167493753
        >>> s1.H
        999.5728716085112

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
                it = 0
                it2 = 0
                while abs(self.T - T_old) > 0.01:
                    T_old = self.T
                    self.T += (H - self.H)/C
                    if it == 5:
                        it = 0
                        it2 += 1
                        C = self.C
                        if it2 > 10:
                            raise SolverError("could not solve temperature "
                                              "given enthalpy")
                    else: it += 1
        except Exception as Error:
            if (self._mol == 0).all():
                raise ValueError(f"cannot set enthalpy to empty stream")
            else:
                raise Error

    @property
    def Hnet(self):
        """Total enthaly flow rate including heats of formation"""
        return self.Hf + self.H

    @property
    def Hf(self):
        """Heat of formation flow rate (kJ/hr).

        >>> from biosteam import *
        >>> Stream.species = Species('Ethanol', 'Water')
        >>> s1 = Stream(Water=2)
        >>> s1.Hf
        -483640.0
        """
        return (self.mol * [i.Hf or 0 for i in self._species._compounds]).sum()

    @property
    def Hc(self):
        """Heat of combustion flow rate (kJ/hr).

        >>> from biosteam import *
        >>> Stream.species = Species('Ethanol', 'Water')
        >>> s1 = Stream(Ethanol=2)
        >>> s1.Hc
        -2819966.0
        """
        return (self.mol * [i.Hc or 0 for i in self._species._compounds]).sum()
        
    # Entropy
    @property
    def S(self):
        """Entropy flow rate as a function of T and P, excluding formation energies (kJ/hr).

        >>> from biosteam import *
        >>> Stream.species = Species('Ethanol', 'Water')
        >>> s1 = Stream(Ethanol=2)
        >>> s1.S # The stream is at the reference state
        0.0
        >>> s1.T = 320
        >>> s1.S
        16.503899351682694
        """
        return self._species._propflow('S', self._mol, self.T, self.P, self._phase)

    # Gibbs free energy
    @property
    def G(self):
        """Gibbs free energy flow rate as a function of T and P, excluding formation energies (kJ/hr).

        >>> from biosteam import *
        >>> Stream.species = Species('Ethanol', 'Water')
        >>> s1 = Stream(Ethanol=2)
        >>> s1.G # The stream is at the reference state
        0.0
        >>> s1.T = 320
        >>> s1.G
        -182.43024800100193
        """
        return self.H - self.S*self.T

    # Internal energy
    @property
    def U(self):
        """Internal energy flow rate as a function of T and P, excluding formation energies (kJ/hr).

        >>> from biosteam import *
        >>> Stream.species = Species('Ethanol', 'Water')
        >>> s1 = Stream(Ethanol=2)
        >>> s1.U # The stream is at the reference state
        0.0
        >>> s1.T = 320
        >>> s1.U
        -7090.732710112934
        """
        return self.H - self.P*self.volnet

    # Helmholtz
    @property
    def A(self):
        """Helmholtz energy flow rate as a function of T and P, excluding formation energies (kJ/hr).

        >>> from biosteam import *
        >>> Stream.species = Species('Ethanol', 'Water')
        >>> s1 = Stream(Ethanol=2)
        >>> s1.A # The stream is at the reference state
        0.0
        >>> s1.T = 320
        >>> s1.A
        -12371.980502651397
        """
        return self.U - self.T*self.S

    # Capacity flow rate
    @property
    def C(self):
        """Heat capacity flow rate as a function of T (kJ/K/hr).

        >>> from biosteam import *
        >>> Stream.species = Species('Ethanol', 'Water')
        >>> s1 = Stream(Ethanol=2)
        >>> s2.C
        243.41704115214097
        """
        return self._species._propflow('Cpm', self._mol, self.T, self.P, self._phase)

    # Material properties
    @property
    def Cp(self):
        """Specific heat capacity as a function of T (kJ/kg/K).

        >>> from biosteam import *
        >>> Stream.species = Species('Ethanol', 'Water')
        >>> s1 = Stream(Ethanol=2)
        >>> s1.Cp
        2.641906706110962
        """
        return self.Cpm*self.molnet/self.massnet

    @property
    def Cpm(self):
        """Molar heat capacity as a function of T (kJ/kmol/K).

        >>> from biosteam import *
        >>> Stream.species = Species('Ethanol', 'Water')
        >>> s1 = Stream(Ethanol=2)
        >>> s1.Cpm # (J/mol/K)
        121.70852057607048
        """
        return self.C/self.molnet

    @property
    def Vm(self):
        """Molar volume as a function of T and P (m^3/mol).

        >>> from biosteam import *
        >>> Stream.species = Species('Ethanol', 'Water')
        >>> s1 = Stream(Ethanol=2)
        >>> s1.Vm
        0.00012030150757118572
        """
        return self._species._prop('Vm', self._mol, self.T, self.P, self._phase)

    @property
    def rho(self):
        """Density as a function of T (kg/m^3).

        >>> from biosteam import *
        >>> Stream.species = Species('Ethanol', 'Water')
        >>> s1 = Stream(Ethanol=2)
        >>> s1.rho
        765.8830039638536
        """
        return self.massnet/self.volnet

    @property
    def rhom(self):
        """Molar density as a function of T and P (mol/m^3).

        >>> from biosteam import *
        >>> Stream.species = Species('Ethanol', 'Water')
        >>> s1 = Stream(Ethanol=2)
        >>> s1.rhom
        16.624895567634884
        """
        return self.molnet/self.volnet

    @property
    def nu(self):
        """Kinematic viscosity as a function of T and P (m^2/s).

        >>> from biosteam import *
        >>> Stream.species = Species('Ethanol', 'Water')
        >>> s1 = Stream(Ethanol=2)
        >>> s1.nu
        5.733020843964377e-06
        """
        return self.mu/self.rho

    @property
    def mu(self):
        """Hydraulic viscosity as a function of T and P (Pa*s).

        >>> from biosteam import *
        >>> Stream.species = Species('Ethanol', 'Water')
        >>> s1 = Stream(Ethanol=2)
        >>> s1.mu
        0.0010780252844121937
        """
        # Katti, P.K.; Chaudhri, M.M. (1964). "Viscosities of Binary Mixtures of Benzyl Acetate with Dioxane, Aniline, and m-Cresol". Journal of Chemical and Engineering Data. 9 (1964): 442–443.
        # molfrac = self.molfrac
        # props = self._species._props
        # mus = np.array(props('mu', self._mol, self.T, self.P, self._phase))
        # Vms = np.array(props('mu', self._mol, self.T, self.P, self._phase))
        # pos = np.where(molfrac != 0)
        # return np.exp((molfrac[pos]*np.log(mus[pos]*Vms[pos])).sum())/self.Vm
        return self._species._prop('mu', self._mol, self.T, self.P, self._phase)

    @property
    def k(self):
        """Thermal conductivity as a function of T and P (W/m/k).

        >>> from biosteam import *
        >>> Stream.species = Species('Ethanol', 'Water')
        >>> s1 = Stream(Ethanol=2)
        >>> s1.k
        0.16476716002011285
        """
        return self._species._prop('k', self._mol, self.T, self.P, self._phase)

    @property
    def alpha(self):
        """Thermal diffusivity as a function of T and P (m^2/s).

        >>> from biosteam import *
        >>> Stream.species = Species('Ethanol', 'Water')
        >>> s1 = Stream(Ethanol=2)
        >>> s1.alpha
        8.614274122585474e-08
        """
        return self._species._prop('alpha', self._mol, self.T, self.P, self._phase)

    @property
    def sigma(self):
        """Surface tension as a function of T (N/m).

        >>> from biosteam import *
        >>> Stream.species = Species('Ethanol', 'Water')
        >>> s1 = Stream(Ethanol=2)
        >>> s1.sigma
        0.02188440412824106
        """
        return self._species._prop('sigma', self._mol, self.T, self.P, self._phase)

    @property
    def Pr(self):
        """Prandtl number as a function of T and P (non-dimensional).

        >>> from biosteam import *
        >>> Stream.species = Species('Ethanol', 'Water')
        >>> s1 = Stream(Ethanol=2)
        >>> s1.Pr
        15.923321696004477
        """
        return self._species._prop('Pr', self._mol, self.T, self.P, self._phase)

    @property
    def P_vapor(self):
        """Vapor pressure (Pa), not taking into account light species always in gas phase.
        
        >>> from biosteam import *
        >>> Stream.species = Species('Ethanol', 'Water')
        >>> s1 = Stream(Ethanol=2)
        >>> s1.P_vapor
        array([7872.1566667784855, 0 ])
        """
        mol = self.mol
        species = self._species
        indices = species._equilibrium_indices(mol>0)
        compounds = species._compounds
        N = len(indices)
        P_vapor = np.zeros_like(mol)
        if N==0: return P_vapor
        species = [compounds[i] for i in indices]
        mol = self.mol[indices]
        x = mol/mol.sum()
        T = self.T
        Psat = [s.VaporPressure(T) for s in species]
        self._gamma.species = species
        P_vapor[indices] = x * Psat * self._gamma(x, T)
        return P_vapor
        
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

        >>> from biosteam import *
        >>> Stream.species = Species('Ethanol', 'Water')
        >>> s1 = Stream(Water=2)
        >>> s1.nonzero_species
        [1], ['Water']
        """
        return nonzero_species(self._species, self.mol)
    
    def quantity(self, prop_ID):
        """Return a property as a Quantity object as described in the `pint package <https://pint.readthedocs.io/en/latest/>`__ 

        **Parameters**

             **prop_ID:** [str] name of the property (e.g. 'mol', 'H', 'k' ...)

        Example:

        >>> from biosteam import *
        >>> Stream.species = Species('Ethanol', 'Water')
        >>> s1 = Stream(Ethanol=2)
        >>> s1.quantity('Cp')
        2.4337467143619698 J/K/g

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

        **Example**

        >>> from biosteam import *
        >>> Stream.species = Species('Ethanol', 'Water')
        >>> s1 = Stream(Ethanol=2)
        >>> s1.derivative('rho', 'T')
        -0.8917376157800969

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
            out = Stream(ID, flow=s._mol, species=s._species,
                         phase=s._phase, T=s.T, P=s.P)
        return out
    
    def copylike(self, stream):
        """Copy flow rates, T, P, and phase of stream to self.

        >>> from biosteam import *
        >>> Stream.species = Species('Ethanol', 'Water')
        >>> s1 = Stream()
        >>> s2 = Stream(Ethanol=1, Water=2, T=400, phase='g')
        >>> s1.copylike(s2)
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
            
    def copyflow(self, stream, species=None, *, remove=False, exclude=False):
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
        if species is None:
            self._mol[:] = stream.mol
            if remove: stream._mol[:] = 0
        else:
            indices = self.indices(species)
            if exclude:
                self._mol[:] = stream.mol
                self._mol[indices] = 0
                if remove:
                    mol = stream._mol
                    if isinstance(stream, MS.MixedStream):
                        mol[:], mol[:, indices] = 0, mol[:, indices]
                    else:
                        mol[:], mol[indices] = 0, mol[indices]
            else:
                self._mol[:] = 0
                self._mol[indices] = stream.mol[indices]
                if remove: 
                    if isinstance(stream, MS.MixedStream):
                        stream._mol[phase_index[self.phase], indices] = 0
                    else:
                        stream._mol[indices] = 0

    def recieve_vent(self, stream, efficiency=1):
        assert self._species is stream._species, 'species must be the same to recieve vent'
        y = stream.P_vapor/self.P
        Y = y.sum()
        vent_mol = self._mol
        mol2vent = vent_mol.sum()*y*Y/(1-Y)*efficiency
        stream._mol[:] -= mol2vent
        vent_mol += mol2vent

    def empty(self):
        """Set flow rates to zero

        >>> from biosteam import *
        >>> Stream.species = Species('Ethanol', 'Water')
        >>> s1 = Stream(Water=1)
        >>> s1.empty()
        >>> s1.mol
        array([0, 0])
        """
        self._mol[:] = 0
    
    def bubble_T(self):
        """Bubble point at current composition and pressure.

        **Returns**

            **T:** [float] Bubble point temperature (T).

            **y:** [numpy ndarray] Vapor phase composition.
            
            **indices:** [list] Indices of species in equilibrium

        >>> from biosteam import *
        >>> stream = Stream(flow=(0.6, 0.4),
        ...                 species=Species('Ethanol', 'Water'))
        >>> stream.bubble_T()
        (352.2820850833474, array([0.703, 0.297]), [0, 1])

        """
        mol = self.mol
        # If just one specie in equilibrium, return Tsat
        indices = self._species._equilibrium_indices(mol>0)
        cmps = self._species._compounds
        N = len(indices)
        if N == 1:
            return (cmps[indices[0]].Tsat(self.P), np.array((1,)), indices)
        elif N == 0:
            raise EquilibriumError('no species available for phase equilibrium')
        mol = mol[indices]
        self._gamma.species = [cmps[i] for i in indices]
        # Solve and return bubble point
        return (*self._bubble_point.solve_Ty(mol/mol.sum(), self.P), indices)
    
    def bubble_P(self):
        """Bubble point at current composition and temperature.

        **Returns**

            **P:** [float] Bubble point pressure (Pa).

            **y:** [numpy ndarray] Vapor phase composition.
            
            **indices:** [list] Indices of species in equilibrium

        >>> from biosteam import *
        >>> s1 = Stream(flow=(0.703, 0.297), T=352.28,
                        species=Species('Ethanol', 'Water'))
        >>> s1.bubble_P()
        (103494.17209657285, array([0.757, 0.243]), [0, 1])

        """
        mol = self.mol
        # If just one specie in equilibrium, return Tsat
        indices = self._species._equilibrium_indices(mol>0)
        cmps = self._species._compounds
        N = len(indices)
        if N == 1:
            return (cmps[indices[0]].VaporPressure(self.T), np.array((1,)), indices)
        elif N == 0:
            raise EquilibriumError('no species available for phase equilibrium')
        mol = mol[indices]
        self._gamma.species = [cmps[i] for i in indices]
        # Solve and return bubble point
        return (*self._bubble_point.solve_Py(mol/mol.sum(), self.T), indices)

    def dew_T(self):
        """Dew point at current composition and pressure.
        
        **Returns**

            **T:** [float] Dew point temperature (K).

            **x:** [numpy array] Liquid phase composition.
            
            **indices:** [list] Indices of species in equilibrium

        >>> stream = Stream(flow=(0.5, 0.5),
        ...                 species=Species('Ethanol', 'Water'))
        >>> stream.dew_T()
        (357.45184742263075, array([0.151, 0.849]), [0, 1])
        
        """
        mol = self.mol
        # If just one specie in equilibrium, return Tsat
        indices = self._species._equilibrium_indices(mol>0)
        cmps = self._species._compounds
        N = len(indices)
        if N == 1:
            return (cmps[0].Tsat(self.P), np.array((1,)), indices)
        elif N == 0:
            raise EquilibriumError('no species available for phase equilibrium')
        mol = mol[indices]
        self._gamma.species = [cmps[i] for i in indices]
        # Solve and return dew point
        return (*self._dew_point.solve_Tx(mol/mol.sum(), self.P), indices)
    
    def dew_P(self):
        """Dew point at current composition and temperature.
        
        **Returns**

            **P:** [float] Dew point pressure (Pa).

            **x:** [numpy array] Liquid phase composition.
            
            **indices:** [list] Indices of species in equilibrium

        >>> from biosteam import *
        >>> stream = Stream(flow=(0.703, 0.297), T=352.28,
        ...                 species=Species('Ethanol', 'Water'))
        >>> stream.dew_P()
        (101328.47030327446, array([0.6, 0.4]), [0, 1])
        
        """
        mol = self.mol
        # If just one specie in equilibrium, return Tsat
        indices = self._species._equilibrium_indices(mol>0)
        cmps = self._species._compounds
        N = len(indices)
        if N == 1:
            return (cmps[0].VaporPressure(self.T), np.array((1,)), indices)
        elif N == 0:
            raise EquilibriumError('no species available for phase equilibrium')        
        mol = mol[indices]    
        self._gamma.species = [cmps[i] for i in indices]
        # Solve and return dew point
        return (*self._dew_point.solve_Px(mol/mol.sum(), self.T), indices)

    # Dimensionless number methods
    def Re(self, L, A=None):
        """Return Reynolds number.

        **Parameters**

             **L:** [float] Characteristic lenght (m)

             **A:** [float] Cross-sectional area (m^2). If A is None, assume flow through a cylindrical pipe.

        Example

        >>> from biosteam import *
        >>> Stream.species = Species('Ethanol', 'Water')
        >>> s1 = Stream(Water=200)
        >>> s1.Re(0.2, 0.031416)
        27639.956372923512
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
    def T_equilibrium(streams, T_guess=None, Q_in=0, approximate=True):
        """Bring all streams to temperature equilibrium.

        **Parameters**

            **streams:** [Stream] All streams
        
            **T_guess:** [float] Equilibrium temperature guess (K)

            **Q_in:** [float] Heat addition to streams (kJ/hr)     
            
            **approximate:** [bool] If True, approximate energy balance with constant heat capacity

        **Return**

             **T:** [float] New temperature (K)

        **Example**

        >>> from biosteam import *
        >>> Stream.species = Species('Water', 'Ethanol')
        >>> s1 = Stream(Water=2, T=300)
        >>> s2 = Stream(Ethanol=1, Water=2, T=340)
        >>> Stream.T_equilibrium([s1, s2])
        >>> (s1.T, s2.T)
        (325.8281231556553, 325.8281231556553)
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

        .. Note:: Ignores contents of `s_sum` as if it was empty

        **Example**

        >>> from biosteam import *
        >>> Stream.species = Species('Water', 'Ethanol')
        >>> s1 = Stream(Water=2, T=300)
        >>> s2 = Stream(Ethanol=1, Water=2, T=340)
        >>> s_sum = Stream('s_sum')
        >>> Stream.sum(s_sum, [s1, s2])
        >>> s_sum.show()
        Stream: s_sum
         phase: 'l', T: 326.29 K, P: 101325 Pa
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
                        s_sum._mol[:] += s._mol
                    else:
                        s_sum._mol[phase_index[s._phase]] += s._mol
            else:
                # For Stream objects
                for s in other_streams:
                    s_sum._mol[:] += s.mol
        except Exception as Error:
            if not inst(s_sum, Stream):
                raise TypeError('s_sum must be a Stream object, not '
                                "'{type(s_sum).__name__}'")
            for s in streams:
                if not inst(s, Stream):
                    raise TypeError('streams must only contain Stream '
                                    f'objects, not {type(s).__name__} objects')
            raise Error

        # Energy Balance
        if energy_balance: s_sum.H = H

    # MixedStream compatibility
    def enable_phases(self):
        """Cast stream into a MixedStream object."""
        mol = self._mol
        self.__class__ = MS.MixedStream
        self._setflows(np.zeros((4, self._species._N)))
        self._mol[phase_index[self._phase]] = mol
        self._lL_split_cached = (None,)
        self._VLE = VLE(self)

    def disable_phases(self, phase):
        """Cast stream into a Stream object.
        
        **Parameters**
        
            **phase:** {'s', 'l', 'g'} desired phase of stream
            
        """
        self._phase = phase
            
    @property
    def VLE(self):
        """A callable VLE object for vapor-liquid equilibrium.

        **Parameters**
        
            **Specify two:**
                * **P:** [float] Operating pressure (Pa)
                * **Q:** [float] Energy input (kJ/hr)
                * **T:** [float] Operating temperature (K)
                * **V:** [float] Molar vapor fraction
                * **x:** [numpy array] Molar composition of liquid (for binary mixture)
                * **y:** [numpy array] Molar composition of vapor (for binary mixture)

            **Optional:**
            
                **species_IDs:** [tuple] IDs of species in equilibrium.
                     
                **LNK:** tuple[str] Light non-keys that remain as a vapor.
        
                **HNK:** tuple[str] Heavy non-keys that remain as a liquid.

        .. Note:
           LNK and HNK are not taken into account for equilibrium. Parameters not specified are None by default.

        """
        self.enable_phases()
        return self._VLE
        
    def LLE(self, species_IDs=(), split=None, lNK=(), LNK=(),
            solvents=(), solvent_split=(),
            P=None, T=None, Q=None):
        self.enable_phases()
        self.LLE(species_IDs, split, lNK, LNK,
                 solvents, solvent_split, P, T, Q)

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
        if hasattr(self, '_mol'):
            nonzero, species = self.nonzero_species
        else:
            return basic_info + f' link: {self._link}'
        T_units, P_units, flow_units, fraction = [(i if i is not None else j) for i, j in
                                                  zip((T, P, flow, fraction), self.display_units)]
        basic_info += self._info_phaseTP(self._phase, T_units, P_units)
        len_ = len(nonzero)
        if len_ == 0:
            return basic_info + ' flow: 0' 
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
        if 'frac' in flow:
            flow = getattr(self, flow)[nonzero]
        else:
            flow = _Q(getattr(self, flow)[nonzero], units[flow]).to(flow_units).magnitude
        
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
              + end.replace('*', new_line_spaces + 'net' + (maxlen-4)*' ' + '  '))

    def show(self, T=None, P=None, flow=None, fraction=None):
        """Print all specifications."""
        print(self._info(T, P, flow, fraction))
    _ipython_display_ = show

    def _disconnect(self):
        outs = self._source and self._source._outs
        if outs: outs[outs.index(self)] = MissingStream
        ins = self._sink and self._sink._ins
        if ins: ins[ins.index(self)] = MissingStream
        self._source = self._sink = None

    def __str__(self):
        if self.ID: return self.ID
        else: return type(self).__name__

    def __repr__(self):
        if self.ID: return f'<{type(self).__name__}: {self.ID}>'
        else: return f'<{type(self).__name__}>'
        
        
from . import _mixed_stream as MS