# -*- coding: utf-8 -*-
"""
Created on Sat Aug 18 13:42:33 2018

@author: yoelr
"""
from . import compounds as comps

Chemical = comps.Chemical
Compound = comps.Compound
__all__ = ('Species',)

# TODO: Fix how Dortmund groups are selected in thermo.Chemical. Glycerol defaults to wrong value.

# %% Managing ChEDL thermodynamic properties

class Species:
    """Create Species object that contains Compound objects as attributes.

    **Parameters**

        ***IDs:** [str] Strings should be one of the following [-]:
                  * Name, in IUPAC form or common form or a synonym registered in PubChem
                  * InChI name, prefixed by 'InChI=1S/' or 'InChI=1/'
                  * InChI key, prefixed by 'InChIKey='
                  * PubChem CID, prefixed by 'PubChem='
                  * SMILES (prefix with 'SMILES=' to ensure smiles parsing)
                  * CAS number
         
        **Optional**
        
        **cls:** [type Compound] The type of objects that will be created and stored as attributes. Defaults to biosteam.Chemical
        
    """
    units = Compound.units
    _immutable = set()
    
    @classmethod
    def tospecies(cls, compounds, copy=True):
        """Return a Species object from an iterable of Compound objects."""
        if not copy and isinstance(compounds, Species): return compounds
        self = cls.__new__(cls)
        setter = super(cls, self).__setattr__
        for c in compounds: 
            if isinstance(c, Compound): setter(c.ID, c)
        return self
    
    def __init__(self, *ID, cls=None):
        # Sphinx docs look ugly if a class is in the default argument
        if cls is None: cls = Chemical
        
        # Make sure ID is a tuple
        if not ID:
            return
        elif isinstance(ID, str):
            ID = (ID,)
        elif not isinstance(ID[0], str):
            ID = ID[0] # For backwards compatibility
        
        # Set Compound object attributes
        attrset = super().__setattr__
        for n in ID:
            try:
                compound = cls(n)
            except:
                raise Exception(f"Compound, '{n}', not defined in databank.")
            attrset(n, compound)

    def __setattr__(self, ID, compound):
        if self in self._immutable:
            raise ValueError('Cannot alter species object attached to Stream objects.')
        if isinstance(compound, Compound):
            compound.ID = ID
            super().__setattr__(ID, compound)
        elif ID == 'ID':
            new_IDs = compound
            new_species = [getattr(self, i) for i in new_IDs]
            self.__dict__.clear()
            for ID, compound in zip(new_IDs, new_species):
                super().__setattr__(ID, compound)
        else:
            raise TypeError('Can only set Compound objects as attributes')
    
    def __delattr__(self, ID):
        if self in self._immutable:
            raise ValueError('Cannot alter species object attached to Stream objects.')
        else:
            super().__delattr__(ID)
    
    def __iter__(self):
        return iter(self.__dict__.values())
    
    @property
    def CAS(self):
        return (s.CAS for s in self)
    
    @property
    def ID(self):
        return iter(self.__dict__.keys())

    def getprops(self, species_IDs, prop_ID, T=None, P=None, phase=None) -> list:
        """Return list of the desired property, prop_ID, for each compound in species_ID."""
        if isinstance(species_IDs, str):
            species = (getattr(self, species_IDs),)
        else:
            species = (getattr(self, i) for i in species_IDs)
        return self._getprops(species, prop_ID, T, P, phase)

    @staticmethod
    def _getprops(species, prop_ID, T, P, phase) -> list:
        """Return list of the desired property, prop_ID, for given species."""
        props = []; sat = setattr; gat = getattr
        for s in species:
            if T: sat(s, 'T', T)
            if P: sat(s, 'P', P)
            if phase: sat(s, 'phase', phase)
            props.append(gat(s, prop_ID))
        return props

    def setprops(self, species_IDs, prop_ID, new_values):
        """Set new values to specie property, prop_ID, for given species_IDs."""
        sat = setattr; gat = getattr
        for ID in species_IDs: 
            sat(gat(self, ID), prop_ID, new_values)

    # Representation
    def _info(self):
        """Information on self"""
        return 'Species: \n ' + str(list(self.ID)).replace("'", '')[1:-1]

    def show(self):
        """print information on self"""
        print(self._info())

    def __repr__(self):
        return '<Species: ' + str(list(self.ID)).replace("'", '')[1:-1] + '>'

    
