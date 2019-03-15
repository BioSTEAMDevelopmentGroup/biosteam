# -*- coding: utf-8 -*-
"""
Created on Sat Aug 18 13:42:33 2018

@author: yoelr
"""
from .chemical import Chemical
from .compound import Compound
from . import np, units_of_measure
from .exceptions import IDconflict
from .utils import missing_method


__all__ = ('Species',)

# TODO: Fix how Dortmund groups are selected in thermo.Chemical. Glycerol defaults to wrong value.

# %% Managing ChEDL thermodynamic properties
inf = np.inf
TPphase = ('T', 'P', 'phase')

class Species(list):
    """Create Species object that contains biosteam 'Chemical' objects as attributes (by ID), or 'Specie' objects if an ID is invalid.

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
    units = units_of_measure

    @classmethod
    def as_species(cls, compounds, copy=True):
        """Return a Species object from an iterable of Compound objects."""
        if not copy and isinstance(compounds, Species): return compounds
        self = cls.__new__(cls)
        for s in compounds: setattr(self, s.ID, s) 
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
        for n in ID:
            try:
                compound = cls(n)
            except:
                raise Exception(f"Compound, '{n}', not defined in databank.")
            setattr(self, n, compound)

    def __setattr__(self, ID, compound):
        spdict = self.__dict__
        if isinstance(compound, Compound):
            compound.ID = ID
            if ID in spdict: super().remove(compound)
            spdict[ID] = compound
            if compound not in self: super().append(compound)
        elif ID == 'ID':
            new_IDs = compound
            old_index = [i.ID for i in self].index
            species = tuple(self)
            new_species = [species[old_index(i)] for i in new_IDs]
            for ID, compound in zip(new_IDs, new_species):
                super().__setattr__(ID, compound)
            super().__init__(new_species)
            
        else:
            raise TypeError('Can only set Compound objects as attributes')
    
    append = missing_method
    insert = missing_method
    extend = missing_method
    __setitem__ = missing_method
    
    def clear(self):
        super().clear()
        self.__dict__.clear()
    
    def sort(self, key=lambda i: i.ID, reverse=False):
        super().sort(key, reverse)
    
    def copy(self):
        copy = super().copy()
        copy.__dict__.update(self)
    
    def pop(self, index):
        compound = super().pop(index)
        del self.__dict__[compound.ID]
        
    def remove(self, compound):
        super().remove(compound)
        super().__delattr__(compound.ID)
    
    def __delattr__(self, ID):
        if hasattr(self, ID):
            compound = getattr(ID)
            super().remove(compound)
            super().__delattr__(ID)
    
    def __delitem__(self, index):
        compound = self.index(index)
        super().__delitem__(compound)
        super().__delattr__(compound.ID)
        
    @property
    def CAS(self):
        return (s.CAS for s in self)
    
    @property
    def ID(self):
        return (s.ID for s in self)

    def get_props(self, species_IDs, prop_ID, T=None, P=None, phase=None) -> list:
        """Return list of the desired property, prop_ID, for each compound in species_ID."""
        if isinstance(species_IDs, str):
            species = (getattr(self, species_IDs),)
        else:
            species = (getattr(self, i) for i in species_IDs)
        return self._get_props(species, prop_ID, T, P, phase)

    @staticmethod
    def _get_props(species, prop_ID, T, P, phase) -> list:
        """Return list of the desired property, prop_ID, for given species."""
        props = []
        for s in species:
            for attr, val in zip(TPphase, (T, P, phase)):
                if val: setattr(s, attr, val)
            props.append(getattr(s, prop_ID))
        return props

    def set_props(self, species_IDs, prop_ID, new_values):
        """Set new values to specie property, prop_ID, for given species_IDs."""
        for ID in species_IDs:
            specie = getattr(self, ID)
            setattr(specie, prop_ID, new_values)

    # Representation
    def _info(self):
        """Information on self"""
        return (f'Species:\n {self.ID}').replace("'", '')

    def show(self):
        """print information on self"""
        print(self._info())

    def __repr__(self):
        return (f'<Species: {list(self.ID)}>').replace("'", '')
