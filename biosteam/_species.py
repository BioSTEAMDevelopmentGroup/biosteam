# -*- coding: utf-8 -*-
"""
Created on Sat Aug 18 13:42:33 2018

@author: yoelr
"""
from . import compounds as cmps
from ._exceptions import UndefinedCompound
import numpy as np

inf = np.inf
__all__ = ('Species', 'WorkingSpecies')


# TODO: Fix how Dortmund groups are selected in thermo.Chemical. Glycerol defaults to wrong value.

# %% Managing ChEDL thermodynamic properties
_setattr = object.__setattr__

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
    @classmethod
    def tospecies(cls, compounds):
        """Return a Species object from an iterable of Compound objects."""
        self = cls.__new__(cls)
        for c in compounds: setattr(self, c.ID, c)
        return self
    
    def __init__(self, *IDs, cls=None):
        # Sphinx docs look ugly if a class is in the default argument
        if not IDs: return
        if cls is None: cls = cmps.Chemical
        
        # Set Compound object attributes
        for n in IDs:
            try:
                compound = cls(n)
            except:
                raise LookupError(f"'{n}', not defined in data bank")
            _setattr(self, n, compound)
    
    @property
    def IDs(self):
        return tuple(self.__dict__)
    
    def __len__(self):
        return len(self.__dict__)
    
    def __contains__(self, compound):
        return compound in self.__dict__.values()
    
    def __iter__(self):
        yield from self.__dict__.values()
    
    def __setattr__(self, ID, compound):
        if isinstance(compound, cmps.Compound):
            if compound in self:
                raise AttributeError(f'compound {repr(compound)} already defined')
            compound.ID = ID
            _setattr(self, ID, compound)
        elif ID == 'IDs':
            _setattr(self, '__dict__', {i: getattr(self, i) for i in compound})
        else:
            raise TypeError('can only set Compound objects as attributes')
      
    def extend(self, compounds):
        """Extend species object with more compounds."""
        for c in compounds: 
            if isinstance(c, cmps.Compound): _setattr(self, c.ID, c)
            else: raise ValueError('compounds must be an iterable of Compound objects')
            
    def getprops(self, IDs, prop_ID, T, P, phase):
        """Return list of the desired property, prop_ID, for each compound in IDz."""
        props = []; getattr_ = getattr
        for ID in IDs:
            s = getattr_(self, ID)
            s.T = T
            s.P = P
            s.phase = phase
            props.append(getattr_(s, prop_ID))
        return props
    
    def setprops(self, IDs, prop_ID, new_values):
        """Set new values to compound property, prop_ID, for given IDs."""
        getattr_ = getattr
        for ID in IDs: 
            _setattr(getattr_(self, ID), prop_ID, new_values)
    
    def __repr__(self):
        return f'<{type(self).__name__}: {", ".join(self.IDs)}>'


class WorkingSpecies:
    """Cast a Species object as a WorkingSpecies object that represents the working species of a process."""
    
    def __new__(cls, self):
        if self.__class__ is cls: return self
        _setattr(self, '__class__', cls)
        me = self.__dict__
        tup = tuple
        compounds = tup(me.values())
        IDs = tup(i for i in compounds)
        CAS = tup(i.CAS for i in compounds)
        Nspecies = len(IDs)
        index = tup(range(Nspecies))
        for i in compounds: me[i.CAS] = i
        me['_indexdct'] = dict((*zip(CAS, index), *zip(IDs, index)))
        me['_immutable'] = True
        me['_index'] = index
        me['_compounds'] = compounds
        me['_IDs'] = IDs
        me['_Nspecies'] = Nspecies
        me['_MW'] = np.array([i.MW for i in compounds])
        return self
    
    def __dir__(self):
        return [i for i in super().__dir__() if i.isalnum()]

    def set_synonym(self, ID, synonym):
        compound = getattr(self, ID)
        dct = self.__dict__
        if synonym in dct and dct[synonym] is not compound:
            raise ValueError(f"synonym '{synonym}' already in use by {repr(dct[synonym])}")
        else:
            self._indexdct[synonym] = self._indexdct[ID]
            dct[synonym] = compound
            
    def kwarray(self, **data):
        """Return an array with entries that correspond to compound IDs.
        
        **Parameters**
        
            ****data:** ID-value pair.
            
        **Examples**
        
            >>> from biosteam import *
            >>> Stream.species = Species('Water', 'Ethanol')
            >>> Stream.species.kwarray(Water=2)
            array([2., 0.])
        
        """
        return self.array(data, [*data.values()])
    
    def array(self, IDs, data):
        """Return an array with entries that correspond to species IDs.
        
        **Parameters**
        
            **IDs:** [iterable] compound IDs.
            
            **data:** [array_like] Data corresponding to IDs.
            
        **Examples**
        
            >>> from biosteam import *
            >>> Stream.species = Species('Water', 'Ethanol')
            >>> Stream.species.array(['Water'], [2])
            array([2., 0.])
        
        """
        array = np.zeros(self._Nspecies)
        array[self.indices(IDs)] = data
        return array

    def indices(self, IDs):
        """Return flow indices of specified species.

        **Parameters**

             **IDs:** [iterable] species IDs or CAS numbers.

        Example

        Indices by ID:
        
        >>> from biosteam import *
        >>> Stream.species = Species(['Ethanol', 'Water'])
        >>> Stream.species.indices(['Water', 'Ethanol'])
        [1, 0]

        Indices by CAS number:
        
        >>> Stream.species.indices(['7732-18-5', '64-17-5']):
        [1, 0]

        """
        try:
            dct = self._indexdct
            return [dct[i] for i in IDs]
        except KeyError:
            for i in IDs:
                if i not in dct: raise UndefinedCompound(i)       
    
    def __delattr__(self, ID):
        AttributeError("'{type(self).__name__} object is read-only")
    
    def __setattr__(self, ID, name):
        AttributeError("'{type(self).__name__} object is read-only")
    
    def __len__(self):
        return self._Nspecies
    
    def __contains__(self, compound):
        return compound in self._compounds
    
    def __iter__(self):
        yield from self._compounds
    
    @property
    def IDs(self):
        """[tuple] IDs of Species object."""
        return self._IDs

    ### Material Properties ###
    getprops = Species.getprops
    setprops = Species.setprops

    # General methods for getting ideal mixture properties
    def _props(self, ID, mol, T, P, phase):
        """Return component property list."""
        out = np.zeros(self._Nspecies)
        getattr_ = getattr
        cmps = self._compounds
        for i in self._index:
            if mol[i]: 
                c = cmps[i]
                c.P = P
                c.T = T
                c.phase = phase
                out[i] = getattr_(c, ID)
        return out

    def _propflows(self, ID, mol, T, P, phase):
        """Return array of component properties * kmol/hr."""
        return self._props(ID, mol, T, P, phase) * mol

    def _propflow(self, ID, mol, T, P, phase):
        """Return flow of component properties."""
        return (self._props(ID, mol, T, P, phase)*mol).sum()

    def _prop(self, ID, mol, T, P, phase):
        """Return molar weighted average property."""
        try:
            return (self._props(ID, mol, T, P, phase)*mol).sum()/mol.sum()
        except:
            # This usually fails due to a missing property.
            # Avoid this problem by ignoring the property.
            # TODO: Maybe include a warning.
            propflow = 0
            molnet = 0
            getattr_ = getattr
            cmps = self._compounds
            for i in self._index:
                if mol[i]:
                    c = cmps[i]
                    c.P = P
                    c.T = T
                    c.phase = phase
                    prop = getattr_(c, ID)
                    if prop:
                        molnet += mol[i]
                        propflow += prop*mol[i]
            return propflow/molnet

    def _mixedpropflows(self, ID, molarray, T, P):
        """Return array of component properties * kmol/hr."""
        pf = self._propflows
        return (pf(ID, molarray[0], T, P, 's')
              + pf(ID, molarray[1], T, P, 'l')
              + pf(ID, molarray[2], T, P, 'l')
              + pf(ID, molarray[3], T, P, 'g'))
    
    def _mixedpropflow(self, ID, molarray, T, P):
        """Return flow of component properties."""
        return self._mixedpropflows(ID, molarray, T, P).sum()
    
    def _mixedprop(self, ID, molarray, T, P):
        """Return molar weighted average property."""
        return self._mixedpropflows(ID, molarray, T, P).sum()/molarray.sum()  
    
    ### Equilibrium ###
    
    def _equilibrium_species(self, mol):
        """Return species and indexes of species in equilibrium."""
        species = []; index = []; cmps = self._compounds
        for i in self._index:
            if mol[i]:
                c = cmps[i]
                if c.UNIFAC_Dortmund_groups and c.Tb not in (inf, None):
                    species.append(c); index.append(i)
        return species, index

    def _heavy_species(self, mol):
        """Return species and indexes of heavy species not in equilibrium."""
        species = []; index = []; cmps = self._compounds
        for i in self._index:
            if mol[i]:
                c = cmps[i]
                if c.Tb in (inf, None) or not c.UNIFAC_Dortmund_groups:
                    species.append(c); index.append(i)
        return species, index

    def _light_species(self, mol):
        """Return species and indexes of light species not in equilibrium."""
        species = []; index = []; cmps = self._compounds
        for i in self._index:
            if mol[i]:
                c = cmps[i]
                if c.Tb is -inf:
                    species.append(c); index.append(i)
        return species, index
    
    __repr__ = Species.__repr__

    
