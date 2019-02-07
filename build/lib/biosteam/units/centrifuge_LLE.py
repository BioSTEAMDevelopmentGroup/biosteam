# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 21:23:56 2018

@author: yoelr
"""
from biosteam import Unit, MixedStream
from biosteam.exceptions import DesignError
from biosteam.units.flash import Flash, EstimateFlash
import copy


class Centrifuge_LLE(Unit):
    r"""Create an equlibrium based centrifuge with the option of having liquid non-keys and LIQUID non-keys completly separate into their respective phases.

    Equation for f.o.b cost:

    :math:`C_{fob}^{2007} = 28100 Q^{0.574} (0.1 < Q < 100 \frac{m^3}{h})` 

    **Parameters**

        **specie_IDs:** *tuple[str]* IDs of equilibrium species
    
        **split:** *tuple[float]* Initial guess split fractions of equilibrium specise to the 'liquid' phase.
    
        **lNK:** *tuple[str]* Species assumed to completely remain in the 'liquid' phase.
    
        **LNK:** *tuple[str]* Species assumed to completely remain in the 'LIQUID' phase.
    
        **solvents:** *tuple[str]* Species corresponding to specified solvent_split.
    
        **solvent_split:** *tuple[float]* Split fractions of each specie to the 'liquid' phase.
         
    **ins**
    
        [0] Input stream
        
    **outs**
    
        [0] 'liquid' phase stream
        
        [1] 'LIQUID' phase stream
    
    **Examples**
    
        :doc:`Centrifuge_LLE Example`
    
    """
    kwargs = {'specie_IDs': None,
              'split': None,
              'lNK': (),
              'LNK': (),
              'solvents': (),
              'solvent_split': ()}

    bounds = {'Flow rate': (0.1, 100)}

    def _setup(self):
        liq, LIQ = self.outs
        self._cached = cached = {}
        liq.phase = 'l'
        LIQ.phase = 'L'
        cached['mixed stream'] =  MixedStream()
        if self.kwargs['specie_IDs'] is None:
            self.kwargs['specie_IDs'] = liq.specie_IDs

    def _run(self):
        liq, LIQ = self.outs
        feed = self.ins[0]

        kwargs = self.kwargs
        LLE_kwargs = copy.copy(kwargs)

        ms = self._cached['mixed stream']
        ms.empty()
        ms.liquid_mol = feed.mol
        ms.LLE(T=feed.T,  **LLE_kwargs)
        liq.mol = ms.liquid_mol
        LIQ.mol = ms.LIQUID_mol
        liq.T = LIQ.T = ms.T
        liq.P = LIQ.P = ms.P

    def _design(self):
        """
        * 'Flow rate': (m^3/hr)
        """
        Design = self.results['Design']
        Design['Flow rate'] = self._volnet_out
        return Design

    def _cost(self):
        """
        * 'Vessel cost': (USD)
        """
        r = self.results
        Design = r['Design']
        Cost = r['Cost']
        Q = Design['Flow rate']
        Cost['Vessel cost'] = 28100 * Q ** 0.574 * self.CEPCI/525.4
        return Cost

class Centrifuge_LLE_Lazy(Centrifuge_LLE):
    _N_heat_util = 0
    kwargs = {'Kspecies': [],  # list of species that correspond to Ks
              'Ks': [],  # list of molar ratio partition coefficinets,
              # Ks = y/x, where y and x are molar ratios of two different phases
              'top_solvents': [],  # list of species that correspond to top_split
              'top_split': [],  # list of splits for top_solvents
              'bot_solvents': [],  # list of species that correspond to bot_split
              'bot_split': []}  # list of splits for bot_solvents

    def _setup(self):
        pass

    _run = EstimateFlash._run
    _simple_run = EstimateFlash._simple_run
    
