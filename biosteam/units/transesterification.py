# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 22:49:58 2018

@author: yoelr
"""
from biosteam.units.reactor import Reactor
from biosteam.exceptions import DesignError
from biosteam import np


class Transesterification(Reactor):
    """Create a transesterification reactor that converts 'Lipid' and 'Methanol' to 'Biodiesel' and 'Glycerol'. Finds the amount of catalyist 'NaOCH3' required and consumes it to 'NaOH' and 'Methanol'.
    
    **Parameters**
    
        **eff:** Efficiency of conversion (on a 'Lipid' basis)
        
        **r:** Methanol to lipid molar ratio
        
        **catalyst_molfrac:** Molar fraction of catalyst in methanol feed 
        
        **T:** Operating temperature (K)
    
    **ins**
    
        [0] Lipid feed
        
        [1] Fresh methanol and catalyst
        
    **outs**
    
        [0] Product
    
    **Examples**
    
        :doc:`Transesterification Example`
    
    """
    _units = {'Volume': 'm^3'}
    _kwargs = {'efficiency': None,  # fraction of theoretical conversion
               'r': None,  # Methanol to lipid molar ratio
               'T': None,
               'catalyst_molfrac': None}  # operating temperature (K)
    
    _bounds = {'Volume': (0.1, 20)}
    
    _tau = 1
    _N_ins = 2
    _N_outs = 1
    _N_heat_utilities = 1
    _has_power_utility = True

    #: kW/m3
    electricity_rate = 1.5

    def _run(self):
        feed, fresh_Methanol = self.ins
        out = self.outs[0]
        eff, r, T, catalyst_molfrac = (self._kwargs[i] for i in (
            'efficiency', 'r', 'T', 'catalyst_molfrac'))
        sp_index = feed._IDs.index

        # Reactant positions
        lipid_pos = sp_index('Lipid')
        Methanol_pos = sp_index('Methanol')
        Glycerol_pos = sp_index('Glycerol')
        biodiesel_pos = sp_index('Biodiesel')
        NaOH_pos = sp_index('NaOH')
        NaOCH3_pos = sp_index('NaOCH3')

        lipid = feed.mol[lipid_pos]

        # Reaction conversions by mol (to add to the feed)
        lipid_ch = -eff*lipid  # - free_lipid_ch/3
        bd_ch = eff*lipid*3
        gly_ch = -lipid_ch
        dummy = r*lipid
        NaOCH3_in = dummy*catalyst_molfrac  # from methanol stream
        Methanol_in = dummy - NaOCH3_in  # from methanol stream
        Methanol_ch = dummy + lipid_ch
        NaOCH3_ch = 0
        NaOH_ch = NaOCH3_in

        change_mol = np.array((lipid_ch, Methanol_ch, gly_ch,
                               bd_ch, NaOCH3_ch, NaOH_ch))

        # Output stream
        fresh_Methanol.mol[ [Methanol_pos, NaOCH3_pos] ] = (Methanol_in, NaOCH3_in)
        index = [lipid_pos, Methanol_pos, Glycerol_pos,
                biodiesel_pos, NaOCH3_pos, NaOH_pos]
        out.phase = 'l'
        out.mol[index] = feed.mol[index] + change_mol
        out.T = T
        out.P = feed.P

    def _design(self):
        self._heat_utilities[0](self._Hnet, self.outs[0].T)
        Design = self._results['Design']
        Design['Volume'] = self._tau * self._volnet_out / 0.8
        return Design
        
    def _cost(self):
        results = self._results
        Design = results['Design']
        Cost = results['Cost']
        Volume = Design['Volume']
        Cost['Reactor'] = self.CEPCI/525.4 * 15000 * Volume ** 0.55
        self._power_utility(Volume*self.electricity_rate)
    
        
        
