# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 14:26:41 2018

@author: yoelr
"""
from .. import Unit, Stream
from .._meta_final import metaFinal

class Mixer(Unit, metaclass=metaFinal):
    """Create a mixer that mixes any number of streams together.
    
    **Parameters**
    
    **ins**
    
        [:] Input streams
        
    **outs**
    
        [0] Mixed stream
    
    **Examples**
    
        Create a Mixer object with an ID and any number of input streams:
            
        .. code-block:: python
        
           >>> from biosteam import Species, Stream, Mixer
           >>> Stream.species = Species('Water', 'Ethanol')
           >>> s1 = Stream('s1', Water=20, T=340)
           >>> s2 = Stream('s2', Ethanol=10, T=300)
           >>> s3 = Stream('s3', Water=3, Ethanol=4)
           >>> M1 = Mixer('M1', ins=(s1, s2, s3), outs='s4')
           >>> M1.simulate()
           >>> M1.show()
           Mixer: M1
           ins...
           [0] s1
               phase: 'l', T: 340 K, P: 101325 Pa
               flow (kmol/hr): Water  20
           [1] s2
               phase: 'l', T: 300 K, P: 101325 Pa
               flow (kmol/hr): Ethanol  10
           [2] s3
               phase: 'l', T: 298.15 K, P: 101325 Pa
               flow (kmol/hr): Water    3
                               Ethanol  4
           outs...
           [0] s4
               phase: 'l', T: 317.54 K, P: 101325 Pa
               flow (kmol/hr): Water    23
               Ethanol  14
           
    """
    results = None
    _N_ins = 2
    _N_outs = 1
    _has_cost = False
    
    def __init__(self, ID='', outs=(), ins=None):
        self.ID = ID
        self._init_ins(ins)
        self._init_outs(outs)

    def _run(self):
        Stream.sum(self.outs[0], self.ins)
        
    simulate = _run
    summary = Unit._cost
