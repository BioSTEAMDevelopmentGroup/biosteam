# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 14:26:41 2018

@author: yoelr
"""
from biosteam import Unit, Stream
from biosteam.meta_classes import metaFinal

class Mixer(Unit, metaclass=metaFinal):
    """Create a mixer that mixes any number of streams together.
    
    **Parameters**
    
    **ins**
    
        [:] Input streams
        
    **outs**
    
        [0] Mixed stream
    
    """
    _N_ins = 2
    _N_outs = 1

    def __init__(self, ID='', outs_ID=(), ins_ID=None):
        """Initialize Mixer object. See help(type(self)) for accurate signature.

        **Parameters**

            ID: [str] Unique identification. If set as '', a default ID will be chosen.

            description: [str] User description of unit

            outs_ID: tuple[str] IDs to initialize output streams. If None, leave streams missing. If empty, default IDs will be given.
        
            ins_ID: tuple[str] IDs to initialize input streams. If None, leave streams missing. If empty, default IDs will be given.
            
        """
        self.ID = ID
        self._init_ins(ins_ID)
        self._init_outs(outs_ID)

    def _run(self):        Stream.sum(self.outs[0], self.ins, self._in_loop)
