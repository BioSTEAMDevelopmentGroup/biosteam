# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena and contributors <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""

from .._unit import Unit

__all__ = ('Duplicator',)


class Duplicator(Unit, isabstract=True):
    """
    Create a Duplicator object that takes in one inlet stream
    and duplicates it to all outlet streams.

    Parameters
    ----------
    ins : stream
        Inlet stream.
    outs : stream sequence
        Duplicated outlet streams.

    Examples
    --------
    Create a Duplicator object with an ID and any number of outlet streams:

    >>> from biosteam import settings, Stream, units
    >>> settings.set_thermo(['Water', 'Ethanol'], cache=True)
    >>> feed = Stream('feed', Water=20, Ethanol=10, T=340)
    >>> D1 = units.Duplicator('D1', ins=feed, outs=('out_a', 'out_b', 'out_c'))
    >>> D1.simulate()
    >>> D1.show()
    Duplicator: D1
    ins...
    [0] feed
        phase: 'l', T: 340 K, P: 101325 Pa
        flow (kmol/hr): Water    20
                        Ethanol  10
    outs...
    [0] out_a
        phase: 'l', T: 340 K, P: 101325 Pa
        flow (kmol/hr): Water    20
                        Ethanol  10
    [1] out_b
        phase: 'l', T: 340 K, P: 101325 Pa
        flow (kmol/hr): Water    20
                        Ethanol  10
    [2] out_c
        phase: 'l', T: 340 K, P: 101325 Pa
        flow (kmol/hr): Water    20
                        Ethanol  10

    """
    _N_outs = 2
    _outs_size_is_fixed = False
    
    def _load_stream_links(self):
        feed, = self.ins
        duplicated_outs = self.outs
        for stream in duplicated_outs:
            stream.link_with(feed)
