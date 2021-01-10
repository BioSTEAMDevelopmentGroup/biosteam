# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import time
import numpy as np

__all__ = ('TicToc',)

class TicToc: # pragma: no coverage
    """Create a TicToc class with tic toc functions that measure elapsed time."""
    __slots__ = ['ID', 'record', '_start']

    def __init__(self, ID=None):
        self.ID = ID
        self.record = [] #: [list] elapsed times from tic toc functions
        self._start = None

    def tic(self):
        """Start timer."""
        # Marks the beginning of a time interval
        self._start = time.perf_counter()

    def toc(self):
        """Record time interval since last 'tic'."""
        # Appends time difference
        try: self.record.append(self.elapsed_time)
        except TypeError:
            if self._start is None:
                raise RuntimeError("Must run 'tic' before 'toc'.")    

    @property
    def elapsed_time(self):
        return time.perf_counter() - self._start

    @property
    def average(self):
        """The mean value of elapsed time"""
        return np.mean(self.record)

    def __repr__(self):
        ID = ': ' + self.ID if self.ID else ''
        return (f"<{type(self).__name__}{ID}, average={self.average:.2g} seconds>")

    def _info(self):
        ID = ': ' + self.ID if self.ID else ''
        record = np.array(self.record)
        average = self.average
        return (f"{type(self).__name__}{ID}\n" + 
                f" average: {average:.2g} second\n"
                f" record: {record}\n" )
    def show(self):
        print(self._info())
        