# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 11:28:33 2019

@author: yoelr
"""
import time
import numpy as np
from .other_utils import factor

__all__ = ('Timer',)

class Timer:
    """Create a Timer class with tic toc functions that measure elapsed time."""
    __slots__ = ['ID', 'tictoc', '_start']

    def __init__(self, ID=None):
        self.ID = ID
        self.tictoc = [] #: [list] elapsed times from tic toc functions
        self._start = None

    def toc(self):
        """Record time interval since last 'tic' in self.tictoc."""
        # Appends time difference
        if self._start is None:
            raise Exception("Must run 'tic' before 'toc'.")
        else:
            self.tictoc.append(time.clock() - self._start)

    def tic(self):
        """Start timer."""
        # Marks the beginning of a time interval
        self._start = time.clock()

    @property
    def average(self):
        """The mean value of elapsed time"""
        return np.mean(self.tictoc)

    def __repr__(self):
        ID = ': ' + self.ID if self.ID else ''
        return (f"<{type(self).__name__}{ID}, average={self.average:.2g} seconds>")

    def _info(self, units='second'):
        ID = ': ' + self.ID if self.ID else ''
        F = factor('second', units)
        tictoc = F * np.array(self.tictoc)
        average = F * self.average
        return (f"{type(self).__name__}{ID}\n" + 
                f" average: {average:.2g} {units}\n"
                f" tictoc: {tictoc}\n" )
    def show(self, units='second'):
        print(self._info(units))
        
        