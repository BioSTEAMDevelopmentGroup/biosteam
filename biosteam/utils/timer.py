# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import time
import numpy as np

__all__ = ('Timer',)

class TimerOffset:
    __slots__ = ['_timer', '_start']
    
    def __init__(self, timer):
        self._timer = timer
        
    def __enter__(self):
        self._start = time.perf_counter()
        
    def __exit__(self, type, exception, traceback):
        self._timer._start += (time.perf_counter() - self._start)
        if exception: raise exception


class Timer: # pragma: no coverage
    """Create a Timer class with functions that measure elapsed time."""
    __slots__ = ['ID', 'record', '_start']

    def __init__(self, ID=None):
        self.ID = ID
        self.record = [] #: [list] elapsed times from tic toc functions
        self._start = None

    def offset(self):
        if self._start is None: 
            raise RuntimeError('timer has not been started yet')
        return TimerOffset(self)

    def start(self):
        """Start timer."""
        # Marks the beginning of a time interval
        self._start = time.perf_counter()
        return self._start

    def mock_start(self):
        self._start = 0
        return 0

    def mock_measure(self):
        record = self.record
        if record: record.append(record[-1] + 1)
        else: record.append(1)

    def measure(self, record=True):
        """Record time interval since last start time."""
        # Appends time difference
        try: elapsed_time = self.elapsed_time
        except TypeError:
            if self._start is None:
                raise RuntimeError("Must run 'tic' before 'toc'.")    
        if record: self.record.append(elapsed_time)
        return elapsed_time

    @property
    def elapsed_time(self):
        return time.perf_counter() - self._start

    @property
    def mean(self):
        """The mean value of elapsed time."""
        return np.mean(self.record)

    @property
    def standard_deviation(self):
        """The standard deviation value of elapsed time."""
        return np.std(self.record)

    std = standard_deviation

    @property
    def summary(self):
        """Return statistical summary of the elapsed time as a string."""
        return f"{self.mean:.2g} \u00B1 {self.std:.2g} s"

    def __repr__(self):
        return f"{type(self).__name__}({self.ID})"

    def _info(self):
        ID = ': ' + self.ID if self.ID else ''
        record = np.array(self.record)
        return (f"{type(self).__name__}{ID}\n" 
                f" mean: {self.mean:.2g} second\n"
                f" std: {self.std:.2g} second\n"
                f" record: {record}\n" )
    
    def show(self):
        print(self._info())
    
    _ipython_display_ = show
        