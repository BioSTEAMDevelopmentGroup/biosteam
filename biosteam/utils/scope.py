# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2024, Joy Zhang <joycheung1994@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

"""
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline as ius
from warnings import warn

__all__ = ('Scope', 'SystemScope')

class Scope():
    """
    A general tracker of attributes of a subject during dynamic simulations.

    Parameters
    ----------
    subject : 
        The subject to scope during dynamic simulation.
    variables : Iterable[str]
        The attributes to track.
    header : list of 2-tuple or :class:`pandas.MultiIndex`, optional
        The header for the tracked time-series data. When none specified, will
        be auto-generated based on the defined variables to track.

    See Also
    --------
    `pandas.MultiIndex <https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.MultiIndex.html>`_

    """    
    def __init__(self, subject, variables, header=None, **kwargs):
        self.subject = subject
        self._ts = []
        self._header = header
        rcd = {}
        for var in variables:
            if hasattr(subject, var): rcd[var] = []
            else: warn(f'Variable {var} ignored in {self.__repr__()} because '
                       f'{self.subject} has no attribute {var}.')
        self._record = rcd
        for k, v in kwargs:
            setattr(self, k, v)

    def getter(self, variable):
        """A function to receive the attribute or variable of interest."""
        return getattr(self.subject, variable).copy()
        
    def __call__(self, t):
        """Tracks the variables at time t."""
        self._ts.append(t)
        for var, log in self._record.items():
            log.append(self.getter(var))
    
    def reset_cache(self):
        """Clears all recorded data."""
        self._ts = []
        self._record = {var:[] for var in self._record.keys()}
    
    def __repr__(self):
        return f'<Scope: {self.subject.ID}>'

    def pop(self):
        """Removes the last tracked time point."""
        self._ts.pop()
        for log in self._record.values():
            log.pop()
    
    def _n_cols(self, make_header=False):
        n = []
        isa = isinstance
        if make_header: 
            names = []
            for var in self._record.keys():
                data = self.getter(var)
                if isa(data, (float, int, str)): ni = 1
                else: ni = len(data)
                n.append(ni)
                names += [f'{var}_{i}' for i in range(ni)]
            return n, names
        else:
            for var in self._record.keys():
                data = self.getter(var)
                if isa(data, (float, int, str)): n.append(1)
                else: n.append(len(data))
            return n
    
    @property
    def header(self):
        """[list of 2-tuple] Header for the tracked time-series data."""
        return self._header
    @header.setter
    def header(self, hd):
        if hd is None:
            ncol, names = self._n_cols(True)
            hd = list(zip([self.subject.ID]*sum(ncol), names))
        else:
            if len(hd) != sum(self._n_cols()):
                raise ValueError(f'Header {hd} has the wrong size {len(hd)}, '
                                 f'it should have length = {sum(self._n_cols)}')
            if len(hd[0]) != 2:
                raise ValueError(f'A list of 2-tuple or a 2-level pandas.MultiIndex '
                                 f'is expected but got an iterable of {len(hd[0])}-tuple')
        self._header = hd
    
    @property
    def record(self):
        """[numpy.ndarray] The tracked time-series data of the variables of interest."""
        data = np.hstack([np.array(v) if len(np.array(v).shape) == 2 \
                          else np.array(v).reshape((len(v), 1))\
                          for v in self._record.values()])
        return data
    
    @property
    def time_series(self):
        """[numpy.1darray] The tracked time points."""
        return np.array(self._ts)
    
    def plot_time_series(self, variable):
        """plot the time series data of a single variable of interest"""
        fig, ax = plt.subplots(figsize=(8, 4.5))
        t = self.time_series
        ys = np.array(self._record[variable])
        if len(ys.shape) == 1:
            ax.plot(t, ys, '-o')
        else:
            for i, y in enumerate(ys.T):
                ax.plot(t, ys, '-o', label=f'#{i}')
            ax.legend(loc='best')
        ax.set(xlabel='Time [d]', ylabel=variable)
        return fig, ax
    
class SystemScope():
    """
    A class for keeping track of time-series data generated during dynamic 
    simulation of a system.

    Parameters
    ----------
    system : :class:`System`
        The `System` object to track.
    *subjects : 
        The subjects of interest to scope during dynamic simulation, e.g., 
        a :class:`Unit` object or a :class:`Stream` object with dynamic 
        state variables.
    interpolator : callable, optional
        An interpolation method that takes in time-series data and returns
        an interpolant. Used to export the data at certain time points. 
        When none specified, will use :class:`scipy.interpolate.InterpolatedUnivariateSpline` 
        with k=1 (i.e., linear) and will raise error when trying to extrapolate.
        
    See Also
    --------
    `scipy.interpolate.InterpolatedUnivariateSpline <https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.InterpolatedUnivariateSpline.html>`_
    """
    def __init__(self, system, *subjects, interpolator=None, **kwargs):
        self.system = system
        self.subjects = subjects
        self._method = interpolator or ius
        self._ts = []
        self.sol = None
        self.sol_header = system._state_header
        for k,v in kwargs:
            setattr(self, k, v)
    
    def __call__(self, t):
        ts = self._ts
        while len(ts) > 0 and t <= ts[-1]:
            ts.pop()
            for s in self.subjects: s.scope.pop()
        ts.append(t)
        for s in self.subjects:
            s.scope(t)
    
    def reset_cache(self):
        '''Clears all recorded data.'''
        self._ts = []
        self.sol = None
        for s in self.subjects:
            s.scope.reset_cache()
    
    def __repr__(self):
        return f'<SystemScope: {self.system.ID}>'

    @property
    def subjects(self):
        """The subjects of interest to scope during dynamic simulation."""
        return self._subjects
    @subjects.setter
    def subjects(self, sjs):
        for s in sjs:
            if not hasattr(s, 'scope'):
                try: # the unit uses `_compile_AE`, but needs to `_init_dynamic` as its state isn't copied from its influent
                    s._init_dynamic()
                except: raise AttributeError(f"{s} has no attribute 'scope'")
            elif not isinstance(s.scope, Scope):
                raise TypeError(f'{s}.scope must be a {Scope} object')
        self._subjects = sjs

    @property
    def sol_header(self):
        """[list of 2-tuple or :class:`pandas.MultiIndex`]The header for the solution data returned by the IVP solver."""
        return self._sol_header
    
    @sol_header.setter
    def sol_header(self, hd):
        units = self.system.units
        if hd is None:
            hd = [('-', 't [d]')]
            for u in units:
                if u.hasode:
                    hd += list(zip([u.ID]*len(u._state_header), u._state_header))                
            self._sol_header = pd.MultiIndex.from_tuples(hd, names=['unit', 'variable'])
        else:
            self._sol_header = hd
            
    @property
    def time_series(self):
        """[numpy.1darray] The tracked time points."""
        return np.array(self._ts)
    
    def _get_records(self):
        data = [s.scope.record for s in self.subjects]
        # each row is one "variable", each column corresponds to one time point
        return np.hstack(data).T
    
    def _get_headers(self):
        headers = [('-', 't [d]')]
        isa = isinstance
        for s in self.subjects:
            if not isa(s.scope.header, list):
                headers += list(s.scope.header)
            else: 
                headers += s.scope.header
        return pd.MultiIndex.from_tuples(headers, names=['ID', 'variable'])

    def _interpolate_eval(self, t_arr, y_arrs, t_eval, **interpolation_kwargs):
        f = self._method
        y_eval = np.empty((y_arrs.shape[0], len(t_eval)))
        if f is ius and 'k' not in interpolation_kwargs.keys(): 
            interpolation_kwargs['k'] = 1
        for i, y in enumerate(y_arrs):            
            intpl = f(t_arr, y, **interpolation_kwargs)
            if max(t_eval) > max(t_arr):
                raise RuntimeError(f'Extrapolation is tempted! t_eval must be '
                                   f'within the range of [{min(t_arr)}, {max(t_arr)}].')
            y_eval[i,:] = intpl(t_eval)
        return np.vstack([t_eval, y_eval])
    
    def export(self, path='', t_eval=None, **interpolation_kwargs):
        """Exports recorded time-series data to given path."""
        ts = self.time_series
        ys = self._get_records()
        if t_eval is None:
            data = np.vstack([ts, ys])
        else:
            data = self._interpolate_eval(ts, ys, t_eval, **interpolation_kwargs)
        df = pd.DataFrame(data.T, columns=self._get_headers())
        if path:
            file, ext = path.rsplit('.', 1)
            if ext == 'npy': 
                np.save(path, data.T)
            elif ext in ('xlsx', 'xls'):
                df.to_excel(path)
            elif ext == 'csv':
                df.to_csv(path)
            elif ext == 'tsv':
                df.to_csv(path, sep='\t')
            else:
                raise ValueError('Only support file extensions of ".npy", '
                                 '".xlsx", ".xls", ".csv", and ".tsv", '
                                 f'not .{ext}.')
        else: return df                