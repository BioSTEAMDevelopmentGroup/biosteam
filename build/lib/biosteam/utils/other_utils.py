# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 17:24:01 2018

This module includes arbitrary classes and functions.

@author: Guest Group
"""

from biosteam import np, Q_
import time

#from multiprocessing import Process
#import threading
#import _thread

__all__ = ('factor', 'checkbounds', 'approx2step', 'copy_attr', 'get_attr', 'Timer', 'run_in_parallel')

# %% Small functions

def factor(base_units, new_units):
    return Q_(1, base_units).to(new_units).magnitude

def checkbounds(x, bounds):
    lb, up = bounds
    return lb < x < up

def approx2step(val, x0, dx):
    """Approximate value, val, to closest increment/step, dx, starting from x0."""
    while True:
        if val < x0:
            break
        x0 += dx
    return x0

# %% For managing objects

def copy_attr(obj0, obj1, *attrs: str):
    """Copy specified attributes 'attrs' from obj1 to obj0."""
    for attr in attrs:
        setattr(obj0, attr, getattr(obj1, attr))


def get_attr(obj, *attrs: str) -> 'list[attributes]':
    """Get a list of attributes, 'attrs', from obj."""
    out = []
    for attr in attrs:
        out.append(getattr(obj, attr))
    return out


# %% Timer

class Timer:
    """Create a Timer class with tic toc functions that measure elapsed time."""
    __slots__ = ['tictoc', '_start']

    def __init__(self):
        self.tictoc = [] #: [list] elapsed times from tic toc functions

    def toc(self):
        """Record time interval since last 'tic' in self.tictoc."""
        # Appends time difference
        if self._start:
            self.tictoc.append(time.time() - self._start)
        else:
            raise Exception("Must run 'tic' before 'toc'.")

    def tic(self):
        """Start timer."""
        # Marks the beginning of a time interval
        self._start = time.time()

    @property
    def average(self):
        """The mean value of elapsed time"""
        return np.mean(self.tictoc)

    def __repr__(self):
        return (f"<{type(self).__name__}, average={self.average:.3g}>")

    def _info(self):
        return (f"{type(self).__name__}: \n" + 
                f" tictoc: {self.tictoc} \n" +
                f" average: {self.average}")
    def show(self):
        print(self._info())
        
        
# %% Iterative solvers

def PI_iteration(process_function, state_function, Kp, x0,
                 invtau, xtol, maxiter=150):
    """Iteratively solves a function using feedback PI control. The process is split in two parts, the process function and the state function. The process function depends on the process state and returns the output that must converge. The state function updates the process state accoring to the error from the process function and the PI settings (gain and integral time constant).


    **Parameters**

        state_function = updates the state of the process according to the input

        process_function = returns iterated output

        Kp = gain

        x0 = initial guess where process_function(x0) = x0

        bal_tau = inverse of integral time constant 

        xtol = margin of allowable error

    Returns:

        it = the number of iterations

    """
    run = True
    it = 0
    int_err = 0
    while run:
        x1 = process_function()
        err = x1 - x0
        if (err < 0 and int_err > 0) or (err > 0 and int_err < 0):
            int_err = 0
        int_err += invtau * err
        state_function(x1 + Kp*(err + int_err))
        x0 = x1
        it += 1
        if it > maxiter:
            raise Exception('PI_iteration did not converge')
            break
        run = abs(err) > xtol
    return it

# def PI_control(function, set_point, x0, Kp = 0.05, invtau = 0.5, xtol = 0.1, maxiter = 100):
#      """
#      PI Control Solver

#      Returns:

#           it = the number of iterations

#      Input arguments:

#           function = the process to be iterated

#           set_point = iteration Stops once function output reaches the set point within xtol

#           x0 = initial condition of function input

#           Kp = gain

#           invtau = inverse of integral time constant

#           xtol = margin of allowable error
#      """
#      err = xtol + 1
#      it = 0
#      int_err = 0
#      while abs(err) > xtol:
#           y = function(x0)
#           err = set_point - y
#           if (err < 0 and int_err > 0) or (err > 0 and int_err < 0):
#                int_err = 0
#           int_err += invtau * err
#           x0 = x0 + Kp*(err + int_err)
#           it += 1
#           if it > maxiter:
#                raise Exception('PI_solver did not converge')
#                break
#      return it

# %% Parallel processing and threading

from multiprocessing import Process
def run_in_parallel(fns, argss):
    """Run functions in parallel.
    
    **Parameters**
    
        fns: iterable[function] Functions to run in parallel
    
        argss: iterable[tuple] Arguments for functions
    
    """
    proc = []
    for fn, args in zip(fns, argss):
          p = Process(target=fn, args=args)
          p.start()
          proc.append(p)
    for p in proc:
          p.join()
#
# def input_with_timeout(prompt, timeout= 5):
#     timer = threading.Timer(timeout, _thread.interrupt_main)
#     astring = None
#     try:
#          timer.start()
#          astring = input(prompt)
#     except KeyboardInterrupt:
#          pass
#     timer.cancel()
#     return astring
