# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 17:24:01 2018

This module includes arbitrary classes and functions.

@author: Guest Group
"""
from biosteam import Q_

#from multiprocessing import Process
#import threading
#import _thread

__all__ = ('factor', 'checkbounds', 'approx2step', 'run_in_parallel', 'strtuple', 'function')


# %% Number functions

def factor(base_units, new_units):
    if base_units == new_units: return 1
    else: return Q_(1, base_units).to(new_units).magnitude

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


# %% String functions

function = type(checkbounds)

def strtuple(iterable):
    """Return string of all items in the tuple""" 
    string = ''
    for i in iterable:
        if isinstance(i , function):
            string += i.__name__ + ', '
        else:
            string += str(i) + ', '
    string = string.rstrip(', ')
    string = '(' + string + ')'
    return string
        
        
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
