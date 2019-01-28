# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 17:24:01 2018

This module includes arbitrary classes and functions.

@author: Guest Group
"""

from biosteam import np
import time
import os

#from multiprocessing import Process
#import threading
#import _thread

__all__ = ('approx2step', 'copy_attr', 'get_attr', 'autodoc', 'TextManager', 'autopep8_all', 'replace_all', 'count_all', 'length_all', 'Timer', 'run_in_parallel')

# %% Small functions

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



# %% For making Sphinx docs

def autodoc(cls, path: str):
    """Make Sphinx documentation txt file in the specified path."""
    # Return if there is no docstring
    if cls.__doc__ is None:
        return
    name = cls.__name__

    # Skeleton doc
    header = f'{name}\n'
    header += '='*len(name) + '\n\n'
    module = f'.. module:: {cls.__module__}\n\n'
    autoclass = ".. autoclass:: {name}\n"
    members = "   :members:"
    skeleton_doc = header + module + autoclass + members

    # Set class name for document
    name = cls.__name__
    doc = skeleton_doc.replace('{name}', name)

    # Write file
    file = open(f"{path}/{name}.txt", 'w')
    file.write(doc)
    file.close()


# %% Text management

class TextManager:
    """Create a TextManager object that can read a document, store its content and update the document.

    **Parameters**

        document: [str] Name of the document

        directory: [str] Relative location of the document

    **Instance variables**

        content: [str] Content of the document

        document: [str] Name of the document
        
        original: [str] Original document

    """

    def __init__(self, document, directory=""):
        self.document = directory + document
        self.refresh()

    def refresh(self):
        """Reads the document and stores its files in content."""
        file = open(self.document, 'r+')
        self.content = file.read()
        self.original = self.content
        file.close()

    def to_utf8(self):
        self.content = self.content.encode('utf-8').strip().decode()

    def replace(self, old=[], new=[]):
        """Replace old strings in self.content with new.

        **Parameters**

            old: list[str] strings to be replaced
    
            new: list[str] new strings that replace old ones
        """
        if type(old) != list:
            old = [old]
        if type(new) != list:
            new = [new]
        if len(old) != len(new):
            raise ValueError('inputs must have the same length')
        for i in range(len(old)):
            self.content = self.content.replace(old[i], new[i])

    def autopep8(self):
        """Fix content to follow PEP 8."""
        from autopep8 import fix_code
        self.content = fix_code(self.content)

    def restart(self):
        """Undo any changes made to content."""
        self.content = self.original

    def update(self):
        """Writes content in document"""
        file = open(self.document, 'w+')
        file.write(self.content)
        file.close()

    def __repr__(self):
        return (f'<Text_manager: {self.document}>')
    
    def show(self):
        print(self)
        print(self.content)



def autopep8_all(directory=".\\"):
    """autopep8 all files in new for every document and subfolder document ending with .py and .ipynb in the specified directory"""
    fileIDs = [os.fsdecode(file) for file in os.listdir(directory)]
    for fileID in fileIDs:
        if '.' not in fileID and '__' not in fileID and fileID != 'Sphinx':
            autopep8_all(directory + fileID + '\\')
        if fileID.endswith('.py') or fileID.endswith('.ipynb'):
            x = TextManager(fileID, directory)
            x.autopep8()
            x.update()

def get_files(directory):
    fileIDs = []
    for file in os.listdir(directory):
        try:
            file = os.fsdecode(file)
            fileIDs.append(file)
        except:
            pass
    return fileIDs

def replace_all(old, new, directory=".\\\\"):
    """Replaces strings in old with strings in new for every document and subfolder document ending with .py and .ipynb in the specified directory"""
    try:
        fileIDs = get_files(directory)
    except:
        return
    
    for fileID in fileIDs:
        if '.' not in fileID and '__' not in fileID and fileID != 'Sphinx':
            replace_all(old, new, directory + fileID + '\\')
        if fileID.endswith('.py') or fileID.endswith('.ipynb') or fileID.endswith('.rst') or fileID.endswith('.txt'):
            try:
                x = TextManager(fileID, directory)
                x.replace(old, new)
                x.update()
            except UnicodeDecodeError:
                print(f'{repr(x)} could not read file.')
                continue
            except NotADirectoryError:
                print(f'{repr(x)} could not read file.')
                continue


def length_all(directory='.\\\\'):
    N = 0
    fileIDs = get_files(directory)
    for fileID in fileIDs:
        if '.' not in fileID and '__' not in fileID and fileID != 'Sphinx':
            N += length_all(directory + fileID + '\\')
        if fileID.endswith('.py') or fileID.endswith('.ipynb'):
            file = open(directory + fileID, 'r+')
            content = file.read()
            N += len(content)
            file.close()
    return N


def count_all(string, directory='.\\\\'):
    N_lines = 0
    fileIDs = get_files(directory)
    for fileID in fileIDs:
        if '.' not in fileID and fileID != 'Sphinx':
            N_lines += count_all(string, directory + fileID + '\\')
        if fileID.endswith('.py') or fileID.endswith('.ipynb'):
            file = open(directory + fileID, 'r+')
            content = file.read()
            N_lines += content.count(string)
            file.close()
    return N_lines

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
