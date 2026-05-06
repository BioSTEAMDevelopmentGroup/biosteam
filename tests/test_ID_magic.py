# -*- coding: utf-8 -*-
"""
Created on Sat Dec  5 13:25:23 2020

@author: yrc2
"""
import biosteam as bst

def test_variable_assignment():
    bst.settings.set_thermo(['Water'], cache=True)
    bst.main_flowsheet.clear()
    feed = bst.Stream()
    assert feed.ID == 'feed'
    mixer = bst.Mixer()
    assert mixer.ID == 'mixer'
    sys = bst.System()
    assert sys.ID == 'sys'
    with bst.System() as sys: pass
    assert sys.ID == 'sys'
    
    @bst.SystemFactory
    def f(ins, outs): pass
    
    sys = f()
    assert sys.ID == 'sys'
    
def test_commas():
    bst.settings.set_thermo(['Water'], cache=True)
    bst.main_flowsheet.clear()
    feed, = bst.Stream(),
    assert feed.ID == 's1'
    (mixer,) = (bst.Mixer(),)
    assert mixer.ID == 'M1'
    sys, = [bst.System()]
    assert sys.ID == 'SYS1'

def test_nested_function():
    bst.settings.set_thermo(['Water'], cache=True)
    bst.main_flowsheet.clear()
    f = lambda obj: obj
    feed = f(bst.Stream())
    assert feed.ID == 's1'
    mixer = f(bst.Mixer())
    assert mixer.ID == 'M1'
    sys = f(bst.System())
    assert sys.ID == 'SYS1'
    with f(bst.System()) as sys: pass
    assert sys.ID == 'SYS2'

def test_comprehension():
    bst.settings.set_thermo(['Water'], cache=True)
    bst.main_flowsheet.clear()
    feeds = [bst.Stream() for i in range(2)]
    assert [i.ID for i in feeds] == ['s1', 's2']
    mixers = [bst.Mixer() for i in range(2)]
    assert [i.ID for i in mixers] == ['M1', 'M2']
    systems = [bst.System() for i in range(2)]
    assert [i.ID for i in systems] == ['SYS1', 'SYS2']

if __name__ == '__main__':
    test_variable_assignment()
    test_commas()
    test_nested_function()
    test_comprehension()