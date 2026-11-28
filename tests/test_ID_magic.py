# -*- coding: utf-8 -*-
"""
Created on Sat Dec  5 13:25:23 2020

@author: yrc2
"""
import biosteam as bst

def test_variable_assignment():
    bst.settings.set_thermo(['Water'], cache=True)
    bst.main_flowsheet.clear()
    for i in (bst.Stream, bst.Unit, bst.System): i.ticket_numbers.clear()
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
    
    s_sum = bst.Stream.sum([bst.Stream() for i in range(4)])
    assert s_sum.ID == 's_sum'
    s_sum_copy = s_sum.copy()
    assert s_sum_copy.ID == 's_sum_copy'
    
def test_commas():
    bst.settings.set_thermo(['Water'], cache=True)
    bst.main_flowsheet.clear()
    for i in (bst.Stream, bst.Unit, bst.System): i.ticket_numbers.clear()
    feed, = bst.Stream(),
    assert feed.ID == 'feed'
    (mixer,) = (bst.Mixer(),)
    assert mixer.ID == 'mixer'
    assert mixer.outs[0].ID == 's1' # Unit outlets take default IDs
    sys, = [bst.System()]
    assert sys.ID == 'sys'
    feed, product = (bst.Stream(), bst.Stream())
    assert (feed.ID, product.ID) == ('feed', 'product')
    water, ethanol, *_ = bst.Stream(), bst.Stream(), bst.Stream()
    other, = _
    assert (water.ID, ethanol.ID) == ('water', 'ethanol')
    assert other.ID == 's2' # Starred targets are not inferred

def test_unpacking_in_function():
    bst.settings.set_thermo(['Water'], cache=True)
    bst.main_flowsheet.clear()
    for i in (bst.Stream, bst.Unit, bst.System): i.ticket_numbers.clear()
    def f():
        feed, product = (
            bst.Stream(),
            bst.Stream(),
        )
        return feed, product
    x, y = f()
    assert (x.ID, y.ID) == ('feed', 'product')

def test_nested_function():
    bst.settings.set_thermo(['Water'], cache=True)
    bst.main_flowsheet.clear()
    for i in (bst.Stream, bst.Unit, bst.System): i.ticket_numbers.clear()
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
    for i in (bst.Stream, bst.Unit, bst.System): i.ticket_numbers.clear()
    feeds = [bst.Stream() for i in range(2)]
    assert [i.ID for i in feeds] == ['s1', 's2']
    mixers = [bst.Mixer() for i in range(2)]
    assert [i.ID for i in mixers] == ['M1', 'M2']
    systems = [bst.System() for i in range(2)]
    assert [i.ID for i in systems] == ['SYS1', 'SYS2']

def test_multiline():
    bst.settings.set_thermo(['Water'], cache=True)
    bst.main_flowsheet.clear()
    for i in (bst.Stream, bst.Unit, bst.System): i.ticket_numbers.clear()
    feed = bst.Stream(
        Water=1,
    )
    assert feed.ID == 'feed'
    mixer = (
        bst.Mixer(ins=feed)
    )
    assert mixer.ID == 'mixer'
    product = \
        mixer.outs[0].copy()
    assert product.ID == 'product'
    with bst.System(
        ) as sys: pass
    assert sys.ID == 'sys'
    def f():
        splitter = bst.Splitter(
            ins=product, split=0.5,
        )
        return splitter
    assert f().ID == 'splitter'

if __name__ == '__main__':
    test_variable_assignment()
    test_commas()
    test_nested_function()
    test_comprehension()
    test_unpacking_in_function()
    test_multiline()