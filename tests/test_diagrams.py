# -*- coding: utf-8 -*-
"""
Created on Sat Dec  5 13:25:23 2020

@author: yrc2
"""
import pytest
import biosteam as bst
import os

folder = os.path.dirname(__file__)
folder = os.path.join(folder, 'Diagram Sources')
DISPLAY = False
COMPARE_CONSISTENCY = False

def save_diagrams():
    import biorefineries.sugarcane as sc
    f = sc.C201.diagram(display=DISPLAY)
    file = os.path.join(folder, 'Clarifier.txt')
    expected_source = open(file, 'w')
    expected_source.write(f.source)
    for kind in ('thorough', 'cluster', 'surface', 'minimal'):
        f = sc.sugarcane_sys.diagram(kind, display=DISPLAY)
        file = os.path.join(folder, f'sugarcane {kind}.txt')
        expected_source = open(file, 'w')
        expected_source.write(f.source)
        expected_source.close()

# TODO: Find out what makes diagrams inconsistent in graphviz.
# For now, just make sure they run.
def test_unit_diagram():
    import biorefineries.sugarcane as sc
    sc.load() # Reload system to make sure all is consistent
    f = sc.C201.diagram(display=DISPLAY)
    if COMPARE_CONSISTENCY:
        file = os.path.join(folder, 'Clarifier.txt')
        with open(file) as f_expected:
            expected_source = f_expected.read()
            assert set(f.source.split()) == set(expected_source.split())
    bst.process_tools.default()

def test_system_thorough_diagram():
    import biorefineries.sugarcane as sc
    f = sc.sugarcane_sys.diagram('thorough', display=DISPLAY)
    file = os.path.join(folder, 'sugarcane thorough.txt')
    if COMPARE_CONSISTENCY:
        with open(file) as f_expected:
            expected_source = f_expected.read()
            assert set(f.source.split()) == set(expected_source.split())
    bst.process_tools.default()

def test_system_cluster_diagram():
    import biorefineries.sugarcane as sc
    f = sc.sugarcane_sys.diagram('cluster', display=DISPLAY)
    file = os.path.join(folder, 'sugarcane cluster.txt')
    if COMPARE_CONSISTENCY:
        with open(file) as f_expected:
            expected_source = f_expected.read()
            assert set(f.source.split()) == set(expected_source.split())
    bst.process_tools.default()

# TODO: Test does not work CI due to FileNotFoundError, but does work locally
# def test_saving():
#     import biorefineries.sugarcane as sc
#     with pytest.raises(TypeError):
#         sc.sugarcane_sys.diagram(file=1.11)
#     sc.sugarcane_sys.diagram(file=os.path.join(folder, 'sugarcane'))
#     os.remove(os.path.join(folder, 'sugarcane.svg'))

# TODO: Find out why this test is not working
# def test_system_surface_diagram():
#     import biorefineries.sugarcane as sc
#     f = sc.sugarcane_sys.diagram('surface', display=DISPLAY)
#     file = os.path.join(folder, 'sugarcane surface.txt')
#     with open(file) as f_expected:
#         expected_source = f_expected.read()
#         assert set(f.source.split()) == set(expected_source.split())
#     bst.process_tools.default()

# def test_system_minimal_diagram():
#     import biorefineries.sugarcane as sc
#     f = sc.sugarcane_sys.diagram('minimal', display=DISPLAY)
#     file = os.path.join(folder, 'sugarcane minimal.txt')
#     with open(file) as f_expected:
#         expected_source = f_expected.read()
#         assert set(f.source.split()) == set(expected_source.split())
#     bst.process_tools.default()

def test_ascii_flowsheet_interface():
    # The layout of ASCII flowsheets is tested in test_ascii_flowsheet.py;
    # here we test the interface with biosteam objects.
    import io
    import tempfile
    from contextlib import redirect_stdout
    bst.settings.set_thermo(['Water'], cache=True)
    with bst.Flowsheet('ascii_flowsheet_interface'):
        recycle = bst.Stream('recycle')
        M1 = bst.Mixer('M1', [bst.Stream('feed', Water=1), recycle])
        S1 = bst.Splitter('S1', M1-0, ['product', recycle], split=0.5)
        sys = bst.System.from_units('sys', [M1, S1])
        expected = (
            "       -------\n"
            "       |     |\n"
            "       v     |\n"
            "feed->[M1]->[S1]->produ"
        )
        assert sys.diagram(format='ascii', display=False) == expected
        for kind in ('cluster', 'surface', 'stage'):
            assert sys.diagram(kind, format='ascii', display=False) == expected
        assert sys.diagram('minimal', format='ascii', display=False) == (
            "feed->[SYS]->produ"
        )
        # Diagram of a unit (streams to other units are feeds and products)
        assert M1.diagram(format='ascii', display=False) == (
            " feed-+>[M1]->s1\n"
            "      |\n"
            "recyc-+"
        )
        # Display and save
        output = io.StringIO()
        with redirect_stdout(output): sys.diagram(format='ascii')
        assert output.getvalue() == expected + '\n'
        with tempfile.TemporaryDirectory() as folder:
            file = os.path.join(folder, 'flowsheet.txt')
            sys.diagram(format='ascii', file=file)
            with open(file, encoding='utf-8') as f: assert f.read() == expected + '\n'
        # Preferences
        with bst.preferences.temporary():
            bst.preferences.graphviz_format = 'ascii'
            assert sys.diagram(display=False) == expected
        # Maximum lengths of unit and stream labels
        assert sys.diagram(format='ascii', display=False,
                           unit_label_length=1, stream_label_length=3) == (
            "      ------\n"
            "      |    |\n"
            "      v    |\n"
            "fee->[M]->[S]->pro"
        )
        assert M1.diagram(format='ascii', display=False,
                          unit_label_length=1, stream_label_length=3) == (
            "fee-+>[M]->s1\n"
            "    |\n"
            "rec-+"
        )

if __name__ == '__main__':
    # save_diagrams()
    test_unit_diagram()
    # test_system_surface_diagram()
    test_system_thorough_diagram()
    test_system_cluster_diagram()
    # test_system_minimal_diagram()
    # test_saving()
    test_ascii_flowsheet_interface()