# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Tests for the layout of ASCII flowsheets (`System.diagram(format='ascii')`).

Tests are in order of increasing complexity. Each test covers one new feature;
later tests only combine features already tested. Flowsheets are given as
nodes and edges, so biosteam is not loaded.

"""
import os
import importlib.util

path = os.path.join(os.path.dirname(__file__), os.pardir, 'biosteam', 'digraph', 'ascii.py')
spec = importlib.util.spec_from_file_location('ascii_flowsheet', path)
ascii_flowsheet = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ascii_flowsheet)

def assert_flowsheet(nodes, edges, expected):
    """
    Assert the ASCII flowsheet of the given nodes (unit IDs in processing
    order) and edges (source, source index, label, sink index, sink; None for
    feeds and products) is as expected.
    """
    labels = [i[:2].upper() for i in nodes]
    assert len(set(labels)) == len(labels), 'each unit must have a different label'
    actual = ascii_flowsheet.ascii_flowsheet(list(zip(nodes, labels)), edges)
    expected = expected.strip('\n')
    assert actual == expected, f"\nactual:\n{actual}\n\nexpected:\n{expected}"

def test_single_unit():
    nodes = ['P1']
    edges = [
        (None, 0, 'feed', 0, 'P1'),
        ('P1', 0, 'product', 0, None),
    ]
    assert_flowsheet(nodes, edges, """
feed->[P1]->produ
""")

def test_stream_labels():
    """Feeds and products are labeled by the first 5 letters of their ID followed by any digits; repeated labels are numbered."""
    nodes = ['M1']
    edges = [
        (None, 0, 'sugarcane_juice', 0, 'M1'),
        (None, 0, 'water2', 1, 'M1'),
        ('M1', 0, 'cooling_water', 0, None),
        ('M1', 1, 'cooling_tower', 0, None),
    ]
    assert_flowsheet(nodes, edges, """
 sugar-+>[M1]+->cooli1
       |     |
water2-+     +->cooli2
""")

def test_units_in_series():
    nodes = ['P1', 'P2', 'P3']
    edges = [
        (None, 0, 'feed', 0, 'P1'),
        ('P1', 0, '', 0, 'P2'),
        ('P2', 0, '', 0, 'P3'),
        ('P3', 0, 'product', 0, None),
    ]
    assert_flowsheet(nodes, edges, """
feed->[P1]->[P2]->[P3]->produ
""")

def test_multiple_feeds():
    """Feeds to the same unit merge through a manifold."""
    nodes = ['M1']
    edges = [
        (None, 0, 'feed1', 0, 'M1'),
        (None, 0, 'feed2', 1, 'M1'),
        (None, 0, 'feed3', 2, 'M1'),
        ('M1', 0, 'product', 0, None),
    ]
    assert_flowsheet(nodes, edges, """
feed1-+
      |
feed2-+>[M1]->produ
      |
feed3-+
""")

def test_two_feeds():
    """Junctions are never next to each other."""
    nodes = ['M1']
    edges = [
        (None, 0, 'feed1', 0, 'M1'),
        (None, 0, 'feed2', 1, 'M1'),
        ('M1', 0, 'product', 0, None),
    ]
    assert_flowsheet(nodes, edges, """
feed1-+>[M1]->produ
      |
feed2-+
""")

def test_multiple_products():
    """Products from the same unit split through a manifold."""
    nodes = ['S1']
    edges = [
        (None, 0, 'feed', 0, 'S1'),
        ('S1', 0, 'product1', 0, None),
        ('S1', 1, 'product2', 0, None),
    ]
    assert_flowsheet(nodes, edges, """
feed->[S1]+->produ1
          |
          +->produ2
""")

def test_unconnected_sections():
    """Sections without connections between them are drawn separately."""
    nodes = ['P1', 'P2', 'P3']
    edges = [
        (None, 0, 'feed1', 0, 'P1'),
        ('P1', 0, 'product1', 0, None),
        (None, 0, 'feed2', 0, 'P2'),
        ('P2', 0, '', 0, 'P3'),
        ('P3', 0, 'product2', 0, None),
    ]
    assert_flowsheet(nodes, edges, """
feed2->[P2]->[P3]->produ2

feed1->[P1]->produ1
""")

def test_feed_to_downstream_unit():
    """Feeds to a unit that already has an inlet from upstream come from below."""
    nodes = ['P1', 'M1']
    edges = [
        (None, 0, 'feed1', 0, 'P1'),
        ('P1', 0, '', 0, 'M1'),
        (None, 0, 'feed2', 1, 'M1'),
        ('M1', 0, 'product', 0, None),
    ]
    assert_flowsheet(nodes, edges, """
feed1->[P1]->[M1]->produ
              ^
              |
          feed2
""")

def test_split_to_two_units():
    """Outlets to units split through a manifold."""
    nodes = ['S1', 'P1', 'P2']
    edges = [
        (None, 0, 'feed', 0, 'S1'),
        ('S1', 0, '', 0, 'P1'),
        ('S1', 1, '', 0, 'P2'),
        ('P1', 0, 'product1', 0, None),
        ('P2', 0, 'product2', 0, None),
    ]
    assert_flowsheet(nodes, edges, """
feed->[S1]+->[P1]->produ1
          |
          -->[P2]->produ2
""")

def test_merge_from_two_units():
    """Inlets from units merge at a junction before the arrow."""
    nodes = ['P1', 'P2', 'M1']
    edges = [
        (None, 0, 'feed1', 0, 'P1'),
        (None, 0, 'feed2', 0, 'P2'),
        ('P1', 0, '', 0, 'M1'),
        ('P2', 0, '', 1, 'M1'),
        ('M1', 0, 'product', 0, None),
    ]
    assert_flowsheet(nodes, edges, """
feed1->[P1]-+>[M1]->produ
            |
feed2->[P2]--
""")

def test_parallel_streams():
    """Parallel streams between the same units are drawn as one line."""
    nodes = ['A1', 'B1']
    edges = [
        (None, 0, 'feed', 0, 'A1'),
        ('A1', 0, '', 0, 'B1'),
        ('A1', 1, '', 1, 'B1'),
        ('B1', 0, 'product', 0, None),
    ]
    assert_flowsheet(nodes, edges, """
feed->[A1]->[B1]->produ
""")

def test_manifolds_in_same_column():
    """Manifolds of units in the same column do not share a bus."""
    nodes = ['X1', 'Y1', 'M1']
    edges = [
        (None, 0, 'feed1', 0, 'X1'),
        (None, 0, 'feed2', 0, 'Y1'),
        ('X1', 0, 'product1', 0, None),
        ('X1', 1, '', 0, 'M1'),
        ('Y1', 0, 'product2', 0, None),
        ('Y1', 1, '', 1, 'M1'),
        ('M1', 0, 'product', 0, None),
    ]
    assert_flowsheet(nodes, edges, """
           +->produ2
           |
feed2->[Y1]+--
             |
feed1->[X1]+-+>[M1]->produ
           |
           +->produ1
""")

def test_manifold_beside_label():
    """Manifold buses pass beside (not through) the labels of other units."""
    nodes = ['S1', 'F1', 'H1', 'G1', 'P1']
    edges = [
        (None, 0, 'feed', 0, 'S1'),
        ('S1', 0, '', 0, 'F1'),
        ('S1', 1, 'fiber_fines', 0, None),
        ('F1', 0, '', 0, 'H1'),
        ('F1', 1, '', 0, 'P1'),
        ('H1', 0, '', 0, 'G1'),
        ('G1', 0, 'fines', 0, None),
        ('G1', 1, 'product', 0, None),
        ('P1', 0, 'side', 0, None),
    ]
    assert_flowsheet(nodes, edges, """
feed->[S1]+---->[F1]+->[H1]->[G1]+->fines
          |         |            |
          +->fiber  |            +->produ
                    -->[P1]->side
""")

def test_bypass_over_manifold():
    """Branches that skip units go above; products go below."""
    nodes = ['B1', 'C1', 'D1']
    edges = [
        (None, 0, 'feed', 0, 'B1'),
        ('B1', 0, '', 0, 'C1'),
        ('B1', 1, '', 1, 'D1'),
        ('C1', 0, '', 0, 'D1'),
        ('C1', 1, 'p1', 0, None),
        ('C1', 2, 'p2', 0, None),
        ('D1', 0, 'product', 0, None),
    ]
    assert_flowsheet(nodes, edges, """
          ----------
          |        |
feed->[B1]+->[C1]+-+>[D1]->produ
                 |
                 +->p1
                 |
                 +->p2
""")

def test_products_on_separate_rows():
    """Products are never on rows right next to each other."""
    nodes = ['X1', 'G1', 'P1']
    edges = [
        (None, 0, 'feed', 0, 'X1'),
        ('X1', 0, '', 0, 'G1'),
        ('X1', 1, '', 0, 'P1'),
        ('G1', 0, 'product', 0, None),
        ('G1', 1, 'waste', 0, None),
        ('P1', 0, 'side', 0, None),
    ]
    assert_flowsheet(nodes, edges, """
feed->[X1]+->[G1]+->produ
          |      |
          |      +->waste
          |
          -->[P1]->side
""")

def test_feeds_on_separate_rows():
    """Feeds are never on rows right next to each other."""
    nodes = ['M1', 'M2', 'P1']
    edges = [
        (None, 0, 'feed1', 0, 'M1'),
        (None, 0, 'feed2', 1, 'M1'),
        (None, 0, 'water', 0, 'M2'),
        ('M1', 0, '', 0, 'P1'),
        ('M2', 0, '', 1, 'P1'),
        ('P1', 0, 'product', 0, None),
    ]
    assert_flowsheet(nodes, edges, """
feed1-+>[M1]-+>[P1]->produ
      |      |
feed2-+      |
             |
 water->[M2]--
""")

def test_simple_recycle():
    """Recycles go around the main line."""
    nodes = ['M1', 'S1']
    edges = [
        (None, 0, 'feed', 0, 'M1'),
        ('M1', 0, '', 0, 'S1'),
        ('S1', 0, 'product', 0, None),
        ('S1', 1, '', 1, 'M1'),
    ]
    assert_flowsheet(nodes, edges, """
       -------
       |     |
       v     |
feed->[M1]->[S1]->produ
""")

def test_recycle_through_unit():
    """Units in the return path of a recycle loop go on their own row."""
    nodes = ['M1', 'S1', 'H1']
    edges = [
        (None, 0, 'feed', 0, 'M1'),
        ('M1', 0, '', 0, 'S1'),
        ('S1', 0, 'product', 0, None),
        ('S1', 1, '', 0, 'H1'),
        ('H1', 0, '', 1, 'M1'),
    ]
    assert_flowsheet(nodes, edges, """
       -----[H1]
       |     ^
       v     |
feed->[M1]->[S1]->produ
""")

def test_recycle_through_many_units():
    """Long return paths of recycle loops wrap around the main line."""
    nodes = ['M1', 'S1', 'H1', 'H2', 'H3']
    edges = [
        (None, 0, 'feed', 0, 'M1'),
        ('M1', 0, '', 0, 'S1'),
        ('S1', 0, 'product', 0, None),
        ('S1', 1, '', 0, 'H1'),
        ('H1', 0, '', 0, 'H2'),
        ('H2', 0, '', 0, 'H3'),
        ('H3', 0, '', 1, 'M1'),
    ]
    assert_flowsheet(nodes, edges, """
      [H3]<-[H2]
       |     ^
       |     |
       |    [H1]
       |     ^
       v     |
feed->[M1]->[S1]->produ
""")

def test_overlapping_recycles():
    """Combines: units in series, simple recycle."""
    nodes = ['A1', 'B1', 'C1', 'D1']
    edges = [
        (None, 0, 'feed', 0, 'A1'),
        ('A1', 0, '', 0, 'B1'),
        ('B1', 0, '', 0, 'C1'),
        ('C1', 0, '', 0, 'D1'),
        ('C1', 1, '', 1, 'A1'),
        ('D1', 0, '', 1, 'B1'),
        ('D1', 1, 'product', 0, None),
    ]
    assert_flowsheet(nodes, edges, """
       -------------
       |           |
       |     ------|------
       |     |     |     |
       v     v     |     |
feed->[A1]->[B1]->[C1]->[D1]->produ
""")

def test_section_in_empty_space():
    """Combines: unconnected sections, simple recycle (sections fill empty space)."""
    nodes = ['M1', 'S1', 'HX']
    edges = [
        (None, 0, 'feed', 0, 'M1'),
        ('M1', 0, '', 0, 'S1'),
        ('S1', 0, 'product', 0, None),
        ('S1', 1, '', 1, 'M1'),
    ]
    assert_flowsheet(nodes, edges, """
[HX]   -------
       |     |
       v     |
feed->[M1]->[S1]->produ
""")

def test_manifold_into_unit_above():
    """Combines: split to two units, simple recycle (the bus goes straight up into units above)."""
    nodes = ['M1', 'C1', 'S1', 'B1']
    edges = [
        (None, 0, 'feed', 0, 'M1'),
        ('M1', 0, '', 0, 'C1'),
        ('C1', 0, '', 0, 'S1'),
        ('C1', 1, '', 0, 'B1'),
        ('B1', 0, 'waste', 0, None),
        ('S1', 0, 'product', 0, None),
        ('S1', 1, '', 1, 'M1'),
    ]
    assert_flowsheet(nodes, edges, """
       ----------------------
       |       [B1]->waste  |
       |        ^   ---------
       v        |   |
feed->[M1]->[C1]+->[S1]->produ
""")

# %% Combined features (each already tested)

def test_recycle_to_same_unit():
    """Combines: multiple products, simple recycle."""
    nodes = ['ST']
    edges = [
        (None, 0, 'feed', 0, 'ST'),
        ('ST', 0, 'vapor', 0, None),
        ('ST', 1, '', 1, 'ST'),
        ('ST', 2, 'liquid', 0, None),
    ]
    assert_flowsheet(nodes, edges, """
       ----
       |  |
       v  |
feed->[ST]+->vapor
          |
          +->liqui
""")

def test_recycle_through_unit_with_side_feed():
    """Combines: units in series, feed to downstream unit, recycle through unit."""
    nodes = ['S1', 'M1', 'O1', 'H1']
    edges = [
        (None, 0, 'feed', 0, 'S1'),
        ('S1', 0, '', 0, 'M1'),
        (None, 0, 'feed2', 1, 'M1'),
        ('M1', 0, '', 0, 'O1'),
        ('O1', 0, 'product', 0, None),
        ('O1', 1, '', 0, 'H1'),
        ('H1', 0, '', 1, 'S1'),
    ]
    assert_flowsheet(nodes, edges, """
       -----------[H1]
       |           ^
       v           |
feed->[S1]->[M1]->[O1]->produ
             ^
             |
         feed2
""")

def test_two_recycle_loops():
    """Combines: units in series, simple recycle, recycle through unit."""
    nodes = ['M1', 'M2', 'S1', 'S2', 'H1']
    edges = [
        (None, 0, 'feed', 0, 'M1'),
        ('M1', 0, '', 0, 'M2'),
        ('M2', 0, '', 0, 'S1'),
        ('S1', 0, '', 0, 'S2'),
        ('S2', 0, 'product', 0, None),
        ('S2', 1, '', 1, 'M2'),
        ('S1', 1, '', 0, 'H1'),
        ('H1', 0, '', 1, 'M1'),
    ]
    assert_flowsheet(nodes, edges, """
       -----------[H1]
       |           ^
       |     ------|------
       |     |     |     |
       v     v     |     |
feed->[M1]->[M2]->[S1]->[S2]->produ
""")

def test_long_recycle_loop_with_inner_recycle():
    """Combines: units in series, feed to downstream unit, simple recycle,
    recycle through many units."""
    nodes = ['S1', 'M1', 'O1', 'R1', 'X1', 'G1', 'A1', 'F1', 'H1', 'B1', 'C1']
    edges = [
        (None, 0, 'feed', 0, 'S1'),
        ('S1', 0, '', 0, 'M1'),
        (None, 0, 'feed2', 1, 'M1'),
        ('M1', 0, '', 0, 'O1'),
        ('O1', 0, '', 0, 'R1'),
        ('R1', 0, '', 0, 'X1'),
        ('X1', 0, '', 0, 'G1'),
        ('G1', 0, 'product', 0, None),
        ('R1', 1, '', 2, 'M1'),
        ('O1', 1, '', 0, 'A1'),
        ('A1', 0, '', 0, 'F1'),
        ('F1', 0, '', 0, 'H1'),
        ('H1', 0, '', 0, 'B1'),
        ('B1', 0, '', 0, 'C1'),
        ('C1', 0, '', 1, 'S1'),
    ]
    assert_flowsheet(nodes, edges, """
      [B1]<-[H1]<-[F1]
       |           ^
       v     ------|------
      [C1]   |    [A1]   |
       |     |     ^     |
       v     v     |     |
feed->[S1]->[M1]->[O1]->[R1]->[X1]->[G1]->produ
             ^
             |
         feed2
""")

def test_sugarcane_biorefinery():
    """Combines all features: the structure of the sugarcane-to-ethanol biorefinery
    in biorefineries.sugarcane (biosteam is not loaded)."""
    nodes = ['U1', 'U2', 'U3', 'U4', 'S1', 'M1', 'T1', 'H1', 'T2', 'P1', 'T3', 'T4', 'P2', 'M2', 'H2', 'T5', 'C1', 'C2', 'P3', 'S2', 'F1', 'P4', 'M3', 'H3', 'R1', 'T6', 'C3', 'S3', 'D1', 'M4', 'P5', 'H4', 'D2', 'P6', 'M5', 'D3', 'H5', 'U5', 'H6', 'T7', 'P7', 'T8', 'P8', 'M6', 'T9', 'P9', 'M7', 'U6', 'M8', 'HXN', 'CWP', 'CT', 'BT', 'PWC']
    edges = [
        (None, None, 'sugarcane', 0, 'U1'),
        ('U1', 0, 's1', 0, 'U2'),
        ('U2', 0, 's2', 0, 'U3'),
        ('U2', 1, 's3', None, None),
        ('U3', 0, 'stream_1', 0, 'U4'),
        ('M1', 0, 's4', 1, 'U4'),
        ('U4', 0, 's5', 0, 'U6'),
        ('U4', 1, 's6', 0, 'S1'),
        ('S1', 0, 's8', 0, 'T1'),
        ('S1', 1, 's7', 0, 'M1'),
        (None, None, 'imbibition_water', 1, 'M1'),
        ('T1', 0, 'stream', 0, 'H1'),
        ('H1', 0, 's9', 0, 'T2'),
        (None, None, 'H3PO4', 1, 'T2'),
        ('T2', 0, 's10', 0, 'P1'),
        ('P1', 0, 's11', 0, 'T4'),
        (None, None, 'lime', 0, 'T3'),
        ('T3', 0, 's12', 1, 'T4'),
        ('T4', 0, 's13', 0, 'P2'),
        ('P2', 0, 's14', 0, 'M2'),
        ('P3', 0, 's15', 1, 'M2'),
        ('M2', 0, 's16', 0, 'H2'),
        ('H2', 0, 's17', 0, 'T5'),
        (None, None, 'polymer', 1, 'T5'),
        ('T5', 0, 's18', 0, 'C1'),
        ('C1', 0, 'stream_2', 0, 'S2'),
        ('C1', 1, 's19', 0, 'C2'),
        (None, None, 'rvf_wash_water', 1, 'C2'),
        ('C2', 0, 'filter_cake', None, None),
        ('C2', 1, 's20', 0, 'P3'),
        ('S2', 0, 'juice', 0, 'F1'),
        ('S2', 1, 'fiber_fines', None, None),
        ('F1', 0, 's21', 0, 'P4'),
        ('F1', 1, 'evaporator_condensate', 1, 'M7'),
        ('P4', 0, 's22', 0, 'M3'),
        (None, None, 'dilution_water', 1, 'M3'),
        ('M3', 0, 's23', 0, 'H3'),
        ('H3', 0, 's24', 0, 'R1'),
        ('S3', 0, 's29', 1, 'R1'),
        ('R1', 0, 'CO2', 1, 'D1'),
        ('R1', 1, 's25', 0, 'T6'),
        ('T6', 0, 's26', 0, 'C3'),
        ('C3', 0, 's27', 0, 'S3'),
        ('C3', 1, 's28', 0, 'M4'),
        ('S3', 1, 'yeast', None, None),
        (None, None, 'stripping_water', 0, 'D1'),
        ('D1', 0, 'vent', None, None),
        ('D1', 1, 's30', 1, 'M4'),
        ('M4', 0, 'beer', 0, 'P5'),
        ('P5', 0, 's31', 0, 'H4'),
        ('P6', 0, 's33', 1, 'H4'),
        ('H4', 0, 's34', 0, 'D2'),
        ('H4', 1, 'vinasse', None, None),
        ('D2', 0, 'distilled_beer', 0, 'M5'),
        ('D2', 1, 's32', 0, 'P6'),
        ('U5', 0, 's39', 1, 'M5'),
        ('M5', 0, 's35', 0, 'D3'),
        ('D3', 0, 's36', 0, 'H5'),
        ('D3', 1, 's37', 0, 'P9'),
        ('H5', 0, 's38', 0, 'U5'),
        ('U5', 1, 's40', 0, 'H6'),
        ('H6', 0, 's41', 0, 'T7'),
        ('T7', 0, 's42', 0, 'P7'),
        ('P7', 0, 's43', 1, 'M6'),
        (None, None, 'denaturant', 0, 'T8'),
        ('T8', 0, 's44', 0, 'P8'),
        ('P8', 0, 's45', 0, 'M6'),
        ('M6', 0, 's46', 0, 'T9'),
        ('T9', 0, 'ethanol', None, None),
        ('P9', 0, 'recycle_process_water', 0, 'M7'),
        ('M7', 0, 's47', 2, 'PWC'),
        ('U6', 0, 'bagasse', 0, 'M8'),
        ('M8', 0, 's48', 0, 'BT'),
        (None, None, 'cooling_tower_makeup_water', 0, 'CT'),
        (None, None, 'cooling_tower_chemicals', 1, 'CT'),
        ('CT', 0, 'cooling_tower_blowdown', None, None),
        ('CT', 1, 'cooling_tower_evaporation', None, None),
        (None, None, 's49', 1, 'BT'),
        (None, None, 'boiler_makeup_water', 2, 'BT'),
        (None, None, 'natural_gas', 3, 'BT'),
        (None, None, 's50', 4, 'BT'),
        (None, None, 's51', 5, 'BT'),
        (None, None, '', 6, 'BT'),
        ('BT', 0, 'emissions', None, None),
        ('BT', 1, 'rejected_water_and_blowdown', None, None),
        ('BT', 2, 'ash_disposal', None, None),
        (None, None, 's52', 0, 'PWC'),
        (None, None, 'makeup_RO_water', 1, 'PWC'),
        (None, None, 's53', 3, 'PWC'),
        ('PWC', 0, 's54', None, None),
        ('PWC', 1, 's55', None, None),
        ('PWC', 2, 's56', None, None),
    ]
    assert_flowsheet(nodes, edges, """
                              imbib    [HX] [CW] cooli1-+>[CT]+->cooli3
                                  |                     |     |                                                                          -----------------
                                  v              cooli2-+     +->cooli4                                                                  |               |
                           ------[M1]                                                                                                    |              [S3]->yeast       ------[P6]   --------------------     denat->[T8]->[P8]-
                           |      ^                lime->[T3]--                                                                          |               ^                |      ^     |                  |                      |
                           v      |                           |                                                                          v               |                v      |     v                  |                      |
sugar->[U1]->[U2]+->[U3]->[U4]+->[S1]->[T1]->[H1]->[T2]->[P1]-+>[T4]->[P2]->[M2]->[H2]->[T5]->[C1]+->[S2]+---->[F1]+->[P4]->[M3]->[H3]->[R1]+->[T6]->[C3]+-+>[M4]->[P5]->[H4]+->[D2]->[M5]->[D3]+->[H5]->[U5]->[H6]->[T7]->[P7]--+>[M6]->[T9]->ethan
                 |            |                     ^                        ^           ^        |      |         |         ^              |              |                 |                  |
                 +->s3        |                     |                        |           |        -----  +->fiber  ----      |              --------       |                 +->vinas           -->[P9]-+>[M7]->[PW]+->s54
                              |                 HPO34                        |       polym            |               |  dilut                     |       |                                            |        ^  |
                              -->[U6]->[M8]->[BT]+->emiss                    |                        |               |                            v       |                                            |        |  +->s55
                                              ^  |                           |                        v               |                    strip->[D1]+-----                                            |    s52-+  |
                                              |  +->rejec                    |                rvfwa->[C2]->filte      |                               |                                                 |        |  +->s56
                                          s49-+  |                           |                        |               |                               +->vent                                           |  makeu-+
                                              |  +->ashdi                    |                        v               -----------------------------------------------------------------------------------        |
                                        boile-+                              ------------------------[P3]                                                                                                    s53-+
                                              |
                                        natur-+
                                              |
                                          s50-+
                                              |
                                          s51-+
                                              |
                                             -+
""")

if __name__ == '__main__':
    test_single_unit()
    test_stream_labels()
    test_units_in_series()
    test_multiple_feeds()
    test_two_feeds()
    test_multiple_products()
    test_unconnected_sections()
    test_feed_to_downstream_unit()
    test_split_to_two_units()
    test_merge_from_two_units()
    test_parallel_streams()
    test_manifolds_in_same_column()
    test_manifold_beside_label()
    test_bypass_over_manifold()
    test_products_on_separate_rows()
    test_feeds_on_separate_rows()
    test_simple_recycle()
    test_recycle_through_unit()
    test_recycle_through_many_units()
    test_overlapping_recycles()
    test_section_in_empty_space()
    test_manifold_into_unit_above()
    test_recycle_to_same_unit()
    test_recycle_through_unit_with_side_feed()
    test_two_recycle_loops()
    test_long_recycle_loop_with_inner_recycle()
    test_sugarcane_biorefinery()
