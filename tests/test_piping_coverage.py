# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import pytest
import biosteam as bst
from numpy.testing import assert_allclose
from biosteam.utils import piping
from biosteam.utils.piping import (
    InletPort, OutletPort, StreamPorts, MissingStream, MockStream,
    as_stream, as_upstream, as_downstream, materialize_connections,
)

def test_ports():
    bst.settings.set_thermo(['Water'], cache=True)
    bst.F.set_flowsheet('piping_ports')
    M = bst.Mixer('M', ins=[bst.Stream('a', Water=1), bst.Stream('b', Water=2)], outs='c')
    inlet_port = InletPort.from_inlet(M.ins[1])
    assert inlet_port == InletPort(M, 1)
    assert str(inlet_port) == '1-M' and 'InletPort' in repr(inlet_port)
    outlet_port = OutletPort.from_outlet(M.outs[0])
    assert str(outlet_port) == 'M-0' and 'OutletPort' in repr(outlet_port)
    outlet_port.set_stream(bst.Stream('d'), 1)
    assert M.outs[0].ID == 'd'
    ports = StreamPorts.from_inlets(M.ins, sort=True)
    assert len(ports) == 2 and ports
    assert [i.ID for i in ports] == ['a', 'b']
    assert ports[0].ID == 'a'
    assert len(ports[0:1]) == 1
    assert len(ports + []) == 2 and len([] + ports) == 2
    assert repr(ports) == '[0-M, 1-M]'
    ports[0] = bst.Stream('e')
    ports[0:2] = [bst.Stream('f'), bst.Stream('g')]
    assert M.ins[0].ID == 'f'
    with pytest.raises(IndexError):
        ports[0:2] = [bst.Stream('h')]
    with pytest.raises(IndexError):
        ports['x'] = bst.Stream('h')

def test_stream_piping_notation():
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    bst.F.set_flowsheet('piping_notation')
    s = bst.Stream('s', Water=1, Ethanol=1)
    assert isinstance(s - 0, piping.InletPipe)
    assert isinstance(0 - s, piping.OutletPipe)
    assert isinstance(s ** 0, piping.InletPipe)
    assert isinstance(0 ** s, piping.OutletPipe)
    assert piping._stream_sub(s, 1.5) is NotImplemented
    assert piping._stream_rsub(s, 1.5) is NotImplemented
    assert piping._stream_pow(s, 1.5) is NotImplemented
    assert piping._stream_rpow(s, 1.5) is NotImplemented
    other = bst.Stream(None, Water=0.5)
    new = s - other
    assert_allclose(new.imol['Water'], 0.5)
    # Connections through junctions
    J = bst.Junction('J', bst.Stream('j_in'), bst.Stream('j_out'))
    P = bst.Pump('P', ins=J-0)
    connection = J.ins[0].get_connection(junction=False)
    assert connection is not None
    # Disjunctions
    piping.mark_disjunction(P.outs[0])
    piping.mark_disjunction(P.outs[0])
    piping.unmark_disjunction(P.outs[0])
    piping.unmark_disjunction(P.outs[0])

def test_missing_and_mock_streams():
    bst.settings.set_thermo(['Water'], cache=True)
    bst.F.set_flowsheet('piping_missing')
    assert isinstance(as_stream(None), MissingStream)
    assert as_stream('x').ID == 'x'
    s = bst.Stream('y')
    assert as_stream(s) is s
    assert as_upstream(None, 'sink').sink == 'sink'
    assert as_downstream(None, 'source').source == 'source'
    missing = MissingStream(None, None)
    assert repr(missing) == '<MissingStream>'
    missing.show()
    assert missing._get_tooltip_string('png', False) == ''
    assert missing._get_tooltip_string('html', False) == '(empty)'
    with pytest.raises(RuntimeError):
        missing.materialize_connection()
    assert missing.get_impact('GWP') == 0. and missing.get_CF('GWP') == 0.
    missing.reset_cache()
    assert missing.get_data() is None
    missing.set_data(None)
    assert missing.get_total_flow('kg/hr') == 0. and missing.get_flow('kg/hr') == 0.
    assert missing.isempty()
    # Materialize a missing stream that is an auxlet
    M = bst.Mixer('M', ins=[None, None])
    materialize_connections(M.ins)
    assert all(M.ins)
    aux = MissingStream(M, M)
    stream = aux.materialize_connection('aux')
    assert stream.source is M and stream.sink is M
    # Mock streams
    mock = MockStream('mock')
    assert not mock.isfeed() and not mock.isproduct()
    # Note: MockStream(1) and MockStream('') fail (no ticket_name)
    mock2 = MockStream(None)
    assert isinstance(mock2.ID, str)

def test_stream_sequences():
    bst.settings.set_thermo(['Water'], cache=True)
    bst.F.set_flowsheet('piping_sequences')
    M = bst.Mixer('M', ins='inlet')
    assert M.ins[0].ID == 'inlet'
    with pytest.raises(TypeError):
        M.ins[0] = 1
    s = bst.Stream('s1')
    M.ins.insert(0, s)
    assert M.ins.index(s) == 0
    M.ins.extend([bst.Stream('s2')])
    assert M.ins.size == 3
    assert M.ins.pop(0) is s
    M.ins.remove(M.ins[0])
    assert 'inlet' not in repr(M.ins)
    with pytest.raises(IndexError):
        M.ins['x'] = s
    # Setting beyond the end appends when size is not fixed
    M.ins[len(M.ins)] = bst.Stream('s3')
    assert M.ins[-1].ID == 's3'
    M.ins.empty()
    M.ins.clear()
    assert len(M.ins) == 0
    # Fixed size sequences
    P = bst.Pump('P', ins='p_in')
    for f in (P.ins.insert, ):
        with pytest.raises(RuntimeError): f(0, bst.Stream())
    with pytest.raises(RuntimeError): P.ins.append(bst.Stream())
    with pytest.raises(RuntimeError): P.ins.extend([bst.Stream()])
    with pytest.raises(IndexError): P.ins[3] = bst.Stream()
    removed = P.ins.pop(0)
    assert removed.ID == 'p_in' and not P.ins[0]
    P.ins.clear()
    assert len(P.ins) == 1
    with pytest.raises(ValueError):
        P.ins.index(bst.Stream('not_here'))
    assert P.ins.sink is P and P.outs.source is P

def test_static_flow_and_functors():
    bst.settings.set_thermo(['Water'], cache=True)
    bst.F.set_flowsheet('piping_static')
    P = bst.Pump('P', ins=bst.Stream('feed', Water=10))
    bst.utils.stream_link_options._run_static_flow(P)
    assert_allclose(P.outs[0].F_mol, 10)
    f = bst.utils.ExponentialFunctor(A=10.0, n=2.0)
    assert repr(f) == 'ExponentialFunctor(A=10.0, n=2.0)'
    assert_allclose(f(2.0), 40.0)

if __name__ == '__main__':
    test_ports()
    test_stream_piping_notation()
    test_missing_and_mock_streams()
    test_stream_sequences()
    test_static_flow_and_functors()
