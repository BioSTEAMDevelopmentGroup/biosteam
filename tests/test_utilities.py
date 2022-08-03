import pytest
import biosteam as bst
from numpy import allclose

def test_heat_util_sum():

    # get all agents
    cool_agents = bst.HeatUtility.cooling_agents
    heat_agents = bst.HeatUtility.heating_agents

    # create two utilities per agent
    Q = 1
    T = 298
    hus = []
    Q_expected = {}
    for agent in cool_agents + heat_agents:
        hu = bst.HeatUtility()
        hu(unit_duty=Q, T_in=T, agent=agent)
        hus += [hu, hu.copy()]
        Q_expected[agent.ID] = 2*Q/agent.heat_transfer_efficiency

    # sum by agent
    hus = sum(hus)

    # check Q_actual == Q_expected
    Q_actual = {hu.agent.ID:hu.duty for hu in hus}
    assert allclose(
        a = [Q_expected[a] for a in Q_actual.keys()],
        b = list(Q_actual.values()),
    )
    pass

def test_power_util_sum():

    # create utilities
    production = [1,0,3,0,5]
    consumption = [0,2,0,4,0]
    pus = [
        bst.PowerUtility(production=p, consumption=c)
        for p, c in zip(production, consumption)
    ]

    # sum utilities
    pus = sum(pus)

    # check result
    expected_production = sum(production)
    expected_consumption = sum(consumption)
    assert allclose(
        a = [pus.production, pus.consumption],
        b = [expected_production, expected_consumption],
    )
    pass

if __name__ == '__main__':
    test_heat_util_sum()
    test_power_util_sum()