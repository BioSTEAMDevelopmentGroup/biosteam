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

    # sum tries to return a HeatUtility object, but cannot due to heat utilities with different agents
    with pytest.raises(ValueError):
        sum(hus) 

    # check Q_actual == Q_expected using both pythonic sum and sum_by_agent method
    hus = bst.HeatUtility.sum_by_agent(hus)
    heat_utilities_by_agent = bst.HeatUtility.heat_utilities_by_agent(hus)
    hus_sum = tuple([sum(i) for i in heat_utilities_by_agent.values()])
    Q_actual_python_sum = {hu.agent.ID: hu.duty for hu in hus_sum}
    Q_actual = {hu.agent.ID: hu.duty for hu in hus}
    assert Q_actual_python_sum == Q_actual == Q_expected

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