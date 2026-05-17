import biosteam as bst
import pytest


def test_aerobic_polishing_filter_satisfies_oxygen_deficit_before_split():
    bst.main_flowsheet.clear()
    chemicals = bst.Chemicals(['CO2', 'AceticAcid', 'H2O', 'O2', 'N2'])
    bst.settings.set_thermo(
        bst.wastewater.high_rate.append_wwt_chemicals(chemicals)
    )
    influent = bst.Stream(
        'influent',
        AceticAcid=1,
        H2O=100,
        units='kg/hr',
    )
    recycle = bst.Stream('recycle')
    air = bst.Stream('air', phase='g')
    unit = bst.wastewater.high_rate.PolishingFilter(
        'PF',
        ins=(influent, recycle, air),
        outs=('biogas', 'effluent', 'waste', 'air_out'),
        filter_type='aerobic',
    )

    unit.simulate()
    assert unit.outs[1].imol['O2'] >= 0
    assert unit.outs[2].imol['O2'] >= 0
    
    # Overall and atomic mass balance
    assert abs(unit.mass_balance_error()) < 1e-6
    for atom, error in unit.atomic_balance_error().items():
        assert abs(error) < 1e-3, f'error in {atom} atomic balance is over {error:.1%}'


def test_polishing_filter_rejects_invalid_filter_type():
    bst.main_flowsheet.clear()
    chemicals = bst.Chemicals(['CO2', 'AceticAcid', 'H2O', 'O2', 'N2'])
    bst.settings.set_thermo(
        bst.wastewater.high_rate.append_wwt_chemicals(chemicals)
    )

    with pytest.raises(ValueError, match="filter_type"):
        bst.wastewater.high_rate.PolishingFilter(
            'PF',
            ins=('influent', 'recycle', 'air'),
            outs=('biogas', 'effluent', 'waste', 'air_out'),
            filter_type='aerobc',
        )

if __name__ == '__main__':
    test_aerobic_polishing_filter_satisfies_oxygen_deficit_before_split()
    test_polishing_filter_rejects_invalid_filter_type()