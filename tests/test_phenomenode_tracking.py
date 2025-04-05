# -*- coding: utf-8 -*-
"""
"""
import biosteam as bst
import thermosteam as tmo
import numpy as np
from numpy.testing import assert_allclose

def test_flash():
    import biosteam as bst
    import numpy as np
    from numpy.testing import assert_allclose
    bst.System.strict_phenomena_decomposition = True
    

def test_2_stage_flash():
    import biosteam as bst
    import numpy as np
    from numpy.testing import assert_allclose
    bst.System.strict_phenomena_decomposition = True

def test_column():
    import biosteam as bst
    import numpy as np
    from numpy.testing import assert_allclose
    bst.System.strict_phenomena_decomposition = True
    # bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    # feed = bst.Stream('feed', Ethanol=80, Water=100, T=80.215 + 273.15)
    # D1 = bst.MESHDistillation('D1', N_stages=5, ins=[feed], feed_stages=[2],
    #     outs=['vapor', 'liquid'],
    #     reflux=0.673, boilup=2.57,
    #     LHK=('Ethanol', 'Water') )
    # D1.set_flow_rates(D1.hot_start())
    # sys = bst.System('sys', path=[D1])
    # sys.track_variables = True
    pass
 
if __name__ == '__main__':
    pass