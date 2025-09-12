# -*- coding: utf-8 -*-
"""
This code tests a mock biorefinery system with integrated carbon capture.
It works as follows:
    
1. Simulate everything except carbon capture to get the emissions from the boiler.
2. Simulate the carbon capture. The emissions may change due to additional heat/power demand.
3. Repeat steps 1 and 2 until the emissions do not change.

"""
import biosteam as bst
from numpy.testing import assert_allclose

def test_integrated_facilities():
    
    class Biorefinery(bst.Unit):
        
        def _cost(self):
            self.add_power_utility(20)
            
    
    class CoHeatAndPower(bst.Facility):
        network_priority = 1
        _N_outs = 1
        @property
        def flue_gas(self):
            return self.outs[0]
        
        def _cost(self):
            power = sum([i.power_utility.consumption for i in self.other_units])
            self.add_power_utility(-power)
            self.flue_gas.imol['CO2'] = power / 10
            
    
    class CarbonCapture(bst.Unit):
        
        @property
        def flue_gas(self):
            return self.ins[0]
        @flue_gas.setter
        def flue_gas(self, flue_gas):
            self.ins[0] = flue_gas
        
        def _cost(self):
            power = self.flue_gas.F_mol * 0.1
            self.add_power_utility(power)
        
    bst.settings.set_thermo(['CO2'], cache=True)
    with bst.System('biorefinery') as biorefinery_sys: 
        B = Biorefinery()
        CHP = CoHeatAndPower()
        
    with bst.System('carbon_capture') as carbon_capture_sys: 
        CC = CarbonCapture()
        CC.flue_gas = CHP.flue_gas
        
    full_system = bst.System('full_system', path=[biorefinery_sys, carbon_capture_sys])
    full_system.set_tolerance(mol=1e-6, rmol=1e-6)
    full_system.simulate()
    
    assert_allclose(CC.flue_gas.mol[0], 2.020202)
    assert_allclose(CHP.power_utility.production, 20.20202)
    
if __name__ == '__main__':
    test_integrated_facilities()