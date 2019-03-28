# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 15:47:26 2018

@author: yoelr
"""

from .. import Unit
from ..exceptions import DesignError
from .mixer import Mixer
from .. import np

class Tank(Unit):
    """Abstract Tank class."""
    _N_ins = 1
    _N_outs = 1
    
    #: Abstract residence time attribute
    _tau = None
    _has_linked_streams = True

    @property
    def tau(self):
        """Residence time (hr)."""
        return self._tau
    @tau.setter
    def tau(self, tau):
        self._tau = tau

    def _run(self):
            self.outs[0].copylike(self.ins[0])


class StorageTank(Tank):
    """Create a storage tank with volume based on residence time [1].

    .. math::

        V &= \\tau Q \\\\
        C_{fob}^{2007} &= 250000 + 94.2 V (2*10^3 < V < 50*10^3 m^3) \\\\

    **References**

        [1] I.K. Kookos, Preliminary Chemical Process Synthesis and Design, Tziolas Publishing,
Thessalonika, Greece, 2008 (book in Greek).

    """
    #: [float] residence time (hr)
    _tau = 4*7*24

    def _design(self):
        """
        * 'Volume': Volume of vessel (m^3)
        * 'N_vessel': Number of vessels ()
        """
        Design = self.results['Design']
        V = self._tau*self._volnet_out
        if V < 50e3:
            Design['N_vessel'] = N = 1
        else:    
            Design['N_vessel'] = N = np.ceil(V/50e3)
        Design['Volume'] = V/N
        return Design

    @staticmethod
    def _calc_vessel_cost(V_vessel: 'm^3', CEPCI):
        if V_vessel < 2e3:
            cost = 65000 + 158.7*V_vessel * CEPCI/525.4
        elif V_vessel < 50e3:
            cost = 250000 + 94.2*V_vessel * CEPCI/525.4
        else:    
            raise DesignError(f"Volume is out of bounds for costing")
        return cost

    def _cost(self):
        """
        * 'Tank': (USD)
        """
        Design = self.results['Design']
        V_vessel, N_vessel = Design['Volume'], Design['N_vessel']
        CEPCI = self.CEPCI
        Cost = self.results['Cost']
        Cost['Tank'] = N_vessel * self._calc_vessel_cost(V_vessel/N_vessel, CEPCI)
        return Cost


class MixTank(Tank):
    """Create a mixing tank with volume based on residence time.

    .. math::

       V &= \\frac{\\tau* Q}{0.8}

       C_{fob}^{2007} &= 12080 * V^{0.525} (0.1 < V < 30 m^3)

    References

         [1] I.K. Kookos, Preliminary Chemical Process Synthesis and Design, Tziolas Publishing, Thessalonika, Greece, 2008 (book in Greek).

    """
    _tau = 1
    _N_ins = 2
    _has_power_utility = True
    _has_linked_streams = False
    _run = Mixer._run
    bounds = {'Volume': (0.1, 30)}
    
    #: Electricity rate (kW/m3)
    electricity_rate = 0.591 # About 10 hp/1000gal
    
    def _design(self):
        """
        * 'Volume': (m^3)
        """
        Design = self.results['Design']
        Design['Volume'] = self._tau * self._volnet_out / 0.8

    def _cost(self):
        """
        * 'Tank': (USD)
        """
        V = self.results['Design']['Volume']
        Cost = self.results['Cost']
        Cost['Tank'] = 12080 * V**0.525 * self.CEPCI/525.4
        self.power_utility(self.electricity_rate*V)


class PCT(MixTank):
    _N_ins = 1
    _N_outs = 2
    kwargs = {'T': 298.15, 'phase': 'l'}

    def _setup(self):
        """Set temperature and phase"""
        self.outs[0].T = self.kwargs['T']
        self.outs[0].phase = self.kwargs['phase']

    def _run(self): 
        out = self.outs[0]
        out.mol = self._mol_in
        out.P = self.ins[0].P
