# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

import flexsolve as flx
import numpy as np
from thermosteam.exceptions import InfeasibleRegion
import thermosteam as tmo

__all__ = ('ReactorSpecification',)

def evaluate_across_TRY(spec, system,
            titer, yield_, metrics, productivities):
    spec.load_specifications(titer=titer, yield_=yield_)
    system.simulate()
    return spec.evaluate_across_productivity(metrics, productivities)
    
evaluate_across_TRY = np.vectorize(
    evaluate_across_TRY, 
    excluded=['spec', 'system', 'metrics', 'productivities'],
    signature='(),(),(),(),(m),(p)->(m,p)'
)

class ReactorSpecification:
    """
    Create a ReactorSpecification object for setting reactor
    process specifications.
    
    Parameters
    ----------
    reactor : Unit
        Reactor unit operation. Must have an Reaction object attribute 
        called "reaction" where all products are produced.
    reaction_name : str
        Name of reaction attribute of the reactor.
    substrates : tuple[str]
        Name of main substrates.
    products : tuple[str]
        Names of main products.
    yield_ : float
        Yield in weight fraction of substrate converted to product 
        over theoretical yield. 
    titer : float
        g products / L effluent
    productivity : float
        g products / L effluent / hr
    production : float
        kg products / hr
    
    Notes
    -----
    Use the `load_specifications` method to load reactor specifications.
    Setting attributes (i.e. yield_, titer, productivity, production) sets 
    defaults for this method, but does not actually load any specifications.
    
    """
    
    __slots__ = ('reactor',
                 'substrates',
                 'products',
                 'yield_',
                 'titer',
                 'productivity',
                 'production',
                 'reaction_name')
    
    def __init__(self, reactor, reaction_name, substrates, products, yield_, titer, 
                 productivity, production):
        self.reactor = reactor #: [Unit] Reactor unit operation
        self.reaction_name = reaction_name #: reaction_name
        self.substrates = substrates #: tuple[str] Name of substrates
        self.products = products #: tuple[str] Names of main products
        self.yield_ = yield_ #: [float] Weight fraction of theoretical yield.
        self.titer = titer #: [float] g products / L effluent
        self.productivity = productivity  #: [float] g products / L effluent / hr
        self.production = float(production) #: [float] kg / hr
        self.reaction_name
      
    def load_specifications(self, yield_=None, titer=None, productivity=None,
                            production=None):
        """
        Load ferementation specifications.

        Parameters
        ----------
        yield_ : float, optional
            Yield in weight fraction of substrates converted to product 
            over theoretical yield. 
        titer : float, optional
            g products / L effluent
        productivity : float, optional
            g products / L effluent / hr
        production : float, optional
            kg products / hr

        """
        self.load_yield(yield_ or self.yield_)
        self.load_productivity(productivity or self.productivity)
        self.load_titer(titer or self.titer)
        self.load_production(production or self.production)

    def evaluate_across_productivity(self, metrics, productivities):
        """
        Evaluate metrics across productivities and return an array with the all
        metric results.
        
        Parameters
        ----------
        metrics : Iterable[Callable; M elements]
            Should return a number given no parameters.
        productivities : array_like[P elements]
            Productivities to evaluate.
        
        Returns
        -------
        results : array[M x P]
            All metric results.
        
        Notes
        -----
        Because setting productivity does not change any parameter associated
        to mass and energy balances, this method only simulates the reactor unit 
        operation at each productivity (as opposed to the whole system).
        
        """
        M = len(metrics)
        P = len(productivities)
        data = np.zeros([M, P])
        for i in range(P):
            self.load_productivity(productivities[i])
            self.reactor._setup()
            self.reactor._summary()
            data[:, i] = [j() for j in metrics]
        return data

    def evaluate_across_TRY(self, system, 
            titer, yield_, metrics, productivities):
        """
        Evaluate metrics at given titer and yield across a set of 
        productivities. Return an array with the all metric results.
            
        Parameters
        ----------
        titer : array_like[shape]
            Titer to evaluate.
        yield_ : array_like[shape]
            Yield to evaluate.
        metrics : Iterable[Callable; M elements]
            Should return a number given no parameters.
        productivities : array_like[P elements]
            Productivities to evaluate.
        
        Returns
        -------
        results : array[shape x M x P]
            All metric results at given titer/yield across productivities.
        
        Notes
        -----
        This method is vectorized along titer and yield. If, for example,
        the parameters had the following dimensions:
            
        titer [Y x T], yield [Y x T], metrics [M], productivities [P]
        
        This method would return an array with the following dimensions:
        
        results [Y x T x M x P]
        
        """
        return evaluate_across_TRY(self, system, 
                                   titer, yield_, 
                                   metrics, productivities)

    @property
    def feed(self):
        """[Stream] Reactor feed."""
        return self.reactor.ins[0]
    
    @property
    def vent(self):
        """[Stream] Reactor vent."""
        return self.reactor.outs[0]    
    
    @property
    def effluent(self):
        """[Stream] Reactor effluent."""
        return self.reactor.outs[1]
    
    @property
    def products_over_substrates(self):
        """[float] g products / g reactants."""
        reaction = self.reaction
        if reaction.basis != 'wt':
            reaction = reaction.copy('wt')
        substrates = self.substrates
        products = self.products
        production = - reaction.istoichiometry[products].sum()
        consumption = reaction.istoichiometry[substrates].sum()
        return production / consumption * reaction.X
    
    @property
    def reaction(self):
        """[Reaction] Stoichiometric reaction for converting substrates to products."""
        return getattr(self.reactor, self.reaction_name)
    
    def load_yield(self, yield_):
        """
        Load yield specification.
        
        Parameters
        ----------
        yield_ : float
            Yield in weight fraction of substrates converted to product 
            over theoretical yield.  
        
        Warnings
        --------
        Changing the yield affects the titer.
        
        """
        self.yield_ = self.reaction.X = yield_
    
    def load_titer(self, titer):
        """
        Load titer specification
        
        Parameters
        ----------
        titer : float
            Titer for fermentors in g products / L effluent.
        
        Notes
        -----
        Substrate concentration in bioreactor feed is adjusted to satisfy this 
        specification. 
        
        Warnings
        --------
        Changing the titer affects the productivity.
        
        """
        f = self._titer_objective_function
        feed = self.feed
        feed.imol[self.products] = 0.
        substrates = feed.imass[self.substrates].sum()
        self.titer = titer
        try:
            flx.aitken_secant(f, substrates, ytol=1e-5)
        except:
            flx.IQ_interpolation(f, 1e-12, 0.50 * feed.F_mass, ytol=1e-5, maxiter=100)
            
        self.reactor.tau = titer / self.productivity
    
    def load_productivity(self, productivity):
        """
        Load productivity specification.
        
        Parameters
        ----------
        productivity : float
            Productivity in g products / effluent L / hr.
        
        Notes
        -----
        Reaction time is adjusted to satisfy titer and productivity 
        specifications.
        
        """
        self.reactor.tau = self.titer / productivity
        self.productivity = productivity
    
    def load_production(self, production):
        """
        Load roduction specification.
        
        Parameters
        ----------
        production : float
            Production in kg product / hr.
        
        Notes
        -----
        The flow rate of feed to bioreactor is adjusted to achieve this
        specification.
        
        """
        F_mass_substrates = production / self.products_over_substrates
        feed = self.feed
        z_mass_substrates = feed.get_mass_composition(self.substrates).sum()
        feed.F_mass = F_mass_substrates / z_mass_substrates
        self.production = production
    
    def _calculate_titer(self):
        """Return titer in g products / effluent L."""
        reactor = self.reactor
        tmo.reaction.CHECK_FEASIBILITY = False
        reactor.run()
        tmo.reaction.CHECK_FEASIBILITY = True
        effluent = self.effluent
        F_mass_products = effluent.imass[self.products].sum()
        if F_mass_products: 
            return F_mass_products / effluent.F_vol
        else:
            return 0.
    
    def _titer_objective_function(self, substrates):
        """
        Return the titer of products given the ratio of substrates over feed 
        water.
        """
        if substrates <= 1e-16: raise InfeasibleRegion('substrate concentration')
        feed = self.feed
        mass_substrates = feed.imass[self.substrates]
        r_substrates = mass_substrates / mass_substrates.sum()
        mass_substrates[:] = substrates * r_substrates
        return self._calculate_titer() - self.titer
