#!/usr/bin/env python
"""
Single-type Galton-Watson processes
"""
from abc import ABC, abstractmethod

from math import exp, log, lgamma, isclose
import numpy as np

class Branching_Process(ABC):
    @abstractmethod
    def __init__(self):
        pass
    @abstractmethod
    def name(self):
        pass
    @abstractmethod
    def probability_generating_function(self, s, n=0):  # n is the derivative #
        pass
    def expected_number_of_offspring(self):
        if not hasattr(self,'_mu'):
            self._mu = self.probability_generating_function(1.0, n=1)
        return self._mu
    # Returns the extinction probability of the Branching_Process.
    def q(self):
        if not hasattr(self,'_q'):
            q = 0.0
            while not isclose(q, self.probability_generating_function(q)):
                q = self.probability_generating_function(q)
            self._q = q
            assert 0.0 <= q <= 1.0
        return self._q
    # Returns the slope of the pgf at the extinction probability q.
    def gamma(self):
        if not hasattr(self,'_gamma'):
            self._gamma = self.probability_generating_function(self.q(), n=1)
        return self._gamma
    # Returns gamma/mu.
    def rho(self):
        if not hasattr(self,'_rho'):
            self._rho = self.gamma()/self.expected_number_of_offspring()
        return self._rho    
    #
    # The following routines are definitely restricted to single-type Galton-Watson processes.
    #
    # Returns the values of the pgfs for the Harris_Sevastyanov transformation of a supercritical Branching_Process.
    def harris_sevastyanov_probability_generating_functions(self, a, b):
        _q = self.q()
        pgf = self.probability_generating_function
        if _q == 1.0:
            return 0.0, pgf(b)
        fa = (pgf((1.0-_q)*a+_q*b)-_q)/(1.0-_q)
        fb = pgf(_q*b)/_q
        return [fa, fb]
    # Returns n-th power of the numpy expectation matrix M for Harris_Sevastyanov transformation of a supercritical Branching_Process.
    def harris_sevastyanov_expectation_matrix(self, n=1): # The matrix m is raised to the power n.
        _q = self.q()
        assert 0.0 < _q < 1.0
        _mu = self.expected_number_of_offspring()
        assert 1.0 < _mu
        _rho = self.rho()
        assert 0.0 < _rho < _mu
        if n == 1:
            return _mu*np.array([[1.0, _q/(1.0-_q)],[0,_rho**n]])
        return _mu**n*np.array([[1.0, (1.0-_rho**n)/(1.0-_rho)*_q/(1.0-_q)],[0,_rho**n]])
    # Returns the expectation matrix M for Harris_Sevastyanov transformation of a supercritical Branching_Process.
    def expectation_z_n_from_immortal(self, n): # The total n-th generation from an immortal.
        _q = self.q()
        assert 0.0 < _q < 1.0
        _mu = self.expected_number_of_offspring()
        assert 1.0 < _mu
        _rho = self.rho()
        assert 0.0 < _rho < _mu
        return (1.0+(1.0-_rho**n)/(1.0-_rho)*_q/(1.0-_q))*_mu**n
# Returns a Galton-Watson process with Negative_Binomial offspring distribution.        
class Negative_Binomial(Branching_Process):
    # Returns the *args for __init__ from epidemic parameters
    def args(r0, dispersion): 
        assert 0.0 < r0
        assert dispersion is not None and 0.0 < dispersion
        (k, p) = (dispersion, dispersion/(r0+dispersion))
        return (k, p)
    
    def __init__(self, *args): 
        (self.k, self.p) = args
    def name(self):
        return 'Negative_Binomial'
    def probability_generating_function(self, s, n=0): # n is the derivative #
        assert isinstance(n,int) and n >= 0
        (k, p) = (self.k, self.p)
        pgf = (p/(1.0-(1.0-p)*s))**k
        if n == 0:
            return pgf
        factor = exp(lgamma(k+n)-lgamma(k)+n*(log(1.0-p)-log(1.0-(1.0-p)*s)))
        return factor*pgf
# Returns a Galton-Watson process with Poisson offspring distribution.        
class Poisson(Branching_Process):
    # Returns the argument for __init__ from epidemic parameters
    def arg(r0):
        assert 0.0 < r0
        mu = r0
        return mu
    
    def __init__(self, mu): 
        self.mu = mu
    def name(self):
        return 'Poisson'
    def probability_generating_function(self, s, n=0): # n is the derivative #
        mu = self.mu
        pgf = exp(mu*(s-1.0))
        if n == 0:
            return pgf
        factor = mu**n
        return factor*pgf
# The interface uses epidemiological notation (r0, dispersion).        
def Branching_Process_Factory(r0, dispersion=None): # epidemic parameters
    if dispersion is None:
        bp = Poisson(r0)
    else:
        bp = Negative_Binomial(*Negative_Binomial.args(r0, dispersion))
    return bp
# Permits arbitrary mixtures of single-type Galton-Watson processes.        
class Branching_Process_Mixture(Branching_Process):
    def __init__(self, gwp0probs):
        assert isinstance(gwp0probs, tuple)
        for t in gwp0probs:
            assert len(t) == 2
            gwp,p = t
            assert isinstance(p, float) and p > 0
        assert isclose(sum([t[1] for t in gwp0probs]), 1.0)
        # Adds probabilities for repeated Galton-Watson processes.
        gwp2prob = {}
        for t in gwp0probs:
            gwp,p = t
            p0 = gwp2prob.get(gwp,0.0)
            p0 += p
            gwp2prob[gwp] = p0
        self.gwp2prob = gwp2prob
    def name(self):
        return 'Galton_Watson_Process_Mixture'
    def probability_generating_function(self, s, n=0):  # n is the derivative #
        assert isclose(sum(self.gwp2prob.values()), 1.0)
        pgf = 0.0
        for k,v in self.gwp2prob.items():
            pgf += v*k.probability_generating_function(s, n)
        return pgf

def test_Branching_Process():
    # test of supercritical Poisson branching process
    bp = Branching_Process_Factory(r0=2.0) # Poisson
    assert bp.name() == 'Poisson'
    assert isclose(bp.probability_generating_function(0.5,n=0), exp(-1.0))
    assert isclose(bp.probability_generating_function(0.5,n=1), exp(-1.0)*2.0)
    assert isclose(bp.probability_generating_function(0.5,n=2), exp(-1.0)*2.0*2.0)
    assert isclose(bp.expected_number_of_offspring(), 2.0)
    assert isclose(bp.q(), 0.20318786982853076)
    assert isclose(bp.gamma(), bp.probability_generating_function(bp.q(),n=1))
    assert isclose(bp.rho(), bp.gamma()/bp.expected_number_of_offspring())
    fa0,fb0 = ((bp.probability_generating_function(0.25+0.5*bp.q())-bp.q())/(1.0-bp.q()),
        bp.probability_generating_function(0.75*bp.q())/bp.q())
    fa,fb = bp.harris_sevastyanov_probability_generating_functions(0.25,0.75)
    assert np.allclose((fa,fb), (fa0,fb0))
    r = bp.expected_number_of_offspring()
    m = bp.harris_sevastyanov_expectation_matrix(n=1)
    q = bp.q()
    rho = bp.rho()
    assert np.allclose(m, [[r, r*q/(1.0-q)],[0,r*rho]])
    m2 = np.matmul(m, m)
    assert np.allclose(bp.harris_sevastyanov_expectation_matrix(n=2),m2)
    assert isclose(bp.expectation_z_n_from_immortal(2), m2[0][0]+m2[0][1])
    # test of supercritical Negative_Binomial branching process
    bp = Branching_Process_Factory(r0=2.0, dispersion=1.0) # Negative Binomial
    assert bp.name() == 'Negative_Binomial'
    assert isclose(bp.probability_generating_function(0.5,n=0), 0.5)
    assert isclose(bp.probability_generating_function(0.5,n=1), 0.5*1.0)
    assert isclose(bp.probability_generating_function(0.5,n=2), 0.5*2.0)
    # additional tests of extinction probability
    assert isclose(bp.expected_number_of_offspring(), 2.0)
    q = bp.q()
    assert isclose(q, 0.4999999990686774)
    bp = Branching_Process_Factory(r0=3.0) # Poisson
    q = bp.q()
    assert isclose(q, 0.0595202092363786)
    bp = Branching_Process_Factory(r0=3.0, dispersion=0.5) # Negative Binomial
    q = bp.q()
    assert isclose(q, 0.4999999996968479)
    # Tests Branching_Process_Factory.
    bp = Branching_Process_Factory(r0=3.0) # Poisson
    q = bp.q()
    assert isclose(q, 0.0595202092363786)
    bp = Branching_Process_Factory(r0=3.0, dispersion=0.5) # Negative Binomial
    q = bp.q()
    assert isclose(q, 0.4999999996968479)
    # test of another supercritical Negative_Binomial branching process
    bp = Branching_Process_Factory(r0=1.5, dispersion=1.0) # Negative Binomial
    assert bp.name() == 'Negative_Binomial'
    assert isclose(bp.probability_generating_function(0.5,n=0), 0.5714285714285715)
    assert isclose(bp.probability_generating_function(0.5,n=1), 0.48979591836734704)
    assert isclose(bp.probability_generating_function(0.5,n=2), 0.8396501457725947)
    q = bp.q()
    assert isclose(q, 0.6666666649022972)
    # Tests Branching_Process_Mixture for a mixture of identical offspring distributions.
    bp0 = Branching_Process_Factory(r0=2.0, dispersion=1.0) # Negative Binomial
    q = bp0.q()
    for p_tilde in list(range(1,10)):
        p = 0.1*p_tilde
        gwp0prob = ((bp0,p),(bp0,1.0-p),)
        bp = Branching_Process_Mixture(gwp0prob)
    assert isclose(bp.expected_number_of_offspring(), bp0.expected_number_of_offspring())
    assert isclose(bp.q(), q)
    # Mixes 0.25*Negative_Binomial(r0=2.0, dispersion=1.0)+0.75*Negative_Binomial(r0=1.5, dispersion=0.4)
    gwp0prob = ((Branching_Process_Factory(2.0,1.0),0.25),
                (Branching_Process_Factory(1.5,0.4),0.75),)
    # Tests Branching_Process_Mixture for reproducibility.
    bp = Branching_Process_Mixture(gwp0prob)
    assert bp.name() == 'Galton_Watson_Process_Mixture'
    assert isclose(bp.expected_number_of_offspring(), 0.25*2.0+0.75*0.75*2.0)
    assert isclose(bp.probability_generating_function(0.5,n=0), 0.6165934908357035)
    assert isclose(bp.probability_generating_function(0.5,n=1), 0.3814835604360193)
    assert isclose(bp.probability_generating_function(0.5,n=2), 0.7183612842744699)
    q = bp.q()
    assert isclose(q, 0.7281434068918151)

def main(): 
    test_Branching_Process()
      
if __name__ == "__main__":
    main()
