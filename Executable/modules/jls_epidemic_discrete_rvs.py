#!/usr/bin/env python
"""
Epidemic successive simulation routines 
"""
from abc import ABC, abstractmethod

from math import isclose, prod, log, exp, ceil
import numpy as np
from scipy.stats import uniform, poisson, binom, bernoulli, chisquare
from collections import Counter
from bisect import bisect_left

class Successive_Simulation(ABC):
    @abstractmethod
    def rvs(self, n, size=None):
        pass  
# Returns a multinomial rv by successive simulation of binomials, while storing probabilities in a list.
#   p_fct:{0,1,2,...}->[0,1], a probability function whose values sum to 1.0
class Multinomial(Successive_Simulation):
    def __init__(self, p_fct):
        self.p_fct = p_fct
        self.ps = {0:p_fct(0)}
        self.pis = {0:1.0}
    def p(self, i):
        if i not in self.ps:
            self.ps[i] = self.p_fct(i)
        return self.ps[i]
    def pi(self, i):
        if i not in self.pis:
            self.pis[i] = self.pi(i-1)-self.p(i-1)
        return self.pis[i]
    def rvs(self, n, size=None):
        if size is not None:
            assert isinstance(size, int) and size > 0
            return [ self.rvs(n) for i in range(size) ]
        rvs = []
        i = 0
        while n > 0:
            rv = binom.rvs(n,self.p(i)/self.pi(i),size=1)[0]
            rvs.append(rv)
            n -= rv
            i += 1
        return rvs
# Simulates the rvs for independent Bernoulli trials P(failure)=qs[i], conditioned on at least one success.
class Zero_Truncated_General_Bernoulli:
    def __init__(self, qs):
        self.qs = qs # probabilities of independent events of no infection
        self.q_prod = prod([q for q in self.qs]) # probability no event occurs
    def rvs(self, size=None): # Simulates outcome, given that at least one event occurs.
        if size is not None:
            assert isinstance(size, int) and size > 0
            return [ self.rvs() for i in range(size) ]
        q_prod = self.q_prod
        assert 0.0 <= q_prod <= 1.0
        if q_prod == 1.0: # No event ever occurs.
            return None
        elif q_prod == 0.0: # Some event is certain.
            return [ bernoulli.rvs(1.0-q) for q in self.qs]
        # No event is certain, but some event can occur.
        has_occurred = False
        rvs = []
        for i,q in enumerate(self.qs):
            if has_occurred:
                occurring = bernoulli.rvs(1.0-q)
            else:
                occurring = bernoulli.rvs(probability((1.0-q)/(1.0-q_prod)))
                if occurring:
                    has_occurred = True
                else:
                    q_prod /= q
            rvs.append(occurring)
        return tuple(rvs)
# Simulates discrete rvs by the inverse cdf method with bisection.
#   Assumes the list cdfs is weakly increasing. Its elements are NOT pmf's.
def discrete_rvs(cdfs, size=None, values=None):
    assert values is None or len(values) == len(cdfs)
    assert isclose(cdfs[-1], 1.0)
    if size is not None:
        rvs = []
        for i in range(size):
            rvs.append(discrete_rvs(cdfs, None, values))
        return rvs
    i = bisect_left(cdfs, np.random.uniform())
    if values is None:
        return i
    return values[i]  
# Simulates 0-truncated poisson rvs by the first-occurrence method.
# https://en.wikipedia.org/wiki/Zero-truncated_Poisson_distribution
# https://stat.ethz.ch/pipermail/r-help/2005-May/070678.html
# https://web.archive.org/web/20180826164029/http://giocc.com/zero_truncated_poisson_sampling_algorithm.html
def zero_truncated_poisson_rvs(mu, size=None):
    if size is not None:
        return [zero_truncated_poisson_rvs(mu) for i in range(size)]
    u = 1.0-uniform.rvs()*(1.0-exp(-mu))
    t = -log(u)
    return 1+poisson.rvs(mu - t)
# Simulates 0-truncated binomial rvs by the first-occurrence method.
# https://www.mathworks.com/matlabcentral/answers/111057-how-can-i-create-a-zero-truncated-binomial-distribution-pf
def zero_truncated_binomial_rvs(n, p, size=None):
    if size is not None:
        return [zero_truncated_binomial_rvs(n, p) for i in range(size)]
    factor = 1.0-uniform.rvs()*(1.0-(1.0-p)**n)
    X = ceil(log(factor, 1.0-p))
    return 1+binom.rvs(n-X,p)
# Truncates to a probability.
def probability(x):
    if 1.0 < x:
        return 1.0
    elif x < 0.0:
        return 0.0
    return x
# 1 DEPRECATED 
# _zero_truncated_binomial_rvs focuses on pmf and cdf, rather than rvs.
#   Simulates discrete binomial rvs by the inverse cdf method.
def _zero_truncated_binomial_rvs(n, p, size=None):
    cdfs = _conditioned_binomial_cdfs(n, p)
    values = tuple(range(n+1))
    return discrete_rvs(cdfs, size, values)
# Simulates discrete binomial rvs by the first-occurrence method.
# Returns binomial cdf conditioned on positive result.
def _conditioned_binomial_cdfs(n, p):
    cdf = 0.0
    cdfs = []
    for k in range(n):
        cdfs.append(cdf)
        cdf += _conditioned_binomial_pmf(k+1, n, p)
    #assert isclose(cdf, 1.0)
    cdfs.append(1.0)
    return cdfs
# Returns binomial pmf conditioned on positive result.
#   The conditioned binomial pmf is numerically accurate even for p near 0.0.
def _conditioned_binomial_pmf(k, n, p):
    assert 0.0 <= p < 1.0
    geom_sum = sum([(1-p)**m for m in range(n)])
    assert isclose(p*geom_sum, 1.0-binom.pmf(0,n,p))
    if k == 0:
        return 0.0
    elif k == n:
        return p**(n-1)/geom_sum
    pmf = binom.pmf(k-1, n-1, p)
    pmf *= n/k
    pmf /= geom_sum
    return pmf
# 0 DEPRECATED
# Returns (chisq, chisq.pvalue, df) comparing dictionary {k:n[k]} from simulations with a discrete non-0 pmf in dictionary0={k:pmf(k)}.
def chisq_discrete_pmf(
        k2probability0, # ideal pmf with 0-elements deleted
        k2count, # rv simulation of the pmf (Counter output as dictionary)
        ddof=0): # deleted degrees of freedom (counts extra constraints)
    N = sum(k2count.values())
    f_exp = [N*v for k,v in k2probability0.items()]
    f_obs = [k2count.get(k,0) for k,v in k2probability0.items()]
    statistic, pvalue = chisquare(f_obs, f_exp, ddof=0, axis=0)
    return {'statistic':statistic, 'pvalue':pvalue, 'df':len(k2probability0)-1-ddof} # -1 matches constraint on sum N
# Tests Multinomial.
def test_chisq_discrete_pmf():
    n,p = 5,0.25
    k2probability0 = {k:binom.pmf(k,n,p) for k in range(0,n+1)}
    N = 1000
    np.random.seed(31415)
    counts = binom.rvs(n,p, size=N)
    k2count = dict(Counter(counts))
    out = chisq_discrete_pmf(k2probability0, k2count)
    outs = [v for v in out.values()]
    assert np.allclose(outs, (0.9515720164609083, 0.9663764897383553, 5))
# Tests Multinomial.
def test_Successive_Multinomial():
    np.random.seed(31415)
    def prob(i):
        p = 0.9
        return p*(1-p)**i
    multi = Multinomial(prob)
    assert isclose(multi.p(2), 0.9*0.1*0.1)
    assert isclose(multi.pi(2), 1.0 - 0.9 - 0.9*0.1)
    N = 10
    SIZE = 100
    rvs = multi.rvs(N,SIZE)
    length = max([len(rv) for rv in rvs])
    adds = [0] * length
    for i in range(length):
        for rv in rvs:
            if i < len(rv):
                adds[i] += rv[i]
        adds[i] /= SIZE*N
    adds0 = [0.902, 0.083, 0.012, 0.003]
    assert np.allclose(adds, adds0)
    chisq = SIZE*N * sum([(adds[i]-multi.p(i))**2/multi.p(i) for i in range(length)]) # df = 4
    assert isclose(chisq, 6.448888888888896)
# Tests the rvs for a Bernoulli trial conditioned on at least a single 1.
def test_Zero_Truncated_General_Bernoulli():
    np.random.seed(31415)
    # 2 trials
    cb = Zero_Truncated_General_Bernoulli((0.9,0.8))
    N = 1600
    cb2count = Counter(cb.rvs(N))
    cb2prob = {k:(v/N) for k,v in cb2count.items()}
    assert cb2prob == {(0, 1): 0.62375, (1, 0): 0.2975, (1, 1): 0.07875}
    cb2prob0 = {(0,1):9/14,(1,0):2/7,(1,1):1/14}
    chisq = N*sum([(cb2prob[k]-v)**2/v for k,v in cb2prob0.items()])
    assert isclose(chisq, 2.887222222222226)
    # 3 trials
    cb = Zero_Truncated_General_Bernoulli((0.5,0.2,0.1))
    assert isclose(cb.q_prod, 0.01)
    N = 1000
    rvs = cb.rvs(size=N)
    cb2count = Counter(rvs)
    cb2prob = {k:(v/N) for k,v in cb2count.items()}
    assert cb2prob == {(1, 1, 1): 0.364, (0, 1, 1): 0.372, (0, 0, 1): 0.097, (1, 0, 1): 0.084, (0, 1, 0): 0.035, (1, 1, 0): 0.044, (1, 0, 0): 0.004}
    cb0 = ((0,0,1),(0,1,0),(0,1,1),(1,0,0),(1,0,1),(1,1,0),(1,1,1))
    cb2prob0 = {}
    for k in cb0:
        v = 1.0
        for i,b in enumerate(k):
            f = (1.0-cb.qs[i])/cb.qs[i]
            v *= cb.qs[i]*(f**b)
        v /= 1.0-cb.q_prod
        cb2prob0[k] = v
    chisq = N*sum([(cb2prob[k]-v)**2/v for k,v in cb2prob0.items()])
    assert isclose(chisq, 5.853749999999998)

def test_discrete_rvs():
    np.random.seed(31415)
    N = 1000
    cdfs = [0.1,0.8,1.0]
    rvs = discrete_rvs(cdfs, size=N)
    i2count = Counter(rvs)
    assert i2count == {1: 695, 2: 202, 0: 103}
    np.random.seed(31415)
    rvs = discrete_rvs(cdfs, size=N, values=('a','b','c'))
    i2count = Counter(rvs)
    assert i2count == {'b': 695, 'c': 202, 'a': 103}
    chisq = (103-100)**2/100+(695-700)**2/700+(202-200)**2/200
    assert isclose(chisq, 0.1457142857142857)

def test_zero_truncated_binomial_pmf():
    n = 3
    p = 0.25
    cdfs = [ binom.cdf(k,n,p) for k in range(n+1) ]
    cdf_no_0 = [ (cdf-cdfs[0])/(1.0-cdfs[0]) for cdf in cdfs ]
    cpmfs = [ binom.pmf(k,n,p)/(1.0-binom.pmf(0,n,p)) for k in range(n+1) ]
    cpmfs[0] = 0.0
    assert np.allclose(cpmfs, [ _conditioned_binomial_pmf(k, n, p) for k in range(n+1) ])
    assert np.allclose(cdf_no_0, _conditioned_binomial_cdfs(n,p))

def test_zero_truncated_poisson_rvs():
    np.random.seed(31415)
    N = 1000
    mu = 3.0
    rvs = zero_truncated_poisson_rvs(mu, size=N)
    assert len(rvs) == N
    dd = Counter(rvs)
    assert 0 not in dd
    obs = [dd[k] for k in sorted(list(dd.keys()))]
    assert len(obs) == (sorted(list(dd.keys())))[-1] # no missing keys
    obs0 = [N*poisson.pmf(i+1, mu)/(1.0-exp(-mu)) for i in range(len(obs))]
    chisq = sum([(obs[i]-obs0[i])**2/obs0[i] for i in range(len(obs))])
    assert isclose(chisq, 9.67135807095378)

def test_zero_truncated_binomial_rvs():
    np.random.seed(31415)
    N = 1000
    n,p = 5,0.25
    obs0 = {(k+1):(N*_conditioned_binomial_pmf(k+1, n, p)) for k in range(n)}
    assert list(obs0.keys()) == list(range(1,6))
    obs0values = list(obs0.values())
    obs0values0 = [518.5659411011524, 345.7106274007683, 115.23687580025612, 19.20614596670936, 1.2804097311139564]
    assert np.allclose(obs0values, obs0values0)
    rvs = zero_truncated_binomial_rvs(n, p, size=N)
    counter = Counter(rvs)
    assert counter == {1: 539, 2: 346, 3: 95, 4: 19, 5: 1}
    chisq = sum([(counter[k]-obs0[k])**2/obs0[k] for k in sorted(list(counter.keys()))])
    assert isclose(chisq, 4.422887654320999)
    rvs = _zero_truncated_binomial_rvs(n, p, size=N)
    counter = Counter(rvs)
    assert counter == {1: 538, 2: 327, 3: 114, 4: 20, 5: 1}
    chisq = sum([(counter[k]-obs0[k])**2/obs0[k] for k in sorted(list(counter.keys()))])
    assert isclose(chisq, 1.848480246913588)

def test_probability():
    assert probability(2.0) == 1.0
    assert probability(0.4) == 0.4
    assert probability(-2.0) == 0.0

def main():
    test_chisq_discrete_pmf()
    test_Successive_Multinomial()
    test_Zero_Truncated_General_Bernoulli()
    test_discrete_rvs()
    test_zero_truncated_binomial_rvs()
    test_zero_truncated_binomial_pmf()
    test_zero_truncated_poisson_rvs()
    test_probability()
    
if __name__ == "__main__":
    main()
