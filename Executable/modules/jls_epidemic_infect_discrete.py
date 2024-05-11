#!/usr/bin/env python
"""
Epidemic Infect=Transit caused by Infectious compartment
"""
from math import isclose, prod, exp
import numpy as np
from scipy.stats import geom
from collections import Counter
from jls_epidemic_discrete_rvs import Zero_Truncated_General_Bernoulli, zero_truncated_binomial_rvs, zero_truncated_poisson_rvs
from jls_epidemic_util import add_count
from jls_epidemic_infect import Infect
# The infection is a discrete-time transit caused by c_infectious with a transmissibility beta to transit[0].
class Infect_Discrete(Infect):
    # Resolves infections with discrete delay to produce one or more earliest time2transition_rvs for the infection.
    #   The function is a class method, but the object appears to select the correct class in jls_epidemic_infect_mgr().
    #   tau_decay is the earliest event on the heap. It is None if the heap is empty.
    #   Returns None, if infection occurs after the earliest event on the heap or if no infection can occur.
    def time2transition_rvs(infectious2infect_binomials, current_time, tau_decay):
        #print('time2transition_rvs')
        i2itbs = infectious2infect_binomials
        assert isinstance(i2itbs, dict)
        for k,vs in i2itbs.items():
            assert isinstance(k, int) and 0 <= k < len(Infect.totals)
            for v in vs:
                assert isinstance(v, Infect_Discrete)
                assert k == v.c_infectious
        # Determines total intensity of infection for each infectious compartment.
        q_no_infections = [1.0] * len(Infect.totals)
        for k,vs in i2itbs.items():
            q_no_infections[k] = prod([v.q_no_infection() for v in vs])
        #print(q_no_infections)
        q_no_infection = prod(q_no_infections) # probability of no infection at an epoch
        assert 0.0 <= q_no_infection <= 1.0
        if q_no_infection == 1.0:
            return None
        cap_q = geom.rvs(1.0-q_no_infection) # geom k>=1, the earliest epoch of some infection among infectious compartments
        if tau_decay is not None and tau_decay < current_time+cap_q: # A decay occurs before the next infection.
            return None
        cb = Zero_Truncated_General_Bernoulli(q_no_infections)
        transition0 = {}
        Infect_Discrete._infectious_rv(transition0, cb, i2itbs)
        assert transition0
        return {current_time+cap_q:transition0} # transit must be a (hashable) tuple.
    # Adds the infection rvs of itbs corresponding to 1s in i_bernoulli_rvs to transition0.
    def _infectious_rv(transition0, cb, i2itbs):
        i_rvs = cb.rvs() # infectious with infection <-> 1
        for i,rv in enumerate(i_rvs): # i indexes the infectious compartments.
            if rv == 1: # The infectious compartment i has an infection.
                assert i2itbs.get(i) is not None
                itbs = i2itbs[i] 
                qs = [infectious.q_no_infection() for infectious in itbs]
                i_cb = Zero_Truncated_General_Bernoulli(qs)
                i_bernoulli_rvs = i_cb.rvs()
                Infect_Discrete._add_transition(transition0, i_bernoulli_rvs, itbs)
    # Adds the infection rvs of itbs corresponding to 1s in i_bernoulli_rvs to transition0.
    def _add_transition(transition0, i_rvs, itbs):
        for j,i_rv in enumerate(i_rvs): # j indexes the itbs, where i_rv=1 if itbs[j] infects.
            if i_rv == 1: # The infectious compartment i had an infection itbs[j].
                itb = itbs[j]
                infect_count = itb.infect_zero_truncated_count_rvs()
                transit = itb.transit
                transition = {transit: infect_count}
                add_count(transition0, transition)
# The infection is a binomial discretized-time transit caused by c_infectious with a transmissibility beta to transit[0].
class Infect_Binomial(Infect_Discrete):
    def __init__(self, c_infectious, beta, transit, delta, mean_recovery_time): 
        super().__init__(c_infectious, beta, transit, delta, mean_recovery_time)
        assert delta is not None 
        assert 0.0 < delta and delta/mean_recovery_time < 1.0 and delta*beta <= 1.0
        mean_recovery_time_delta = mean_recovery_time-delta
        self.b = mean_recovery_time_delta/mean_recovery_time*beta*delta
    def p_infection(self): # the probability that a single infectious infects a susceptible
        return self.b*Infect.totals[self.transit[0]]/Infect.total*Infect.nu
    def q_no_infection(self): # that probability that no infectious infects a susceptible
        return (1.0-self.p_infection())**Infect.totals[self.c_infectious]
    def infect_zero_truncated_count_rvs(self): # Returns a count of infections, given that an infection occurred.
        return zero_truncated_binomial_rvs(Infect.totals[self.c_infectious], self.p_infection(), size=None) 

    def test():
        # single infection
        Infect.nu = 1.0
        Infect.total = 10
        Infect.totals = [2,5,3]
        ite = Infect_Binomial(1, 2.0, (0,1), 0.1, 10.0)
        N = 1000
        np.random.seed(31415)
        counter = Counter([ite.infect_zero_truncated_count_rvs() for i in range(N)])
        assert counter == {1: 932, 2: 64, 3: 4}  
        # [ _conditioned_binomial_pmf(k, Infect.totals[ite.c_infectious], ite.p_infection()) for k in range(4) ]
        #     [0.0, 0.9208640640303245, 0.07593964376426666, 0.0031312056362608922 ... ]
        assert isclose(ite.b, 0.198) # 0.198 = (10-0.1)/10*2.0*1.0
        i2itbs = {1:[Infect_Binomial(1, 2.0, (0,1), 0.1, 10.0)]}
        current_time = 100
        tau_decay = 200
        np.random.seed(31415)
        assert Infect_Discrete.time2transition_rvs(i2itbs, current_time, tau_decay) == {106: {(0, 1): 1}}
        tau_decay = 105
        np.random.seed(31415)
        assert Infect_Discrete.time2transition_rvs(i2itbs, current_time, tau_decay) is None
        # 2 infectious, each with 1 infection
        Infect.nu = 1.0
        Infect.total = 10
        Infect.totals = [2,3,4,1]
        i2itbs = {2: [Infect_Binomial(2,2.0,(0,2), 0.1, 10.0)], 3: [Infect_Binomial(3,4.0,(1,3), 0.1, 10.0)]}
        assert isclose(i2itbs[2][0].p_infection(), 0.0396) # 0.99 * 2.0 * 0.1 * 2 / 10 * 1.0 
        assert isclose(i2itbs[3][0].p_infection(), 0.1188) # 0.99 * 4.0 * 0.1 * 3 / 10 * 1.0 
        assert isclose(i2itbs[2][0].q_no_infection(), 0.8507630225817857)
        assert isclose(i2itbs[3][0].q_no_infection(), 0.8812)
        current_time = 100
        tau_decay = 200
        np.random.seed(31415)
        assert Infect_Discrete.time2transition_rvs(i2itbs, current_time, tau_decay) == {104: {(1, 3): 1}}
        tau_decay = 101
        np.random.seed(31415)
        assert Infect_Discrete.time2transition_rvs(i2itbs, current_time, tau_decay) is None
        tau_decay = 200
        rvs = [Infect_Discrete.time2transition_rvs(i2itbs, current_time, tau_decay) for i in range(N)]
        mean_delay,counter_transits,counter_counts = Infect._stats(rvs, N, current_time, True)
        assert isclose(mean_delay, 3.975999999999999) # exact = 3.9950840570431074
        assert counter_transits == {(0, 2): 582, (1, 3): 489} # 71 infects beyond single: exponential gave {(0, 2): 557, (1, 3): 443}.
        obs = 582
        p2 = 1-i2itbs[2][0].q_no_infection()
        p3 = 1-i2itbs[3][0].q_no_infection()
        obs0 = p2/(p2+p3)*(582+489)
        chisq = (obs-obs0)**2*(1/obs0+1/(582+489-obs0))
        assert isclose(chisq, 0.7746619751795442)
        assert counter_counts == {1: 1036, 2: 33, 3: 2} # 71 infects beyond single: exponential gave {(0, 2): 557, (1, 3): 443}.
        # 1 infectious, with 2 infections
        Infect.nu = 1.0
        Infect.total = 10
        Infect.totals = [2,3,5]
        i2itbs = {2: [Infect_Binomial(2,2.0,(0,2), 0.1, 10.0), Infect_Binomial(2,4.0,(1,2), 0.1, 10.0)]}
        assert isclose(i2itbs[2][0].p_infection(), 0.0396)
        assert isclose(i2itbs[2][1].p_infection(), 0.1188)
        assert isclose(i2itbs[2][0].q_no_infection(), 0.817072806887547)
        assert isclose(i2itbs[2][1].q_no_infection(), 0.5313399155475583)
        current_time = 100
        tau_decay = 200
        np.random.seed(31415)
        assert Infect_Discrete.time2transition_rvs(i2itbs, current_time, tau_decay) == {102: {(0, 2): 1, (1, 2): 1}}
        tau_decay = 101
        np.random.seed(31415)
        assert Infect_Discrete.time2transition_rvs(i2itbs, current_time, tau_decay) is None
        tau_decay = 200
        rvs = [Infect_Discrete.time2transition_rvs(i2itbs, current_time, tau_decay) for i in range(N)]
        mean_delay,counter_transits,counter_counts = Infect._stats(rvs, N, current_time, True)
        assert isclose(mean_delay, 1.7309999999999945)
        assert counter_transits == {(1, 2): 802, (0, 2): 337}
        assert counter_counts == {1: 928, 2: 186, 3: 24, 4: 1}
        i2itbs = {2: [Infect_Binomial(2,2.0,(0,2), 0.001, 10.0), Infect_Binomial(2,4.0,(1,2), 0.001, 10.0)]}
        current_time = 0
        tau_decay = 20000
        rvs = [Infect_Discrete.time2transition_rvs(i2itbs, current_time, tau_decay) for i in range(N)]
        mean_delay,counter_transits,counter_counts = Infect._stats(rvs, N, current_time, True)
        assert isclose(mean_delay, 131.486) # cf exponential 0.125
        assert counter_transits == {(1, 2): 748, (0, 2): 254}
        assert counter_counts == {1: 1001, 2: 1}

# The infection is a Poisson discretized-time transit caused by c_infectious with a transmissibility beta to transit[0].
#   A deterministic decay permits an exact Poisson lattice infection probability, a discrete-time transit caused by c_infectious with a transmissibility beta to transit[0].
class Infect_Poisson(Infect_Discrete):
    def __init__(self, c_infectious, beta, transit, delta, mean_recovery_time): 
        super().__init__(c_infectious, beta, transit, delta, mean_recovery_time)
        assert delta is not None 
        assert isclose(mean_recovery_time, round(mean_recovery_time/delta)*delta)
        assert 0.0 < delta and delta/mean_recovery_time < 1.0 and delta*beta <= 1.0
    def single_mu(self): # the probability that a single infectious infects a susceptible
        return self.beta*self.delta*Infect.nu*Infect.totals[self.transit[0]]/Infect.total
    def mu(self): # the probability that a single infectious infects a susceptible
        return self.single_mu()*Infect.totals[self.c_infectious]
    def p_infection(self): # the probability that a single infectious infects a susceptible
        single_mu = self.beta*self.delta*Infect.nu*Infect.totals[self.transit[0]]/Infect.total
        return 1.0-exp(-single_mu)
    def q_no_infection(self): # that probability that no infectious infects a susceptible
        return exp(-self.mu())
    def infect_zero_truncated_count_rvs(self): # Returns a count of infections, given that an infection occurred.
        return zero_truncated_poisson_rvs(self.mu(), size=None) 

    def test():
        # single infection
        Infect.nu = 1.0
        Infect.total = 10
        Infect.totals = [2,5,3]
        ite = Infect_Poisson(1, 2.0, (0,1), 0.1, 10.0)
        N = 1000
        np.random.seed(31415)
        counter = Counter([ite.infect_zero_truncated_count_rvs() for i in range(N)])
        assert counter == {1: 919, 2: 75, 3: 5, 4: 1}  
        assert isclose(ite.q_no_infection(), 0.8187307530779818)
        assert isclose(exp(-ite.mu()), 0.8187307530779818)
        i2itbs = {1:[Infect_Poisson(1, 2.0, (0,1), 0.1, 10.0)]}
        current_time = 100
        tau_decay = 200
        np.random.seed(31415)
        assert Infect_Discrete.time2transition_rvs(i2itbs, current_time, tau_decay) == {106: {(0, 1): 1}}
        tau_decay = 105
        np.random.seed(31415)
        assert Infect_Discrete.time2transition_rvs(i2itbs, current_time, tau_decay) is None
        # 2 infectious, each with 1 infection
        Infect.nu = 1.0
        Infect.total = 10
        Infect.totals = [2,3,4,1]
        i2itbs = {2: [Infect_Poisson(2, 2.0, (0,2), 0.1, 10.0)], 3: [Infect_Poisson(3, 4.0, (1,3), 0.1, 10.0)]}
        assert isclose(i2itbs[2][0].p_infection(), 0.03921056084767682) 
        assert isclose(i2itbs[3][0].p_infection(), 0.11307956328284252) 
        assert isclose(i2itbs[2][0].q_no_infection(), 0.8521437889662113)
        assert isclose(i2itbs[3][0].q_no_infection(), 0.8869204367171575)
        current_time = 100
        tau_decay = None
        np.random.seed(31415)
        assert Infect_Discrete.time2transition_rvs(i2itbs, current_time, tau_decay) == {104: {(1, 3): 1}}
        tau_decay = 101
        np.random.seed(31415)
        assert Infect_Discrete.time2transition_rvs(i2itbs, current_time, tau_decay) is None
        tau_decay = None
        rvs = [Infect_Discrete.time2transition_rvs(i2itbs, current_time, tau_decay) for i in range(N)]
        mean_delay,counter_transits,counter_counts = Infect._stats(rvs, N, current_time, True)
        assert isclose(mean_delay, 4.046999999999997) # exact = 4.0947314726742805
        counter_transits = {(0, 2): 595, (1, 3): 465}
        obs = 595
        p2 = 1-i2itbs[2][0].q_no_infection()
        p3 = 1-i2itbs[3][0].q_no_infection()
        obs0 = p2/(p2+p3)*(595+465)
        chisq = (obs-obs0)**2/obs0
        assert isclose(chisq, 0.05289650583175232)
        assert counter_counts == {1: 982, 2: 75, 3: 3} # Cannot check, so for reproducibility only.
        # 1 infectious, with 2 infections
        Infect.nu = 1.0
        Infect.total = 10
        Infect.totals = [2,3,5]
        i2itbs = {2: [Infect_Poisson(2, 2.0, (0,2), 0.1, 10.0), Infect_Poisson(2, 4.0, (1,2), 0.1, 10.0)]}
        assert isclose(i2itbs[2][0].p_infection(), 0.03921056084767682)
        assert isclose(i2itbs[2][1].p_infection(), 0.11307956328284252)
        assert isclose(i2itbs[2][0].q_no_infection(), 0.8187307530779818)
        assert isclose(i2itbs[2][1].q_no_infection(), 0.5488116360940265)
        current_time = 100
        tau_decay = None
        np.random.seed(31415)
        assert Infect_Discrete.time2transition_rvs(i2itbs, current_time, tau_decay) == {102: {(0, 2): 1, (1, 2): 1}}
        tau_decay = 101
        np.random.seed(31415)
        assert Infect_Discrete.time2transition_rvs(i2itbs, current_time, tau_decay) is None
        tau_decay = None
        rvs = [Infect_Discrete.time2transition_rvs(i2itbs, current_time, tau_decay) for i in range(N)]
        mean_delay,counter_transits,counter_counts = Infect._stats(rvs, N, current_time, True)
        assert isclose(mean_delay, 1.7860000000000014) # exact = 1.815966220916094
        assert counter_transits == {(1, 2): 825, (0, 2): 331}
        obs = 331 # [2][0]
        p0 = 1-i2itbs[2][0].q_no_infection()
        p1 = 1-i2itbs[2][1].q_no_infection()
        obs0 = p0/(p0+p1)*(331+825)
        chisq = (obs-obs0)**2/obs0
        assert isclose(chisq, 0.0003133363837647393)
        assert counter_counts == {1: 909, 2: 202, 3: 38, 4: 6, 5: 1} # Cannot check, so for reproducibility only. p[0] Poisson > binomial.
        # 1 infectious, with 2 infections
        i2itbs = {2: [Infect_Poisson(2, 2.0, (0,2), 0.001, 10.0), Infect_Poisson(2, 4.0, (1,2), 0.001, 10.0)]}
        assert isclose(i2itbs[2][0].q_no_infection(), 0.9980019986673331)
        assert isclose(i2itbs[2][1].q_no_infection(), 0.9940179640539353)
        current_time = 0
        tau_decay = None
        rvs = [Infect_Discrete.time2transition_rvs(i2itbs, current_time, tau_decay) for i in range(N)]
        mean_delay,counter_transits,counter_counts = Infect._stats(rvs, N, current_time, True)
        #127.915 Counter({(1, 2): 737, (0, 2): 266}) Counter({1: 1000, 2: 3})
        assert isclose(mean_delay, 127.915) # exact = 125.50066666595599
        assert counter_transits == {(1, 2): 737, (0, 2): 266} # cf exponential {(1, 2): 724, (0, 2): 276}
        obs = 266 # [2][0]
        p0 = 1-i2itbs[2][0].q_no_infection()
        p1 = 1-i2itbs[2][1].q_no_infection()
        obs0 = p0/(p0+p1)*(737+266)
        chisq = (obs-obs0)**2/obs0
        assert isclose(chisq, 0.8809680243061018)
        assert counter_counts == {1: 1000, 2: 3} # Cannot check, so for reproducibility only. p[0] Poisson > binomial.

def main():
    Infect_Binomial.test()
    Infect_Poisson.test()
    
if __name__ == "__main__":
    main()
