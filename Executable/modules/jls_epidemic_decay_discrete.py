#!/usr/bin/env python
"""
Epidemic Discrete Decay
"""
from math import isclose, floor
import numpy as np
from scipy.stats import binom, nbinom
from jls_epidemic_infect import Infect
from jls_epidemic_decay import Decay
from jls_epidemic_discrete_rvs import Multinomial

# Returns a (possibly trivial) mixture of two deltas, matching a deterministic mean recovery time. 
class Discrete_Deterministic(Decay):    
    def __init__(self, transit, mean_recovery_time, delta): 
        super().__init__(transit, mean_recovery_time, None, delta)
        x = self.mean_recovery_time/self.delta
        self.m = m = floor(x)
        self.p_m = m+1-x
    # Returns the time2transition for the time of decay. 
    def time2transition_rvs(self, current_time, count):
        if self.p_m == 1.0:
            time2transition = {current_time+self.m: {self.transit: count}}
        else:
            count_m = binom.rvs(count, self.p_m)
            time2transition = {current_time+self.m: {self.transit: count_m}, current_time+self.m+1: {self.transit: count-count_m} }
        return time2transition
# Returns the discrete RV M_delta, matching the mean gamma recovery time. 
class Discrete_Shifted_Neg_Binomial(Decay):    
    def __init__(self, transit, mean_recovery_time, dispersion, delta): 
        super().__init__(transit, mean_recovery_time, dispersion, delta)
        mean_recovery_time_delta = mean_recovery_time-delta
        factor = dispersion/mean_recovery_time_delta*delta
        self.a = factor/(1.0+factor)
        self.multinomial = Multinomial(self._pmf)
    def _pmf(self, i):
        return nbinom.pmf(i,self.dispersion,self.a)
    # Returns the time2transition for the time of decay. 
    def time2transition_rvs(self, current_time, count):
        rvs=self.multinomial.rvs(count)
        time2transition = {}
        for i,rv in enumerate(rvs):
            if rv > 0:
               time2transition[current_time+1+i] = {self.transit: rv} # Shifts the negative binomial rv by 1.
        return time2transition

# Tests probabilities of transits to c_tos = probability0delay0c_tos = [probability, delay, c_to].
def test_Discrete_Decay():
    Infect.totals = [0,0] # Permits Infect.assert_transit()
    N = 1000
   # discrete deterministic decays
    np.random.seed(31415)
    decay = Discrete_Deterministic(transit=(0,1), mean_recovery_time=2.0, delta=0.5)
    time2transition = decay.time2transition_rvs(current_time=10, count=N)
    assert time2transition == {14: {(0, 1): N}}
   # discrete deterministic decays
    np.random.seed(31415)
    decay = Discrete_Deterministic(transit=(0,1), mean_recovery_time=2.2, delta=0.5)
    time2transition = decay.time2transition_rvs(current_time=10, count=N)
    assert time2transition == {14: {(0,1): 592}, 15: {(0,1): 408}}
    # mean = (4 * 592 + 5 * 408) / 1000 = 4.408
    # mean0 = 4 + 2/5 = 4.4
    # discretized gamma decays = discrete shifted negative binomial decay
    decay = Discrete_Shifted_Neg_Binomial(transit=(0,1), mean_recovery_time=2.0, dispersion=10.0, delta=0.5)
    time2transition = decay.time2transition_rvs(current_time=10, count=N)
    # a = factor/(1.0+factor) where factor = 10.0 * (2.0 - 0.5) * 0.5 = 7.5, so a = 7.5/8.5 = 0.8823529411764706
    # mean = k * (1.0-a)/a / delta = 8/3 = 2.666666666666667
    time2count = {}
    for k,v in time2transition.items():
        v0 = list(v.values())
        time2count[k - 10] = v0[0]
    assert time2count == {1: 86, 2: 169, 3: 186, 4: 198, 5: 143, 6: 93, 7: 58, 8: 40, 9: 23, 10: 1, 11: 1, 12: 1, 13: 1}
    assert sum(time2count.values()) == N
    s = 0
    for k,v in time2count.items():
        s += k*v/1000
    assert isclose(s, 4.025999999999999) # mean = 4 = 2.0/0.5 = mean_recovery_time/delta
    time2count0 = {i+1:int(decay._pmf(i)*N) for i in range(len(time2count))}
    assert time2count0 == {1: 72, 2: 167, 3: 212, 4: 196, 5: 147, 6: 95, 7: 54, 8: 28, 9: 14, 10: 6, 11: 2, 12: 1, 13: 0}
    chisq = 0.0
    for i in range(8):
        chisq += (time2count[i+1]-time2count0[i+1])**2/time2count0[i+1]
    assert isclose(chisq, 11.545363966305231) # 7 to 8 df

def main():
    test_Discrete_Decay()

if __name__ == "__main__":
    main()
