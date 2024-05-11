#!/usr/bin/env python
"""
Epidemic Continuous Decay
"""
from math import isclose
import numpy as np
from scipy.stats import gamma
from jls_epidemic_decay import Decay
from jls_epidemic_infect import Infect

class Continuous_Deterministic(Decay):    
    def __init__(self, transit, mean_recovery_time): 
        super().__init__(transit, mean_recovery_time)
    # Returns the time2transition for the time of decay. 
    def time2transition_rvs(self, current_time, count):
        time2transition = {}
        time2transition = {current_time+self.mean_recovery_time: {self.transit:count}}
        return time2transition

class Continuous_Gamma(Decay):    
    def __init__(self, transit, mean_recovery_time, dispersion): 
        super().__init__(transit, mean_recovery_time, dispersion)
    # Returns the time2transition for the time of decay. 
    def time2transition_rvs(self, current_time, count):
        rvs = gamma.rvs(self.dispersion,scale=1.0/(self.dispersion/self.mean_recovery_time),size=count) # parametrization: a=shape, scale=1.0/rate
        time2transition = {current_time+k: {self.transit:1} for k in rvs}
        return time2transition
    
def test_Continuous_Decay():
    Infect.totals = [0,0] # Permits Infect.assert_transit()
    N = 1000
   # continuous deterministic decays
    np.random.seed(31415)
    decay = Continuous_Deterministic((0,1), 2.0)
    time2transition = decay.time2transition_rvs(current_time=10.0, count=N)
    assert time2transition == {12.0: {(0, 1): N}}
   # continuous gamma decays (k, lambda) mean = k/lambda and var = k/lambda**2
    np.random.seed(31415)
    decay = Continuous_Gamma((0,1), 2.0, dispersion=10.0) # mean = 2.0 = k/lambda and var = k/lambda**2 = mean**2/k = 4.0/10.0
    time2transition = decay.time2transition_rvs(current_time=10.0, count=5)
    assert time2transition == {12.910292992655329: {(0, 1): 1}, 12.727772511304437: {(0, 1): 1}, 13.806465876564971: {(0, 1): 1}, 11.776378119458881: {(0, 1): 1}, 12.110276516211126: {(0, 1): 1}}
    np.random.seed(31415)
    time2transition = decay.time2transition_rvs(current_time=10.0, count=N)
    assert isclose(np.average(list(time2transition.keys())) - 10.0, 1.994022803994522) # mean = 2.0 = k/lambda and var = k/lambda**2 = mean**2/k = 4.0/10.0
    assert isclose(np.var(list(time2transition.keys())), 0.4000517218961863) # mean = 2.0 = k/lambda and var = k/lambda**2 = mean**2/k = 4.0/10.0

def main():
    test_Continuous_Decay()

if __name__ == "__main__":
    main()
