#!/usr/bin/env python
"""
Epidemic Infect=Transit caused by Infectious compartment
"""
from math import isclose
import numpy as np
from scipy.stats import expon
from collections import Counter
from jls_epidemic_infect import Infect
# The infection is a continuous-time transit caused by c_infectious with a transmissibility beta to transit[0].
class Infect_Continuous(Infect):    

    def __init__(self, c_infectious, beta, transit): 
        super().__init__(c_infectious, beta, transit)
    def _weight_of_infection(self): # Given c_infectious, weights for choosing the transit proportional to intensity.
        return self.beta*Infect.totals[self.transit[0]]
    def _factor_from_weight_to_intensity(self):
        return Infect.totals[self.c_infectious]/Infect.total*Infect.nu
    def intensity_of_infection(self): # Poisson intensity of infection
        return self._weight_of_infection()*self._factor_from_weight_to_intensity()
    def infect_count_rvs(self): # Returns a count of infections, given that an infection occurred.
        return 1
    # Resolves infections with continuous delay to produce a single earliest time2transition_rvs for the infection.
    #   The function is a class method, but the object appears to select the correct class in jls_epidemic_infect_mgr().
    #   tau_decay is the earliest event on the heap. It is None if the heap is empty.
    #   Returns None, if infection occurs after the earliest event on the heap or if no infection can occur.
    def time2transition_rvs(index2infect_countinuous, current_time, tau_decay): 
        assert tau_decay is None or current_time <= tau_decay
        i2ites = index2infect_countinuous
        assert isinstance(i2ites, dict)
        for k,vs in i2ites.items():
            assert isinstance(k, int) and 0 <= k < len(Infect.totals)
            for v in vs:
                assert isinstance(v, Infect_Continuous)
                assert k == v.c_infectious
        # Determines total intensity of infection for each infectious compartment.
        i2intensity = {}
        for k,vs in i2ites.items():
            i2intensity[k] = sum([e.intensity_of_infection() for e in vs])
        sum_intensity = sum(i2intensity.values()) # total intensity of infection
        if sum_intensity == 0.0:
            return None
        cap_q = expon.rvs()/sum_intensity # delay until next infection
        if tau_decay is not None and tau_decay < current_time+cap_q: # A decay occurs before the next infection.
            return None
        ps = [v/sum_intensity for k,v in i2intensity.items()]
        infectious = np.random.choice(list(i2intensity.keys()), p=ps) # Chooses an infectious compartment.
        ps = [e.intensity_of_infection()/i2intensity[infectious] for e in i2ites[infectious]]
        infect_transit = np.random.choice(i2ites[infectious], p=ps) # Chooses a transition for the infectious compartment.
        return {current_time+cap_q:{infect_transit.transit:1}} # transit must be a (hashable) tuple.

    def test():
        # single infection
        Infect.nu = 1.0
        Infect.total = 10
        Infect.totals = [2,5,3]
        ite = Infect_Continuous(1,2.0,(0,1))
        assert ite.intensity_of_infection() == 2.0*2*5/10*1.0 # 2.0
        assert ite.infect_count_rvs() == 1
        N = 1000
        rvs = [ite.infect_count_rvs() for i in range(N)]
        assert Counter(rvs) == {1:N}
        i2ites = {1:[Infect_Continuous(1,2.0,(0,1))]}
        current_time = 10.0
        tau_decay = 20.0
        np.random.seed(31415)
        assert Infect_Continuous.time2transition_rvs(i2ites, current_time, tau_decay) == {10.51897304219989: {(0, 1): 1}}
        Infect.totals = [0,5,5]
        assert Infect_Continuous.time2transition_rvs(i2ites, current_time, tau_decay) is None
        np.random.seed(31415)
        tau_decay = 10.5
        np.random.seed(31415)
        assert Infect_Continuous.time2transition_rvs(i2ites, current_time, tau_decay) is None
        # 2 infectious, each with 1 infection
        Infect.nu = 1.0
        Infect.total = 10
        Infect.totals = [2,3,4,1]
        i2ites = {2: [Infect_Continuous(2,2.0,(0,2))], 3: [Infect_Continuous(3,4.0,(1,3))]}
        ite = i2ites[2][0]
        assert isclose(i2ites[2][0].intensity_of_infection(), 2.0*2*4/10*1.0) # 1.6
        assert isclose(i2ites[3][0].intensity_of_infection(), 4.0*3*1/10*1.0) # 1.2
        current_time = 10.0
        tau_decay = 20.0
        np.random.seed(31415)
        assert Infect_Continuous.time2transition_rvs(i2ites, current_time, tau_decay) == {10.370695030142778: {(1, 3): 1}}
        tau_decay = 10.1
        np.random.seed(31415)
        assert Infect_Continuous.time2transition_rvs(i2ites, current_time, tau_decay) is None
        tau_decay = 20.0
        rvs = [Infect_Continuous.time2transition_rvs(i2ites, current_time, tau_decay) for i in range(N)]
        mean_delay,counter_transits,counter_counts = Infect._stats(rvs, N, current_time)
        assert np.allclose((mean_delay, 1/2.8), (0.36611130176902584, 0.35714285714285715))
        assert counter_transits == {(0, 2): 557, (1, 3): 443} # frequencies = (1.6/2.8,1.2/2.8) = (0.5714285714285715, 0.4285714285714286)
        assert counter_counts == {1: 1000} 
        chisq = (557-571.4)**2/571.4+(557-571.4)**2/428.6
        assert isclose(chisq, 0.8467058910574993)
        # 1 infectious, with 2 infections
        Infect.nu = 1.0
        Infect.total = 10
        Infect.totals = [2,3,5]
        i2ites = {2: [Infect_Continuous(2,2.0,(0,2)), Infect_Continuous(2,4.0,(1,2))]}
        ite = i2ites[2][0]
        assert isclose(i2ites[2][0].intensity_of_infection(), 2.0*2*5/10*1.0) # 2.0
        assert isclose(i2ites[2][1].intensity_of_infection(), 4.0*3*5/10*1.0) # 6.0
        current_time = 10.0
        tau_decay = 20.0
        np.random.seed(31415)
        assert Infect_Continuous.time2transition_rvs(i2ites, current_time, tau_decay) == {10.129743260549972: {(1, 2): 1}}
        tau_decay = 10.1
        np.random.seed(31415)
        assert Infect_Continuous.time2transition_rvs(i2ites, current_time, tau_decay) is None
        tau_decay = 20.0
        rvs = [Infect_Continuous.time2transition_rvs(i2ites, current_time, tau_decay) for i in range(N)]
        mean_delay,counter_transits,counter_counts = Infect._stats(rvs, N, current_time)
        assert np.allclose((mean_delay, 1/8.0), (0.12813895561915878, 0.125))
        assert counter_transits == {(1, 2): 724, (0, 2): 276} # (3/4,1/4)
        assert counter_counts == {1: 1000} 
        chisq = (724-750)**2/750+(724-750)**2/250
        assert isclose(chisq, 3.6053333333333333)
    
def main():
    Infect_Continuous.test()
    
if __name__ == "__main__":
    main()
