#!/usr/bin/env python
"""
Epidemic Infect=Transit caused by Infectious compartment
"""
from abc import ABC, abstractmethod

from collections import Counter

# base class for Infect_Discrete and Infect_Continuous.
class Infect(ABC):
    delta = None # delta is the time-step for the simulation (None for continuous)
    nu = 1.0 # 0.0 < nu <= 1.0 is the effect of non-pharmaceutical interventions
    totals = [] # totals[i] counts individuals in compartment i and requires refreshing after every transit.
    total = 0 # initial total of all compartments, unchanging throughout computation
    index2infects = {}
    # Sets Infect class variables.
    def set_delta(delta):
        Infect.delta = delta
    def set_nu(nu):
        assert 0 < nu <= 1.0 # non-pharmaceutical interventions do not increase infection.
        Infect.nu = nu
    def set_totals(totals):
        Infect.totals = totals
        Infect.total = sum(Infect.totals) # The computation leaves Compartment.total unchanged.
        assert 0 < Infect.total
    def clear_index2infects():
        Infect.index2infects = {}
    def set_index2infect(infects, i):
        Infect.index2infects[i] = infects
    # Updates the Infect.totals, where transit==(int, int)
    def transit(transition):
        #print('transit')
        for k,v in transition.items():
            Infect.assert_transit(k)
            c_from, c_to = k
            v = int(v) # Removes the possibility of numpy type.
            if v > Infect.totals[c_from]:
                v = Infect.totals[c_from]
            Infect.totals[c_from] -= v
            Infect.totals[c_to] += v
    # Returns the intensity of infection.
    def beta(mean_recovery_time, r0):
        return r0/mean_recovery_time
    @abstractmethod
    def __init__(self, 
            c_infectious, # index of the infectious compartment
            beta, # infection intensity from infectious compartment on transit[0] compartment
            transit, # (susceptible, infected) compartment indexes
            delta=None, # grid for the computation
            mean_recovery_time=None): # decay time for infectious compartment
        Infect.assert_infect(c_infectious, beta, transit, delta, mean_recovery_time)
        self.c_infectious = c_infectious
        self.beta = beta
        self.transit = transit
        self.delta = delta
        self.mean_recovery_time = mean_recovery_time
    # Returns an arbitrary infect object to determine the subclass in time2transition_rvs() after calling Compartment.initialize(). 
    def infect_object():
        assert Infect.index2infects
        return list(Infect.index2infects.values())[0][0]
    # Checks the input for a transit.
    def assert_transit(transit):
        assert isinstance(transit, tuple) and len(transit) == 2
        for i in transit:
            assert isinstance(i, int) and 0 <= i < len(Infect.totals)
    # Checks the input for the Infect object.
    def assert_infect(c_infectious, beta, transit, delta, mean_recovery_time):
        assert isinstance(c_infectious, int) and 0 <= c_infectious
        assert isinstance(beta, float) and 0.0 < beta
        Infect.assert_transit(transit)
        assert delta is None or isinstance(delta, float)
        if delta is None:
            assert mean_recovery_time is None
        else:
            assert isinstance(mean_recovery_time, float) and 0.0 < mean_recovery_time
            assert 0.0 < delta
            assert beta*delta < 1.0
            assert delta/mean_recovery_time < 1.0
    # helper routine for debugging derived classes
    def _stats(rvs, N, current_time, is_no_assert_count=False):
        mean_delay = 0
        transits = []
        counts = []
        for rv in rvs:
            for delay,v in rv.items():
                mean_delay += delay
                for transit,count in v.items():
                    assert is_no_assert_count or count == 1
                    transits.append(transit)
                    counts.append(count)
        mean_delay /= N
        mean_delay -= current_time
        counter_transits = Counter(transits)
        counter_counts = Counter(counts)
        return mean_delay, counter_transits, counter_counts

def test_Infect():
    Infect.totals = [20,2,0]
    # The test tests changes to the totals.
    transition = {(0,1):2}
    Infect.transit(transition)
    assert Infect.totals == [18,4,0]
    # The next test challenges the restriction to weakly positive totals.
    transition = {(1,2):5}
    Infect.transit(transition)
    assert Infect.totals == [18,0,4]
    
def main():
    test_Infect() # The tests are in the derived classes: Infect_Exponential and Infect_Binomial
    
if __name__ == "__main__":
    main()
