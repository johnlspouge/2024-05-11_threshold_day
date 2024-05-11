#!/usr/bin/env python
"""
Epidemic Decay=Transit occurring spontaneously (must occur to an infected compartment)
    The Prob0Decay class exists for extensibility but is not presently used except with probs = [1.0] and a single decay.
"""
from math import isclose
import numpy as np
from scipy.stats import multinomial
from functools import reduce

from jls_epidemic_util import resolve_collision
from jls_epidemic_infect import Infect
from jls_epidemic_decay import Decay
from jls_epidemic_decay_continuous import Continuous_Deterministic
# Constructs an object containing probabilities of different decays for a single compartment.
#   The transit[0] and the delta for all decays must be the same.
class Prob0Decay:
    # Raises a ValueError if sum of decay probabilities is not close to 1.0.
    def error(probs, decays):
        and_fct = lambda x,y : x and y
        checks = [
            isinstance(probs, list) or isinstance(probs, tuple),
            isinstance(decays, list) or isinstance(decays, tuple),
            len(probs) == len(decays) ]
        is_alright = reduce(and_fct, checks) 
        if not is_alright:
            raise ValueError(f'The probs and decays must be matching lists.')
        for prob in probs:
            if not isinstance(prob, float):
                raise ValueError(f'prob = {prob} must be a float.')
            if not 0.0 <= prob <= 1.0:
                raise ValueError(f'0.0 <= prob = {prob} <= 1.0')
        if not isclose(sum(probs), 1.0):
            raise ValueError(f'Sum of probabilities = {sum(probs)} must be close to 1.0.')
        for decay in decays:
            if not isinstance(decay, Decay):
                raise ValueError(f'decay = {decay} must be a Decay.')
            if not decay.delta == decays[0].delta:
                raise ValueError(f'decay.delta = {decay.delta} == {decays[0].delta} = decays[0].delta')
            if not decay.transit[0] == decays[0].transit[0]:
                raise ValueError(f'decay.transit[0] = {decay.transit[0]} == {decays[0].transit[0]} = decays[0].transit[0]')
            
    def __init__(self, probs, decays):
        Prob0Decay.error(probs, decays)
        self.probs = probs
        self.decays = decays
    def time2transition_rvs(self, current_time, count):
        counts = multinomial.rvs(count, self.probs)
        time2transition = {}
        for i,cnt in enumerate(counts):
            decay = self.decays[i]
            time2transition_rvs = decay.time2transition_rvs(current_time, cnt)
            for k,v in time2transition_rvs.items():
                resolve_collision(k, time2transition, time2transition_rvs)
        return time2transition
        
def test_Prob0Decay():
    Infect.totals = [0,0,0] # Permits Infect.assert_transit()
    np.random.seed(31415)
    N = 1000
    # Tests decay_rvs_factory.
    decay2 = Continuous_Deterministic((2,0), 2.0)
    decay3 = Continuous_Deterministic((2,1), 3.0)
    decays = Prob0Decay((0.33,0.67),(decay2, decay3))
    rvs = decays.time2transition_rvs(0.0, N)
    for k,v in rvs.items():
        if k == 2.0:
            assert list(rvs[k].keys()) == [(2, 0)]
            assert list(rvs[k].values()) == [338]
        elif k == 3.0:
            assert list(rvs[k].keys()) == [(2, 1)]
            assert list(rvs[k].values()) == [662]
        else:
            assert False
    chisq = N*((0.338-0.33)**2/0.33 + (0.338-0.33)**2/0.67)
    assert isclose(chisq, 0.289461781999096)

def main(): 
    test_Prob0Decay()
    
if __name__ == "__main__":
    main()
