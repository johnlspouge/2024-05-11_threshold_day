#!/usr/bin/env python
"""
Epidemic Decay
"""
from abc import ABC, abstractmethod

from jls_epidemic_infect import Infect

class Decay(ABC):
    @abstractmethod
    def __init__(self, 
            transit, # (self, next_compartment) indexes
            mean_recovery_time, # mean decay time
            dispersion=None, # dispersion of the decay time : the shape parameter
            delta=None): # grid for the computation
        Decay.assert_decay(transit, mean_recovery_time, dispersion, delta)
        self.transit = transit
        self.mean_recovery_time = mean_recovery_time
        self.dispersion = dispersion
        self.delta = delta
    def assert_decay(transit, mean_recovery_time, dispersion, delta):
        Infect.assert_transit(transit)
        assert isinstance(mean_recovery_time, (float,int)) and 0.0 < mean_recovery_time
        assert dispersion is None or isinstance(dispersion, (float,int)) and 0.0 < dispersion
        assert delta is None or isinstance(delta, (float,int)) and 0.0 < delta and delta/mean_recovery_time < 1.0
    @abstractmethod
    # Returns the time2transition for the next time of infection, or None if no infection is possible or tau_decay is earlier. 
    def time2transition_rvs(self, current_time, count):
        pass

def main(): 
    pass 

if __name__ == "__main__":
    main()
