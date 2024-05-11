#!/usr/bin/env python
"""
Epidemic Infect=Transit caused by Infectious compartment
"""
from jls_epidemic_infect import Infect
from jls_epidemic_infect_continuous import Infect_Continuous
from jls_epidemic_infect_discrete import Infect_Binomial, Infect_Poisson

class Infect_Mgr:
    # The argument delta determines whether the Infect is discrete or continuous.
    def factory(c_infectious, beta, transit, delta=None, mean_recovery_time=None, is_poisson=False):
        if delta is None:
            return Infect_Continuous(c_infectious, beta, transit)
        elif is_poisson:
            return Infect_Poisson(c_infectious, beta, transit, delta, mean_recovery_time)
        else:
            return Infect_Binomial(c_infectious, beta, transit, delta, mean_recovery_time)
    # Returns the time2transition for the next time of infection, or None if no infection is possible or tau_decay is earlier. 
    #   The function is a class method, but the object appears to select the correct class in jls_epidemic_infect_mgr().
    # The argument index2infects determines which subclass of Infect pertains.
    def time2transition_rvs(index2infects, current_time, tau_decay):
        if Infect.delta is None:
            return Infect_Continuous.time2transition_rvs(index2infects, current_time, tau_decay)
        else:
            return Infect_Binomial.time2transition_rvs(index2infects, current_time, tau_decay)

def main():
    pass # Other tests are in the derived classes, Infect_Continuous and Infect_Discrete.
    
if __name__ == "__main__":
    main()
