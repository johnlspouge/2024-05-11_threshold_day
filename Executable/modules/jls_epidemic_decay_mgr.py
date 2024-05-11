#!/usr/bin/env python
"""
Epidemic Decay=Transit occurring spontaneously (never to an infected compartment)
"""
from jls_epidemic_decay_continuous import Continuous_Deterministic, Continuous_Gamma
from jls_epidemic_decay_discrete import Discrete_Deterministic, Discrete_Shifted_Neg_Binomial

class Decay_Mgr:
    # The argument delta determines whether the Infect is discrete or continuous.
    def factory(transit, mean_recovery_time, dispersion=None, delta=None):
        if delta is None:
            if dispersion is None:
                return Continuous_Deterministic(transit, mean_recovery_time)
            else:
                return Continuous_Gamma(transit, mean_recovery_time, dispersion)
        else:
            if dispersion is None:
                return Discrete_Deterministic(transit, mean_recovery_time, delta)
            else:
                return Discrete_Shifted_Neg_Binomial(transit, mean_recovery_time, dispersion, delta)

def main(): 
    pass
    
if __name__ == "__main__":
    main()
