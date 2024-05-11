#!/usr/bin/env python
"""
Epidemic Infect=Transit caused by Infectious compartment
"""
from math import isclose
from scipy.stats import poisson, nbinom

# The argument dispersion determines whether the individual_r0_pmf is Poisson or negative binomial.
def individual_r0_pmf(k, beta, mean_recovery_time, dispersion=None):
    if dispersion is None:
        return poisson.pmf(k, beta*mean_recovery_time)
    else: # continuous gamma recovery time
        p = dispersion/(beta*mean_recovery_time+dispersion) # dispersion/(r0+dispersion)
        return nbinom.pmf(k, dispersion, p)

def test_individual_r0_pmf():
    k = 2
    beta = 3.0
    mean_recovery_time = 5.0
    assert isclose(poisson.pmf(k, beta*mean_recovery_time), 3.4414011056455366e-05)
    assert isclose(individual_r0_pmf(k, beta, mean_recovery_time, dispersion=None), 3.4414011056455366e-05)
    dispersion = 10.0
    p = dispersion/(beta*mean_recovery_time+dispersion) 
    assert isclose(p, 0.4)
    assert isclose(nbinom.pmf(k, dispersion, p), 0.0020761804800000097)
    assert isclose(individual_r0_pmf(k, beta, mean_recovery_time, dispersion), 0.0020761804800000097)
    
def main():
    test_individual_r0_pmf() # Other tests are in the derived classes, Infect_Continuous and Infect_Discrete.
    
if __name__ == "__main__":
    main()
