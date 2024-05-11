#!/usr/bin/env python
"""
Computes percentiles and confidence intervals. 
"""
from math import isclose
import numpy as np
from scipy.stats import norm, bootstrap
from collections import Counter

# https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.bootstrap.html
class Quantile:
    quantile_ = None # quantile static variable
    
    #   The standard quantile definition always rounds up.
    def set(quantile):
        Quantile.quantile_ = quantile
    # The callable statistic calculates the quantile in scipy.stats.bootstrap.
    #   https://numpy.org/doc/stable/reference/generated/numpy.percentile.html
    #   The quantile follows the linear numpy default
    #   https://numpy.org/doc/stable/reference/generated/numpy.quantile.html
    def quantile(samples):
        return np.quantile(samples, Quantile.quantile_)

def test_quantile():
    samples = np.array([2,3,1,5,4])
    quantiles2value = {0.2:1.8,0.3:2.2,0.7:3.8,0.8:4.2}
    for p,v in quantiles2value.items():
        Quantile.set(p)
        assert isclose(Quantile.quantile(samples), v)
    samples = np.array([2,3,0,6,4])
    quantiles2value = {0.4:2.6,0.5:3.0,0.6:3.4}
    for p,v in quantiles2value.items():
        Quantile.set(p)
        assert isclose(Quantile.quantile(samples), v)
# Converts a counter to the corresponding list of samples.
def counter2samples(counter):
    d = dict(counter)
    samples = []
    for k in sorted(d.keys()):
        for i in range(d[k]):
            samples.append(k)
    return samples
        
def test_counter2samples():
    samples0 = [0, 0, 1, 0, 1, 2, 0, 1, 2, 3]
    counter = Counter(samples0)
    assert counter2samples(counter) == [0,0,0,0,1,1,1,2,2,3]
    # Reverses the samples
    samples0.reverse()
    counter = Counter(samples0)
    assert counter2samples(counter) == [0,0,0,0,1,1,1,2,2,3]
    # Has an element with 0 count.
    counter = {2:2,0:0,3:3,1:1,4:4}
    assert counter2samples(counter) == [1,2,2,3,3,3,4,4,4,4]

def test_bootstrap():
    np.random.seed(31415)
    dist = norm(loc=2, scale=4)  # our "unknown" distribution
    assert isclose(dist.std(), 4.0)
    data = dist.rvs(size=100)
    assert isclose(np.std(data), 3.5744389405719175)
    data = (data,)  # samples must be in a sequence
    res = bootstrap(data, np.std, confidence_level=0.9)
    assert isclose(res.standard_error, 0.23866111637883353)
    assert isclose(res.confidence_interval.low, 3.2343776431541524)
    assert isclose(res.confidence_interval.high, 4.0389891526854385)
    data = (np.array(list(range(20))),)
    res = bootstrap(data, np.std, confidence_level=0.9)
    assert isclose(res.confidence_interval.low, 4.972675336275233)
    assert isclose(res.confidence_interval.high, 6.877043555568017)
    Quantile.set(0.4)
    res = bootstrap(data, Quantile.quantile, confidence_level=0.9)
    assert isclose(res.confidence_interval.low, 4.0)
    assert isclose(res.confidence_interval.high, 11.0)

def main(): 
    test_quantile()
    test_counter2samples()
    test_bootstrap()
    
if __name__ == "__main__":
    main()
