#!/usr/bin/env python
"""
Epidemic utility routines 
"""
from copy import copy
from math import isclose, sqrt
import numpy as np
from statistics import fmean

# Asserts transit type.
def assert_transit(transit, length=None):
    assert isinstance(transit, tuple) 
    assert isinstance(transit[0], int) and isinstance(transit[1], int)
    assert 0 <= transit[0] and 0 <= transit[1] 
    if length is not None:
        assert transit[0] < length and transit[1] < length
# Returns a random element from a list of choices whose first elements are probabilities.
#   For consistent debugging avoids call to random if the list has one element.
def one_list_rv(p0lists):
    if len(p0lists) == 1: 
        p0list = p0lists[0]
        assert p0list[0] == 1.0
        index = 0
    else:
        ps = [k[0] for k in p0lists]
        assert isclose(sum(ps), 1.0)
        index = np.random.choice(range(len(p0lists)), p=ps)
    p0list_copy = copy(p0lists[index])
    p0list_copy.pop(0)
    list_copy = p0list_copy
    return list_copy
# Resolves collision for dictionaries transit2count by adding count when two transits are the same.
def resolve_collision(k, dictionary0, dictionary):
    assert dictionary.get(k) is not None
    if dictionary0.get(k) is None:
        dictionary0[k] = dictionary[k]
    else:
        add_count(dictionary0[k], dictionary[k])
# Resolves collision for dictionaries transition=transit2count by adding count when two transits are the same.
def add_count(transition0, transition): # event = {transit:count}
    for k,v in transition.items():
        if k in transition0:
            transition0[k] += v
        else:
            transition0[k] = v

def test_resolve_collision():
    dictionary0 = {1:{2:3,4:5},2:{2:-3,4:1},3:{5:-3,6:1}}
    dictionary = {1:{2:3,4:5},2:{2:-3,4:1},4:{5:-3,6:1}}
    resolve_collision(1, dictionary0, dictionary)
    assert dictionary0 == {1: {2: 6, 4: 10}, 2: {2: -3, 4: 1}, 3: {5: -3, 6: 1}}
    resolve_collision(2, dictionary0, dictionary)
    assert dictionary0 == {1: {2: 6, 4: 10}, 2: {2: -6, 4: 2}, 3: {5: -3, 6: 1}}
    resolve_collision(4, dictionary0, dictionary)
    assert dictionary0 == {1: {2: 6, 4: 10}, 2: {2: -6, 4: 2}, 3: {5: -3, 6: 1}, 4: {5: -3, 6: 1}}

# Normalizes population sizes to fractions of the original population size = Compartment.total.
def normalize(totals): # totals = pts[1] = [pts[1][i] for i,c in enumerate(Compartment.compartments)]
    # Checks totals is a matrix = rows of equal length with constant column sum
    def check(totals):
        length = len(totals[0])
        for ts in totals:
            assert len(ts) == length
        size = sum([ts[0] for ts in totals])
        for i,t in enumerate(totals[0]):
            assert sum([ts[i] for ts in totals]) == size
        return size
    size = check(totals)
    fractions = []
    for ts in totals:
        fractions.append([t/size for t in ts])
    return fractions
# Tests the normalization of totals by the population size = Compartment.total.
def test_one_list():
    np.random.seed(31415)
    p0lists = [[1.0,2]]
    oc = one_list_rv(p0lists)
    assert oc == [2]
    p0lists = [[0.0,0],[0.1,1],[0.2,2],[0.3,3],[0.4,4]]
    N = 1000
    assert fmean([one_list_rv(p0lists)[0] for i in range(N)]) == 2.998
# Tests the normalization of totals by the population size = Compartment.total.
def test_normalize():
    totals0 = [
        [8, 6, 2, 0], 
        [2, 4, 6, 9], 
        [0, 0, 2, 1]
    ]
    fractions0 = [
        [0.8, 0.6, 0.2, 0.0], 
        [0.2, 0.4, 0.6, 0.9], 
        [0.0, 0.0, 0.2, 0.1]
    ]
    results = normalize(totals0)
    for i,t in enumerate(results):
        assert np.allclose(results[i],fractions0[i])

# https://www.mikulskibartosz.name/wilson-score-in-python-example/
# Returns the Wilson score confidence interval.
def wilson(k, n, z = 1.96): # The default is 0.95 confidence interval. k = successes observed, n = trials
    p = k/n
    denominator = 1 + z**2/n
    centre_adjusted_probability = p + z*z / (2*n)
    adjusted_standard_deviation = sqrt((p*(1 - p) + z*z / (4*n)) / n)
    # bounds
    lower_bound = (centre_adjusted_probability - z*adjusted_standard_deviation) / denominator
    upper_bound = (centre_adjusted_probability + z*adjusted_standard_deviation) / denominator
    return (lower_bound, upper_bound)

def test_wilson(): 
    assert np.allclose(wilson(3, 12), (0.08894003962896137, 0.53231033912066))

def main():
    test_one_list()
    test_resolve_collision()
    test_normalize()
    test_wilson()
    
if __name__ == "__main__":
    main()
