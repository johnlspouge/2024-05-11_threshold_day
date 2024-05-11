#!/usr/bin/env python
"""
Epidemic Compartments
"""
from math import isclose
import numpy as np

from jls_epidemic_infect import Infect
from jls_epidemic_infect_mgr import Infect_Mgr
from jls_epidemic_decay_mgr import Decay_Mgr
from jls_epidemic_prob0decay import Prob0Decay

# The present version assumes all random times for infection are distinct, being from a Poisson process.
# transit=(compartment_from.index, compartment_to.index)
# transition=transit2count
class Compartment:
    compartments = None
    name2index = None
    # Sets compartments from a tuple of pair tuples (name0total).
    #   jls_epidemic_watcher.py provides all output from the epidemic.
    #   It ensures that the output accounts for the time-grid delta, if delta > 0.0.
    #   The internal representation in the computation is integer time i, but the output is i*delta.
    def initialize_from(delta, *args): # *args=name0totals, delta is the time-step (None for continuous)
        compartments = [Compartment(*name0total) for name0total in args]
        Compartment.initialize(compartments)
        Infect.set_delta(delta)
    # Adds a decay from the compartment name_from to name_to.
    def add_decay_to(name_from, name_to, mean_recovery_time, dispersion=None, prob=1.0): 
        Compartment.compartments[Compartment.name2index[name_from]].add_decay(name_to, mean_recovery_time, dispersion, prob)
    # Adds the ability to the compartment name_infectious to infect at rate beta compartment name_from, converting it to name_to.
    def add_infect_to(name_infectious, beta, name_from, name_to, mean_recovery_time=None, is_poisson=False): # Adds an infect.
        Compartment.compartments[Compartment.name2index[name_infectious]].add_infect(beta, name_from, name_to, mean_recovery_time, is_poisson)
    # Sets compartments, compartment.name, compartment.index, Compartment.name2index, Infect.total, and Infect.total.
    def initialize(compartments):
        Compartment.totals = [0] * len(compartments)
        Compartment.compartments = compartments
        for i,c in enumerate(compartments):
            c.index = i
        Compartment.name2index = {c.name:c.index for c in compartments}
        if len(Compartment.name2index) < len(compartments):
            raise ValueError(f'Each compartment must have a unique name.')
        if 't' in Compartment.name2index:
            raise ValueError(f'The name "t" is reserved for time in the output.')
        Compartment.initialize_infect(compartments)
    # Initializes Infect.total and Infect.total.
    def initialize_infect(compartments):
        for c in compartments:
            assert isinstance(c.total, int) and 0 <= c.total
        Infect.set_totals([c.total for c in compartments])
    # Wraps up the initiation of the compartments.
    def wrapup(nu=1.0):
        Infect.set_nu(nu)
        Infect.clear_index2infects()
        delta = None
        for i,c in enumerate(Compartment.compartments):
            if hasattr(c, 'infects'):
                if not hasattr(c, 'decays'):
                    raise ValueError(f'If compartment "{c.name}" infects, it must decay.')
                Infect.set_index2infect(c.infects, i)  # Registers infectious compartments.
            if hasattr(c, 'decays'): 
                if delta is None:
                    delta = c.decays[0].delta
                if c.decays[0].delta != delta:
                    raise ValueError(f'The grid delta "{delta}" for compartment "{c.name}" does not match the grids for other compartments.')
                Prob0Decay.error(c.probs, c.decays) # Construction ensures the len(self.probs) == len(self.decays).
                c.prob0decay = Prob0Decay(c.probs, c.decays)
                c.decay_time2transition_rvs = c.prob0decay.time2transition_rvs
        assert Infect.index2infects # An epidemic requires an infectious compartment.
    # Constructs an Epidemic Compartment.
    def __init__(self, name, total=0):
        self.name = name
        self.total = total # initial total - Infect handles updated compartment totals during the epidemic.
    def name(self):
        return self.name
    def get_index(self):
        return self.index
    def add_decay(self, name_to, mean_recovery_time, dispersion=None, prob=1.0):
        # Delays the following until Compartment.wrapup():
        #   self.decay_time2transition_rvs = self.prob0decays.time2transition_rvs
        #   assert isclose(sum(self.probs), 1.0)
        index_to = Compartment.name2index[name_to]
        transit = (self.index, index_to)
        decay = Decay_Mgr.factory(transit, mean_recovery_time, dispersion, Infect.delta)
        if hasattr(self, 'decays'):
            self.probs.append(prob)
            self.decays.append(decay)
        else:
            self.probs = [prob]
            self.decays = [decay]
    def decay_time2transition_rvs(self, current_time, count): # decay
        return None
    def add_infect(self, beta, name_from, name_to, mean_recovery_time=None, is_poisson=False):
        if not hasattr(self, 'decays'):
            raise Exception('Before becoming infectious, a compartment must be given a decay mode.')
        Prob0Decay.error(self.probs, self.decays)
        index_from = Compartment.name2index[name_from]
        index_to = Compartment.name2index[name_to]
        transit = (index_from, index_to)
        # To match continuous individual R0, deterministic decay uses Poisson infection; gamma, NegBinomial.
        infect = Infect_Mgr.factory(self.index, beta, transit, Infect.delta, mean_recovery_time, is_poisson)
        if hasattr(self, 'infects'):
            assert infect not in self.infects
            self.infects.append(infect)
        else:
            self.infects = [infect]    
# Tests initialization with variables.
def test_Compartment():
    Infect.delta = None
    susceptible = Compartment(name='Susceptible', total=20)
    infected = Compartment(name='Infected', total=2)
    recovered = Compartment(name='Recovered', total=0)
    Compartment.initialize((susceptible, infected, recovered))
    infected.add_decay('Recovered', mean_recovery_time=5.0) # Adds a decay.
    infected.add_infect(beta=2.0, name_from='Susceptible', name_to='Infected') # Adds an infect.
    Compartment.wrapup()
    # Checks construction.
    tuples = (('Susceptible', 20), ('Infected', 2), ('Recovered', 0))
    assert len(Compartment.compartments) == 3
    for i,c in enumerate(Compartment.compartments):
        assert (c.name, c.total) == tuples[i] and c.index == i
# Tests initialization with variables.
def test_Name_Interface():
    # s,i,r for branching process
    Compartment.initialize_from(None, ('s',10000000000), ('i',1), ('r',0))
    Compartment.add_decay_to(name_from='i', name_to='r', mean_recovery_time=5.0) # Adds a decay.
    Compartment.add_infect_to(name_infectious='i', beta=2.0, name_from='s', name_to='i') # Adds an infect.
    Compartment.wrapup()
    assert Infect.total == 10000000001
    assert Infect.totals == [10000000000,1,0]
    assert len(Infect.index2infects) == 1 and 1 in Infect.index2infects
    # s,i,i0,r for variant
    Compartment.initialize_from(None, ('s',20), ('i',2), ('i0',2), ('r',0))
    Compartment.add_decay_to(name_from='i', name_to='r', mean_recovery_time=5.0, dispersion=1.0) # Adds a decay.
    Compartment.add_infect_to(name_infectious='i', beta=2.0, name_from='s', name_to='i') # Adds an infect.
    Compartment.add_decay_to(name_from='i0', name_to='r', mean_recovery_time=5.0, dispersion=1.0) # Adds a decay.
    Compartment.add_infect_to(name_infectious='i0', beta=2.0, name_from='s', name_to='i0') # Adds an infect.
    Compartment.wrapup()
    assert Infect.total == 24
    assert Infect.totals == [20,2,2,0]
    assert len(Infect.index2infects) == 2
    assert 1 in Infect.index2infects and 2 in Infect.index2infects
# Tests the infectious transits and compartment_to r0s that compartment c_index causes in compartment_froms.
def test_Add_Decay_Infect():
    np.random.seed(31415)
    # Checks infection with one pathway.
    susceptible = Compartment(name='Susceptible', total=8)
    infected = Compartment(name='Infected', total=2)
    recovered = Compartment(name='Recovered', total=0)
    Compartment.initialize((susceptible, infected, recovered))
    for i,c in enumerate(Compartment.compartments):
        assert c.index == i
    infected.add_decay('Recovered', mean_recovery_time=5.0) # Adds a continuous deterministic decay.
    infected.add_infect(beta=1.25, name_from='Susceptible', name_to='Infected') # Adds a continuous Infect.
    Compartment.wrapup(nu=0.5)
    s_index, i_index, r_index = (0,1,2)
    intensity = infected.infects[0].intensity_of_infection()
    intensity0 = Infect.nu*infected.infects[0].beta*Infect.totals[i_index]*Infect.totals[s_index]/Infect.total
    assert isclose(intensity, 1.0) and isclose(intensity0, 1.0)
    N = 1000
    time2transition = infected.decay_time2transition_rvs(current_time=10.0, count=N)
    assert time2transition == {15.0: {(1, 2): N}}
    # Checks infection with two pathways.
    susceptible = Compartment(name='Susceptible', total=4)
    infected = Compartment(name='Infected', total=2)
    recovered = Compartment(name='Recovered', total=4)
    Compartment.initialize((susceptible, infected, recovered))
    infected.add_decay('Recovered', mean_recovery_time=5.0) # Adds a continuous deterministic decay.
    infected.add_infect(beta=1.25, name_from='Susceptible', name_to='Infected') # Adds a continuous Infect.
    infected.add_infect(beta=0.625, name_from='Recovered', name_to='Infected') # Adds a second pathway for a continuous Infect.
    Compartment.wrapup(nu=0.5)
    intensities = [infected.infects[0].intensity_of_infection(), 
                   infected.infects[1].intensity_of_infection()]
    intensities0 = [Infect.nu*infected.infects[0].beta*Infect.totals[i_index]*Infect.totals[s_index]/Infect.total,
                    Infect.nu*infected.infects[1].beta*Infect.totals[i_index]*Infect.totals[r_index]/Infect.total]
    assert np.allclose(intensities, intensities0) and np.allclose(intensities, (0.5,0.25))
    # Checks infection with alternative pathways.
    N = 1000
    count_s = times_s = times_r = 0
    current_time=10.0
    for i in range(N):
        time2transition = Infect_Mgr.time2transition_rvs(Infect.index2infects, current_time=10.0, tau_decay=None)
        assert len(time2transition) == 1
        for time,transition in time2transition.items():
            assert transition == {(0, 1): 1} or transition == {(2, 1): 1} 
            if transition == {(0, 1): 1}:
                count_s += 1
                times_s += time
            else:
                times_r += time
    assert isclose(times_s/count_s-current_time, 1.250855614079953) # about 4/3=1(0.5+0.25)
    assert isclose(times_r/(N-count_s)-current_time, 1.336102245825359) # about 4/3=1(0.5+0.25)
    assert isclose(count_s/N, 2/3, abs_tol=0.005)
    # Adds an exposed compartment tp infection with alternative pathways.
    susceptible = Compartment(name='Susceptible', total=20)
    exposed = Compartment(name='Exposed', total=0)
    infected = Compartment(name='Infected', total=2)
    recovered = Compartment(name='Recovered', total=0)
    Compartment.initialize((susceptible, exposed, infected, recovered))
    infected.add_decay('Recovered', mean_recovery_time=5.0) # Adds a continuous deterministic decay.
    infected.add_infect(beta=2.0, name_from='Susceptible', name_to='Exposed') # Adds a continuous Infect.
    infected.add_infect(beta=4.0, name_from='Susceptible', name_to='Infected') # Adds a second pathway for a continuous Infect.
    Compartment.wrapup(nu=0.5)
    # Checks an exposed compartment tp infection with alternative pathways.
    N = 1000
    count_e = times_e = times_i = 0
    current_time=10.0
    for i in range(N):
        time2transition = Infect_Mgr.time2transition_rvs(Infect.index2infects, current_time=10.0, tau_decay=None)
        assert len(time2transition) == 1
        for time,transition in time2transition.items():
            assert transition == {(0, 1): 1} or transition == {(0, 2): 1} 
            if transition  == {(0, 1): 1}:
                count_e += 1
                times_e += time
            else:
                times_i += time
    assert isclose(times_e/count_e-current_time, 0.18443597412959) # about 1/6=1(2.0+4.0)
    assert isclose(times_i/(N-count_e)-current_time, 0.176183633944174) # about 1/6=1(2.0+4.0)
    assert isclose(count_e/N, 1/3, abs_tol=0.001)

def main():
    test_Compartment()
    test_Add_Decay_Infect()
    test_Name_Interface()
    
if __name__ == "__main__":
    main()