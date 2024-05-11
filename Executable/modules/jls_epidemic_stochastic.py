#!/usr/bin/env python
"""
Simulates a stochastic Epidemic with watchers that observe after each epidemic event (infection, decay, etc.).
"""
from jls_heap_dictionary import Heap_Dictionary # A heap maintains the non-Markovian events.
from jls_epidemic_infect import Infect
from jls_epidemic_infect_mgr import Infect_Mgr
from jls_epidemic_compartment import Compartment
from jls_epidemic_util import resolve_collision
# transit = (compartment_from.index, compartment_to.index)
# transition = transit2count = {(compartment_from.index, compartment_to.index): count}
# event = time2_transit2count = {time: {(compartment_from.index, compartment_to.index): count}}
class Epidemic():
    def __init__(self, compartments): # initial counts in Markov compartments (with class count)
        assert compartments == Compartment.compartments
        self.events = Heap_Dictionary() # Holds (possibly non-Markov) {time:transition}.
        self.events.set_resolve_collision(resolve_collision) # Resolves collisions when a new event has the same time as an event in self.events.
    # Callback objects have watchers for epidemic output (see jls_callback.py).
    #   The default callbock, [], does nothing.
    def initialize(self, callback=[]):
        for c in callback.watchers:
            if '_obj' in c:
                raise ValueError('Epidemic watchers reserve the keyword "_obj" for the epidemic itself.')
            c['_obj'] = self # Each epidemic watcher has a key ['_obj'] that holds epidemic.
        self.callback = callback # The callback observes the events.
        self.current_time = 0
        for i,c in enumerate(Compartment.compartments):
            event = c.decay_time2transition_rvs(self.current_time, Infect.totals[i])
            if event is not None:
                self.events.push(event)
        self.callback.execute()
    # Updates the epidemic until every callback has seen enough or the epidemic has terminated.
    def evolve(self):
        while self.update():
            pass
    # Updates structures with the next transition and returns True if simulation should continue.
    #   Execution terminates if: (1) callback returns False; (2) no compartment has an event and events.is_empty().
    #   event = transit2count
    def update(self):
        tau_decay = self.events.peek()
        time2transition = Infect_Mgr.time2transition_rvs(Infect.index2infects, self.current_time, tau_decay)
        # Computes the consequences of the time2transition from Infect.
        if time2transition is None:
            if self.events.is_empty(): # No decay event can restore infectiousness.
                return False # Ends the simulation regardless of callbacks.
        else:
            assert len(time2transition) == 1
            self.events.push(time2transition)
        assert not self.events.is_empty()
        self.current_time = self.events.peek()
        transition0 = self.events.pop()
        assert transition0
        for transit,count in transition0.items():
            Infect.assert_transit(transit)
            i_from, i_to = transit
            time2transition = Compartment.compartments[i_to].decay_time2transition_rvs(self.current_time, count)
            if time2transition is not None:
                self.events.push(time2transition)
        Infect.transit(transition0)
        return self.callback.execute()

def main():
    pass
    
if __name__ == "__main__":
    main()
