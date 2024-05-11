#!/usr/bin/env python
"""
Adds watchers (dictionaries) for an epidemic simulation. 
    All watchers have a 'name' key.
    Initially, the watchers are either passive or active (['is_continue'] == None or True).
    An active watcher either watches a combination of specific compartment counts (either totals or differences of totals).
    After the combination exceeds a threshold eta the active watcher becomes inactive (['is_continue'] == False).
    The epidemic ceases if all individuals are in stable states, or if all active watchers have become inactive. 
"""
from jls_epidemic_compartment import Compartment
from jls_epidemic_infect import Infect
from jls_epidemic_watcher import totals, constant_times, ts_integer
from jls_callback import Callback # only here to facilitate debugging
#
# The following routine is necessary to connect the epidemic to the passive watchers.
#
# Returns the epidemic watchers, a list of dictionaries.
def to_watchers(epidemic): # Epidemic
    return epidemic.callback.watchers
# Returns the two standard passive_watchers.
def to_passive_watchers():
    return [
         {'callback':totals,'name':'all_times','is_continue':None}, # passive watcher at event times
         {'callback':constant_times,'name':'integer_times','ts':ts_integer,'is_continue':None} # passive watcher at integer times (e.g., days)
    ]
#
# The following routines access active_watchers.
#
# Returns an active watcher for the total of the last element of a passive watcher.
def active_for_total(
        passive_watcher, # Tracks totals in all compartments, possibly at integer days.
        eta, # threshold for stopping the computation: total > eta
        compartments): # Sums differences over these compartments, e.g., to find new cases from 'i'+'r'.
    k2v = {'callback':_is_total_less, 'name':'total',
           'passive_watcher':passive_watcher, 'eta':eta, 'compartments':compartments,
           'is_continue':True,'out':None}
    return k2v
# Stops watching if a sum of counts over the watcher['compartments'] > eta.
def _is_total_less(watcher):
    passive_watcher = watcher['passive_watcher']
    eta = watcher['eta']
    compartments = watcher['compartments']
    if ('out' not in passive_watcher or
        't' not in passive_watcher['out'] or
        len(passive_watcher['out']['t']) == 0):
        return True
    total = sum([passive_watcher['out'][c][-1] for c in compartments])
    if total > eta:
        watcher['out'] = {k:[v[-1]] for k,v in passive_watcher['out'].items()}
        return False
    return True
# Appends an active watcher of the difference of the last two elements of a passive watcher.
#   With passive_integer_times(), the active watcher can watch for new cases over a threshold eta.
# Returns the watcher.
def active_for_difference(
        passive_watcher, # Tracks totals in all compartments, possibly at integer days.
        eta, # threshold for stopping the computation: difference > eta
        compartments): # Sums differences over these compartments, e.g., to find new cases from 'i'+'r'.
    k2v = {'callback':_is_difference_less, 'name':'difference',
           'passive_watcher':passive_watcher, 'eta':eta, 'compartments':compartments,
           'is_continue':True,'out':None}
    return k2v
# Stops watching if the consecutive difference of a sum of counts over watcher['compartments'] > eta.
# Stops computation if the consecutive difference of a sum of counts over certain compartments > eta.
def _is_difference_less(watcher):
    passive_watcher = watcher['passive_watcher']
    eta = watcher['eta']
    compartments = watcher['compartments']
    if ('out' not in passive_watcher or
        't' not in passive_watcher['out'] or
        len(passive_watcher['out']['t']) <= 1):
        return True
    difference = sum([passive_watcher['out'][c][-1]-passive_watcher['out'][c][-2] for c in compartments])
    if difference > eta:
        watcher['out'] = {k:[v[-1]-v[-2]] for k,v in passive_watcher['out'].items()}
        watcher['out']['t'] = [passive_watcher['out']['t'][-2]] # Overwrites 't' difference.
        return False
    return True
#
# The following routines determine the state of the watchers.
#
# Returns True is the epidemic watcher is a passive watcher.
def is_passive(watcher): 
    return watcher['is_continue'] is None
# Returns True is the epidemic watcher is active.
def is_active(watcher): 
    return watcher['is_continue'] == True
# Returns True is the epidemic watcher is inactive.
def is_inactive(watcher): 
    return watcher['is_continue'] == False
#
# The following routines test the watchers.
#
# The watchers observe during an epidemic, but they report independently.    
def test_total_all_times():
    Compartment.initialize_from(None, ('s',20),('i',2))
    Compartment.initialize_infect(Compartment.compartments)
    epidemic = Epidemic0()
    # Tests stopping_time with constant_times.
    watchers = to_passive_watchers()
    passive_0 = watchers[0]
    assert passive_0 == {'callback':totals,'name':'all_times','is_continue':None}
    passive_1 = watchers[1]
    assert passive_1 == {'callback':constant_times,'name':'integer_times','ts':ts_integer,'is_continue':None}
    ETA = 5
    compartments = ['i']
    active_2 = active_for_total(passive_0, ETA, compartments)
    watchers.append(active_2)
    assert watchers == [passive_0, passive_1, active_2] 
    # Watches epidemic.
    epidemic.initialize(watchers, is_add_obj=True)
    assert [is_passive(w) for w in watchers] == [True, True, False] # The passive status is permanent.
    assert [is_inactive(w) for w in watchers] == [False, False, False] # The passive status is permanent.
    epidemic.watch()
    assert [w['out'] for w in watchers] == [
        {"t": [0.0], "s": [20], "i": [2]}, 
        {"t": [0.0], "s": [20], "i": [2]}, 
        None] 
    
    epidemic.current_time += 0.5
    Infect.transit({(0,1):1})
    epidemic.watch()
    assert [w['out'] for w in watchers] == [
        {"t": [0.0, 0.5], "s": [20, 19], "i": [2, 3]},
        {"t": [0.0], "s": [20], "i": [2]}, 
        None] 

    epidemic.current_time += 0.6
    Infect.transit({(0,1):2})
    epidemic.watch()
    assert [w['out'] for w in watchers] == [
        {"t": [0.0, 0.5, 1.1], "s": [20, 19, 17], "i": [2, 3, 5]},
        {"t": [0.0, 1.0], "s": [20, 19], "i": [2, 3]}, 
        None] 

    epidemic.current_time += 1.0
    Infect.transit({(0,1):4})
    epidemic.watch()
    assert [w['out'] for w in watchers] == [
        {"t": [0.0, 0.5, 1.1, 2.1], "s": [20, 19, 17, 13], "i": [2, 3, 5, 9]},
        {"t": [0.0, 1.0, 2.0], "s": [20, 19, 17], "i": [2, 3, 5]}, 
        {"t": [2.1], "s": [13], "i": [9]}] # 'i' > ETA=5, so the computation ends and the watcher gives the last total.
    assert not active_2['is_continue']
    assert [is_passive(w) for w in watchers] == [True, True, False] # The passive status is permanent.
    assert [is_inactive(w) for w in watchers] == [False, False, True] # The passive status is permanent.
# The watchers observe during an interaction with the epidemic, but they report independently.    
def test_difference_integer_times():
    Compartment.initialize_from(None, ('s',20),('i',2))
    Compartment.initialize_infect(Compartment.compartments)
    epidemic = Epidemic0()
    # Tests stopping_time with constant_times.
    watchers = to_passive_watchers()
    passive_0 = watchers[0]
    assert passive_0 == {'callback':totals,'name':'all_times','is_continue':None}
    passive_1 = watchers[1]
    assert passive_1 == {'callback':constant_times,'name':'integer_times','ts':ts_integer,'is_continue':None}
    ETA = 1
    compartments = ['i']
    active_2 = active_for_difference(passive_1, ETA, compartments)
    watchers.append(active_2)
    assert watchers == [passive_0, passive_1, active_2] 
    # Watches epidemic.
    epidemic.initialize(watchers, is_add_obj=True)
    epidemic.watch()
    assert [w['out'] for w in watchers] == [
        {"t": [0.0], "s": [20], "i": [2]}, 
        {"t": [0.0], "s": [20], "i": [2]}, 
        None] 
    epidemic.current_time += 0.5
    Infect.transit({(0,1):1})
    epidemic.watch()
    assert [w['out'] for w in watchers] == [
        {"t": [0.0, 0.5], "s": [20, 19], "i": [2, 3]},
        {"t": [0.0], "s": [20], "i": [2]}, 
        None] 

    epidemic.current_time += 0.6
    Infect.transit({(0,1):2})
    epidemic.watch()
    assert [w['out'] for w in watchers] == [
        {"t": [0.0, 0.5, 1.1], "s": [20, 19, 17], "i": [2, 3, 5]},
        {"t": [0.0, 1.0], "s": [20, 19], "i": [2, 3]}, 
        None] 

    epidemic.current_time += 1.0
    Infect.transit({(0,1):4})
    epidemic.watch()
    assert [w['out'] for w in watchers] == [
        {"t": [0.0, 0.5, 1.1, 2.1], "s": [20, 19, 17, 13], "i": [2, 3, 5, 9]},
        {"t": [0.0, 1.0, 2.0], "s": [20, 19, 17], "i": [2, 3, 5]}, 
        {'t': [1.0], 's': [-2], 'i': [2]}] # 'i' > ETA=1, so the computation ends and the watcher gives the integer day of the final difference i[day+1]-i[day].
    assert not active_2['is_continue']
# a crippled Epidemic class for testing purposes only
class Epidemic0:
    def __init__(self):
        self.current_time = 0.0
    def initialize(self, callback, is_add_obj=False):
        if is_add_obj:
            for c in callback:
                c['_obj'] = self
        self.callback = Callback(callback)
    def watch(self):
        return self.callback.execute()
    
def main(): 
    test_total_all_times()
    test_difference_integer_times()
      
if __name__ == "__main__":
    main()
