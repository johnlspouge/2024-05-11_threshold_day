#!/usr/bin/env python
"""
Callbacks for watching epidemic states
"""
import pandas as pd

from jls_epidemic_compartment import Compartment
from jls_epidemic_infect import Infect
from jls_callback import Callback # only here to facilitate debugging

# Returns a df containing the epidemic event times and the corresponding compartment counts.
def to_df(out):
    for k,v in out.items():
        if len(v) != len(out['t']):
            raise ValueError(f'len(out[{k}]) != {len(out["t"])} = len(out["t"])')
    df = pd.DataFrame(list(zip(*list(out.values()))), columns=list(out.keys())) 
    return df
    
def test_to_df():
    ts = list(range(10,70,10))
    js = list(range(6,0,-1))
    ks = list(range(0,12,2))
    out = {'t':ts,'j':js,'k':ks}
    df = to_df(out)
    df0 = pd.DataFrame(list(zip(ts, js, ks)), columns = ['t','j','k']) 
    assert df.equals(df0)
#
# The following routines watch the epidemic and update watchers as aappropriate.
#
# Watches & appends pt = {compartment: compartmental count} to watcher['out'].
def totals(watcher):
    pt = total(watcher)
    total_append(pt, watcher['out'])
# Appends pt = {compartment: compartmental count} to watcher['out'].
def total_append(pt, out):
    for k,v in pt.items():
        out[k].append(pt[k])
# Watches time t for each epidemic event & determines pt = {compartment: compartmental count}.
#   The time t is Infect.delta ? epidemic.current_time : epidemic.current_time*Infect.delta. 
#   All adjustments for 0.0 < Infect.delta should be done by total(watcher).
#   Then, the time-grid can be stored as integers.
#   The compartment counts are in Infect.totals.
def total(watcher):
    if 'out' not in watcher: 
        out = watcher['out'] = {}
        out['t'] = []
        for c in Compartment.compartments:
            out[c.name] = []
    epidemic = watcher['_obj']
    # Multiplies t by Infect.delta if 0.0 < Infect.delta:.
    def _delta_adjust(t):
        if Infect.delta is not None:
            assert isinstance(Infect.delta, float) and 0.0 < Infect.delta
            return t*Infect.delta
        return t
    pt = {}
    pt['t'] = _delta_adjust(epidemic.current_time)
    for i,c in enumerate(Compartment.compartments):
        pt[c.name] = Infect.totals[i]
    return pt

def test_totals():
    Compartment.initialize_from(None, ('s',20),('i',0))
    Compartment.initialize_infect(Compartment.compartments)
    epidemic = Epidemic0()
    # Tests totals.
    watchers = [{'callback':totals,'_obj':epidemic}] 
    epidemic.initialize(watchers)
    epidemic.watch()
    data0 = del_ptrs(watchers[0],('callback', '_obj'))
    assert data0 == {'is_continue': None, 'out': {'t': [0.0], 's': [20], 'i': [0]}}
    epidemic.current_time += 1.0
    Infect.transit({(0,1):2})
    epidemic.watch()
    data0 = del_ptrs(watchers[0],('callback', '_obj'))
    assert data0 == {'is_continue': None, 'out': {'t': [0.0, 1.0], 's': [20, 18], 'i': [0, 2]}}
# The watcher has watcher['ts'](index) to map {0,1,2,...} to a set of strictly increasing positive reals.
#   Behaves like total(watcher), except it only records the time t & compartmental counts at the times watcher['ts'](index). 
def constant_times(watcher):
    if 'index' not in watcher: 
        assert 'ts' in watcher # The function watcher['ts'] maps 'index' to watching time.
        watcher['index'] = 0 # the next time index in ts
        watcher['pt0'] = None # previous point
    pt0 = watcher['pt0']
    pt = total(watcher) # adjusted for Infect.delta
    t = watcher['ts'](watcher['index']) # not adjusted for Infect.delta
    while t is not None:
        if pt['t'] < t or pt0 is None and t < pt['t']:
            watcher['pt0'] = pt
            return
        elif t == pt['t']:
            pt_star = pt.copy()
        else:
            pt_star = pt0.copy()
        watcher['pt0'] = pt
        watcher['index'] += 1
        pt_star['t'] = t
        t = watcher['ts'](watcher['index'])
        total_append(pt_star, watcher['out'])
# Maps integer indexes to corresponding integer times for index set {0,1,2,...}
def ts_integer(i):
    return i
    
def test_constant_times():
    Compartment.initialize_from(None, ('s',20),('i',2))
    Compartment.initialize_infect(Compartment.compartments)
    epidemic = Epidemic0()
    # Tests totals.
    watchers = [{'callback':constant_times,'_obj':epidemic,'ts':ts_integer}] 
    epidemic.initialize(watchers)
    epidemic.watch()
    data0 = del_ptrs(watchers[0],('callback', '_obj', 'ts', 'pt0'))
    assert data0 == {'is_continue': None, 'index': 1, 'out': {'t': [0.0], 's': [20], 'i': [2]}}
    epidemic.current_time += 0.5
    Infect.transit({(0,1):1})
    epidemic.watch()
    data0 = del_ptrs(watchers[0],('callback', '_obj', 'ts', 'pt0'))
    assert data0 == {'is_continue': None, 'index': 1, 'out': {'t': [0.0], 's': [20], 'i': [2]}}
    epidemic.current_time += 0.6
    Infect.transit({(0,1):2})
    epidemic.watch()
    data0 = del_ptrs(watchers[0],('callback', '_obj', 'ts', 'pt0'))
    assert data0 == {'is_continue': None, 'index': 2, 'out': {'t': [0.0, 1.0], 's': [20, 19], 'i': [2, 3]}}
    epidemic.current_time += 1.0
    Infect.transit({(0,1):4})
    epidemic.watch()
    data0 = del_ptrs(watchers[0],('callback', '_obj', 'ts', 'pt0'))
    assert data0 == {'is_continue': None, 'index': 3, 'out': {'t': [0.0, 1.0, 2.0], 's': [20, 19, 17], 'i': [2, 3, 5]}}

### 1 DEPRECATED - For backward compatibility only.

# active watcher (out['is_continue']=True)
#   Requires out['condition'], a lambda function of the current epidemic counts that returns True or False.
#   Terminates watching with False at a stopping_time, and watchers['out'] counts the epidemic compartments at the stopping time.
#   Adjustment for 0.0 < Infect.delta is done by total(watcher).
def stopping_time(watcher):
    if 'condition' not in watcher:
        ValueError('watcher = {watcher} needs a stopping condition.')
    if 'out' not in watcher: 
        out = watcher['out'] = {}
        out['t'] = []
        for c in Compartment.compartments:
            out[c.name] = []
    pt = total(watcher)
    if watcher['condition'](pt):
        total_append(pt, watcher['out'])
        return False
    return True

def test_stopping_time():
    Compartment.initialize_from(None, ('s',20),('i',2))
    Compartment.initialize_infect(Compartment.compartments)
    epidemic = Epidemic0()
    # Tests stopping_time with constant_times.
    THRESHOLD = 6
    watchers = [
            {'callback':totals,'_obj':epidemic},
            {'callback':constant_times,'_obj':epidemic,'ts':ts_integer},
            {'callback':stopping_time,'_obj':epidemic,'condition':lambda pt: pt['i'] > THRESHOLD,'is_continue':True}
    ] 
    epidemic.initialize(watchers)
    epidemic.watch()
    assert epidemic.callback.out() == [
        {"t": [0.0], "s": [20], "i": [2]},
        {"t": [0.0], "s": [20], "i": [2]}, 
        {"t": [], "s": [], "i": []}] 

    epidemic.current_time += 0.5
    Infect.transit({(0,1):1})
    epidemic.watch()
    assert epidemic.callback.out() == [
        {"t": [0.0, 0.5], "s": [20, 19], "i": [2, 3]},
        {"t": [0.0], "s": [20], "i": [2]}, 
        {"t": [], "s": [], "i": []}] 

    epidemic.current_time += 0.6
    Infect.transit({(0,1):2})
    epidemic.watch()
    assert epidemic.callback.out() == [
        {"t": [0.0, 0.5, 1.1], "s": [20, 19, 17], "i": [2, 3, 5]},
        {"t": [0.0, 1.0], "s": [20, 19], "i": [2, 3]}, 
        {"t": [], "s": [], "i": []}] 

    epidemic.current_time += 1.0
    Infect.transit({(0,1):4})
    epidemic.watch()
    assert epidemic.callback.out() == [
        {"t": [0.0, 0.5, 1.1, 2.1], "s": [20, 19, 17, 13], "i": [2, 3, 5, 9]},
        {"t": [0.0, 1.0, 2.0], "s": [20, 19, 17], "i": [2, 3, 5]}, 
        {"t": [2.1], "s": [13], "i": [9]}] 

### 0 DEPRECATED

# a crippled Epidemic class
class Epidemic0:
    def __init__(self):
        self.current_time = 0.0
    def initialize(self, watchers):
        self.callback = Callback(watchers)
    def watch(self):
        return self.callback.execute()
# Deletes keys from watchers to promote brevity of debugging comparisons.
def del_ptrs(watchers, keys):
    data0 = watchers.copy()
    for k in keys:
        del data0[k]
    return data0
    
def main(): 
    test_to_df()
    test_totals()
    test_constant_times()
    test_stopping_time()
    
if __name__ == "__main__":
    main()
