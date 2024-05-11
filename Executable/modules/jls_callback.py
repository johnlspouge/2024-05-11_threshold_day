#!/usr/bin/env python
"""
A Callback object has a list of dictionaries, self.watchers.
    Each watcher must have a key ['callback'] whose value is a callback function.
    The callback function takes the watcher as an argument and returns True, False, or None.
    Each watcher has other keys whose values monitor and store the progress of a computation.
    Each watcher must have a key ['is_continue'] to indicate the watcher's state (initial state defaults to None).
    Each watcher is passive, active, or inactive (['is_continue'] == None, True, or False).
    Events in a computation trigger an update with the ['callback'] function for all watchers simultaneously.
    If a watcher is passive (['is_continue'] == None), it remains passive but its ['callback'] function is updated.
    If a watcher is active, its ['callback'] function is updated.
    If an active ['callback'] function returns False, the watcher becomes inactive (['is_continue'] = False)
    If all watchers have a falsy ['is_continue'] (None or False) then self.execute() returns False (which is useful for terminating computations).
"""
# The Callback class permits a class to monitor the progress or state of a second by adding a watcher. 
#   The watcher code can be altered without disturbing the code for the monitored class.
class Callback:
    # Initializes watcher['is_continue']=None, if it is absent.
    def __init__(self, watchers=[]): 
        self.watchers = watchers
        for watcher in self.watchers:
            if 'callback' not in watcher:
                raise KeyError("'callback' not in watcher")
            if 'is_continue' not in watcher:
                watcher['is_continue'] = None
    # Callback execution halts after not watcher['is_continue'], but continues passively if watcher['is_continue'] is None.
    def execute(self):
        if not self.watchers: # The computation does not use watchers, so Callback does not terminate it.
            return True 
        for watcher in self.watchers:
            if watcher['is_continue'] is None: # passive watcher
                watcher['callback'](watcher) # A passive watcher remains a possive watcher.
        is_continue = False
        for watcher in self.watchers:
            if watcher['is_continue']: # The watcher still awaits termination.
                watcher['is_continue'] = watcher['callback'](watcher)
                if watcher['is_continue']:
                    is_continue = True
        return is_continue # Return False only if all active watchers have been terminated.
    # Returns watchers[i][out], singly or all together, as a convenience to iterate over key=out in the watchers.
    #   Returns None if [out] not in watchers[i].
    def out(self, i=None, key='out'):
        if i is None:
            return [w.get(key, None) for w in self.watchers]
        elif len(self.watchers) <= i:
            raise ValueError(f'len(self.watchers) = {len(self.watchers)} <= {i} = i')
        else:
            return self.watchers[i].get(key, None)
#
# The following routines test the watcher functionality.
#
# watcher function for testing with a 'condition' for termination (see test_Callback below).
def add(watcher):
    if 'out' not in watcher: # Initialize watchers['out'] before first watcher updates it.
        watcher['out'] = [2]
    pt = watcher['out'][-1]+watcher['n'] # next number
    # watcher['condition'] is a finesse permitting a stopping 'condition' in the watcher to inactivate the watcher.
    if watcher['is_continue']: # If active, updates watcher['is_continue'].
        watcher['is_continue'] = watcher['condition'](pt) 
    if watcher['is_continue'] is None or watcher['is_continue']: # passive or active watcher
        watcher['out'].append(pt)
    return watcher['is_continue']
# The Append class has a function self.append() that triggers watchers to two functions.
#   Each watcher add(k2v) adds 'n' to the last element of watchers['out'].
#   Note that Append.append() could be a much more complicated subroutine.
class Append:
    def __init__(self, watcher):
        self.watcher = watcher
    def append(self):
        return self.watcher.execute()
    
def test_Callback():
    # Tests 3 watchers, all active.
    watchers = [
        {'callback':add,'n':3,'is_continue':True,'condition':(lambda pt: pt < 12)},
        {'callback':add,'n':5,'is_continue':True,'condition':(lambda pt: pt < 13)},
        {'callback':add,'n':1,'is_continue':True,'condition':(lambda pt: pt < 16)}
    ] 
    a = Append(Callback(watchers))
    while a.append(): 
        pass
    # Execution halst when every watcher['is_continue'] has become False.
    # Every watchers is then inactive.
    out = [[2, 5, 8, 11], [2, 7, 12], [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]]
    for j in range(3):
        assert a.watcher.out(j) == out[j]
    assert a.watcher.out() == out
    # Tests 3 watchers, 2 active and 1 passive.
    watchers = [
        {'callback':add,'n':3,'is_continue':True,'condition':(lambda pt: pt < 12)},
        {'callback':add,'n':5,'is_continue':True,'condition':(lambda pt: pt < 13)},
        {'callback':add,'n':1}
    ] 
    a = Append(Callback(watchers))
    # Execution halts when both watchers have returned False and the 2 active watchers are inactive.
    while a.append(): 
        pass
    # Execution halst when every watcher['is_continue'] has become falsy (False or None).
    # The passive watcher records one more event than any active watcher.
    out = [[2, 5, 8, 11], [2, 7, 12], [2, 3, 4, 5, 6]]
    assert a.watcher.out() == out
    # Tests 3 watchers, 1 active and 2 passive.
    #   Passive watchers do not have watcher['condition'] to make watcher['is_continue'] == False.
    watchers = [
        {'callback':add,'n':3,},
        {'callback':add,'n':5,'is_continue':True,'condition':(lambda pt: pt < 13)},
        {'callback':add,'n':1}
    ] 
    a = Append(Callback(watchers))
    # Execution halts when both watchers have returned False and the 2 active watchers are inactive.
    while a.append(): 
        pass
    # Execution halst when every watcher['is_continue'] has become falsy (False or None).
    # The 2 passive watchers records one more event than any active watcher.
    out = [[2, 5, 8, 11], [2, 7, 12], [2, 3, 4, 5]]
    assert a.watcher.out() == out
    
def main(): 
      test_Callback()
      
if __name__ == "__main__":
    main()
