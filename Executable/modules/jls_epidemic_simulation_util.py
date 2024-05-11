#!/usr/bin/env python
"""
Utility routines for epidemic simulation. 
"""

# Converts the dictionary 'out' in a passive watcher to CSV for a DataFrame. 
#   Assumes the dictionary 'out' has a key 't' for time.
def to_csv(out):
    if 't' not in out:
        raise ValueError(f'Error : The dictionary is not the output from an epidemic simulation.')
    if not all(isinstance(x, (int,float)) for x in out['t']):
        raise ValueError(f'Error : Not all elements of the list out["t"] have type float.')
    for k,vs in out.items():
        if k == 't':
            continue
        if len(vs) != len(out['t']):
            raise ValueError(f'Error : The lists in out have different lengths len("{k}") = {len(k)} != {len(out["t"])} = len(out["t"]).')
        if not all(isinstance(x, int) for x in vs):
            raise ValueError(f'Error : Not all elements of the list "{k}" have type int.')
    cols = [k for k in out.keys() if k != 't']
    cols.insert(0, 't')
    s = ','.join(cols)+'\n'
    lists = [out[c] for c in cols]
    tuples = list(zip(*lists))
    for i,t in enumerate(tuples):
        tuples[i] = list(map(str,t))
    lines = [','.join(t) for t in tuples]
    s = s+'\n'.join(lines)
    return s
# Returns the integer times for a watcher.
def ts(i):
    return float(i)
# Returns the first day of exceedance of eta.
#   Returns None if no exceedance
def day_exceeds(eta, ts, i0, r0):
    if ts != list(range(len(ts))) and ts != tuple(range(len(ts))):
        raise ValueError(f'ts != list(range(len(ts))) and ts != tuple(range(len(ts)))')
    if len(ts) != len(i0):
        raise ValueError(f'len(ts) = {len(ts)} != {len(i0)} = len(i0)')
    if len(ts) != len(r0):
        raise ValueError(f'len(ts) = {len(ts)} != {len(r0)} = len(r0)')
    day = 0
    if len(ts) <= day+1: # no more days with new cases
        return None # is_continue
    cumulative_cases0 = i0[day]+r0[day]
    while day+1 < len(ts):
        day += 1
        cumulative_cases = i0[day]+r0[day]
        cases = cumulative_cases-cumulative_cases0
        cumulative_cases0 = cumulative_cases
        if cases > eta:
            return day # The day contains the day of exceedance.
    return None
    
def main(): 
    assert day_exceeds(3, (0,1,2,3,4,5), (1,1,3,2,1,4), (0,0,2,2,1,4)) == 2
    assert day_exceeds(3, (0,1,2,3,4,5), (1,1,3,2,1,4), (0,0,1,1,0,1)) == 5
    assert day_exceeds(3, (0,1,2,3,4,5), (1,1,3,2,1,4), (0,0,1,1,1,0)) is None
    out = {'t':[0,1,1.2,1.3],'s':[10,10,8,8],'i':[1,1,2,1],'r':[0,0,1,2]}
    csv0 = """\
t,s,i,r
0,10,1,0
1,10,1,0
1.2,8,2,1
1.3,8,1,2\
"""
    assert csv0 == to_csv(out)
      
if __name__ == "__main__":
    main()
