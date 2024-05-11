#!/usr/bin/env python
"""
Contains auxiliary subroutines for output. 
"""
import sys
sys.path.insert(0,"../modules")

from math import isclose
import pandas as pd
import numpy as np
from collections import Counter

from jls_epidemic_io import Io

def _dfs_for_tests():
    # initialize_from
    compartment_names = ['s', 'e', 'i', 'r', 'i0', 'r0']
    counts = [100,1,0,0,0,0]
    df_initialize_from = pd.DataFrame([counts], columns=compartment_names)
    # 
    cols = ['from','to','mean','dispersion','prob']
    lst = [['e','i',2.5,4.0,0.25],['e','i0',2.5,4.0,0.75],['i','r',5.5,0.3,1.0],['i0','r0',5.5,0.3,1.0]]
    df_add_decay_to = pd.DataFrame(lst, columns=cols)
    #
    cols = ['from','to','infecting','R_0']
    lst = [['s','e','i',2.0],['s','e','i0',1.5]] # i0 is 0.75 as infectious as i
    df_add_infect_to = pd.DataFrame(lst, columns=cols)
    #
    cols = ['passive_watcher','is_difference','compartments','eta_threshold','outfile']
    lsts = [
        ['all',0,'i+r',29,'pmf'],
        ['integer',0,'i+r',39,'pmf'],
        ['integer',1,'i+r',49,'pmf']
    ]
    df_watchers = pd.DataFrame(lsts, columns=cols) # The df has a single row
    #
    df_epidemic = [df_initialize_from, df_add_decay_to, df_add_infect_to]
    return df_epidemic, df_watchers

# Returns the columns of an output DataFrame for the active watchers.    
#   Assumes that the df_epidemic have been checked, as follows.
#   All inputs to the simulation have the same columns, so the infectious network is the same.
#   df_epidemic is a list of pandas DataFrame-s 
#       [df_initialize_from, df_add_decay_to, df_add_infect_to].
#   Inputs to the simulation may have different parameters in the columns.
def to_df_out_columns(df_epidemic, df_watchers):
    # columns for initial counts
    dfs = df_epidemic
    compartment_names = list(dfs[0].columns) 
    # columns for decays
    decays = []
    for i in range(len(dfs[1])):
        decays.extend(list([dfs[1].loc[i,'from']+'->'+dfs[1].loc[i,'to']+'_'+param for param in Io.DECAY_PARAMS])) 
    # columns for infects
    infects = []
    for i in range(len(dfs[2])):
        infects.extend(list([dfs[2].loc[i,'from']+'->'+dfs[2].loc[i,'to']+':'+dfs[2].loc[i,'infecting']+'_'+param for param in Io.INFECT_PARAMS])) 
    # columns for the watcher
    watchers = list(df_watchers.columns[:-1]) # ['passive_watcher', 'is_difference', 'compartments',	'eta_threshold']
    cols = compartment_names+decays+infects+watchers    
    cols.append('q_estimated')
    cols.append('realizations')
    cols.append('0')
    return pd.DataFrame(columns=cols)

def test_to_df_out_columns():
    df_epidemic, df_watchers = _dfs_for_tests()
    cols0 = ['s', 'e', 'i', 'r', 'i0', 'r0', 'e->i_mean', 'e->i_dispersion', 'e->i_prob', 'e->i0_mean', 'e->i0_dispersion', 'e->i0_prob', 'i->r_mean', 'i->r_dispersion', 'i->r_prob', 'i0->r0_mean', 'i0->r0_dispersion', 'i0->r0_prob', 's->e:i_R_0', 's->e:i0_R_0', 'passive_watcher', 'is_difference', 'compartments', 'eta_threshold', 'q_estimated', 'realizations', '0']
    assert list(to_df_out_columns(df_epidemic, df_watchers)) == cols0
# Returns epidemic parameters for part of one row of the output DataFrame.
def to_epidemic_params(df_epidemic):
    # columns for initial counts
    dfs = df_epidemic
    compartment_counts = list(dfs[0].iloc[0]) 
    # columns for decays
    decays = []
    for row_index in range(len(dfs[1])):
        decays.extend(list(dfs[1].loc[row_index,Io.DECAY_PARAMS]))
    # columns for infects
    infects = []
    for row_index in range(len(dfs[2])):
        infects.extend(list(dfs[2].loc[row_index,Io.INFECT_PARAMS]))
    epidemic_params = compartment_counts+decays+infects
    return epidemic_params

def test_to_epidemic_params():
    df_epidemic, df_watchers = _dfs_for_tests()
    epidemic_params0 = [
            100, 1, 0, 0, 0, 0, # counts
            2.5, 4.0, 0.25, # e->i
            2.5, 4.0, 0.75, # e->i0
            5.5, 0.3, 1.0, # i->r
            5.5, 0.3, 1.0, # i0->r0
            2.0, 1.5] # s->e:i s->e:i0
    epidemic_params = to_epidemic_params(df_epidemic)
    assert epidemic_params == epidemic_params0
# Returns one set of watcher parameters for part of one row of the output DataFrame.
def to_watcher_params(df_watchers, row_index):
    watcher_params = []
    watcher_params.append(df_watchers.loc[row_index,'passive_watcher'])
    watcher_params.append(int(df_watchers.loc[row_index,'is_difference']))
    watcher_params.append(df_watchers.loc[row_index,'compartments'])
    watcher_params.append(int(df_watchers.loc[row_index,'eta_threshold']))
    return watcher_params

def test_to_watcher_params():
    df_epidemic, df_watchers = _dfs_for_tests()
    watcher_params0 = ['all', 0, 'i+r',	29]
    assert to_watcher_params(df_watchers, 0) == watcher_params0
    watcher_params0 = ['integer', 0, 'i+r',	39]
    assert to_watcher_params(df_watchers, 1) == watcher_params0
    watcher_params0 = ['integer', 1, 'i+r',	49]
    assert to_watcher_params(df_watchers, 2) == watcher_params0
# Returns one empirical extinction probability and the count of successful realizations for part of one row of the output DataFrame.
def to_success_params(all_realization, successful_realization):
    q = 1.0-successful_realization/all_realization
    return [q, successful_realization]
# Returns a list counting the 'integer' stopping times in integers (e.g., days). 
def to_counts(times): # timess[0:#(active_watchers))[0:REALIZATION) ; times = timess[i] for some i in [0:REALIZATION)
    assert all([isinstance(time, int) for time in times])
    counter = dict(Counter(times))
    max_time = max(counter.keys())
    counts = []
    for t in range(max_time+1):
        counts.append(counter.get(t,0)) # Outputs {days, counts}.
    assert sum(counts) == len(times)
    return counts

def test_to_counts():
    times = [0,3,2,2,1,1,1]
    counts0 = [1,3,2,1]
    assert to_counts(times) == counts0
# Adjusts the old df to the current counts, or the current counts to the old df.
def _adjust(counts, df): # The counts have max_count = len(counts)-1.
    max_pmf = int(df.columns[-1])
    max_count = len(counts)-1
    if max_pmf < max_count: # The DataFrame must be padded to the right with 0s.
        cols = [str(i) for i in range(max_pmf+1, max_count+1)]
        if len(df):
            df0 = pd.DataFrame(0, index=np.arange(len(df)), columns=cols) # new columns
        else:
            df0 = pd.DataFrame(columns=cols)
        df = pd.concat([df,df0],axis=1)
    elif max_count < max_pmf:
        counts += [0 for i in range(max_pmf-max_count)]
    return df, counts
# Adds row to the old df.
def add_to(df_out, # old output DataFrame
        df_epidemic, df_watchers, index_watcher, # epidemic parameters and watcher
        all_realization, successful_realization, counts): # counters for realizations, all and successful, current counts 
    row = to_epidemic_params(df_epidemic)
    row += to_watcher_params(df_watchers, index_watcher)
    row += to_success_params(all_realization, successful_realization)
    df_out, counts = _adjust(counts, df_out)
    row += counts
    df_out.loc[len(df_out.index)] = row
    return df_out

def test_add():
    df_epidemic, df_watchers = _dfs_for_tests()
    df_out = pd.DataFrame(columns=list(to_df_out_columns(df_epidemic, df_watchers))) 
    # row 0 of data
    times = [0,3,2,2,1,1,1]
    counts = to_counts(times)
    assert counts == [1,3,2,1]
    index_watcher = 0
    all_realization = 20
    successful_realization = 10
    df_out = add_to(df_out, 
                 df_epidemic, df_watchers, index_watcher, 
                 all_realization, successful_realization, counts)
    assert isclose(df_out.loc[0,'q_estimated'], 0.5)
    assert df_out.loc[0,'realizations'] == 10
    assert list(df_out.loc[0, to_strings(4)]) == [1,3,2,1]
    # row 1 of data
    times = [0,2,2,2,1,1,1]
    counts = to_counts(times)
    assert counts == [1,3,3]
    index_watcher = 0
    all_realization = 20
    successful_realization = 9
    df_out = add_to(df_out, 
                 df_epidemic, df_watchers, index_watcher, 
                 all_realization, successful_realization, counts)
    assert isclose(df_out.loc[1,'q_estimated'], 0.55)
    assert df_out.loc[1,'realizations'] == 9
    assert list(df_out.loc[1, to_strings(4)]) == [1,3,3,0]
    # row 2 of data
    times = [0,2,3,2,1,4,1,4]
    counts = to_counts(times)
    assert counts == [1,2,2,1,2]
    index_watcher = 0
    all_realization = 20
    successful_realization = 8
    df_out = add_to(df_out, 
                 df_epidemic, df_watchers, index_watcher, 
                 all_realization, successful_realization, counts)
    assert isclose(df_out.loc[2,'q_estimated'], 0.6)
    assert df_out.loc[2,'realizations'] == 8
    assert list(df_out.loc[2, to_strings(5)]) == [1,2,2,1,2]
# Converts range to list of strings.
def to_strings(*range_args): # arguments for range function
    return list(map(str, list(range(*range_args))))

def test_to_strings():
    assert to_strings(4) == ['0','1','2','3']
    assert to_strings(1,4) == ['1','2','3']
    assert to_strings(1,4,2) == ['1','3']

def main():
    test_to_df_out_columns()
    test_to_epidemic_params()
    test_to_watcher_params()
    test_to_counts()
    test_add()
    test_to_strings()
    
if __name__ == "__main__":
    main()
