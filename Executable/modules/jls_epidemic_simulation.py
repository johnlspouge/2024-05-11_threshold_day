#!/usr/bin/env python
"""
Runs simulations. 
"""
import sys
sys.path.insert(0,"../modules")

import time
from os.path import exists
from os import mkdir
from itertools import product

from jls_epidemic_simulation_active_out import add_to, to_df_out_columns, to_counts
from jls_epidemic_simulation_util import to_csv
from jls_epidemic_io import Io
from jls_epidemic_watcher_threshold import is_passive, to_passive_watchers

CSV = '.csv'
PASSIVE = 'Passive/'
# Runs different simulations with a common epidemic topology but with different rates on the edges.
def product_simulate(
        dfs_epidemic, # [dfs_initialize_from, dfs_add_decay_to, dfs_add_infect_to]
        df_watchers, # all watchers, which will be preserved with deepcopy
        realization_count, # successful realizations (i.e., simulated epidemics not going extinct)
        odir, # output directory
        is_passive_output): # ? output the passive watchers ?
    # Creates passive output directory, if necessary.
    if is_passive_output and not exists(odir+PASSIVE):
        mkdir(odir+PASSIVE)
    # Organizes the output. 
    #   Assumes the driver program (e.g., run_ui_threshold_counts.py) has checked dfs_epidemics.
    #   All infectious networks in dfs_epidemic must be the same, although rates may differ.
    row_indexes = list(range(len(df_watchers)))
    if row_indexes:
        outfiles = list(df_watchers['outfile']) 
        df_epidemic = [df[0] for df in dfs_epidemic]
        out_dfs = [to_df_out_columns(df_epidemic, df_watchers) for out in outfiles] 
    # Lists all possible epidemic parameter combinations.
    dfs_initialize_from, dfs_add_decay_to, dfs_add_infect_to = dfs_epidemic
    dfs_epidemics = product(dfs_initialize_from, 
                            dfs_add_decay_to, 
                            dfs_add_infect_to)
    watchers = Io.to_watchers(df_watchers) # Initializes fresh watchers for each pass through simulate.
    print(watchers, flush=True)
    for dummy,df_epidemic in enumerate(dfs_epidemics): # Avoids memory issues with enumerate.
        print('epidemic [', dummy, ']')
        # Prepares final output from the watchers.
        d = simulate(df_epidemic, df_watchers, realization_count, odir, is_passive_output)
        assert len(df_watchers) == len(d['all_realizations'])
        # Accumulates output for the active watchers.
        for out in outfiles:
            df = df_watchers.loc[df_watchers['outfile'] == out]
            for index_watcher in list(df.index):
                out_dfs[index_watcher] = add_to(out_dfs[index_watcher], # old output DataFrame
                    df_epidemic, df_watchers, index_watcher, # epidemic parameters and watcher
                    d['all_realizations'][index_watcher], realization_count, to_counts(d['timess'][index_watcher]))
    # Completes the output for the active watchers.
    for i,out in enumerate(outfiles): 
        out_dfs[i].to_csv(odir+out+CSV, index=False)
# Returns stopping time distribution
def simulate(
        df_epidemic, # [df_initialize_from, df_add_decay_to, df_add_infect_to]
        df_watchers, # all watchers
        realization_count, # realizations in the simulation
        odir, # output directory
        is_passive_output): # ? output the passive watchers ?
    df_initialize_from, df_add_decay_to, df_add_infect_to = df_epidemic
    # constants of the computation
    DELTA = None # continuous-time simulation
    NU = 1 # full infectivity, no non-pharmaceutical interventions in place
    # Logs the simulation parameters.
    print(df_initialize_from, flush=True)
    print(df_add_decay_to, flush=True)
    print(df_add_infect_to, flush=True)
    watchers = Io.to_watchers(df_watchers) # Initializes fresh watchers for each pass through simulate.
    # Counts active watchers for output.
    len_active_watchers = len([w for w in watchers if not is_passive(w)])
    timess = [[] for i in range(len_active_watchers)] # stopping times in the watchers' realizations
    successful_realizations = [0]*len_active_watchers # count of successful realizations for watchers[i]
    all_realizations = [0]*len_active_watchers # total count of realizations, successful and failed, for watchers[i]
    # Simulates the epidemic.
    all_rs = 0 # Counts the total realizations.
    tic = time.perf_counter()
    # Simulates realizations until all active watchers have realization_count realizations.
    while any([r < realization_count for r in successful_realizations]):
        # Sets the parameters for epidemic simulation.
        Io.initialize_from(df_initialize_from, DELTA=DELTA) # experiment later to remove
        Io.add_decay_to(df_add_decay_to) # experiment later to remove
        Io.add_infect_to(df_add_infect_to, NU=NU)
        # Sets all watchers and partitions them into passive and active.
        watchers = Io.to_watchers(df_watchers) # Initializes fresh watchers for each pass through simulate.        
        passive_watchers = watchers[:len(to_passive_watchers())] # passive watchers 
        active_watchers = watchers[len(to_passive_watchers()):] # active watchers, defaults to [] if empty.
        # Simulates epidemic.
        epidemic = Io.to_epidemic_with_watchers(watchers)
        epidemic.evolve()
        # All active watchers are now inactive.
        # Stores the results from active watchers.
        for i,w in enumerate(active_watchers):
            if successful_realizations[i] == realization_count: # The i-th watcher has reached the required realization_count.
                continue
            all_realizations[i] += 1 # count of epidemic realizations, successful and unsuccessful
            if w['out'] is not None:
                timess[i].append(w['out']['t'][0])
                successful_realizations[i] += 1
        # Outputs the realizations from the passive watchers, if desired.
        if is_passive_output:
            for w in passive_watchers:
                out(w, all_rs, odir) 
        all_rs += 1
    assert all([r == realization_count for r in successful_realizations])
    toc = time.perf_counter()
    print(f'The simulation for {realization_count} successful realizations took {toc - tic:0.4f} seconds.', flush=True)
    # Returns the results from the active watchers.
    #   timess[0:#(active_watchers))[0:REALIZATION) ; all_realizations[0:#(active_watchers)
    d = {'all_realizations':all_realizations,'timess':timess}
    assert len(d['all_realizations']) == len(d['timess'])
    return d
# Outputs an entire passive watcher w at the end of a realization_count.
def out(w, all_rs, odir, passive=PASSIVE): # The argument all_rs differentiates the realizations.
    int_len = 5 # integer string length for left padding with 0s
    name = w['name']
    if name is None:
        raise ValueError('Error : w is not a passive watcher.')
    csv = to_csv(w['out'])
    fmt = f'0{int_len}d'
    txt = '{0:{1}}'.format(all_rs, fmt)
    ofn = odir+passive+name+'_'+txt+CSV
    with open(ofn, "w") as ofh:
        ofh.write(csv)

def main():
    pass
    
if __name__ == "__main__":
    main()
