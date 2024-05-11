#!/usr/bin/env python
"""
Runs simulations for the first day where the new cases exceeded a threshold. 
    For reproducible Monte Carlo, e.g., the program sets np.random.seed(seed=argument.seed).
    See the subroutine check for restrictions on the input.
"""
import sys
sys.path.insert(0,"../modules")

import argparse
from os.path import isfile, exists, split
from os import mkdir, listdir

import numpy as np

from jls_epidemic_io import Io
from jls_epidemic_simulation import product_simulate

def main():
    parser = getArguments()
    argument = parser.parse_args()
    check(argument)  
    test_epidemic(argument) 

def test_epidemic(argument):
    # seed for the simulations
    np.random.seed(seed=argument.seed)
    # a list of lists of pd.DataFrame-s
    dfs_epidemic = [argument.dfs_initialize_from, 
                    argument.dfs_add_decay_to,
                    argument.dfs_add_infect_to]
    product_simulate(dfs_epidemic,
            argument.df_watchers,
            argument.realization_count,
            argument.odir, argument.is_passive_output)
# Returns all files starting with the prefix = folder+ifn_prefix, lexicographically sorted for reproducibility.    
def ifns_with(prefix):
    folder, ifn_prefix = split(prefix)
    if not exists(folder):
        raise ValueError('Error: "{folder}" does not exist.')
    ifns = sorted([folder+'/'+ifn for ifn in listdir(folder) if ifn.startswith(ifn_prefix)])
    if not ifns:
        raise ValueError('Error: "No files in {folder}" have the prefix "{prefix}".')
    return ifns
# Checks arguments.    
def check(argument): 
    # argument.odir is the output directory.
    if not argument.odir.endswith('/'):
        argument.odir += '/'
    if not exists(argument.odir):
        mkdir(argument.odir)
    # argument.realization_count contains the number of successful realizations in each epidemic simulation.
    if not 0 < argument.realization_count:
        raise ValueError(f'Error : not 0 < {argument.realization_count} = realization_count')
    # argument.seed contains the random seed, for reproducibility of simulations.
    if not 0 < argument.seed:
        raise ValueError(f'Error : not 0 < {argument.seed} = random number seed')
    # argument.dfs_initialize_from contains a list of DataFrames with the epidemic compartments and set(s) of initial counts.
    #   argument.initialize_from is a CSV file with compartments and column(s) of initial counts. 
    if not isfile(argument.initialize_from):
        raise ValueError(f'Input file "{argument.initialize_from}" does not exist.')
    argument.dfs_initialize_from = Io.dfs_initialize_from(argument.initialize_from)
    compartments = list(argument.dfs_initialize_from[0].columns) # compartments common to all simulated epidemics
    # argument.dfs_add_decay_to contains DataFrames from every file that starts.with(argument.add_decay_to).
    ifns = ifns_with(argument.add_decay_to)
    dfs = [Io.df_add_decay_to(compartments, ifn) for ifn in ifns] # a list of DataFrames containing the decays
    # argument.dfs_add_decay_to contains a list of DataFrames describing decays from one compartment to another.
    #   The columns ['from','to'] in the DataFrames must be identical, so decay topology is consistent across the DataFrames.
    df0 = dfs[0]
    for df in dfs:
        if not df[['from','to']].equals(df0[['from','to']]): # Decay topology is consistent across the DataFrames.
            raise ValueError(f'Error: Columns ["from","to"] are not equal in df_add_decay_to.\n{df}\n{dfs[0]}')
    argument.dfs_add_decay_to = dfs
    decay_froms = list((dfs[0])['from']) # Maintains a list of unstable compartments, because every infected compartment must decay.
    # argument.dfs_add_infect_to contains a list of DataFrames with the possible infections.
    #   add_infect_to contains all files starts.with(argument.add_infect_to).
    #   The columns ['from','to','infecting'] are identical, so infection topology is consistent across the files.
    ifns = ifns_with(argument.add_infect_to)
    dfs = [Io.df_add_infect_to(compartments, decay_froms, ifn) for ifn in ifns] 
    # argument.dfs_add_infect_to contains a list of DataFrames describing infections.
    #   The columns ['from','to','infecting'] in the DataFrames must be identical, so infection topology is consistent across the DataFrames.
    df0 = dfs[0]
    for df in dfs:
        if not df[['from','to','infecting']].equals(df0[['from','to','infecting']]): # Infection topology is consistent across the DataFrames.
            raise ValueError(f'Error: Columns ["from","to","infecting"] are not equal in df_add_decay_to.\n{df}\n{dfs[0]}')
    argument.dfs_add_infect_to = dfs
    # argument.df_watchers contains the watchers.
    if not isfile(argument.add_watchers):
        raise ValueError(f'Input file "{argument.add_watchers}" does not exist.')
    argument.df_watchers = Io.df_watchers(compartments, argument.add_watchers)
    if len(set(argument.df_watchers['outfile'])) != len(argument.df_watchers['outfile']): # The output for the watchers goes to distinct files.
        raise ValueError('Each active watcher must have its own outfile.')
        
def getArguments():
    parser = argparse.ArgumentParser(description='Runs continuous stochastic epidemic to determine the exceedance time for cases over cutoff.\n')
    parser.add_argument("-o", "--output_dir", dest="odir", type=str, default='../../Output/Threshold/', 
                        help="OUTPUT_DIR contains the output from the watchers, default='../../Output/Threshold/'.", metavar="OUTPUT_DIR")
    parser.add_argument("-r", "--realization_count", dest="realization_count", type=int, 
                        help="REALIZATION_COUNT is the number of realizations in the simulation.", metavar="REALIZATION_COUNT")
    parser.add_argument("-s", "--seed", dest="seed", type=int, default=31415, 
                        help="SEED is the random number seed, default=31415.", metavar="SEED")
    parser.add_argument("-f", "--initialize_from", dest="initialize_from", type=str, required=True,
                        help="INITIALIZE_FROM gives the CSV file of the DataFrame of compartment counts for initialize_from.", metavar="INITIALIZE_FROM")
    parser.add_argument("-d", "--add_decay_to", dest="add_decay_to", type=str, required=True, 
                        help="ADD_DECAY_TO gives the CSV file prefix of (possibly several) DataFrames for add_decay_to.", metavar="ADD_DECAY_TO")
    parser.add_argument("-i", "--add_infect_to", dest="add_infect_to", type=str, required=True,  
                        help="ADD_INFECT_TO gives the CSV file prefix of (possibly several) DataFrames for add_infect_to.", metavar="ADD_INFECT_TO")
    parser.add_argument("-w", "--add_watchers", dest="add_watchers", type=str, default=None,  
                        help="ADD_WATCHERS gives the CSV file of DataFrames for to_epidemic_with_watchers.", metavar="ADD_WATCHERS")
    parser.add_argument("-p", "--is_passive_output", dest="is_passive_output", default=False, action='store_true',  
                        help="IS_PASSIVE_OUTPUT outputs passive watchers for all events and for integer times.")
    return parser
    
if __name__ == "__main__":
    main()
