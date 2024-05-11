#!/usr/bin/env python
"""
Calculates the extinction probability for (sufficiently simple) branching processes. 
    See the subroutine check for restrictions on the input.
"""
import sys
sys.path.insert(0,"../modules")

import argparse
from os.path import isfile, exists, dirname
from os import mkdir
from io import StringIO

import pandas as pd

from jls_branching_process import Branching_Process_Factory

def main():
    parser = getArguments()
    argument = parser.parse_args()
    check(argument)  
    df = argument.df
    (qs, gammas, lineage_doomeds) = ([], [], [])
    for index, row in df.iterrows(): 
        n = row['ancestors']
        r0 = row['R0']
        dispersion = row['dispersion']
        if pd.isna(dispersion):
            bp = Branching_Process_Factory(r0)
        else:
            bp = Branching_Process_Factory(r0, dispersion)
        qs.append(bp.q()**n)
        gammas.append(bp.gamma())
        lineage_doomeds.append(1.0/(1.0-bp.gamma()))
    df['q_exact'] = qs
    df['gamma'] = gammas
    df['lineage_doomed'] = lineage_doomeds
    df.to_csv(argument.ofn, index=False)

# Reads the string from argument.ifn, a DataFrames with columns 'ancestors', 'R0', and 'dispersion'.
def to_df(string):
    ANCESTORS = 1
    #DISPERSION = None
    ifh = StringIO(string)
    df = pd.read_csv(ifh, sep=',', dtype={'ancestors':int,'R0':float,'dispersion':float})
    df['ancestors'] = df['ancestors'].fillna(ANCESTORS)
    #df['dispersion'] = df['dispersion'].fillna(DISPERSION)
    #df['ancestors'] = df['ancestors'].astype(int)
    #df['r0'] = df['r0'].astype(float)
    #df['dispersion'] = df['dispersion'].astype(float)
    for c in df['ancestors']:
        if not isinstance(c, int) or c <= 0:
            raise ValueError(f'Error: Each entry {c} in the column "ancestors" must be a positive integer.')
    for c in df['R0']:
        if not isinstance(c, float) or c <= 0.0:
            raise ValueError(f'Error: Each entry {c} in the column "R0" must be a positive float.')
    for c in df['dispersion']:
        if pd.isna(c):
            c = None
        if c is not None and not (isinstance(c, float) and c > 0.0):
            raise ValueError(f'Error: Each entry in the column "dispersion" must be a positive float or None.')
    return df
# Checks arguments.    
def check(argument): 
    # argument.odir is the output directory.
    odir = dirname(argument.ofn)
    if not exists(odir):
        mkdir(odir)
    # argument.ifn contains a DataFrames with columns 'ancestors', 'R0', and 'dispersion'.
    #   The 'ancestors' column count the ancestors of the population (default 1).
    #   The 'R0', and 'dispersion' columns define a Negative Binomial or Poisson offspring distribution. 
    if not isfile(argument.ifn):
        raise ValueError(f'Input file "{argument.ifn}" does not exist.')
    try:
        with open(argument.ifn, 'r') as ifh:
            string = ifh.read()
    except:
        raise ValueError(f'The read of the input file "{argument.ifn}" failed.')
    try:
        argument.df = to_df(string)
    except:
        raise ValueError(f'The input file "{argument.ifn}" contained bad values.')
        
def getArguments():
    parser = argparse.ArgumentParser(description='Calculates the extinction probability for (sufficiently simple) single-type Galton-Watson processes.\n')
    parser.add_argument("-o", "--ofn", dest="ofn", type=str, default='../../Output/Extinction_Probability/q_exact.csv', 
                        help="OFN contains the output DataFrame with the extinction probability.", metavar="OFN")
    parser.add_argument("-i", "--ifn", dest="ifn", type=str, required=True,  
                        help="IFN is a CSV with DataFrame whose columns define Negative Binomial or Poisson offspring distributions.", metavar="IFN")
    return parser
    
if __name__ == "__main__":
    main()
