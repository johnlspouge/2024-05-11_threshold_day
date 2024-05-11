#!/usr/bin/env python
"""
Calculates the exponent for gamma-distributed latent and infectious periods. 
"""
import argparse
from os.path import isfile, exists, dirname
from os import mkdir
from io import StringIO

import numpy as np
import pandas as pd

from math import log, isclose
from scipy.optimize import fsolve

# names of the relevant columns
COLS = ['e->i_mean',
        'e->i_dispersion',
        'i->r_mean',
        'i->r_dispersion',
        's->e:i_R_0']

def main():
    parser = getArguments()
    argument = parser.parse_args()
    check(argument) 
    DEBUG = argument.debug
    df = argument.df[COLS]
    lambdas = []
    doubling_times = []
    for index, row in df.iterrows(): 
        (e_mu, e_kappa, i_mu, i_kappa, r0) = row.to_list()
        def _laplace_exposed(theta, mu, kappa):
            laplace = (1.0+theta*mu/kappa)**(-kappa)
            return laplace
        def _laplace_infectious(theta, mu, kappa):
            if theta <= 0.0:
                return 1.0
            laplace = 1.0-(1.0+theta*mu/kappa)**(-(kappa-1.0))
            laplace /= kappa-1.0
            laplace /= theta*mu/kappa
            return laplace
        def laplace_generation(theta):
            return _laplace_exposed(theta, e_mu, e_kappa)*_laplace_infectious(theta, i_mu, i_kappa)-1.0/r0
        theta0 = 0.1
        theta = fsolve(laplace_generation, theta0) # exponent lambda
        if DEBUG:
            assert np.allclose((e_mu, e_kappa, i_mu, i_kappa, r0), (e_mu, e_kappa, i_mu, i_kappa, r0))
            assert isclose(theta, 0.14156035, abs_tol=1.0e-06)
            assert isclose(_laplace_exposed(theta, e_mu, e_kappa), 0.62682015, abs_tol=1.0e-06)
            assert isclose(_laplace_infectious(theta, i_mu, i_kappa), 0.79767698, abs_tol=1.0e-06)
            assert isclose(1.0/r0, 0.5, abs_tol=1.0e-06)
            assert isclose(laplace_generation(theta), 0.0, abs_tol=1.0e-09)
        lambdas.append(theta)
        doubling_times.append(log(2.0)/theta)
    df['lambda'] = np.array(lambdas) # Adds values only not arrays of dimension 1 with np.array,
    df['doubling_time'] = np.array(doubling_times) # Adds values only not arrays of dimension 1 with np.array,
    if not DEBUG:
        df.to_csv(argument.ofn, index=False)
# Reads the string from test data, with columns COLS.
def to_test(argument):
    data = [[3.5, 4, 5.5, 0.3, 2]]
    df = pd.DataFrame(data, columns=COLS)
    return df
# Reads the string from argument.ifn, a CSV DataFrame with columns COLS.
def to_df(string):
    ifh = StringIO(string)
    df = pd.read_csv(ifh, sep=',', dtype={COLS[0]:float,COLS[1]:float,COLS[2]:float,COLS[3]:float,COLS[4]:float})
    arr = df.to_numpy()
    if (np.abs(arr) <= 0.0).any():
        raise ValueError('The dataframe elements must be positive.')
    return df
# Checks arguments.    
def check(argument): 
    if argument.debug:
        argument.df = to_test(argument)
        return
    # argument.odir is the output directory.
    odir = dirname(argument.ofn)
    if not exists(odir):
        mkdir(odir)
    # argument.ifn contains a DataFrames with columns 
    #    'e->i_mean' : E gamma distribution mean
    #    'e->i_dispersion' : E gamma distribution dispersion
    #    'i->r_mean' : I gamma distribution mean
    #    'i->r_dispersion' : I gamma distribution mean
    #    's->e:i_R_0' : I mean basic reproduction number
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
    parser = argparse.ArgumentParser(description='Calculates the exponential growth lambda for SEIR model, where E and I are gamma-distributed.\n')
    parser.add_argument("-o", "--ofn", dest="ofn", type=str, default='../../Output/Exponent/exponent.csv', 
                        help="OFN contains the output DataFrame with the exponential growth lambda.", metavar="OFN")
    parser.add_argument("-i", "--ifn", dest="ifn", type=str, required=True,  
                        help="IFN is a CSV with DataFrame whose columns define the parameters of the gamma distributions.", metavar="IFN")
    parser.add_argument("-d", "--debug", dest="debug", default=False, action='store_true', 
                        help="flag for debugging.")
    return parser
    
if __name__ == "__main__":
    main()
