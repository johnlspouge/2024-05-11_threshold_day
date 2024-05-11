#!/usr/bin/env python
"""
Calculates confidence region for threshold time from input file and parameters.
"""
import sys
sys.path.insert(0,"../modules")

from argparse import ArgumentParser, RawTextHelpFormatter
from os.path import isfile
import numpy as np
import pandas as pd
from scipy.stats import bootstrap
from collections import namedtuple

from jls_quantile_util import Quantile, counter2samples

def main(): 
    parser = getArguments()
    argument = parser.parse_args()
    df0 = check(argument) 
    # Reads input file.
    ts = sorted([ int(t) for t in df0.columns if t is not None and t.isnumeric() ])
    ts_c = list(map(str, ts))
    cols = [ col for col in df0.columns if col not in ts_c ]
    cols.pop(0)
    cols.insert(0, '')
    quantiles = []
    for q in argument.quantiles:
        quantiles.extend([f'{q} left', f'{q} sample', f'{q} right'])
    cols.extend(quantiles)
    #print(cols, flush=True)
    ofh = open(argument.ofn, 'w')
    ofh.write(','.join(map(str, cols)))
    ofh.write('\n')
    ofh.flush()
    index2row_dict = df0.to_dict(orient='index')
    for i,col2value in index2row_dict.items():
        print('Row', i, flush=True)
        row = index2row_dict[i]
        if argument.n_resamples == 0:
            n_resamples = int(row['realizations'])
        else:
            n_resamples = argument.n_resamples 
        counter = {int(t):int(row[t]) for t in ts_c}
        samples = (np.array(counter2samples(counter)),)        
        cendpts = []
        for q in argument.quantiles:
            print(f'Bootstrap {q}', flush=True)
            cendpts.extend(estimates(samples, q, n_resamples, argument.confidence))
        # Outputs the row.
        values = [ row[col] for col in cols if col != '' and col not in quantiles ]
        values.insert(0, i)
        values.extend(cendpts)
        ofh.write(','.join(map(str, values)))
        ofh.write('\n')
        ofh.flush()
        print('', flush=True)
    ofh.close()
# Returns the estimates as a list to extend cendpts.
def estimates(samples, q, n_resamples, confidence_level): # The quantile is a.
    Quantile.set(q)
    sample_a = Quantile.quantile(samples)
    try:
        res = bootstrap(samples, statistic=Quantile.quantile, n_resamples=n_resamples, confidence_level=confidence_level, method='basic')
        ci_a = res.confidence_interval
    except:
        # https://docs.python.org/3/library/collections.html#collections.namedtuple
        ci_a = namedtuple('Point', ['low', 'high'])
        ci_a.low = ci_a.high = sample_a
    return [ci_a.low, sample_a, ci_a.high]
# Check and fixes arguments if possible.    
def check( argument ):
    ifn = argument.ifn
    if not isfile( f'{ifn}' ):
        raise ValueError( f'Error: a valid INPUT_FILE_CSV "{ifn}" is required.' )
    df = pd.read_csv(ifn)
    integers = [ int(i) for i in df.columns if i is not None and i.isnumeric() ]
    if sorted(integers) != list(range(min(integers), max(integers)+1)):
        raise ValueError('df integer headers are not consecutive.')
    argument.df = df # Adds dataframe to argument to avoid re-reading it.
    for q in argument.quantiles:
        if not 0.0 <= q <= 1.0:
            raise ValueError( f'Error: not 0.0 <= q = {q} <= 1.0.' )
    if not 0.0 < argument.confidence < 1.0:
        raise ValueError( f'Error: not 0.0 < confidence = {argument.confidence} < 1.0.' )
    if argument.n_resamples < 0:
        raise ValueError( f'Error: n_resamples = {argument.n_resamples} < 0.' )
    return df
        
def getArguments():
    parser = ArgumentParser(description='Calculates confidence intervals probability mass functions.\n',
                            formatter_class=RawTextHelpFormatter)
    parser.add_argument("-i", "--input_file_csv", dest="ifn", type=str, default='../../Output/Quantile/pmfs.csv', 
                        help="INPUT_FILE_CSV contains the data with probability mass functions.", metavar="INPUT_FILE_CSV")
    parser.add_argument("-o", "--output_file_csv", dest="ofn", type=str, default='../../Output/Quantile/quantile_ci.csv', 
                        help="OUTPUT_FILE_CSV contains the data with the confidence regions.", metavar="OUTPUT_FILE_CSV")
    parser.add_argument("-q", "--quantiles", dest="quantiles", type=float, nargs='+', # sizes of the confidence interval 
                        help="QUANTILES estimates the quantile of the threshold day.", metavar="QUANTILES")
    parser.add_argument("-c", "--confidence", dest="confidence", type=float, default=0.95, # size of the confidence interval 
                        help="CONFIDENCE estimates the probability that the threshold day lies within the bootstrap confidence interval.", metavar="CONFIDENCE")
    parser.add_argument("-n", "--n_resamples", dest="n_resamples", type=int, default=10000, #(bootstrap resamples) (0: defaults to the sample size)
                        help="N_RESAMPLES counts the resamples in the bootstrap.", metavar="N_RESAMPLES")
    return parser
    
if __name__ == "__main__":
    main()
