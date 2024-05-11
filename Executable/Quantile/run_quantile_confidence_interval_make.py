#!/usr/bin/env python

from os import system

log = f'run_quantile_confidence_interval.log'

I = ' -i ../../Output/Threshold/pmfs.csv'
O = ' -o ../../Output/Quantile/quantile_confidence_interval.csv'
Q = ' -q 0.025 0.05 0.5 0.95 0.975'
C = ' -c 0.95'
N = ' -n 10000'

system( f'python run_quantile_confidence_interval.py {I} {O} {Q} {C} {N} > {log}' )

# To examine, e.g., columns 1 and 4 of the output in UNIX:
# cat ../../Output/pmfs_top_20_slopes.csv | cut -d ',' -f1,4