#!/usr/bin/env python

from os import system

log = f'run_ui_threshold.log'

# Explanations follow each option as a comment after '#', with (defaults) in parentheses. 
O = ' -o ../../Output/Threshold/Final' # -o output directory
F = ' -f ../../Data/Threshold/0_exposed_123.csv'
D = ' -d ../../Data/Threshold/1_add_decay_to_' # prefix of all files used
I = ' -i ../../Data/Threshold/2_add_infect_to_' # prefix of all files used
W = ' -w ../../Data/Threshold/3_add_watchers_1.csv'
R = ' -r 10000' # -r 10000 REALIZATION is the number of realizations in the simulation.
P = ' -p' # ' -p' # flag for output from passive watchers.

system( f'python run_ui_threshold.py {O} {F} {D} {I} {W} {R} {P} > {log}' )

