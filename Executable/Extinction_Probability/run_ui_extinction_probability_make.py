#!/usr/bin/env python

from os import system

log = f'run_ui_extinction_probability.log'

# Explanations follow each option as a comment after '#', with (defaults) in parentheses. 
O = ' -o ../../Output/Extinction_Probability/q_exact.csv' # -o output file
I = ' -i ../../Data/Extinction_Probability/q_exact0.csv' # -i input file

system( f'python run_ui_extinction_probability.py {O} {I} > {log}' )
