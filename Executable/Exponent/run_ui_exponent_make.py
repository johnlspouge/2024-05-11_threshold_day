#!/usr/bin/env python

from os import system

log = f'run_ui_exponent.log'

# Explanations follow each option as a comment after '#', with (defaults) in parentheses. 
O = ' -o ../../Output/Exponent/exponent.csv' # -o output file
I = ' -i ../../Data/Exponent/exponent0.csv' # -i input file
D = '' # DEBUG flag ' -d '

system( f'python run_ui_exponent.py {O} {I} {D} > {log}' )
