#!/usr/bin/env python

from os import system, listdir

dir = './'

# Runs the files jls*.py. 
for file in listdir(dir):
    if file.startswith("jls") and file.endswith(".py"):
        print(file)
        system( f'python {file}')
