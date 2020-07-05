#!/usr/bin/env python3
# Use GNU Parallel to sweep across parameters
import numpy as np
import subprocess
from distutils.util import strtobool
import os
import sys

jobfile = 'parallelJobs.txt'
# params = np.linspace(0,0.1, num=11)
# params= [1E-6, 1E-5, 1E-4, 1E-3, 1E-2]
params = np.linspace(1E-5,1E-3,num=20)

if not strtobool(input('\nRun jobs? (y/n): ')):
    print('Quitting...')
    sys.exit()

with open(jobfile, 'w') as runfile:
    for param in params:
        # runfile.write(f'python3 monteCarlo.py -f windowLoss_{param:.2f} -n 500000 -wl {param}\n')
        runfile.write(f'python3 monteCarlo.py -f lpb_{param:.2e}_noWL -n 500000 -lpb {param} -wl 0\n')

process = subprocess.Popen([f'parallel < {jobfile}'], shell=True)
print(f'Running {len(params)} jobs...')
