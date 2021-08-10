#!/usr/bin/env python3
# Use GNU Parallel to sweep across parameters
import numpy as np
import subprocess
from distutils.util import strtobool
import os
import sys

jobfile = 'parallelJobs.txt'

params = [ 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2]

if not strtobool(input('\nRun jobs? (y/n): ')):
    print('Quitting...')
    sys.exit()

with open(jobfile, 'w') as runfile:
    for i, param in enumerate(params):

        runfile.write(f'./randomWalk_t.x --f cellExitMfpScan/mfp_{param:.2e}.h5 --n 1000000 --mfp2 {param} --progress false\n')
        runfile.write(f'./randomWalk_t.x --f cellExitMfpScan/mfp_{param:.2e}_noWL.h5 --n 1000000 --wl 0 --mfp2 {param} --progress false\n')

process = subprocess.Popen([f'parallel < {jobfile} --bar'], shell=True)
print(f'Running {len(params)} jobs...')
