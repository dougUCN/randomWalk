#!/usr/bin/env python3
import numpy as np
import tables
import os
import argparse
from collections import namedtuple

from randomWalkParticle import neutron1D
from tqdm import trange

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-n', '--n', type=int, required=True, help='Number of particles to simulate')
    parser.add_argument('-f', '--file', type=str, required=True, help='Stores end particle states to h5 file')
    parser.add_argument('-lpb', '--lossPerBounce', type=float, default=1E-4, help='Loss per bounce on pipe')
    parser.add_argument('-ns', '--nonSpec', type=float, default=0.05, help='Chance for nonspecular bounce')
    parser.add_argument('-wl', '--windowLoss', type=float, default=0.03, help='Chance of neutron loss for single window pass')
    args = parser.parse_args()

    ## parameters ##
    pipeID = 3 * 0.0254 # Inner diameter of the pipe [meters]
    cellEntranceID = 0.0745 # [m]
    pipeL  = 12         # Total length of the pipe [meters]
    start  = 6          # starting position [meters]
    window = 9.1          # PPM window location
    nonspec = args.nonSpec      # Chance for a nonspecular bounce
    lossPerBounce = args.lossPerBounce # Loss per bounce (NiPh)
    windowLoss = args.windowLoss   # Chance for loss when passing through window
    ################

    if nonspec<=0:
        print('ERROR: non specularity cannot be 0 for this calculation!!')
        return

    # Prep file IO
    outputFormat = {'particleNum':     tables.Int64Col(pos=0),
                    'windowHits':      tables.Int64Col(),
                    'totalSteps':      tables.Int64Col(),
                    'location':        tables.Int64Col(),
                    'status':          tables.StringCol(16),
                    'cellRejections':  tables.Int64Col()}
    filename = noExt(noExt(args.file, '.hdf'),'.h5') + '.h5'
    if os.path.exists(filename):
        print(f'{filename} already exists')
        return
    else:
        h5file = tables.open_file(filename, mode='w')
        group = h5file.create_group('/', 'branch')
        table = h5file.create_table(group, 'neutron1D', description=outputFormat, expectedrows=args.n)

    # Set up random walk 1D line
    params = calcEnvironment(pipeID, pipeL, cellEntranceID, start, window, nonspec, lossPerBounce, windowLoss)
    print('### Input parameters ###')
    for field in params._fields:
        print("{:<20}{:<20}".format(field, getattr(params, field)) )
    # Save params to file
    table.attrs.params = params._asdict()

    neutron = neutron1D(params)

    for i in trange(args.n):
        neutron.walk()
        endState = neutron.getState()
        for key in outputFormat.keys():
            table.row[key] = endState[key]
        table.row.append()
        neutron.resetState(label=i+1)

    # flush buffer and write to disc
    table.flush()
    return

def calcEnvironment(pipeID, pipeL, cellEntranceID, start, window, nonspec, lossPerBounce, windowLoss):
    '''
    Cleans up parameters to pass to monte carlo
    '''
    parameters = namedtuple('parameters',['cellChance','cell','source',
                            'gateValve', 'window','stepSize','lossPerStep', 'windowLoss', 'lossPerBounce'])
    mfp = pipeID * np.sqrt( 2*(2-nonspec)/nonspec/3 )
    cell = pipeL
    window = window
    gatevalve = start
    cellChance = (1 - ( 1-(cellEntranceID/2)**2/(pipeID*mfp) ) ** (1/nonspec) )**2
    lossPerStep = (1/nonspec) * lossPerBounce

    return parameters(cellChance=cellChance, cell=cell, window=window, stepSize=mfp,
                        gateValve=gatevalve, source=0, lossPerStep=lossPerStep, windowLoss=windowLoss, lossPerBounce=lossPerBounce)

def noExt(filename, extensionName):
    '''Removes the extension name from a filename if present'''
    if filename[-len(extensionName):].find(extensionName) != -1:
        return filename[:-len(extensionName)]
    return filename


if ( __name__ == '__main__' ):
    main()
