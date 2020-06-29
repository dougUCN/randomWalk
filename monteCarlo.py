#!/usr/bin/env python3
import numpy as np
import tables
import os
import argparse
from collections import namedtuple

from randomWalkParticle import neutron1D
from tqdm import trange

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--n', type=int, required=True, help='Number of particles to simulate')
    parser.add_argument('-f', '--file', type=str, required=True, help='Stores end particle states to h5 file')
    args = parser.parse_args()

    ## parameters ##
    pipeID = 3 * 0.0254 # Inner diameter of the pipe [meters]
    cellEntranceID = 2* 0.0254 # [m]
    pipeL  = 12         # Total length of the pipe [meters]
    start  = 6          # starting position [meters]
    window = 9          # PPM window location
    nonspec = 0.05      # Chance for a nonspecular bounce
    lossPerBounce = 1E-4 # Loss per bounce (NiPh)
    windowLoss = 0.06   # Chance for loss when passing through window
    # windowLoss = 0   # Chance for loss when passing through window
    ################

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
    print('### Converted input parameters ###')
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
    Calculates mean free path between nonspecular bounces (For a cylinder = 4V/A * number of specular bounces)
    Converts starting position and window position from meters to integer multiples of this mfp
    D2 source defined as 0, cell entrance defined at the other end of the pipe
    Odds of getting into the cell = 1 -  [1 - (opening area)/(total inner surface of last segment)]^N_bounces
    Loss Per Bounce converted to loss per random walk step
    '''
    parameters = namedtuple('parameters',['cellChance','cell','source',
                            'gateValve', 'window','mfp','lossPerStep', 'windowLoss'])
    mfp = 1/nonspec * 4 * (pipeL * (pipeID/2)**2 * np.pi) / (pipeID * np.pi * pipeL + 2* (pipeID/2)**2 * np.pi)
    cell = np.ceil(pipeL/mfp)
    window = window/mfp
    gatevalve = np.ceil(start/mfp)
    cellChance = 2 * (1 - ( 1-(cellEntranceID/2)**2/(pipeID*mfp) ) ** (1/nonspec))
    lossPerStep = (1/nonspec) * lossPerBounce

    return parameters(cellChance=cellChance, cell=cell, window=window, mfp=mfp,
                        gateValve=gatevalve, source=0, lossPerStep=lossPerStep, windowLoss=windowLoss)

def noExt(filename, extensionName):
    '''Removes the extension name from a filename if present'''
    if filename[-len(extensionName):].find(extensionName) != -1:
        return filename[:-len(extensionName)]
    return filename


if ( __name__ == '__main__' ):
    main()
