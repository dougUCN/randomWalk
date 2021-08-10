import numpy as np
import sys

class neutron1D:
    maxSteps = 1000 # Stops by force after this many steps
    rng = np.random.default_rng() # PCG64

    def __init__(self, params):
        '''
        PARAMS must have fields [gateValve, window, cell, source, cellChance, lossPerStep, windowLoss, stepSize]
        '''
        self.particleNum = 0
        try:
            self.gateValve = getattr(params, 'gateValve')
            self.window = getattr(params, 'window')
            self.cell = getattr(params, 'cell')
            self.source = getattr(params, 'source')
            self.cellChance = getattr(params, 'cellChance')     # Chance for particle to enter cell
            self.lossPerStep = getattr(params, 'lossPerStep')   # Loss per random walk step
            self.windowLoss = getattr(params, 'windowLoss')
            self.stepSize = getattr(params, 'stepSize')

        except AttributeError as error:
            sys.stderr.write(f'\n{error}')
            sys.stderr.write('\nERROR: neutron1D.__init__(self, params)\n' + self.__init__.__doc__ + '\n')
            sys.exit()

        self.resetState()
        if self.source < self.cell: # Shouldn't matter where the source/cell are relative to each other
            self.sourceLeftCellRight = True
        else:
            self.sourceLeftCellRight = False

    def resetState(self, label = 0):
        self.particleNum = label
        self.windowHits = 0
        self.totalSteps = 0
        self.location = self.gateValve
        self.prevLocation = None
        self.status = 'alive'
        self.cellRejections = 0

    def getState(self):
        return {'particleNum':self.particleNum,
                'windowHits':self.windowHits,
                'totalSteps':self.totalSteps,
                'location':self.location,
                'status':self.status,
                'cellRejections':self.cellRejections}

    def step(self):
        self.prevLocation = self.location
        self.location += self.rng.choice([-1,1]) * self.stepSize
        self.totalSteps += 1
        # Check collisions
        if self.crossedWindow():
            self.windowHits += 1
            if self.rng.random() < self.windowLoss: # chance for loss on the window
                self.status = 'window'
        elif (self.location <= self.source and self.sourceLeftCellRight) \
              or (self.location >= self.source and not self.sourceLeftCellRight): # neutrons get absorbed by the source
            self.status = 'source'
        elif (self.location >= self.cell and self.sourceLeftCellRight) \
              or (self.location <= self.cell and not self.sourceLeftCellRight):   # chance for neutrons to get into cell (we assume they won't leave)
            if self.rng.random() < self.cellChance:
                self.status = 'cell'
            else: # neutron rejected from entrance
                self.location = self.prevLocation
                self.totalSteps += 1
                self.cellRejections += 1
        elif self.rng.random() < self.lossPerStep:  # chance to be absorbed by pipe
            self.status = 'pipe'

    def walk(self):
        while (self.totalSteps < self.maxSteps) and (self.status=='alive'):
            self.step()

    def crossedWindow(self):
        if self.prevLocation == None:
            return False
        elif self.prevLocation < self.window < self.location:
            return True
        elif self.prevLocation > self.window > self.location:
            return True
        elif self.location == self.window:
            return True
        else:
            return False
