import numpy as np
import sys

class neutron1D:
    maxSteps = 1000 # Stops by force after this many steps
    rng = np.random.default_rng() # PCG64

    def __init__(self, params):
        '''
        PARAMS must have fields [gateValve, window, cell, source, cellChance, lossPerStep, windowLoss]
        Random walk step size = 1
        All fields other than 'window' should be an int
        '''
        self.particleNum = 0
        try:
            self.gateValve = params.gateValve
            self.window = params.window
            self.cell = params.cell
            self.source = params.source
            self.cellChance = params.cellChance     # Chance for particle to enter cell
            self.lossPerStep = params.lossPerStep   # Loss per random walk step
            self.windowLoss = params.windowLoss
            if self.window.is_integer():
                raise TypeError('For window hit logging, please make params.window a non-integer')
        except AttributeError as error:
            sys.stderr.write('\nERROR: neutron1D.__init__(self, params)\n' + self.__init__.__doc__ + '\n')
            sys.exit()
        except TypeError as error:
            sys.stderr.write('\nERROR: neutron1D.__init__(self, params)\n' + self.__init__.__doc__ + '\n')
            sys.stderr.write(str(error) + '\n\n')
            sys.exit()

        self.resetState()
        if self.source < self.cell: # bounds of the 1D random walk line
            self.leftBound = self.source
            self.rightBound = self.cell
        else:
            self.leftBound = self.cell
            self.rightBound = self.source

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
        self.location += self.rng.choice([-1,1])
        self.totalSteps += 1
        # Check collisions
        if self.crossedWindow():
            self.windowHits += 1
            if self.rng.random() < self.windowLoss: # chance for loss on the window
                self.status = 'window'
        elif self.location == self.source: # neutrons get absorbed by the source
            self.status = 'source'
        elif self.location == self.cell:   # chance for neutrons to get into cell (we assume they won't leave)
            if self.rng.random() < self.cellChance:
                self.status = 'cell'
            else: # neutron rejected from entrance
                self.prevLocation = self.location
                self.location -= 1
                self.totalSteps += 1
                self.cellRejections += 1
        elif self.rng.random() < self.lossPerStep:  # chance to be absorbed by pipe
            self.status = 'pipe'
        elif (self.location > self.rightBound) or (self.location < self.leftBound):
            self.status = 'error'
            raise Exception('Something went wrong. Neutron out of bounds')


    def walk(self):
        while (self.totalSteps < self.maxSteps) and (self.status=='alive'):
            self.step()

    def crossedWindow(self):
        if self.prevLocation < self.window < self.location:
            return True
        elif self.prevLocation > self.window > self.location:
            return True
        else:
            return False
