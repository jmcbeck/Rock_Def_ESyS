#WallLoader.py: A Runnable for moving walls in ESyS-Particle simulations
# Author: D. Weatherley
# Date: 28 December 2008
# Organisation: ESSCC, University of Queensland
# (C) All rights reserved, 2008.
#
#
#import the division module for compatibility between Python 2 and Python 3
from __future__ import division
#import the appropriate ESyS-Particle modules:
from esys.lsm import *
from esys.lsm.util import *

#This script implements a Runnable designed to move a wall at a specified
#speed. The Runnable also implements initial acceleration of the wall
#from zero to the desired speed as well as an optional initial idle
#period during which the wall does not move.
class WallLoaderRunnable (Runnable):
    def __init__ (self,
                  LsmMpi=None,
                  wallName=None,
                  vPlate=Vec3(0,0,0),
                  startTime=0,
                  rampTime = 200):
        """
        Subroutine to initialise the Runnable and store parameter values.
        """
        Runnable.__init__(self)
        self.sim = LsmMpi
        self.wallName = wallName
        self.Vplate = vPlate
        self.dt = self.sim.getTimeStepSize()
        self.rampTime = rampTime
        self.startTime = startTime
        self.Nt = 0

    def run (self):
        """
        Subroutine to move the specified wall. After self.startTime
        timesteps, the speed of the wall increases linearly over
        self.rampTime timesteps until the desired wall speed is achieved.
        Thereafter the wall is moved at that speed.
        """
        if (self.Nt >= self.startTime):
            #compute the slowdown factor if still accelerating the wall:
            if (self.Nt < (self.startTime + self.rampTime)):
                f = float(self.Nt - self.startTime) / float(self.rampTime)
            else:
                f = 1.0
            #compute the amount by which to move the wall this timestep:
            Dplate = Vec3(
                f*self.Vplate[0]*self.dt,
                f*self.Vplate[1]*self.dt,
                f*self.Vplate[2]*self.dt
            )
            #instruct the simulation to move the wall:
            self.sim.moveWallBy (self.wallName, Dplate)

        #count the number of timesteps completed thus far:
        self.Nt += 1

