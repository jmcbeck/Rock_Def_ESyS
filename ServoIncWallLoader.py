#import the division module for compatibility between Python 2 and Python 3
from __future__ import division
#import the appropriate ESyS-Particle modules:
from esys.lsm import *
from esys.lsm.util import *

class ServoIncWallLoaderRunnable (Runnable):
    def __init__ (self,
                  LsmMpi=None,
                  interactionName=None,
                  maxForce=Vec3(0,0,0),
                  startTime=0,
                  maxTime = 200
              ):
        """
        Subroutine to initialise the Runnable and store parameter values.
        """
        Runnable.__init__(self)
        self.sim = LsmMpi
        self.interactionName = interactionName
        self.maxForce = maxForce
        self.dt = self.sim.getTimeStepSize()
        self.maxTime = maxTime
        self.startTime = startTime
        self.Nt = 0

    def run (self):
        """
        Subroutine to apply the force to a wall interaction. After self.startTime
        timesteps, the force on the wall increases linearly over
        self.rampTime timesteps until the desired wall force is achieved.
        Thereafter the wall force is kept fixed.
        """
        if (self.Nt > self.startTime):
            #compute the slowdown factor if still accelerating the wall:
            if (self.Nt < (self.startTime + self.maxTime)):
                f = float(self.Nt - self.startTime) / float(self.maxTime)
            else:
                f = 1.0

            #compute the amount by which to move the wall this timestep:
            Dforce = Vec3(
                f*self.maxForce[0],
                f*self.maxForce[1],
                f*self.maxForce[2]
            )

            #instruct the simulation to apply the prescribed force to the wall:
            self.sim.applyForceToWall (self.interactionName, Dforce)

        self.Nt += 1
