#MeshLoader.py: A Runnable for moving meshes in esysparticle

#import the division module for compatibility between Python 2 and Python 3
from __future__ import division
#import the appropriate ESyS-Particle modules:
from esys.lsm import *
from esys.lsm.util import *

class MeshLoaderRunnable (Runnable):
    def __init__ (self,
                  LsmMpi=None,
                  meshName=None,
                  vPlate=Vec3(0,0,0),
                  startTime=0,
                  rampTime = 200,
                  endTime = 100000000000
              ):
        """
        Subroutine to initialise the Runnable and store parameter values.
        """
        Runnable.__init__(self)
        self.sim = LsmMpi
        self.meshName = meshName
        self.vPlate = vPlate
        self.dt = self.sim.getTimeStepSize()
        self.rampTime = rampTime
        self.startTime = startTime
        self.endTime = endTime
        self.Nt = 0

    def run (self):
        """
        Subroutine to apply the force to a wall interaction. After self.startTime
        timesteps, the force on the wall increases linearly over
        self.rampTime timesteps until the desired wall force is achieved.
        Thereafter the wall force is kept fixed.
        """
        if ((self.Nt > self.startTime) and (self.Nt < self.endTime)):
            #compute the slowdown factor if still accelerating the wall:
            if (self.Nt < (self.startTime + self.rampTime)):
                f = float(self.Nt - self.startTime) / float(self.rampTime)
            else:
                f = 1.0

            # calculate distance to move the meshwall
            Dplate = Vec3(
                f*self.vPlate[0]*self.dt,
                f*self.vPlate[1]*self.dt,
                f*self.vPlate[2]*self.dt
            )

            #instruct the simulation to move the mesh
            self.sim.translateMeshBy (self.meshName, Dplate)

        #count the number of timesteps completed thus far:
        self.Nt += 1

