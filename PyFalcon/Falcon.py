from __future__ import print_function
import numpy as np
import os
from optparse import OptionParser
from ctypes import cdll, c_int, c_double, byref, c_bool
from particle import Particle

class Falcon(object):
    __instance__ == None

    def __init__(self):
        Falcon.__instance__ = self

        self.library = cdll.LoadLibrary('./libfalcon.so')

        # input parameters
        # N = number of particles
        # ETA = accuracy parameter
        # OUT_N = number of outputs
        # FINAL_T = duration of integration

        # particle parameters
        # X , Y, Z = particle positions
        # VX, VY, VZ = particle velocities
        # M = particle mass

        # TODO:
        # include a parameter for a debug mode


        # default values

        self.N = 0
        self.eta = 0.02
        self.out_n = 100
        self.final_t = 1

        self.Particles = None


    def create_particles(self, N):
        __error__ = 0

        __particles__ = np.empty(N,dtype=Particle)

        for i in xrange(0,N):
            __particles__[i] = Particle()

        self.Particles = __particles__
        self.N = N
        if(__particles__.size == N):
            __error__ = 0
        else:
            __error__ = 1

        return __error__


    def __set_falcon_values__(self):

        if(self.N == None or self.N < 2):
            print("Error in the number of particles! Check...")

            # add error handler

        else:
            for i in range(0,self.N):
                if(self.Particles[i].isNone == True):
                    print("Error in initialization of particle {}".format(i))

                    # add error handler




