from __future__ import print_function
import numpy as np
import os
from optparse import OptionParser
from ctypes import *
from particle import Particle

class Falcon(object):
    __instance__ = None

    def __init__(self):
        Falcon.__instance__ = self

        self.library = cdll.LoadLibrary('./libfalcon.so')
        #self.library = cdll.LoadLibrary('./libtest.so')

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
        self.dt = 0.5
        self.final_t = 1

        self.Particles = None

        self.__ctypes_input_params__ = None
        self.__ctypes_particle_params__ = None

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
                    break
                    # add error handler


        # convert to ctypes


        P_x = (c_double * self.N)(*[self.Particles[i].x for i in range(0,self.N)])
        P_y = (c_double * self.N)(*[self.Particles[i].y for i in range(0,self.N)])
        P_z = (c_double * self.N)(*[self.Particles[i].z for i in range(0,self.N)])
        P_vx = (c_double * self.N)(*[self.Particles[i].vx for i in range(0,self.N)])
        P_vy = (c_double * self.N)(*[self.Particles[i].vy for i in range(0,self.N)])
        P_vz = (c_double * self.N)(*[self.Particles[i].vz for i in range(0,self.N)])
        P_m = (c_double * self.N)(*[self.Particles[i].m for i in range(0,self.N)])

        cast(P_x,POINTER(c_double))
        cast(P_y,POINTER(c_double))
        cast(P_z,POINTER(c_double))
        cast(P_vx,POINTER(c_double))
        cast(P_vy,POINTER(c_double))
        cast(P_vz,POINTER(c_double))
        cast(P_m,POINTER(c_double))


        # now pass them to the code via the C interface

        # C interface returns pointers of type input_params and particle_params

        class __input_params__(Structure):
            _fields_=[("N",c_int),
                    ("eta",c_double),
                    ("out_n",c_double),
                    ("final_t",c_double)]


        class __particle_params__(Structure):
            _fields_=[("x",POINTER(c_double)),
                    ("y",POINTER(c_double)),
                    ("z",POINTER(c_double)),
                    ("vx",POINTER(c_double)),
                    ("vy",POINTER(c_double)),
                    ("vz",POINTER(c_double)),
                    ("m",POINTER(c_double))]

        self.__ctypes_input_params__ = __input_params__(self.N,self.eta,self.dt,self.final_t)
        self.__ctypes_particle_params__ = __particle_params__(P_x,P_y,P_z,P_vx,P_vy,P_vz,P_m)


    def integrate(self):
        self.__set_falcon_values__()
        print("starting integration...")
        self.library.start(byref(self.__ctypes_input_params__),byref(self.__ctypes_particle_params__))

