#ifndef INTERFACE_H
#define INTERFACE_H

#include "common.h"

input_params* __set_input_params__(int N, double eta, int out_n, double final_t);

particle_params* __set_particle_params__(double *x, double* y, double* z,
        double* vx, double* vy, double* vz,
        double *m);
#endif
