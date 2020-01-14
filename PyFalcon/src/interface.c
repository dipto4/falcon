#include "common.h"
// TODO:
// add functions to set values from the python interface to the C interface
// add functions to get values from the C interface to the python interface
//


input_params* __set_input_params__(int N, double eta, int out_n, double final_t) {
    input_params i = {N,eta,out_n,final_n};
    return &i;
}


particle_params* __set_particle_params__(double *x, double* y, double* z,
        double* vx, double* vy, double* vz,
        double *m) {

    particle_params p = {x,y,z,vx,vy,vz,m};
    return &p;

}


