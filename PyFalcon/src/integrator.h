#ifndef INTEGRATOR_H
#define INTEGRATOR_H
#include<stddef.h>
#include "globals.h"


size_t integrate(double **pos, double **vel, double **acc, double **jerk, double *mass, 
        double *step, double* t_last, double *t, double t_end, 
        size_t N, double eta); 

#endif
