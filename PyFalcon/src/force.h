#ifndef FORCE_H
#define FORCE_H
#include<stddef.h>
#include "globals.h"
void calculate_potential(double **pos, double *mass, size_t N, double pot[NMAX]);
void calculate_acceleration_and_jerk(double **pos,double **vel, double *mass, size_t N, double **acc, double **jerk);

void calculate_acceleration_and_jerk_single(double **pos, double **vel, double *mass, size_t N,
        size_t req_i, double acc_i[3], double jerk_i[3]);
#endif
