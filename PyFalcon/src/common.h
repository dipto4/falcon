#ifndef COMMON_H
#define COMMON_H

typedef struct input_params {
    int N;
    double eta;
    int out_n;
    double final_t;
} input_params;

typedef struct particle_params {
    double *x;
    double *y;
    double *z;
    double *vx;
    double *vy;
    double *vz;
    double *m;
} particle_params;

#endif
