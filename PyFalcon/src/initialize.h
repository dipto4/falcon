#ifndef INITIALIZE_H
#define INITIALIZE_H
#include<stddef.h>
void initialize(double **pos, double **vel, double **acc, double **jerk, 
        double *step, double *t_last, double *mass, size_t N, double eta); 
 
void scale(double **pos, double **vel, double *mass, size_t N,
        double* mass_scale, double *pos_scale, double *vel_scale, double *t_scale);


void move_to_com(double **pos, double **vel, double *mass, size_t N); 




#endif

