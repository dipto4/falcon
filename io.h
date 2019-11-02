#ifndef IO_H
#define IO_H
#include<stddef.h>
//void read_from_input_file()

double get_potential_energy(double **pos, double *mass, size_t N);

double get_kinetic_energy(double **vel, double *mass, size_t N);

size_t get_N_from_file(const char* filename);

void read_inputs_from_file(const char* filename, double **pos, double **vel, double *mass, size_t N);

void write_output_to_hdf5(size_t output_num, double *mass, double **pos, double **vel, double **acc, double **jerk, 
        size_t N, double t); 
 
#endif
