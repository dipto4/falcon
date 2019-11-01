#include<stddef.h>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<assert.h>
#include<string.h>
#include "hdf5.h"
#include "globals.h"
#include "io.h"

double get_potential_energy(double **pos, double *mass, size_t N) {
    size_t i;
    double potential_energy = 0;
    double pot[NMAX];

    calculate_potential(pos,mass,N,pot);

    for (i=0;i<N;i++) {
        potential_energy += mass[i] *pot[i];
    }
    
    return potential_energy;

}

double get_kinetic_energy(double **vel, double *mass, size_t N) {
    size_t i;
    
    double kinetic_energy = 0;
    for(i=0;i<N;i++) {

        vel_sq = vel[i][0] * vel[i][0] + vel[i][1] * vel[i][1] + vel[i][2] * vel[i][2];

        kinetic_energy += 0.5*mass[i]*vel_sq;
    }
    
    return kinetic_energy;
}

size_t get_N_from_file(const char* filename) {
    FILE *input;

    char *line = NULL;
    size_t len = 0;
    ssize_t read;

    input = fopen(filename, "r");
    size_t line_count = 0;

    while ((read = getline(&line, &len, input)) != -1) {
        line_count++;
    
 
    }

    fclose(input);
    assert(line_count > 1);
    return line_count;
}


void read_inputs_from_file(const char* filename, double **pos, double **vel, double *mass, size_t N) {
    FILE *input;

    char *line = NULL;
    size_t len = 0;
    ssize_t read;

    input = fopen(filename, "r");
    size_t line_count = 0;
    
    double params[7];

    while ((read = getline(&line, &len, input)) != -1) {
        char *l = line;
        char *token;
        
        size_t param_count = 0;
        
        while ((token = strsep(&l, " ")) != NULL) {
            params[param_count] = atof(token);
            param_count++;
        }

        mass[line_count] = params[0];
        pos[line_count][0] = params[1];
        pos[line_count][1] = params[2];
        pos[line_count][2] = params[3];
        vel[line_count][0] = params[4];
        vel[line_count][1] = params[5];
        vel[line_count][2] = params[6];

        line_count++;

    }

    fclose(input);

}

// Used from Taichi io.h with appropriate modifications

void write_output_to_hdf5(size_t output_num, double *mass, double **pos, double **vel, double **acc, double **jerk, 
        size_t N, double t) {
    hsize_t dims[1] = { N };
    herr_t status;   hid_t file_id, space_id, dset_id, memspace, handle = 0;
    hid_t file_id, space_id, dset_id, memspace, handle = 0;
    char buf[500];
    sprintf(buf, "output_%d.hdf5", snapshot_num);

    double *data = (double *) malloc(N*sizeof(double));
    double pot[NMAX];
    
    file_id = H5Fcreate(buf, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    // first write the header information
    handle = H5Gcreate(file_id, "/info", 0);

    hid_t hdf5_dataspace, hdf5_attribute;

    printf("Writing %d particles at t=%g \n", N, t);
    fflush(stdout);

    hdf5_dataspace = H5Screate(H5S_SCALAR);
    hdf5_attribute = H5Acreate(handle, "Time", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
    H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &t);
    H5Aclose(hdf5_attribute);
    H5Sclose(hdf5_dataspace);
    //Header writing done 

 
    memspace = H5Screate_simple(1, dims, NULL);
    space_id = H5Screate_simple(1, dims, NULL);
    for(unsigned int b = 0; b < s.n; b++)
      data[b] = pos[b][0];
    dset_id = H5Dcreate(file_id, "Posx", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT);
    status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, space_id, H5P_DEFAULT, data);
    status = H5Sclose(space_id);
    status = H5Sclose(memspace);
    status = H5Dclose(dset_id);

    memspace = H5Screate_simple(1, dims, NULL);
    space_id = H5Screate_simple(1, dims, NULL);
    for(unsigned int b = 0; b < s.n; b++)
      data[b] = pos[b][1];
    dset_id = H5Dcreate(file_id, "Posy", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT);
    status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, space_id, H5P_DEFAULT, data);
    status = H5Sclose(space_id);
    status = H5Sclose(memspace);
    status = H5Dclose(dset_id);

    memspace = H5Screate_simple(1, dims, NULL);
    space_id = H5Screate_simple(1, dims, NULL);
    for(unsigned int b = 0; b < s.n; b++)
      data[b] = pos[b][2];
    dset_id = H5Dcreate(file_id, "Posz", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT);
    status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, space_id, H5P_DEFAULT, data);
    status = H5Sclose(space_id);
    status = H5Sclose(memspace);
    status = H5Dclose(dset_id);

    memspace = H5Screate_simple(1, dims, NULL);
    space_id = H5Screate_simple(1, dims, NULL);
    for(unsigned int b = 0; b < s.n; b++)
      data[b] = vel[b][0];
    dset_id = H5Dcreate(file_id, "Velx", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT);
    status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, space_id, H5P_DEFAULT, data);
    status = H5Sclose(space_id);
    status = H5Sclose(memspace);
    status = H5Dclose(dset_id);

    memspace = H5Screate_simple(1, dims, NULL);
    space_id = H5Screate_simple(1, dims, NULL);
    for(unsigned int b = 0; b < s.n; b++)
      data[b] = vel[b][1];
    dset_id = H5Dcreate(file_id, "Vely", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT);
    status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, space_id, H5P_DEFAULT, data);
    status = H5Sclose(space_id);
    status = H5Sclose(memspace);
    status = H5Dclose(dset_id);

    memspace = H5Screate_simple(1, dims, NULL);
    space_id = H5Screate_simple(1, dims, NULL);
    for(unsigned int b = 0; b < s.n; b++)
      data[b] = vel[b][2];
    dset_id = H5Dcreate(file_id, "Velz", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT);
    status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, space_id, H5P_DEFAULT, data);
    status = H5Sclose(space_id);
    status = H5Sclose(memspace);
    status = H5Dclose(dset_id);

    memspace = H5Screate_simple(1, dims, NULL);
    space_id = H5Screate_simple(1, dims, NULL);
    for(unsigned int b = 0; b < s.n; b++)
      data[b] = mass[b];
    dset_id = H5Dcreate(file_id, "Mass", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT);
    status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, space_id, H5P_DEFAULT, data);
    status = H5Sclose(space_id);
    status = H5Sclose(memspace);
    status = H5Dclose(dset_id);


    // output the potentials
    calculate_potential(pos,mass,pot);
    
    for(unsigned int b = 0; b < s.n; b++)
        data[b] = pot[b];
    dset_id = H5Dcreate(file_id, "Potential", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT);
    status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, space_id, H5P_DEFAULT, data);
    status = H5Sclose(space_id);
    status = H5Sclose(memspace);
    status = H5Dclose(dset_id);
    
    free(data);
    
    // close the file
    status = H5Gclose(handle);
    status = H5Fclose(file_id);
    if(status < 0)
        printf("Writing snapshot error\n");

}
