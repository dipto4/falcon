#include<stddef.h>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<assert.h>
#include<string.h>
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


void read_stars_from_file(const char* filename, double **pos, double **vel, double *mass, size_t N) {
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
