#include <math.h>
#include "force.h"


// calculate_acceleration_and_jerk takes following inputs:
// double **pos: positions of N particles
// double **vel: velocities of N particles
// double **mass: masses of N particles
// size_t N: number of particles in the simulation
//
// calculate_acceleration_and_jerk produces the following outputs
// double **acc: acceleration of N particles
// double **jerk: time derivative of acceleration of N particles
// calculate_acceleration_and_jerk calculates acceleration using the following formulae:
// F_ij_x = -m_j*R_x/|R|^3
// J_ij_x = -m_j*V_x/|R|^3 -3*(R_x*V_x+R_y*V_y+R_z*V_z)*F_ij_x
//
// acc[i][0] = acc_i_x
// acc[i][1] = acc_i_y
// acc[i][2] = acc_i_z

void calculate_potential(double **pos, double *mass, size_t N, double pot[NMAX]) {
    size_t i,j;

    for(i = 0; i < N;i++) {
        pot[i] = 0;
    }


    for(i=0;i<N;i++) {
        for(j=0;j<N;j++) {
            if(i != j) {
                double dx = pos[i][0] - pos[j][0];
                double dy = pos[i][1] - pos[j][1];
                double dz = pos[i][2] - pos[j][2];

                double dist = sqrt(dx*dx + dy*dy + dz*dz);
                
                pot[i] += mass[j] / dist;

            } else {
                continue;
            }
        }
    }
}


void calculate_acceleration_and_jerk(double **pos,double **vel, double *mass, size_t N, double **acc, double **jerk,
        size_t req_i, double *acc_i, double *jerk_i) {

    size_t i,j;

    for(i = 0; i <N; i++) {
        for( j = 0; j < N; j++) {
            if(i != j) {
                // calculate |R|
                double dx = pos[i][0] - pos[j][0];
                double dy = pos[i][1] - pos[j][1];
                double dz = pos[i][2] - pos[j][2];

                double dist = sqrt(dx*dx + dy*dy + dz*dz);
                double dist2 = dist*dist;
                double dist3 = dist*dist*dist;

                double mass_over_dist3 = -mass[j]/dist3;

                // calculate F_ij
                acc[i][0] += -mass_over_dist3*dx;
                acc[i][1] += -mass_over_dist3*dy;
                acc[i][2] += -mass_over_dist3*dz;

                // calculate J_ij
                double dvx = vel[i][0] - vel[j][0];
                double dvy = vel[i][1] - vel[j][1];
                double dvz = vel[i][2] - vel[j][2];
                // a = R.V/|R|^2
                double a = dx*dvx+dy*dvy+dz*dvz/dist2;


                jerk[i][0] += -mass_over_dist3*dvx-3*a*acc[i][0];
                jerk[i][1] += -mass_over_dist3*dvy-3*a*acc[i][1];
                jerk[i][2] += -mass_over_dist3*dvz-3*a*acc[i][2];

            } else {
                continue;
            }
        }

    }
    acc_i = acc[req_i];
    jerk_i = jerk[req_i];

}

void calculate_acceleration_and_jerk_single(double **pos,double **vel, double *mass, size_t N, 
        size_t req_i, double acc_i[3], double jerk_i[3]) {

    size_t i,j;
    i = req_i;

    for( j = 0; j < N; j++) {
        if(i != j) {
            // calculate |R|
            double dx = pos[i][0] - pos[j][0];
            double dy = pos[i][1] - pos[j][1];
            double dz = pos[i][2] - pos[j][2];

            double dist = sqrt(dx*dx + dy*dy + dz*dz);
            double dist2 = dist*dist;
            double dist3 = dist*dist*dist;

            double mass_over_dist3 = -mass[j]/dist3;

            // calculate F_ij
            acc_i[0] += -mass_over_dist3*dx;
            acc_i[1] += -mass_over_dist3*dy;
            acc_i[2] += -mass_over_dist3*dz;

            // calculate J_ij
            double dvx = vel[i][0] - vel[j][0];
            double dvy = vel[i][1] - vel[j][1];
            double dvz = vel[i][2] - vel[j][2];
            // a = R.V/|R|^2
            double a = dx*dvx+dy*dvy+dz*dvz/dist2;


            jerk_i[0] += -mass_over_dist3*dvx-3*a*acc_i[0];
            jerk_i[1] += -mass_over_dist3*dvy-3*a*acc_i[1];
            jerk_i[2] += -mass_over_dist3*dvz-3*a*acc_i[2];

        } else {
            continue;
        }
    }


}

