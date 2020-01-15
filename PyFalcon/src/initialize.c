#include<math.h>
#include<stdio.h>
#include "globals.h"
#include "force.h"
#include "io.h"

// TODO: Add scaling at some point to get energy scaling
// Input is in arbitrary code units that get transformed to N-Body units

void scale(double **pos, double **vel, double *mass, size_t N,
        double* mass_scale, double *pos_scale, double *vel_scale, double *t_scale) {
    // mass scaling
    // scaling: total mass = 1
    
    double initial_total_mass = 0;
    double kinetic_total, potential_total;
    double Q_v;
    double beta;


    size_t i = 0;
    for(i = 0; i<N;i++) {
        initial_total_mass+=mass[i];
    }

    for(i = 0; i<N;i++) {
        mass[i] = mass[i]/initial_total_mass;
    }
    
    *mass_scale = initial_total_mass;

    // scaling of masses complete
    

    //scaling of energy, position and velocites

    kinetic_total = get_kinetic_energy(vel, mass, N);   
    potential_total = get_potential_energy(pos, mass, N);

    potential_total = fabs(potential_total);
    
    Q_v = sqrt(QVIR*potential_total/kinetic_total);

    for(i = 0; i < N; i++) {
        vel[i][0] = vel[i][0] * Q_v;
        vel[i][1] = vel[i][1] * Q_v;
        vel[i][2] = vel[i][2] * Q_v;
    }
    //scaling total energy to -0.25
    beta = -0.25/((QVIR-1.0)*potential_total);
    
    for(i = 0; i<N;i++) {
        pos[i][0] = pos[i][0]/beta;
        pos[i][1] = pos[i][1]/beta;
        pos[i][2] = pos[i][2]/beta;

        vel[i][0] = vel[i][0]*sqrt(beta);
        vel[i][1] = vel[i][1]*sqrt(beta);
        vel[i][2] = vel[i][2]*sqrt(beta);
    }

    *pos_scale = beta;
    *vel_scale = 1./(sqrt(beta)*Q_v);
    *t_scale = *pos_scale/ (*vel_scale);

}

void move_to_com(double **pos, double **vel, double *mass, size_t N) {
    //calculate center of mass for the system
    double pos_com[3], vel_com[3];
    double tot_mass = 0;
    
    pos_com[0] = 0;
    pos_com[1] = 0;
    pos_com[2] = 0;
    
    vel_com[0] = 0;
    vel_com[1] = 0;
    vel_com[2] = 0;

    int i;
    
    for(i = 0; i < N; i++) {
        pos_com[0] += mass[i]*pos[i][0];    
        pos_com[1] += mass[i]*pos[i][1];    
        pos_com[2] += mass[i]*pos[i][2];    
        
        vel_com[0] += mass[i]*vel[i][0];
        vel_com[1] += mass[i]*vel[i][1];
        vel_com[2] += mass[i]*vel[i][2];

        tot_mass += mass[i];
    }

    pos_com[0] = pos_com[0]/tot_mass;
    pos_com[1] = pos_com[1]/tot_mass;
    pos_com[2] = pos_com[2]/tot_mass;
    
    vel_com[0] = vel_com[0]/tot_mass;
    vel_com[1] = vel_com[1]/tot_mass;
    vel_com[2] = vel_com[2]/tot_mass;
    
    for(i = 0; i<N; i++) {
        pos[i][0]-=pos_com[0];
        pos[i][1]-=pos_com[1];
        pos[i][2]-=pos_com[2];

        vel[i][0] -= vel_com[0];
        vel[i][1] -= vel_com[1];
        vel[i][2] -= vel_com[2];
    }
}


void initialize(double **pos, double **vel, double **acc, double **jerk, 
        double *step, double *t_last, double *mass, size_t N, double eta) {
    double acc_i[3], jerk_i[3];
    double f_mod, f1_mod;
    int i;
    double eta_ini = 0.1*eta;    
    // initialize t_last and step
    for(i = 0; i<N;i++) {
        t_last[i] = 0.0;
        
        calculate_acceleration_and_jerk_single(pos,vel,mass, N, i, acc_i, jerk_i);
        

        // set the acceleration and jerks
        acc[i][0] = acc_i[0];
        acc[i][1] = acc_i[1];
        acc[i][2] = acc_i[2];
        
        jerk[i][0] = jerk_i[0];
        jerk[i][1] = jerk_i[1];
        jerk[i][2] = jerk_i[2];
        
        // find the initial timestep
        //

        f_mod = acc_i[0] * acc_i[0] + acc_i[1] * acc_i[1] + acc_i[2] * acc_i[2];
        f1_mod = jerk_i[0] * jerk_i[0] + jerk_i[1] * jerk_i[1] + jerk_i[2] * jerk_i[2];
        

#ifdef DEBUG
        double step_i = eta_ini*sqrt(f_mod/f1_mod);
        printf("step_i %0.16e for i = %d", step_i, i);
#endif
        step[i] = eta_ini*sqrt(f_mod/f1_mod); 
  
        if(f1_mod < TINY) {
            step[i] = 1e-6;
        }
   
    
    }
    
}
