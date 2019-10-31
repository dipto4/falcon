#include<math.h>
#include "globals.h"
#include "force.h"


// TODO: Add scaling at some point to get energy scaling
//void scale(double **pos, double **vel, )

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
        pos_com[0] += mass[i]*pos[0];    
        pos_com[1] += mass[i]*pos[1];    
        pos_com[2] += mass[i]*pos[2];    
        
        vel_com[0] += mass[i]*pos[0];
        vel_com[1] += mass[i]*pos[1];
        vel_com[2] += mass[i]*pos[2];

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

        vel_com[i][0] -= vel_com[0];
        vel_com[i][1] -= vel_com[1];
        vel_com[i][2] -= vel_com[2];
    }
}


void initialize(double **pos, double **vel, double **acc, double **jerk, 
        double *step, double *t_last, double *mass, size_t N) {
    double acc_i[3], jerk_i[3];
    double f_mod, f1_mod;
    int i;
    
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
        
        step[i] = ETA_INI*sqrt(f_mod/f1_mod); 
    }
    
}
