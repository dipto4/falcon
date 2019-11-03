#include<stdio.h>
#include<math.h>
#include<assert.h>
#include "globals.h"
#include "force.h"
#include "integrator.h"

void update_timestep(double acc_i[3], double jerk_i[3], double a2[3], double a3[3], double *step, size_t i) {
    double mod_f, mod_f1, mod_f2, mod_f3;
    double temp_t_i;

    mod_f = sqrt(acc_i[0] * acc_i[0] + acc_i[1] * acc_i[1] + acc_i[2] * acc_i[2]);
    mod_f1 = sqrt(jerk_i[0] * jerk_i[0] + jerk_i[1] * jerk_i[1] + jerk_i[2] * jerk_i[2]);
    mod_f2 = sqrt(a2[0] * a2[0] + a2[1] * a2[1] + a2[2] * a2[2]);
    mod_f3 = sqrt(a3[0] * a3[0] + a3[1] * a3[1] + a3[2] * a3[2]);

    //find an optimal timestep
    temp_t_i = sqrt(ETA* ((mod_f * mod_f2) + mod_f1*mod_f1 ) / (mod_f1*mod_f3 + mod_f2*mod_f2));
    
    //check if timestep is not zero
    assert(temp_t_i != 0);

    step[i] = min((1.2*step[i]), temp_t_i);
    
}


void predictor(double **pos, double **vel, double **acc, double **jerk, double t,size_t N, 
        double temp_pos[NMAX][3], double temp_vel[NMAX][3], double *t_last) {
    int j;
    
    double dt;
    double dt2_over_2;
    double dt3_over_6;

    for(j = 0;j < N;j++) {
        dt = t - t_last[j];
        dt2_over_2 = dt*dt/2.0;
        dt3_over_6 = dt*dt*dt/6.0;
        
        temp_pos[j][0] = pos[j][0] + vel[j][0]*dt + acc[j][0]*dt2_over_2 + jerk[j][0]*dt3_over_6; 
        temp_pos[j][1] = pos[j][1] + vel[j][1]*dt + acc[j][1]*dt2_over_2 + jerk[j][1]*dt3_over_6; 
        temp_pos[j][2] = pos[j][2] + vel[j][2]*dt + acc[j][2]*dt2_over_2 + jerk[j][2]*dt3_over_6;

        temp_vel[j][0] = vel[j][0] + acc[j][0]*dt + jerk[j][0]*dt2_over_2;
        temp_vel[j][1] = vel[j][1] + acc[j][1]*dt + jerk[j][1]*dt2_over_2;
        temp_vel[j][2] = vel[j][2] + acc[j][2]*dt + jerk[j][2]*dt2_over_2;
    }
}

void corrector(double **pos, double **vel, double **acc, double **jerk, double *step, size_t i,
        double temp_pos[NMAX][3], double temp_vel[NMAX][3], double acc_i[3], double jerk_i[3]) {

    double a2[3], a3[3];
    double step_i2 = step[i] * step[i];
    double step_i3 = step[i] * step[i] * step[i];

    //find the second and third derivatives of acceleration from prediced and original positions
    a2[0] = (-6*(acc[i][0]-acc_i[0]) - step[i]*(4*jerk[i][0] + 2*jerk_i[0]))/(step_i2);
    a2[1] = (-6*(acc[i][1]-acc_i[1]) - step[i]*(4*jerk[i][1] + 2*jerk_i[1]))/(step_i2);
    a2[2] = (-6*(acc[i][2]-acc_i[2]) - step[i]*(4*jerk[i][2] + 2*jerk_i[2]))/(step_i2);
    
    a3[0] = (12*(acc[i][0]-acc_i[0]) + 6*step[i]*(jerk[i][0] + jerk_i[0]))/(step_i3);
    a3[1] = (12*(acc[i][1]-acc_i[1]) + 6*step[i]*(jerk[i][1] + jerk_i[1]))/(step_i3);
    a3[2] = (12*(acc[i][2]-acc_i[2]) + 6*step[i]*(jerk[i][2] + jerk_i[2]))/(step_i3);


    //correct the predicted positions for particle i
    pos[i][0] = temp_pos[i][0] + (step_i2*step_i2)*a2[0]/24 + (step_i2*step_i3)*a3[0]/120;
    pos[i][1] = temp_pos[i][1] + (step_i2*step_i2)*a2[1]/24 + (step_i2*step_i3)*a3[1]/120;
    pos[i][2] = temp_pos[i][2] + (step_i2*step_i2)*a2[2]/24 + (step_i2*step_i3)*a3[2]/120;
    
    vel[i][0] = temp_vel[i][0] + (step_i3)*a2[0]/6 + (step_i2*step_i2)*a3[0]/24;
    vel[i][1] = temp_vel[i][1] + (step_i3)*a2[1]/6 + (step_i2*step_i2)*a3[1]/24;
    vel[i][2] = temp_vel[i][2] + (step_i3)*a2[2]/6 + (step_i2*step_i2)*a3[2]/24;
    
    //update the force and jerk for particle i

    acc[i][0] = acc_i[0];
    acc[i][1] = acc_i[1];
    acc[i][2] = acc_i[2];
    
    jerk[i][0] = jerk_i[0];
    jerk[i][1] = jerk_i[1];
    jerk[i][2] = jerk_i[2];

    update_timestep(acc_i,jerk_i,a2,a3,step,i);
}



size_t integrate(double **pos, double **vel, double **acc, double **jerk,double *mass, 
        double *step, double* t_last, double *t, double t_end, 
        size_t N) {
    
    double temp_pos[NMAX][3], temp_vel[NMAX][3], acc_i[3], jerk_i[3];
    
    int i, j,i_min;
    size_t nsteps;
    double t_min;

    assert(N < NMAX);

    nsteps = 0;
    do {
        t_min = HUGE;
        
        //find the next smallest timestep 

        for(i=0; i < N; i++) {
            if((step[i] + t_last[i]) < t_min) { 
                t_min = step[i] + t_last[i];    
                i_min = i;
            }

        } // end for

        assert(i_min < N);

        *t = t_min;
        i = i_min;
        
        // loop over all particles to predict their new positions
        
        // O(N) steps
        predictor(pos, vel, acc, jerk, *t, N, temp_pos, temp_vel, t_last);

        // O(N) steps
        calculate_acceleration_and_jerk_single(pos, vel, mass, N, i, acc_i, jerk_i);

        // O(1) steps
        corrector(pos, vel, acc, jerk, step, i,temp_pos, temp_vel, acc_i, jerk_i);

        t_last[i] = *t;

        nsteps++;
#ifdef DEBUG
        printf("Integration step completed. Current time: %f \n",*t);
#endif
    } while(*t < t_end);

    return nsteps;
}
