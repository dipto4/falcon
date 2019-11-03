#include<stdio.h>
#include<stddef.h>
#include<stdlib.h>
#include<math.h>
#include<assert.h>
#include<string.h>
#include "io.h"
#include "initialize.h"
#include "integrator.h"


void main(int argc, char **argv) {
    double **pos = NULL, **vel = NULL, *mass = NULL;
    double **acc = NULL, **jerk = NULL;
    size_t N;
    double time = 0, dt;
    double *step, *t_last;
    size_t nsteps = 0;
    double final_t;
    double out_frequency;
    double mass_scale, pos_scale, vel_scale, t_scale;
    size_t integrator_steps;
    
    assert(argc==4);
    
    const char* filename = argv[1];
    final_t = atof(argv[2]);
    out_frequency = atof(argc[3]);

    /* a cool code needs a cooler logo*/ 
    printf("    *                                                                                                                            \n");
    printf("     (                                                                                                                           \n");
    printf("      (                                                                                                                          \n");
    printf("       (                                                                                                                         \n");
    printf("       ((                                                                                                                        \n");
    printf("        ((                                                                                                                       \n");
    printf("         ((                                                                                                                      \n");
    printf("          ((                                                                                                                     \n");
    printf("           ((                                                                                                                    \n");
    printf("            ((*                                                                                                                  \n");
    printf("             (((                                                                                                                 \n");
    printf("              (((                                                                                                                \n");
    printf("               ((((                                                                                                              \n");
    printf("                ((((                                                                                                             \n");
    printf("                 (((((                                                                                                           \n");
    printf("                  (((((     .                                                                                                    \n");
    printf("                   /(((((     (.                                                                                                 \n");
    printf("                     (((((,     (((                                                                                              \n");
    printf("                      ((((((      ((((,                                                                                          \n");
    printf("                       (((((((      (((((((                                                                                      \n");
    printf("                        .(((((((      .((((((((,                                                                                 \n");
    printf("                          ((((((((       ((((((((((((                                                                            \n");
    printf("                           (((((((((       ((((((((((((((((*                                                                     \n");
    printf("                             (((((((((       .(((((((((((((((((((((,                                                             \n");
    printf("                              ((((((((((        ((((((((((((((((((((((((((((((*                                                  \n");
    printf("                                ((((((((((         ((((((((((((((((((((((((((((((((((((((((((((((((((((                             \n");
    printf("                                 ((((((((((((        .(((((((((((((((((((((((((((/.                  .((((*                           \n");
    printf("                                   ((((((((((((         .((((((((((((((.                       *((,     /(((                          \n");
    printf("                                     (((((((((((((          (((((                                  ((((( (((*                         \n");
    printf("                                      .((((((((((((((.                                            (((((((((((                         \n");
    printf("                                        *((((((((((((((((                                              ((((((((((((                   \n");
    printf("                                          .(((((((((((                                                    ((((((((((((                \n");
    printf("                                             ((((((                                                     ((((((((((((((((              \n");
    printf("                                               /                                                                  *((((((             \n");
    printf("                                                                                                                      ((((            \n");
    printf("                                                                                                                        (((           \n");
    printf("                                                                                                                         ((           \n");
    printf("                                                                                                                          (           \n");
    printf("\n Falcon v0.1: N-Body code for small N systems\n");
    printf("\n Reading input file...");


    
    
    /* Flow:
     
       read inputs -> move inputs to CoM -> scale the energy to N-Body units ->
       get initial jerk, acc, timestep -> print scaling constants -> start integration
        
    */
    N = get_N_from_file(filename);
    assert(N > 1);
    assert(N<NMAX);
    
    printf("\nNumber of particles: %i\n", N);
    
    //set up the dataset arrays
    
    mass = (double *) malloc(N*sizeof(double));
    step = (double *) malloc(sizeof(double)*N);
    t_last = (double *) malloc(sizeof(double)*N);
    
    pos = (double **) malloc(N*sizeof(double*));
    vel = (double **) malloc(N*sizeof(double*));
    acc = (double **) malloc(N*sizeof(double*));
    jerk = (double **) malloc(N*sizeof(double*));
    int row;
    
    for(row=0;row<N;row++) {
        pos[row] = (double *) malloc(3*sizeof(double));
        vel[row] = (double *) malloc(3*sizeof(double));
        acc[row] = (double *) malloc(3*sizeof(double));
        jerk[row] = (double *) malloc(3*sizeof(double));
    }

    // get the inputs from the file now

    read_inputs_from_file(filename,pos,vel,mass,N);
    

    // move the system to the center of mass
    
    printf("intializing the system...\n");
    move_to_com(pos,vel,mass,N);
    
    scale(pos,vel,mass,N,&mass_scale, &pos_scale, &vel_scale, &t_scale);

    printf("mass_scale = %f , pos_scale = %f, vel_scale= %f, t_scale=%f \n", mass_scale, pos_scale, vel_scale, t_scale);

    initialize(pos,vel,acc,jerk,step,t_last,mass,N);
    
    // dump the initial conditions to an output file

    write_output_to_hdf5(nsteps,mass,pos,vel,acc,jerk,N,time); 
    
    // start the simulation 
    printf("starting integration...");
    
    while(time <= final_t) {

        integrator_steps = integrate(pos,vel,acc,jerk,mass,step,t_last,&time,(t+out_frequency),N);
        nsteps++;
        write_output_to_hdf5(nsteps,mass,pos,vel,acc,jerk,N,time); 
    }
    printf("Integration completed!\n");
}
