#include<stdio.h>
#include<stddef.h>
#include<stdlib.h>
#include<math.h>
#include<assert.h>
#include<string.h>
#include "common.h"
#include "io.h"
#include "initialize.h"
#include "integrator.h"


void start(input_params* in, particle_params* pa) {
    // TODO: add extra stuff to transfer input and particle params
    
    double **pos = NULL, **vel = NULL, *mass = NULL;
    double **acc = NULL, **jerk = NULL;
    size_t N;
    double eta;
    double time = 0, dt;
    double *step, *t_last;
    size_t nsteps = 0;
    double final_t;
    double out_frequency;
    double mass_scale, pos_scale, vel_scale, t_scale;
    size_t integrator_steps = 0;
    size_t total_steps = 0;
    double e0, kin, pot;
    double energy_current;


    // input parameters from file
    N = in->N;
    eta = in->eta;
    out_frequency = in->dt;
    final_t = in->final_t;
    
    printf("N = %d\n",N); 
    printf("eta = %f\n",eta); 
    printf("out_fre = %d\n",out_frequency); 
    printf("final_t = %f\n",final_t); 
    
    
    
    //assert(argc==4);
    
    //const char* filename = argv[1];
    //final_t = atof(argv[2]);
    //out_frequency = atof(argv[3]);

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
    //N = get_N_from_file(filename);
    //printf("N = %d\n",N);
    
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



    // TODO: change the following line to get the data from python
    //read_inputs_from_file(filename,pos,vel,mass,N);
    
    for(row = 0 ; row < N ; row++) {
        pos[row][0] = pa->x[row];
        pos[row][1] = pa->y[row];
        pos[row][2] = pa->z[row];

        vel[row][0] = pa->vx[row];
        vel[row][1] = pa->vy[row];
        vel[row][2] = pa->vz[row];

        mass[row] = pa->m[row];
    }

    // move the system to the center of mass
    
    printf("intializing the system...\n");
    move_to_com(pos,vel,mass,N);
    
    //skip scaling for now
    //scale(pos,vel,mass,N,&mass_scale, &pos_scale, &vel_scale, &t_scale);

    printf("mass_scale = %f , pos_scale = %f, vel_scale= %f, t_scale=%f \n", mass_scale, pos_scale, vel_scale, t_scale);

    initialize(pos,vel,acc,jerk,step,t_last,mass,N,eta);
    
    // dump the initial conditions to an output file

    write_output_to_hdf5(nsteps,mass,pos,vel,acc,jerk,N,time); 
    
    kin = get_kinetic_energy(vel,mass,N);
    pot = get_potential_energy(pos, mass, N);
    e0 = kin - pot;


    printf("kinetic energy = %f, potential_energy = %f, total energy = %f\n\n",kin,pot,e0);   
    // start the simulation 
    printf("starting integration...");
    
    while(time <= final_t) {

        integrator_steps = integrate(pos,vel,acc,jerk,mass,step,t_last,&time,(time+out_frequency),N,eta);
        total_steps += integrator_steps;
        
        nsteps++;
            
        kin = get_kinetic_energy(vel,mass,N);
        pot = get_potential_energy(pos, mass, N);
        energy_current = kin - pot;
        
        printf("steps taken by integrator: %d\n",integrator_steps);
        printf("time = %f\n",time);
        printf("kinetic energy = %18.12f, potential_energy = %18.12f, total energy = %18.12f\n",kin,pot,energy_current);
        printf("relative energy error = %.16e\n\n", fabs(e0-energy_current)/e0);

        write_output_to_hdf5(nsteps,mass,pos,vel,acc,jerk,N,time); 
    }
    printf("Total steps taken %lu\n",total_steps);
    printf("Integration completed!\n");
    free(pos); free(vel); free(mass); free(acc); free(jerk);
}
