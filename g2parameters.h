#ifndef PARAMETERS_INCLUDED
#define PARAMETERS_INCLUDED
/****************************
 PARAMETERS FILE
 ****************************/

/*
 This headr file contains the adjustable and non adjustable global parameters for the GABE2.0 program.
 
 Copyright (2013): 
 Kenyon College
 John T. Giblin, Jr
 Tate Deskins and Hillary Child
 Last Updated: 06.27.2013
*/

/***************
 model constants
 ***************/
#define num_flds 2    // number of fields
extern real_t mphi;   // mass of phi field
extern real_t phi0;   // initial avg phi field value
extern real_t gsq;    // g^2 value for phi chi coupling
extern real_t f0[2];  // array storing intial phi and chi field values
extern real_t df0[2]; // array storing intial phi and chi field derivative values
extern real_t grav;   // the gravitational constant
extern real_t c;      // speed of light, shouldn't change, if you do, fix all equations

/***************************
 model independent parameters
 ***************************/
const int randseed = 44463132; // seed for rand number generator
const int N = 32;              // number of points along one side of grid
const int NX = N;              // number of points along one side of grid
const int NY = N;              // number of points along one side of grid
const int NZ = N;              // number of points along one side of grid
const int POINTS = NX*NY*NZ;   // total number of points
extern real_t L;               // length of one side of box in prgm units
extern real_t starttime;       // start time of simulation
extern real_t endtime;         // end time of simulations
extern real_t dt;              // time step size
// (you will need to adjust functions file (adot and such) and type two
// evolution in the step() function fnd g2init.cpp initexpansion() for
// user defined expansion )

/*************************
 model dependent parameters
 *************************/

const real_t rescale_B = mphi; // rescalings

#define rand_init 1               // 1 to have random initialization 0 to not (see model file)
#define field_full_rand 1         // 1 to have full random 0 to have symmetric kspace initilaization

/*********************************
 These are important DO NOT CHANGE
 *********************************/
const int nflds = num_flds;                 // stores number of fields for looping
const real_t dx = L/((real_t) N); // stores the change in x from point to point
const real_t gridsize = NX*NY*NZ;      // stores size of grid for averaging

#endif
