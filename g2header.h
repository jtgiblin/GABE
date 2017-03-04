/**********************************
 Header File (declarations)
 **********************************/

/*
 This header file contains all the declarations for program functions.
 
 
 Copyright (2013): 
 Kenyon College
 John T. Giblin, Jr
 Tate Deskins and Hillary Child
*/

#include <stdlib.h> // needed for malloc
#include <math.h>   // needed for fabs, tanh
#include <time.h>   // needed for time and ctime
#include <string>
#include <unistd.h>

typedef float real_t; // type used for simulation

//this is a directive to declare all the common indecies for the fields
#define DECLARE_INDEX int fld,i,j,k; 

/****************
LOOP DEFINITIONS
 ****************/
//the following are definitions for looping over different field indicies

// this loops over the fld and all three indicies
#define fldLOOP for(fld=0;fld<nflds;fld++) for(i=0;i<NX;i++) for(j=0;j<NY;j++) for(k=0;k<NZ;k++)

// this loops over just the indicies (3D)
#define LOOP for(i=0;i<NX;i++) for(j=0;j<NY;j++) for(k=0;k<NZ;k++)

// this loops over just the j and k indicies (2d)
#define LOOP2 for(j=0;j<NY;j++) for(k=0;k<NZ;k++)

// Indexing macro for arrays
#define IDX(i,j,k) ((i)*NY*NZ + (j)*NZ + k)
#define FIELD(s,fld) ((fld)*2 + (s))

#include "g2parameters.h"
//this is the include statement for the parameters file

/******************************
 global parameters DO NOT CHANGE
 ******************************/
extern real_t t;                          // this is the variable that stores the evolution time

extern real_t ** field;  // this stores the field values for each step along the grid
extern real_t ** dfield; // this stores the derivative of the field for each step along the grid

//The following are all arrays of length two, one for each step of the RK2 integration
extern real_t a[2];      // this stores the scale facator for each step
extern real_t adot[2];   // this stores the time derivative of the scale factor
extern real_t edpot[2];  // this stores the average potential energy
extern real_t edkin[2];  // this stores the average kinetic energy
extern real_t edgrad[2]; // this stores the average gradient energy
extern real_t edrho[2];  // this stores the avg. energy density over the box


/***********************
Main file Header
 ***********************/

// (re-)allocate memory for and initialize fields
void init(int nx, int ny, int nz);

/***********************
Initialization Header
 ***********************/
//These are the declerations of the initialization functions whose definitions are found in g2init.cpp

// initializes all of the energy densities and scale factor as appropriet for the begining of the run
void initexpansion();

// function which initializes the random conditions for fields f and df
void zeroInit( real_t * f, real_t * df );


/************
Model Header
 ************/
 //These are the declerations of the functions defined in g2model.cpp

// function to evaluate the potential of the field(s)
real_t potential(int s, int idx);

//function to store derivative wrt field of the potential
real_t dVdf(int s, int fld, int idx);

// function holds the effective mass of the fields, returns 1. if none is stored.
inline real_t effMass(int s, int fld);

// function to initialize the fields (and anything else)
void initfields();


/************
 Functions Header
 ************/
 
// There are the declerations of the functions defined in g2functions.cpp

//This function squares doubles
real_t pw2(real_t x);

//this is the function to call for the 7pt laplacian
real_t laplacian(real_t * f, int i, int j, int k);

//spatial derivative of a field in the i (x) direction
real_t dfdi(real_t * f, int i, int j, int k);

//spatial derivative of a field in the j (y) direction
real_t dfdj(real_t * f, int i, int j, int k);

//spatial derivative of a field in the j (y) direction
real_t dfdk(real_t * f, int i, int j, int k);

//spatial derivative of the field f in the "x" direction (stores the three functions above).
real_t dfdx(real_t * f, int x, int i, int j, int k);

//takes the gradient of the field at a point and squares it
real_t gradF2(real_t * f,int i,int j,int k);

//calculates the average gradient energy over the box
real_t avgGrad(int s);

//calculates the average potential energy over the box
real_t avgPot(int s);

// calculates the avereage kinetic energy over the box
real_t avgKin(int s);

// calculates the total average energy  over the box
void calcEnergy(int s);

// calculates adot from the average energy density
real_t adf(int s);

// equation of motion for the fields (klein gordon)
real_t ddfield(int s, int fld, int i, int j, int k);

//performs the full RK2 integration
void step();
