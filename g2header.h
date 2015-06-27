/**********************************
 Header File (declarations)
 **********************************/

/*
 This header file contains all the declarations for program functions.
 
 
 Copyright (2013):
 Kenyon College
 John T. Giblin, Jr
 Tate Deskins and Hillary Child
 Last Updated: 06.27.2013
 */


#include <stdio.h>//needed for printf fprintf
#include <stdlib.h>//needed for malloc
#include <math.h>//needed for fabs, tanh
#include <time.h>//needed for time and ctime
#include <fftw3.h>//needed for the rand initial conditions
#include <unistd.h>

#include <omp.h>// needed for omp

#define omp_set_nested(x)
#define omp_get_num_threads() 1

//this is a directive to declare all the common indices for the fields
#define DECLARE_INDEX gIdx fld,i,j,k;

//Indecies for function deceleration
#define INDECIES_f gIdx s, gIdx fld, gIdx i, gIdx j, gIdx k

//indecies for function call
#define INDECIES s, fld, i, j, k


/****** TYPE DEFINITION **********/
typedef long double gNum;//number (for fields) type
typedef long int gIdx;//index type

/******************************
 Global Preprocessor Directives
 ******************************/
//the following are definitions for looping over different field indices

#define fldLOOP for(fld=0;fld<nflds;fld++) for(i=0;i<N;i++) for(j=0;j<N;j++) for(k=0;k<N;k++)
//this loops over the fld and all three indices

#define LOOP for(i=0;i<N;i++) for(j=0;j<N;j++) for(k=0;k<N;k++)
//this loops over just the indices (3D)

#define LOOP2 for(j=0;j<N;j++) for(k=0;k<N;k++)
//this loops over just the j and k indices (2D)

#define INDEX(s,fld,i,j,k) N*(N*(N*(nflds*s+fld)+i)+j)+k//index function to make array format easier to adjust

#include "g2parameters.h"
//this is the include statement for the parameters file



/******************************
 global parameters DO NOT CHANGE
 ******************************/
extern gNum t;//this is the variable that stores the evolution time

extern gNum field[int_order*N*N*N*nflds];//this stores the field values for each step along the grid
extern gNum dfield[int_order*N*N*N*nflds];//this stores the derivative of the field for each step along the grid

//The following are all arrays of length two, one for each step of the RK2 integration
extern gNum a[int_order];//this stores the scale factor for each step
extern gNum adot[int_order];// this stores the time derivative of the scale factor
extern gNum edpot[int_order]; //this stores the average potential energy
extern gNum edkin[int_order]; //this stores the average kinetic energy
extern gNum edgrad[int_order]; //this stores the average gradient energy
extern gNum edrho[int_order]; // this stores the avg. energy density over the box



/***********************
 Initialization Header
 ***********************/
//These are the declerations of the initialization functions whose definitions are found in g2init.cpp


void initexpansion();// initializes all of the energy densities and scale factor as appropriate for the beginning of the run

void dftMemAlloc();// Allocates the memory for the dft for moving the random initialization in momentum space to configuration space

void randInit( gNum f[],gNum df[],gNum d2vdf2);//function which initializes the random conditions for fields f and df

void initDestroy();//destroys the fftw extra stuffs needed only during initialization


/************
 Model Header
 ************/
//These are the decelerations of the functions defined in g2model.cpp


void modelinfo(FILE *info);//function which prints model dependent information to info.txt

gNum potential(INDECIES_f);// function to evaluate the potential of the field(s)

gNum profile(gNum tin); //function to set how the initial energy density grows

gNum radius(gNum tin, gNum tdone, gNum rate); //function to set what the radius of the orbit is.

gNum massGaussian1 (INDECIES_f,gNum tin);//the mass Gaussian for one of the stars

gNum massGaussian2 (INDECIES_f,gNum tin);//the mass Gaussian for the other star

gNum energyDensity(INDECIES_f,gNum tin); //function defining the energy density of the source.

gNum galileon3(INDECIES_f,gNum tin); // function to numerically calculate pi-double-dot for the cubic Galileon

gNum galileon4(INDECIES_f,gNum tin); // function to numerically calculate pi-double-dot for the quartic Galileon

inline gNum effMass(INDECIES_f);// function holds the effective mass of the fields, returns 1. if none is stored.

void initfields();//function to initialize the fields (and anything else)


/************
 Functions Header
 ************/

// There are the decelerations of the functions defined in g2functions.cpp

gNum pw2(gNum x);//This function squares doubles

gIdx incr(gIdx i);///for incrementing with periodic boundary conditions

gIdx decr(gIdx i);//for decrementing with periodic boundary conditions

gNum laplacian(gNum f[], INDECIES_f);//this is the function to call for the 7pt laplacian

gNum dfdi(gNum f[], INDECIES_f);//spatial derivative of a field in the i (x) direction

gNum dfdj(gNum f[], INDECIES_f);//spatial derivative of a field in the j (y) direction

gNum dfdk(gNum f[], INDECIES_f);//spatial derivative of a field in the j (y) direction

gNum dfdi_s(gNum f[], INDECIES_f);// spatial derivative of the field f in i (x) direction direction one sided (check sign)

gNum dfdj_s(gNum f[], INDECIES_f);// spatial derivative of the field f in j (y) direction one sided (check sign)

gNum dfdk_s(gNum f[], INDECIES_f);// spatial derivative of the field f in k (z) direction one sided (check sign)

gNum dfdx(gNum f[], int x, INDECIES_f);//spatial derivative of the field f in the "x" direction (stores the three functions above).

gNum dfdii(gNum f[], INDECIES_f); //double spatial derivative in x or i direction.

gNum dfdjj(gNum f[], INDECIES_f); // double spatial derivative in y or j direction.

gNum dfdkk(gNum f[], INDECIES_f); // double spatial derivative in z or k direction.

gNum dfdij (gNum f[], INDECIES_f); //mixed partial wrt x and y (i and j).

gNum dfdjk(gNum f[], INDECIES_f); //mixed partial wrt y and z (j and k).

gNum dfdik(gNum f[], INDECIES_f); //mixed partial wrt x and z (i and k).

gNum dfdr_bound(gNum f[], INDECIES_f);// radial derivative for the boundary

gNum dfdr_analytic(INDECIES_f,gNum tin);// full analytic derivative.

gNum dfdr_pert(gNum f[], INDECIES_f, gNum tin);// radial derivative for the perturbation (on boundary)

gNum dfdr_bulk(gNum f[], INDECIES_f);//radial derivative for the bulk

gNum dfdr(gNum f[], INDECIES_f);//radial derivative

gIdx sign(gIdx x); //returns 1 for x>0, 0 for x=0 -1 for x<0 (signum)

gNum gradF2(gNum f[], INDECIES_f);//takes the gradient of the field at a point and squares it

gNum avgGrad(gIdx s);//calculates the average gradient energy over the box

gNum avgPot(gIdx s);//calculates the average potential energy over the box

gNum avgKin(gIdx s);// calculates the average kinetic energy over the box

void calcEnergy(gIdx s);// calculates the total average energy  over the box

gNum adf(gIdx s);// calculates adot from the average energy density

gNum ddfield(INDECIES_f,gNum tin);//equation of motion for the fields (Klein-Gordon)

void step();//performs the full RK2 integration

void steprk4();// performs the full RK4 integration

/**************
 Output Header
 *************/
//Declaration for all output functions found in g2output.cpp

gNum pi_powerout(gNum *bkgf);//Calculates the power of pi at the sides of he box.

void outputfield(int first);//outputs the field values over box (dimension and sampling determined in g2parameters.h

void screenout();//determines when there should be screen output

//gNum ENERGYDENISTY(int i, int j, int k);//calculates the energy density at each grid point

void meansvars(); //outputs means and variances of all fields

int slicewaitf();//evaluates the slicewait value (how long to wait between output slices) for the outputslice function

void outputslice();//this function determines what is output how and when

void output_parameters();//this creates the info.txt and populates it

void nancheck(); //this checks field[index] for nan.

void readable_time(int tt, FILE *info);//prints meaningful time to the info.txt


/**************
 Spectra Header
 **************/
//Deceleration for all functions found in g2spectra.cpp
void specOut(int first);//the calculates and prints to file the spectra of the fields

void specClear();//clears memory from the dft's used in specOut

