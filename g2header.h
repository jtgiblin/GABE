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
#include <unistd.h>

#include <omp.h>// needed for omp

#define omp_set_nested(x)
#define omp_get_num_threads() 1

//this is a directive to declare all the common indecies for the fields
#define DECLARE_INDEX int fld,i,j,k; 

/****************
LOOP DEFINITIONS
 ****************/
//the following are definitions for looping over different field indicies

#define fldLOOP for(fld=0;fld<nflds;fld++) for(i=0;i<N;i++) for(j=0;j<N;j++) for(k=0;k<N;k++)
//this loops over the fld and all three indicies

#define LOOP for(i=0;i<N;i++) for(j=0;j<N;j++) for(k=0;k<N;k++)
//this loops over just the indicies (3D)

#define LOOP2 for(j=0;j<N;j++) for(k=0;k<N;k++)
//this loops over just the j and k indicies (2d)

#include "g2parameters.h"
//this is the include statement for the parameters file

/******************************
 global parameters DO NOT CHANGE
 ******************************/
extern long double t; //this is the variable that stores the evolution time

extern long double (* field)[nflds][N][N][N]; //this stores the field values for each step along the grid
extern long double (* dfield)[nflds][N][N][N]; //this stores the derivative of the field for each step along the grid


/************
Model Header
 ************/
 //These are the declerations of the functions defined in g2model.cpp

void modelinfo(FILE *info); //function which prints model dependent information to info.txt

void initfields(); //function to initialize the fields (and anything else)


/************
 Functions Header
 ************/
 
// There are the declerations of the functions defined in g2functions.cpp

long double pw2(long double x); //This function squares doubles

int incr(int i); ///for incremiting with periodic boundary conditions

int decr(int i); //for decremiting with periodic boundary conditions

long double laplacian(long double f[][N][N], int i, int j, int k); //this is the function to call for the 7pt laplacian

long double dfdi(long double f[][N][N], int i, int j, int k); //spatial derivative of a field in the i (x) direction

long double dfdj(long double f[][N][N], int i, int j, int k); //spatial derivative of a field in the j (y) direction

long double dfdk(long double f[][N][N], int i, int j, int k); //spatial derivative of a field in the j (y) direction

long double dfdx(long double f[][N][N], int x, int i, int j, int k); //spatial derivative of the field f in the "x" direction (stores the three functions above).

long double ddfield(int s, int fld, int i, int j, int k); //equation of motion for the fields (klein gordon)

void step(); //performs the full RK2 integration



/**************
 Output Header
 *************/
 //Decleartion for all output functions found in g2output.cpp

void outputfield(int first); //outputs the field values over box (dimension and sampling determined in g2parameters.h

int slicewaitf(); //evaluates the slicewait value (how long to wait between output slices) for the outputslice function

void outputslice(); //this function determines what is output how and when

void output_parameters(); //this creates the info.txt and populates it

void readable_time(int tt, FILE *info); //prints meanigful time to the info.txt

void screenout(); //determines when there should be screen output

void meansvars(); //outputs means and variances of all fields
