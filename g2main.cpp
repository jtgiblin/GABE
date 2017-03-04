/****************************
 LatticeSolve Main File
 ****************************/

/* 
 Copyright (2013): 
 Kenyon College
 John T. Giblin, Jr
 Tate Deskins and Hillary Child
 Last Updated: 06.27.2013
*/

/*********************
 GABE Header
 *********************/

//contains all the parameters, changable and unchangeable
#include "g2header.h"
#include <iostream>

real_t (*field)[nflds][NX][NY][NZ];  // this stores the field values for each step along the grid
real_t (*dfield)[nflds][NX][NY][NZ]; // this stores the derivative of the field for each step along the grid
real_t a[2];                         // this stores the scale facator for each step
real_t adot[2];                      // this stores the time derivative of the scale factor
real_t edpot[2];                     // this stores the average potential energy
real_t edkin[2];                     // this stores the average kinetic energy
real_t edgrad[2];                    // this stores the average gradient energy
real_t edrho[2];                     // this stores the avg. energy density over the box

// values of model constants from g2parameters.h
real_t mphi = 1.e-6;            // mass of phi field
real_t phi0 = 0.193;            // initial avg phi field value
real_t gsq = 2.5e-5;            // g^2 value for phi chi coupling
real_t f0[2] = {phi0,0.};       // array storing intial phi and chi field values
real_t df0[2] = {-0.142231,0.}; // array storing intial phi and chi field derivative values
real_t grav = 1.;               // the gravitational constant
real_t c = 1.;                  // speed of light, shouldn't change, if you do, fix all equations

// values of model-independent parameters from g2parameters.h
real_t L = 20.;         //  length of one side of box in prgm units
real_t starttime = 0.;  // start time of simulation
real_t endtime = 10.;   // end time of simulations
real_t dt = 0.1;        // time step size


// Allocate space for fields. Other functions will be called via javascript.
void alloc()
{
    //allocates memory for the fields
    field = (real_t(*)[nflds][NX][NY][NZ]) malloc(sizeof(real_t)*nflds*2*NX*NY*NZ);
    
    //allocates memory for the fields' time derivatives
    dfield = (real_t(*)[nflds][NX][NY][NZ]) malloc(sizeof(real_t)*nflds*2*NX*NY*NZ);
}

// Just initialize things. 
// Steps will be called separately (by javascript routines).
void init()
{
    initfields(); //initializes fields as defined in model.h
    initexpansion(); //start expansion;
    initfields(); //do fluctuations
}

int main()
{
    alloc();
    init();
    
    // This is the main for-loop over time.
    // Note that the current t value is the time for the currently evaluated fields
    for(real_t t=starttime; t<=endtime; t+=dt)
    {
        // evolves the fields one step
        step();
    }

    return 0;
}

#ifdef EMSCRIPTEN_BINDINGS
EMSCRIPTEN_BINDINGS(my_module) {
    function("init", &init);
    function("alloc", &alloc);
}
#endif
