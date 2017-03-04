/**********************************
 Initialization FILE
 **********************************/

/*
 This header file contains all the functions which are independent of the model needed to initialize the fields and random initial conditions.
 
 Copyright Kenyon College
 John T. Giblin, Jr
 Tate Deskins and Hillary Child
*/

#include "g2header.h" //contains declerations for program functions.

/** expansion initialization **/
void initexpansion()
{
    //for power law
    a[0]=1;
    adot[0]=0;
    a[1]=1;
    adot[1]=0;
    calcEnergy(0);
}

void randInit(real_t f[][N][N], real_t df[][N][N], real_t d2vdf2)
{
    // TODO: initialize something
}
