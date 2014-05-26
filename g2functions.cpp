/**********************************
 FUNCTIONS FILE
 **********************************/

/*
 This header file contains all the functions which are independent of the model needed to evolve the fields by the Runge-Kutta Second order method (for first order finite derivatives). The incr and decr commands set periodic boundary conditions on the lattice durring evolution.
 
 Copyright (2013): 
 Kenyon College
 John T. Giblin, Jr
 Tate Deskins and Hillary Child
 Last Updated: 06.27.2013
*/

#include "g2header.h" //contains declerations for program functions.

long double pw2(long double x)//squares long doubles
{
    return x*x;
}

inline int incr(int i)//for periodic boundaries
{
    return i == N-1? 0: i+1;
}

inline int decr(int i)//for periodic boundaries
{
    return i == 0? N-1: i-1;
}

/** Laplacian Functions **/


long double laplacian(long double f[][N][N], int i, int j, int k)//this calculates the seven point laplacian 
{
    
    return (
            f[incr(i)][j][k] + f[decr(i)][j][k]
            + f[i][incr(j)][k] + f[i][decr(j)][k]
            + f[i][j][incr(k)] + f[i][j][decr(k)]
            - 6.*f[i][j][k]
        )/(dx*dx);
    
}

/** Spatial Derivative Functions **/
/**
 * The following fucntions are used to calculate the
 * gradient energy of the field only; but they can also
 * be implemented if derivative couplings are desired
 **/

//If you know which partial derivative you need

long double dfdi(long double f[][N][N], int i, int j, int k)// spatial derivative of the field f in i (x) direction direction
{
    return (f[incr(i)][j][k] - f[decr(i)][j][k])/2./dx;
}

long double dfdj(long double f[][N][N], int i, int j, int k)// spatial derivative of the field f in j (y) direction
{
    return (f[i][incr(j)][k] - f[i][decr(j)][k])/2./dx;
}

long double dfdk(long double f[][N][N], int i, int j, int k)// spatial derivative of the field f in k (z) direction
{
    return (f[i][j][incr(k)] - f[i][j][decr(k)])/2./dx;
}

// If you want to loop over spatial derivatives this form is somewhat more convienent
long double dfdx(long double f[][N][N], int x, int i, int j, int k)//spatial derivative of the field f in the "x" direction.
{
    switch (x)
    {
        case 0:
            return dfdi(f,i,j,k);
        case 1:
            return dfdj(f,i,j,k);
        case 2:
            return dfdk(f,i,j,k);
        default:
            return 0;
    }
}

long double ddfield( int s, int fld, int i, int j, int k)//evaluates the double time derivative of the field fld (s) at i,j,k.
{
   return laplacian(field[s][fld],i,j,k);
}


/** RK2 function **/

void step()//this steps (integrates) the field and it time derivative via the rk2 meathod.
{
    
    DECLARE_INDEX

    //the first part of the RK2 step
    for(fld = 0; fld<nflds; fld++)
    {
        //paralleleizes over the index i
        #pragma omp parallel for private (j,k) num_threads (tot_num_thrds)
        LOOP//loops over fld i,j,k
        {
            field[1][fld][i][j][k] = field[0][fld][i][j][k] + .5*dt*dfield[0][fld][i][j][k];
            dfield[1][fld][i][j][k] = dfield[0][fld][i][j][k] + .5*dt*ddfield(0,fld,i,j,k);
        }
    }

    // the second part of the RK2 step
    for(fld = 0; fld<nflds; fld++)
    {
        //paralleleizes over the index i
        #pragma omp parallel for private (j,k) num_threads (tot_num_thrds)
        LOOP//This returns the actuall value of the field and derivative at t
        {
            field[0][fld][i][j][k] = field[0][fld][i][j][k] + dt*dfield[1][fld][i][j][k];
            dfield[0][fld][i][j][k] = dfield[0][fld][i][j][k] + dt*ddfield(1,fld,i,j,k);
        }
    }	

}
