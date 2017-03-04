/**********************************
 FUNCTIONS FILE
 **********************************/

/*
 This header file contains all the functions which are independent of the
 model needed to evolve the fields by the Runge-Kutta Second order method
 (for first order finite derivatives). The modular arithmetic sets periodic
 boundary conditions on the lattice durring evolution.
 
 Copyright (2013): 
 Kenyon College
 John T. Giblin, Jr
 Tate Deskins and Hillary Child
 Last Updated: 06.27.2013
*/

#include "g2header.h" //contains declerations for program functions.

real_t pw2(real_t x) // squares real_ts
{
    return x*x;
}

/** Laplacian Function **/
real_t laplacian(real_t * f, int i, int j, int k)//this calculates the seven point laplacian 
{
    return (
        f[IDX((i+1)%NX, j, k)] + f[IDX((i+NX-1)%NX, j, k)]
        + f[IDX(i, (j+1)%NY, k)] + f[IDX(i, (j+NY-1)%NY, k)]
        + f[IDX(i, j, (k+1)%NZ)] + f[IDX(i, j, (k+NZ-1)%NZ)]
        - 6.*f[IDX(i,j,k)]
    )/(dx*dx);
}

/** Spatial Derivative Functions
 * The following fucntions are used to calculate the gradient energy of the
 * field only; but they can also be used if derivative couplings are desired
 */

// If you know which partial derivative you need
real_t dfdi(real_t * f, int i, int j, int k)// spatial derivative of the field f in i (x) direction direction
{
    return (f[IDX((i+1)%NX,j,k)] - f[IDX((i+NX-1)%NX,j,k)])/2./dx;
}

real_t dfdj(real_t * f, int i, int j, int k)// spatial derivative of the field f in j (y) direction
{
    return (f[IDX(i,(j+1)%NY,k)] - f[IDX(i, (j+NY-1)%NY, k)])/2./dx;
}

real_t dfdk(real_t * f, int i, int j, int k)// spatial derivative of the field f in k (z) direction
{
    return (f[IDX(i,j,(k+1)%NZ)] - f[IDX(i, j, (k+NZ-1)%NZ)])/2./dx;
}

// If you want to loop over spatial derivatives this form is somewhat more convienent
// spatial derivative of the field f in the "x" direction.
real_t dfdx(real_t * f, int x, int i, int j, int k)
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

/**Functions needed for self-consistent expansion**/

// this is the unscaled gradient of the field at the point i,j,k
real_t gradF2(real_t * f, int i, int j, int k)
{
    return dfdi(f,i,j,k)*dfdi(f,i,j,k) + dfdj(f,i,j,k)*dfdj(f,i,j,k)
     + dfdk(f,i,j,k)*dfdk(f,i,j,k);
}

// Find the average gradient energy
real_t avgGrad(int s)
{
    real_t grad=0;
    DECLARE_INDEX
    // loops over fld, i,j,k
    fldLOOP
    {
        // sums the gradient energy at each point
        grad+=gradF2(field[FIELD(s,fld)],i,j,k);
    }

    // divides by the gridsize (to normalize) and 1/(2a^2) to get the gradient energy density
    return grad/gridsize/2./a[s]/a[s];
}

// Find the average potential energy
real_t avgPot(int s)
{
    real_t pot=0;
    DECLARE_INDEX

    for(int p=0; p<POINTS; p++)
    {
        // sums the potential at every point
        pot += potential(s,p);
    }
    
    // averages over the grid
    return pot/gridsize;
}


// Find the average kinetic energy
real_t avgKin(int s)
{
    real_t kin=0;
    DECLARE_INDEX
    //loops over i,j,k
    for(int fld=0; fld<nflds; fld++)
        for(int p=0; p<POINTS; p++)
    {
        //sums the square field derivative at every point
        kin += dfield[FIELD(s,fld)][p]*dfield[FIELD(s,fld)][p];
    }
    // divide by the grid size to get the average and 2
    return kin/gridsize/2.;
}

// Calculate the total energy
void calcEnergy(int s)
{
    edkin[s] = avgKin(s);
    edpot[s] = avgPot(s);
    edgrad[s] = avgGrad(s);
    edrho[s] = edkin[s]+edpot[s]+edgrad[s];
}

// the friedman equation
real_t adf(int s)
{
    return sqrt(8.*M_PI*grav/3.*edrho[s])*a[s];
}

// evaluates the double time derivative of the field fld (s) at i,j,k.
real_t ddfield( int s, int fld, int i, int j, int k)
{
    int idx = IDX(i,j,k);
    return laplacian(field[FIELD(s,fld)],i,j,k)/a[s]/a[s]
        - dVdf(s,fld,idx) - 3*adot[s]/a[s]*dfield[FIELD(s,fld)][idx];
}


/** RK2 function **/

// this steps (integrates) the field and it time derivative via the rk2 method.
void step()
{
    
    DECLARE_INDEX

    // the first part of the RK2 step
    for(fld=0; fld<nflds; fld++)
    {
        int flds0 = FIELD(0,fld);
        int flds1 = FIELD(1,fld);
        for(int p=0; p<POINTS; p++)
        {
            field[flds1][p] = field[flds0][p] + .5*dt*dfield[flds0][p];
        }
        // loops over fld i,j,k
        LOOP
        {
            int idx = IDX(i,j,k);
            dfield[flds1][idx] = dfield[flds0][idx] + .5*dt*ddfield(0,fld,i,j,k);
        }
    }


    // the second part of the RK2 step
    // This computes the actual value of the field and derivative at t
    for(fld=0; fld<nflds; fld++)
    {
        int flds0 = FIELD(0,fld);
        int flds1 = FIELD(1,fld);
        for(int p=0; p<POINTS; p++)
        {
            field[flds0][p] = field[flds0][p] + dt*dfield[flds1][p];
        }
        LOOP
        {
            int idx = IDX(i,j,k);
            dfield[flds0][idx] = dfield[flds0][idx] + dt*ddfield(1,fld,i,j,k);
        }
    }

    // calcualtes the energy at the end of the step.
    calcEnergy(0);
}
