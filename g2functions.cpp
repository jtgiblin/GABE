/**********************************
 FUNCTIONS FILE
 **********************************/

/*
 This header file contains all the functions which are independent of the model needed to evolve the fields by the Runge-Kutta Second order method (for first order finite derivatives). The incr and decr commands set periodic boundary conditions on the lattice during evolution.
 
 Copyright (2013): 
 Kenyon College
 John T. Giblin, Jr
 Tate Deskins and Hillary Child
 Last Updated: 06.27.2013
*/

#include "g2header.h" //contains declarations for program functions.

gNum pw2(gNum x)//squares gNums
{
    return x*x;
}

inline int incr(int i)//for periodic boundaries
{
    return i==N-1? 0: i+1;
}

inline int decr(int i)//for periodic boundaries
{
    return i==0? N-1: i-1;
}

/** Laplacian Functions **/


gNum laplacian(gNum f[][N][N], int i, int j, int k)//this calculates the seven point laplacian 
{
    
    return (f[incr(i)][j][k]+f[decr(i)][j][k]
            +f[i][incr(j)][k]+f[i][decr(j)][k]
            +f[i][j][incr(k)]+f[i][j][decr(k)]
            -6.*f[i][j][k])/(dx*dx);
    
}

/** Spatial Derivative Functions **/
/*
The following functions are used to calculate the gradient energy of the field only; but they can also be implemented if derivative couplings are desired*/

//If you know which partial derivative you need

/*
//First-order spatial derivatives
gNum dfdi(gNum f[][N][N], int i, int j, int k)// spatial derivative of the field f in i (x) direction direction
{
    return (f[incr(i)][j][k]-f[decr(i)][j][k])/2./dx;
}

gNum dfdj(gNum f[][N][N], int i, int j, int k)// spatial derivative of the field f in j (y) direction
{
    return (f[i][incr(j)][k]-f[i][decr(j)][k])/2./dx;
}

gNum dfdk(gNum f[][N][N], int i, int j, int k)// spatial derivative of the field f in k (z) direction
{
    return (f[i][j][incr(k)]-f[i][j][decr(k)])/2./dx;
}
 */

//Second-order spatial derivatives
gNum dfdi(gNum f[][N][N], int i, int j, int k)
{
    return (
        (1. / 12.) * f[decr(decr(i))][j][k]
        + ((-2.) / 3.) * f[decr(i)][j][k]
        + (2. / 3.) * f[incr(i)][j][k]
        + ((-1.) / 12.) * f[incr(incr(i))][j][k]
    ) / dx;
}

gNum dfdj(gNum f[][N][N], int i, int j, int k)
{
    return (
        (1. / 12.) * f[i][decr(decr(j))][k]
        + ((-2.) / 3.) * f[i][decr(j)][k]
        + (2. / 3.) * f[i][incr(j)][k]
        + ((-1.) / 12.) * f[i][incr(incr(j))][k]
    ) / dx;
}

gNum dfdk(gNum f[][N][N], int i, int j, int k)
{
    return (
        (1. / 12.) * f[i][j][decr(decr(k))]
        + ((-2.) / 3.) * f[i][j][decr(k)]
        + (2. / 3.) * f[i][j][incr(k)]
        + ((-1.) / 12.) * f[i][j][incr(incr(k))]
    ) / dx;
}

//If you want to loop over spatial derivatives this form is somewhat more convenient
gNum dfdx(gNum f[][N][N], int x, int i, int j, int k)//spatial derivative of the field f in the "x" direction.
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

gNum gradF2(gNum f[][N][N],int i,int j,int k){
    
    return  dfdi(f,i,j,k)*dfdi(f,i,j,k)+dfdj(f,i,j,k)*dfdj(f,i,j,k)+dfdk(f,i,j,k)*dfdk(f,i,j,k);//this is the unscaled gradient fo the field at the point i,j,k
    
}

gNum avgGrad(int s) //Find the average gradient energy
{
    gNum grad=0;
    DECLARE_INDEX
    for(fld=0; fld<nflds; fld++){
#pragma omp parallel for private (j,k) reduction(+:grad) num_threads (tot_num_thrds)
        LOOP//loops over i,j,k
        {
            grad-=field[s][fld][i][j][k]*laplacian(field[s][fld],i,j,k);//sums the gradient energy at each point
        }
    }
    return grad/gridsize/2./a[s]/a[s];//divides by the gridsize (to normalize) and 1/(2a^2) to get the gradient energy density
}


gNum avgPot(int s) //Find the average potential energy
{
    gNum pot=0;
    DECLARE_INDEX
    #pragma omp parallel for private (j,k) reduction(+:pot) num_threads (tot_num_thrds)
        LOOP//loops over i,j,k
        {
            pot+=potential(s,i,j,k);//sums the potential at every point
        }
    
    return pot/gridsize;//averages over the grid
}


gNum avgKin(int s) //Find the average kinetic energy
{
    gNum kin=0;
    DECLARE_INDEX
    for(fld=0; fld<nflds; fld++){
#pragma omp parallel for private (j,k) reduction(+:kin) num_threads (tot_num_thrds)
        LOOP//loops over i,j,k
        {
            kin+=dfield[s][fld][i][j][k]*dfield[s][fld][i][j][k];//sums the square field derivative at every point
        }
    }
    return kin/gridsize/2.;//divide by the grid size to get the average and 2
}


void calcEnergy(int s) //Calculate the total energy
{
    edkin[s]=avgKin(s);
    edpot[s]=avgPot(s);
    edgrad[s]=avgGrad(s);
    edrho[s]=edkin[s]+edpot[s]+edgrad[s];
}

gNum calcrho(int s, int i, int j, int k){ //calculates rho at a point
    
    gNum rho = potential(s,i,j,k);
    for (int fld = 0; fld<nflds; fld++){
        rho += (.5*dfield[s][fld][i][j][k]*dfield[s][fld][i][j][k] +
                .5*gradF2(field[0][fld],i,j,k)/a[s]/a[s]);
    }
    
    return rho;
            
}

gNum adf(int s)//the friedman equation
{
    return sqrt(8.*M_PI/3.*edrho[s])*a[s];
}

gNum ddfield( int s, int fld, int i, int j, int k)//evaluates the double time derivative of the field fld (s) at i,j,k.
{
   return laplacian(field[s][fld],i,j,k)/a[s]/a[s] - dVdf(s,fld,i,j,k) - 3*adot[s]/a[s]*dfield[s][fld][i][j][k];
}


/** RK2 function **/

void step()//this steps (integrates) the field and it time derivative via the rk2 meathod.
{
    
    DECLARE_INDEX

	
#if expansion_type==0
    //no expansion note that this only calculates the energies at the end of the full step

    for(fld=0;fld<nflds;fld++)//the first part of the RK2 step 
    {
    //paralleleizes over the index i
#pragma omp parallel for private (j,k) num_threads (tot_num_thrds)
        LOOP//loops over fld i,j,k
        {
            field[1][fld][i][j][k]=field[0][fld][i][j][k]+.5*dt*dfield[0][fld][i][j][k];
            dfield[1][fld][i][j][k]=dfield[0][fld][i][j][k]+.5*dt*ddfield(0,fld,i,j,k);
        }
    }
	
	for(fld=0;fld<nflds;fld++)//the second part of the RK2 step 
    {
    //paralleleizes over the index i
#pragma omp parallel for private (j,k) num_threads (tot_num_thrds)
        LOOP//This returns the actual value of the field and derivative at t
        {
            field[0][fld][i][j][k]=field[0][fld][i][j][k]+dt*dfield[1][fld][i][j][k];
            dfield[0][fld][i][j][k]=dfield[0][fld][i][j][k]+dt*ddfield(1,fld,i,j,k);
        }
    }	
    
    calcEnergy(0);//calculates the energy at the end of the step.
	
#elif expansion_type==1
	
    for(fld=0;fld<nflds;fld++)//first step of the Rk2 integration
    {
#pragma omp parallel for private (j,k) num_threads (tot_num_thrds)
        LOOP//loops over fld i,j,k
        {
            field[1][fld][i][j][k]=field[0][fld][i][j][k]+.5*dt*dfield[0][fld][i][j][k];
            dfield[1][fld][i][j][k]=dfield[0][fld][i][j][k]+.5*dt*ddfield(0,fld,i,j,k);
        }
    }
	
#if calc_gws == 1
    evolve_perts(0, 1, 0.5 * dt); //evolve tensor modes
#endif
    
    a[1]=a[0]+.5*dt*adot[0];//this does the first step of the RK2 for the scale factor
    calcEnergy(1);//this calculates the energy based on this half step
    adot[1]=adf(1);//this updates adot based off of the energy at this step
    
    	for(fld=0;fld<nflds;fld++)//second step of the Rk2 integration
    {
#pragma omp parallel for private (j,k) num_threads (tot_num_thrds)
        LOOP//This returns the actuall value of the field and derivative at t
        {
            field[0][fld][i][j][k]=field[0][fld][i][j][k]+dt*dfield[1][fld][i][j][k];
            dfield[0][fld][i][j][k]=dfield[0][fld][i][j][k]+dt*ddfield(1,fld,i,j,k);
        }
    }
    
#if calc_gws == 1
    evolve_perts(1, 0, dt); //evolve tensor modes
#endif
    
	a[0]=a[0]+dt*adot[1];//this calculates the full step scale factor
    calcEnergy(0);//calculates the energy at the full step
    adot[0]=adf(0);//then calculates adot based off of the full step
    
#elif expansion_type==2
	
    for(fld=0;fld<nflds;fld++)//first step of the Rk2 integration
    {
#pragma omp parallel for private (j,k) num_threads (tot_num_thrds)
        LOOP//loops over fld i,j,k
        {
            field[1][fld][i][j][k]=field[0][fld][i][j][k]+.5*dt*dfield[0][fld][i][j][k];
            dfield[1][fld][i][j][k]=dfield[0][fld][i][j][k]+.5*dt*ddfield(0,fld,i,j,k);
        }
    }
	
    /* this may need to chage based off of user defined expansion*/
    a[1]=a[0]+.5*dt*adot[0];//this does the first step of the RK2 for the scale factor
    calcEnergy(1);//this calculates the energy based on this half step
    adot[1]=adf(1);//this updates adot based off of the energy at this step
    
    	for(fld=0;fld<nflds;fld++)//second step of the Rk2 integration
    {
#pragma omp parallel for private (j,k) num_threads (tot_num_thrds)
        LOOP//This returns the actuall value of the field and derivative at t
        {
            field[0][fld][i][j][k]=field[0][fld][i][j][k]+dt*dfield[1][fld][i][j][k];
            dfield[0][fld][i][j][k]=dfield[0][fld][i][j][k]+dt*ddfield(1,fld,i,j,k);
        }
    }	
    /* this may need to change based off of user defined expansion*/
	a[0]=a[0]+dt*adot[1];//this calculates the full step scale factor
    calcEnergy(0);//calculates the energy at the full step
    adot[0]=adf(0);//then calculates adot based off of the full step
#endif
    
}
