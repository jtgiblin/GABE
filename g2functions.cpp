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

#include "g2header.h" //contains decelerations for program functions.

gNum pw2(gNum x)//squares long doubles
{
    return x*x;
}

inline gIdx incr(gIdx i)//for damping BC
{
    return i==N-1? N-1: i+1;
}

inline gIdx decr(gIdx i)//for damping BC
{
    return i==0? 0: i-1;
}

/** Laplacian Functions **/


gNum laplacian(gNum f[], INDECIES_f)//this calculates the seven point laplacian
{
    
    return (f[INDEX(s,fld,incr(i),j,k)]+f[INDEX(s,fld,decr(i),j,k)]
            +f[INDEX(s,fld,i,incr(j),k)]+f[INDEX(s,fld,i,decr(j),k)]
            +f[INDEX(s,fld,i,j,incr(k))]+f[INDEX(s,fld,i,j,decr(k))]
            -6.*f[INDEX(s,fld,i,j,k)])/(dx*dx);
    
}

/** Spatial Derivative Functions **/
/*
 The following functions are used to calculate the gradient energy of the field only; but they can also be implemented if derivative couplings are desired*/

//If you know which partial derivative you need

gNum dfdi(gNum f[], INDECIES_f)// spatial derivative of the field f in i (x) direction direction
{
    return (f[INDEX(s,fld,incr(i),j,k)]-f[INDEX(s,fld,decr(i),j,k)])/2./dx;
}

gNum dfdj(gNum f[], INDECIES_f)// spatial derivative of the field f in j (y) direction
{
    return (f[INDEX(s,fld,i,incr(j),k)]-f[INDEX(s,fld,i,decr(j),k)])/2./dx;
}

gNum dfdk(gNum f[], INDECIES_f)// spatial derivative of the field f in k (z) direction
{
    return (f[INDEX(s,fld,i,j,incr(k))]-f[INDEX(s,fld,i,j,decr(k))])/2./dx;
}


//one sided derivative (the sign on these may or may not be right so either square or think hard before use)

gNum dfdi_s(gNum f[], INDECIES_f)// spatial derivative of the field f in i (x) direction direction one sided (check sign)
{
    return sign(i-N/2)*(f[INDEX(s,fld,i,j,k)]-f[INDEX(s,fld,i-sign(i-N/2),j,k)])/dx;
}

gNum dfdj_s(gNum f[], INDECIES_f)// spatial derivative of the field f in j (y) direction one sided (check sign)
{
    return sign(j-N/2)*(f[INDEX(s,fld,i,j,k)]-f[INDEX(s,fld,i,j-sign(j-N/2),k)])/dx;
}

gNum dfdk_s(gNum f[], INDECIES_f)// spatial derivative of the field f in k (z) direction one sided (check sign)
{
    return sign(k-N/2)*(f[INDEX(s,fld,i,j,k)]-f[INDEX(s,fld,i,j,k-sign(k-N/2))])/dx;
}

//If you want to loop over spatial derivatives this form is somewhat more convenient
gNum dfdx(gNum f[], int x, INDECIES_f)//spatial derivative of the field f in the "x" direction.
{
    switch (x)
    {
        case 0:
            return dfdi(f,s,fld,i,j,k);
        case 1:
            return dfdj(f,s,fld,i,j,k);
        case 2:
            return dfdk(f,s,fld,i,j,k);
        default:
            return 0;
    }
}

gNum dfdii(gNum f[], INDECIES_f)
{
    return (-2.*f[INDEX(s,fld,i,j,k)] + f[INDEX(s,fld,incr(i),j,k)] +f[INDEX(s,fld,decr(i),j,k)])/(dx*dx);
}

gNum dfdjj(gNum f[], INDECIES_f)
{
    return (-2.*f[INDEX(s,fld,i,j,k)] +f[INDEX(s,fld,i,incr(j),k)] +f[INDEX(s,fld,i,decr(j),k)])/(dx*dx);
}

gNum dfdkk(gNum f[], INDECIES_f)
{
    return (-2.*f[INDEX(s,fld,i,j,k)] +f[INDEX(s,fld,i,j,incr(k))] +f[INDEX(s,fld,i,j,decr(k))])/(dx*dx);
}

/** Mixed Partial Spatial Derivatives **/
/*
 The following functions will be used in the equations of motion of the Galileon (or anywhere else if one would need them) */

//Realize that in these functions the error is proportional to dx*dx

gNum dfdij (gNum f[], INDECIES_f) // mixed partial wrt x and y (i and j)
{
    return (f[INDEX(s,fld,incr(i),incr(j),k)] - f[INDEX(s,fld,incr(i),decr(j),k)] - f[INDEX(s,fld,decr(i),incr(j),k)] + f[INDEX(s,fld,decr(i),decr(j),k)])/(4.*dx*dx);
}

gNum dfdjk(gNum f[], INDECIES_f)
{
    return (f[INDEX(s,fld,i,incr(j),incr(k))]-f[INDEX(s,fld,i,incr(j),decr(k))]-f[INDEX(s,fld,i,decr(j),incr(k))]+f[INDEX(s,fld,i,decr(j),decr(k))])/(4.*dx*dx);
}

gNum dfdik(gNum f[], INDECIES_f)
{
    return (f[INDEX(s,fld,incr(i),j,incr(k))]-f[INDEX(s,fld,incr(i),j,decr(k))]-f[INDEX(s,fld,decr(i),j,incr(k))]+f[INDEX(s,fld,decr(i),j,decr(k))])/(4.*dx*dx);
}


/*
 one-sided r deriv for boundaries above monopole background
 NOTE the analytic term may/will change with rescaling.
 NOTE the analytic term may need to be different for GR.
 */
gNum dfdr_bound(gNum f[], INDECIES_f)
{
    return (abs(i-N/2.)*(f[INDEX(s,fld,i,j,k)]-f[INDEX(s,fld,i-sign(i-N/2),j,k)])
            +abs(j-N/2.)*(f[INDEX(s,fld,i,j,k)]-f[INDEX(s,fld,i,j-sign(j-N/2),k)])
            +abs(k-N/2.)*(f[INDEX(s,fld,i,j,k)]-f[INDEX(s,fld,i,j,k-sign(k-N/2))]))
    /dx/sqrt(pw2(i-N/2.)+pw2(j-N/2.)+pw2(k-N/2.));//normal one-sided derivative
}//discrete derivative problems?

gNum dfdr_analytic(INDECIES_f, gNum tin){
    
    gNum r2=dx*dx*(pw2(i-N/2.)+pw2(j-N/2.)+pw2(k-N/2.));
    
    return 1./(4*sqrtl(r2)*kappa)*(-3.*r2+sqrtl(32.*pw2(rstar)*rstar*profile(tin)*sqrtl(r2)/M_PI+9.*pw2(r2)));//profile is because rs\propto m which is tanhed.
}

gNum dfdr_pert(gNum f[], INDECIES_f, gNum tin)//need to double check the form of this with the new solution
{
    
    return dfdr_bound(f,s,fld,i,j,k) - dfdr_analytic(s,fld,i,j,k,tin);//profile(t)/(4*M_PI*dx*dx*(pw2(i-N/2.)+pw2(j-N/2.)+pw2(k-N/2.)));//analytic term goes as M/4/Pi/R^2 where M is half the mass (double check with new)
    
}

gNum dfdr_bulk(gNum f[], INDECIES_f)
{
    return ((i-N/2.)*dfdi(f,s,fld,i,j,k)+ (j-N/2.)*dfdj(f,s,fld,i,j,k)+ (k-N/2.)*dfdk(f,s,fld,i,j,k))/sqrtl(pw2(i-N/2.)+pw2(j-N/2.)+pw2(k-N/2.));
}

gNum dfdr(gNum f[], INDECIES_f)
{
    if(i==0||j==0||k==0||i==N-1||j==N-1||k==N-1)
        return dfdr_bound(f,s,fld,i,j,k);
    else
        return dfdr_bulk(f,s,fld,i,j,k);
}


gIdx sign(gIdx x)
{
    return (0<x)-(x<0);
}



/**Functions needed for self-consistent expansion**/

gNum gradF2(gNum f[], INDECIES_f){
    
    return  dfdi(f,s,fld,i,j,k)*dfdi(f,s,fld,i,j,k)+dfdj(f,s,fld,i,j,k)*dfdj(f,s,fld,i,j,k)+dfdk(f,s,fld,i,j,k)*dfdk(f,s,fld,i,j,k);//this is the unscaled gradient of the field at the point i,j,k
    
}

gNum avgGrad(gIdx s) //Find the average gradient energy
{
    gNum grad=0.;
    DECLARE_INDEX
    for(fld=0; fld<nflds; fld++){
#pragma omp parallel for private (j,k) reduction(+:grad) num_threads (tot_num_thrds)
        LOOP//loops over i,j,k
        {
            grad+=gradF2(field,INDECIES);//sums the gradient energy at each point
        }
    }
    return grad/gridsize/2./a[s]/a[s];//divides by the grid-size (to normalize) and 1/(2a^2) to get the gradient energy density
}


gNum avgPot(gIdx s) //Find the average potential energy
{
    gNum pot=0;
    DECLARE_INDEX
#pragma omp parallel for private (j,k) reduction(+:pot) num_threads (tot_num_thrds)
    LOOP//loops over i,j,k
    {
        pot+=potential(INDECIES);//sums the potential at every point
    }
    
    return pot/gridsize;//averages over the grid
}


gNum avgKin(gIdx s) //Find the average kinetic energy
{
    gNum kin=0.;
    DECLARE_INDEX
    for(fld=0; fld<nflds; fld++){
#pragma omp parallel for private (j,k) reduction(+:kin) num_threads (tot_num_thrds)
        LOOP//loops over i,j,k
        {
            kin+=dfield[INDEX(s,fld,i,j,k)]*dfield[INDEX(s,fld,i,j,k)];//sums the square field derivative at every point
        }
    }
    return kin/gridsize/2.;//divide by the grid size to get the average and 2
}


void calcEnergy(gIdx s) //Calculate the total energy
{
    edkin[s]=avgKin(s);
    edpot[s]=avgPot(s);
    edgrad[s]=avgGrad(s);
    edrho[s]=edkin[s]+edpot[s]+edgrad[s];
}

gNum adf(gIdx s)//the Friedman equation
{
    return 0;//return sqrt(8.*M_PI*grav/3.*edrho[s])*a[s];
}



gNum ddfield(INDECIES_f,gNum tin)//evaluates the double time derivative of the field fld (s) at i,j,k.
{
      return galileon3(s,fld,i,j,k,tin);
}



/** RK2 function **/

void step()//this steps (integrates) the field and it time derivative via the rk2 method.
{
    
    DECLARE_INDEX
    
    gNum tin=t;
#if expansion_type==0
    //no expansion note that this only calculates the energies at the end of the full step
    
    for(fld=0;fld<nflds;fld++)//the first part of the RK2 step
    {
        //paralleleizes over the index i
#pragma omp parallel for private (j,k) num_threads (tot_num_thrds)
        LOOP//loops over fld i,j,k
        {
            field[INDEX(1,fld,i,j,k)]=field[INDEX(0,fld,i,j,k)]+.5*dt*dfield[INDEX(0,fld,i,j,k)];
            dfield[INDEX(1,fld,i,j,k)]=dfield[INDEX(0,fld,i,j,k)]+.5*dt*ddfield(0,fld,i,j,k,tin);
        }
    }
    
    
    for (fld=0;fld<nflds;fld++)
    {
        
        
        for(k=0;k<N;k++)
        {
            for(i=0;i<N;i++)
            {
                dfield[INDEX(1,fld,N-1,i,k)]=-dfdr_pert(field,0,fld,N-1,i,k,tin);
                dfield[INDEX(1,fld,i,k,N-1)]=-dfdr_pert(field,0,fld,i,k,N-1,tin);
                dfield[INDEX(1,fld,k,N-1,i)]=-dfdr_pert(field,0,fld,k,N-1,i,tin);
                
                dfield[INDEX(1,fld,0,i,k)]=-dfdr_pert(field,0,fld,0,i,k,tin);
                dfield[INDEX(1,fld,i,k,0)]=-dfdr_pert(field,0,fld,i,k,0,tin);
                dfield[INDEX(1,fld,k,0,i)]=-dfdr_pert(field,0,fld,k,0,i,tin);
            }
        }
    }
    
    tin=t+.5*dt;
    for(fld=0;fld<nflds;fld++)//the second part of the RK2 step
    {
        //paralleleizes over the index i
#pragma omp parallel for private (j,k) num_threads (tot_num_thrds)
        LOOP//This returns the actual value of the field and derivative at t
        {
            field[INDEX(0,fld,i,j,k)]=field[INDEX(0,fld,i,j,k)]+dt*dfield[INDEX(1,fld,i,j,k)];
            dfield[INDEX(0,fld,i,j,k)]=dfield[INDEX(0,fld,i,j,k)]+dt*ddfield(1,fld,i,j,k,tin);
            
        }
    }
    
    for (fld=0;fld<nflds;fld++)
    {
        //tates way of doing BC   phidot=-dphidr where dphidr is dpidr-analytic dpidr
        for(k=0;k<N;k++)
        {
            for(i=0;i<N;i++)
            {
                dfield[INDEX(0,fld,N-1,i,k)]=-dfdr_pert(field,1,fld,N-1,i,k,tin);
                dfield[INDEX(0,fld,i,k,N-1)]=-dfdr_pert(field,1,fld,i,k,N-1,tin);
                dfield[INDEX(0,fld,k,N-1,i)]=-dfdr_pert(field,1,fld,k,N-1,i,tin);
                
                dfield[INDEX(0,fld,0,i,k)]=-dfdr_pert(field,1,fld,0,i,k,tin);
                dfield[INDEX(0,fld,i,k,0)]=-dfdr_pert(field,1,fld,i,k,0,tin);
                dfield[INDEX(0,fld,k,0,i)]=-dfdr_pert(field,1,fld,k,0,i,tin);
            }
        }
    }
    
    /***********NOT IMPLEMNTED FOR THIS VERSION***************/
    
    
#elif expansion_type==1
    
    for(fld=0;fld<nflds;fld++)//first step of the Rk2 integration
    {
#pragma omp parallel for private (j,k) num_threads (tot_num_thrds)
        LOOP//loops over fld i,j,k
        {
            field[INDEX(1,fld,i,j,k)]=field[INDEX(0,fld,i,j,k)]+.5*dt*dfield[INDEX(0,fld,i,j,k)];
            dfield[INDEX(1,fld,i,j,k)]=dfield[INDEX(0,fld,i,j,k)]+.5*dt*ddfield(0,fld,i,j,k);
        }
    }
    
    
    a[1]=a[0]+.5*dt*adot[0];//this does the first step of the RK2 for the scale factor
    calcEnergy(1);//this calculates the energy based on this half step
    adot[1]=adf(1);//this updates adot based off of the energy at this step
    
    for(fld=0;fld<nflds;fld++)//second step of the Rk2 integration
    {
#pragma omp parallel for private (j,k) num_threads (tot_num_thrds)
        LOOP//This returns the actual value of the field and derivative at t
        {
            field[INDEX(0,fld,i,j,k)]=field[INDEX(0,fld,i,j,k)]+dt*dfield[INDEX(1,fld,i,j,k)];
            dfield[INDEX(0,fld,i,j,k)]=dfield[INDEX(0,fld,i,j,k)]+dt*ddfield(1,fld,i,j,k);
        }
    }
    
    a[0]=a[0]+dt*adot[1];//this calculates the full step scale factor
    calcEnergy(0);//calculates the energy at the full step
    adot[0]=adf(0);//then calculates adot based off of the full step
    
    
    
    
#elif expansion_type==2
    
    for(fld=0;fld<nflds;fld++)//first step of the Rk2 integration
    {
#pragma omp parallel for private (j,k) num_threads (tot_num_thrds)
        LOOP//loops over fld i,j,k
        {
            field[INDEX(1,fld,i,j,k)]=field[INDEX(0,fld,i,j,k)]+.5*dt*dfield[INDEX(0,fld,i,j,k)];
            dfield[INDEX(1,fld,i,j,k)]=dfield[INDEX(0,fld,i,j,k)]+.5*dt*ddfield(0,fld,i,j,k);
        }
    }
    
    /* this may need to chage based off of user defined expansion*/
    a[1]=a[0]+.5*dt*adot[0];//this does the first step of the RK2 for the scale factor
    calcEnergy(1);//this calculates the energy based on this half step
    adot[1]=adf(1);//this updates adot based off of the energy at this step
    
    for(fld=0;fld<nflds;fld++)//second step of the Rk2 integration
    {
#pragma omp parallel for private (j,k) num_threads (tot_num_thrds)
        LOOP//This returns the actual value of the field and derivative at t
        {
            field[INDEX(0,fld,i,j,k)]=field[INDEX(0,fld,i,j,k)]+dt*dfield[INDEX(1,fld,i,j,k)];
            dfield[INDEX(0,fld,i,j,k)]=dfield[INDEX(0,fld,i,j,k)]+dt*ddfield(1,fld,i,j,k);
        }
    }
    /* this may need to change based off of user defined expansion*/
    a[0]=a[0]+dt*adot[1];//this calculates the full step scale factor
    calcEnergy(0);//calculates the energy at the full step
    adot[0]=adf(0);//then calculates adot based off of the full step
#endif
  
    
    /**************END BROKEN CODE**************/
    
}


void steprk4()//this steps (integrates) the field and it time derivative via the rk2 method.
{
    
    DECLARE_INDEX
    gNum tin=t;
    /** k1 step **/
    fld=0;
    //paralleleizes over the index i
#pragma omp parallel for private (j,k) num_threads (tot_num_thrds)
    LOOP//loops over fld i,j,k
    {
        field[INDEX(1,fld,i,j,k)]=field[INDEX(0,fld,i,j,k)]+.5*dt*dfield[INDEX(0,fld,i,j,k)];
        dfield[INDEX(1,fld,i,j,k)]=dfield[INDEX(0,fld,i,j,k)]+.5*dt*ddfield(0,fld,i,j,k,tin);
    }
    
    for(k=0;k<N;k++)
    {
        for(i=0;i<N;i++)
        {
            dfield[INDEX(1,fld,N-1,i,k)]=-dfdr_pert(field,1,fld,N-1,i,k,tin);
            dfield[INDEX(1,fld,i,k,N-1)]=-dfdr_pert(field,1,fld,i,k,N-1,tin);
            dfield[INDEX(1,fld,k,N-1,i)]=-dfdr_pert(field,1,fld,k,N-1,i,tin);
            
            dfield[INDEX(1,fld,0,i,k)]=-dfdr_pert(field,1,fld,0,i,k,tin);
            dfield[INDEX(1,fld,i,k,0)]=-dfdr_pert(field,1,fld,i,k,0,tin);
            dfield[INDEX(1,fld,k,0,i)]=-dfdr_pert(field,1,fld,k,0,i,tin);
        }
    }
    
    /** k2 step **/
    tin=t+.5*dt;
#pragma omp parallel for private (j,k) num_threads (tot_num_thrds)
    LOOP
    {
        field[INDEX(2,fld,i,j,k)]=field[INDEX(0,fld,i,j,k)]+.5*dt*dfield[INDEX(1,fld,i,j,k)];
        dfield[INDEX(2,fld,i,j,k)]=dfield[INDEX(0,fld,i,j,k)]+.5*dt*ddfield(1,fld,i,j,k,tin);
        
    }

    for(k=0;k<N;k++)
    {
        for(i=0;i<N;i++)
        {
            dfield[INDEX(2,fld,N-1,i,k)]=-dfdr_pert(field,2,fld,N-1,i,k,tin);
            dfield[INDEX(2,fld,i,k,N-1)]=-dfdr_pert(field,2,fld,i,k,N-1,tin);
            dfield[INDEX(2,fld,k,N-1,i)]=-dfdr_pert(field,2,fld,k,N-1,i,tin);
            
            dfield[INDEX(2,fld,0,i,k)]=-dfdr_pert(field,2,fld,0,i,k,tin);
            dfield[INDEX(2,fld,i,k,0)]=-dfdr_pert(field,2,fld,i,k,0,tin);
            dfield[INDEX(2,fld,k,0,i)]=-dfdr_pert(field,2,fld,k,0,i,tin);
        }
    }
    
    
    /** k3 step **/
    tin=t+.5*dt;
#pragma omp parallel for private (j,k) num_threads (tot_num_thrds)
    LOOP
    {
        field[INDEX(3,fld,i,j,k)]=field[INDEX(0,fld,i,j,k)]+dt*dfield[INDEX(2,fld,i,j,k)];
        dfield[INDEX(3,fld,i,j,k)]=dfield[INDEX(0,fld,i,j,k)]+dt*ddfield(2,fld,i,j,k,tin);
        
    }
    
    for(k=0;k<N;k++)
    {
        for(i=0;i<N;i++)
        {
            dfield[INDEX(3,fld,N-1,i,k)]=-dfdr_pert(field,3,fld,N-1,i,k,tin);
            dfield[INDEX(3,fld,i,k,N-1)]=-dfdr_pert(field,3,fld,i,k,N-1,tin);
            dfield[INDEX(3,fld,k,N-1,i)]=-dfdr_pert(field,3,fld,k,N-1,i,tin);
            
            dfield[INDEX(3,fld,0,i,k)]=-dfdr_pert(field,3,fld,0,i,k,tin);
            dfield[INDEX(3,fld,i,k,0)]=-dfdr_pert(field,3,fld,i,k,0,tin);
            dfield[INDEX(3,fld,k,0,i)]=-dfdr_pert(field,3,fld,k,0,i,tin);
        }
    }
    
    
    
    /** k4 and actual value **/
    tin=t+dt;
#pragma omp parallel for private (j,k) num_threads (tot_num_thrds)
    LOOP
    {
        field[INDEX(0,fld,i,j,k)]=(-field[INDEX(0,fld,i,j,k)]+field[INDEX(1,fld,i,j,k)]+2.*field[INDEX(2,fld,i,j,k)]+field[INDEX(3,fld,i,j,k)]+.5*dt*dfield[INDEX(3,fld,i,j,k)])/3.;
        dfield[INDEX(0,fld,i,j,k)]=(-dfield[INDEX(0,fld,i,j,k)]+dfield[INDEX(1,fld,i,j,k)]+2.*dfield[INDEX(2,fld,i,j,k)]+dfield[INDEX(3,fld,i,j,k)]+.5*dt*ddfield(3,fld,i,j,k,tin))/3.;
    }
    
    
    for(k=0;k<N;k++)
    {
        for(i=0;i<N;i++)
        {
            dfield[INDEX(0,fld,N-1,i,k)]=-dfdr_pert(field,0,fld,N-1,i,k,tin);
            dfield[INDEX(0,fld,i,k,N-1)]=-dfdr_pert(field,0,fld,i,k,N-1,tin);
            dfield[INDEX(0,fld,k,N-1,i)]=-dfdr_pert(field,0,fld,k,N-1,i,tin);
            
            dfield[INDEX(0,fld,0,i,k)]=-dfdr_pert(field,0,fld,0,i,k,tin);
            dfield[INDEX(0,fld,i,k,0)]=-dfdr_pert(field,0,fld,i,k,0,tin);
            dfield[INDEX(0,fld,k,0,i)]=-dfdr_pert(field,0,fld,k,0,i,tin);
        }
    }
    
}


