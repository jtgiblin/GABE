/**************************************
 MODEL FILE
 **************************************/


/*
 This header file contains contains all the model dependet functions necessary for the program.
 Note that model dependent parameters should go in g2parameters.h.
 Note that all functions must be in program units which are given by the following rescallings (pr denotes quantity used in program).
 Any other model functions should be added here (and may need to tweek g2function.cpp and g2output.cpp).
 
 
 B=mphi //Note that B may change depending on your model
 dt_(pr)=dt*B
 x_(pr)=x*B
 f_(pr)=f
 V_(pr)=1/B^2*V
 dV_(pr)/df_(pr)=1/B^2*dV/df
 
 
 Copyright (2013):
 Kenyon College
 John T. Giblin, Jr
 Tate Deskins and Hillary Child
 Last Updated: 06.27.2013
 */

#include "g2header.h" //contains declerations for program functions.

// Model specific details about the run to be output to an information file
//change as needed
void modelinfo(FILE *info)
{
    // Name and description of model
    fprintf(info,"Testing GR and Galileons from Binary System.\n");
    fprintf(info,"Using Hulse-Taylor Binary Pulsar as Base Model\n Model Parametrs\n");
    
    // Model specific parameter values
    fprintf(info,"sigma2 = %Le",sigma2);
    fprintf(info,"alpha = %Le",alpha);
    fprintf(info,"Omega = %Le\n",omega);
    fprintf(info,"Gallileon Order = %d", gallileon_order);
    fprintf(info,"Kappa = %Le\n",kappa);
    fprintf(info,"Kappa^2 = %Le\n",kappa2);
    
    
    
#if rand_init==1
    const char initType[6]="U";
#else
    const char initType[6]="Not u";
#endif
    fprintf(info,"%ssing random initial conditions",initType);
}


/***************************
 user defined model functions
 ***************************/
#define PI field,s,nflds-1,i,j,k
#define dPI dfield,s,nflds-1,i,j,k




gNum potential(INDECIES_f)//user defined potential
{
    return 0;
}

gNum profile(gNum r)
{
    return 0.49330714907571516+0.5*tanh((.1*r -2.5));
}




gNum massGaussian1 (INDECIES_f,gNum tin)
{
    return profile(tin)*(PNorm*exp(-1.*(pw2(dx*(i-N/2.)-cos(omega*tin)) + pw2(dx*(j-N/2.)-sin(omega*tin)) + pw2(dx*(k-N/2.)))/(sigma2)));
    
    // the last term is so that the avg value of the source is zero. note that there is something odd with normalization going on. or not....
}

gNum massGaussian2 (INDECIES_f,gNum tin)
{
    return profile(tin)*alpha*PNorm*exp(-1.*(pw2(dx*(i-N/2.)+cos(omega*tin)) + pw2(dx*(j-N/2.)+sin(omega*tin)) + pw2(dx*(k-N/2.)))/(sigma2));
}


gNum energyDensity(INDECIES_f,gNum tin)
{
    return - massGaussian1(INDECIES,tin)-massGaussian2(INDECIES,tin);
}


gNum galileon2(INDECIES_f,gNum tin)
{
    
    
    return (energyDensity(INDECIES,tin)*lambda
            +9.*(
              dfdii(PI)
              +dfdjj(PI)
              +dfdkk(PI)
              )
#if gallileon_order>=3
            +2.*(
                 dfdi(dPI)*dfdi(dPI)
                 +dfdj(dPI)*dfdj(dPI)
                 +dfdk(dPI)*dfdk(dPI)
                 -dfdij(PI)*dfdij(PI)
                 -dfdjk(PI)*dfdjk(PI)
                 -dfdik(PI)*dfdik(PI)
                 +dfdii(PI)*dfdjj(PI)
                 +dfdii(PI)*dfdkk(PI)
                 +dfdjj(PI)*dfdkk(PI)
                 )*kappa
#if gallileon_order>=4 //need to check that it agrees
            +6.*(
                 dfdi(dPI)*dfdi(dPI)*dfdkk(PI)
                 +dfdi(dPI)*dfdi(dPI)*dfdjj(PI)
                 +dfdj(dPI)*dfdj(dPI)*dfdii(PI)
                 +dfdj(dPI)*dfdj(dPI)*dfdkk(PI)
                 +dfdk(dPI)*dfdk(dPI)*dfdii(PI)
                 +dfdk(dPI)*dfdk(dPI)*dfdjj(PI)
                 -2.*dfdi(dPI)*dfdj(dPI)*dfdij(PI)
                 -2.*dfdi(dPI)*dfdk(dPI)*dfdik(PI)
                 -2.*dfdj(dPI)*dfdj(dPI)*dfdjk(PI)
                 +2.*dfdij(PI)*dfdik(PI)*dfdjk(PI)
                 -dfdii(PI)*dfdjk(PI)*dfdjk(PI)
                 -dfdjj(PI)*dfdik(PI)*dfdik(PI)
                 -dfdkk(PI)*dfdij(PI)*dfdij(PI)
                 +dfdii(PI)*dfdjj(PI)*dfdkk(PI)
                 )*kappa2
#endif
#endif
            )/(9.
#if gallileon_order>=3
    
        +2.*(
              dfdii(PI)
              +dfdjj(PI)
              +dfdkk(PI)
              )*kappa
#if gallileon_order>=4 //need to check that it agrees
       +6.*(
           dfdij(PI)*dfdij(PI)
           +dfdik(PI)*dfdik(PI)
           +dfdjk(PI)*dfdjk(PI)
           -dfdii(PI)*dfdjj(PI)
           -dfdii(PI)*dfdkk(PI)
           -dfdjj(PI)*dfdkk(PI)
           )*kappa2
#endif
#endif
       )
    ;
    //the new comment
}


inline gNum effMass(INDECIES_f)//the effective mass used for random inital conditions
{
    return 0.;
}
/*******************
 field initialization
 *******************/

void initfields()//here the user may decide how the fields will be initialized
{
    static int first=0,s=0;
    DECLARE_INDEX
    
    if(first==0)
    {
        fldLOOP//loops over fld i,j,k
        {
            field[INDEX(s,fld,i,j,k)]=f0[fld];//initialize each field as its initial value
            dfield[INDEX(s,fld,i,j,k)]=df0[fld];// initialize each field derivative as its initial value
        }
        
    }
    
    if(first==1)
    {
#if rand_init==1
        for(fld=0; fld<nflds; fld++){
            randInit(field[s][fld],dfield[s][fld],effMass(s,fld));//adds random intial conditions ontop of mean value above
        }
        initDestroy();
        printf("Fields fluctuated\n");
#endif
        
        //Any other model specific initialization can go here -- i.e. Bubbles, etc	
    }
    
    calcEnergy(0); //This is important -- needed for first step of evolution
    
    first++;
}
#undef PI
#undef dPI




