/**************************************
 MODEL FILE
 **************************************/

/*
 This header file contains contains all the model dependent functions
 necessary for the program. 
 Note that model dependent parameters should go in g2parameters.h. 
 Note that all functions must be in program units which are given by the
 following rescallings (pr denotes quantity used in program). 
 Any other model functions should be added here (and may need to tweak
 g2function.cpp and g2output.cpp).
 
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
s*/

#include "g2header.h" //contains declerations for program functions.


/***************************
 user defined model functions
 ***************************/
#define PHI field[FIELD(s,0)]
#define CHI field[FIELD(s,1)]
#define PHIDOT dfield[FIELD(s,0)]
#define CHIDOT dfield[FIELD(s,1)]

#define COUP gsq/mphi/mphi //coupling term

// user defined potential
real_t potential(int s, int idx)
{
    return (0.5*PHI[idx]*PHI[idx]
        + 0.5*COUP*PHI[idx]*PHI[idx]*CHI[idx]*CHI[idx]);
}

// user defined derivative of the potential
real_t dVdf(int s, int fld, int idx)

{
    switch (fld) {
        case 0:
            //derivative with respect to phi
            return (PHI[idx] + COUP*CHI[idx]*CHI[idx]*PHI[idx]);
        case 1:
            //derivative with respect to chi
            return (COUP*PHI[idx]*PHI[idx]*CHI[idx]);
        default:
            return 0;
    }
}

// the effective mass used for random inital conditions
inline real_t effMass(int s, int fld)
{
    real_t avemass=0.;
    switch (fld) {
        case 0:
            for(int p=0; p<POINTS; p++)
            {
                avemass += (1. + COUP*CHI[p]*CHI[p]);
            }
            return avemass/gridsize;
        case 1:
            for(int p=0; p<POINTS; p++)
            {
                avemass += (COUP*PHI[p]*PHI[p]);
            }
            return avemass/gridsize;
        default:
            // set mass as 1 (rescaled) if there is no case structure
            return 1.;
    }

    return 0;
}

/*******************
 field initialization
 *******************/

// here the user may decide how the fields will be initialized
void initfields()
{
    static int first=0,s=0;
    DECLARE_INDEX

    //loops over fld i,j,k
    for(int fld=0; fld<nflds; fld++)
    LOOP
    {
        int idx = IDX(i,j,k);
        field[FIELD(s,fld)][idx] = f0[fld] + 0.01*sin(2.0*3.14159*j/NY); // initialize each field as its initial value
        dfield[FIELD(s,fld)][idx] = df0[fld]; // initialize each field derivative as its initial value
    }
        
    calcEnergy(0); //This is important -- needed for first step of evolution

    first++;
}

#undef PHI
#undef CHI
#undef PHIDOT
#undef CHIDOT
