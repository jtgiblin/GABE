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
#define PHI field[s][0]
#define CHI field[s][1]
#define PHIDOT dfield[s][0]
#define CHIDOT dfield[s][1]

#define COUP gsq/mphi/mphi //coupling term

// user defined potential
real_t potential(int s, int i, int j, int k)
{
    return (0.5*PHI[i][j][k]*PHI[i][j][k]
        + 0.5*COUP*PHI[i][j][k]*PHI[i][j][k]*CHI[i][j][k]*CHI[i][j][k]);
}

// user defined derivative of the potential
real_t dVdf(int s, int fld, int i, int j, int k)

{
    switch (fld) {
        case 0:
            //derivative with respect to phi
            return (PHI[i][j][k] + COUP*CHI[i][j][k]*CHI[i][j][k]*PHI[i][j][k]);
        case 1:
            //derivative with respect to chi
            return (COUP*PHI[i][j][k]*PHI[i][j][k]*CHI[i][j][k]);
        default:
            return 0;
    }

}

// the effective mass used for random inital conditions
inline real_t effMass(int s, int fld)
{
    real_t avemass=0.;
    int i,j,k;
    switch (fld) {
        case 0:
            LOOP{
                avemass += (1. + COUP*CHI[i][j][k]*CHI[i][j][k]);
            }
            return avemass/gridsize;
        case 1:
            LOOP{
                avemass += (COUP*PHI[i][j][k]*PHI[i][j][k]);
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
    
    if(first==0)
    {
        //loops over fld i,j,k
        fldLOOP
        {
            field[s][fld][i][j][k]=f0[fld]; // initialize each field as its initial value
            dfield[s][fld][i][j][k]=df0[fld]; // initialize each field derivative as its initial value
        }
    }
    
    if(first==1)
    {
        // Any other model specific initialization can go here -- i.e. Bubbles, etc
    }
    
    calcEnergy(0); //This is important -- needed for first step of evolution

    first++;
}

#undef PHI
#undef CHI
#undef PHIDOT
#undef CHIDOT
