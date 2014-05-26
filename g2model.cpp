/**************************************
 MODEL FILE
 **************************************/


/*
 This header file contains contains all the model dependet functions necessary for the program. 
 Note that model dependent parameters should go in g2parameters.h. 
 Note that all functions must be in program units which are given by the following rescallings (pr denotes quantity used in program). 
 Any other model functions should be added here (and may need to tweek g2function.cpp and g2output.cpp). 
 
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
    fprintf(info,"Basic 3-d wave equation solver.\n");

}


/***************************
 user defined model functions
 ***************************/
#define PHI field[s][0]
#define PHIDOT dfield[s][0]


/********************
 field initialization
 ********************/

void initfields()//here the user may decide how the fields will be initialized
{
    static int first=0,s=0;
    DECLARE_INDEX
    
    if(first==0)
    {
        fldLOOP//loops over fld i,j,k
        {
			field[s][fld][i][j][k]=f0[fld];//initialize each field as its initial value
            dfield[s][fld][i][j][k]=df0[fld];// initialize each field derivative as its initial value
        }
   	
    }
    	
    calcEnergy(0); //This is important -- needed for first step of evolution

}
#undef PHI
#undef PHIDOT
