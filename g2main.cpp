/****************************
 LatticeSolve Main File
 ****************************/

/*
 This is the .cpp file that contains the main function for the GABE2.0 program. 
 In this function the looping over time occurs aswell as the calls to the output functions and the initialize function. 
 Also all the header calls and macro definitions are made in this file.
 
 Copyright (2013): 
 Kenyon College
 John T. Giblin, Jr
 Tate Deskins and Hillary Child
 Last Updated: 06.27.2013
*/

/*********************
 GABE Header
 *********************/
 //contains all the parameters, changeable and unchangeable
#include "g2header.h" //contains declarations for program functions and parameters file


gNum t=starttime;//this is the variable that stores the time
gNum (*field)[nflds][N][N][N];//this stores the field values for each step along the grid
gNum (*dfield)[nflds][N][N][N];//this stores the derivative of the field for each step along the grid
gNum a[2];//this stores the scale factor for each step
gNum adot[2];// this stores the time derivative of the scale factor
gNum edpot[2]; //this stores the average potential energy
gNum edkin[2]; //this stores the average kinetic energy
gNum edgrad[2]; //this stores the average gradient energy
gNum edrho[2]; // this stores the avg. energy density over the box

#if calc_gws ==1 //if we have gravitational waves
gNum (*h)[12][N][N][N/2+1];
gNum (*l)[12][N][N][N/2+1];
gNum hub = 0;
gNum (*T_gw)[12][N][N][N/2+1];
fftw_plan plan_fft_gw;
#endif

int main()
{
    omp_set_nested(1);//allows for nested parallelization
    printf("\n\nLattice evolution program started\n\n");
    output_parameters(); //Outputs general run information (info.txt)
    printf("Info file made\n");
    
    field=(gNum(*)[nflds][N][N][N]) malloc(sizeof(gNum)*nflds*2*N*N*N);//allocates memory for the fields
    printf("'field' memory allocated\n");
    
#if calc_gws ==1
    h = (gNum (*)[12][N][N][N/2+1]) malloc(sizeof(gNum)*2.*(N/2+1)*N*N*12);
    l = (gNum (*)[12][N][N][N/2+1]) malloc(sizeof(gNum)*2.*(N/2+1)*N*N*12);
    T_gw = (gNum (*)[12][N][N][N/2+1]) malloc(sizeof(gNum)*2.*(N/2+1)*N*N*12);
    
    //zero out metric perturbations, just in case
    for(int temps=0; temps<2; temps++){
        int j;
        int k;
#pragma omp parallel for private (j,k) num_threads (tot_num_thrds)
        for(int i=0;i<N;i++){
            for(j=0;j<N;j++){
                for(k=0;k<N/2+1;k++){
                    
                    h[temps][0][i][j][k]=0.;
                    h[temps][1][i][j][k]=0.;
                    h[temps][2][i][j][k]=0.;
                    h[temps][3][i][j][k]=0.;
                    h[temps][4][i][j][k]=0.;
                    h[temps][5][i][j][k]=0.;
                    h[temps][6][i][j][k]=0.;
                    h[temps][7][i][j][k]=0.;
                    h[temps][8][i][j][k]=0.;
                    h[temps][9][i][j][k]=0.;
                    h[temps][10][i][j][k]=0.;
                    h[temps][11][i][j][k]=0.;
                    
                    l[temps][0][i][j][k]=0.;
                    l[temps][1][i][j][k]=0.;
                    l[temps][2][i][j][k]=0.;
                    l[temps][3][i][j][k]=0.;
                    l[temps][4][i][j][k]=0.;
                    l[temps][5][i][j][k]=0.;
                    l[temps][6][i][j][k]=0.;
                    l[temps][7][i][j][k]=0.;
                    l[temps][8][i][j][k]=0.;
                    l[temps][9][i][j][k]=0.;
                    l[temps][10][i][j][k]=0.;
                    l[temps][11][i][j][k]=0.;
                    
                }
            }
        }
    }
    printf("Gravitational Wave calculation memory allocated\n");
#endif
    
    dfield=(gNum(*)[nflds][N][N][N]) malloc(sizeof(gNum)*nflds*2*N*N*N);//allocates memory for the fields' time derivatives
    printf("'dfield' memory allocated\n");
    
    printf("Memory allocated for all arrays\n\n");
    printf("Starting run...\n\n");
    
    
    for(t=starttime;t<=endtime;t+=dt)//This is the main for-loop over time. Note that the current t value is the time for the currently evaluated fields
    {
		
        if(t==starttime)
        {
            initfields();//initializes fields as defined in model.h
            printf("Fields initialized (no fluctuations)\n");
            initexpansion();//start expansion;
            printf("Expansion started\n");
			
            initfields();//do fluctuations
            printf("Model Specific Initialization Completed\n");

            printf("Time evolution begining\n");
            
			screenout();//outputs current pr time to screen (see g2output.cpp)
			outputslice();//outputs slice values to file

        }

		step();//evolves the fields one step
        
		screenout();//outputs current pr time to screen
        outputslice();//outputs slice values to file
        
	}
#if spec_output==1
    specClear();//clears fftw stuff needed for spectra at the end of the run
#endif
    printf("\nRun complete\n\n");
    output_parameters();//finalizes any necessary info for the run to info.txt)
    printf("GABE2.2 finished\n\n\n");
}
