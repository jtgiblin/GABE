/**********************************
 SPECTRA FILE
 **********************************/

/*
 This file contains the functions needed to output spectra.
 
 Copyright (2013): 
 Kenyon College
 John T. Giblin, Jr
 Tate Deskins and Hillary Child
 Last Updated: 06.27.2013
*/

#include "g2header.h" //contains declarations for program functions.

#if spec_output==1

static gNum *indft; //declares input field for dft
static fftw_complex *outdft; //declares output field for dft
static fftw_plan spec_plan; //declares dft plan for fftw

void specOut(int first)
{
    
    DECLARE_INDEX
    const gNum spec_norm=pow(L/(gNum)N,3.); //normalization factor for the spectra
    int px,py,pz;//tracks real place in momentum grid
    gNum pmagnitude;//stores the magnitude of the momentum vector for current point
    const int numbins=((int)(sqrt(3.)*N/2+1));//number of spectra bins based off of size of the grid
    gNum spec_power[nflds][numbins], spec_numpts[nflds][numbins];// holds the power and number of points in each spectra bin
       
    if(first==0)
    {
        
        fftw_plan_with_nthreads(tot_num_thrds);//lets fftw know the number of available threads (tot_num_thrds)
        indft = (gNum *) fftw_malloc(sizeof(gNum) * N * N * N);//allocates memory for the input array
        outdft =  (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N * (N/2+1));//allocates memory for the output array
        spec_plan = fftw_plan_dft_r2c_3d(N, N, N, indft, outdft, FFTW_MEASURE);// defines the plan for fftw (moving from configuration space to momentum space
        
    }
    
    
    for(fld=0;fld<nflds;fld++){
        LOOP//this loop asigns the field values to the input array
        {
                indft[k + N*(j + N*i)]=field[0][fld][i][j][k];
        }
        
        fftw_execute(spec_plan);//performs the dft
        
        for(i=0;i<numbins;i++){//initialize the power and number of points in each bin to 0
            spec_power[fld][i]=0.0;
            spec_numpts[fld][i]=0;
        }
        
        for(i=0;i<N;i++) {//this is the loop over momentum space
            px=(i<=N/2 ? i: i-N);//calculates the real momentum space position for the x direction
            for(j=0; j<N;j++) {
                py = (j<=N/2 ? j : j-N); //calculates the real momentum space position for the y direction
                for(k=1;k<N/2;k++) { //since the k=0 and k=n/2+1 are un matched we treat the separately
                    pz=k;//calculates the real momentum space position for the z direction
                    pmagnitude = sqrt(pw2(px)+pw2(py)+pw2(pz));//calculates the magnitude of the momentum of the point
                    spec_power[fld][(int) pmagnitude]+=2.*(pw2(outdft[k + (N/2+1)*(j + N*i)][0])+pw2(outdft[k + (N/2+1)*(j + N*i)][1])); //adds the power in this mode to the proper bin times two 
                    //since the symmetry in the dft we only need half of the z direction box (see fftw documentation)
                    spec_numpts[fld][(int) pmagnitude]+=2;//adds two points for the 
                }
                //this is the same procedure as seen above with out the multiplication by 2 since the zero mode and the k=N/2 modes are unmatched
                k=0;
                pz=k;
                pmagnitude = sqrt(pw2(px)+pw2(py)+pw2(pz));
                spec_power[fld][(int) pmagnitude]+=(pw2(outdft[k + (N/2+1)*(j + N*i)][0])+pw2(outdft[k + (N/2+1)*(j + N*i)][1]));
                spec_numpts[fld][(int) pmagnitude]+=1;
                k=N/2;
                pz=k;
                pmagnitude = sqrt(pw2(px)+pw2(py)+pw2(pz));
                spec_power[fld][(int) pmagnitude]+=(pw2(outdft[k + (N/2+1)*(j + N*i)][0])+pw2(outdft[k + (N/2+1)*(j + N*i)][1]));
                spec_numpts[fld][(int) pmagnitude]+=1;
                
            }
        }
        
        for(i=0;i<numbins;i++){ //convert the sums in each bin to averages
            if(spec_numpts[fld][i]>0)
                spec_power[fld][i]=spec_power[fld][i]*spec_norm/((gNum) spec_numpts[fld][i]);
            else
                spec_power[fld][i]=0.;
        }
		
	}
    
    
    //Print the spectra to a file
    
    static FILE *slicespectra;
    static char name[30];
    
    sprintf(name,"./slices/slices_spectra_%d.dat", first);
    slicespectra=fopen(name,"w");
    
    for(i=0;i<numbins;i++)
    {
#if fftw_flag==1
        fprintf(slicespectra,"%Le", 2.*M_PI*i/L);//this prints the mode
        for(fld=0;fld<nflds;fld++){
            fprintf(slicespectra," %Le", spec_power[fld][i]);//and the power associated with it for each field
        }
#else
        fprintf(slicespectra,"%e", 2.*M_PI*i/L);//this prints the mode
        for(fld=0;fld<nflds;fld++){
            fprintf(slicespectra," %e", spec_power[fld][i]);//and the power associated with it for each field
        }

#endif
        fprintf(slicespectra,"\n");
    }

    fclose(slicespectra);
    
}

void specClear(){//clears the memory used in the dft's
    fftw_destroy_plan(spec_plan);
    fftw_free(indft);
    fftw_free(outdft);
}

#endif
