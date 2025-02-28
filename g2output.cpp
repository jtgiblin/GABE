/****************************
 OUTPUT FILE
 ****************************/

/*
 This header file contains all the functions for the output of the data. Note that the outputslice() function needs to be generalized for n-fields and right now only gives meaning full output for 1 field. The only information outputted is the value of the field for the slice, the times at which the slices were output and the info.dat file.
 
 Copyright (2013):
 Kenyon College
 John T. Giblin, Jr
 Tate Deskins and Hillary Child
 Last Updated: 06.27.2013
 */

#include "g2header.h" //contains declarations for program functions.

//readable time variables

time_t tStart,tCurrent,tLast,tFinish; // Keep track of elapsed clock time

/*****************
 output stuffs
 *****************/

void outputfield(int first)//outputs the field values
{
    static FILE *slicefield;
    //static int fld,i,j,k;
    static char name[30];
    
    sprintf(name,"./slices/slices_fields_%d.dat", first);
    //sprintf(name,"slices_field%d_%d.dat", fld, first);
    slicefield=fopen(name,"w");
    
#if field_outdim>2//outputs slice for 3dimensions
    for(int i=0;i<N;i+=field_sliceskip)
#else
    int i=0;
#endif
    {
#if field_outdim>1//outputs slice for 3dimensions
        for(int j=0;j<N;j+=field_sliceskip)
#else
        int j=0;
#endif
        {
            for(int k=0;k<N;k+=field_sliceskip){
                for(int fld=0;fld<nflds;fld++){
#if fftw_flag==1
                    fprintf(slicefield,"%Le ", field[0][fld][i][j][k]);//,rescale_B*dfield[0][fld][i][j][k]);
#else
                    fprintf(slicefield,"%e ", field[0][fld][i][j][k]);//,rescale_B*dfield[0][fld][i][j][k]);
#endif
                }
#if fftw_flag==1
                fprintf(slicefield,"%Le ", calcrho(0,i,j,k));//,rescale_B*dfield[0][fld][i][j][k]);
#else
                fprintf(slicefield,"%e ", calcrho(0,i,j,k));//,rescale_B*dfield[0][fld][i][j][k]);
#endif
                
                fprintf(slicefield,"\n");
            }
        }
    }


    fclose(slicefield);
}

void meansvars()//calculates the mean and variance of each field
{
    char name_[500];
    static FILE *meansvarsout_;
    int i, j, k, fld;
    gNum av,av_sq,var;
    
    sprintf(name_,"./slices/meansvariances.dat");
    meansvarsout_=fopen(name_,"a");
    
#if fftw_flag==1
    fprintf(meansvarsout_,"%Lf",t);
#else
    fprintf(meansvarsout_,"%f",t);
#endif
    
    for(fld=0;fld<nflds;fld++)
    {
        av=0.;
        var=0.;
        // Calculate field mean
        LOOP
        av += field[0][fld][i][j][k];
        av = av/(gNum)gridsize; // Convert sum to average
        av_sq = av*av; // Calculate mean squared for finding variance
        // Calculate variance
        LOOP
        var += field[0][fld][i][j][k]*field[0][fld][i][j][k] - av_sq;
        var = var/(gNum)gridsize; // Convert sum to variance
        // Output means and variances to files
#if fftw_flag==1
        fprintf(meansvarsout_," %Le %Le",av,var);
#else
        fprintf(meansvarsout_," %e %e",av,var);
#endif
        // Check for instability. See if the field has grown exponentially and become non-numerical at any point.
        if((av!=0. && av/av!=1.)) // These two separate checks between them work on all the compilers I've tested
        {
#if fftw_flag==1
            printf("Unstable solution developed. Field %d not numerical at t=%Le\n",fld,t);
#else
            printf("Unstable solution developed. Field %d not numerical at t=%e\n",fld,t);
#endif
            
            output_parameters();
            fflush(meansvarsout_);
            exit(1);
        }
    } // End of loop over fields
    fprintf(meansvarsout_,"\n");
    fclose(meansvarsout_);
}


int slicewaitf()
{
    return (int) (endtime/dt/slicenumber);//calculates number of timesteps to wait for slicenumber slices by time endtime
}

void outputslice()//externally called function for outputing the data from the run
{
    static FILE  *slicetime;
    static int first=0;
    static int lastslice,lastspec;
    static int slicewaitm;
    
    if(first==0)
    {
        if(slicewait==0)//determins number of steps to wait between each slice output based off of g2parameters.h
            slicewaitm=slicewaitf();
        else
            slicewaitm=slicewait;
        
        lastslice=slicewaitm;
        lastspec=specnumber;
    }
    
    lastslice++;
    
    if(lastslice>=slicewaitm)//this statement delays the output every slicewait*dt
    {
        
        
        
#if field_outdim!=0
        outputfield(first);//outputs field slice
#endif
        
#if spec_output==1//outputs spectra if last spectra output was not too soon (see g2paramters.h)
        lastspec++;
        if (lastspec>=specnumber) {
            specOut(first/specnumber);
            lastspec=0;
        }
#endif
        
#if var_output!=0
        meansvars();//outputs mean and variance of fields
#endif
        
        
        
        /*times file output */
        //this routine outputs time, scale factor, scale factor derivative, hubble constant, and energy components at each slice output
        slicetime=fopen("./slices/slices_time.dat","a");
        //slicetime=fopen("slices_time.dat","a");
#if fftw_flag==1
        fprintf(slicetime,"%Le %Le %Le %Le %Le %Le %Le %Le\n", t,a[0],adot[0],adot[0]/a[0],edkin[0],edpot[0],edgrad[0],edrho[0]);
#else
        fprintf(slicetime,"%e %e %e %e %e %e %e %e\n", t,a[0],adot[0],adot[0]/a[0],edkin[0],edpot[0],edgrad[0],edrho[0]);
#endif
        fclose(slicetime);
        
        first++;
        lastslice=0;
    }
}



void output_parameters()//this creates info.dat which contains information about run parameters and statistics
{
    static FILE *info;
    
    static int first=1;
    if(first) // At beginning of run output run parameters
    {
        char expanType[13];
#if expansion_type==1
        sprintf(expanType,"a dot");
#elif expansion_type==2
        sprintf(expanType,"a duble dot");
#else
        sprintf(expanType,"no");
#endif
        
        info=fopen("info.txt","a");
        
        fprintf(info,"--------------------------\n");
        fprintf(info,"Model Specific Information\n");
        fprintf(info,"--------------------------\n");
        modelinfo(info);
        
        fprintf(info,"\n--------------------------\n");
        fprintf(info,"General Program Information\n");
        fprintf(info,"-----------------------------\n");
        fprintf(info,"Grid size=%d^3\n",N);
#if fftw_flag==1
        fprintf(info,"L=%Le\n",L);
        fprintf(info,"dt=%Le, dx=%Le\n",dt,dx);
        fprintf(info,"end time=%Le\n",endtime);
        fprintf(info,"B=%Le\n",rescale_B);
#else
        fprintf(info,"L=%e\n",L);
        fprintf(info,"dt=%e, dx=%e\n",dt,dx);
        fprintf(info,"end time=%e\n",endtime);
        fprintf(info,"B=%e\n",rescale_B);
#endif
        fprintf(info, "\nUsing %s expansion\n",expanType);
        fprintf(info, "\nUsing %d cores\n",tot_num_thrds);
        fprintf(info, "%d momentum bins for spectra\n",((int)(sqrt(3.)*N/2+1)));
        time(&tStart);
        fprintf(info,"\nRun began at %s",ctime(&tStart)); // Output date in readable form
        first=0;
    }
    else // If not at beginning, record elapsed time for run
    {
        time(&tFinish);
        fprintf(info,"Run ended at %s",ctime(&tFinish)); // Output ending date
#if fftw_flag==1
        fprintf(info,"\nRun from t=%Le to t=%Le took ",starttime,t);
        readable_time((int)(tFinish-tStart),info);
        fprintf(info,"\n");
        fprintf(info, "\nFinal scale factor is %Le\n",a[0]);
#else
        fprintf(info,"\nRun from t=%e to t=%e took ",starttime,t);
        readable_time((int)(tFinish-tStart),info);
        fprintf(info,"\n");
        fprintf(info, "\nFinal scale factor is %e\n",a[0]);
        
#endif
        
    }
    
    fflush(info);
    return;
}

// Convert a time in seconds to a more readable form and print the results
void readable_time(int tt, FILE *info)
{
    int tminutes=60,thours=60*tminutes,tdays=24*thours;
    
    if(tt==0)
    {
        fprintf(info,"less than 1 second\n");
        return;
    }
    
    // Days
    if(tt>tdays)
    {
        fprintf(info,"%d days",tt/tdays);
        tt = tt%tdays;
        if(tt>0)
            fprintf(info,", ");
    }
    // Hours
    if(tt>thours)
    {
        fprintf(info,"%d hours",tt/thours);
        tt = tt%thours;
        if(tt>0)
            fprintf(info,", ");
    }
    // Minutes
    if(tt>tminutes)
    {
        fprintf(info,"%d minutes",tt/tminutes);
        tt = tt%tminutes;
        if(tt>0)
            fprintf(info,", ");
    }
    // Seconds
    if(tt>0)
        fprintf(info,"%d seconds",tt);
    fprintf(info,"\n");
    return;
}



void screenout()//this calculates the time elapsed from last screen output before outputting current program time.
{
    
    time(&tCurrent);
    
    if(screentime==0)
        return;
    
    if(tCurrent-tLast>screentime)
    {
#if fftw_flag==1
        printf("%Lf\n",t);
#else
        printf("%f\n",t);
#endif
        time(&tLast);
    }
}
