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

#include "g2header.h" //contains decelerations for program functions.
gNum gravrho1=0.;  // this stores the energy of the gravitational waves
gNum gravrho2=0.;

//readable time variables

time_t tStart,tCurrent,tLast,tFinish; // Keep track of elapsed clock time




/*****************
 output stuffs
 *****************/

gNum piPowerOut(gNum *bkgf){
    *bkgf=0;
    gNum mbkg0=0, mbkgN=0, flux=0;
    for (gIdx i=0; i<N; i++) {
        for (gIdx j=0; j<N; j++) {
            
            mbkg0=profile(t)*(-N/2.)/sqrtl(pw2(i-N/2.)+pw2(j-N/2.)+pw2(-N/2.))*dfdr_analytic(0,0,i,j,0,t);
            
            //1./pw2(4*M_PI*dx*dx*(pw2(i-N/2.)+pw2(j-N/2.)+pw2(-N/2.)));
            
            mbkgN=profile(t)*(N/2.-1.)/sqrtl(pw2(i-N/2.)+pw2(j-N/2.)+pw2(N/2.-1.))*dfdr_analytic(0,0,i,j,N-1,t);
            //1./pw2(4*M_PI*dx*dx*(pw2(i-N/2.)+pw2(j-N/2.)+pw2(N/2.-1.)));
            
            *bkgf+=(mbkgN)*dfield[INDEX(0, nflds-1,N-1,i,j)];
            
            *bkgf+=(mbkg0)*dfield[INDEX(0, nflds-1,0,i,j)];
            
            *bkgf+=(mbkgN)*dfield[INDEX(0, nflds-1,i,N-1,j)];
            
            *bkgf+=(mbkg0)*dfield[INDEX(0, nflds-1,i,0,j)];
            
            *bkgf+=(mbkgN)*dfield[INDEX(0, nflds-1,i,j,N-1)];
            
            *bkgf+=(mbkg0)*dfield[INDEX(0, nflds-1,i,j,0)];
            
            
            flux+=(dfdi_s(field,0,nflds-1,N-1,i,j)-mbkgN)*dfield[INDEX(0, nflds-1,N-1,i,j)];
            
            flux+=(dfdi_s(field,0,nflds-1,0,i,j)-mbkg0)*dfield[INDEX(0, nflds-1,0,i,j)];
            
            flux+=(dfdj_s(field,0,nflds-1,i,N-1,j)-mbkgN)*dfield[INDEX(0, nflds-1,i,N-1,j)];
            
            flux+=(dfdj_s(field,0,nflds-1,i,0,j)-mbkg0)*dfield[INDEX(0, nflds-1,i,0,j)];
            
            flux+=(dfdk_s(field,0,nflds-1,i,j,N-1)-mbkgN)*dfield[INDEX(0, nflds-1,i,j,N-1)];
            
            flux+=(dfdk_s(field,0,nflds-1,i,j,0)-mbkg0)*dfield[INDEX(0, nflds-1,i,j,0)];
        }
    }

    return flux;
}

   
void outputfield(int first)//outputs the field values
{
    static FILE *slicefield;
    static FILE *energyslice;
    static FILE *gravrhoslice;
    static FILE *sourceslice;
    
    static gIdx fld,i,j,k;
    
    static char name[5000];
    static char energyname[5000];
    static char gravname[5000];
    static char sourcename[5000];
    
    
    sprintf(name,"./slices/slices_fields_%d.dat", first);
    //sprintf(name,"slices_field%d_%d.dat", fld, first);
    slicefield=fopen(name,"w");
    
    sprintf(energyname, "./energy/energydensity_%d.dat", first);
    energyslice = fopen(energyname, "w");
    
    sprintf(gravname, "./gravrhoslices/slices_gravrho_%d.dat", first);
    //gravrhoslice = fopen(gravname, "w");
    
    sprintf(sourcename, "./source/slices_source_%d.dat", first);
    sourceslice = fopen(sourcename, "w");
    
#if field_outdim==3//outputs slice for 3 dimensions
    fld=0;
    for(i=1;i<N-1;i+=field_sliceskip){
        for(j=1;j<N-1;j+=field_sliceskip){
            for(k=1;k<N-1;k+=field_sliceskip){
                fprintf(slicefield,"%Le ", field[INDEX(0,fld,i,j,k)]);
                fprintf(sourceslice,"%Le ", dfield[INDEX(0, fld, i, j, k)]);
                //fprintf(sourceslice,"%Le ", energyDensity(0,fld,i,j,k,t));
                fprintf(energyslice,"%Le ", energyDensity(0,fld,i,j,k,t));
            }
            fprintf(slicefield,"\n");
            fprintf(sourceslice,"\n");
            fprintf(energyslice,"\n");
        }
        fprintf(slicefield,"\n");
        fprintf(sourceslice,"\n");
        fprintf(energyslice,"\n");
    }

//THE BELOW NEEDS TO BE UPDATED WITH THE CURRENT OUTPUT AVAILABLE THINGS
#elif field_outdim==2//outputs slice for 2dimensions
    
#if slice_orient==0  //for xy-slice
    for(int i=0;i<N;i+=field_sliceskip)
    {
        for(int j=0;j<N;j+=field_sliceskip)
        {
            fprintf(slicefield,"%Le %Le ",(i-N/2.)*dx, (j-N/2.)*dx);
            //fprintf(sourceslice,"%d %d ",i, j);
            //fprintf(energyslice,"%d %d %d %15.15Lf \n", i, j, N/2, energyDen(0,fld,i,j,N/2)+energyDen2(0,fld,i,j,N/2));

            for(int fld=0;fld<nflds;fld++)
            {
                fprintf(slicefield,"%Le %Le ", field[INDEX(0,fld,i,j,N/2)], pw2(dfield[INDEX(0,fld,i,j,N/2)])+gradF2(field,0,fld,i,j,N/2));
                // fprintf(sourceslice,"%Le ", source(0,fld,i,j,N/2)+source2(0,fld,i,j,N/2));
                
            }
            fprintf(slicefield, "\n");
            // fprintf(sourceslice, "\n");
            
        }
        //fprintf(slicefield, "\n");
    }
    
#elif slice_orient==1//for yz-slice
    
    for(int i=0;i<N;i+=field_sliceskip)
    {
        for(int j=0;j<N;j+=field_sliceskip)
        {
            gravenergy1(N/2,i,j);
            gravenergy2(N/2,i,j);
            fprintf(gravrhoslice,"%d %d %d %15.15Lf %15.15Lf \n",N/2, i, j, gravrho1,gravrho2);
            fprintf(slicefield,"%d %d %d ",N/2, i, j);
            fprintf(sourceslice,"%d %d %d ",N/2, i, j);
            fprintf(energyslice,"%d %d %d %15.15Lf \n", N/2, i, j, energyDen(0,fld,N/2,i,j)+energyDen2(0,fld,N/2,i,j));
            
            for(int fld=0;fld<nflds;fld++)
            {
                fprintf(slicefield,"%15.15Lf ", field[INDEX(0,fld,N/2,i,j)]);
                fprintf(sourceslice,"%15.15Lf ", source(0,fld,N/2,i,j)+source2(0,fld,N/2,i,j));
                
            }
            fprintf(slicefield, "\n");
            fprintf(sourceslice, "\n");
            
        }
        // fprintf(slicefield, "\n");
    }
    
#elif slice_orient==2//for xz-slice
    
    for(int i=0;i<N;i+=field_sliceskip)
    {
        for(int j=0;j<N;j+=field_sliceskip)
        {
            gravenergy1(i,N/2,j);
            gravenergy2(i,N/2,j);
            fprintf(gravrhoslice,"%d %d %d %15.15Lf %15.15Lf \n",i, N/2, j, gravrho1,gravrho2);
            fprintf(slicefield,"%d %d %d ",i, N/2, j);
            fprintf(sourceslice,"%d %d %d ",i, N/2, j);
            fprintf(energyslice,"%d %d %d %15.15Lf \n", i, N/2, j, energyDen(0,fld,i,N/2,j)+energyDen2(0,fld,i,N/2,j));
            
            for(int fld=0;fld<nflds;fld++)
            {
                fprintf(slicefield,"%15.15Lf ", field[INDEX(0,fld,i,N/2,j)]);
                fprintf(sourceslice,"%15.15Lf ", source(0,fld,i,N/2,j)+source2(0,fld,i,N/2,j));
                
            }
            fprintf(slicefield, "\n");
            fprintf(sourceslice, "\n");
            
        }
        // fprintf(slicefield, "\n");
    }
    
#endif
    
    
    
#elif field_outdim==1//outputs slice for 1dimension
    
    for(int k=0;k<N;k+=field_sliceskip)
    {
        //int fld=3;
        for(int fld=0;fld<nflds;fld++)
        {
            fprintf(slicefield,"%Le ", field[INDEX(0,fld,N/2,N/2,k)]);
            fprintf(sourceslice,"%Le ", source(0,fld,N/2,N/2,k)+source2(0,fld,N/2,N/2,k));
        }
        fprintf(slicefield, "\n");
        fprintf(sourceslice, "\n");
        
        
    }
    
#endif
    fclose(slicefield);
    fclose(energyslice);
     //fclose(gravrhoslice);
    fclose(sourceslice);
}


void meansvars()//calculates the mean and variance of each field
{
    char name_[5000];
    static FILE *meansvarsout_;
    int i, j, k, fld;
    gNum av,av_sq,var;
    
    sprintf(name_,"./slices/meansvariances.dat");
    meansvarsout_=fopen(name_,"a");
    
    fprintf(meansvarsout_,"%Lf",t);
    for(fld=0;fld<nflds;fld++)
    {
        av=0.;
        var=0.;
        // Calculate field mean
        LOOP
        av += field[INDEX(0,fld,i,j,k)];
        av = av/(gNum)gridsize; // Convert sum to average
        av_sq = av*av; // Calculate mean squared for finding variance
        // Calculate variance
        LOOP
        var += field[INDEX(0,fld,i,j,k)]*field[INDEX(0,fld,i,j,k)] - av_sq;
        var = var/(gNum)gridsize; // Convert sum to variance
        // Output means and variances to files
        fprintf(meansvarsout_," %Le %Le",av,var);
        // Check for instability. See if the field has grown exponentially and become non-numerical at any point.
        if((av!=0. && av/av!=1.)) // These two separate checks between them work on all the compilers I've tested
        {
            printf("Unstable solution developed. Field %d not numerical at t=%Le\n",fld,t);
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
    return (int) (endtime/dt/slicenumber);//calculates number of time-steps to wait for slicenumber slices by time endtime
}




void outputslice()//externally called function for outputting the data from the run
{
    static FILE  *slicetime;
    static int first=0;
    static int lastslice,lastspec;
    static int slicewaitm;
    
    if(first==0)
    {
        if(slicewait==0)//determines number of steps to wait between each slice output based off of g2parameters.h
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
        //meansvars();//outputs mean and variance of fields
#endif

#if pow_output!=0
        modePowerOut(t,first);
#endif     
        
        /*times file output */
        //this routine outputs time, scale-factor, scale-factor derivative, Hubble constant, and energy components at each slice output
        slicetime=fopen("./slices/slices_time.dat","a");
        fprintf(slicetime,"%Lf %Le\n", t, profile(t));
#if pi_powerout!=0
        gNum bkgf, piPow;
        piPow=piPowerOut(&bkgf);
        if((piPow!=0. && piPow/piPow!=1.)) // These two separate checks between them work on all the compilers I've tested
        {
            printf("Unstable solution developed. Pi not numerical at t=%Le\n",t);
            output_parameters();
            exit(1);
        }
        
        fprintf(slicetime,"%Lf %Le ",piPow, bkgf);//hij_powerout(),
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
        sprintf(expanType,"a double dot");
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
        fprintf(info,"Grid size=%ld^3\n",N);
        fprintf(info,"L=%Le\n",L);
        fprintf(info,"dt=%Le, dx=%Le\n",dt,dx);
        fprintf(info,"End time=%Le\n",endtime);
        //fprintf(info,"B = %Le\n",rescale_B);
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
        fprintf(info,"\nRun from t=%Le to t=%Le took ",starttime,t);
        readable_time((int)(tFinish-tStart),info);
        fprintf(info,"\n");
        fprintf(info, "\nFinal scale factor is %Le\n",a[0]);
    }
    
    fflush(info);
    return;
}

void nancheck()
{
    gNum av=0.;
    
    DECLARE_INDEX
    
    fldLOOP
        av+=field[INDEX(0,0,i,j,k)];
    
    if((av!=0. && av/av!=1.)) // These two separate checks between them work on all the compilers I've tested
    {
        printf("Unstable solution developed. Field %d not numerical at t=%Le\n",fld,t);
        output_parameters();
        exit(1);
    }
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
        printf("%Lf\n",t);
        time(&tLast);
    }
}
