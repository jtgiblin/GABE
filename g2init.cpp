/**********************************
 Initialization FILE
 **********************************/

/*
 This header file contains all the functions which are independent of the model needed to initialize the fields and random initial conditions.
 
 Copyright Kenyon College
 John T. Giblin, Jr
 Tate Deskins and Hillary Child
*/

#include "g2header.h" //contains declerations for program functions.

/** expansion initialization **/
void initexpansion()
{
#if expansion_type==0
    //for power law
    a[0]=1;
    adot[0]=0;
    //addot[1]=0;
    a[1]=1;
    adot[1]=0;
    addot[1]=0;
    calcEnergy(0);
#elif expansion_type==1
    
    a[0]=1;
    calcEnergy(0);
    edrhot=edrho[0];
    adot[0]=adf(0);
    //addot[0]=0;
    //addot[1]=0;
#elif expansion_type==2
    /*
     a[0]=1;
     edkin[0]=avgKin(0);
     edpot[0]=avgPot(0);
     edgrad[0]=avgGrad(0);
     edrho[0]=edkin[0]+edpot[0]+edgrad[0];
     adot[0]=adf(0);
     addot[0]=addf(0);
     */
#endif
}

/** global parameters needed for fftwl**/
#if rand_init==1

long double *inic;
fftwl_complex *fkf, *fkd;
fftwl_plan picf, picd;


void dftMemAlloc()
{
    fftwl_plan_with_nthreads(tot_num_thrds);
    
    inic = (long double *) fftwl_malloc(sizeof(long double) * N * N * N);
    
    fkf   = (fftwl_complex*) fftwl_malloc(sizeof(fftwl_complex) * N * N * (N/2+1));
    fkd   = (fftwl_complex*) fftwl_malloc(sizeof(fftwl_complex) * N * N * (N/2+1));
    
    picf = fftwl_plan_dft_c2r_3d(N, N, N, fkf, inic, FFTW_MEASURE);
    picd = fftwl_plan_dft_c2r_3d(N, N, N, fkd, inic, FFTW_MEASURE);
}

/**function to initialize random initial conditions**/

void randInit(long double f[][N][N],long double df[][N][N],long double d2vdf2)
{
    static int first=0;
    int i,j,k,tester;
    long double i2,j2,k2, wk, randVar, randAngle;
    long double dk2 = 4.*M_PI*M_PI/L/L;
    long double dk = 2*M_PI/L;
    long double fkRI1[2], fkRI2[2];
    
    
    long double H0=sqrt(8.*M_PI*grav/3.*edrho[0]);
    //printf("%Lf\n",H0);
    if(first==0)
    {
        dftMemAlloc();
        first++;
        srand(randseed);
    }
    
    /* first we calculate phi(0, x) and phidot(0, x)	*/
#pragma omp parallel for private (j,k,i2,j2,k2,tester,wk,randVar,randAngle, fkRI1,fkRI2) num_threads (tot_num_thrds)
    for(i=0; i<N; i++){
        i2 = (i<(N/2+1) ? (i*i) : ((N-i)*(N-i)));
        
        for(j=0; j<N; j++){
            j2 = (j<(N/2+1) ? (j*j) : ((N-j)*(N-j)));
            
            for(k=0; k<N/2+1; k++){
                k2 = (k*k);
                tester=i+j+k;
                switch (tester) {//this initializes the zero mode to 0
                    case 0:
                        
                        fkf[k + (N/2+1)*(j + N*i)][0] = 0.;
                        fkf[k + (N/2+1)*(j + N*i)][1] = 0.;
                        
                        
                        fkd[k + (N/2+1)*(j + N*i)][0] = 0.;
                        fkd[k + (N/2+1)*(j + N*i)][1]= 0.;
                        break;
                        
                    default:
                        
#if field_full_rand==1
                        
                        wk = sqrt(dk2*((long double)i2 + (long double)j2 + (long double)k2) + d2vdf2);//dk*i\doti+d^2v/df^2
                                               
                       
                        
# ifdef spec_cut_off
                        
                        randVar = (long double)rand()/(long double)RAND_MAX;
                        randAngle = 2.*M_PI*(long double)rand()/(long double)RAND_MAX;
                        
                        fkRI1[0] = rescale_B*cos(randAngle)*sqrt(-log(randVar)/(2.*wk*L*L*L))*(1.-tanh(spec_smooth*(sqrt(i2+j2+k2)-N*spec_cut_off)))/2.;//long double
                        fkRI1[1] = rescale_B*sin(randAngle)*sqrt(-log(randVar)/(2.*wk*L*L*L))*(1.-tanh(spec_smooth*(sqrt(i2+j2+k2)-N*spec_cut_off)))/2.;//imaginary
                        
                        randVar = (long double)rand()/(long double)RAND_MAX;//unlike lattice easy (becuse that was an error
                        randAngle = 2.*M_PI*(long double)rand()/(long double)RAND_MAX;

                        fkRI2[0] = rescale_B*cos(randAngle)*sqrt(-log(randVar)/(2.*wk*L*L*L))*(1.-tanh(spec_smooth*(sqrt(i2+j2+k2)-N*spec_cut_off)))/2.;//long double
                        fkRI2[1] = rescale_B*sin(randAngle)*sqrt(-log(randVar)/(2.*wk*L*L*L))*(1.-tanh(spec_smooth*(sqrt(i2+j2+k2)-N*spec_cut_off)))/2.;//imaginary
                        
# else
                        
                        randVar = (long double)rand()/(long double)RAND_MAX;
                        randAngle = 2.*M_PI*(long double)rand()/(long double)RAND_MAX;
                        
                        fkRI1[0] = rescale_B*cos(randAngle)*sqrt(-log(randVar)/(2.*wk*L*L*L));//long double
                        fkRI1[1] = rescale_B*sin(randAngle)*sqrt(-log(randVar)/(2.*wk*L*L*L));//imaginary
                        
                        randVar = (long double)rand()/(long double)RAND_MAX;//unlike lattice easy (becuse that was an error
                        randAngle = 2.*M_PI*(long double)rand()/(long double)RAND_MAX;

                        fkRI2[0] = rescale_B*cos(randAngle)*sqrt(-log(randVar)/(2.*wk*L*L*L));//long double
                        fkRI2[1] = rescale_B*sin(randAngle)*sqrt(-log(randVar)/(2.*wk*L*L*L));//imaginary
# endif
                        

                        fkf[k + (N/2+1)*(j + N*i)][0] = (fkRI1[0] + fkRI2[0])/sqrt(2.);//long double
                        fkf[k + (N/2+1)*(j + N*i)][1] = (fkRI1[1] + fkRI2[1])/sqrt(2.);//imaginary

                        fkd[k + (N/2+1)*(j + N*i)][0] = (fkRI2[1] - fkRI1[1])*wk/sqrt(2.) - H0*fkf[k + (N/2+1)*(j + N*i)][0];//long double
                        fkd[k + (N/2+1)*(j + N*i)][1]= (fkRI1[0] - fkRI2[0])*wk/sqrt(2.) -   H0*fkf[k + (N/2+1)*(j + N*i)][1];//imaginary
#endif
#if field_full_rand==0
                        
                        wk = sqrt(dk2*((long double)i2 + (long double)j2 + (long double)k2) + d2vdf2);//dk*i\doti+d^2v/df^2
                        
                        
                        fkf[k + (N/2+1)*(j + N*i)][0] = rescale_B*sqrt(1/(2.*wk*L*L*L))*(1.-tanh(spec_smooth*(sqrt(i2+j2+k2)-N*spec_cut_off)))/2.;//long double
                        fkf[k + (N/2+1)*(j + N*i)][1] = 0.;//(fkRI1[1])*sqrt(2.);//imaginary
                        
                        
                        fkd[k + (N/2+1)*(j + N*i)][0] = - H0*fkf[k + (N/2+1)*(j + N*i)][0];//long double
                        fkd[k + (N/2+1)*(j + N*i)][1]= 0.;//(fkRI1[0])*wk*sqrt(2.) -   H0*fkf[k + (N/2+1)*(j + N*i)][1];//imaginary
                        
#endif
                        
                        
                        break;
                }
            }
        }
    }
    
    fftwl_execute(picf);
#pragma omp parallel for private (j,k) num_threads (tot_num_thrds)
    LOOP
    {
        f[i][j][k] += inic[k + N*(j + N*i)];
        
    }
    
    fftwl_execute(picd);
#pragma omp parallel for private (j,k) num_threads (tot_num_thrds)
    LOOP
    {
        df[i][j][k] += inic[k + N*(j + N*i)];
    }
    
}


/**function to free up memory used for the dft**/

void initDestroy()
{
    fftwl_destroy_plan(picf);
    fftwl_destroy_plan(picd);
    fftwl_free(inic);
    fftwl_free(fkf);
    fftwl_free(fkd);
}
#endif
