/**********************************
 Initialization FILE
 **********************************/

/*
 This header file contains all the functions which are independent of the model needed to initialize the fields and random initial conditions.
 
 Copyright Kenyon College
 John T. Giblin, Jr
 Tate Deskins and Hillary Child
*/

#include "g2header.h" //contains declarations for program functions.

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

/** global parameters needed for fftw**/
#if rand_init==1

gNum *inic;
fftw_complex *fkf, *fkd;
fftw_plan picf, picd;


void dftMemAlloc()
{
    fftw_plan_with_nthreads(tot_num_thrds);
    
    inic = (gNum *) fftw_malloc(sizeof(gNum) * N * N * N);
    
    fkf   = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N * (N/2+1));
    fkd   = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N * (N/2+1));
    
    picf = fftw_plan_dft_c2r_3d(N, N, N, fkf, inic, FFTW_MEASURE);
    picd = fftw_plan_dft_c2r_3d(N, N, N, fkd, inic, FFTW_MEASURE);
}

/**function to initialize random initial conditions**/

void randInit(gNum f[][N][N],gNum df[][N][N],gNum d2vdf2)
{
    static int first=0;
    int i,j,k,tester;
    gNum i2,j2,k2, wk, randVar, randAngle;
    gNum dk2 = 4.*M_PI*M_PI/L/L;
    gNum dk = 2*M_PI/L;
    gNum fkRI1[2], fkRI2[2];
    
    
    gNum H0=sqrt(8.*M_PI/3.*edrho[0]);
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
                        
                        wk = sqrt(dk2*((gNum)i2 + (gNum)j2 + (gNum)k2) + d2vdf2);//dk*i\doti+d^2v/df^2
                                               
                       
                        
# ifdef spec_cut_off
                        
                        randVar = (gNum)rand()/(gNum)RAND_MAX;
                        randAngle = 2.*M_PI*(gNum)rand()/(gNum)RAND_MAX;
                        
                        fkRI1[0] = rescale_B*cos(randAngle)*sqrt(-log(randVar)/(2.*wk*L*L*L))*(1.-tanh(spec_smooth*(sqrt(i2+j2+k2)-N*spec_cut_off)))/2.;//gNum
                        fkRI1[1] = rescale_B*sin(randAngle)*sqrt(-log(randVar)/(2.*wk*L*L*L))*(1.-tanh(spec_smooth*(sqrt(i2+j2+k2)-N*spec_cut_off)))/2.;//imaginary
                        
                        randVar = (gNum)rand()/(gNum)RAND_MAX;//unlike lattice easy (becuse that was an error
                        randAngle = 2.*M_PI*(gNum)rand()/(gNum)RAND_MAX;

                        fkRI2[0] = rescale_B*cos(randAngle)*sqrt(-log(randVar)/(2.*wk*L*L*L))*(1.-tanh(spec_smooth*(sqrt(i2+j2+k2)-N*spec_cut_off)))/2.;//gNum
                        fkRI2[1] = rescale_B*sin(randAngle)*sqrt(-log(randVar)/(2.*wk*L*L*L))*(1.-tanh(spec_smooth*(sqrt(i2+j2+k2)-N*spec_cut_off)))/2.;//imaginary
                        
# else
                        
                        randVar = (gNum)rand()/(gNum)RAND_MAX;
                        randAngle = 2.*M_PI*(gNum)rand()/(gNum)RAND_MAX;
                        
                        fkRI1[0] = rescale_B*cos(randAngle)*sqrt(-log(randVar)/(2.*wk*L*L*L));//gNum
                        fkRI1[1] = rescale_B*sin(randAngle)*sqrt(-log(randVar)/(2.*wk*L*L*L));//imaginary
                        
                        randVar = (gNum)rand()/(gNum)RAND_MAX;//unlike lattice easy (becuse that was an error
                        randAngle = 2.*M_PI*(gNum)rand()/(gNum)RAND_MAX;

                        fkRI2[0] = rescale_B*cos(randAngle)*sqrt(-log(randVar)/(2.*wk*L*L*L));//gNum
                        fkRI2[1] = rescale_B*sin(randAngle)*sqrt(-log(randVar)/(2.*wk*L*L*L));//imaginary
# endif
                        

                        fkf[k + (N/2+1)*(j + N*i)][0] = (fkRI1[0] + fkRI2[0])/sqrt(2.);//gNum
                        fkf[k + (N/2+1)*(j + N*i)][1] = (fkRI1[1] + fkRI2[1])/sqrt(2.);//imaginary

                        fkd[k + (N/2+1)*(j + N*i)][0] = (fkRI2[1] - fkRI1[1])*wk/sqrt(2.) - H0*fkf[k + (N/2+1)*(j + N*i)][0];//gNum
                        fkd[k + (N/2+1)*(j + N*i)][1]= (fkRI1[0] - fkRI2[0])*wk/sqrt(2.) -   H0*fkf[k + (N/2+1)*(j + N*i)][1];//imaginary
#endif
#if field_full_rand==0
                        
                        wk = sqrt(dk2*((gNum)i2 + (gNum)j2 + (gNum)k2) + d2vdf2);//dk*i\doti+d^2v/df^2
                        
                        
                        fkf[k + (N/2+1)*(j + N*i)][0] = rescale_B*sqrt(1/(2.*wk*L*L*L))*(1.-tanh(spec_smooth*(sqrt(i2+j2+k2)-N*spec_cut_off)))/2.;//gNum
                        fkf[k + (N/2+1)*(j + N*i)][1] = 0.;//(fkRI1[1])*sqrt(2.);//imaginary
                        
                        
                        fkd[k + (N/2+1)*(j + N*i)][0] = - H0*fkf[k + (N/2+1)*(j + N*i)][0];//gNum
                        fkd[k + (N/2+1)*(j + N*i)][1]= 0.;//(fkRI1[0])*wk*sqrt(2.) -   H0*fkf[k + (N/2+1)*(j + N*i)][1];//imaginary
                        
#endif
                        
                        
                        break;
                }
            }
        }
    }
    
    fftw_execute(picf);
#pragma omp parallel for private (j,k) num_threads (tot_num_thrds)
    LOOP
    {
        f[i][j][k] += inic[k + N*(j + N*i)];
        
    }
    
    fftw_execute(picd);
#pragma omp parallel for private (j,k) num_threads (tot_num_thrds)
    LOOP
    {
        df[i][j][k] += inic[k + N*(j + N*i)];
    }
    
}


/**function to free up memory used for the dft**/

void initDestroy()
{
    fftw_destroy_plan(picf);
    fftw_destroy_plan(picd);
    fftw_free(inic);
    fftw_free(fkf);
    fftw_free(fkd);
}
#endif
