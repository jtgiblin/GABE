/****************************
 PARAMETERS FILE
 ****************************/

/*
 This headr file contains the adjustable and non adjustable global parameters for the GABE2.0 program.
 
 Copyright (2013): 
 Kenyon College
 John T. Giblin, Jr
 Tate Deskins and Hillary Child
 Last Updated: 03.15.2024
*/

/***************
 model constants
 ***************/

#define fftw_flag 0  // 0 for double 1 for long double
#if fftw_flag == 1
    #define gNum long double
    #define exp expl
    #define sqrt sqrtl
    #define cbrt cbrtl
    #define log logl
    #define fftw_plan_dft_c2r_3d fftwl_plan_dft_c2r_3d
    #define fftw_plan_dft_r2c_3d fftwl_plan_dft_r2c_3d
    #define fftw_execute fftwl_execute
    #define fftw_plan fftwl_plan
    #define fftw_destroy_plan fftwl_destroy_plan
    #define fftw_complex fftwl_complex
    #define fftw_alloc_complex fftwl_alloc_complex
    #define fftw_alloc_real fftwl_alloc_real
    #define fftw_execute_dft_r2c fftwl_execute_dft_r2c
    #define fftw_execute_dft_c2r fftwl_execute_dft_c2r
    #define fftw_malloc fftwl_malloc
    #define fftw_free fftwl_free
    #define fftw_plan_with_nthreads fftwl_plan_with_nthreads
#else
    #define gNum double
#endif

#define num_flds 2// number of fields
const gNum mphi=1.e-6;//mass of phi field
const gNum phi0=0.193;//initial avg phi field value
const gNum gsq=2.5e-7;//g^2 value for phi chi coupling
const gNum f0[2]={phi0,0.};//array storing initial phi and chi field values
const gNum df0[2]={-0.142231,0.};//array storing initial phi and chi field derivative values

/***************************
 model independent parameters
 ***************************/
#define parallelize 1// for parallelization set to 1 and set other variables set to 0 for no parallelization
#define tot_num_thrds 4//total (max) number of threads to run during program
#defite calc_gws 1//0 for no gravitational waves, 1 for gravitational waves
const int randseed=44463132;//seed for rand number generator
const int N=128;//number of points along one side of grid
const gNum L=20.;// length of one side of box in prgm units
const gNum starttime=0.;//start time of simulation
const gNum endtime=300.;//end time of simulations
const gNum dt=L/(gNum)N/20.;//time step size
#define expansion_type 1//(0 for no expansion 1 for evolving from adot 2 for user defined expansion
//(will need to adjust functions file (adot and such) and type two evolution in the step() function and g2init.cpp initexpansion() for user defined expansion )

/*************************
 model dependent parameters
 *************************/

const gNum rescale_B=mphi;//rescaling

#define rand_init 1//1 to have random initialization 0 to not (see model file)
#define field_full_rand 1// 1 to have full random 0 to have symmetric kspace initialization
// #define spec_cut_off sqrt(3.)/8.// do not define for no cut off or smoothing
// #define spec_smooth 0.5// determines the sharpness of the tanh window function

/****************
 output parameters
 *****************/
const gNum screentime=60;// in seconds how frequently output prgm time to screen
const int slicewait=10;//how many dt's to wait between outputs (1 for no waiting) if 0 then slicenumber will be used.
const int slicenumber=20;//approx number of slices to output (only used if slicewait=0)
const int field_sliceskip=2;//how many points to print in field profile (1 is every, 2 every two, 3 every three...)
const int specnumber=1; //how many spectra to out put (1= every output slice 2 every two....)
#define field_outdim 2// number of dimensions of output in field profile (0 for no output)
#define spec_output 1// 1 to output spectra zero for no spectra output
#define var_output 1// 1 to output mean and variance zero for no variance output


/*********************************
 These are important DO NOT CHANGE
 *********************************/
const int nflds=num_flds; //stores number of fields for looping
const gNum dx=L/((gNum) N);//stores the change in x from point to point
const gNum gridsize=N*N*N;//stores size of grid for averaging


#if parallelize!=1
#undef tot_num_thrds
#define tot_num_thrds 1
#endif

