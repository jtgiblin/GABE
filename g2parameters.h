/****************************
 PARAMETERS FILE
 ****************************/

/*
 This headr file contains the adjustable and non adjustable global parameters for the GABE2.0 program.
 
 Copyright (2013): 
 Kenyon College
 John T. Giblin, Jr
 Tate Deskins and Hillary Child
 Last Updated: 06.27.2013
*/

/***************
 model constants
 ***************/
 #define num_flds 1// number of fields
//const gNum mphi=1.e-6;//mass of phi field
//const gNum phi0=0.20123;//initial avg phi field value
//const gNum gsq=2.5e-7;//g^2 value for phi chi coupling
const gNum f0[num_flds]={0.};//array storing intial h_{ij} field values
const gNum df0[num_flds]={0.};//array storing intial h_{ij} field derivative values

const gNum c=1.; //speed of light, shouldn't change, if you do, fix all equations
const gNum omega = 0.142; // Radial velocity of source
const gNum rstar =10.;//This is the size of the vainshtein radius (where rb=2r_o)
const gNum sigma2 = 0.2; // Value of sigma-squared*2 in the gaussian energy density. May not need this value if energy density is not gaussian or some other form.
const gNum alpha = 1.; // This is the proportionality constant relating the energy densities of the two masses in binary. alpha = roh_2 / roh_1
//R0 shall be assumed to stay constant

#define gallileon_order 3// 2 for linear 3 for cubic 4 for quartic eventually 5 for quintic
//the following assume A=rb/2
const gNum kappa = rstar*rstar*rstar*2.*sqrt((1.+alpha)/(M_PI*omega));; //used in the Galileon EOM
const gNum kappa2 = kappa*kappa; //Kappa-squared
const gNum lambda = 8.*sqrt(M_PI*omega/(1.+alpha));
/***************************
 model independent parameters
 ***************************/
#define parallelize 1// for parallelization set to 1 and set other variables set to 0 for no parallelization
#define tot_num_thrds 6//total (max) number of threads to run durring program
const int randseed=44463132;//seed for rand number generator
const gIdx N=256;//number of points along one side of grid
const gNum L=50.;// length of one side of box in prgm units
const gNum starttime=0.;//start time of simulation
const gNum endtime=100.;//end time of simulations
const gNum dt=0.002;//time step size
#define int_order 2//integration order (2 or 4)
#define expansion_type 0//(0 for no expansion 1 for evolving from adot 2 for user defined expansion 
//(will need to adjust functions file (adot and such) and type two evolution in the step() function fnd g2init.cpp initexpansion() for user defined expansion )

/*************************
 model dependent parameters
 *************************/

//const gNum rescale_B=lambda;//rescallings

#define rand_init 0//1 to have random initialization 0 to not (see model file)
#define field_full_rand 1// 1 to have full random 0 to have symmetric kspace initilaization
// #define spec_cut_off sqrt(3.)/8.// do not define for no cut off or smoothing
// #define spec_smooth 0.5// determines the sharpness of the tanh window function

/****************
 output parameters
 *****************/
const gNum screentime=60;// in seconds how frequently output prgm time to screen
const int slicewait=0;//how many dt's to wait between outputs (1 for no waiting) if 0 then slicenumber will be used.
const int slicenumber=300;//approx number of slices to output (only used if slicewait=0)
const int field_sliceskip=4;//howmany points to print in field profile (1 is every, 2 every two, 3 every three...)
const int specnumber=1; //how many spectra to out put (1= every output slice 2 every two....)
#define field_outdim 1// number of dimensions of output in field profile (0 for no output)
#define spec_output  0// 1 to ouptut spectra, 0 for no spectra output
#define var_output   0// 1 to output mean and varraince, 0 for no varriance output
#define slice_orient 0// 0 for xy-slice; 1 for yz-slice; 2 for xz-slice
#define nan_check 0//1 to have ind nan-check

/*********************************
 These are important DO NOT CHANGE
 *********************************/
const int nflds=num_flds; //stores number of fields for looping
const gNum dx=L/((gNum) N);//stores the change in x from point to point
const gNum gridsize=N*N*N;//stores size of grid for averaging


const gNum PNorm = /*dx*dx*dx*/ 1./sqrt(sigma2*sigma2*sigma2*M_PI*M_PI*M_PI); //Normalization constant for source terms

#if parallelize!=1
#undef tot_num_thrds
#define tot_num_thrds 1
#endif

