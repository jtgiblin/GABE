/*
 This file has the functions needed to evolve the metric perturbations.  It uses a 4th order runge-kutta integrator that ties into the GABE program.
 Written by John T. Giblin, jr .. Last modified 4.22.2007/7.7.2016
 */

#include "g2header.h"
#include<fftw3.h>

//const int Block_size = N/8;


//start by defining the same wrapping as in evolution.cpp
#if calc_gws == 1

// Increments a grid location accounting for periodic wrapping
inline int INCREMENT(int i)
{
    return( (i==N-1) ? 0 : i+1 );
}

// Decrements a grid location accounting for periodic wrapping
inline int DECREMENT(int i)
{
    return( (i==0) ? N-1 : i-1 );
}

// delfld calculates the components of the gradient of the field.
// fld corresponds to which field, i,j,k are placement on the lattice and
// n is the direction
// *** needs to be divided by dx to be a real gradient
inline gNum delfld(int fld, int s, int i, int j, int k, int n) {
    if (n==1) {
        return (gNum)(field[s][fld][INCREMENT(i)][j][k]-field[s][fld][DECREMENT(i)][j][k])/2;
    }
    else if(n==2) {
        return (gNum)(field[s][fld][i][INCREMENT(j)][k]-field[s][fld][i][DECREMENT(j)][k])/2;
    }
    else if(n==3) {
        return (gNum)(field[s][fld][i][j][INCREMENT(k)]-field[s][fld][i][j][DECREMENT(k)])/2;
    }
    else {
        return 0.;
    }
}


// jtg_stresstensor calculates the T_ab for any point (i,j,k). It's not really the SET, but appropriate combinations of the SET to get S_{ij}--the source term of the metric perturbations
inline gNum jtg_stresstensor(int aa, int bb, int s, int i, int j, int k) {
    
    gNum tempfld[nflds];  //to temp. store the field values for calling in
    
    //the potential energy function
    gNum deln2[3] = {(gNum)(1/pw2(dx)), (gNum)(1/pw2(dx)), (gNum)(1/pw2(dx))};  //normalization for the gradient
    //since they appear twice, only
    //one factor is saved - but we need to add the factor for chi
    
    for(int dd=0; dd<nflds; dd++)
    {
        tempfld[dd] = field[s][dd][i][j][k];   //save field values in temp.
    }
    
    gNum tempout=0.0;           //define output variable
    
    
    if (aa==1) {
        if (bb==1) {
            for(int m=0; m < nflds; m++) {
                tempout += (2./3.)*deln2[m]*delfld(m,s,i,j,k,1)*delfld(m,s,i,j,k,1)
                - (1./3.)*deln2[m]*delfld(m,s,i,j,k,2)*delfld(m,s,i,j,k,2)
                - (1./3.)*deln2[m]*delfld(m,s,i,j,k,3)*delfld(m,s,i,j,k,3);
            }
            return tempout;
        }
        else if (bb==2) {
            for(int m=0; m < nflds; m++) {
                tempout += deln2[m]*delfld(m,s,i,j,k,1)*delfld(m,s,i,j,k,2);
            }
            return tempout;
        }
        else if (bb==3) {
            for(int m=0; m < nflds; m++) {
                tempout += (deln2[m]*delfld(m,s,i,j,k,1)*delfld(m,s,i,j,k,3));
            }
            return tempout;
        }
        else {
            printf("calling wrong SET %i and %i\n", aa, bb);
            return 0.0;
        }
    }
    else if (aa==2) {
        if (bb==2) {
            for(int m=0; m < nflds; m++) {
                tempout += -(1./3.)*deln2[m]*delfld(m,s,i,j,k,1)*delfld(m,s,i,j,k,1)
                + (2./3.)*deln2[m]*delfld(m,s,i,j,k,2)*delfld(m,s,i,j,k,2)
                - (1./3.)*deln2[m]*delfld(m,s,i,j,k,3)*delfld(m,s,i,j,k,3);
            }
            return tempout;
        }
        else if (bb==3) {
            for(int m=0; m < nflds; m++) {
                tempout += (deln2[m]*delfld(m,s,i,j,k,2)*delfld(m,s,i,j,k,3));
            }
            return tempout;
        }
        else {
            printf("calling wrong SET %i and %i\n", aa, bb);
            return 0.0;
        }
    }
    else if (aa==3) {
        if (bb==3) {
            for(int m=0; m < nflds; m++) {
                tempout += -(1./3.)*deln2[m]*delfld(m,s,i,j,k,1)*delfld(m,s,i,j,k,1)
                - (1./3.)*deln2[m]*delfld(m,s,i,j,k,2)*delfld(m,s,i,j,k,2)
                + (2./3.)*deln2[m]*delfld(m,s,i,j,k,3)*delfld(m,s,i,j,k,3);
            }
            return tempout;
        }
        else
            printf("calling wrong SET %i and %i\n", aa, bb);
        return 0.0;
    }
    else {
        printf("calling wrong SET %i and %i\n", aa, bb);
        return 0.0;
    }
    
    printf("calling wrong SET NOT DOING ANYTHING %i and %i\n", aa, bb);
    return 0.;
    
}

//put the sourceterm in momentum space
void fft_stresstensor(int s) {
    
    int pz, px, py;
    gNum dpx, dpz, dpy;
    gNum momentum2;
    
    //need to declare the objects to use FFTW
    
    gNum (*tij0)[N][N], (*tij1)[N][N], (*tij2)[N][N];
    gNum (*tij3)[N][N], (*tij4)[N][N], (*tij5)[N][N];
    
    fftw_complex (*tij0_out)[N][N/2+1], (*tij1_out)[N][N/2+1], (*tij2_out)[N][N/2+1];
    fftw_complex (*tij3_out)[N][N/2+1], (*tij4_out)[N][N/2+1], (*tij5_out)[N][N/2+1];
    
    tij0 = (gNum (*)[N][N]) malloc(sizeof(gNum)*N*N*N);
    tij1 = (gNum (*)[N][N]) malloc(sizeof(gNum)*N*N*N);
    tij2 = (gNum (*)[N][N]) malloc(sizeof(gNum)*N*N*N);
    tij3 = (gNum (*)[N][N]) malloc(sizeof(gNum)*N*N*N);
    tij4 = (gNum (*)[N][N]) malloc(sizeof(gNum)*N*N*N);
    tij5 = (gNum (*)[N][N]) malloc(sizeof(gNum)*N*N*N);
    
    tij0_out = (fftw_complex (*)[N][N/2+1]) malloc(sizeof(fftw_complex)*N*N*(N/2+1));
    tij1_out = (fftw_complex (*)[N][N/2+1]) malloc(sizeof(fftw_complex)*N*N*(N/2+1));
    tij2_out = (fftw_complex (*)[N][N/2+1]) malloc(sizeof(fftw_complex)*N*N*(N/2+1));
    tij3_out = (fftw_complex (*)[N][N/2+1]) malloc(sizeof(fftw_complex)*N*N*(N/2+1));
    tij4_out = (fftw_complex (*)[N][N/2+1]) malloc(sizeof(fftw_complex)*N*N*(N/2+1));
    tij5_out = (fftw_complex (*)[N][N/2+1]) malloc(sizeof(fftw_complex)*N*N*(N/2+1));
    
    int i,j,k;
    
#pragma omp parallel for private (j,k) num_threads (tot_num_thrds)
    for(i=0; i<N; i++)  // Fills the arrays with the SET
    {
        for (j=0; j<N; j++)            //  adds up diagonal components
        {
            for (k=0; k<N; k++)
            {
                tij0[i][j][k] = jtg_stresstensor(1,1,s,i,j,k);
                tij1[i][j][k] = jtg_stresstensor(1,2,s,i,j,k);
                tij2[i][j][k] = jtg_stresstensor(1,3,s,i,j,k);
                tij3[i][j][k] = jtg_stresstensor(2,2,s,i,j,k);
                tij4[i][j][k] = jtg_stresstensor(2,3,s,i,j,k);
                tij5[i][j][k] = jtg_stresstensor(3,3,s,i,j,k);
            }
        }
    }
    
    
    if (!plan_fft_gw){
        int fftw_init_threads(void);
        fftw_plan_with_nthreads(tot_num_thrds);//tells fftw that tot_num_thrds are available for use durring dft
        plan_fft_gw = fftw_plan_dft_r2c_3d(N, N, N, &tij0[0][0][0], &tij0_out[0][0][0], FFTW_MEASURE);
    }
    
    
    fftw_execute_dft_r2c(plan_fft_gw, &tij0[0][0][0], &tij0_out[0][0][0]);
    fftw_execute_dft_r2c(plan_fft_gw, &tij1[0][0][0], &tij1_out[0][0][0]);
    fftw_execute_dft_r2c(plan_fft_gw, &tij2[0][0][0], &tij2_out[0][0][0]);
    fftw_execute_dft_r2c(plan_fft_gw, &tij3[0][0][0], &tij3_out[0][0][0]);
    fftw_execute_dft_r2c(plan_fft_gw, &tij4[0][0][0], &tij4_out[0][0][0]);
    fftw_execute_dft_r2c(plan_fft_gw, &tij5[0][0][0], &tij5_out[0][0][0]);
    //calculates the FFT
    
    gNum normmm = pow(dx,3);//*pow(a[s],3);
    //coefficient in front of the SET in the equations of motion (without the
    //  8\p/3)
    
    
#pragma omp parallel for private (j,k,px,py,pz,momentum2,dpz,dpx,dpy) num_threads (tot_num_thrds)
    for(i=0; i<N; i++) {              // stores them in the correct
        px = (i<=N/2 ? i : i-N);
        for (j=0; j<N; j++) {         // arrays
            py = (j<=N/2 ? j : j-N);
            for (k=0; k<N/2+1; k++) {
                pz = k;
                momentum2=pw2((gNum)px)+pw2((gNum)py)+pw2((gNum)pz);
                dpz = ((gNum)pz)/sqrt(momentum2);
                dpx = ((gNum)px)/sqrt(momentum2);
                dpy = ((gNum)py)/sqrt(momentum2);
                
                //T_{11}
                T_gw[s][0][i][j][k] = .5*pw2(1-dpx*dpx)*normmm*(tij0_out[i][j][k][0])
                -(1-dpx*dpx)*dpx*dpy*normmm*(tij1_out[i][j][k][0])
                -(1-dpx*dpx)*dpx*dpz*normmm*(tij2_out[i][j][k][0])
                +pw2(dpx*dpy)*normmm*(tij3_out[i][j][k][0])
                -.5*(1-dpx*dpx)*(1-dpy*dpy)*normmm*(tij3_out[i][j][k][0])
                +(1+dpx*dpx)*dpy*dpz*normmm*(tij4_out[i][j][k][0])
                +pw2(dpx*dpz)*normmm*(tij5_out[i][j][k][0])
                -.5*(1-dpx*dpx)*(1-dpz*dpz)*normmm*(tij5_out[i][j][k][0]);
                //T_{11}(i)
                T_gw[s][1][i][j][k] = .5*pw2(1-dpx*dpx)*normmm*(tij0_out[i][j][k][1])
                -(1-dpx*dpx)*dpx*dpy*normmm*(tij1_out[i][j][k][1])
                -(1-dpx*dpx)*dpx*dpz*normmm*(tij2_out[i][j][k][1])
                +pw2(dpx*dpy)*normmm*(tij3_out[i][j][k][1])
                -.5*(1-dpx*dpx)*(1-dpy*dpy)*normmm*(tij3_out[i][j][k][1])
                +(1+dpx*dpx)*dpy*dpz*normmm*(tij4_out[i][j][k][1])
                +pw2(dpx*dpz)*normmm*(tij5_out[i][j][k][1])
                -.5*(1-dpx*dpx)*(1-dpz*dpz)*normmm*(tij5_out[i][j][k][1]);
                //T_{22}
                T_gw[s][6][i][j][k] = .5*pw2(1-dpy*dpy)*normmm*(tij3_out[i][j][k][0])
                -(1-dpy*dpy)*dpx*dpy*normmm*(tij1_out[i][j][k][0])
                -(1-dpy*dpy)*dpy*dpz*normmm*(tij4_out[i][j][k][0])
                +pw2(dpx*dpy)*normmm*(tij0_out[i][j][k][0])
                -.5*(1-dpx*dpx)*(1-dpy*dpy)*normmm*(tij0_out[i][j][k][0])
                +(1+dpy*dpy)*dpx*dpz*normmm*(tij2_out[i][j][k][0])
                +pw2(dpy*dpz)*normmm*(tij5_out[i][j][k][0])
                -.5*(1-dpy*dpy)*(1-dpz*dpz)*normmm*(tij5_out[i][j][k][0]);
                //T_{22}(i)
                T_gw[s][7][i][j][k] = .5*pw2(1-dpy*dpy)*normmm*(tij3_out[i][j][k][1])
                -(1-dpy*dpy)*dpx*dpy*normmm*(tij1_out[i][j][k][1])
                -(1-dpy*dpy)*dpy*dpz*normmm*(tij4_out[i][j][k][1])
                +pw2(dpx*dpy)*normmm*(tij0_out[i][j][k][1])
                -.5*(1-dpx*dpx)*(1-dpy*dpy)*normmm*(tij0_out[i][j][k][1])
                +(1+dpy*dpy)*dpx*dpz*normmm*(tij2_out[i][j][k][1])
                +pw2(dpy*dpz)*normmm*(tij5_out[i][j][k][1])
                -.5*(1-dpy*dpy)*(1-dpz*dpz)*normmm*(tij5_out[i][j][k][1]);
                //T_{33}
                T_gw[s][10][i][j][k] = .5*pw2(1-dpz*dpz)*normmm*(tij5_out[i][j][k][0])
                -(1-dpz*dpz)*dpx*dpz*normmm*(tij2_out[i][j][k][0])
                -(1-dpz*dpz)*dpy*dpz*normmm*(tij4_out[i][j][k][0])
                +pw2(dpx*dpz)*normmm*(tij0_out[i][j][k][0])
                -.5*(1-dpx*dpx)*(1-dpz*dpz)*normmm*(tij0_out[i][j][k][0])
                +(1+dpz*dpz)*dpx*dpy*normmm*(tij1_out[i][j][k][0])
                +pw2(dpy*dpz)*normmm*(tij3_out[i][j][k][0])
                -.5*(1-dpy*dpy)*(1-dpz*dpz)*normmm*(tij3_out[i][j][k][0]);
                //T_{33}(i)
                T_gw[s][11][i][j][k] = .5*pw2(1-dpz*dpz)*normmm*(tij5_out[i][j][k][1])
                -(1-dpz*dpz)*dpx*dpz*normmm*(tij2_out[i][j][k][1])
                -(1-dpz*dpz)*dpy*dpz*normmm*(tij4_out[i][j][k][1])
                +pw2(dpx*dpz)*normmm*(tij0_out[i][j][k][1])
                -.5*(1-dpx*dpx)*(1-dpz*dpz)*normmm*(tij0_out[i][j][k][1])
                +(1+dpz*dpz)*dpx*dpy*normmm*(tij1_out[i][j][k][1])
                +pw2(dpy*dpz)*normmm*(tij3_out[i][j][k][1])
                -.5*(1-dpy*dpy)*(1-dpz*dpz)*normmm*(tij3_out[i][j][k][1]);
                //T_{12}
                T_gw[s][2][i][j][k] = -.5*dpx*dpy*(1-dpx*dpx)*normmm*(tij0_out[i][j][k][0])
                -.5*dpx*dpy*(1-dpy*dpy)*normmm*(tij3_out[i][j][k][0])
                +.5*dpx*dpy*(1+dpz*dpz)*normmm*(tij5_out[i][j][k][0])
                +(1-dpx*dpx)*(1-dpy*dpy)*normmm*(tij1_out[i][j][k][0])
                -dpy*dpz*(1-dpx*dpx)*normmm*(tij2_out[i][j][k][0])
                -dpx*dpz*(1-dpy*dpy)*normmm*(tij4_out[i][j][k][0]);
                //T_{12}(i)
                T_gw[s][3][i][j][k] = -.5*dpx*dpy*(1-dpx*dpx)*normmm*(tij0_out[i][j][k][1])
                -.5*dpx*dpy*(1-dpy*dpy)*normmm*(tij3_out[i][j][k][1])
                +.5*dpx*dpy*(1+dpz*dpz)*normmm*(tij5_out[i][j][k][1])
                +(1-dpx*dpx)*(1-dpy*dpy)*normmm*(tij1_out[i][j][k][1])
                -dpy*dpz*(1-dpx*dpx)*normmm*(tij2_out[i][j][k][1])
                -dpx*dpz*(1-dpy*dpy)*normmm*(tij4_out[i][j][k][1]);
                //T_{13}
                T_gw[s][4][i][j][k] = -.5*dpx*dpz*(1-dpx*dpx)*normmm*(tij0_out[i][j][k][0])
                +.5*dpx*dpz*(1+dpy*dpy)*normmm*(tij3_out[i][j][k][0])
                -.5*dpx*dpz*(1-dpz*dpz)*normmm*(tij5_out[i][j][k][0])
                +(1-dpx*dpx)*(1-dpz*dpz)*normmm*(tij2_out[i][j][k][0])
                -dpy*dpz*(1-dpx*dpx)*normmm*(tij1_out[i][j][k][0])
                -dpx*dpy*(1-dpz*dpz)*normmm*(tij4_out[i][j][k][0]);
                //T_{13}
                T_gw[s][5][i][j][k] = -.5*dpx*dpz*(1-dpx*dpx)*normmm*(tij0_out[i][j][k][1])
                +.5*dpx*dpz*(1+dpy*dpy)*normmm*(tij3_out[i][j][k][1])
                -.5*dpx*dpz*(1-dpz*dpz)*normmm*(tij5_out[i][j][k][1])
                +(1-dpx*dpx)*(1-dpz*dpz)*normmm*(tij2_out[i][j][k][1])
                -dpy*dpz*(1-dpx*dpx)*normmm*(tij1_out[i][j][k][1])
                -dpx*dpy*(1-dpz*dpz)*normmm*(tij4_out[i][j][k][1]);
                //T_{23}
                T_gw[s][8][i][j][k] = +.5*dpy*dpz*(1+dpx*dpx)*normmm*(tij0_out[i][j][k][0])
                -.5*dpy*dpz*(1-dpy*dpy)*normmm*(tij3_out[i][j][k][0])
                -.5*dpy*dpz*(1-dpz*dpz)*normmm*(tij5_out[i][j][k][0])
                +(1-dpy*dpy)*(1-dpz*dpz)*normmm*(tij4_out[i][j][k][0])
                -dpx*dpy*(1-dpz*dpz)*normmm*(tij2_out[i][j][k][0])
                -dpx*dpz*(1-dpy*dpy)*normmm*(tij1_out[i][j][k][0]);
                //T_{23}(i)
                T_gw[s][9][i][j][k] = +.5*dpy*dpz*(1+dpx*dpx)*normmm*(tij0_out[i][j][k][1])
                -.5*dpy*dpz*(1-dpy*dpy)*normmm*(tij3_out[i][j][k][1])
                -.5*dpy*dpz*(1-dpz*dpz)*normmm*(tij5_out[i][j][k][1])
                +(1-dpy*dpy)*(1-dpz*dpz)*normmm*(tij4_out[i][j][k][1])
                -dpx*dpy*(1-dpz*dpz)*normmm*(tij2_out[i][j][k][1])
                -dpx*dpz*(1-dpy*dpy)*normmm*(tij1_out[i][j][k][1]);
            }
        }
    }
    
    
    
    
    //subtract zero mode (just in case)
    for(int i=0; i<12; i++){
        T_gw[s][i][0][0][0] = 0;
    }
    
    //    FILTER(kmin_gw, kmax_gw, T_gw);
    // allows for a filter on the set
    
    free(tij0);
    free(tij1);
    free(tij2);
    free(tij3);
    free(tij4);
    free(tij5);
    
    free(tij0_out);
    free(tij1_out);
    free(tij2_out);
    free(tij3_out);
    free(tij4_out);
    free(tij5_out);
}


//  this updates the metric perturbations
void evolve_perts(int s, int snew, gNum deet) {
    
    fft_stresstensor(s);   //calculate the source terms in momentum space
    
    //printf("called fft_stresstensor \n");
    
    int i,j,k;
    
#pragma omp parallel for private (j,k) num_threads (tot_num_thrds)
    for(i=0;i<N;i++){        // EOM's
        int ii=(i<=N/2 ? i : i-N);
        for(j=0;j<N;j++){
            int jj = (j<=N/2 ? j : j-N);
            for(k=0;k<N/2+1;k++){ //only need to do half the array
                
                gNum source[12];
                
                source[0] = T_gw[s][0][i][j][k];      //take the correct parts
                source[1] = T_gw[s][1][i][j][k];
                source[2] = T_gw[s][2][i][j][k];
                source[3] = T_gw[s][3][i][j][k];
                source[4] = T_gw[s][4][i][j][k];
                source[5] = T_gw[s][5][i][j][k];
                source[6] = T_gw[s][6][i][j][k];
                source[7] = T_gw[s][7][i][j][k];
                source[8] = T_gw[s][8][i][j][k];
                source[9] = T_gw[s][9][i][j][k];
                source[10] = T_gw[s][10][i][j][k];
                source[11] = T_gw[s][11][i][j][k];
                
                
                gNum hd1[12];
                gNum ld1[12];
                
                derivs(h[s][0][i][j][k],
                       h[s][1][i][j][k],
                       h[s][2][i][j][k],
                       h[s][3][i][j][k],
                       h[s][4][i][j][k],
                       h[s][5][i][j][k],
                       h[s][6][i][j][k],
                       h[s][7][i][j][k],
                       h[s][8][i][j][k],
                       h[s][9][i][j][k],
                       h[s][10][i][j][k],
                       h[s][11][i][j][k],
                       l[s][0][i][j][k],
                       l[s][1][i][j][k],
                       l[s][2][i][j][k],
                       l[s][3][i][j][k],
                       l[s][4][i][j][k],
                       l[s][5][i][j][k],
                       l[s][6][i][j][k],
                       l[s][7][i][j][k],
                       l[s][8][i][j][k],
                       l[s][9][i][j][k],
                       l[s][10][i][j][k],
                       l[s][11][i][j][k],
                       ii, jj, k,
                       hd1, ld1,
                       source, s);
                
                
                h[snew][0][i][j][k] = h[0][0][i][j][k] + deet*hd1[0];
                h[snew][1][i][j][k] = h[0][1][i][j][k] + deet*hd1[1];
                h[snew][2][i][j][k] = h[0][2][i][j][k] + deet*hd1[2];
                h[snew][3][i][j][k] = h[0][3][i][j][k] + deet*hd1[3];
                h[snew][4][i][j][k] = h[0][4][i][j][k] + deet*hd1[4];
                h[snew][5][i][j][k] = h[0][5][i][j][k] + deet*hd1[5];
                h[snew][6][i][j][k] = h[0][6][i][j][k] + deet*hd1[6];
                h[snew][7][i][j][k] = h[0][7][i][j][k] + deet*hd1[7];
                h[snew][8][i][j][k] = h[0][8][i][j][k] + deet*hd1[8];
                h[snew][9][i][j][k] = h[0][9][i][j][k] + deet*hd1[9];
                h[snew][10][i][j][k] = h[0][10][i][j][k] + deet*hd1[10];
                h[snew][11][i][j][k] = h[0][11][i][j][k] + deet*hd1[11];
                
                
                l[snew][0][i][j][k] = l[0][0][i][j][k] + deet*ld1[0];
                l[snew][1][i][j][k] = l[0][1][i][j][k] + deet*ld1[1];
                l[snew][2][i][j][k] = l[0][2][i][j][k] + deet*ld1[2];
                l[snew][3][i][j][k] = l[0][3][i][j][k] + deet*ld1[3];
                l[snew][4][i][j][k] = l[0][4][i][j][k] + deet*ld1[4];
                l[snew][5][i][j][k] = l[0][5][i][j][k] + deet*ld1[5];
                l[snew][6][i][j][k] = l[0][6][i][j][k] + deet*ld1[6];
                l[snew][7][i][j][k] = l[0][7][i][j][k] + deet*ld1[7];
                l[snew][8][i][j][k] = l[0][8][i][j][k] + deet*ld1[8];
                l[snew][9][i][j][k] = l[0][9][i][j][k] + deet*ld1[9];
                l[snew][10][i][j][k] = l[0][10][i][j][k] + deet*ld1[10];
                l[snew][11][i][j][k] = l[0][11][i][j][k] + deet*ld1[11];
                
                
            }
            
        }
        
    }
    
    
    
    
    
    // file for testing in case you want to output any variables--these are some i found useful
    /*
    FILE *tomout;
    
    tomout = fopen("tomout.txt", "a");
    
    fprintf(tomout, "%10.10lf   %10.10lf \n", t, h[0][0][1][5][6]);
    
    
    fclose (tomout);
    */
}


//function for calculating the derivatives
void derivs(gNum h11,
            gNum hi11,
            gNum h12,
            gNum hi12,
            gNum h13,
            gNum hi13,
            gNum h22,
            gNum hi22,
            gNum h23,
            gNum hi23,
            gNum h33,
            gNum hi33,
            gNum l11,
            gNum li11,
            gNum l12,
            gNum li12,
            gNum l13,
            gNum li13,
            gNum l22,
            gNum li22,
            gNum l23,
            gNum li23,
            gNum l33,
            gNum li33,	    int ii, int jj, int k,
            gNum hd[12], gNum ld[12],
            gNum source_gw[12], int s) {
    
    gNum norm = (4.*M_PI*M_PI/L/L/a[s]/a[s]);//this includes the 1/a^2 in the EOM
    gNum omega = (8.*M_PI)/a[s]/a[s];
    gNum kk = (pw2((gNum)ii) + pw2((gNum)jj) + pw2((gNum)k))*norm;
    
    ld[0] = - kk*h11 - 3.*adot[s]/a[s]*l11 + 2.*omega*source_gw[0];
    ld[1] = - kk*hi11 - 3.*adot[s]/a[s]*li11 + 2.*omega*source_gw[1];
    ld[2] = - kk*h12 - 3.*adot[s]/a[s]*l12 + 2.*omega*source_gw[2];
    ld[3] = - kk*hi12 - 3.*adot[s]/a[s]*li12 + 2.*omega*source_gw[3];
    ld[4] = - kk*h13 - 3.*adot[s]/a[s]*l13 + 2.*omega*source_gw[4];
    ld[5] = - kk*hi13 - 3.*adot[s]/a[s]*li13 + 2.*omega*source_gw[5];
    ld[6] = - kk*h22 - 3.*adot[s]/a[s]*l22 + 2.*omega*source_gw[6];
    ld[7] = - kk*hi22 - 3.*adot[s]/a[s]*li22 + 2.*omega*source_gw[7];
    ld[8] = - kk*h23 - 3.*adot[s]/a[s]*l23 + 2.*omega*source_gw[8];
    ld[9] = - kk*hi23 - 3.*adot[s]/a[s]*li23 + 2.*omega*source_gw[9];
    ld[10] = - kk*h33 - 3.*adot[s]/a[s]*l33 + 2.*omega*source_gw[10];
    ld[11] = - kk*hi33 - 3.*adot[s]/a[s]*li33 + 2.*omega*source_gw[11];
    
    hd[0]=l11;
    hd[1]=li11;
    hd[2]=l12;
    hd[3]=li12;
    hd[4]=l13;
    hd[5]=li13;
    hd[6]=l22;
    hd[7]=li22;
    hd[8]=l23;
    hd[9]=li23;
    hd[10]=l33;
    hd[11]=li33;
    
}

#endif
