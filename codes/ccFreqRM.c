#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include "fftw3.h"

/* return the cross correlation of array x (length nx) and array y (length ny) 
 * the length of the cross correlation is nx + ny */
/* Note that the cross-correlation is computed in frequency domain */
/* Typically, x is the template event, and nx << ny */
double* ccFreq_( double *x, int nx,  double *y, int ny){
    int i, n, nHalf;
    double *tx, *ty, *xfft;
    double xNorm=0.0, *yWinNorm, *yWinMean;
    fftw_complex *xWkspace = NULL;
    fftw_complex *yWkspace = NULL;
    fftw_plan pForward, pInverse;

    n = nx + ny;
    tx = (double *)malloc(sizeof(double)*n);
    ty = (double *)malloc(sizeof(double)*n);
    memset(tx, 0, n*sizeof(double)); /* set tx to zero */
    memset(ty, 0, n*sizeof(double)); /* set ty to zero */
    for(i=0;i<nx;i++) {
	tx[i] = x[i];
	xNorm += x[i]*x[i];
    }//for
    xNorm = sqrt( xNorm  );
    for(i=0;i<ny;i++) ty[i] = y[i];

    nHalf = n/2 + 1;
    xWkspace = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nHalf);
    yWkspace = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nHalf);

    xfft = (double*) malloc(sizeof(double) * n);
    memset(xfft, 0, n*sizeof(double)); /* set xfft to zero */
    pForward = fftw_plan_dft_r2c_1d(n, xfft, xWkspace, FFTW_ESTIMATE);
    pInverse = fftw_plan_dft_c2r_1d(n, xWkspace, xfft, FFTW_ESTIMATE);

    /* sliding window norm of array y*/
    for (i=0; i<n; i++) xfft[i] = ty[i]*ty[i];
    fftw_execute(pForward);
    for(i=0; i<nHalf; i++) yWkspace[i] = xWkspace[i];

    /* fft(win) */
    for(i=0; i<nx; i++) xfft[i] = 1.0;
    for(i=nx; i<n; i++) xfft[i] = 0.0;
    fftw_execute(pForward);

    for(i=0; i<nHalf; i++) xWkspace[i] = yWkspace[i] * conj(xWkspace[i]);
    fftw_execute(pInverse);
    yWinNorm = (double *)malloc(sizeof(double)*n);
    memset(yWinNorm, 0, n*sizeof(double)); /* set tx to zero */
    for(i=0; i<n; i++) yWinNorm[i] = fabs(xfft[i]/n);

    for (i=0; i<n; i++) xfft[i] = ty[i];
    fftw_execute(pForward);
    for(i=0; i<nHalf; i++) yWkspace[i] = xWkspace[i];
    for(i=0; i<nx; i++) xfft[i] = 1.0;
    for(i=nx; i<n; i++) xfft[i] = 0.0;
    fftw_execute(pForward);

    for(i=0; i<nHalf; i++) xWkspace[i] = yWkspace[i] * conj(xWkspace[i]);
    for(i=0; i<nHalf; i++) xWkspace[i] = xWkspace[i];
    fftw_execute(pInverse);
    yWinMean = (double *)malloc(sizeof(double)*n);
    memset(yWinMean, 0, n*sizeof(double)); /* set ty to zero */
    for(i=0; i<n; i++) yWinMean[i] = xfft[i]/n;
    for(i=0; i<n; i++) yWinNorm[i] = sqrt(yWinNorm[i] - (yWinMean[i]*yWinMean[i])/nx);
   // for(i=0; i<n; i++) if(yWinMean[i]>=0.000000001) printf("%d  %d %d  %f  %f \n",i,nx,n,yWinNorm[i],yWinMean[i]);

    /* cross-correlation */
    for (i=0; i<n; i++) xfft[i] = tx[i];
    fftw_execute(pForward);
    for (i=0; i<nHalf; i++) yWkspace[i] = xWkspace[i];

    for (i=0; i<n; i++) xfft[i] = ty[i];
    fftw_execute(pForward);

    for (i=0; i<nHalf; i++) xWkspace[i] = xWkspace[i] * conj(yWkspace[i]);
    fftw_execute(pInverse);

    //Norm the cross correlation
    for (i=0; i<nx; i++) tx[i] = xfft[ny+i]/(yWinNorm[ny+i]*xNorm)/n;
    for (i=nx; i<n; i++) tx[i] = xfft[i-nx]/(yWinNorm[i-nx]*xNorm)/n;
    for (i=0; i<n; i++) {
      if (isinf(tx[i])) tx[i] = 0.;
    }
    //Free memory
    fftw_destroy_plan(pForward);
    fftw_destroy_plan(pInverse);
    fftw_free(yWkspace);
    fftw_free(xWkspace);
    free(xfft); 
    free(ty);
    free(yWinNorm);
    free(yWinMean);

    return tx;
}//ccFreq_
