/***********************************************
Written by Xiaobo Meng
This program allows data shift when stacking and output 9times the MAD.
************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#define     max(a,b) ((a) > (b) ? (a):(b))
int     compare (const void *a, const void *b) {
        const float *da = (const float *) a;
        const float *db = (const float *) b;
        return  (*da > *db) - (*da < *db);
}


float comThreshold(float  *stack,int npts) 
{
    
    float *StackSort,median,MAD;
    int i=0;
    MAD=0.0;
    median=0.0;
    StackSort=(float*)malloc(npts*sizeof(float));
    memset(StackSort, 0x0, sizeof(float) * npts);
//    memcpy(StackSort,stack,sizeof(float) * npts);
    for(i=0;i<npts;i++){StackSort[i]=stack[i];}
    qsort(StackSort, npts, sizeof(float), compare);
    if ( npts%2 == 0 ) {
         median=(StackSort[npts/2]+StackSort[npts/2-1])/2;
    }

    if ( npts%2 == 1 ) {
         median=StackSort[(npts-1)/2];
    }    
   
    for (i=0; i<npts; i++) {
         StackSort[i]=fabsf(StackSort[i]-median);
    }
    
    qsort(StackSort, npts, sizeof(float), compare);   

    if ( npts%2 == 0 ) {
         MAD=(StackSort[npts/2]+StackSort[npts/2-1])/2;
    }

    if ( npts%2 == 1 ) {
         MAD=StackSort[(npts-1)/2];
    }
    float threshold=0.0;
    threshold=median+9*MAD;
    printf("%f %f %f\n",median,MAD,threshold);
    return threshold;
    free(StackSort);
    StackSort = NULL;
}
